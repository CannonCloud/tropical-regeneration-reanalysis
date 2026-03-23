#!/usr/bin/env python3
"""
Step 2: Sample Environmental Layers
====================================

This script samples environmental variables at your master spine points.

Usage:
    # Run ALL 10 layers (default)
    python step2_sample_layers.py --samples 200000
    
    # Run just ONE layer
    python step2_sample_layers.py --samples 200000 --layer soil_biome
    python step2_sample_layers.py --samples 200000 --layer climate
"""

import ee
import argparse
import logging
import time
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    PROJECT_ID = '# replace with your GEE project ID'
    SCALE = 30
    OUTPUT_DIR = Path('./gee_outputs')
    LOG_FILE = OUTPUT_DIR / 'step2_sample_layers.log'
    
    def get_asset_id(self, n_samples):
        return f'projects/# replace with your GEE project ID/assets/master_sample_points_{n_samples}'

# All layer groups (soil split into 3 for lighter computation)
LAYER_GROUPS = [
    'climate',
    'topography', 
    'landcover',
    'npp_fire',
    'soil_chemical',
    'soil_texture',
    'soil_physical',
    'socioeconomic',
    'highres',
    'pasture',
    'opportunity_cost',
    'plantation_raw'
]

# =============================================================================
# LAYER DEFINITIONS
# =============================================================================

class LayerDefinitions:
    
    @staticmethod
    def get_climate():
        """Climate (6 bands): WorldClim bioclimatic variables"""
        clim = ee.Image("WORLDCLIM/V1/BIO")
        return clim.select(['bio01', 'bio04', 'bio06', 'bio12', 'bio14', 'bio15']) \
            .rename(['bio_temp_mean', 'bio_temp_seasonality', 'bio_min_temp_coldest', 
                     'bio_precip_annual', 'bio_precip_driest', 'bio_precip_seasonality']) \
            .unmask(0)
    
    @staticmethod
    def get_topography():
        """Topography (3 bands): Elevation, slope, forest density 2018"""
        topo = ee.Image("USGS/SRTMGL1_003")
        slope = ee.Terrain.slope(topo).rename('bio_slope').unmask(0)
        elev = topo.select('elevation').rename('bio_elevation').unmask(0)
        
        hansen = ee.Image("UMD/hansen/global_forest_change_2018_v1_6")
        tc2000 = hansen.select('treecover2000')
        lossyear = hansen.select('lossyear')
        gain = hansen.select('gain')
        
        cover_post_loss = tc2000.where(lossyear.gte(1).And(lossyear.lte(17)), 0)
        is_stable_gain = gain.eq(1).And(lossyear.eq(0).Or(lossyear.gt(12)))
        forest_2018 = cover_post_loss.where(is_stable_gain, 100).rename('bio_forest_density_2018').unmask(0)
        
        return elev.addBands(slope).addBands(forest_2018)
    
    @staticmethod
    def get_landcover():
        """Land Cover (3 bands): ESA CCI classification, cropland weights"""
        lc_raw = ee.Image("projects/# replace with your GEE project ID/assets/lc_300m_2015").select('b1')
        from_vals = [10, 11, 12, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220]
        to_vals = [1, 1, 1, 1, 3, 3, 4, 5, 6, 6, 6, 7, 7, 8, 9, 2, 10, 11, 11, 11, 10, 10, 10, 10]
        lc_11class = lc_raw.remap(from_vals, to_vals, 0).rename('bio_landcover_class').unmask(0)
        
        crop_weights = lc_raw.remap(
            [10, 11, 12, 20, 30, 40], 
            [1.0, 1.0, 1.0, 1.0, 0.75, 0.25], 
            0
        ).rename('cropland_weight').unmask(0)
        
        return lc_11class.addBands(crop_weights).addBands(lc_raw.rename('lc_raw_value'))
    
    @staticmethod
    def get_npp_fire():
        """NPP & Fire (2 bands): MODIS productivity and fire frequency
        
        bio_npp_mean: Mean annual NPP (kg C/m²/yr scaled), 2001-2016
        bio_fire_freq: Count of years with any burning, 2001-2016 (0-16 integer)
        
        NOTE: No .clip() used - filterBounds() is sufficient and .clip() on a
        global rectangle causes GEE masking issues that produce all-zero outputs.
        """
        tropics = ee.Geometry.Rectangle([-180, -25, 180, 25])
        
        # NPP: mean annual net primary productivity
        npp = ee.ImageCollection("MODIS/061/MOD17A3HGF") \
            .filterDate('2001-01-01', '2017-12-31') \
            .filterBounds(tropics) \
            .select('Npp') \
            .mean() \
            .rename('bio_npp_mean').unmask(0)
        
        # FIRE: count of years (2001-2016) with any detected burning
        # BurnDate encodes Julian day of burn (1-365), or 0=unburned, negative=water/missing
        # .gt(0) correctly identifies burned pixels; .max() per year = burned at least once
        # Result is a clean 0-16 integer (number of years with fire)
        def year_burned(year):
            year = ee.Number(year)
            start = ee.Date.fromYMD(year, 1, 1)
            end = ee.Date.fromYMD(year, 12, 31)
            return ee.ImageCollection("MODIS/061/MCD64A1") \
                .filterDate(start, end) \
                .filterBounds(tropics) \
                .select('BurnDate') \
                .map(lambda img: img.gt(0)) \
                .max() \
                .unmask(0)  # months with no data treated as unburned
        
        years = ee.List.sequence(2001, 2016)
        fire_freq = ee.ImageCollection(years.map(year_burned)) \
            .sum() \
            .rename('bio_fire_freq') \
            .unmask(0)
        
        return npp.addBands(fire_freq)
    
    @staticmethod
    def get_soil_chemical():
        """Soil Chemical & Biome (3 bands): pH, organic carbon, biome
        
        Returns:
        - bio_soil_organic_carbon: Organic carbon (g/kg), top 30cm weighted avg
        - bio_soil_ph: pH × 10, top 30cm weighted avg
        - bio_biome_id: Biome identifier (categorical)
        """
        
        # Biome
        biomes_fc = ee.FeatureCollection("RESOLVE/ECOREGIONS/2017")
        biomes = ee.Image().byte().paint(biomes_fc, 'BIOME_NUM').rename('bio_biome_id')
        
        # Helper: Calculate weighted average for top 30cm
        def top30cm_average(image_path, band_prefix='b'):
            img = ee.Image(image_path)
            avg = img.select(f'{band_prefix}0').multiply(0.167) \
                .add(img.select(f'{band_prefix}10').multiply(0.333)) \
                .add(img.select(f'{band_prefix}30').multiply(0.500))
            return avg
        
        # Chemical properties
        soil_oc = top30cm_average("OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02") \
            .rename('bio_soil_organic_carbon').unmask(0)
        
        soil_ph = top30cm_average("OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02") \
            .rename('bio_soil_ph').unmask(0)
        
        return soil_oc.addBands(soil_ph).addBands(biomes)
    
    @staticmethod
    def get_soil_texture():
        """Soil Texture (2 bands): Sand and clay content
        
        Returns:
        - bio_soil_sand: Sand content (g/kg), top 30cm weighted avg
        - bio_soil_clay: Clay content (g/kg), top 30cm weighted avg
        """
        
        # Helper: Calculate weighted average for top 30cm
        def top30cm_average(image_path, band_prefix='b'):
            img = ee.Image(image_path)
            avg = img.select(f'{band_prefix}0').multiply(0.167) \
                .add(img.select(f'{band_prefix}10').multiply(0.333)) \
                .add(img.select(f'{band_prefix}30').multiply(0.500))
            return avg
        
        # Texture variables
        soil_sand = top30cm_average("OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02") \
            .rename('bio_soil_sand').unmask(0)
        
        soil_clay = top30cm_average("OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02") \
            .rename('bio_soil_clay').unmask(0)
        
        return soil_sand.addBands(soil_clay)
    
    @staticmethod
    def get_soil_physical():
        """Soil Physical Properties (2 bands): Bulk density and water content
        
        Returns:
        - bio_soil_bulk_density_fine: Bulk density fine earth (kg/dm³), top 30cm weighted avg
        - bio_soil_water_33kpa: Water content at field capacity (vol %), top 30cm weighted avg
        """
        
        # Helper: Calculate weighted average for top 30cm
        def top30cm_average(image_path, band_prefix='b'):
            img = ee.Image(image_path)
            avg = img.select(f'{band_prefix}0').multiply(0.167) \
                .add(img.select(f'{band_prefix}10').multiply(0.333)) \
                .add(img.select(f'{band_prefix}30').multiply(0.500))
            return avg
        
        # Physical properties
        soil_bd_fine = top30cm_average("OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02") \
            .rename('bio_soil_bulk_density_fine').unmask(0)
        
        soil_water_33 = top30cm_average("OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01") \
            .rename('bio_soil_water_33kpa').unmask(0)
        
        return soil_bd_fine.addBands(soil_water_33)
    
    @staticmethod
    def get_socioeconomic():
        """Socioeconomic (5 bands): Population, GDP, roads, water, protected areas"""
        wdpa = ee.FeatureCollection("WCMC/WDPA/current/polygons")
        protected_mask = ee.Image(0).byte().paint(wdpa, 1).rename('paper_protected_binary')
        
        roads = ee.Image("projects/# replace with your GEE project ID/assets/grip4_density").rename('paper_roads_density')
        dist_water = ee.Image("projects/# replace with your GEE project ID/assets/distance2water_30arcsec").rename('paper_dist_water')
        
        gdp_img = ee.Image("projects/sat-io/open-datasets/GRIDDED_HDI_GDP/adm2_gdp_perCapita_1990_2022")
        gdp = gdp_img.select(['PPP_2015']).rename('paper_gdp_2015')
        
        ghs_pop = ee.ImageCollection("JRC/GHSL/P2023A/GHS_POP") \
            .filterDate('2015-01-01', '2015-12-31') \
            .first().select([0]).rename('paper_ghs_pop')
        
        return ghs_pop.addBands(gdp).addBands(roads).addBands(dist_water).addBands(protected_mask)
    
    @staticmethod
    def get_highres():
        """High-Resolution (4 bands): Nightlights, travel time, human modification, WorldPop"""
        lights = ee.ImageCollection("NOAA/VIIRS/DNB/ANNUAL_V21") \
            .filterDate('2014-01-01', '2016-12-31') \
            .select('median_masked') \
            .median() \
            .rename('highres_lights_2014_2016')
        
        travel = ee.Image("Oxford/MAP/accessibility_to_cities_2015_v1_0") \
            .select('accessibility').rename('highres_travel')
        
        ghm = ee.ImageCollection("CSP/HM/GlobalHumanModification") \
            .first().select('gHM').rename('highres_ghm')
        
        worldpop = ee.ImageCollection("WorldPop/GP/100m/pop") \
            .filterDate('2015-01-01', '2016-01-01') \
            .select('population') \
            .mosaic() \
            .rename('highres_worldpop_density') \
            .unmask(0)
        
        return worldpop.addBands(lights).addBands(travel).addBands(ghm)
    
    @staticmethod
    def get_pasture():
        """Pasture (1 band): Global Pasture Watch classifications"""
        gpw_class = ee.ImageCollection("projects/global-pasture-watch/assets/ggc-30m/v1/grassland_c") \
            .filterDate('2015-01-01', '2015-12-31') \
            .mosaic() \
            .unmask(0) \
            .rename('gpw_raw_class')
        
        return gpw_class
    
    @staticmethod
    def get_opportunity_cost():
        """Opportunity Cost (1 band): Agricultural value (SPAM)"""
        opp_cost = ee.Image("projects/# replace with your GEE project ID/assets/Opp_Cost_USD_per_Ha") \
            .rename('econ_opp_cost_usd_ha') \
            .unmask(0)
        
        return opp_cost
    
    @staticmethod
    def get_plantation_raw():
        """Plantation (1 band): Xiao et al. planted forests"""
        plantation_rgb = ee.ImageCollection("projects/sat-io/open-datasets/GLOBAL-NATURAL-PLANTED-FORESTS") \
            .mosaic() \
            .unmask(0)
        
        return plantation_rgb.select('b1').rename('plantation_raw_b1')

# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Sample Environmental Layers')
    parser.add_argument('--samples', type=int, required=True,
                        help='Sample size (must match asset, e.g., 200000)')
    parser.add_argument('--layer', type=str, default='all',
                        choices=['all'] + LAYER_GROUPS,
                        help='Layer to sample (default: all)')
    parser.add_argument('--delay', type=int, default=30,
                        help='Seconds between task submissions when running all (default: 30)')
    
    args = parser.parse_args()
    
    # Setup
    config = Config()
    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
        handlers=[
            logging.FileHandler(config.LOG_FILE),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)
    
    # Determine which layers to run
    if args.layer == 'all':
        layers_to_run = LAYER_GROUPS
        logger.info("="*60)
        logger.info(f"STEP 2: SAMPLE ALL {len(LAYER_GROUPS)} LAYERS")
        logger.info("="*60)
    else:
        layers_to_run = [args.layer]
        logger.info("="*60)
        logger.info(f"STEP 2: SAMPLE LAYER '{args.layer.upper()}'")
        logger.info("="*60)
    
    logger.info(f"Sample size: {args.samples:,}")
    
    # Initialize GEE
    try:
        ee.Initialize(project=config.PROJECT_ID)
        logger.info(f"\n✓ GEE initialized: {config.PROJECT_ID}")
    except:
        logger.info("\nAuthenticating...")
        ee.Authenticate()
        ee.Initialize(project=config.PROJECT_ID)
        logger.info("✓ Authentication successful")
    
    # Load master spine
    asset_id = config.get_asset_id(args.samples)
    logger.info(f"\nLoading master spine from asset...")
    logger.info(f"Asset: {asset_id}")
    
    try:
        master_spine = ee.FeatureCollection(asset_id)
        logger.info("✓ Master spine loaded")
    except Exception as e:
        logger.error(f"✗ Failed to load asset: {e}")
        exit(1)
    
    # Submit layer tasks
    logger.info(f"\n{'='*60}")
    logger.info(f"SUBMITTING {len(layers_to_run)} LAYER TASK(S)")
    logger.info(f"{'='*60}\n")
    
    submitted_tasks = []
    
    for i, layer_name in enumerate(layers_to_run, 1):
        logger.info(f"[{i}/{len(layers_to_run)}] Processing: {layer_name}")
        
        try:
            # Get layer function
            layer_function = getattr(LayerDefinitions, f'get_{layer_name}')
            layer_stack = layer_function()
            
            band_names = layer_stack.bandNames().getInfo()
            logger.info(f"  ✓ Built stack: {len(band_names)} bands")
            logger.info(f"    Bands: {', '.join(band_names)}")
            
            # Sample at master spine points
            sampled_data = layer_stack.sampleRegions(
                collection=master_spine,
                scale=config.SCALE,
                geometries=True
            )
            
            # Export
            description = f"PNV_Layer_{layer_name}_N{args.samples}"
            
            task = ee.batch.Export.table.toDrive(
                collection=sampled_data,
                description=description,
                fileFormat='CSV'
            )
            
            task.start()
            task_id = task.id
            
            submitted_tasks.append({
                'layer': layer_name,
                'task_id': task_id,
                'description': description
            })
            
            logger.info(f"  ✓ Task submitted: {task_id}")
            logger.info(f"  ✓ Will export: {description}.csv")
            
            # Wait between submissions (only if running multiple layers)
            if len(layers_to_run) > 1 and i < len(layers_to_run):
                logger.info(f"  Waiting {args.delay}s before next task...\n")
                time.sleep(args.delay)
            else:
                logger.info("")
            
        except Exception as e:
            logger.error(f"  ✗ Failed to submit {layer_name}: {e}\n")
            continue
    
    # Summary
    logger.info("="*60)
    logger.info("TASKS SUBMITTED")
    logger.info("="*60 + "\n")
    
    logger.info(f"Successfully submitted {len(submitted_tasks)}/{len(layers_to_run)} task(s):\n")
    for task_info in submitted_tasks:
        logger.info(f"  ✓ {task_info['layer']}: {task_info['task_id']}")
    
    logger.info(f"\n{'='*60}")
    logger.info("WHAT HAPPENS NEXT")
    logger.info(f"{'='*60}\n")
    
    logger.info("Tasks are running on Google's servers.")
    logger.info("Monitor: https://code.earthengine.google.com/tasks\n")
    
    if len(layers_to_run) == 1:
        logger.info("Estimated time: ~30-60 minutes")
    else:
        logger.info("Estimated time:")
        if args.samples <= 1000:
            logger.info("  ~1-2 hours for all layers")
        elif args.samples <= 50000:
            logger.info("  ~2-4 hours for all layers")
        else:
            logger.info("  ~4-6 hours for all layers")
    
    logger.info("\nWhen complete:")
    logger.info(f"1. Download CSV(s) from Google Drive")
    logger.info(f"2. Place in: {config.OUTPUT_DIR}/")
    logger.info("3. Rename to: layer_<name>.csv")
    
    if len(layers_to_run) == len(LAYER_GROUPS):
        logger.info("4. Run Step 3 (distance calculations)")
    
    logger.info(f"\n{'='*60}\n")

if __name__ == '__main__':
    main()
