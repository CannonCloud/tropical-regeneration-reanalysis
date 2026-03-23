#!/usr/bin/env python3
"""
Brazil Step 2: Sample Environmental Layers
===========================================

Samples environmental variables at Brazil regrowth or non-regrowth points.

Input assets (merged in GEE before running this step):
    Regrowth:     brazil_regrowth_points
    Non-regrowth: brazil_nonregrowth_points_Cleaned
    Annulus:      brazil_nonregrowth_annulus_batch_00_Cleaned

Key changes vs the global pipeline:
    - climate:      WorldClim V2 (v2.1), all 19 raw bioclim bands (PCA done in R)
    - landcover:    ESA CCI 2000 (not 2015) — matches Williams training-era data
    - soil_chemical: includes bio_biome_id (RESOLVE BIOME_NUM) as categorical
    - forest density and distance to forest are NOT extracted here —
      those require focal/distance-transform operations done in step3

Layer groups and the variables they produce:
    climate_part1   bio01–bio04
    climate_part2   bio05–bio08
    climate_part3   bio09–bio12
    climate_part4   bio13–bio16
    climate_part5   bio17–bio19
    topography      bio_elevation, bio_slope
    landcover       bio_landcover_class, cropland_weight, lc_raw_value
    npp_fire        bio_npp_mean, bio_fire_freq
    soil_chemical   bio_soil_organic_carbon, bio_soil_ph, bio_biome_id
    soil_texture    bio_soil_sand, bio_soil_clay
    soil_physical   bio_soil_bulk_density_fine, bio_soil_water_33kpa
    socioeconomic   econ_dist_water, econ_road_density, paper_gdp_2015,
                    paper_ghs_pop, paper_protected_binary
    highres_lights  highres_lights_2012_2014   (sampled alone — GHM has mask gaps)
    highres_travel  highres_travel
    highres_ghm     highres_ghm                (sat-io year 2000, not CSP)
    highres_worldpop highres_worldpop_density
    pasture         gpw_raw_class
    opportunity_cost econ_opp_cost_usd_ha

Each highres band is sampled as its own layer so that a mask gap in one
(e.g. GHM over unmodified forest) does not cause sampleRegions to drop
the entire feature, which would null out all four columns at that point.

Outputs go to Google Drive as CSVs, one per layer group.
Merge on point_id in R.

Usage:
    # Run all layers:
    python brazil_step2_sample_layers.py --sample_type regrowth
    python brazil_step2_sample_layers.py --sample_type nonregrowth

    # Run a specific layer only:
    python brazil_step2_sample_layers.py --sample_type regrowth --layer socioeconomic
    python brazil_step2_sample_layers.py --sample_type regrowth --layer highres_lights
    python brazil_step2_sample_layers.py --sample_type regrowth --layer highres_ghm
    python brazil_step2_sample_layers.py --sample_type regrowth --layer pasture
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
    SCALE      = 30
    OUTPUT_DIR = Path('./gee_outputs')
    LOG_FILE   = OUTPUT_DIR / 'brazil_step2_sample_layers.log'

    # ── Sample asset mapping ──────────────────────────────────────────────────
    SAMPLE_ASSETS = {
        'regrowth':         'projects/# replace with your GEE project ID/assets/brazil_regrowth_points',
        'nonregrowth':      'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_points_Cleaned',
        'nonregrowth_b02':  'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_batch_02',
        'nonregrowth_b03':  'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_batch_03',
        'nonregrowth_b04':  'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_batch_04',
        'nonregrowth_b05':  'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_batch_05',
        'annulus':          'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_annulus_batch_00_Cleaned',
    }


LAYER_GROUPS = [
    'climate_part1',   # bio01-bio04
    'climate_part2',   # bio05-bio08
    'climate_part3',   # bio09-bio12
    'climate_part4',   # bio13-bio16
    'climate_part5',   # bio17-bio19
    'topography',
    'landcover',
    'npp_fire',
    'soil_chemical',
    'soil_texture',
    'soil_physical',
    'socioeconomic',   # roads, dist_water, GDP, GHS pop, protected areas
    'highres_lights',  # VIIRS nightlights 2012-2014 median
    'highres_travel',  # travel time to cities 2015
    'highres_ghm',     # Global Human Modification year 2000
    'highres_worldpop',# WorldPop 100m population 2000
    'pasture',
    'opportunity_cost',
]


# =============================================================================
# LAYER DEFINITIONS
# =============================================================================

class LayerDefinitions:

    @staticmethod
    def get_climate_part1():
        """bio01-bio04"""
        clim = ee.Image("projects/# replace with your GEE project ID/assets/brazil_worldclim_v2")
        return clim.select(['b1','b2','b3','b4']).rename(['bio01','bio02','bio03','bio04']).unmask(0)

    @staticmethod
    def get_climate_part2():
        """bio05-bio08"""
        clim = ee.Image("projects/# replace with your GEE project ID/assets/brazil_worldclim_v2")
        return clim.select(['b5','b6','b7','b8']).rename(['bio05','bio06','bio07','bio08']).unmask(0)

    @staticmethod
    def get_climate_part3():
        """bio09-bio12"""
        clim = ee.Image("projects/# replace with your GEE project ID/assets/brazil_worldclim_v2")
        return clim.select(['b9','b10','b11','b12']).rename(['bio09','bio10','bio11','bio12']).unmask(0)

    @staticmethod
    def get_climate_part4():
        """bio13-bio16"""
        clim = ee.Image("projects/# replace with your GEE project ID/assets/brazil_worldclim_v2")
        return clim.select(['b13','b14','b15','b16']).rename(['bio13','bio14','bio15','bio16']).unmask(0)

    @staticmethod
    def get_climate_part5():
        """bio17-bio19"""
        clim = ee.Image("projects/# replace with your GEE project ID/assets/brazil_worldclim_v2")
        return clim.select(['b17','b18','b19']).rename(['bio17','bio18','bio19']).unmask(0)

    @staticmethod
    def get_topography():
        """Topography (2 bands): elevation and slope from SRTM."""
        topo  = ee.Image("USGS/SRTMGL1_003")
        slope = ee.Terrain.slope(topo).rename('bio_slope').unmask(0)
        elev  = topo.select('elevation').rename('bio_elevation').unmask(0)
        return elev.addBands(slope)

    @staticmethod
    def get_landcover():
        """Land Cover — ESA CCI 2000 reclassified to 11 simplified classes."""
        lc_raw    = ee.Image("projects/# replace with your GEE project ID/assets/lc_300m_2000").select('b1')
        from_vals = [10, 11, 12, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220]
        to_vals   = [1,   1,  1,  1,  3,  3,  4,  5,  6,  6,  6,   7,   7,   8,   9,   2,  10,  11,  11,  11,  10,  10,  10,  10]
        lc_11class = lc_raw.remap(from_vals, to_vals, 0).rename('bio_landcover_class').unmask(0)
        crop_weights = lc_raw.remap(
            [10, 11, 12, 20, 30, 40],
            [1.0, 1.0, 1.0, 1.0, 0.75, 0.25], 0
        ).rename('cropland_weight').unmask(0)
        return lc_11class.addBands(crop_weights).addBands(lc_raw.rename('lc_raw_value'))

    @staticmethod
    def get_npp_fire():
        """NPP & Fire (2 bands): MODIS productivity and fire frequency"""
        tropics = ee.Geometry.Rectangle([-180, -25, 180, 25])
        npp = ee.ImageCollection("MODIS/061/MOD17A3HGF") \
            .filterDate('2001-01-01', '2017-12-31') \
            .filterBounds(tropics) \
            .select('Npp') \
            .mean() \
            .rename('bio_npp_mean').unmask(0)

        def year_burned(year):
            year  = ee.Number(year)
            start = ee.Date.fromYMD(year, 1, 1)
            end   = ee.Date.fromYMD(year, 12, 31)
            return ee.ImageCollection("MODIS/061/MCD64A1") \
                .filterDate(start, end) \
                .filterBounds(tropics) \
                .select('BurnDate') \
                .map(lambda img: img.gt(0)) \
                .max() \
                .unmask(0)

        years    = ee.List.sequence(2001, 2016)
        fire_freq = ee.ImageCollection(years.map(year_burned)) \
            .sum() \
            .rename('bio_fire_freq') \
            .unmask(0)
        return npp.addBands(fire_freq)

    @staticmethod
    def get_soil_chemical():
        """Soil chemical properties + biome (3 bands)."""
        biomes_fc = ee.FeatureCollection("RESOLVE/ECOREGIONS/2017")
        biomes    = ee.Image().byte().paint(biomes_fc, 'BIOME_NUM').rename('bio_biome_id')

        def top30cm_average(image_path, band_prefix='b'):
            img = ee.Image(image_path)
            return img.select(f'{band_prefix}0').multiply(0.167) \
                .add(img.select(f'{band_prefix}10').multiply(0.333)) \
                .add(img.select(f'{band_prefix}30').multiply(0.500))

        soil_oc = top30cm_average("OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02") \
            .rename('bio_soil_organic_carbon').unmask(0)
        soil_ph = top30cm_average("OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02") \
            .rename('bio_soil_ph').unmask(0)
        return soil_oc.addBands(soil_ph).addBands(biomes)

    @staticmethod
    def get_soil_texture():
        """Soil Texture (2 bands): Sand and clay content"""
        def top30cm_average(image_path, band_prefix='b'):
            img = ee.Image(image_path)
            return img.select(f'{band_prefix}0').multiply(0.167) \
                .add(img.select(f'{band_prefix}10').multiply(0.333)) \
                .add(img.select(f'{band_prefix}30').multiply(0.500))

        soil_sand = top30cm_average("OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02") \
            .rename('bio_soil_sand').unmask(0)
        soil_clay = top30cm_average("OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02") \
            .rename('bio_soil_clay').unmask(0)
        return soil_sand.addBands(soil_clay)

    @staticmethod
    def get_soil_physical():
        """Soil Physical Properties (2 bands): Bulk density and water content"""
        def top30cm_average(image_path, band_prefix='b'):
            img = ee.Image(image_path)
            return img.select(f'{band_prefix}0').multiply(0.167) \
                .add(img.select(f'{band_prefix}10').multiply(0.333)) \
                .add(img.select(f'{band_prefix}30').multiply(0.500))

        soil_bd = top30cm_average("OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02") \
            .rename('bio_soil_bulk_density_fine').unmask(0)
        soil_wc = top30cm_average("OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v02") \
            .rename('bio_soil_water_33kpa').unmask(0)
        return soil_bd.addBands(soil_wc)

    @staticmethod
    def get_socioeconomic():
        """
        Socioeconomic (5 bands):
            econ_road_density      — GRIP4 road density
            econ_dist_water        — distance to water (30 arc-sec)
            paper_gdp_2015         — GDP per capita PPP 2015 (adm2-gridded)
            paper_ghs_pop          — GHS population 2015
            paper_protected_binary — WDPA protected area mask
        """
        dist_water = ee.Image("projects/# replace with your GEE project ID/assets/distance2water_30arcsec") \
            .rename('econ_dist_water').unmask(0)
        road_density = ee.Image("projects/# replace with your GEE project ID/assets/grip4_density") \
            .rename('econ_road_density').unmask(0)

        gdp_img = ee.Image("projects/sat-io/open-datasets/GRIDDED_HDI_GDP/adm2_gdp_perCapita_1990_2022")
        gdp = gdp_img.select(['PPP_2015']).rename('paper_gdp_2015')

        ghs_pop = ee.ImageCollection("JRC/GHSL/P2023A/GHS_POP") \
            .filterDate('2015-01-01', '2015-12-31') \
            .first().select([0]).rename('paper_ghs_pop')

        wdpa = ee.FeatureCollection("WCMC/WDPA/current/polygons")
        protected_mask = ee.Image(0).byte().paint(wdpa, 1).rename('paper_protected_binary')

        return dist_water.addBands(road_density).addBands(gdp).addBands(ghs_pop).addBands(protected_mask)

    @staticmethod
    def get_highres_lights():
        """
        VIIRS nightlights (1 band):
            highres_lights_2012_2014 — annual median 2012-2014 (earliest available)

        Sampled alone so VIIRS mask gaps (common in remote Amazon in early
        years) don't drag down the other highres variables via sampleRegions'
        all-bands-must-be-valid behavior.
        """
        return ee.ImageCollection("NOAA/VIIRS/DNB/ANNUAL_V21") \
            .filterDate('2012-01-01', '2014-12-31') \
            .select('median_masked') \
            .median() \
            .rename('highres_lights_2012_2014')

    @staticmethod
    def get_highres_travel():
        """Travel time to cities (1 band): Oxford/MAP accessibility 2015"""
        return ee.Image("Oxford/MAP/accessibility_to_cities_2015_v1_0") \
            .select('accessibility').rename('highres_travel')

    @staticmethod
    def get_highres_ghm():
        """
        Global Human Modification (1 band):
            highres_ghm — year-2000 overall HM at 300m (Kennedy et al. v3)

        Uses sat-io 1990-2020 time series filtered to year 2000 to match
        training-era land cover. ~1% of points in remote forest have no
        data (genuine zero HM areas not encoded in the dataset); these
        are left as NA rather than forced to 0.
        """
        return ee.ImageCollection("projects/sat-io/open-datasets/GHM/HM_1990_2020_OVERALL_300M") \
            .filter(ee.Filter.eq('year', 2000)) \
            .first().select('constant').rename('highres_ghm')

    @staticmethod
    def get_highres_worldpop():
        """Population density (1 band): WorldPop 100m, year 2000"""
        return ee.ImageCollection("WorldPop/GP/100m/pop") \
            .filterDate('2000-01-01', '2001-01-01') \
            .select('population') \
            .mosaic() \
            .rename('highres_worldpop_density') \
            .unmask(0)

    @staticmethod
    def get_pasture():
        """Pasture (1 band): Global Pasture Watch grassland class 2000"""
        gpw_class = ee.ImageCollection("projects/global-pasture-watch/assets/ggc-30m/v1/grassland_c") \
            .filterDate('2000-01-01', '2000-12-31') \
            .mosaic() \
            .unmask(0) \
            .rename('gpw_raw_class')
        return gpw_class

    @staticmethod
    def get_opportunity_cost():
        """Opportunity Cost (1 band): Agricultural value (SPAM)"""
        return ee.Image("projects/# replace with your GEE project ID/assets/Opp_Cost_USD_per_Ha") \
            .rename('econ_opp_cost_usd_ha') \
            .unmask(0)

    @staticmethod
    def get_plantation_raw():
        """Plantation (1 band): Xiao et al. planted forests"""
        plantation_rgb = ee.ImageCollection("projects/sat-io/open-datasets/GLOBAL-NATURAL-PLANTED-FORESTS") \
            .mosaic().unmask(0)
        return plantation_rgb.select('b1').rename('plantation_raw_b1')


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Brazil Step 2: Sample Environmental Layers')
    parser.add_argument('--sample_type', type=str, required=True,
                        choices=list(Config.SAMPLE_ASSETS.keys()),
                        help='Which points to sample: regrowth, nonregrowth, or annulus')
    parser.add_argument('--layer', type=str, default='all',
                        choices=['all'] + LAYER_GROUPS,
                        help='Layer group to sample (default: all)')
    parser.add_argument('--delay', type=int, default=2,
                        help='Seconds between task submissions (default: 2)')
    parser.add_argument('--point-batch', type=int, default=None,
                        help='Filter points to this batch_id value (e.g. 2,3,4,5 for nonregrowth). '
                             'Uses server-side ee.Filter.eq("batch_id") — no materialization overhead.')
    args = parser.parse_args()

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

    layers_to_run = LAYER_GROUPS if args.layer == 'all' else [args.layer]

    logger.info("=" * 60)
    logger.info("BRAZIL STEP 2: SAMPLE ENVIRONMENTAL LAYERS")
    logger.info("=" * 60)
    logger.info(f"Sample type: {args.sample_type}")
    logger.info(f"Layers:      {layers_to_run}")

    # ── Initialize GEE ────────────────────────────────────────────────────────
    try:
        ee.Initialize(project=config.PROJECT_ID)
        logger.info(f"✓ GEE initialized: {config.PROJECT_ID}")
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=config.PROJECT_ID)
        logger.info("✓ GEE authenticated and initialized")

    # ── Load points ───────────────────────────────────────────────────────────
    asset_id = config.SAMPLE_ASSETS[args.sample_type]
    logger.info(f"\nLoading points from: {asset_id}")
    points = ee.FeatureCollection(asset_id)

    if args.point_batch is not None:
        points = points.filter(ee.Filter.eq('batch_id', args.point_batch))
        logger.info(f"✓ Filtered to batch_id={args.point_batch}")
    else:
        logger.info("✓ Points loaded (full collection)")

    # ── Submit layer tasks ────────────────────────────────────────────────────
    logger.info(f"\n{'='*60}")
    logger.info(f"SUBMITTING {len(layers_to_run)} LAYER TASK(S)")
    logger.info(f"{'='*60}\n")

    submitted_tasks = []

    for i, layer_name in enumerate(layers_to_run, 1):
        logger.info(f"[{i}/{len(layers_to_run)}] {layer_name}")

        try:
            layer_fn    = getattr(LayerDefinitions, f'get_{layer_name}')
            layer_stack = layer_fn()

            band_names = layer_stack.bandNames().getInfo()
            logger.info(f"  ✓ Built stack: {len(band_names)} bands — {', '.join(band_names)}")

            sampled = layer_stack.sampleRegions(
                collection = points,
                scale      = config.SCALE,
                geometries = True
            )

            batch_suffix = f"_pb{args.point_batch}" if args.point_batch is not None else ""
            description  = f"Brazil_{args.sample_type}_{layer_name}{batch_suffix}"

            task = ee.batch.Export.table.toDrive(
                collection  = sampled,
                description = description,
                fileFormat  = 'CSV'
            )
            task.start()

            submitted_tasks.append({
                'batch':       args.sample_type,
                'layer':       layer_name,
                'task_id':     task.id,
                'description': description
            })
            logger.info(f"  ✓ Task submitted: {task.id}")
            logger.info(f"  ✓ Output: {description}.csv\n")

            if len(layers_to_run) > 1 and i < len(layers_to_run):
                time.sleep(args.delay)

        except Exception as e:
            logger.error(f"  ✗ Failed {layer_name}: {e}\n")
            continue

    # ── Summary ───────────────────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info(f"SUBMITTED {len(submitted_tasks)}/{len(layers_to_run)} TASKS")
    logger.info("=" * 60)
    for t in submitted_tasks:
        logger.info(f"  {t['layer']:<20} | {t['task_id']}")

    logger.info("\nMonitor: https://code.earthengine.google.com/tasks")
    logger.info("Download CSVs from Google Drive when complete")


if __name__ == '__main__':
    main()
