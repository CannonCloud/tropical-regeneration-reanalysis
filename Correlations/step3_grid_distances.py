#!/usr/bin/env python3
"""
Step 3: Grid-Based Distance Calculations (Optimized)
====================================================

This script uses a grid-tiling approach to eliminate redundant calculations.
Instead of computing distance maps 200K times (once per point), it:
1. Divides the tropics into 3°×3° tiles (~330km squares)
2. Calculates distance map ONCE per tile
3. Samples all points in that tile from the same map

Massive speedup for large datasets and future reruns.

Usage:
    # Run all distance types sequentially
    python step3_grid_distances.py --samples 200000
    
    # Or run one at a time
    python step3_grid_distances.py --samples 200000 --distance plantation
    python step3_grid_distances.py --samples 200000 --distance forest
    ... etc
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
    OUTPUT_DIR = Path('./gee_outputs')
    LOG_FILE = OUTPUT_DIR / 'step3_grid_distances.log'
    
    # Grid settings
    GRID_SIZE_DEG = 3.0  # 3 degrees (~330km) - proven stable, no task failures
    BUFFER_SIZE_M = 50000  # 50km buffer for edge effects
    
    # Tropics bounds
    LAT_MIN, LAT_MAX = -25, 25
    LON_MIN, LON_MAX = -180, 180
    
    def get_asset_id(self, n_samples):
        return f'projects/{self.PROJECT_ID}/assets/master_sample_points_{n_samples}'

# All distance/density types
CALCULATION_TYPES = [
    'cropland_density',         # Special case: 5km buffer density, not distance (old ESA CCI)
    'forest',                   # Distance to forest + 1km density
    'urban',                    # Distance to urban areas
    'settlement',               # Distance to populated areas + 2km population density
    'agriculture',              # Distance to cropland OR cultivated grassland + 1km density
    'plantation',               # Distance to planted forests + 1km density
    'nightlight_density',       # 2km nightlight density (smoothed neighborhood context)
    'cropland',                 # Distance to cropland only + 2km density
    'cultivated_grass'          # Distance to cultivated grassland only + 2km density
]

# =============================================================================
# DISTANCE/DENSITY IMAGE BUILDERS
# =============================================================================

def get_calculation_image(calc_type, geometry, buffer_5km=False):
    """
    Returns the calculation image (distance or density) for a buffered tile.
    
    Args:
        calc_type: Type of calculation (forest, urban, plantation, etc.)
        geometry: Buffered tile geometry
        buffer_5km: For cropland_density, use 5km buffer instead of fastDistanceTransform
    
    Returns:
        ee.Image with the calculated values
    """
    
    # =========================================================================
    # CROPLAND DENSITY (special case - not a distance)
    # =========================================================================
    if calc_type == 'cropland_density':
        lc_raw = ee.Image("projects/# replace with your GEE project ID/assets/lc_300m_2015").select('b1')
        crop_weights = lc_raw.remap(
            [10, 11, 12, 20, 30, 40], 
            [1.0, 1.0, 1.0, 1.0, 0.75, 0.25], 
            0
        ).unmask(0).clip(geometry)
        
        # For density, we'll compute focal_mean with 5km radius when sampling
        # Return the weights image directly
        return crop_weights.rename('bio_cropland_density')
    
    # =========================================================================
    # FOREST DISTANCE (uses fastDistanceTransform for speed and reliability)
    # =========================================================================
    elif calc_type == 'forest':
        hansen = ee.Image("UMD/hansen/global_forest_change_2018_v1_6")
        tc2000 = hansen.select('treecover2000')
        lossyear = hansen.select('lossyear')
        gain = hansen.select('gain')
        
        # Forest = >30% cover in 2000, no loss through 2017, OR gained
        cover_post_loss = tc2000.where(lossyear.gte(1).And(lossyear.lte(17)), 0)
        is_stable_gain = gain.eq(1).And(lossyear.eq(0).Or(lossyear.gt(12)))
        forest_2018 = cover_post_loss.where(is_stable_gain, 100)
        
        # BAND 1: Distance to forest
        # Create mask: 1=forest, 0=not forest, ensuring NULLs become 0
        forest_mask = forest_2018.gt(30).unmask(0).clip(geometry)
        dist = forest_mask.fastDistanceTransform(1024).sqrt().multiply(30)  # 30m pixels
        
        # BAND 2: Local forest density in 1km² area (Williams et al. key variable!)
        # 1km² = 1,000,000 m² = circle with radius ~564m
        # We'll use a 564m radius circular kernel to approximate 1km² area
        # CRITICAL: Unmask NULLs to 0 so they're counted as 0% forest (not ignored)
        forest_percent = forest_2018.unmask(0).clip(geometry)  # 0-100 scale, NULLs=0
        
        local_density = forest_percent.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(564, 'meters')  # 564m radius ≈ 1km² area
        )
        
        # Return both as a multi-band image
        return dist.rename('bio_dist_forest_2018').addBands(
            local_density.rename('bio_forest_density_1km2')
        )
    
    # =========================================================================
    # URBAN DISTANCE
    # =========================================================================
    elif calc_type == 'urban':
        lc_raw = ee.Image("projects/# replace with your GEE project ID/assets/lc_300m_2015").select('b1')
        # Compare FIRST (NULL.eq(190) → NULL), then unmask (NULL → 0)
        urban_mask = lc_raw.eq(190).unmask(0).clip(geometry)  # Class 190 = urban
        
        dist = urban_mask.fastDistanceTransform(1024).sqrt().multiply(300)  # 300m pixels
        
        return dist.rename('paper_dist_urban')
    
    # =========================================================================
    # SETTLEMENT DISTANCE (based on WorldPop)
    # =========================================================================
    elif calc_type == 'settlement':
        # WorldPop 2015 (100m resolution)
        worldpop = ee.ImageCollection("WorldPop/GP/100m/pop") \
            .filterDate('2015-01-01', '2016-01-01') \
            .select('population') \
            .mosaic()
        
        # BAND 1: Distance to any settlement (population ≥ 1)
        # Compare FIRST (NULL.gte(1) → NULL), then unmask (NULL → 0)
        settlement_mask = worldpop.gte(1).unmask(0).clip(geometry)
        dist = settlement_mask.fastDistanceTransform(1024).sqrt().multiply(100)  # 100m pixels
        
        # BAND 2: Population density in 2km radius
        # 2km provides neighborhood context while remaining computationally stable
        worldpop_filled = worldpop.unmask(0).clip(geometry)
        pop_density = worldpop_filled.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(2000, 'meters')  # 2km radius
        )
        
        return dist.rename('highres_dist_settlement').addBands(
            pop_density.rename('highres_worldpop_density_2km')
        )
    
    # =========================================================================
    # AGRICULTURE DISTANCE (cropland OR cultivated grassland) + 1km density
    # =========================================================================
    elif calc_type == 'agriculture':
        worldcover = ee.ImageCollection("ESA/WorldCover/v100").first()
        gpw_class = ee.ImageCollection("projects/global-pasture-watch/assets/ggc-30m/v1/grassland_c") \
            .filterDate('2015-01-01', '2015-12-31') \
            .mosaic()
        
        # Compare FIRST for both layers, then unmask
        crop = worldcover.select('Map').eq(40).unmask(0).clip(geometry)       # ESA cropland
        cultivated = gpw_class.eq(1).unmask(0).clip(geometry)                 # GPW cultivated grassland
        ag_mask = crop.Or(cultivated)  # Both already 0/1 binary
        
        # BAND 1: Distance to agriculture
        dist = ag_mask.fastDistanceTransform(1024).sqrt().multiply(30)  # 30m pixels
        
        # BAND 2: Local agriculture density in 1km radius (0-100%)
        ag_density = ag_mask.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(564, 'meters')  # 564m radius ≈ 1km² area
        ).multiply(100)
        
        return dist.rename('highres_dist_agriculture').addBands(
            ag_density.rename('highres_agriculture_density_1km')
        )
    
    # =========================================================================
    # PLANTATION DISTANCE (using Xiao et al. Global Planted Forests)
    # =========================================================================
    elif calc_type == 'plantation':
        # 1. Load Xiao et al. Global Natural & Planted Forests
        xiao_collection = ee.ImageCollection("projects/sat-io/open-datasets/GLOBAL-NATURAL-PLANTED-FORESTS")
        xiao = xiao_collection.mosaic()
        
        # 2. Identify plantation pixels (Yellow = high R + high G + low B)
        is_red_high = xiao.select('b1').gte(100)
        is_green_high = xiao.select('b2').gte(100)
        is_blue_low = xiao.select('b3').lt(100)
        
        raw_plantation = is_red_high.And(is_green_high).And(is_blue_low).unmask(0).clip(geometry)
        
        # 3. Reproject BEFORE connectedPixelCount (force 30m scale)
        raw_plantation_30m = raw_plantation.reproject(crs='EPSG:3857', scale=30)
        
        # 4. Filter by patch size - only keep plantations ≥ 4 hectares
        connected_count = raw_plantation_30m.connectedPixelCount(maxSize=256, eightConnected=True)
        is_plantation = raw_plantation_30m.updateMask(connected_count.gte(45)).unmask(0)
        
        # 5. Distance calculation
        # fastDistanceTransform measures distance TO nearest mask=True (plantation=1)
        # NO .Not() needed - this was the bug the whole time!
        dist = is_plantation \
            .reproject(crs='EPSG:3857', scale=30) \
            .fastDistanceTransform(2048).sqrt().multiply(30)
        
        # 6. Density calculation (1km radius)
        local_density = is_plantation.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(564, 'meters')
        ).multiply(100)
        
        # 7. Mean canopy height
        glad_collection = ee.ImageCollection("projects/sat-io/open-datasets/GLAD/GEDI_V27")
        glad_height = glad_collection.mosaic().unmask(0)
        
        mean_height = glad_height.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(564, 'meters')
        )
        
        return dist.rename('bio_dist_plantation') \
            .addBands(local_density.rename('bio_plantation_density_1km')) \
            .addBands(mean_height.rename('bio_canopy_height_1km_mean'))
    
    # =========================================================================
    # NIGHTLIGHT DENSITY (2km radius mean)
    # =========================================================================
    elif calc_type == 'nightlight_density':
        # VIIRS DNB Monthly composite - 2014-2016 median (matches existing point values)
        viirs = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG') \
            .filterDate('2014-01-01', '2016-12-31') \
            .select('avg_rad') \
            .median()  # Median across 2014-2016 period
        
        # Unmask NULLs to 0 (no light = 0, not missing data)
        viirs_filled = viirs.unmask(0).clip(geometry)
        
        # Calculate mean nightlight in 2km radius
        # This provides smoothed neighborhood context vs point-level values
        nightlight_density = viirs_filled.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(2000, 'meters')  # 2km radius
        )
        
        return nightlight_density.rename('highres_nightlight_density_2km')
    
    # =========================================================================
    # CROPLAND DISTANCE + 2KM DENSITY (balance between impact scale and computation)
    # =========================================================================
    elif calc_type == 'cropland':
        worldcover = ee.ImageCollection("ESA/WorldCover/v100").first()
        
        # Compare FIRST, then unmask (NULL.eq(40) → NULL → unmask(0) → 0)
        crop_mask = worldcover.select('Map').eq(40).unmask(0).clip(geometry)
        
        # BAND 1: Distance to cropland
        dist = crop_mask.fastDistanceTransform(1024).sqrt().multiply(10)  # 10m pixels
        
        # BAND 2: Local cropland density in 2km radius (0-100%)
        # 2km balances ecological relevance with computational feasibility
        crop_density = crop_mask.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(2000, 'meters')  # 2km radius
        ).multiply(100)
        
        return dist.rename('highres_dist_cropland').addBands(
            crop_density.rename('highres_cropland_density_2km')
        )
    
    # =========================================================================
    # CULTIVATED GRASSLAND DISTANCE (Global Pasture Watch class 1 only) + 2km density
    # =========================================================================
    elif calc_type == 'cultivated_grass':
        gpw_class = ee.ImageCollection("projects/global-pasture-watch/assets/ggc-30m/v1/grassland_c") \
            .filterDate('2015-01-01', '2015-12-31') \
            .mosaic()
        
        # Compare FIRST, then unmask (NULL.eq(1) → NULL → unmask(0) → 0)
        cultivated_mask = gpw_class.eq(1).unmask(0).clip(geometry)
        
        # BAND 1: Distance to cultivated grassland
        dist = cultivated_mask.fastDistanceTransform(1024).sqrt().multiply(30)  # 30m pixels
        
        # BAND 2: Local cultivated grassland density in 2km radius (0-100%)
        # 2km provides landscape context while remaining computationally stable
        cultivated_density = cultivated_mask.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(2000, 'meters')  # 2km radius
        ).multiply(100)
        
        return dist.rename('highres_dist_cultivated_grass').addBands(
            cultivated_density.rename('highres_cultivated_grass_density_2km')
        )
    
    else:
        raise ValueError(f"Unknown calculation type: {calc_type}")

# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Grid-Based Distance Calculations')
    parser.add_argument('--samples', type=int, required=True,
                        help='Sample size (must match asset, e.g., 200000)')
    parser.add_argument('--distance', type=str, default='all',
                        choices=['all'] + CALCULATION_TYPES,
                        help='Which calculation to run (default: all)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be submitted without actually submitting')
    parser.add_argument('--test-one', action='store_true',
                        help='Submit only ONE task (first valid tile) as a smoke test')
    parser.add_argument('--tile', type=str, default=None,
                        help='Submit only this specific tile, e.g. --tile "Lon138_Lat-13". '
                             'Grid size is halved automatically for reliability. '
                             'Use for debugging failed tiles or verifying logic.')
    parser.add_argument('--tile-scale', type=int, default=4,
                        help='GEE tileScale for computation (default: 4). '
                             'Increase to 8 or 16 for tiles with computational timeout. '
                             'Higher values = more memory but slower.')
    parser.add_argument('--delay', type=int, default=2,
                        help='Seconds between task submissions (default: 2)')
    
    args = parser.parse_args()
    
    # Parse --tile argument if provided
    # Expected format: "Lon138_Lat-13" (matching GEE task name convention)
    # Snaps user input to nearest grid corner so you don't need to know exact grid coords
    target_tile = None
    if args.tile:
        try:
            parts = args.tile.replace('Lon', '').split('_Lat')
            target_lon = int(parts[0])
            target_lat = int(parts[1])
            
            # Snap to grid: find which tile's bottom-left corner contains this coordinate
            snapped_lon = int(Config.LON_MIN + int((target_lon - Config.LON_MIN) // Config.GRID_SIZE_DEG) * Config.GRID_SIZE_DEG)
            snapped_lat = int(Config.LAT_MIN + int((target_lat - Config.LAT_MIN) // Config.GRID_SIZE_DEG) * Config.GRID_SIZE_DEG)
            target_tile = (snapped_lon, snapped_lat)
            
            if snapped_lon != target_lon or snapped_lat != target_lat:
                print(f"📍 Snapped Lon{target_lon}_Lat{target_lat} → grid tile Lon{snapped_lon}_Lat{snapped_lat}")
            else:
                print(f"📍 Tile Lon{snapped_lon}_Lat{snapped_lat} is exactly on the grid")
        except Exception:
            print(f"ERROR: Could not parse --tile '{args.tile}'. "
                  f"Expected format: 'Lon138_Lat-13'")
            exit(1)
    
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
    
    logger.info("="*70)
    logger.info("STEP 3: GRID-BASED DISTANCE CALCULATIONS")
    logger.info("="*70)
    logger.info(f"Sample size: {args.samples:,}")
    logger.info(f"Grid size: {config.GRID_SIZE_DEG}° (~{int(config.GRID_SIZE_DEG * 111)}km)")
    logger.info(f"Buffer: {config.BUFFER_SIZE_M/1000:.0f}km")
    
    if args.dry_run:
        logger.info("\n⚠️  DRY RUN MODE - No tasks will be submitted")
    
    if args.test_one:
        logger.info("\n🧪 TEST MODE - Will submit only ONE task (first valid tile)")
        logger.info("   Use this to verify the task completes successfully before full run")
    
    if target_tile:
        logger.info(f"\n🔍 TILE MODE - Submitting ONLY tile: Lon{target_tile[0]}_Lat{target_tile[1]}")
        logger.info("   Use this to resubmit failed tiles or debug specific locations")
    
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
    
    # Determine which calculations to run
    if args.distance == 'all':
        calc_types = CALCULATION_TYPES
    else:
        calc_types = [args.distance]
    
    logger.info(f"\nCalculations to run: {', '.join(calc_types)}")
    
    # Generate grid
    lat_steps = int((config.LAT_MAX - config.LAT_MIN) / config.GRID_SIZE_DEG) + 1
    lon_steps = int((config.LON_MAX - config.LON_MIN) / config.GRID_SIZE_DEG)
    total_cells = lat_steps * lon_steps
    
    logger.info(f"\nGrid dimensions: {lon_steps} × {lat_steps} = {total_cells:,} potential tiles")
    logger.info("Scanning for tiles with points (skipping oceans/deserts)...\n")
    
    # For each calculation type
    for calc_type in calc_types:
        logger.info("="*70)
        logger.info(f"PROCESSING: {calc_type.upper()}")
        logger.info("="*70)
        
        submitted_count = 0
        skipped_count = 0
        scanned_count = 0
        
        # Scan grid
        for i in range(lon_steps):
            for j in range(lat_steps):
                scanned_count += 1
                
                # Progress update every 50 tiles
                if scanned_count % 50 == 0:
                    logger.info(f"  Progress: Scanned {scanned_count}/{total_cells} tiles, " +
                              f"submitted {submitted_count}, skipped {skipped_count}")
                
                # Define tile boundaries
                w = config.LON_MIN + (i * config.GRID_SIZE_DEG)
                e = w + config.GRID_SIZE_DEG
                s = config.LAT_MIN + (j * config.GRID_SIZE_DEG)
                n = s + config.GRID_SIZE_DEG
                
                # If --tile specified, skip all tiles except the matching one
                if target_tile is not None:
                    if round(w) != target_tile[0] or round(s) != target_tile[1]:
                        skipped_count += 1
                        continue
                
                # Create tile geometry (STRICT boundaries - no overlap)
                tile_geom = ee.Geometry.Rectangle([w, s, e, n], 'EPSG:4326', False)
                
                # Check if points exist in this tile
                points_in_tile = master_spine.filterBounds(tile_geom)
                
                # Quick existence check (don't count all points, just see if any exist)
                try:
                    count = points_in_tile.limit(1).size().getInfo()
                except:
                    # If check fails, assume tile might have points (be conservative)
                    count = 1
                
                if count == 0:
                    skipped_count += 1
                    continue  # Skip empty tiles
                
                # Process this tile
                submitted_count += 1
                # Use round() not int() to avoid floating point drift
                # e.g. -25 + 4*3.0 = -13.000000000000002 → int = -14 (WRONG)
                tile_id = f"Lon{round(w)}_Lat{round(s)}"
                
                logger.info(f"\nTile {submitted_count}: {tile_id}")
                logger.info(f"  Bounds: ({w:.1f}°, {s:.1f}°) to ({e:.1f}°, {n:.1f}°)")
                logger.info(f"  Center: ({(w+e)/2:.1f}°E, {(s+n)/2:.1f}°N)")
                logger.info(f"  (Open in maps: https://www.google.com/maps/@{(s+n)/2:.2f},{(w+e)/2:.2f},8z)")
                
                try:
                    # Buffer tile for calculation (to see nearby features)
                    calc_geom = tile_geom.buffer(config.BUFFER_SIZE_M)
                    
                    # Get calculation image
                    if calc_type == 'cropland_density':
                        # Special case: density within 5km buffer
                        # Use RASTER approach (not vector) for efficiency
                        base_img = get_calculation_image(calc_type, calc_geom)
                        
                        # Calculate density map for entire tile at once using focal mean
                        # This is MUCH faster than buffering each point individually
                        density_map = base_img.reduceNeighborhood(
                            reducer=ee.Reducer.mean(),
                            kernel=ee.Kernel.circle(5000, 'meters')  # 5km radius
                        ).rename('bio_cropland_density')
                        
                        # Always include geometries so lon/lat appear in output CSV
                        sampled = density_map.sampleRegions(
                            collection=points_in_tile,
                            scale=300,
                            geometries=True
                        )
                    else:
                        # Standard distance calculation
                        calc_img = get_calculation_image(calc_type, calc_geom)
                        
                        # Use appropriate scale per dataset
                        # plantation/forest: 30m (Hansen/DeepMind native)
                        # urban: 300m (ESA CCI native)  
                        # settlement: 100m (WorldPop native)
                        # agriculture/cropland: 30m (WorldCover/GPW native)
                        scale_map = {
                            'forest': 30,
                            'plantation': 30,
                            'urban': 300,
                            'settlement': 100,
                            'agriculture': 30,
                            'cropland': 10,
                            'cultivated_grass': 30,
                            'nightlight_density': 500,
                            'population_density': 100,
                            'cropland_coverage': 10,
                        }
                        sample_scale = scale_map.get(calc_type, 30)
                        
                        # Always include geometries so lon/lat appear in output CSV
                        sampled = calc_img.sampleRegions(
                            collection=points_in_tile,
                            scale=sample_scale,
                            geometries=True,
                            tileScale=args.tile_scale  # Configurable via --tile-scale
                        )
                    
                    # Export
                    task_desc = f'Grid_{calc_type}_{tile_id}_N{args.samples}'
                    
                    if args.dry_run:
                        logger.info(f"  [DRY RUN] Would submit: {task_desc}")
                        
                        # In dry-run + test-one mode, exit after first tile
                        if args.test_one:
                            logger.info("\n" + "="*70)
                            logger.info("🧪 TEST MODE (DRY RUN): Found first valid tile")
                            logger.info("="*70)
                            logger.info(f"\nTile: {tile_id}")
                            logger.info(f"Would submit: {task_desc}")
                            logger.info("\nRemove --dry-run to submit this task for real")
                            logger.info("="*70)
                            return  # Exit script
                    else:
                        task = ee.batch.Export.table.toDrive(
                            collection=sampled,
                            description=task_desc,
                            folder=f'GEE_Grid_{calc_type}',
                            fileFormat='CSV'
                        )
                        
                        task.start()
                        task_id = task.id
                        
                        logger.info(f"  ✓ Task submitted: {task_id}")
                        logger.info(f"  ✓ Output folder: GEE_Grid_{calc_type}/")
                        
                        # If test-one mode, exit immediately after first task
                        if args.test_one:
                            logger.info("\n" + "="*70)
                            logger.info("🧪 TEST MODE: First task submitted successfully!")
                            logger.info("="*70)
                            logger.info(f"\nTask ID: {task_id}")
                            logger.info(f"Description: {task_desc}")
                            logger.info(f"\nMonitor at: https://code.earthengine.google.com/tasks")
                            logger.info("\nWait for this task to complete (~5-10 min)")
                            logger.info("If it succeeds, run the full batch without --test-one flag")
                            logger.info("\n" + "="*70)
                            return  # Exit the entire script
                        
                        # Brief delay to avoid API rate limits
                        time.sleep(args.delay)
                
                except Exception as e:
                    logger.error(f"  ✗ Failed to submit tile {tile_id}: {e}")
                    continue
        
        # Summary for this calculation type
        logger.info(f"\n{'='*70}")
        logger.info(f"SUMMARY: {calc_type.upper()}")
        logger.info(f"{'='*70}")
        logger.info(f"Submitted: {submitted_count} tiles")
        logger.info(f"Skipped: {skipped_count} empty tiles")
        logger.info(f"Total scanned: {total_cells} tiles")
        logger.info(f"Coverage: {(submitted_count/total_cells)*100:.1f}% of grid has points")
        
        if not args.dry_run:
            logger.info(f"\nTasks running on GEE servers")
            logger.info(f"Monitor: https://code.earthengine.google.com/tasks")
            logger.info(f"Output folder: Google Drive → GEE_Grid_{calc_type}/")
        
        logger.info("")
    
    # Final summary
    logger.info("="*70)
    logger.info("ALL CALCULATIONS SUBMITTED")
    logger.info("="*70)
    
    if args.dry_run:
        logger.info("\nThis was a dry run. No tasks were actually submitted.")
        logger.info("Remove --dry-run flag to submit for real.")
    else:
        logger.info(f"\nSubmitted tasks for: {', '.join(calc_types)}")
        logger.info("\nNext steps:")
        logger.info("1. Monitor tasks: https://code.earthengine.google.com/tasks")
        logger.info("2. When complete, download all CSVs from Google Drive")
        logger.info("3. Organize into folders:")
        for calc_type in calc_types:
            logger.info(f"   gee_outputs/distance_{calc_type}/")
        logger.info("4. Run Step 4 to merge all data")
    
    logger.info("\n" + "="*70)
    
    if not args.dry_run:
        logger.info("\n💡 Pro tip: You can close this terminal now.")
        logger.info("   All tasks are running on Google's servers.\n")

if __name__ == '__main__':
    main()
