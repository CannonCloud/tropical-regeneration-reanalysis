#!/usr/bin/env python3
"""
Brazil Step 3: Distance & Density Variables
============================================

Computes spatially-intensive variables (distance transforms, neighborhood
densities) that cannot be done as simple pixel lookups in step2.

All variables use 1.5×1.5 degree grid tiling with a 50km buffer to keep each
GEE task within memory limits. Points are filtered per tile so there is
no double-counting.

Variables produced:
    forest            bio_dist_forest_2000, bio_forest_density_1km2
                      (Hansen treecover2000)
    urban             paper_dist_urban
                      (ESA CCI 2000, class 190)
    settlement        highres_dist_settlement, highres_worldpop_density_2km
                      (WorldPop GP/100m V2, year 2000)
    plantation        bio_dist_plantation, bio_plantation_density_1km,
                      bio_canopy_height_1km_mean
    nightlight_density  highres_nightlight_density_2km
                      (VIIRS DNB annual V21, median 2012-2014)
    cropland          highres_dist_cropland, highres_cropland_density_5km
                      (ESA CCI 2000, weighted classes 10/11/12/20/30/40)
    cultivated_grass  highres_dist_cultivated_grass,
                      highres_cultivated_grass_density_2km
                      (Global Pasture Watch 2000 — earliest available in v1)

Inputs (run separately per sample_type):
    regrowth:    brazil_regrowth_points
    nonregrowth: brazil_nonregrowth_points_Cleaned
    annulus:     brazil_nonregrowth_annulus_batch_00_Cleaned

Outputs:
    Google Drive folder: GEE_Brazil_{variable}/
    One CSV per tile:    Brazil_{variable}_{sample_type}_Lon{w}_Lat{s}.csv
    Merge all tiles on point_id in R.

Usage:
    # Run all variables:
    python brazil_step3_forest_variables.py --sample_type regrowth
    python brazil_step3_forest_variables.py --sample_type nonregrowth

    # Run a single variable:
    python brazil_step3_forest_variables.py --sample_type regrowth --variable forest
    python brazil_step3_forest_variables.py --sample_type regrowth --variable urban
    python brazil_step3_forest_variables.py --sample_type regrowth --variable nightlight_density

    # Debugging:
    python brazil_step3_forest_variables.py --sample_type regrowth --variable urban --dry-run
    python brazil_step3_forest_variables.py --sample_type regrowth --variable urban --test-one
    python brazil_step3_forest_variables.py --sample_type regrowth --variable urban --tile -63 -6

    # Fix "Reprojection output too large" on failed plantation tiles:
    python brazil_step3_forest_variables.py --sample_type annulus --variable plantation --tile -36.5 -5.5 --buffer 25000
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
    LOG_FILE   = OUTPUT_DIR / 'brazil_step3_forest_variables.log'

    # Grid settings — 3 degree tiles proven stable in global step3
    GRID_SIZE_DEG = 1.5
    BUFFER_SIZE_M = 50_000   # 50km handles edge effects for all kernels used here

    # Brazil bounding box (slightly padded beyond country border)
    LAT_MIN, LAT_MAX = -34, 6
    LON_MIN, LON_MAX = -74, -28

    # Sample assets
    SAMPLE_ASSETS = {
        'regrowth':    'projects/# replace with your GEE project ID/assets/brazil_regrowth_points',
        'nonregrowth': 'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_points_Cleaned',
        'annulus':     'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_annulus_batch_00_Cleaned',
    }


ALL_VARIABLES = [
    'forest',
    'urban',
    'settlement',
    'plantation',
    'nightlight_density',
    'cropland',
    'cultivated_grass',
]


# =============================================================================
# IMAGE BUILDERS
# Each function returns an ee.Image clipped to the buffered tile geometry.
# =============================================================================

def build_forest_image(geometry):
    """
    Band 1: bio_dist_forest_2000
        Euclidean distance in metres to nearest treecover2000 > 30%.
    Band 2: bio_forest_density_1km2
        Mean treecover2000 (0-100) within 564m radius (≈ 1km²).
    """
    hansen = ee.Image("UMD/hansen/global_forest_change_2018_v1_6")
    tc2000 = hansen.select('treecover2000')

    forest_mask = tc2000.gt(30).unmask(0).clip(geometry)

    dist_forest = forest_mask \
        .fastDistanceTransform(256) \
        .sqrt() \
        .multiply(30) \
        .rename('bio_dist_forest_2000')

    tc2000_filled = tc2000.unmask(0).clip(geometry)
    forest_density = tc2000_filled \
        .reduceNeighborhood(
            reducer = ee.Reducer.mean(),
            kernel  = ee.Kernel.circle(564, 'meters')
        ) \
        .rename('bio_forest_density_1km2')

    return dist_forest.addBands(forest_density)


def build_urban_image(geometry):
    """
    Band: paper_dist_urban
        Distance (m) to nearest ESA CCI class 190 (urban) pixel.
        Uses year-2000 land cover to match Williams training-era data.
    """
    lc_raw = ee.Image("projects/# replace with your GEE project ID/assets/lc_300m_2000").select('b1')
    urban_mask = lc_raw.eq(190).unmask(0).clip(geometry)
    dist = urban_mask.fastDistanceTransform(1024).sqrt().multiply(300)  # 300m pixels
    return dist.rename('paper_dist_urban')


def build_settlement_image(geometry):
    """
    Band 1: highres_dist_settlement
        Distance (m) to nearest WorldPop pixel with ≥ 1 person.
        Uses year-2000 (earliest available in WorldPop GP/100m V2).
    Band 2: highres_worldpop_density_2km
        Mean population in 2km radius.
    """
    worldpop = ee.ImageCollection("WorldPop/GP/100m/pop") \
        .filterDate('2000-01-01', '2001-01-01') \
        .select('population') \
        .mosaic()

    settlement_mask = worldpop.gte(1).unmask(0).clip(geometry)
    dist = settlement_mask.fastDistanceTransform(1024).sqrt().multiply(100)  # 100m pixels

    worldpop_filled = worldpop.unmask(0).clip(geometry)
    pop_density = worldpop_filled.reduceNeighborhood(
        reducer = ee.Reducer.mean(),
        kernel  = ee.Kernel.circle(2000, 'meters')
    )

    return dist.rename('highres_dist_settlement') \
               .addBands(pop_density.rename('highres_worldpop_density_2km'))


def build_plantation_image(geometry):
    """
    Band 1: bio_dist_plantation
        Distance (m) to nearest Xiao et al. planted forest patch (≥ 4 ha).
    Band 2: bio_plantation_density_1km
        Plantation fraction (0-100%) in 1km radius.
    Band 3: bio_canopy_height_1km_mean
        Mean GEDI canopy height in 1km radius.

    Note: no .reproject() calls. The Xiao dataset is natively 30m, so
    connectedPixelCount and fastDistanceTransform already operate at
    that resolution through lazy evaluation. Forcing .reproject() caused
    "Reprojection output too large" errors on buffered 1.5° tiles
    (~10,600 × 10,600 pixels exceeds GEE's eager materialization limit).
    The forest builder uses the same pattern (no reproject) successfully.
    """
    xiao = ee.ImageCollection("projects/sat-io/open-datasets/GLOBAL-NATURAL-PLANTED-FORESTS") \
        .mosaic()

    # Identify plantation pixels by RGB signature (yellow = high R, high G, low B)
    raw_plantation = xiao.select('b1').gte(100) \
        .And(xiao.select('b2').gte(100)) \
        .And(xiao.select('b3').lt(100)) \
        .unmask(0).clip(geometry)

    # Filter to patches ≥ 4 ha (45 pixels at 30m)
    connected_count = raw_plantation.connectedPixelCount(maxSize=256, eightConnected=True)
    is_plantation   = raw_plantation.updateMask(connected_count.gte(45)).unmask(0)

    dist = is_plantation \
        .fastDistanceTransform(2048).sqrt().multiply(30)

    local_density = is_plantation.reduceNeighborhood(
        reducer = ee.Reducer.mean(),
        kernel  = ee.Kernel.circle(564, 'meters')
    ).multiply(100)

    glad_height = ee.ImageCollection("projects/sat-io/open-datasets/GLAD/GEDI_V27") \
        .mosaic().unmask(0)
    mean_height = glad_height.reduceNeighborhood(
        reducer = ee.Reducer.mean(),
        kernel  = ee.Kernel.circle(564, 'meters')
    )

    return dist.rename('bio_dist_plantation') \
               .addBands(local_density.rename('bio_plantation_density_1km')) \
               .addBands(mean_height.rename('bio_canopy_height_1km_mean'))


def build_nightlight_density_image(geometry):
    """
    Band: highres_nightlight_density_2km
        VIIRS DNB median 2012-2014 (earliest available period), smoothed to
        mean within 2km radius. Provides landscape-level light context.
    """
    viirs = ee.ImageCollection('NOAA/VIIRS/DNB/ANNUAL_V21') \
        .filterDate('2012-01-01', '2014-12-31') \
        .select('median_masked') \
        .median()

    viirs_filled = viirs.unmask(0).clip(geometry)

    nightlight_density = viirs_filled.reduceNeighborhood(
        reducer = ee.Reducer.mean(),
        kernel  = ee.Kernel.circle(2000, 'meters')
    )

    return nightlight_density.rename('highres_nightlight_density_2km')


def build_cropland_image(geometry):
    """
    Band 1: highres_dist_cropland
        Distance (m) to nearest ESA CCI year-2000 cropland pixel.
        Any pixel with cropland weight > 0 counts as cropland for distance.
    Band 2: highres_cropland_density_5km
        Weighted cropland density (0-100%) in 5km radius.
        Weights match step2 cropland_weight: classes 10/11/12/20 = 1.0,
        class 30 (mosaic >50% cropland) = 0.75, class 40 (mosaic <50%) = 0.25.
    """
    lc_raw = ee.Image("projects/# replace with your GEE project ID/assets/lc_300m_2000").select('b1')

    # Weighted cropland mask — matches step2 cropland_weight exactly
    crop_weights = lc_raw.remap(
        [10, 11, 12, 20, 30, 40],
        [1.0, 1.0, 1.0, 1.0, 0.75, 0.25], 0
    ).unmask(0).clip(geometry)

    # Distance: based on presence of any cropland (weight > 0)
    dist = crop_weights.gt(0).fastDistanceTransform(1024).sqrt().multiply(300)  # 300m pixels

    # Density: weighted mean in 5km radius
    crop_density = crop_weights.reduceNeighborhood(
        reducer = ee.Reducer.mean(),
        kernel  = ee.Kernel.circle(5000, 'meters')
    ).multiply(100)

    return dist.rename('highres_dist_cropland') \
               .addBands(crop_density.rename('highres_cropland_density_5km'))


def build_cultivated_grass_image(geometry):
    """
    Band 1: highres_dist_cultivated_grass
        Distance (m) to nearest Global Pasture Watch cultivated grassland pixel (class 1).
    Band 2: highres_cultivated_grass_density_2km
        Cultivated grassland fraction (0-100%) in 2km radius.
    Uses year 2000 (earliest available in GPW v1, which spans 2000-2022).
    """
    gpw_class = ee.ImageCollection("projects/global-pasture-watch/assets/ggc-30m/v1/grassland_c") \
        .filterDate('2000-01-01', '2000-12-31') \
        .mosaic()

    cultivated_mask = gpw_class.eq(1).unmask(0).clip(geometry)

    dist = cultivated_mask.fastDistanceTransform(1024).sqrt().multiply(30)  # 30m pixels

    cultivated_density = cultivated_mask.reduceNeighborhood(
        reducer = ee.Reducer.mean(),
        kernel  = ee.Kernel.circle(2000, 'meters')
    ).multiply(100)

    return dist.rename('highres_dist_cultivated_grass') \
               .addBands(cultivated_density.rename('highres_cultivated_grass_density_2km'))


# ── Dispatch table ────────────────────────────────────────────────────────────
IMAGE_BUILDERS = {
    'forest':             build_forest_image,
    'urban':              build_urban_image,
    'settlement':         build_settlement_image,
    'plantation':         build_plantation_image,
    'nightlight_density': build_nightlight_density_image,
    'cropland':           build_cropland_image,
    'cultivated_grass':   build_cultivated_grass_image,
}


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Brazil Step 3: Distance & density variables (tiled GEE exports)')
    parser.add_argument('--sample_type', required=True,
                        choices=['regrowth', 'nonregrowth', 'annulus'],
                        help='Which point asset to process')
    parser.add_argument('--variable', type=str, default='all',
                        choices=['all'] + ALL_VARIABLES,
                        help='Which variable to compute (default: all)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print tiles without submitting tasks')
    parser.add_argument('--test-one', action='store_true',
                        help='Submit only the first valid tile then exit')
    parser.add_argument('--tile', nargs=2, type=float, metavar=('LON', 'LAT'),
                        help='Resubmit a single tile by SW corner e.g. --tile -63 -6 or --tile -72.5 -4.5')
    parser.add_argument('--delay', type=int, default=0,
                        help='Seconds between submissions (default: 1)')
    parser.add_argument('--tile-scale', type=int, default=4,
                        help='GEE tileScale for sampleRegions (default: 4). '
                             'Increase to 8 or 16 for computationally heavy tiles.')
    parser.add_argument('--buffer', type=int, default=None,
                        help='Override buffer size in metres (default: 50000). '
                             'Use 25000 or 30000 to fix "Reprojection output too large" '
                             'errors on plantation tiles near the pixel-count limit.')
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

    variables_to_run = ALL_VARIABLES if args.variable == 'all' else [args.variable]

    logger.info("=" * 70)
    logger.info("BRAZIL STEP 3: DISTANCE & DENSITY VARIABLES")
    logger.info("=" * 70)
    logger.info(f"Sample type: {args.sample_type}")
    logger.info(f"Variables:   {variables_to_run}")
    logger.info(f"Grid:        {config.GRID_SIZE_DEG} deg tiles + {config.BUFFER_SIZE_M//1000}km buffer (~{int(config.GRID_SIZE_DEG * 111)}km squares)")

    buffer_m = args.buffer if args.buffer is not None else config.BUFFER_SIZE_M
    if args.buffer is not None:
        logger.info(f"Buffer:      OVERRIDDEN to {buffer_m//1000}km (default {config.BUFFER_SIZE_M//1000}km)")

    if args.dry_run:    logger.info("\n  DRY RUN — no tasks will be submitted")
    if args.test_one:   logger.info("\n  TEST MODE — first valid tile only")
    if args.tile:       logger.info(f"\n  TILE MODE — Lon{args.tile[0]}_Lat{args.tile[1]} only")

    # ── Initialize GEE ────────────────────────────────────────────────────────
    try:
        ee.Initialize(project=config.PROJECT_ID)
        logger.info(f"\n✓ GEE initialized: {config.PROJECT_ID}")
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=config.PROJECT_ID)
        logger.info("✓ GEE authenticated and initialized")

    # ── Load points ───────────────────────────────────────────────────────────
    asset_id = config.SAMPLE_ASSETS[args.sample_type]
    logger.info(f"\nPoints: {asset_id}")
    points = ee.FeatureCollection(asset_id)
    logger.info("✓ Points loaded")

    # ── Grid ──────────────────────────────────────────────────────────────────
    lon_steps   = int((config.LON_MAX - config.LON_MIN) / config.GRID_SIZE_DEG) + 1
    lat_steps   = int((config.LAT_MAX - config.LAT_MIN) / config.GRID_SIZE_DEG) + 1
    total_cells = lon_steps * lat_steps

    logger.info(f"\nGrid: {lon_steps} x {lat_steps} = {total_cells} potential tiles")

    # ── Loop over variables, then tiles ───────────────────────────────────────
    for var_name in variables_to_run:
        builder = IMAGE_BUILDERS[var_name]
        folder  = f'GEE_Brazil_{var_name}'
        logger.info(f"\n{'='*70}")
        logger.info(f"VARIABLE: {var_name}  →  Drive folder: {folder}/")
        logger.info(f"{'='*70}")

        submitted = 0
        skipped   = 0
        scanned   = 0

        for i in range(lon_steps):
            for j in range(lat_steps):
                scanned += 1

                if scanned % 20 == 0:
                    logger.info(f"  Progress: {scanned}/{total_cells} scanned | "
                                f"{submitted} submitted | {skipped} skipped")

                w = round(config.LON_MIN + i * config.GRID_SIZE_DEG, 6)
                e = round(w + config.GRID_SIZE_DEG, 6)
                s = round(config.LAT_MIN + j * config.GRID_SIZE_DEG, 6)
                n = round(s + config.GRID_SIZE_DEG, 6)

                tile_geom = ee.Geometry.Rectangle([w, s, e, n], 'EPSG:4326', False)
                tile_id   = f"Lon{w:g}_Lat{s:g}"

                # Tile filter (--tile mode)
                if args.tile is not None:
                    if w != args.tile[0] or s != args.tile[1]:
                        skipped += 1
                        continue

                # Skip empty tiles
                try:
                    has_points = points.filterBounds(tile_geom).limit(1).size().getInfo()
                except Exception:
                    has_points = 1  # conservative

                if has_points == 0:
                    skipped += 1
                    continue

                submitted += 1
                logger.info(f"\n  Tile {submitted}: {tile_id}  "
                            f"({w:g} to {e:g}E, {s:g} to {n:g}N)")

                try:
                    calc_geom      = tile_geom.buffer(buffer_m)
                    img            = builder(calc_geom)
                    points_in_tile = points.filterBounds(tile_geom)

                    sampled = img.sampleRegions(
                        collection = points_in_tile,
                        scale      = 30,
                        geometries = True,
                        tileScale  = args.tile_scale
                    )

                    task_desc = f"Brazil_{var_name}_{args.sample_type}_{tile_id}"

                    if args.dry_run:
                        logger.info(f"  [DRY RUN] Would submit: {task_desc}")
                        if args.test_one:
                            logger.info("  First tile found. Remove --dry-run to submit.")
                            return
                    else:
                        task = ee.batch.Export.table.toDrive(
                            collection  = sampled,
                            description = task_desc,
                            folder      = folder,
                            fileFormat  = 'CSV'
                        )
                        task.start()
                        logger.info(f"  ✓ {task.id}")

                        if args.test_one:
                            logger.info("\n  TEST MODE: first task submitted. "
                                        "Check it runs cleanly, then re-run without --test-one.")
                            return

                        time.sleep(args.delay)

                except Exception as e:
                    logger.error(f"  ✗ Failed {tile_id}: {e}")
                    continue

        logger.info(f"\n  {var_name}: {submitted} tiles submitted, "
                    f"{skipped} skipped, {scanned} scanned")

    # ── Final summary ─────────────────────────────────────────────────────────
    logger.info(f"\n{'=' * 70}")
    logger.info(f"DONE")
    logger.info(f"{'=' * 70}")
    if not args.dry_run:
        logger.info("\nMonitor: https://code.earthengine.google.com/tasks")
        logger.info("Output folders per variable: GEE_Brazil_{variable}/")
        logger.info("\nMerge in R (example for one variable):")
        logger.info("  files <- list.files('GEE_Brazil_urban', '*.csv', full.names=TRUE)")
        logger.info("  df    <- dplyr::bind_rows(lapply(files, readr::read_csv))")
        logger.info("  training <- dplyr::left_join(training, df, by='point_id')")


if __name__ == '__main__':
    main()
