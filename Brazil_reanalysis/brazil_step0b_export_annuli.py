#!/usr/bin/env python3
"""
Brazil Step 0b: Export Near-Regrowth Nonregrowth Mask
======================================================

Creates a binary raster mask of eligible non-regrowth pixels within
a specified distance of any regrowth pixel. Uses fastDistanceTransform()
— a pure raster operation — to compute per-pixel distance to the nearest
regrowth pixel. This is fast, scalable, and avoids the vector buffering
memory errors that occur when buffering 100K individual points.

Why raster distance transform instead of vector buffering:
    Buffering 100K regrowth points and taking their union produces a
    polygon with millions of vertices that GEE cannot handle for inner
    bands. fastDistanceTransform() computes the same result as a single
    raster pass with no geometry construction at all.

Why use the full regrowth raster (not sampled points):
    Using only sampled regrowth points to define distance bands would
    miss unsampled regrowth pixels, causing some non-regrowth points
    to be incorrectly classified as "far from regrowth" when they aren't.
    The full raster mask guarantees every regrowth pixel contributes.

Output:
    brazil_nonregrowth_annulus_0_3km — binary mask, 1 = eligible
    non-regrowth pixel within 3km of any regrowth pixel. Sample this
    heavily in step1 Task C and do matching/pairing in R.

Usage:
    python brazil_step0b_export_annuli.py           # full Brazil
    python brazil_step0b_export_annuli.py --test    # Manaus box (~minutes)
    python brazil_step0b_export_annuli.py --dist 5000  # change distance threshold

After task completes, verify in GEE console then run:
    python brazil_step1_sample_points.py --task C --test
    python brazil_step1_sample_points.py --task C
"""

import ee
import argparse
import logging
import time
import json
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    PROJECT_ID = '# replace with your GEE project ID'

    REGROWTH_MASK_ASSET    = 'projects/# replace with your GEE project ID/assets/brazil_regrowth_mask_v1'
    NONREGROWTH_MASK_ASSET = 'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_mask_v4'

    # Distance threshold in metres — non-regrowth pixels beyond this are excluded
    DISTANCE_M = 3000

    # Hansen grid — output raster locked to this grid
    HANSEN_CRS_TRANSFORM = [0.00025, 0, -180, 0, -0.00025, 80]

    MAX_PIXELS = 1e13

    # Test mode — small box around Manaus
    TEST_BOUNDS = [-60.5, -3.6, -59.5, -2.6]

    OUTPUT_DIR = Path('./gee_outputs')
    LOG_FILE   = OUTPUT_DIR / 'brazil_step0b.log'


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Brazil Step 0b: Export near-regrowth nonregrowth mask')
    parser.add_argument('--test', action='store_true',
                        help='Run on small Manaus test region (~minutes)')
    parser.add_argument('--dist', type=int, default=None,
                        help=f'Distance threshold in metres (default: {Config.DISTANCE_M})')
    args = parser.parse_args()

    config = Config()
    if args.dist:
        config.DISTANCE_M = args.dist
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

    logger.info("=" * 60)
    logger.info("BRAZIL STEP 0b: EXPORT NEAR-REGROWTH NONREGROWTH MASK")
    logger.info("=" * 60)
    logger.info(f"Distance threshold: {config.DISTANCE_M}m")

    # ── Initialize GEE ────────────────────────────────────────────────────────
    try:
        ee.Initialize(project=config.PROJECT_ID)
        logger.info(f"✓ GEE initialized: {config.PROJECT_ID}")
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=config.PROJECT_ID)
        logger.info("✓ GEE authenticated and initialized")

    # ── Region ────────────────────────────────────────────────────────────────
    if args.test:
        b      = config.TEST_BOUNDS
        region = ee.Geometry.Rectangle([b[0], b[1], b[2], b[3]])
        asset_suffix = '_test'
        logger.info(f"✓ TEST MODE: Manaus box {config.TEST_BOUNDS}")
    else:
        region = (ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
                  .filter(ee.Filter.eq('country_na', 'Brazil'))
                  .geometry())
        asset_suffix = ''
        logger.info("✓ Full Brazil region loaded")

    # ── Build annulus mask ────────────────────────────────────────────────────
    logger.info("Loading regrowth and nonregrowth masks...")
    regrowth_mask    = ee.Image(config.REGROWTH_MASK_ASSET)
    nonregrowth_mask = ee.Image(config.NONREGROWTH_MASK_ASSET)

    # Raster distance transform — distance in metres to nearest regrowth pixel
    # fastDistanceTransform() returns distance in pixels, so we convert:
    #   pixels -> metres via sqrt(pixelArea)
    logger.info("Computing distance transform...")
    dist_to_regrowth = (regrowth_mask
                        .fastDistanceTransform()
                        .sqrt()
                        .multiply(ee.Image.pixelArea().sqrt())
                        .rename('dist_m'))

    # Eligible non-regrowth pixels within distance threshold
    dist_m = config.DISTANCE_M
    annulus_mask = (nonregrowth_mask
                    .where(dist_to_regrowth.gt(dist_m), 0)
                    .selfMask()
                    .rename('mask')
                    .toByte())

    logger.info(f"✓ Annulus mask built: nonregrowth within {dist_m}m of regrowth")

    # ── Export ────────────────────────────────────────────────────────────────
    dist_km    = dist_m // 1000
    asset_name = f'brazil_nonregrowth_annulus_0_{dist_km}km{asset_suffix}'
    asset_id   = f'projects/{config.PROJECT_ID}/assets/{asset_name}'

    task = ee.batch.Export.image.toAsset(
        image            = annulus_mask,
        description      = f'BrazilNonRegrowthAnnulus_0_{dist_km}km{asset_suffix}',
        assetId          = asset_id,
        region           = region,
        crs              = 'EPSG:4326',
        crsTransform     = config.HANSEN_CRS_TRANSFORM,
        maxPixels        = config.MAX_PIXELS,
        pyramidingPolicy = {'.default': 'mode'}
    )
    task.start()

    logger.info(f"✓ Task submitted: {task.id}")
    logger.info(f"✓ Asset: {asset_id}")

    logger.info("\n" + "=" * 60)
    logger.info("TASK SUBMITTED")
    logger.info("=" * 60)
    logger.info("\nMonitor: https://code.earthengine.google.com/tasks")
    logger.info("Once complete, run:")
    if args.test:
        logger.info("  python brazil_step1_sample_points.py --task C --test")
    else:
        logger.info("  python brazil_step1_sample_points.py --task C")

    task_info = {
        'task_id':    task.id,
        'asset_id':   asset_id,
        'distance_m': dist_m,
        'test_mode':  args.test,
    }
    task_file = config.OUTPUT_DIR / f'brazil_step0b_tasks_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
    with open(task_file, 'w') as f:
        json.dump(task_info, f, indent=2)
    logger.info(f"✓ Task info saved: {task_file}")


if __name__ == '__main__':
    main()
