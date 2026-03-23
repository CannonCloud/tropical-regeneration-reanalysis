#!/usr/bin/env python3
"""
Brazil Step 0: Export Binary Mask Assets
=========================================

Exports two binary raster assets used by all downstream sampling steps:

  Mask A — brazil_regrowth_mask_v1
    Binary (1 = natural regrowth), derived from the Hansen-validated Fagan
    mosaic. Pixels where Fagan class == 1 (natural regrowth).

  Mask B — brazil_nonregrowth_mask_v1
    Binary (1 = eligible non-regrowth), pixels satisfying ALL of:
      - NOT forest in 2000 (Hansen treecover2000 >= 30%, i.e. strict forest definition)
      - NOT water / bare / urban / sparse veg (ESA CCI 2000 classes 150, 190, 200, 210)
      - NOT any Fagan pixel (classes 1, 2, 3 — regrowth, plantation, open)
    Clipped to Brazil boundary.

Why pre-export these:
    Each sampling batch in Steps A and B calls randomPoints() on a mask.
    If the mask is a live expression (multiple layers combined on the fly),
    GEE recomputes it for every batch. Pre-baking as assets means each
    batch loads a single binary image — much faster and avoids redundant
    computation across batches.

Usage:
    python brazil_step0_export_masks.py           # export both masks
    python brazil_step0_export_masks.py --mask A  # regrowth mask only
    python brazil_step0_export_masks.py --mask B  # non-regrowth mask only

After these tasks complete (~10-20 min each), verify both masks visually
in the GEE code editor before proceeding to Steps A and B.

Assets:
    Fagan validated mosaic:  projects/# replace with your GEE project ID/assets/brazil_plantations_v2_0_hansen_validated
      Band b1: 1=natural regrowth, 2=plantation, 3=open
    ESA CCI 2000:            projects/# replace with your GEE project ID/assets/lc_300m_2000
      Band b1: uint8, values 10-220
    Hansen:                  UMD/hansen/global_forest_change_2021_v1_9
      treecover2000: canopy cover % in year 2000
"""

import ee
import argparse
import logging
import json
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    PROJECT_ID = '# replace with your GEE project ID'

    # ── Input assets ──────────────────────────────────────────────────────────
    FAGAN_ASSET      = 'projects/# replace with your GEE project ID/assets/brazil_plantations_v2_0_hansen_validated'
    FAGAN_BAND       = 'b1'
    FAGAN_ALL_VALUES = [1, 2, 3]       # all classes excluded from non-regrowth domain
    FAGAN_REGROWTH   = 1               # natural regrowth class

    LC_2000_ASSET    = 'projects/# replace with your GEE project ID/assets/lc_300m_2000'
    LC_2000_BAND     = 'b1'
    LC_EXCLUDE       = [150, 190, 200, 210]  # sparse, urban, bare, water

    HANSEN_ASSET     = 'UMD/hansen/global_forest_change_2021_v1_9'
    FOREST_THRESHOLD = 30              # treecover2000 >= this = forest (strict, uses .lt())

    # Hansen grid parameters — all layers are reprojected to this grid before combining.
    # This ensures pixel-perfect alignment and the output asset is anchored to Hansen's grid.
    # Origin: (-180, 80), pixel size: 0.00025 degrees (~25m at equator)
    HANSEN_CRS_TRANSFORM = [0.00025, 0, -180, 0, -0.00025, 80]

    # ── Output assets ─────────────────────────────────────────────────────────
    ASSET_REGROWTH    = 'brazil_regrowth_mask_v2'
    ASSET_NONREGROWTH = 'brazil_nonregrowth_mask_v3'

    # ── Export parameters ─────────────────────────────────────────────────────
    SCALE      = 30        # export at 30m (Fagan/Hansen native resolution)
                           # Step A/B will sample at this resolution too
    MAX_PIXELS = 1e13

    # ── Test mode ─────────────────────────────────────────────────────────────
    # A small patch around Manaus for quick validation before full Brazil run.
    # Run with --test to export just this region (~minutes instead of hours).
    # Manaus centre: -60.025, -3.117
    TEST_MODE   = False
    TEST_BOUNDS = [-60.5, -3.6, -59.5, -2.6]  # ~100km box around Manaus

    # ── Output ────────────────────────────────────────────────────────────────
    OUTPUT_DIR = Path('./gee_outputs')
    LOG_FILE   = OUTPUT_DIR / 'brazil_step0.log'


# =============================================================================
# BRAZIL BOUNDARY
# =============================================================================

def get_brazil_geometry(test_mode=False, test_bounds=None):
    if test_mode:
        # Small box around Manaus for fast validation
        b = test_bounds
        return ee.Geometry.Rectangle([b[0], b[1], b[2], b[3]])
    return (ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
            .filter(ee.Filter.eq('country_na', 'Brazil'))
            .geometry())


# =============================================================================
# MASK BUILDERS
# =============================================================================

def build_regrowth_mask(config):
    """
    Binary mask: 1 where Fagan class == 1 (natural regrowth), 0 elsewhere.
    The Fagan asset is already Hansen-validated so no further intersection needed.
    """
    fagan = ee.Image(config.FAGAN_ASSET)
    mask  = fagan.select(config.FAGAN_BAND).eq(config.FAGAN_REGROWTH)
    return mask.rename('regrowth').toByte()


def build_nonregrowth_mask(config, brazil_geom):
    """
    Binary mask: 1 where pixel is eligible for non-regrowth sampling.

    Eligible = inside Brazil AND not forest AND not excluded LC class AND not any Fagan pixel.

    Key implementation notes:
    (1) All inputs are reprojected to the Hansen grid (crsTransform) before combining.
        Without this, subtle pixel grid misalignments between Fagan (~27m) and Hansen
        (25m) cause pixels at layer boundaries to be compared against the wrong
        neighboring pixel, producing ghost overlaps and missed exclusions.
    (2) All inputs are unmasked to 0 within Brazil so that every pixel has an explicit
        value. Without unmask(0), NoData outside Fagan/ESA footprints propagates through
        And()/Not() operations, producing NoData instead of correct 1s.
    (3) Forest threshold is strict: treecover2000 >= 30% = forest (excluded).
        Using .lt() not .lte() — pixels at exactly 30% are forest by standard definition.
    (4) Export uses explicit crsTransform to lock output to Hansen grid, preventing
        GEE from defaulting to Fagan's grid origin for the output asset.
    """
    hansen = ee.Image(config.HANSEN_ASSET)
    lc2000 = ee.Image(config.LC_2000_ASSET)
    fagan  = ee.Image(config.FAGAN_ASSET)

    # Reproject all inputs to Hansen grid before any operations.
    # This ensures pixel-perfect alignment across all layers.
    hansen_crs = {'crs': 'EPSG:4326', 'crsTransform': config.HANSEN_CRS_TRANSFORM}

    tc2000     = (hansen.select('treecover2000')
                  .reproject(**hansen_crs)
                  .clip(brazil_geom)
                  .unmask(0))
    lc_band    = (lc2000.select(config.LC_2000_BAND)
                  .reproject(**hansen_crs)
                  .clip(brazil_geom)
                  .unmask(0))
    fagan_band = (fagan.select(config.FAGAN_BAND)
                  .reproject(**hansen_crs)
                  .clip(brazil_geom)
                  .unmask(0))

    # Not forest in 2000 — strict definition: treecover2000 >= 30% is forest
    not_forest = tc2000.lt(config.FOREST_THRESHOLD)

    # Not excluded ESA CCI land cover class (150=sparse, 190=urban, 200=bare, 210=water)
    lc_exclude = ee.Image(0)
    for cls in config.LC_EXCLUDE:
        lc_exclude = lc_exclude.Or(lc_band.eq(cls))
    not_lc_excluded = lc_exclude.Not()

    # Not any Fagan pixel (regrowth=1, plantation=2, open=3)
    fagan_any = ee.Image(0)
    for val in config.FAGAN_ALL_VALUES:
        fagan_any = fagan_any.Or(fagan_band.eq(val))
    not_fagan = fagan_any.Not()

    # Combine: eligible where all conditions true, then mask to Brazil
    eligible = (not_forest
                .And(not_lc_excluded)
                .And(not_fagan)
                .clip(brazil_geom)
                .selfMask())

    return eligible.rename('nonregrowth').toByte()


# =============================================================================
# EXPORT HELPER
# =============================================================================

def export_mask(image, asset_name, description, brazil_geom, config, logger):
    """Export a binary mask image to a GEE asset."""
    asset_id = f'projects/{config.PROJECT_ID}/assets/{asset_name}'

    task = ee.batch.Export.image.toAsset(
        image            = image,
        description      = description,
        assetId          = asset_id,
        region           = brazil_geom,
        crs              = 'EPSG:4326',
        crsTransform     = config.HANSEN_CRS_TRANSFORM,  # lock output to Hansen grid
        maxPixels        = config.MAX_PIXELS,
        pyramidingPolicy = {'.default': 'mode'}  # correct for binary/categorical
    )
    task.start()

    info = {
        'task_id':    task.id,
        'asset_id':   asset_id,
        'asset_name': asset_name,
        'description': description,
    }
    logger.info(f"  ✓ Submitted: {description}")
    logger.info(f"    Task ID:  {task.id}")
    logger.info(f"    Asset:    {asset_id}")
    return info


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Brazil Step 0: Export binary mask assets')
    parser.add_argument('--mask', choices=['A', 'B', 'both'], default='both',
                        help='Which mask to export: A=regrowth, B=non-regrowth, both (default)')
    parser.add_argument('--test', action='store_true',
                        help='Run on small Manaus test region instead of all Brazil (~minutes not hours)')
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

    logger.info("=" * 60)
    logger.info("BRAZIL STEP 0: EXPORT BINARY MASKS")
    logger.info("=" * 60)

    # ── Initialize GEE ────────────────────────────────────────────────────────
    try:
        ee.Initialize(project=config.PROJECT_ID)
        logger.info(f"✓ GEE initialized: {config.PROJECT_ID}")
    except Exception:
        logger.info("Authenticating with GEE...")
        ee.Authenticate()
        ee.Initialize(project=config.PROJECT_ID)
        logger.info("✓ Authentication successful")

    config.TEST_MODE = args.test
    brazil_geom = get_brazil_geometry(test_mode=args.test, test_bounds=config.TEST_BOUNDS)
    if args.test:
        logger.info("✓ TEST MODE: using Manaus region (~100km box)")
        # Use distinct asset names for test outputs so they don't overwrite real masks
        config.ASSET_REGROWTH    = 'brazil_regrowth_mask_test'
        config.ASSET_NONREGROWTH = 'brazil_nonregrowth_mask_test'
    else:
        logger.info("✓ Brazil boundary loaded (full country)")

    all_task_infos = []

    # ── Mask A: Regrowth ──────────────────────────────────────────────────────
    if args.mask in ('A', 'both'):
        logger.info("\n── Mask A: Regrowth ──")
        logger.info(f"  Source: {config.FAGAN_ASSET}")
        logger.info(f"  Logic:  b1 == 1 (natural regrowth)")

        mask_a = build_regrowth_mask(config)
        info   = export_mask(
            image       = mask_a,
            asset_name  = config.ASSET_REGROWTH,
            description = 'BrazilRegrowthMask',
            brazil_geom = brazil_geom,
            config      = config,
            logger      = logger
        )
        all_task_infos.append(info)

    # ── Mask B: Non-regrowth eligible ─────────────────────────────────────────
    if args.mask in ('B', 'both'):
        logger.info("\n── Mask B: Non-regrowth eligible ──")
        logger.info(f"  Hansen:  treecover2000 >= {config.FOREST_THRESHOLD}% excluded (strict forest definition)")
        logger.info(f"  ESA CCI: exclude classes {config.LC_EXCLUDE}")
        logger.info(f"  Fagan:   exclude all classes {config.FAGAN_ALL_VALUES}")

        mask_b = build_nonregrowth_mask(config, brazil_geom)
        info   = export_mask(
            image       = mask_b,
            asset_name  = config.ASSET_NONREGROWTH,
            description = 'BrazilNonRegrowthMask',
            brazil_geom = brazil_geom,
            config      = config,
            logger      = logger
        )
        all_task_infos.append(info)

    # ── Summary ───────────────────────────────────────────────────────────────
    logger.info("\n" + "=" * 60)
    logger.info("TASKS SUBMITTED")
    logger.info("=" * 60)
    for info in all_task_infos:
        logger.info(f"  {info['description']:<30} | Task: {info['task_id']}")
        logger.info(f"    → {info['asset_id']}")

    logger.info("\nExpected runtime: 10-20 min per mask")
    logger.info("Monitor: https://code.earthengine.google.com/tasks")
    logger.info("\nBefore running Step 1, verify masks in GEE code editor:")
    logger.info("  - Regrowth mask should show green patches across Amazon/Atlantic Forest")
    logger.info("  - Non-regrowth mask should be absent over dense forest, water, urban areas")
    logger.info("  - Check boundaries look clean at Brazil border")
    logger.info("\nOnce verified, proceed to:")
    logger.info("  python brazil_step1_sample_points.py --task A")
    logger.info("  python brazil_step1_sample_points.py --task B")

    task_file = config.OUTPUT_DIR / f'brazil_step0_tasks_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
    with open(task_file, 'w') as f:
        json.dump(all_task_infos, f, indent=2)
    logger.info(f"\n✓ Task info saved: {task_file}")


if __name__ == '__main__':
    main()
