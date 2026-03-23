#!/usr/bin/env python3
"""
Brazil Step 1: Sample Points
==============================

Samples points using pre-baked binary mask assets from Step 0.
Uses randomPoints() — no tileScale, so sampling is truly random
across the full Brazil extent rather than tile-by-tile.

Tasks:
  A — Regrowth points (Y=1): from brazil_regrowth_mask_v2
  B — Non-regrowth points (Y=0): from brazil_nonregrowth_mask_v6
  C — Near-regrowth annulus points (Y=0): from brazil_nonregrowth_annulus_0_3km_v2
        Batched (4 × 50K) with tileScale=8 to handle the thin/complex mask at 30m.
        Spatial distance-band pairing (30m–1km, 1–5km, 5–25km, 25–200km) is done in R.

Usage:
    python brazil_step1_sample_points.py --task A                      # regrowth only
    python brazil_step1_sample_points.py --task B                      # non-regrowth, start at batch 0
    python brazil_step1_sample_points.py --task B --batch-start 2      # start at batch 2 (avoids overwrite)
    python brazil_step1_sample_points.py --task C                      # annulus points (4 batches of 50K)
    python brazil_step1_sample_points.py --task C --test               # test on Manaus box, 1K points
    python brazil_step1_sample_points.py --task all                    # A + B + C

Assets used:
    brazil_regrowth_mask_v2              — binary, biome + 25°S lat filtered
    brazil_nonregrowth_mask_v6           — binary, biome + 25°S lat filtered
    brazil_nonregrowth_annulus_0_3km_v2  — binary, nonregrowth within 3km of regrowth, biome filtered
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

    # ── Input mask assets (from Step 0) ───────────────────────────────────────
    REGROWTH_MASK_ASSET    = 'projects/# replace with your GEE project ID/assets/brazil_regrowth_mask_v2'
    NONREGROWTH_MASK_ASSET = 'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_mask_v6'

    # ── Sampling targets ──────────────────────────────────────────────────────
    N_REGROWTH    = 100_000
    N_NONREGROWTH = 200_000
    BATCH_SIZE    = 200_000

    # Seeds: batch i uses SEED_BASE + i * 10
    SEED_BASE_A = 42
    SEED_BASE_B = 100

    # ── Output asset names ────────────────────────────────────────────────────
    ASSET_A_BATCH = 'brazil_regrowth_batch_{i:02d}'
    ASSET_B_BATCH = 'brazil_nonregrowth_batch_{i:02d}'

    # ── Annulus sampling (Task C) ─────────────────────────────────────────────
    # Near-regrowth mask (0-3km) built by step0b, biome filtered.
    # Batched into 50K chunks — the thin/complex geometry of this mask causes
    # OOM errors at 200K in a single task at 30m scale. tileScale=8 subdivides
    # the sampling grid to avoid holding the full mask in memory at once.
    # Spatial matching into distance bands (30m–1km, 1–5km, 5–25km, 25–200km)
    # is done post-hoc in R using sf nearest-neighbour matching.
    ANNULUS_MASK_ASSET  = 'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_annulus_0_3km_v2'
    N_ANNULUS           = 200_000
    ANNULUS_BATCH_SIZE  = 200_000    # 4 batches of 50K
    SEED_BASE_C         = 300
    ASSET_C_BATCH       = 'brazil_nonregrowth_annulus_batch_{i:02d}'

    # ── Test mode ─────────────────────────────────────────────────────────────
    # Small box around Manaus — fast validation before full Brazil run
    TEST_BOUNDS    = [-60.5, -3.6, -59.5, -2.6]
    N_ANNULUS_TEST = 1_000   # points in test mode

    # ── Output ────────────────────────────────────────────────────────────────
    OUTPUT_DIR = Path('./gee_outputs')
    LOG_FILE   = OUTPUT_DIR / 'brazil_step1.log'


# =============================================================================
# HELPERS
# =============================================================================

def get_brazil_geometry():
    return (ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
            .filter(ee.Filter.eq('country_na', 'Brazil'))
            .geometry())


def add_metadata(fc, sample_type, batch_id, seed):
    """Add point_id, sample_type, batch_id, seed to each feature."""
    prefix = f"{sample_type}_b{batch_id:02d}_"

    def _add(feature):
        return feature.set({
            'sample_type': sample_type,
            'point_id':    ee.String(prefix).cat(feature.id()),
            'batch_id':    batch_id,
            'seed':        seed,
        })
    return fc.map(_add)


def submit_export(fc, asset_name, description, project_id, logger):
    """Submit Export.table.toAsset and return task info dict."""
    asset_id = f'projects/{project_id}/assets/{asset_name}'
    task = ee.batch.Export.table.toAsset(
        collection  = fc,
        description = description,
        assetId     = asset_id
    )
    task.start()
    info = {
        'task_id':     task.id,
        'asset_id':    asset_id,
        'asset_name':  asset_name,
        'description': description,
    }
    logger.info(f"    ✓ Task: {task.id} → {asset_name}")
    return info


# =============================================================================
# SAMPLING
# =============================================================================

def run_batches(logger, config, mask_asset, sample_type, n_total,
                seed_base, asset_template, brazil_geom, batch_start=0, scale=30):
    """
    Sample n_total points from mask_asset in batches of BATCH_SIZE.
    Each batch uses a different seed for reproducibility.
    Uses stratifiedSample without tileScale for true spatial randomness.
    batch_start: offset the batch index to avoid overwriting existing assets.
    """
    n_batches = (n_total + config.BATCH_SIZE - 1) // config.BATCH_SIZE + 1
    logger.info(f"  Target:      {n_total:,}")
    logger.info(f"  Batches:     {n_batches} x {config.BATCH_SIZE:,}")
    logger.info(f"  Batch start: {batch_start}")
    logger.info(f"  Mask:        {mask_asset}")

    mask = ee.Image(mask_asset).rename('mask').toInt()

    task_infos = []
    for i in range(n_batches):
        batch_i = i + batch_start
        seed = seed_base + batch_i * 10
        logger.info(f"\n  Batch {batch_i+1} (loop {i+1}/{n_batches}) | seed={seed}")

        points = mask.stratifiedSample(
            numPoints   = 0,
            classBand   = 'mask',
            region      = brazil_geom,
            scale       = scale,
            seed        = seed,
            geometries  = True,
            dropNulls   = True,
            classValues = [1],
            classPoints = [config.BATCH_SIZE]
            # No tileScale — samples randomly across the full region
        )

        points = add_metadata(points, sample_type, batch_id=batch_i, seed=seed)

        asset_name = asset_template.format(i=batch_i)
        info = submit_export(
            fc          = points,
            asset_name  = asset_name,
            description = f'Brazil_{sample_type}_batch{batch_i:02d}',
            project_id  = config.PROJECT_ID,
            logger      = logger
        )
        task_infos.append(info)
        time.sleep(1.5)   # avoid API rate limits

    return task_infos


# =============================================================================
# ANNULUS SAMPLING (Task C)
# =============================================================================

def run_annulus_bands(logger, config, brazil_geom, test_mode=False):
    """
    Sample non-regrowth points from the 0-3km near-regrowth annulus mask
    in batches of ANNULUS_BATCH_SIZE.

    Why batched:
        The annulus mask is geometrically thin and complex (narrow corridors
        of eligible pixels threading through the landscape). A single 200K
        stratifiedSample causes OOM/timeout errors. Splitting into 50K
        batches keeps each task within GEE memory limits.

    Why .bounds() instead of brazil_geom:
        At 90m scale, GEE still has to check every pixel against the region
        polygon. brazil_geom (from LSIB) has tens of thousands of vertices,
        making this a very expensive point-in-polygon check. Using
        brazil_geom.bounds() reduces this to a trivial min/max coordinate
        check. The annulus mask was already clipped to Brazil during Step 0b,
        so no points outside Brazil can be sampled regardless — the complex
        border is pure overhead here.

    Why scale=90 not 30:
        Distribution bias from coarser scale does not matter here because
        these points are paired to individual regrowth points by distance
        band in R. A 90m grid still resolves the 30m-1km inner band
        adequately and reduces pixel count 9x vs 30m.

    Why tileScale=4 not 8:
        With .bounds() removing geometric complexity, tileScale=8 adds
        unnecessary scheduling overhead (64 sub-tiles vs 16). tileScale=4
        is sufficient for the remaining memory pressure from the sparse mask.
    """
    if test_mode:
        n_total      = config.N_ANNULUS_TEST
        asset_suffix = '_test'
        logger.info(f"  TEST MODE: {n_total:,} points (1 batch)")
    else:
        n_total      = config.N_ANNULUS
        asset_suffix = ''

    n_batches     = (n_total + config.ANNULUS_BATCH_SIZE - 1) // config.ANNULUS_BATCH_SIZE
    batch_size    = config.N_ANNULUS_TEST if test_mode else config.ANNULUS_BATCH_SIZE

    # Bounding box — avoids expensive point-in-polygon checks against the
    # complex LSIB border. The mask asset handles the actual spatial filtering.
    simple_region = brazil_geom.bounds()

    logger.info(f"  Mask:      {config.ANNULUS_MASK_ASSET}")
    logger.info(f"  Target:    {n_total:,} points")
    logger.info(f"  Batches:   {n_batches} x {batch_size:,}")
    logger.info(f"  Scale:     90m | tileScale=4 | region=bounds()")

    annulus_mask = ee.Image(config.ANNULUS_MASK_ASSET).rename('mask').toInt()
    task_infos   = []

    for i in range(n_batches):
        seed = config.SEED_BASE_C + i * 10
        logger.info(f"\n  Batch {i+1}/{n_batches} | seed={seed}")

        points = annulus_mask.stratifiedSample(
            numPoints   = 0,
            classBand   = 'mask',
            region      = simple_region,  # bounding box — fast rejection, mask handles eligibility
            scale       = 90,
            seed        = seed,
            geometries  = True,
            dropNulls   = True,
            classValues = [1],
            classPoints = [batch_size],
            tileScale   = 1
        )

        points = add_metadata(points, 'nonregrowth_annulus_0_3km', batch_id=i, seed=seed)

        asset_name = config.ASSET_C_BATCH.format(i=i)
        if test_mode:
            asset_name = asset_name.replace('annulus_batch', 'annulus_batch_test')

        info = submit_export(
            fc          = points,
            asset_name  = asset_name,
            description = f'Brazil_annulus_batch{i:02d}',
            project_id  = config.PROJECT_ID,
            logger      = logger
        )
        task_infos.append(info)
        time.sleep(1.5)

    return task_infos


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Brazil Step 1: Sample Points')
    parser.add_argument('--task', choices=['A', 'B', 'C', 'all'], default='A',
                        help='A=regrowth, B=non-regrowth, C=annulus bands, all=A+B+C (default: A)')
    parser.add_argument('--batch-start', type=int, default=0,
                        help='Batch index to start from (default: 0). Use to avoid overwriting existing batches.')
    parser.add_argument('--test', action='store_true',
                        help='Task C only: run on small Manaus box with 1K points per band (~minutes)')
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
    logger.info("BRAZIL STEP 1: SAMPLE POINTS")
    logger.info("=" * 60)
    logger.info(f"Task(s):    {args.task}")
    logger.info(f"Batch size: {config.BATCH_SIZE:,}")

    # ── Initialize GEE ────────────────────────────────────────────────────────
    try:
        ee.Initialize(project=config.PROJECT_ID)
        logger.info(f"✓ GEE initialized: {config.PROJECT_ID}")
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=config.PROJECT_ID)
        logger.info("✓ GEE authenticated and initialized")

    brazil_geom = get_brazil_geometry()
    logger.info("✓ Brazil boundary loaded")

    all_task_infos = []

    batch_start = args.batch_start

    # ── Task A: Regrowth ──────────────────────────────────────────────────────
    if args.task in ('A', 'all'):
        logger.info(f"\n── Task A: Regrowth points ──")
        infos = run_batches(
            logger         = logger,
            config         = config,
            mask_asset     = config.REGROWTH_MASK_ASSET,
            sample_type    = 'regrowth',
            n_total        = config.N_REGROWTH,
            seed_base      = config.SEED_BASE_A,
            asset_template = config.ASSET_A_BATCH,
            brazil_geom    = brazil_geom,
            batch_start    = batch_start
        )
        all_task_infos.extend(infos)

    # ── Task B: Non-regrowth ──────────────────────────────────────────────────
    if args.task in ('B', 'all'):
        logger.info(f"\n── Task B: Non-regrowth points ──")
        infos = run_batches(
            logger         = logger,
            config         = config,
            mask_asset     = config.NONREGROWTH_MASK_ASSET,
            sample_type    = 'nonregrowth_brazil',
            n_total        = config.N_NONREGROWTH,
            seed_base      = config.SEED_BASE_B,
            asset_template = config.ASSET_B_BATCH,
            brazil_geom    = brazil_geom,
            batch_start    = batch_start,
            scale          = 120   # coarser scale for dense mask — much faster
        )
        all_task_infos.extend(infos)

    # ── Task C: Annulus bands ─────────────────────────────────────────────────
    if args.task in ('C', 'all'):
        logger.info(f"\n── Task C: Annulus points (4 × 50K batches) ──")
        if args.test:
            logger.info("  TEST MODE: Manaus box, 1K points")
        infos = run_annulus_bands(
            logger      = logger,
            config      = config,
            brazil_geom = brazil_geom,
            test_mode   = args.test,
        )
        all_task_infos.extend(infos)

    # ── Summary ───────────────────────────────────────────────────────────────
    logger.info("\n" + "=" * 60)
    logger.info("ALL TASKS SUBMITTED")
    logger.info("=" * 60)
    for info in all_task_infos:
        logger.info(f"  {info['description']:<40} | {info['task_id']}")
        logger.info(f"    → {info['asset_id']}")

    logger.info("\nMonitor: https://code.earthengine.google.com/tasks")
    logger.info("Expected: ~10-20 min per batch (batches run in parallel)")

    task_file = config.OUTPUT_DIR / f'brazil_step1_tasks_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
    with open(task_file, 'w') as f:
        json.dump(all_task_infos, f, indent=2)
    logger.info(f"\n✓ Task info saved: {task_file}")


if __name__ == '__main__':
    main()
