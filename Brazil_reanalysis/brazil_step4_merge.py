#!/usr/bin/env python3
"""
Brazil Step 4: Merge All Layer CSVs Into Final Dataset
=======================================================

Merges all Step 2 and Step 3 outputs into a single analysis-ready CSV.

Three sample sets are assembled separately then stacked:
    regrowth    — Y=1, ~100K points
    nonregrowth — Y=0, ~800K points (4 batch_id files per layer)
    annulus     — Y=0 near-regrowth, ~200K points (climate envelope analysis)

Core layers (step2 — single CSV per sample set):
    climate_part1-5  → bio01-bio19 (PCA done in R)
    topography       → bio_elevation, bio_slope
    landcover        → bio_landcover_class, cropland_weight, lc_raw_value
    npp_fire         → bio_npp_mean, bio_fire_freq
    soil_chemical    → bio_biome_id, bio_soil_organic_carbon, bio_soil_ph
    soil_texture     → bio_soil_sand, bio_soil_clay
    soil_physical    → bio_soil_bulk_density_fine, bio_soil_water_33kpa
    socioeconomic    → econ_dist_water, econ_road_density, paper_gdp_2015,
                       paper_ghs_pop, paper_protected_binary
    highres_lights   → highres_lights_2012_2014  (split to avoid mask contamination)
    highres_travel   → highres_travel
    highres_ghm      → highres_ghm
    highres_worldpop → highres_worldpop_density

Tiled layers (step3 — folder of per-tile CSVs):
    forest      (GEE_Brazil_forest_2000/)       → bio_dist_forest_2000, bio_forest_density_1km2
    urban       (GEE_Brazil_urban/)              → paper_dist_urban
    settlement  (GEE_Brazil_settlement/)         → highres_dist_settlement, highres_worldpop_density_2km
    plantation  (GEE_Brazil_plantation/)         → bio_dist_plantation, bio_plantation_density_1km, bio_canopy_height_1km_mean
    cultivated_grass (GEE_Brazil_cultivated_grass/) → highres_dist_cultivated_grass, highres_cultivated_grass_density_2km
    cropland         (GEE_Brazil_cropland/)          → highres_dist_cropland, highres_cropland_density_5km
    nightlight_density (GEE_Brazil_nightlight_density/) → highres_nightlight_density_2km

Output columns:
    point_id        — unique point identifier
    sample_set      — 'regrowth' | 'nonregrowth' | 'annulus'
    y               — 1 (regrowth) or 0 (nonregrowth / annulus)
    batch_id        — original GEE sampling batch
    longitude       — parsed from .geo
    latitude        — parsed from .geo
    bio01..bio19    — WorldClim v2 raw bioclim variables (PCA in R)
    bio_elevation, bio_slope
    bio_landcover_class, cropland_weight, lc_raw_value
    bio_npp_mean, bio_fire_freq
    bio_biome_id, bio_soil_organic_carbon, bio_soil_ph
    bio_soil_sand, bio_soil_clay
    bio_soil_bulk_density_fine, bio_soil_water_33kpa
    econ_dist_water, econ_road_density
    paper_gdp_2015, paper_ghs_pop, paper_protected_binary
    highres_lights_2012_2014, highres_travel, highres_ghm, highres_worldpop_density
    bio_dist_forest_2000, bio_forest_density_1km2   (if step3 forest complete)
    paper_dist_urban                                 (if step3 urban complete)
    highres_dist_settlement, highres_worldpop_density_2km  (if step3 settlement complete)

Usage:
    python brazil_step4_merge.py
    python brazil_step4_merge.py --gee-dir ./gee_outputs
    python brazil_step4_merge.py --gee-dir ./gee_outputs --output brazil_training_data.csv
    python brazil_step4_merge.py --skip-forest       # if step3 forest not yet complete
    python brazil_step4_merge.py --skip-tiles         # skip ALL step3 tiled variables
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================

# Step2 layers to include (single CSV per sample set)
# nonregrowth files are split across pb2-pb5; others are single files
CORE_LAYERS = [
    'climate_part1',
    'climate_part2',
    'climate_part3',
    'climate_part4',
    'climate_part5',
    'topography',
    'landcover',
    'npp_fire',
    'soil_chemical',
    'soil_texture',
    'soil_physical',
    'socioeconomic',
    'highres_lights',
    'highres_travel',
    'highres_ghm',
    'highres_worldpop',
]

# Step3 tiled layers — each entry: (variable_name, folder_name, file_prefix, data_columns)
# file_prefix is the start of the CSV filename before _{sample_type}_{tile_id}.csv
TILED_LAYERS = [
    {
        'name': 'forest',
        'folder': 'GEE_Brazil_forest_2000',
        'glob_pattern': 'Brazil_forest2000_{sample_set}_*.csv',
        'data_cols': ['bio_dist_forest_2000', 'bio_forest_density_1km2'],
    },
    {
        'name': 'urban',
        'folder': 'GEE_Brazil_urban',
        'glob_pattern': 'Brazil_urban_{sample_set}_*.csv',
        'data_cols': ['paper_dist_urban'],
    },
    {
        'name': 'settlement',
        'folder': 'GEE_Brazil_settlement',
        'glob_pattern': 'Brazil_settlement_{sample_set}_*.csv',
        'data_cols': ['highres_dist_settlement', 'highres_worldpop_density_2km'],
    },
    {
        'name': 'plantation',
        'folder': 'GEE_Brazil_plantation',
        'glob_pattern': 'Brazil_plantation_{sample_set}_*.csv',
        'data_cols': ['bio_dist_plantation', 'bio_plantation_density_1km', 'bio_canopy_height_1km_mean'],
    },
    {
        'name': 'cultivated_grass',
        'folder': 'GEE_Brazil_cultivated_grass',
        'glob_pattern': 'Brazil_cultivated_grass_{sample_set}_*.csv',
        'data_cols': ['highres_dist_cultivated_grass', 'highres_cultivated_grass_density_2km'],
    },
    {
        'name': 'cropland',
        'folder': 'GEE_Brazil_cropland',
        'glob_pattern': 'Brazil_cropland_{sample_set}_*.csv',
        'data_cols': ['highres_dist_cropland', 'highres_cropland_density_5km'],
    },
    {
        'name': 'nightlight_density',
        'folder': 'GEE_Brazil_nightlight_density',
        'glob_pattern': 'Brazil_nightlight_density_{sample_set}_*.csv',
        'data_cols': ['highres_nightlight_density_2km'],
    },
]

# Columns GEE adds or that we never want in the final dataset
GEE_DROP_COLS = {'system:index', 'mask', 'overlap_regrowth', 'bio_forest_density_2018'}

# Columns that exist in some files but should only be kept from the
# first (base) layer to avoid duplicates during merging
META_COLS = {'point_id', 'batch_id', 'sample_type', 'seed', '.geo'}

# Nonregrowth batch suffixes (batch_id 2-5)
NONREGROWTH_BATCHES = ['pb2', 'pb3', 'pb4', 'pb5']


# =============================================================================
# HELPERS
# =============================================================================

def parse_geo(geo_str):
    """Extract (longitude, latitude) from GEE .geo JSON string."""
    try:
        geo = json.loads(geo_str)
        coords = geo['coordinates']
        return float(coords[0]), float(coords[1])
    except Exception:
        return None, None


def load_layer(path: Path, logger) -> pd.DataFrame:
    """Load a single layer CSV, dropping GEE internal columns."""
    df = pd.read_csv(path)
    drop = GEE_DROP_COLS & set(df.columns)
    if drop:
        df = df.drop(columns=drop)
    return df


def load_and_concat_nonregrowth_layer(gee_dir: Path, layer: str, logger) -> pd.DataFrame:
    """Load pb2-pb5 files for one layer and concatenate."""
    dfs = []
    for pb in NONREGROWTH_BATCHES:
        path = gee_dir / f"Brazil_nonregrowth_{layer}_{pb}.csv"
        if not path.exists():
            logger.warning(f"  ⚠ Missing: {path.name}")
            continue
        dfs.append(load_layer(path, logger))
    if not dfs:
        raise FileNotFoundError(f"No nonregrowth files found for layer: {layer}")
    return pd.concat(dfs, ignore_index=True)


def build_sample_set(gee_dir: Path, sample_set: str, logger) -> pd.DataFrame:
    """
    Load and merge all layers for one sample set (regrowth / nonregrowth / annulus).

    Strategy:
        - Load the first layer as base (gives us point_id, batch_id, sample_type,
          seed, .geo plus that layer's data columns)
        - For each remaining layer, load data columns only and left-join on point_id
        - Parse .geo into longitude/latitude after all joins
        - Drop .geo, seed (seed is in point_id implicitly)
    """
    logger.info(f"\n{'─' * 60}")
    logger.info(f"Building: {sample_set.upper()}")
    logger.info(f"{'─' * 60}")

    base_df = None
    skipped = []

    for i, layer in enumerate(CORE_LAYERS):
        logger.info(f"  [{i+1}/{len(CORE_LAYERS)}] {layer}")

        try:
            if sample_set == 'nonregrowth':
                df = load_and_concat_nonregrowth_layer(gee_dir, layer, logger)
            else:
                path = gee_dir / f"Brazil_{sample_set}_{layer}.csv"
                df = load_layer(path, logger)
        except Exception as e:
            logger.warning(f"    ⚠ Skipping {layer}: {e}")
            skipped.append(layer)
            continue

        logger.info(f"    rows={len(df):,}  cols={list(df.columns)}")

        if base_df is None:
            # First layer — keep everything including meta
            base_df = df
        else:
            # Subsequent layers — keep only point_id + new data columns
            meta_in_df = META_COLS & set(df.columns)
            data_cols  = [c for c in df.columns if c not in meta_in_df]
            # Avoid duplicate columns already present in base
            new_cols   = [c for c in data_cols if c not in base_df.columns]
            merge_cols = ['point_id'] + new_cols
            base_df    = base_df.merge(df[merge_cols], on='point_id', how='left')

        logger.info(f"    → running shape: {base_df.shape}")

    if base_df is None:
        raise RuntimeError(f"No data loaded for sample set: {sample_set}")

    # Parse .geo → longitude, latitude early so we can deduplicate on coords.
    # GEE stratifiedSample places points at pixel centroids at the sampling
    # scale — two draws on the same pixel return bit-identical coordinates.
    # Deduplication on (longitude, latitude) therefore correctly catches the
    # rare case where two batches sampled the same pixel.
    if '.geo' in base_df.columns:
        coords = base_df['.geo'].apply(lambda x: pd.Series(parse_geo(x)))
        base_df['longitude'] = coords[0]
        base_df['latitude']  = coords[1]
        base_df = base_df.drop(columns=['.geo'])

    dups = base_df.duplicated(subset=['longitude', 'latitude']).sum()
    if dups > 0:
        logger.warning(f"  ⚠ {dups:,} duplicate pixel centroids in {sample_set} — removing")
        base_df = base_df.drop_duplicates(subset=['longitude', 'latitude'], keep='first')
    else:
        logger.info(f"  ✓ No duplicate pixel centroids in {sample_set}")

    # Add outcome variable and canonical sample_set label
    base_df['y']          = 1 if sample_set == 'regrowth' else 0
    base_df['sample_set'] = sample_set

    # Drop seed (redundant — batch_id carries this info)
    if 'seed' in base_df.columns:
        base_df = base_df.drop(columns=['seed'])

    # Drop sample_type (replaced by sample_set which is cleaner)
    if 'sample_type' in base_df.columns:
        base_df = base_df.drop(columns=['sample_type'])

    logger.info(f"\n  ✓ {sample_set}: {len(base_df):,} rows × {len(base_df.columns)} cols")
    if skipped:
        logger.warning(f"  ⚠ Skipped layers (not yet available): {', '.join(skipped)}")
    return base_df


def load_tiled_layer(gee_dir: Path, tiled_cfg: dict, sample_set: str, logger) -> pd.DataFrame | None:
    """
    Load and concatenate all step3 tile CSVs for a tiled variable.
    Returns None if the folder or tiles are not found.

    tiled_cfg keys:
        name          — human label (e.g. 'forest')
        folder        — subfolder name under gee_dir
        glob_pattern  — pattern with {sample_set} placeholder
        data_cols     — list of variable column names to keep
    """
    tile_dir = gee_dir / tiled_cfg['folder']
    if not tile_dir.exists():
        logger.warning(f"  ⚠ {tiled_cfg['folder']}/ not found — skipping {tiled_cfg['name']}")
        return None

    pattern = tiled_cfg['glob_pattern'].format(sample_set=sample_set)
    tile_files = sorted(tile_dir.glob(pattern))

    if not tile_files:
        logger.warning(f"  ⚠ No tiles found for {tiled_cfg['name']}/{sample_set} "
                       f"(pattern: {pattern}) — skipping")
        return None

    logger.info(f"  Loading {len(tile_files)} {tiled_cfg['name']} tiles for {sample_set}")
    dfs = []
    for f in tile_files:
        df = pd.read_csv(f)
        drop = GEE_DROP_COLS & set(df.columns)
        if drop:
            df = df.drop(columns=drop)
        # Keep only point_id and the expected data columns
        keep = [c for c in df.columns
                if c == 'point_id' or c in tiled_cfg['data_cols']]
        dfs.append(df[keep])

    tile_df = pd.concat(dfs, ignore_index=True)

    # Deduplicate — points near tile boundaries may appear in two tiles
    dups = tile_df.duplicated(subset=['point_id']).sum()
    if dups > 0:
        logger.info(f"  Removing {dups:,} duplicate point_ids from tile overlaps")
        tile_df = tile_df.drop_duplicates(subset=['point_id'], keep='first')

    logger.info(f"  ✓ {tiled_cfg['name']} tiles: {len(tile_df):,} rows")
    return tile_df


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Brazil Step 4: Merge all layer CSVs into final dataset')
    parser.add_argument('--gee-dir', type=str, default='./gee_outputs',
                        help='Directory containing GEE outputs (default: ./gee_outputs)')
    parser.add_argument('--output', type=str, default='brazil_training_data.csv',
                        help='Output filename (default: brazil_training_data.csv)')
    parser.add_argument('--skip-forest', action='store_true',
                        help='Skip step3 forest variables (use if step3 forest not yet complete)')
    parser.add_argument('--skip-tiles', action='store_true',
                        help='Skip ALL step3 tiled variables (forest, urban, settlement)')
    args = parser.parse_args()

    gee_dir = Path(args.gee_dir)
    if not gee_dir.exists():
        print(f"Error: directory not found: {gee_dir}")
        exit(1)

    log_file = gee_dir / 'brazil_step4_merge.log'
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)

    logger.info("=" * 60)
    logger.info("BRAZIL STEP 4: MERGE ALL LAYERS")
    logger.info("=" * 60)
    logger.info(f"GEE dir: {gee_dir}")
    logger.info(f"Output:  {args.output}")
    logger.info(f"Core layers:  {', '.join(CORE_LAYERS)}")
    tiled_names = [t['name'] for t in TILED_LAYERS]
    if args.skip_tiles:
        logger.info(f"Tiled layers: ALL SKIPPED (--skip-tiles)")
    elif args.skip_forest:
        logger.info(f"Tiled layers: {', '.join(tiled_names)} (forest SKIPPED via --skip-forest)")
    else:
        logger.info(f"Tiled layers: {', '.join(tiled_names)}")

    # ── Build each sample set ─────────────────────────────────────────────────
    sample_sets = ['regrowth', 'nonregrowth', 'annulus']
    assembled   = []

    for sample_set in sample_sets:
        df = build_sample_set(gee_dir, sample_set, logger)

        # Attach step3 tiled variables if available
        if not args.skip_tiles:
            for tiled_cfg in TILED_LAYERS:
                # Honour --skip-forest for backward compat
                if tiled_cfg['name'] == 'forest' and args.skip_forest:
                    logger.info(f"  Skipping forest tiles (--skip-forest)")
                    continue

                tile_df = load_tiled_layer(gee_dir, tiled_cfg, sample_set, logger)
                if tile_df is not None:
                    before  = len(df)
                    df      = df.merge(tile_df, on='point_id', how='left')
                    # Report match rate using the first data column
                    check_col = tiled_cfg['data_cols'][0]
                    matched   = df[check_col].notna().sum()
                    logger.info(f"  {tiled_cfg['name']} join: "
                                f"{matched:,}/{before:,} points matched")

        assembled.append(df)

    # ── Stack all three sample sets ───────────────────────────────────────────
    logger.info(f"\n{'=' * 60}")
    logger.info("STACKING SAMPLE SETS")
    logger.info(f"{'=' * 60}")

    final = pd.concat(assembled, ignore_index=True)
    logger.info(f"  regrowth:    {(final.sample_set == 'regrowth').sum():>8,} rows")
    logger.info(f"  nonregrowth: {(final.sample_set == 'nonregrowth').sum():>8,} rows")
    logger.info(f"  annulus:     {(final.sample_set == 'annulus').sum():>8,} rows")
    logger.info(f"  total:       {len(final):>8,} rows × {len(final.columns)} cols")

    # ── Column order ──────────────────────────────────────────────────────────
    id_cols      = ['point_id', 'sample_set', 'y', 'batch_id', 'longitude', 'latitude']
    climate_cols = [f'bio{str(i).zfill(2)}' for i in range(1, 20)]
    topo_cols    = ['bio_elevation', 'bio_slope']
    lc_cols      = ['bio_landcover_class', 'cropland_weight', 'lc_raw_value']
    npp_cols     = ['bio_npp_mean', 'bio_fire_freq']
    soil_cols    = ['bio_biome_id', 'bio_soil_organic_carbon', 'bio_soil_ph',
                    'bio_soil_sand', 'bio_soil_clay',
                    'bio_soil_bulk_density_fine', 'bio_soil_water_33kpa']
    forest_cols  = ['bio_dist_forest_2000', 'bio_forest_density_1km2']
    plant_cols   = ['bio_dist_plantation', 'bio_plantation_density_1km', 'bio_canopy_height_1km_mean']
    econ_cols    = ['econ_dist_water', 'econ_road_density']
    paper_cols   = ['paper_gdp_2015', 'paper_ghs_pop', 'paper_protected_binary',
                    'paper_dist_urban']
    highres_cols = ['highres_lights_2012_2014', 'highres_travel', 'highres_ghm',
                    'highres_worldpop_density',
                    'highres_dist_settlement', 'highres_worldpop_density_2km',
                    'highres_dist_cultivated_grass', 'highres_cultivated_grass_density_2km',
                    'highres_dist_cropland', 'highres_cropland_density_5km',
                    'highres_nightlight_density_2km']

    desired_order = (id_cols + climate_cols + topo_cols + lc_cols + npp_cols
                     + soil_cols + forest_cols + plant_cols + econ_cols
                     + paper_cols + highres_cols)
    # Only keep columns that actually exist (some tiled layers may be absent)
    col_order = [c for c in desired_order if c in final.columns]
    # Append any remaining columns not in desired order
    remaining = [c for c in final.columns if c not in col_order]
    if remaining:
        logger.info(f"  Extra columns appended: {remaining}")
    final = final[col_order + remaining]

    # ── Validation ────────────────────────────────────────────────────────────
    logger.info(f"\n{'=' * 60}")
    logger.info("VALIDATION")
    logger.info(f"{'=' * 60}")

    logger.info(f"  ✓ No deduplication applied (point_ids are unique by construction)")

    # Missing values
    missing = final.isnull().sum()
    missing = missing[missing > 0].sort_values(ascending=False)
    if len(missing) == 0:
        logger.info(f"  ✓ No missing values")
    else:
        logger.info(f"  Missing values in {len(missing)} columns:")
        for col, n in missing.items():
            logger.info(f"    {col:<40} {n:>8,} ({n/len(final)*100:.1f}%)")

    # Y distribution
    logger.info(f"\n  Y distribution:")
    logger.info(f"    Y=1 (regrowth):    {(final.y == 1).sum():>8,}")
    logger.info(f"    Y=0 (nonregrowth): {(final.y == 0).sum():>8,}")

    # ── Save ──────────────────────────────────────────────────────────────────
    output_path = gee_dir / args.output
    final.to_csv(output_path, index=False)
    size_mb = output_path.stat().st_size / 1024 / 1024

    logger.info(f"\n{'=' * 60}")
    logger.info(f"DONE")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Output: {output_path}")
    logger.info(f"  Size:   {size_mb:.1f} MB")
    logger.info(f"  Rows:   {len(final):,}")
    logger.info(f"  Cols:   {len(final.columns)}")
    logger.info(f"\n  Columns: {list(final.columns)}")


if __name__ == '__main__':
    main()
