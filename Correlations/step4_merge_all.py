#!/usr/bin/env python3
"""
Step 4: Smart Merge All Layers Into Final Dataset
==================================================

This script AUTOMATICALLY discovers and merges:
1. All PNV_Layer_*.csv files (from Step 2)
2. All GEE_Grid_* directories (from Step 3)

No hardcoded file lists - adapts to whatever you have!

Usage:
    python step4_merge_all_smart.py
    python step4_merge_all_smart.py --output final_dataset.csv
    python step4_merge_all_smart.py --gee-dir ./gee_outputs
"""

import argparse
import logging
from pathlib import Path
import pandas as pd
from datetime import datetime
import json

# =============================================================================
# AUTO-DISCOVERY FUNCTIONS
# =============================================================================

def discover_layer_files(output_dir: Path, logger):
    """Find all PNV_Layer_*.csv files automatically"""
    pattern = 'PNV_Layer_*.csv'
    layer_files = sorted(output_dir.glob(pattern))
    
    logger.info(f"Discovered {len(layer_files)} layer files:")
    for f in layer_files:
        # Extract layer name from filename: PNV_Layer_XXXXX_N200000.csv -> XXXXX
        name = f.stem.replace('PNV_Layer_', '').rsplit('_', 1)[0]
        logger.info(f"  - {name}: {f.name}")
    
    return layer_files

def discover_distance_dirs(output_dir: Path, logger):
    """Find all GEE_Grid_* directories automatically"""
    pattern = 'GEE_Grid_*'
    distance_dirs = sorted([d for d in output_dir.glob(pattern) if d.is_dir()])
    
    logger.info(f"Discovered {len(distance_dirs)} distance directories:")
    for d in distance_dirs:
        # Count CSV files in each directory
        csv_count = len(list(d.glob('Grid_*.csv')))
        name = d.name.replace('GEE_Grid_', '')
        logger.info(f"  - {name}: {csv_count} grid files")
    
    return distance_dirs

def extract_coords_from_geo(geo_str):
    """Extract longitude, latitude from GeoJSON string"""
    try:
        geo = json.loads(geo_str)
        coords = geo['coordinates']
        return coords[0], coords[1]  # longitude, latitude
    except:
        return None, None

# =============================================================================
# MERGE LOGIC
# =============================================================================

def merge_all_data(output_dir: Path, logger, output_file: Path):
    """Merge all CSVs into final dataset using auto-discovery"""
    
    # =========================================================================
    # STEP 1: Auto-discover all files
    # =========================================================================
    logger.info("="*60)
    logger.info("STEP 1: Auto-Discover Available Data")
    logger.info("="*60)
    
    layer_files = discover_layer_files(output_dir, logger)
    distance_dirs = discover_distance_dirs(output_dir, logger)
    
    if not layer_files:
        logger.error("✗ No PNV_Layer_*.csv files found!")
        logger.error(f"  Searched in: {output_dir}")
        exit(1)
    
    # =========================================================================
    # STEP 2: Load base point locations from first layer
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("STEP 2: Load Base Point Locations")
    logger.info("="*60)
    
    base_file = layer_files[0]
    base_name = base_file.stem.replace('PNV_Layer_', '').rsplit('_', 1)[0]
    logger.info(f"Using {base_name} as base layer: {base_file.name}")
    
    base_df = pd.read_csv(base_file)
    
    # Extract lat/long from .geo column if present
    if '.geo' in base_df.columns:
        logger.info("Extracting coordinates from .geo column...")
        base_df[['longitude', 'latitude']] = base_df['.geo'].apply(
            lambda x: pd.Series(extract_coords_from_geo(x))
        )
        logger.info(f"✓ Extracted lat/long for {len(base_df):,} points")
    
    # Keep only point_id, longitude, latitude, pnv_probability as base
    base_cols = ['point_id', 'longitude', 'latitude']
    if 'pnv_probability' in base_df.columns:
        base_cols.append('pnv_probability')
        logger.info("  Including pnv_probability from base layer")
    
    df = base_df[base_cols].copy()
    logger.info(f"✓ Base dataset: {len(df):,} points with coordinates")
    
    # Now merge the data columns from the base layer (excluding metadata)
    exclude_cols = ['.geo', 'system:index', 'longitude', 'latitude', 
                   'batch_id', 'batch_seed', 'pnv_probability']  # Exclude pnv_probability here - already added above
    data_cols = [c for c in base_df.columns if c not in exclude_cols and c != 'point_id']
    
    if data_cols:
        logger.info(f"  Adding {len(data_cols)} columns from {base_name}")
        for col in data_cols:
            df[col] = base_df[col]
    
    # =========================================================================
    # STEP 3: Merge remaining layer files
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("STEP 3: Merge Remaining Layer Files")
    logger.info("="*60)
    
    for layer_file in layer_files[1:]:  # Skip first (already loaded as base)
        layer_name = layer_file.stem.replace('PNV_Layer_', '').rsplit('_', 1)[0]
        
        logger.info(f"\nMerging {layer_name}...")
        layer_df = pd.read_csv(layer_file)
        
        # Extract only data columns (exclude geometry and GEE metadata)
        exclude_cols = ['.geo', 'system:index', 'longitude', 'latitude', 
                       'batch_id', 'batch_seed', 'pnv_probability']
        data_cols = [c for c in layer_df.columns if c not in exclude_cols and c != 'point_id']
        
        logger.info(f"  Rows: {len(layer_df):,}")
        logger.info(f"  New columns: {len(data_cols)} - {', '.join(data_cols[:5])}")
        if len(data_cols) > 5:
            logger.info(f"               ... and {len(data_cols) - 5} more")
        
        # Merge on point_id
        merge_cols = ['point_id'] + data_cols
        df = df.merge(layer_df[merge_cols], on='point_id', how='left')
        
        logger.info(f"✓ Merged {layer_name}")
        logger.info(f"  Total columns now: {len(df.columns)}")
    
    # =========================================================================
    # STEP 4: Merge distance/density calculations
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("STEP 4: Merge Distance/Density Calculations")
    logger.info("="*60)
    
    for dist_dir in distance_dirs:
        dist_name = dist_dir.name.replace('GEE_Grid_', '')
        
        # Find all grid files
        batch_files = sorted(dist_dir.glob('Grid_*.csv'))
        
        if not batch_files:
            logger.warning(f"⚠ No Grid_*.csv files found in {dist_dir.name}")
            continue
        
        logger.info(f"\nMerging {dist_name} ({len(batch_files)} grid files)...")
        
        # Concatenate all batches
        batch_dfs = []
        for batch_file in batch_files:
            batch_df = pd.read_csv(batch_file)
            batch_dfs.append(batch_df)
        
        dist_df = pd.concat(batch_dfs, ignore_index=True)
        logger.info(f"  Total rows: {len(dist_df):,}")
        
        # Extract data columns (exclude geometry and GEE metadata)
        exclude_cols = ['.geo', 'system:index', 'longitude', 'latitude',
                       'batch_id', 'batch_seed', 'pnv_probability']
        data_cols = [c for c in dist_df.columns if c not in exclude_cols and c != 'point_id']
        
        logger.info(f"  New columns: {len(data_cols)} - {', '.join(data_cols)}")
        
        # Merge on point_id
        merge_cols = ['point_id'] + data_cols
        df = df.merge(dist_df[merge_cols], on='point_id', how='left')
        
        logger.info(f"✓ Merged {dist_name}")
        logger.info(f"  Total columns now: {len(df.columns)}")
    
    # =========================================================================
    # STEP 5: Final validation
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("STEP 5: Final Validation")
    logger.info("="*60)
    
    logger.info(f"\nFinal dataset statistics:")
    logger.info(f"  Rows: {len(df):,}")
    logger.info(f"  Columns: {len(df.columns)}")
    
    # Show column names organized by type
    logger.info(f"\nColumn inventory:")
    
    # ID columns
    id_cols = ['point_id', 'longitude', 'latitude']
    logger.info(f"  ID columns ({len(id_cols)}): {', '.join(id_cols)}")
    
    # Distance columns
    dist_cols = [c for c in df.columns if c.startswith(('highres_dist_', 'bio_dist_'))]
    if dist_cols:
        logger.info(f"  Distance columns ({len(dist_cols)}): {', '.join(dist_cols)}")
    
    # Density columns
    density_cols = [c for c in df.columns if 'density' in c.lower()]
    if density_cols:
        logger.info(f"  Density columns ({len(density_cols)}): {', '.join(density_cols)}")
    
    # Other columns
    other_cols = [c for c in df.columns if c not in id_cols + dist_cols + density_cols]
    if other_cols:
        logger.info(f"  Other columns ({len(other_cols)}): {', '.join(other_cols[:10])}")
        if len(other_cols) > 10:
            logger.info(f"                     ... and {len(other_cols) - 10} more")
    
    # Check for missing data
    logger.info(f"\nMissing value analysis:")
    missing = df.isnull().sum()
    missing = missing[missing > 0].sort_values(ascending=False)
    
    if len(missing) == 0:
        logger.info("  ✓ No missing values!")
    else:
        logger.info(f"  Found missing values in {len(missing)} columns:")
        for col, count in list(missing.items())[:10]:  # Show top 10
            pct = (count / len(df)) * 100
            logger.info(f"    {col}: {count:,} ({pct:.1f}%)")
        if len(missing) > 10:
            logger.info(f"    ... and {len(missing) - 10} more columns with missing data")
    
    # Check for duplicates
    dup_count = df.duplicated(subset=['point_id']).sum()
    if dup_count > 0:
        logger.warning(f"⚠ Found {dup_count} duplicate point_ids")
        logger.info("  Removing duplicates...")
        df = df.drop_duplicates(subset=['point_id'], keep='first')
        logger.info(f"  ✓ Removed duplicates, {len(df):,} rows remaining")
    else:
        logger.info("  ✓ No duplicates found")
    
    # =========================================================================
    # STEP 6: Save
    # =========================================================================
    logger.info("\n" + "="*60)
    logger.info("STEP 6: Save Final Dataset")
    logger.info("="*60)
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)
    
    file_size = output_file.stat().st_size / (1024 * 1024)  # MB
    logger.info(f"\n✓ Final dataset saved to: {output_file}")
    logger.info(f"  File size: {file_size:.1f} MB")
    logger.info(f"  Rows: {len(df):,}")
    logger.info(f"  Columns: {len(df.columns)}")
    
    logger.info("\n" + "="*60)
    logger.info("MERGE COMPLETE!")
    logger.info("="*60)
    logger.info(f"\nYour final dataset is ready: {output_file}")
    logger.info("You can now use this for analysis.")
    logger.info("\n" + "="*60)
    
    return df

# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Smart Merge - Auto-discovers and merges all layer and distance files'
    )
    parser.add_argument('--output', type=str, default='pnv_final_dataset.csv',
                        help='Output filename (default: pnv_final_dataset.csv)')
    parser.add_argument('--gee-dir', type=str, default='./gee_outputs',
                        help='Directory containing GEE outputs (default: ./gee_outputs)')
    
    args = parser.parse_args()
    
    # Setup
    output_dir = Path(args.gee_dir)
    if not output_dir.exists():
        print(f"Error: Directory not found: {output_dir}")
        exit(1)
    
    log_file = output_dir / 'step4_merge_smart.log'
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)
    
    logger.info("="*60)
    logger.info("STEP 4: SMART MERGE ALL DATA")
    logger.info("="*60)
    logger.info(f"Input directory: {output_dir}")
    logger.info(f"Output file: {args.output}")
    logger.info("")
    
    output_file = output_dir / args.output
    
    # Run merge
    df = merge_all_data(output_dir, logger, output_file)

if __name__ == '__main__':
    main()
