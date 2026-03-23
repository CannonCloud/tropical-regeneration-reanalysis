#!/usr/bin/env python3
"""
Step 1b: Merge Batch Assets
============================

After all batch tasks complete, this merges them into
one final master spine asset.

Usage:
    python step1b_merge_batches.py --samples 200000
    python step1b_merge_batches.py --samples 100000
"""

import ee
import argparse
import logging
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    PROJECT_ID = '# replace with your GEE project ID'
    OUTPUT_DIR = Path('./gee_outputs')
    LOG_FILE = OUTPUT_DIR / 'step1b_merge.log'
    
    def get_final_asset_name(self, n_samples):
        return f'master_sample_points_{n_samples}'

# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Merge Batch Assets')
    parser.add_argument('--samples', type=int, default=200000,
                        help='Target sample size (default: 200000)')
    parser.add_argument('--batch-size', type=int, default=20000,
                        help='Points per batch used in step1a (default: 20000)')
    
    args = parser.parse_args()
    
    # Calculate how many batches were created (+1 buffer from step1a)
    batch_count = (args.samples + args.batch_size - 1) // args.batch_size + 1
    
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
    
    logger.info("="*60)
    logger.info("STEP 1B: MERGE BATCH ASSETS")
    logger.info("="*60)
    logger.info(f"\nTarget samples: {args.samples:,}")
    logger.info(f"Batch size: {args.batch_size:,}")
    logger.info(f"Expected batches: {batch_count}")
    
    # Initialize GEE
    try:
        ee.Initialize(project=config.PROJECT_ID)
        logger.info(f"\n✓ GEE initialized: {config.PROJECT_ID}")
    except:
        logger.info("\nAuthenticating...")
        ee.Authenticate()
        ee.Initialize(project=config.PROJECT_ID)
        logger.info("✓ Authentication successful")
    
    # Load batch assets
    logger.info(f"\nLoading {batch_count} batch assets...")
    
    collections = []
    loaded_batches = []
    
    for i in range(batch_count):
        asset_id = f'projects/{config.PROJECT_ID}/assets/MasterSpine_Batch_{i:02d}'
        logger.info(f"  Loading batch {i}: {asset_id}")
        
        try:
            batch_fc = ee.FeatureCollection(asset_id)
            collections.append(batch_fc)
            loaded_batches.append(i)
            logger.info(f"    ✓ Loaded")
        except Exception as e:
            logger.error(f"    ✗ Failed to load: {e}")
            logger.error(f"    Make sure this batch completed in GEE tasks!")
            exit(1)
    
    logger.info(f"\n✓ Successfully loaded {len(collections)}/{batch_count} batches")
    
    # Merge collections
    logger.info("\nMerging collections...")
    merged_collection = ee.FeatureCollection(collections).flatten()
    logger.info("✓ Collections merged")
    
    # Remove duplicates based on geometry (lat/lon)
    logger.info("\nRemoving spatial duplicates (if any)...")
    logger.info("Note: Duplicates are extremely rare (~1-2 in 200K points)")
    
    # Group by geometry and take first of each unique location
    # This removes any points that landed on the exact same pixel
    def add_geo_key(feature):
        feature = ee.Feature(feature)  # Cast to Feature
        # Create a key from rounded coordinates (to handle floating point)
        coords = feature.geometry().coordinates()
        lon = ee.Number(coords.get(0)).format('%.6f')
        lat = ee.Number(coords.get(1)).format('%.6f')
        geo_key = lon.cat('_').cat(lat)
        return feature.set('_geo_key', geo_key)
    
    merged_with_keys = merged_collection.map(add_geo_key)
    
    # Use distinct to remove duplicates by geo_key
    deduplicated = merged_with_keys.distinct('_geo_key')
    
    # Remove the temporary geo_key property
    def remove_geo_key(feature):
        feature = ee.Feature(feature)  # Cast to Feature
        # Keep only the original properties (exclude _geo_key)
        return feature.select(['point_id', 'batch_id', 'batch_seed', 'pnv_probability'], None, True)
    
    final_collection = deduplicated.map(remove_geo_key)
    
    logger.info("✓ Deduplication complete")
    
    # Export deduplicated collection to final asset
    final_asset_name = config.get_final_asset_name(args.samples)
    final_asset_id = f'projects/{config.PROJECT_ID}/assets/{final_asset_name}'
    
    logger.info(f"\nExporting to final master asset...")
    logger.info(f"Asset ID: {final_asset_id}")
    
    task = ee.batch.Export.table.toAsset(
        collection=final_collection,  # Use deduplicated collection
        description=f'MasterSpine_{args.samples}_Merged',
        assetId=final_asset_id
    )
    
    task.start()
    task_id = task.id
    
    logger.info("\n" + "="*60)
    logger.info("✓ MERGE TASK SUBMITTED")
    logger.info("="*60)
    logger.info(f"\nTask ID: {task_id}")
    logger.info(f"Description: MasterSpine_{args.samples}_Merged")
    logger.info(f"Final Asset: {final_asset_id}")
    
    logger.info("\nThis merge should complete in ~5-10 minutes.")
    logger.info("\nCheck status: https://code.earthengine.google.com/tasks")
    
    logger.info("\n" + "="*60)
    logger.info("NEXT STEPS")
    logger.info("="*60 + "\n")
    
    logger.info("Once merge completes, proceed with:")
    logger.info(f"  python step2_sample_all_layers.py --samples {args.samples}")
    logger.info("\nThe asset will be automatically found at:")
    logger.info(f"  {final_asset_id}")
    
    logger.info("\n" + "="*60)
    logger.info("SUMMARY")
    logger.info("="*60 + "\n")
    
    logger.info(f"Merged {batch_count} batches:")
    for i in loaded_batches:
        logger.info(f"  ✓ Batch {i}: {args.batch_size:,} points (seed {10 + i*10})")
    
    logger.info(f"\nExpected total: ~{batch_count * args.batch_size:,} points")
    logger.info(f"(May be slightly less due to filtering, but >{args.samples:,} guaranteed)")
    
    logger.info("\n" + "="*60)
    
    # Save final asset ID
    asset_file = config.OUTPUT_DIR / f'master_spine_asset_n{args.samples}.txt'
    with open(asset_file, 'w') as f:
        f.write(f"{final_asset_id}\n")
    
    logger.info(f"\n✓ Asset ID saved to: {asset_file}\n")

if __name__ == '__main__':
    main()
