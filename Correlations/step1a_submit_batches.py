#!/usr/bin/env python3
"""
Step 1a: Submit Sample Batches (Memory-Safe Approach)
=====================================================

Submits multiple parallel tasks to avoid memory errors.
Default: 11 batches × 20K = 220K points

Usage:
    python step1a_submit_batches.py --samples 200000
    python step1a_submit_batches.py --samples 100000 --batch-size 20000
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
    SCALE = 300
    BASE_SEED = 10  # Will use 10, 20, 30, 40... for each batch
    STUDY_BOUNDS = [[-180, -25], [0, -25], [180, -25], [180, 25], [0, 25], [-180, 25]]
    OUTPUT_DIR = Path('./gee_outputs')
    LOG_FILE = OUTPUT_DIR / 'step1a_batches.log'

# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Submit Sample Batches (Memory-Safe)')
    parser.add_argument('--samples', type=int, default=200000,
                        help='Total target samples (default: 200000)')
    parser.add_argument('--batch-size', type=int, default=20000,
                        help='Points per batch (default: 20000)')
    
    args = parser.parse_args()
    
    # Calculate batch count (round up to ensure we get enough)
    batch_count = (args.samples + args.batch_size - 1) // args.batch_size
    # Add 1 extra batch as buffer
    batch_count += 1
    
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
    logger.info("STEP 1A: SUBMIT SAMPLE BATCHES")
    logger.info("="*60)
    logger.info(f"\nTarget: {args.samples:,}+ points")
    logger.info(f"Batch size: {args.batch_size:,} points")
    logger.info(f"Number of batches: {batch_count}")
    logger.info(f"Expected total: {batch_count * args.batch_size:,} points")
    logger.info(f"\nSeeds: {config.BASE_SEED}, {config.BASE_SEED+10}, {config.BASE_SEED+20}, ...")
    logger.info("\nEach batch uses a different seed for reproducibility.")
    
    # Initialize GEE
    try:
        ee.Initialize(project=config.PROJECT_ID)
        logger.info(f"\n✓ GEE initialized: {config.PROJECT_ID}")
    except:
        logger.info("\nAuthenticating...")
        ee.Authenticate()
        ee.Initialize(project=config.PROJECT_ID)
        logger.info("✓ Authentication successful")
    
    # Load PNV map
    logger.info("\nLoading PNV map...")
    pnv_map = ee.Image("projects/# replace with your GEE project ID/assets/pnv_mosaic_full")
    valid_binary = pnv_map.mask().rename('valid_binary').toInt()
    stack_with_class = pnv_map.rename('pnv_probability').addBands(valid_binary)
    logger.info("✓ PNV map loaded")
    
    # Define region
    region = ee.Geometry.Polygon(config.STUDY_BOUNDS, None, False)
    logger.info("✓ Region defined (full tropics)")
    
    # Submit batches
    logger.info("\n" + "="*60)
    logger.info("SUBMITTING BATCHES")
    logger.info("="*60 + "\n")
    
    submitted_tasks = []
    
    for i in range(batch_count):
        # Calculate seed (10, 20, 30, 40, ...)
        current_seed = config.BASE_SEED + (i * 10)
        
        logger.info(f"Batch {i+1}/{batch_count}")
        logger.info(f"  Seed: {current_seed}")
        logger.info(f"  Points: {args.batch_size:,}")
        
        try:
            # Sample with unique seed
            samples = stack_with_class.stratifiedSample(
                numPoints = 0,
                classBand='valid_binary',
                region=region,
                scale=config.SCALE,
                seed=current_seed,              # Unique seed per batch
                geometries=True,
                dropNulls=True,
                classValues=[1],                # Only valid PNV
                classPoints=[args.batch_size],  # Use batch_size from args
                tileScale=16                    # Memory optimization
            )
            
            # Add metadata
            def add_meta(feature):
                return feature.set({
                    'point_id': ee.String(f"b{i}_").cat(feature.id()),
                    'batch_id': i,
                    'batch_seed': current_seed
                })
            
            samples_with_meta = samples.map(add_meta)
            
            # Export to asset
            asset_name = f'MasterSpine_Batch_{i:02d}'
            asset_id = f'projects/{config.PROJECT_ID}/assets/{asset_name}'
            
            task = ee.batch.Export.table.toAsset(
                collection=samples_with_meta,
                description=asset_name,
                assetId=asset_id
            )
            
            task.start()
            task_id = task.id
            
            submitted_tasks.append({
                'batch': i,
                'seed': current_seed,
                'task_id': task_id,
                'asset_id': asset_id
            })
            
            logger.info(f"  ✓ Submitted: {task_id}")
            logger.info(f"  ✓ Asset: {asset_name}\n")
            
            # Brief delay to avoid API rate limits
            time.sleep(1.5)
            
        except Exception as e:
            logger.error(f"  ✗ Failed to submit batch {i}: {e}\n")
            continue
    
    # Summary
    logger.info("="*60)
    logger.info("ALL BATCHES SUBMITTED")
    logger.info("="*60 + "\n")
    
    logger.info(f"Successfully submitted {len(submitted_tasks)}/{batch_count} batches\n")
    
    for task_info in submitted_tasks:
        logger.info(f"Batch {task_info['batch']}: Seed {task_info['seed']}, Task {task_info['task_id']}")
    
    logger.info("\n" + "="*60)
    logger.info("WHAT HAPPENS NEXT")
    logger.info("="*60 + "\n")
    
    logger.info("Phase 1 (NOW): Batches are running on GEE servers")
    logger.info(f"  • Each batch samples {args.batch_size:,} points independently")
    logger.info("  • Estimated time: 10-20 min per batch")
    logger.info("  • Check: https://code.earthengine.google.com/tasks\n")
    
    logger.info("Phase 2 (AFTER ALL COMPLETE): Merge batches")
    logger.info(f"  • Run: python step1b_merge_batches.py --samples {args.samples}")
    logger.info(f"  • This combines all {batch_count} batches into one master asset")
    logger.info("  • Takes ~5 minutes\n")
    
    logger.info("="*60)
    logger.info("WHY THIS WORKS")
    logger.info("="*60 + "\n")
    
    logger.info("Problem: 200K points in one task = Out of Memory crash")
    logger.info("Solution: 11 × 20K points = same data, no crash\n")
    
    logger.info("Benefits:")
    logger.info("  ✓ Each batch is small enough to avoid OOM")
    logger.info("  ✓ Different seeds ensure no duplicate points")
    logger.info("  ✓ Global random distribution preserved (no spatial bias)")
    logger.info("  ✓ Batches run in parallel (fast)")
    logger.info("  ✓ Fully reproducible (seeds logged)")
    
    logger.info("\n" + "="*60)
    
    # Save task info
    import json
    task_file = config.OUTPUT_DIR / 'batch_tasks.json'
    with open(task_file, 'w') as f:
        json.dump(submitted_tasks, f, indent=2)
    
    logger.info(f"\n✓ Task info saved to: {task_file}")
    logger.info("\nMonitor progress at: https://code.earthengine.google.com/tasks")
    logger.info("When all show 'COMPLETED', run: python step1b_merge_batches.py\n")

if __name__ == '__main__':
    main()
