# GEE Tropical Forest Sampling - Complete Workflow Guide

## Overview

This workflow samples 200K+ points from tropical forests and extracts ~53 environmental variables for correlation analysis. It uses a **batch-based approach** to avoid GEE memory limits.

## The Complete Pipeline

```
Step 1a: Sample 220K points in 11 batches (20K each)
         ↓
Step 1b: Merge batches + deduplicate → Final asset
         ↓
Step 2: Sample 10 layer groups from asset
         ↓
Step 3: Calculate 6 distance variables in batches
         ↓
Step 4: Merge everything into final dataset
```

## Setup (One-Time)

```bash
# Create environment
cd ~/Dropbox/Regrowth/pnv_sampling
python3 -m venv venv
source venv/bin/activate
pip install earthengine-api pandas tabulate

# Authenticate with GEE (one time only)
earthengine authenticate
```

## Step 1a: Sample Point Locations in Batches

**What it does:** Creates 11 separate GEE assets, each with 20K points (220K total). Uses different random seeds to avoid duplicates.

**Why batches?** Sampling 200K points at once causes "Out of Memory" errors. 11× 20K = no problem!

```bash
python step1a_submit_batches.py --samples 200000
```

**Parameters:**
- `--samples 200000` - Target sample size
- `--batch-size 20000` - Points per batch (default: 20000)

**What happens:**
- Submits 11 parallel tasks to GEE
- Seeds: 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110
- Each task takes ~10-20 minutes
- They run in parallel on Google's servers

**Monitor progress:**
- https://code.earthengine.google.com/tasks
- Wait until ALL 11 tasks show "COMPLETED"

**Time:** ~20-30 minutes total (parallel execution)

## Step 1b: Merge Batches into Master Asset

**After all 11 batches complete:**

```bash
python step1b_merge_batches.py --samples 200000
```

**What it does:**
- Loads all 11 batch assets
- Removes spatial duplicates (if any)
- Merges into final asset: `master_sample_points_200000`

**Time:** ~5-10 minutes

**Output asset:**
```
projects/# replace with your GEE project ID/assets/master_sample_points_200000
```

This is your permanent master spine - use it for all subsequent steps!

## Step 2: Sample All Layer Groups

**What it does:** Extracts 12 environmental layer groups at each point location.

```bash
python step2_sample_layers.py --samples 200000
```

**Layer Groups (12 total):**
1. **Climate** (6 bands): WorldClim bioclimatic variables
2. **Topography** (3 bands): Elevation, slope, forest density 2018
3. **Land Cover** (3 bands): ESA CCI classification, cropland weights
4. **NPP & Fire** (2 bands): MODIS productivity and fire frequency
5. **Soil Chemical** (3 bands): Organic carbon, pH, biome ID
6. **Soil Texture** (2 bands): Sand and clay content
7. **Soil Physical** (2 bands): Bulk density, water content at field capacity
8. **Socioeconomic** (5 bands): Protected areas, roads, water, GDP, population
9. **High-Res** (4 bands): Nightlights, travel time, human modification, WorldPop
10. **Pasture** (1 band): Global Pasture Watch classifications
11. **Opportunity Cost** (1 band): Agricultural value (SPAM)
12. **Plantation Raw** (1 band): Xiao et al. planted forests

**Note:** Soil was split into 3 separate tasks (chemical, texture, physical) to avoid "computed value too large" errors.

**What happens:**
- Submits 12 export tasks to GEE (one per layer group)
- Each task exports a CSV to Google Drive
- Tasks run in parallel

**Time:** ~4-6 hours total for all 12 layers

**After completion:**
1. Download all 12 CSVs from Google Drive
2. Place in `gee_outputs/`
3. Files will be named:
   - `PNV_Layer_climate_N200000.csv`
   - `PNV_Layer_topography_N200000.csv`
   - `PNV_Layer_soil_chemical_N200000.csv`
   - `PNV_Layer_soil_texture_N200000.csv`
   - `PNV_Layer_soil_physical_N200000.csv`
   - ... etc for all 12 layers

## Step 3: Calculate Distance Variables (Grid-Based)

**What it does:** Computes 11 distance/density metrics using a 6°×6° grid approach. Much faster than per-point batches!

**IMPORTANT - Updated Variables:** Several bugs were fixed (unmask order, ocean handling) and 5 new variables added. See section below for which folders to delete/regenerate.

```bash
# Option 1: Run all 11 variables in parallel (RECOMMENDED)
python step3_grid_distances.py --samples 200000 --distance all

# Option 2: Run individually (if keeping some existing correct data)
python step3_grid_distances.py --samples 200000 --distance cropland_density &
python step3_grid_distances.py --samples 200000 --distance forest &
python step3_grid_distances.py --samples 200000 --distance urban &
python step3_grid_distances.py --samples 200000 --distance settlement &
python step3_grid_distances.py --samples 200000 --distance agriculture &
python step3_grid_distances.py --samples 200000 --distance plantation &
python step3_grid_distances.py --samples 200000 --distance nightlight_density &
python step3_grid_distances.py --samples 200000 --distance population_density &
python step3_grid_distances.py --samples 200000 --distance cropland_coverage &
python step3_grid_distances.py --samples 200000 --distance cropland &
python step3_grid_distances.py --samples 200000 --distance cultivated_grass &
wait
```

**Distance/Density Variables (13 columns from 11 calculation types):**
1. **bio_cropland_density** - Cropland density within 5km buffer (mean %)
2. **bio_dist_forest_2018** - Distance to nearest forest edge (meters)
3. **bio_forest_density_1km2** - Local forest density in 1km² area (0-100%)
4. **paper_dist_urban** - Distance to nearest urban area (meters)
5. **highres_dist_settlement** - Distance to nearest settlement (meters)
6. **highres_dist_agriculture** - Distance to agriculture (cropland OR cultivated grass)
7. **bio_dist_plantation** - Distance to planted forest (meters)
8. **bio_plantation_density_1km** - Local plantation density in 1km radius (0-100%)
9. **highres_nightlight_density_2km** - Smoothed nightlight intensity (2014-2016 median, 2km radius)
10. **highres_worldpop_density_2km** - Population density (2km radius mean)
11. **highres_cropland_coverage_1km** - Cropland coverage % (1km radius, WorldCover 10m)
12. **highres_dist_cropland** - Distance to cropland only (meters)
13. **highres_dist_cultivated_grass** - Distance to cultivated grassland only (meters)

**Grid Strategy:**
- Divides tropics (-25° to +25°) into 6°×6° tiles (~660km × 660km)
- Grid includes +1 to latitude steps to ensure full 23°N to 25°N coverage
- Adds 50km buffer around each tile (prevents edge effects)
- Only processes tiles with data points (skips ocean/desert)
- Result: ~540 tiles total (60 longitude × 9 latitude steps)
- Per variable: ~90 tiles with data (rest are ocean/desert)

**What happens:**
- For 200K points:
  - ~90 tiles per variable type (only tiles with data)
  - 11 variable types
  - = ~990 total tasks submitted
- Tasks run on GEE servers (you can close your laptop!)
- Each tile processes in 3-10 minutes

**Time:** ~2-3 hours total (parallel execution)

**Bug Fixes Applied (Feb 2026):**
The following bugs were identified and fixed:
- **Unmask order bug:** 5 variables had `.unmask()` before comparison instead of after, causing NULL handling issues
- **Ocean/water bug:** Distance/density calculations failed across water (coastal points got NULL values)
- **Plantation inversion bug:** Natural forest was being treated as plantation (logic was backwards)
- **Density ignore-NULL bug:** Forest density was ignoring ocean pixels instead of counting them as 0%

**Which folders to regenerate:**
If you ran Step 3 before Feb 2026, DELETE and regenerate these folders:
- `GEE_Grid_plantation/` - Had inverted logic bug
- `GEE_Grid_forest/` - Had ocean-ignore bug in density
- `GEE_Grid_urban/` - Had unmask order bug
- `GEE_Grid_settlement/` - Had unmask order bug
- `GEE_Grid_agriculture/` - Had unmask order bug

**New folders (didn't exist before):**
- `GEE_Grid_nightlight_density/` - NEW
- `GEE_Grid_population_density/` - NEW
- `GEE_Grid_cropland_coverage/` - NEW
- `GEE_Grid_cropland/` - NEW
- `GEE_Grid_cultivated_grass/` - NEW

**Unchanged (safe to keep):**
- `GEE_Grid_cropland_density/` - No bugs, can keep existing data

**After completion:**
1. Download all grid CSVs from Google Drive (11 folders)
2. Organize into folders:
   ```
   gee_outputs/GEE_Grid_cropland_density/
   gee_outputs/GEE_Grid_forest/
   gee_outputs/GEE_Grid_urban/
   gee_outputs/GEE_Grid_settlement/
   gee_outputs/GEE_Grid_agriculture/
   gee_outputs/GEE_Grid_plantation/
   gee_outputs/GEE_Grid_nightlight_density/
   gee_outputs/GEE_Grid_population_density/
   gee_outputs/GEE_Grid_cropland_coverage/
   gee_outputs/GEE_Grid_cropland/
   gee_outputs/GEE_Grid_cultivated_grass/
   ```
3. Each folder contains ~90 Grid_*.csv files (tiles with data)

**Pro tip:** Use `--distance all` to submit all 11 types at once. Total submission time: ~10-15 minutes.

## Step 4: Merge Everything

**What it does:** Combines all layer CSVs + grid distance CSVs into one final dataset.

```bash
python step4_merge_all.py
```

**What it does:**
- Uses first layer (climate) as base and extracts coordinates from `.geo` column
- Merges all 12 layer CSVs (keeps point_id, lon, lat)
- Merges all 11 grid distance folders (~90 CSVs each with data)
- Automatically excludes GEE metadata columns (batch_id, batch_seed, pnv_probability)
- Removes duplicate point_ids (if any)
- Reports missing data statistics

**Output:**
```
gee_outputs/pnv_final_dataset.csv
```

Your complete dataset: ~219,482 rows × ~59 columns, ready for analysis!

**Time:** ~30-60 seconds

**What you'll see:**
- Point count: 219,482 (not exactly 200K due to ocean/invalid area filtering in Step 1)
- Missing data: Should be minimal (<1%) with bug fixes applied
- Columns: point_id, longitude, latitude, + all environmental variables

## Complete Timeline

| Day | Tasks | Time |
|-----|-------|------|
| **Day 1** | Step 1a (submit batches) | 30 min |
| | Step 1b (merge) | 10 min |
| | Step 2 (layer sampling) | 4-6 hours |
| | Step 3 (grid distance calculations) | 2-3 hours |
| | Step 4 (merge) | 1 min |
| **Total** | | ~7-10 hours (mostly automated) |

**Note:** Step 3 now uses 6°×6° grid (faster than old 3°×3° approach)!

## File Organization

```
pnv_sampling/
├── venv/                           # Python environment
├── gee_outputs/                    # All outputs
│   ├── PNV_Layer_climate_N200000.csv           # Step 2 outputs (12 files)
│   ├── PNV_Layer_topography_N200000.csv
│   ├── PNV_Layer_soil_chemical_N200000.csv
│   ├── PNV_Layer_soil_texture_N200000.csv
│   ├── PNV_Layer_soil_physical_N200000.csv
│   ├── ...
│   ├── GEE_Grid_cropland_density/              # Step 3 outputs (11 folders)
│   │   ├── Grid_cropland_density_Lon-84_Lat-7_N200000.csv
│   │   ├── Grid_cropland_density_Lon120_Lat15_N200000.csv
│   │   └── ... (~90 grid tiles with data)
│   ├── GEE_Grid_forest/
│   ├── GEE_Grid_urban/
│   ├── GEE_Grid_settlement/
│   ├── GEE_Grid_agriculture/
│   ├── GEE_Grid_plantation/
│   ├── GEE_Grid_nightlight_density/
│   ├── GEE_Grid_population_density/
│   ├── GEE_Grid_cropland_coverage/
│   ├── GEE_Grid_cropland/
│   ├── GEE_Grid_cultivated_grass/
│   ├── pnv_final_dataset.csv      # Step 4 output (FINAL!)
│   └── *.log                       # Execution logs
└── *.py                            # Scripts
```

## Key Advantages of This Approach

✅ **Memory-safe**: No OOM crashes (batches are small enough)  
✅ **Resumable**: Can restart from any failed step  
✅ **Modular**: Each layer/distance is independent  
✅ **Server-side**: All heavy lifting on GEE, not your laptop  
✅ **Reproducible**: Fixed seeds (10, 20, 30...) for each batch  
✅ **Flexible**: Easy to add new layers later

## Troubleshooting

### "Asset not found" in Step 2/3
**Problem:** Step 1b hasn't completed yet  
**Solution:** Check https://code.earthengine.google.com/tasks - wait for merge to finish

### "Too many concurrent tasks" 
**Problem:** GEE rate limiting  
**Solution:** Wait 10 minutes, tasks will resume automatically

### "Computation timed out" in Step 3
**Problem:** Batch size too large  
**Solution:** Reduce `--batch-size` to 500:
```bash
python step3_calc_all_distances.py --samples 200000 --batch-size 500
```

### Missing CSVs after Step 2/3
**Problem:** Task failed silently  
**Solution:** Check GEE tasks page for failed tasks, rerun specific layer:
```bash
python step2_sample_layers_from_asset.py --layer climate --samples 200000
```

### Step 4 merge fails with "KeyError"
**Problem:** Missing a required CSV file  
**Solution:** Check which file is missing:
```bash
ls gee_outputs/layer_*.csv  # Should show 10 files
ls gee_outputs/distance_*/  # Should show 6 folders
```

## Quick Commands Reference

```bash
# SETUP (once)
python3 -m venv venv
source venv/bin/activate
pip install earthengine-api pandas
earthengine authenticate

# STEP 1: Point sampling (2 stages)
python step1a_submit_batches.py --samples 200000
# Wait for all 11 batches to complete
python step1b_merge_batches.py --samples 200000

# STEP 2: Layer sampling (12 layers including 3 soil splits)
python step2_sample_layers.py --samples 200000

# STEP 3: Grid distance calculations (11 variables, 6° grid)
# Option 1: All at once (RECOMMENDED)
python step3_grid_distances.py --samples 200000 --distance all

# Option 2: Run individually (parallel submission)
python step3_grid_distances.py --samples 200000 --distance cropland_density &
python step3_grid_distances.py --samples 200000 --distance forest &
python step3_grid_distances.py --samples 200000 --distance urban &
python step3_grid_distances.py --samples 200000 --distance settlement &
python step3_grid_distances.py --samples 200000 --distance agriculture &
python step3_grid_distances.py --samples 200000 --distance plantation &
python step3_grid_distances.py --samples 200000 --distance nightlight_density &
python step3_grid_distances.py --samples 200000 --distance population_density &
python step3_grid_distances.py --samples 200000 --distance cropland_coverage &
python step3_grid_distances.py --samples 200000 --distance cropland &
python step3_grid_distances.py --samples 200000 --distance cultivated_grass &
wait

# STEP 4: Final merge
python step4_merge_all.py

# DONE! Your data: gee_outputs/pnv_final_dataset.csv
```

## Variables in Final Dataset (~59 total)

**Location (3):** point_id, longitude, latitude  
**Climate (6):** temp_mean, temp_seasonality, min_temp_coldest, precip_annual, precip_driest, precip_seasonality  
**Topography (3):** elevation, slope, forest_density_2018  
**Land Cover (3):** landcover_class, cropland_weight, lc_raw_value  
**Productivity (2):** npp_mean, fire_freq  
**Soil Chemical (3):** soil_organic_carbon, soil_ph, biome_id  
**Soil Texture (2):** soil_sand, soil_clay  
**Soil Physical (2):** soil_bulk_density_fine, soil_water_33kpa  
**Socioeconomic (5):** ghs_pop, gdp_2015, roads_density, dist_water, protected_binary  
**High-Res (4):** lights_2014_2016, travel_time, ghm, worldpop_density  
**Land Use (3):** gpw_raw_class, opp_cost_usd_ha, plantation_raw_b1  
**Distances (7):** dist_forest_2018, dist_urban, dist_settlement, dist_agriculture, dist_plantation, dist_cropland, dist_cultivated_grass  
**Densities (6):** cropland_density_5km, forest_density_1km2, plantation_density_1km, nightlight_density_2km, worldpop_density_2km, cropland_coverage_1km

**Note:** Soil silt and available water capacity (AWC) were removed due to OpenLandMap asset availability issues.

## Notes

- **Scale:** Points sampled at 300m resolution (faster, still robust for correlations)
- **Deduplication:** Spatial duplicates removed in Step 1b (typically <10 duplicates)
- **Sample size:** ~219,482 points (target was 200K, some filtered due to ocean/invalid areas)
- **Seeds:** Batches use seeds 10, 20, 30...110 for full reproducibility
- **Asset permanence:** Master spine asset never needs to be recreated - reuse forever!
- **Grid efficiency:** 6°×6° grid (updated from 3°×3°) with 50km buffer processes ~540 tiles instead of 200K+ individual points
- **Grid coverage:** Latitude calculation includes +1 to ensure full 23°N to 25°N coverage (no missing strip!)
- **Soil split:** Soil variables split into 3 tasks (chemical, texture, physical) to avoid GEE memory limits
- **Missing variables:** Silt and AWC removed from soil suite due to OpenLandMap asset inconsistencies
- **Bug fixes (Feb 2026):** Fixed unmask order issues, ocean handling, and plantation inversion logic
- **New variables (Feb 2026):** Added nightlight density, population density, cropland coverage, and separated cropland/cultivated grass distances

## Advanced: Custom Sample Sizes

Want 100K instead of 200K?

```bash
# Step 1a: Fewer batches
python step1a_submit_batches.py --samples 100000 --batch-size 20000
# Creates 6 batches (5 needed + 1 buffer)

# Step 1b: Merge with matching args
python step1b_merge_batches.py --samples 100000 --batch-size 20000

# Steps 2-4: Use --samples 100000
python step2_sample_layers.py --samples 100000

# Step 3: All distance variables
python step3_grid_distances.py --samples 100000 --distance all
# Or individually with matching --samples flag

python step4_merge_all.py
```

The scripts automatically calculate the right number of batches/tiles!
