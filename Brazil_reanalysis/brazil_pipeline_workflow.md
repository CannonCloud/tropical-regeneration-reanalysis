# Brazil Natural Regeneration — Sampling Pipeline

Replication and extension of Williams et al. (2024) *Nature* for Brazil, using a
Brazilian-specific training dataset to test sensitivity of predicted regeneration
potential to (1) prevalence assumptions and (2) climate envelope confounding.

---

## Overview

The pipeline produces a training dataset of ~1.2M georeferenced points across three
sample sets, each with a full suite of biophysical and socioeconomic covariates.

| Sample set | Y | N (approx) | Purpose |
|---|---|---|---|
| `regrowth` | 1 | 200,000 | Locations where natural regeneration occurred 2000–2016 (Fagan et al.) |
| `nonregrowth` | 0 | 800,000 | Brazil-wide eligible non-regrowth pixels |
| `annulus` | 0 | 200,000 | Non-regrowth pixels within 3km of any regrowth pixel (climate envelope analysis) |

**GEE project:** `# set your GEE ID here`  
**Working directory:** `# set your working directory here`  
**GEE outputs:** `./gee_outputs/`

---

## Prerequisites

### Python environment
```bash
cd ~/Dropbox/Regrowth/pnv_sampling
source venv/bin/activate
pip install earthengine-api pandas
earthengine authenticate
```

### GEE assets required (pre-existing)
| Asset | Description |
|---|---|
| `brazil_plantations_v2_0_hansen_validated` | Fagan et al. regrowth/plantation/open classification |
| `lc_300m_2000` | ESA CCI land cover, year 2000 |
| `brazil_worldclim_v2` | WorldClim v2.1 bioclim variables 1–19 (uploaded manually) |
| `UMD/hansen/global_forest_change_2021_v1_9` | Hansen Global Forest Change (public) |
| `RESOLVE/ECOREGIONS/2017` | RESOLVE biomes (public) |
| `distance2water_30arcsec` | Distance to water (uploaded asset) |
| `grip4_density` | GRIP4 road density (uploaded asset) |
| `Opp_Cost_USD_per_Ha` | Agricultural opportunity cost — SPAM (uploaded asset) |

---

## Step 0 — Export Binary Mask Assets

**Script:** `brazil_step0_export_masks.py`

Exports two binary raster assets to GEE that define the sampling domains for all
downstream steps. Pre-baking these as assets avoids recomputing the same expensive
multi-layer expressions for every sampling batch.

**Mask A — `brazil_regrowth_mask_v2`**  
Binary (1 = natural regrowth). Pixels where Fagan class == 1. Filtered to RESOLVE
tropical forest biomes (BIOME_NUM 1–3) and ≤25°S latitude.

**Mask B — `brazil_nonregrowth_mask_v6`**  
Binary (1 = eligible non-regrowth). Pixels satisfying all of:
- Not forest in 2000 (Hansen treecover2000 < 30%)
- Not water / bare / urban / sparse vegetation (ESA CCI 2000 classes 150, 190, 200, 210)
- Not any Fagan pixel (classes 1, 2, 3)
- Within RESOLVE tropical forest biomes and ≤25°S latitude

All inputs are reprojected to the Hansen grid (`crsTransform = [0.00025, 0, -180, 0, -0.00025, 80]`)
before combining to ensure pixel-perfect alignment.

```bash
python brazil_step0_export_masks.py           # both masks
python brazil_step0_export_masks.py --mask A  # regrowth only
python brazil_step0_export_masks.py --mask B  # non-regrowth only
python brazil_step0_export_masks.py --test    # Manaus test box (~minutes)
```

**Runtime:** 10–20 min per mask. Verify both in GEE code editor before proceeding.

---

## Step 0b — Export Near-Regrowth Annulus Mask

**Script:** `brazil_step0b_export_annuli.py`

Exports a third mask defining eligible non-regrowth pixels within 3km of any regrowth
pixel. Used to sample the `annulus` control group for the climate envelope analysis.

Uses `fastDistanceTransform()` on the regrowth raster — a pure raster operation that
computes per-pixel distance to the nearest regrowth pixel. This avoids the vector
buffering approach (buffering 100K points) which produces geometry with millions of
vertices that GEE cannot handle.

**Output asset:** `brazil_nonregrowth_annulus_0_3km_v2`

```bash
python brazil_step0b_export_annuli.py           # full Brazil, 3km threshold
python brazil_step0b_export_annuli.py --test    # Manaus test box
python brazil_step0b_export_annuli.py --dist 5000  # change distance threshold
```

**Runtime:** ~30 min.

---

## Step 1 — Sample Random Points

**Script:** `brazil_step1_sample_points.py`

Samples random points from each mask using `stratifiedSample()`. Points are exported
as GEE FeatureCollection assets with metadata columns: `point_id`, `sample_type`,
`batch_id`, `seed`.

### Task A — Regrowth points (Y=1)

Samples 100,000 points per batch × 2 batches = 200,000 total from `brazil_regrowth_mask_v2`.
Each batch uses a different seed. Scale: 30m.

```bash
python brazil_step1_sample_points.py --task A
```

**Output assets:** `brazil_regrowth_batch_00`, `brazil_regrowth_batch_01`

### Task B — Non-regrowth points (Y=0)

Samples 200,000 points per batch × 4 batches = 800,000 total from `brazil_nonregrowth_mask_v6`.
Scale: 30m.

```bash
python brazil_step1_sample_points.py --task B
```

**Output assets:** `brazil_nonregrowth_batch_02` through `brazil_nonregrowth_batch_05`

> **Note:** Batch numbering starts at 02 to avoid overwriting earlier test runs.
> Controlled via `--batch-start 2`.

### Task C — Annulus points (Y=0, near-regrowth)

Samples 200,000 points from `brazil_nonregrowth_annulus_0_3km_v2`. Uses:
- `scale=90` (vs 30m for Tasks A/B) — reduces pixel count 9× and avoids OOM. Distribution
  bias from coarser scale is irrelevant since points are used for climate control, not
  spatial representativeness.
- `region=brazil_geom.bounds()` — bounding box instead of complex LSIB polygon.
  The annulus mask was already clipped to Brazil in Step 0b so points outside Brazil
  cannot be sampled. The `.bounds()` call eliminates expensive point-in-polygon checks
  against the LSIB polygon's tens of thousands of vertices.
- `tileScale=4` — subdivides sampling grid to manage memory on the thin/complex annulus mask.

```bash
python brazil_step1_sample_points.py --task C
python brazil_step1_sample_points.py --task C --test  # Manaus test, 1K points
```

**Output asset:** `brazil_nonregrowth_annulus_batch_00`

---

## ⚠ Manual Step — Merge Point Batches in GEE

After Step 1 tasks complete, the individual batch assets must be merged into combined
assets before Step 2 can run. This is done in the **GEE Code Editor** (JavaScript).

### Merge regrowth batches (A)

```javascript
var batch00 = ee.FeatureCollection('projects/# replace with your GEE project ID/assets/brazil_regrowth_batch_00');
var batch01 = ee.FeatureCollection('projects/# replace with your GEE project ID/assets/brazil_regrowth_batch_01');

var merged = batch00.merge(batch01);
print('Total regrowth points:', merged.size());

Export.table.toAsset({
  collection:  merged,
  description: 'BrazilRegrowthPointsMerged',
  assetId:     'projects/# replace with your GEE project ID/assets/brazil_regrowth_points'
});
```

### Post-sampling filter — remove points overlapping regrowth mask

Before merging non-regrowth batches, run a 30m pixel-level check to remove any
non-regrowth points that land on actual regrowth pixels (can occur due to resolution
mismatch between 120m sampling scale and 30m mask). This produced ~3.7% overlap in
testing.

Run the following in GEE Code Editor for each asset to be cleaned, then export as
`brazil_nonregrowth_points_Cleaned` and `brazil_nonregrowth_annulus_batch_00_Cleaned`:

```javascript
var points       = ee.FeatureCollection('projects/# replace with your GEE project ID/assets/brazil_nonregrowth_points');
var regrowthMask = ee.Image('projects/# replace with your GEE project ID/assets/brazil_regrowth_mask_v2').rename('regrowth_flag');

var checked = points.map(function(feature) {
  var sample = regrowthMask.reduceRegion({
    reducer:   ee.Reducer.max(),
    geometry:  feature.geometry(),
    scale:     30,
    maxPixels: 1e9
  });
  var isRegrowth = ee.Number(sample.get('regrowth_flag'));
  var flag = ee.Algorithms.If(isRegrowth, isRegrowth, 0);
  return feature.set('overlap_regrowth', flag);
});

var cleaned = checked.filter(ee.Filter.eq('overlap_regrowth', 0));

Export.table.toAsset({
  collection:  cleaned,
  description: 'Clean_Brazil_NonRegrowth_Points',
  assetId:     'projects/# replace with your GEE project ID/assets/brazil_nonregrowth_points_Cleaned'
});
```

**Cleaned assets required for Step 2:**
- `brazil_regrowth_points` (merged, no cleaning needed)
- `brazil_nonregrowth_points_Cleaned`
- `brazil_nonregrowth_annulus_batch_00_Cleaned`

---

## Step 2 — Sample Environmental Layers

**Script:** `brazil_step2_sample_layers.py`

Extracts raster values at each point location for all biophysical and socioeconomic
covariates. One GEE export task is submitted per layer group per sample set. Outputs
are CSVs exported to Google Drive.

### Layer groups

| Layer group | Columns | Source |
|---|---|---|
| `climate_part1` | bio01–bio04 | WorldClim v2.1 (uploaded asset) |
| `climate_part2` | bio05–bio08 | WorldClim v2.1 |
| `climate_part3` | bio09–bio12 | WorldClim v2.1 |
| `climate_part4` | bio13–bio16 | WorldClim v2.1 |
| `climate_part5` | bio17–bio19 | WorldClim v2.1 |
| `topography` | bio_elevation, bio_slope | SRTM 30m |
| `landcover` | bio_landcover_class, cropland_weight, lc_raw_value | ESA CCI 2000 |
| `npp_fire` | bio_npp_mean, bio_fire_freq | MODIS MOD17A3HGF, MODIS MCD64A1 burned area |
| `soil_chemical` | bio_biome_id, bio_soil_organic_carbon, bio_soil_ph | OpenLandMap SoilGrids, RESOLVE |
| `soil_texture` | bio_soil_sand, bio_soil_clay | OpenLandMap SoilGrids |
| `soil_physical` | bio_soil_bulk_density_fine, bio_soil_water_33kpa | OpenLandMap SoilGrids |
| `socioeconomic` | econ_dist_water, econ_road_density, paper_gdp_2015, paper_ghs_pop, paper_protected_binary | GRIP4, sat-io GDP, GHS-POP, WDPA |
| `highres_lights` | highres_lights_2012_2014 | VIIRS DNB annual V21 (median 2012–2014) |
| `highres_travel` | highres_travel | Oxford/MAP accessibility to cities 2015 |
| `highres_ghm` | highres_ghm | sat-io GHM 1990–2020, year 2000 (not CSP) |
| `highres_worldpop` | highres_worldpop_density | WorldPop GP/100m, year 2000 |
| `pasture` | gpw_raw_class | Global Pasture Watch v1, year 2000 |
| `opportunity_cost` | econ_opp_cost_usd_ha | SPAM agricultural value (uploaded asset) |

> **Climate:** All 19 raw bioclim variables are extracted here. PCA to PC1–PC4
> (matching Williams et al.) is computed in R using `prcomp()`.

> **Forest density and distance to forest** are NOT extracted in Step 2 — they require
> `fastDistanceTransform()` and `reduceNeighborhood()` which cannot be run as simple
> pixel lookups. These are computed in Step 3.

> **High-resolution layers** (`highres_lights`, `highres_ghm`, `highres_travel`,
> `highres_worldpop`) are each sampled as separate layer groups so that a mask gap
> in one (e.g., GHM over unmodified forest) does not cause `sampleRegions` to drop
> the entire feature, which would null out all columns at that point.

### Batching strategy

The nonregrowth collection (~800K points) causes GEE memory errors if processed as a
single collection. Instead, Step 2 uses `ee.Filter.eq('batch_id', n)` to filter to
one batch at a time server-side (a true GEE predicate, not a `toList()` slice which
would materialise the full collection).

```bash
# Regrowth — single collection, all layers
python brazil_step2_sample_layers.py --sample_type regrowth

# Nonregrowth — run one batch_id at a time (batch_id values 2-5)
python brazil_step2_sample_layers.py --sample_type nonregrowth --point-batch 2
python brazil_step2_sample_layers.py --sample_type nonregrowth --point-batch 3
python brazil_step2_sample_layers.py --sample_type nonregrowth --point-batch 4
python brazil_step2_sample_layers.py --sample_type nonregrowth --point-batch 5

# Annulus — single collection, all layers
python brazil_step2_sample_layers.py --sample_type annulus

# Run a single layer group only (useful for reruns)
python brazil_step2_sample_layers.py --sample_type regrowth --layer socioeconomic
python brazil_step2_sample_layers.py --sample_type regrowth --layer highres_ghm

# Climate was split further due to GEE memory limits — run parts separately
python brazil_step2_sample_layers.py --sample_type nonregrowth --layer climate_part1 --point-batch 2
```

### Output file naming

| Sample set | Pattern |
|---|---|
| regrowth | `Brazil_regrowth_{layer}.csv` |
| nonregrowth batch | `Brazil_nonregrowth_{layer}_pb{n}.csv` |
| annulus | `Brazil_annulus_{layer}.csv` |

All CSVs land in Google Drive root. Download to `./gee_outputs/` before Step 4.

---

## Step 3 — Distance & Density Variables (Tiled)

**Script:** `brazil_step3_forest_variables.py`

Computes spatially intensive variables (distance transforms, neighbourhood densities)
that cannot be done as simple pixel lookups in Step 2. All variables use 1.5°×1.5°
grid tiling with a 50km buffer for edge effects.

### Variables

| Variable key | Columns produced | Source |
|---|---|---|
| `forest` | bio_dist_forest_2000, bio_forest_density_1km2 | Hansen treecover2000 |
| `urban` | paper_dist_urban | ESA CCI 2000 class 190 |
| `settlement` | highres_dist_settlement, highres_worldpop_density_2km | WorldPop GP/100m year 2000 |
| `plantation` | bio_dist_plantation, bio_plantation_density_1km, bio_canopy_height_1km_mean | Xiao et al. planted forests, GEDI canopy height |
| `nightlight_density` | highres_nightlight_density_2km | VIIRS DNB median 2012–2014, 2km focal mean |
| `cropland` | highres_dist_cropland, highres_cropland_density_5km | ESA CCI 2000 (weighted classes 10/11/12/20/30/40) |
| `cultivated_grass` | highres_dist_cultivated_grass, highres_cultivated_grass_density_2km | Global Pasture Watch v1 year 2000 |

### Grid tile approach

`fastDistanceTransform()` and `reduceNeighborhood()` cannot run over all of Brazil
in a single GEE task. The script divides Brazil into 1.5°×1.5° tiles (~167km) with a
50km buffer for edge effects. One export task is submitted per tile containing points.
Each task computes the distance transform and focal mean over its buffered tile, then
samples all points within the strict tile boundary.

```bash
# Run all variables for each sample type
python brazil_step3_forest_variables.py --sample_type regrowth
python brazil_step3_forest_variables.py --sample_type nonregrowth
python brazil_step3_forest_variables.py --sample_type annulus

# Run a single variable
python brazil_step3_forest_variables.py --sample_type regrowth --variable urban
python brazil_step3_forest_variables.py --sample_type regrowth --variable nightlight_density
python brazil_step3_forest_variables.py --sample_type regrowth --variable plantation

# Debugging
python brazil_step3_forest_variables.py --sample_type regrowth --variable urban --dry-run
python brazil_step3_forest_variables.py --sample_type regrowth --variable forest --test-one
python brazil_step3_forest_variables.py --sample_type regrowth --variable forest --tile -63 -6

# Fix "Reprojection output too large" on failed plantation tiles
python brazil_step3_forest_variables.py --sample_type annulus --variable plantation --tile -36.5 -5.5 --buffer 25000
```

> **Plantation note:** The Xiao et al. dataset is natively 30m. No `.reproject()` calls
> are used — `connectedPixelCount` and `fastDistanceTransform` operate at native
> resolution through lazy evaluation. Forcing `.reproject()` caused "Reprojection output
> too large" errors. Use `--buffer 25000` to shrink the buffer on failing tiles.

**Output folders:** `GEE_Brazil_{variable}/` in Google Drive  
**File naming:** `Brazil_{variable}_{sample_type}_Lon{w}_Lat{s}.csv`

Download all tile CSVs to `./gee_outputs/GEE_Brazil_{variable}/` before Step 4.

---

## Step 4 — Merge Into Final Dataset

**Script:** `brazil_step4_merge.py`

Merges all Step 2 and Step 3 CSVs into a single analysis-ready file.

### What it does

1. For each sample set (regrowth, nonregrowth, annulus), loads all Step 2 layer CSVs
   and left-joins them on `point_id`
2. Parses `.geo` JSON column → `longitude`, `latitude`
3. Deduplicates on `(longitude, latitude)` — GEE `stratifiedSample` places points at
   pixel centroids, so two draws on the same pixel produce bit-identical coordinates
4. Adds `y` (1=regrowth, 0=nonregrowth/annulus) and `sample_set` columns
5. Loads and joins all Step 3 tiled variable folders (deduplicating tile-boundary overlaps)
6. Stacks all three sample sets and writes a single CSV

### Step 3 tiled layers joined

| Variable | Folder | Glob pattern | Columns |
|---|---|---|---|
| forest | `GEE_Brazil_forest_2000/` | `Brazil_forest2000_{sample_set}_*.csv` | bio_dist_forest_2000, bio_forest_density_1km2 |
| urban | `GEE_Brazil_urban/` | `Brazil_urban_{sample_set}_*.csv` | paper_dist_urban |
| settlement | `GEE_Brazil_settlement/` | `Brazil_settlement_{sample_set}_*.csv` | highres_dist_settlement, highres_worldpop_density_2km |
| plantation | `GEE_Brazil_plantation/` | `Brazil_plantation_{sample_set}_*.csv` | bio_dist_plantation, bio_plantation_density_1km, bio_canopy_height_1km_mean |
| cultivated_grass | `GEE_Brazil_cultivated_grass/` | `Brazil_cultivated_grass_{sample_set}_*.csv` | highres_dist_cultivated_grass, highres_cultivated_grass_density_2km |
| cropland | `GEE_Brazil_cropland/` | `Brazil_cropland_{sample_set}_*.csv` | highres_dist_cropland, highres_cropland_density_5km |
| nightlight_density | `GEE_Brazil_nightlight_density/` | `Brazil_nightlight_density_{sample_set}_*.csv` | highres_nightlight_density_2km |

```bash
# Full merge including all tiled variables
python brazil_step4_merge.py --gee-dir ./gee_outputs

# Skip all tiled (step3) variables
python brazil_step4_merge.py --gee-dir ./gee_outputs --skip-tiles

# Skip only forest tiles (backward compat)
python brazil_step4_merge.py --gee-dir ./gee_outputs --skip-forest

# Custom output path
python brazil_step4_merge.py --gee-dir ./gee_outputs --output brazil_training_data.csv
```

### Output columns

| Group | Columns | Description |
|---|---|---|
| **ID** | point_id, sample_set, y, batch_id, longitude, latitude | Identifiers and outcome |
| **Climate** | bio01–bio19 | WorldClim v2.1 raw bioclim variables (PCA in R) |
| **Topography** | bio_elevation, bio_slope | SRTM elevation (m) and slope (degrees) |
| **Land cover** | bio_landcover_class, cropland_weight, lc_raw_value | ESA CCI 2000 reclassified to 11 classes + raw |
| **Productivity** | bio_npp_mean, bio_fire_freq | MODIS NPP 2001–2017 mean, fire frequency 2001–2016 |
| **Soil** | bio_biome_id, bio_soil_organic_carbon, bio_soil_ph, bio_soil_sand, bio_soil_clay, bio_soil_bulk_density_fine, bio_soil_water_33kpa | RESOLVE biome + OpenLandMap top 30cm |
| **Forest** | bio_dist_forest_2000, bio_forest_density_1km2 | Hansen treecover2000 distance (m) and 1km² mean density |
| **Plantation** | bio_dist_plantation, bio_plantation_density_1km, bio_canopy_height_1km_mean | Xiao planted forests + GEDI canopy height |
| **Econ (base)** | econ_dist_water, econ_road_density | Distance to water, GRIP4 road density |
| **Econ (paper)** | paper_gdp_2015, paper_ghs_pop, paper_protected_binary, paper_dist_urban | GDP, GHS population, WDPA protected areas, distance to urban |
| **High-res** | highres_lights_2012_2014, highres_travel, highres_ghm, highres_worldpop_density, highres_dist_settlement, highres_worldpop_density_2km, highres_dist_cultivated_grass, highres_cultivated_grass_density_2km, highres_dist_cropland, highres_cropland_density_5km, highres_nightlight_density_2km | VIIRS nightlights, travel time, GHM year 2000, WorldPop, settlement distance, pasture distance/density, cropland distance/density, nightlight density |
| **Pasture** | gpw_raw_class | Global Pasture Watch grassland class 2000 |
| **Opp. cost** | econ_opp_cost_usd_ha | Agricultural opportunity cost (SPAM) |

### Expected output size

| Sample set | Rows |
|---|---|
| regrowth | ~199,000 |
| nonregrowth | ~797,000 |
| annulus | ~199,700 |
| **total** | **~1,196,000** |

---

## Next Steps — R Analysis

### 1. Bioclim PCA

Williams et al. derive PC1–PC4 from all 19 bioclim variables using `prcomp()` on a
1M random tropical point sample, then apply the rotation to training points.

```r
library(tidyverse)

df <- read_csv('gee_outputs/brazil_training_data.csv')

# Drop legacy variable if present
df <- df %>% select(-any_of('bio_forest_density_2018'))

# Fit PCA on a sample (or full dataset)
bio_cols <- paste0('bio', str_pad(1:19, 2, pad='0'))
pca_fit  <- prcomp(df[, bio_cols], center=TRUE, scale.=TRUE)

# Apply to all points
pcs      <- predict(pca_fit, df[, bio_cols])
df       <- bind_cols(df, as_tibble(pcs[, 1:4]) %>%
              rename(bioclim_pc1=PC1, bioclim_pc2=PC2, bioclim_pc3=PC3, bioclim_pc4=PC4))
```

### 2. Analysis 1 — Prevalence sensitivity

Vary the Y=1:Y=0 ratio from 1:1 to 1:8 using the Brazil-wide nonregrowth points,
fit a random forest for each ratio, and compare predicted regeneration area estimates.

### 3. Analysis 2 — Climate envelope test

**Run 1 (Williams baseline):** 50/50 split, nonregrowth 100% Brazil-wide.  
**Run 2 (Climate-controlled):** 50/50 split, nonregrowth 80% annulus + 20% Brazil-wide.

If predicted area or spatial distribution shifts substantially between runs, this
provides evidence that the original model was partly learning climate envelope rather
than land-use-driven regeneration signal.

---

## File Structure

```
brazil_predict/
├── brazil_step0_export_masks.py
├── brazil_step0b_export_annuli.py
├── brazil_step1_sample_points.py
├── brazil_step2_sample_layers.py
├── brazil_step3_forest_variables.py
├── brazil_step4_merge.py
└── gee_outputs/
    ├── Brazil_regrowth_*.csv                 # Step 2 regrowth layers
    ├── Brazil_nonregrowth_*_pb*.csv          # Step 2 nonregrowth layers (batched)
    ├── Brazil_annulus_*.csv                   # Step 2 annulus layers
    ├── GEE_Brazil_forest_2000/               # Step 3 forest tiles
    │   └── Brazil_forest2000_*.csv
    ├── GEE_Brazil_urban/                     # Step 3 urban distance tiles
    │   └── Brazil_urban_*.csv
    ├── GEE_Brazil_settlement/                # Step 3 settlement tiles
    │   └── Brazil_settlement_*.csv
    ├── GEE_Brazil_plantation/                # Step 3 plantation tiles
    │   └── Brazil_plantation_*.csv
    ├── GEE_Brazil_nightlight_density/        # Step 3 nightlight density tiles
    │   └── Brazil_nightlight_density_*.csv
    ├── GEE_Brazil_cropland/                  # Step 3 cropland tiles
    │   └── Brazil_cropland_*.csv
    ├── GEE_Brazil_cultivated_grass/          # Step 3 pasture tiles
    │   └── Brazil_cultivated_grass_*.csv
    └── brazil_training_data.csv              # Final merged dataset (Step 4 output)
```
