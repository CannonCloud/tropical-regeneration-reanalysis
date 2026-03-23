# ============================================================
# 04_socioeconomic_models.R
# Brazil Natural Regeneration — Socioeconomic Model Extension
#
# Standalone script. Loads saved objects from 01_train_models.R
# pipeline and produces:
#   rf_socio_50_50.rds         — baseline (50/50) with socioeconomic vars
#   rf_socio_paired.rds        — paired (50/40/10) with socioeconomic vars
#   socio_vars.rds             — variable list
#   fig_supp_importance_socio  — variable importance shift plot
#   fig_supp_pdp_socio         — PDPs for key socioeconomic variables
#
# Prereqs on disk (from 01_train_models.R):
#   gee_outputs/brazil_training_data.csv
#   rf_outputs/pca_fit.rds
#   rf_outputs/williams_vars.rds
#   rf_outputs/extended_vars.rds
#   rf_outputs/negatives_paired_p1.rds
#   rf_outputs/pred_sample.rds
#
# Prereqs on disk (for sensitivity/specificity comparison):
#   rf_outputs/preds_rf_prev_50_50.rds
#   rf_outputs/preds_rf_ce.rds
#   rf_outputs/preds_rf_ext_50_50.rds
#   rf_outputs/preds_rf_ext_paired.rds
#
# Run time: ~5 min (data load + 2 RF fits + PDPs + predictions)
# ============================================================

library(data.table)
library(ranger)
options(ranger.num.threads = parallel::detectCores() - 3)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("# set your working directory here")

output_dir <- "rf_outputs"
dir.create(output_dir, showWarnings = FALSE)


# ============================================================
# S1: Load data & project PCA
# ============================================================
# Mirrors R1+R2 of 01_train_models.R but loads the SAVED PCA
# (no refitting) so PC scores are identical to the main pipeline.

message("S1: Loading data...")
dt <- fread("gee_outputs/brazil_training_data.csv")
message(sprintf("Loaded: %d rows x %d cols", nrow(dt), ncol(dt)))

# Clean factors (same as R1)
if ("bio_forest_density_2018" %in% names(dt))
  dt[, bio_forest_density_2018 := NULL]
dt[, bio_landcover_class := as.factor(bio_landcover_class)]
dt[, bio_biome_id        := as.factor(bio_biome_id)]
dt[, y                   := as.factor(y)]

print(dt[, .N, by = sample_set])

# Project PCA (same as R2, using saved fit)
pca_fit  <- readRDS(file.path(output_dir, "pca_fit.rds"))
bio_cols <- paste0("bio", sprintf("%02d", 1:19))
pcs      <- predict(pca_fit, dt[, ..bio_cols])
dt[, bioclim_pc1 := pcs[, 1]]
dt[, bioclim_pc2 := pcs[, 2]]
dt[, bioclim_pc3 := pcs[, 3]]
dt[, bioclim_pc4 := pcs[, 4]]
dt[, bioclim_pc5 := pcs[, 5]]
rm(pcs); gc()

# Load variable lists from pipeline
williams_vars <- readRDS(file.path(output_dir, "williams_vars.rds"))
extended_vars <- readRDS(file.path(output_dir, "extended_vars.rds"))

# Drop rows missing Williams vars (same filter as 01_train_models.R R2)
dt <- dt[complete.cases(dt[, ..williams_vars])]
message(sprintf("Rows after Williams NA drop: %d", nrow(dt)))

# Training budget (same as 01_train_models.R)
n_regrowth_all <- dt[sample_set == "regrowth", .N]
n_total        <- 2 * n_regrowth_all
message(sprintf("n_total = %d (n_regrowth = %d)", n_total, n_regrowth_all))


# ============================================================
# S2: Define socioeconomic variable set
# ============================================================
# extended_vars (18 biophysical) + socioeconomic additions.
#
# Design principles:
#   - Keep bio_landcover_class: it's in the Williams baseline and
#     captures pixel-level land state, complementary to neighborhood
#     distance/density variables.
#   - Drop plantation vars: Xiao et al. dataset maps 2021 extent,
#     so it leaks future information into a year-2000 training model.
#     PDPs confirmed it mirrors forest density (updated tree line).
#   - Drop bio_canopy_height_1km_mean: ETH/GEDI product is c.2020,
#     same future-leak problem. Also largely biophysical (Amazon proxy).
#   - Hardcode highres_ghm NAs → 0: missings are deep forest in 2000
#     with no human modification signal.
#   - Variables are grouped by concept so importance can be mentally
#     aggregated (e.g., settlement distance + population density both
#     measure "proximity to people").
#
# VARIABLE STATUS:
#   [ACTIVE]  = in the model now
#   [PENDING] = GEE job running, uncomment when available
#   [FUTURE]  = needs re-run (e.g. pasture for year 2000)
# ============================================================

# ---- Hardcode GHM missings to 0 before defining vars ---------
if ("highres_ghm" %in% names(dt)) {
  n_ghm_na <- sum(is.na(dt$highres_ghm))
  dt[is.na(highres_ghm), highres_ghm := 0]
  message(sprintf("  GHM: filled %d NAs with 0 (%.1f%% of rows)",
                  n_ghm_na, 100 * n_ghm_na / nrow(dt)))
}

socio_vars <- c(
  # ── Extended biophysical (18 vars, from 01_train_models.R R6) ──
  extended_vars,
  
  # ── Human settlement & population ──────────────────────────────
  "highres_worldpop_density_2km",   # [ACTIVE]  WorldPop 2015, 2km focal mean
  "highres_dist_settlement",        # [ACTIVE]  distance to settlement (WorldPop)
  #"paper_dist_urban",               # [Redundant]  distance to ESA CCI urban (2000)
  
  # ── Infrastructure & economy ───────────────────────────────────
  "econ_road_density",              # [ACTIVE]  GRIP4 road network density
  "econ_dist_water",                # [ACTIVE]  distance to water bodies
  "paper_gdp_2015",                 # [ACTIVE]  GDP (Kummu et al. 2018)
  "paper_protected_binary",         # [ACTIVE]  WDPA protected area (binary)
  
  # ── Agriculture & land use ─────────────────────────────────────
  "highres_dist_cropland",          # [ACTIVE]  distance to ESA CCI cropland (2000)
  # "cropland_weight",              # [CONSIDER] pixel-level ag intensity (ESA CCI 2000)
  # "highres_cropland_density_5km", # [PENDING]  weighted cropland fraction in 5km radius
  "highres_dist_cultivated_grass",       # [ACTIVE] Global Pasture Watch — needs 2000 rerun
  "highres_cultivated_grass_density_2km",# [ACTIVE] Global Pasture Watch — needs 2000 rerun
  
  # ── Accessibility & human footprint ────────────────────────────
  "highres_ghm",                  # [ACTIVE]  Global Human Modification index (NAs → 0)
  "highres_travel",               # [ACTIVE]  travel time to nearest city (MAP 2015)
  "highres_nightlight_density_2km",     # [ACTIVE]  VIIRS nightlights density
  
  # ── DROPPED ────────────────────────────────────────────────────
  # highres_lights_2012_2014      — dropped due to overlap with density
  # bio_plantation_density_1km    — Xiao et al. maps 2021 extent (future leak)
  # bio_dist_plantation           — same dataset, same problem
  # bio_canopy_height_1km_mean    — ETH/GEDI c.2020 (future leak, also biophysical)
  # paper_ghs_pop                 — lower-res duplicate of worldpop density
  
  NULL  # trailing NULL allows clean commenting
)

# Remove the trailing NULL
socio_vars <- socio_vars[!is.na(socio_vars) & socio_vars != ""]

# ---- Verify all exist in dt ----------------------------------
missing_socio <- socio_vars[!socio_vars %in% names(dt)]
if (length(missing_socio) > 0) {
  warning("Missing socio vars (dropped): ", paste(missing_socio, collapse = ", "))
  socio_vars <- socio_vars[socio_vars %in% names(dt)]
}

cat("\nSocioeconomic variable set:\n")
cat(sprintf("  %d total (%d extended biophysical + %d socioeconomic)\n",
            length(socio_vars), length(extended_vars),
            length(socio_vars) - length(extended_vars)))
cat(paste("  ", socio_vars, collapse = "\n"), "\n")

# ---- Drop NAs ------------------------------------------------
dt_socio <- dt[complete.cases(dt[, ..socio_vars])]
message(sprintf("Rows after socio NA drop: %d (lost %d)",
                nrow(dt_socio), nrow(dt) - nrow(dt_socio)))

saveRDS(socio_vars, file.path(output_dir, "socio_vars.rds"))


# ============================================================
# S3: Train models
# ============================================================

# ---- S3a: Socioeconomic baseline (50/50) ---------------------

message("\n--- S3a: Fitting socioeconomic RF: 50/50 baseline ---")

set.seed(42)
n_pos_s <- min(round(n_total * 0.50), dt_socio[sample_set == "regrowth",    .N])
n_neg_s <- min(n_total - round(n_total * 0.50), dt_socio[sample_set == "nonregrowth", .N])

train_pos_s <- dt_socio[sample_set == "regrowth"]    |> _[sample(.N, n_pos_s)]
train_neg_s <- dt_socio[sample_set == "nonregrowth"] |> _[sample(.N, n_neg_s)]
train_s     <- rbind(train_pos_s, train_neg_s)
rm(train_pos_s, train_neg_s); gc()

message(sprintf("  Training: %d pos, %d neg (%.1f%% positive)",
                n_pos_s, n_neg_s, 100 * n_pos_s / (n_pos_s + n_neg_s)))

t0 <- proc.time()
rf_socio_50_50 <- ranger(
  formula     = y ~ .,
  data        = train_s[, c("y", socio_vars), with = FALSE],
  num.trees   = 500,
  mtry        = floor(sqrt(length(socio_vars))),
  probability = TRUE,
  importance  = "permutation",
  seed        = 42
)
message(sprintf("  Fitted in %.1f seconds", (proc.time() - t0)["elapsed"]))
rm(train_s); gc()

saveRDS(rf_socio_50_50, file.path(output_dir, "rf_socio_50_50.rds"))
message("  Saved: rf_socio_50_50.rds")
rm(rf_socio_50_50); gc()

# ---- S3b: Socioeconomic paired (50/40/10, pass 1) -----------

message("\n--- S3b: Fitting socioeconomic RF: paired (50/40/10) ---")

# Load pass-1 matched negatives and re-hydrate with current columns
negatives_paired_raw <- readRDS(file.path(output_dir, "negatives_paired_p1.rds"))
matched_ids <- negatives_paired_raw[, .(point_id)]
negatives_paired <- dt[matched_ids, on = "point_id", nomatch = NULL]
message(sprintf("  Re-hydrated matched negatives: %d / %d rows joined",
                nrow(negatives_paired), nrow(matched_ids)))
rm(negatives_paired_raw, matched_ids); gc()

negatives_paired_socio <- negatives_paired[
  complete.cases(negatives_paired[, ..socio_vars])
]
message(sprintf("  Matched negatives with complete socio vars: %d / %d",
                nrow(negatives_paired_socio), nrow(negatives_paired)))

set.seed(42)
n_pos_sp  <- min(round(n_total * 0.50), dt_socio[sample_set == "regrowth",    .N])
n_pair_sp <- min(round(n_total * 0.40), nrow(negatives_paired_socio))
n_wide_sp <- min(n_total - round(n_total * 0.50) - round(n_total * 0.40),
                 dt_socio[sample_set == "nonregrowth", .N])

train_pos_sp  <- dt_socio[sample_set == "regrowth"]    |> _[sample(.N, n_pos_sp)]
train_pair_sp <- negatives_paired_socio[sample(.N, n_pair_sp)]
train_wide_sp <- dt_socio[sample_set == "nonregrowth"] |> _[sample(.N, n_wide_sp)]
train_sp      <- rbind(train_pos_sp, train_pair_sp, train_wide_sp)
rm(train_pos_sp, train_pair_sp, train_wide_sp, negatives_paired_socio); gc()

message(sprintf("  Training: %d regrowth, %d matched, %d wide",
                n_pos_sp, n_pair_sp, n_wide_sp))

t0 <- proc.time()
rf_socio_paired <- ranger(
  formula     = y ~ .,
  data        = train_sp[, c("y", socio_vars), with = FALSE],
  num.trees   = 500,
  mtry        = floor(sqrt(length(socio_vars))),
  probability = TRUE,
  importance  = "permutation",
  seed        = 42
)
message(sprintf("  Fitted in %.1f seconds", (proc.time() - t0)["elapsed"]))
rm(train_sp); gc()

saveRDS(rf_socio_paired, file.path(output_dir, "rf_socio_paired.rds"))
message("  Saved: rf_socio_paired.rds")
rm(rf_socio_paired); gc()

message("\nS3 complete. Both socioeconomic models saved.")


# ============================================================
# S4: Variable importance figure
# ============================================================

message("\n--- S4: Variable importance ---")

compute_rel_imp <- function(rf, model_label) {
  vi    <- rf$variable.importance
  total <- sum(vi)
  data.table(
    model    = model_label,
    variable = names(vi),
    rel_imp  = 100 * vi / total
  )
}

nature_theme <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(0.35, "cm"),
    plot.title       = element_text(face = "bold", size = 9),
    axis.title       = element_text(size = 8),
    legend.text      = element_text(size = 7)
  )

col_baseline <- "#E41A1C"
col_ce       <- "#377EB8"
col_socio    <- "#4DAF4A"

rf_tmp          <- readRDS(file.path(output_dir, "rf_socio_50_50.rds"))
vi_socio_base   <- compute_rel_imp(rf_tmp, "50/50 baseline")
rm(rf_tmp); gc()

rf_tmp           <- readRDS(file.path(output_dir, "rf_socio_paired.rds"))
vi_socio_paired  <- compute_rel_imp(rf_tmp, "Paired")
rm(rf_tmp); gc()

# ---- Clean display names -------------------------------------
# Add entries here as new variables are activated

clean_names_socio <- c(
  # Williams vars
  bio_forest_density_1km2         = "Forest density",
  bioclim_pc1                     = "Bioclim PC1",
  bio_dist_forest_2000            = "Dist. to forest",
  bio_soil_ph                     = "Soil pH",
  bioclim_pc2                     = "Bioclim PC2",
  bioclim_pc3                     = "Bioclim PC3",
  bio_landcover_class             = "Landcover",
  bioclim_pc4                     = "Bioclim PC4",
  bio_soil_organic_carbon         = "Soil org. C",
  bio_biome_id                    = "Biome",
  # Extended biophysical additions
  bio_npp_mean                    = "NPP",
  bio_fire_freq                   = "Fire frequency",
  bio_slope                       = "Slope",
  bio_elevation                   = "Elevation",
  bioclim_pc5                     = "Bioclim PC5",
  bio_soil_sand                   = "Soil sand",
  bio_soil_clay                   = "Soil clay",
  bio_soil_bulk_density_fine      = "Soil bulk density",
  bio_soil_water_33kpa            = "Soil water capacity",
  # Socioeconomic — settlement & population
  highres_worldpop_density_2km    = "Population density",
  highres_dist_settlement         = "Dist. to settlement",
  paper_dist_urban                = "Dist. to urban",
  # Socioeconomic — infrastructure & economy
  econ_road_density               = "Road density",
  econ_dist_water                 = "Dist. to water",
  paper_gdp_2015                  = "GDP (2015)",
  paper_protected_binary          = "Protected area",
  # Socioeconomic — agriculture & land use
  highres_dist_cropland           = "Dist. to cropland",
  cropland_weight                 = "Cropland intensity",
  highres_cropland_density_2km    = "Cropland density",
  highres_dist_cultivated_grass   = "Dist. to pasture",
  highres_cultivated_grass_density_2km = "Pasture density",
  # Socioeconomic — accessibility & human footprint
  highres_ghm                     = "Human modification",
  highres_travel                  = "Travel time",
  highres_nightlight_density_2km  = "Nightlight density"
)

vi_socio_both <- rbind(vi_socio_base, vi_socio_paired)
vi_socio_both[, var_label := clean_names_socio[variable]]
vi_socio_both[, model     := factor(model, levels = c("50/50 baseline", "Paired"))]

var_order_socio <- vi_socio_base[order(-rel_imp), clean_names_socio[variable]]
vi_socio_both[, var_label := factor(var_label, levels = rev(var_order_socio))]

# ---- Direction flags -----------------------------------------
#   "down"    = between-climate classifiers (expect ↓ under paired)
#   "up"      = local biophysical (expect ↑ under paired)
#   "socio"   = socioeconomic (the new test)
#   "neutral" = other (PCs 2-5, landcover)

socio_labels <- c(
  "Population density", "Dist. to settlement", "Dist. to urban",
  "Road density", "Dist. to water", "GDP (2015)", "Protected area",
  "Dist. to cropland", "Cropland intensity", "Cropland density",
  "Dist. to pasture", "Pasture density",
  "Human modification", "Travel time", "Nightlight density"
)

vi_socio_both[, direction := fifelse(
  var_label %in% c("Bioclim PC1", "Biome", "Soil pH", "Elevation"), "down",
  fifelse(var_label %in% c("Forest density", "Dist. to forest",
                           "Soil org. C", "NPP", "Soil bulk density",
                           "Soil sand", "Fire frequency", "Slope"),
          "up",
          fifelse(var_label %in% socio_labels,
                  "socio", "neutral"))
)]

vi_socio_labels <- vi_socio_both[model == "Paired"]

# Dynamic figure height: ~6mm per variable + 30mm overhead
n_vars   <- length(unique(vi_socio_both$variable))
fig_h_mm <- max(120, n_vars * 6 + 30)

pB_socio <- ggplot(vi_socio_both, aes(x = model, y = rel_imp, group = var_label)) +
  geom_line(aes(colour = direction), linewidth = 0.7, alpha = 0.85) +
  geom_point(aes(colour = direction), size = 2) +
  scale_colour_manual(
    values = c("down"    = col_baseline,
               "up"      = col_ce,
               "socio"   = col_socio,
               "neutral" = "grey70"),
    labels = c("down"    = "Climate classifiers (\u2193)",
               "up"      = "Local biophysical (\u2191)",
               "socio"   = "Socioeconomic",
               "neutral" = "Other"),
    name   = NULL
  ) +
  geom_text_repel(
    data               = vi_socio_labels,
    aes(label = var_label, colour = direction),
    hjust              = 0,
    direction          = "y",
    nudge_x            = 0.08,
    segment.size       = 0.3,
    segment.alpha      = 0.5,
    size               = 2.2,
    show.legend        = FALSE,
    force              = 3,
    min.segment.length = 0,
    max.overlaps       = 30
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.55))) +
  labs(
    title = "Variable importance: socioeconomic model (baseline vs paired)",
    x     = NULL,
    y     = "Relative importance (% of total)"
  ) +
  nature_theme +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

ggsave(file.path(output_dir, "fig_supp_importance_socio.pdf"),
       pB_socio, width = 135, height = fig_h_mm, units = "mm",
       device = cairo_pdf)
ggsave(file.path(output_dir, "fig_supp_importance_socio.png"),
       pB_socio, width = 135, height = fig_h_mm, units = "mm", dpi = 300)
message(sprintf("Saved: fig_supp_importance_socio.pdf/.png (%d mm tall)", fig_h_mm))

# ---- Print importance shift table ----------------------------
cat("\nSocioeconomic variable importance shifts (baseline \u2192 paired):\n")
imp_shift <- merge(
  vi_socio_base[,  .(variable, var_label = clean_names_socio[variable], base = rel_imp)],
  vi_socio_paired[, .(variable, paired = rel_imp)],
  by = "variable"
)
imp_shift[, shift := paired - base]
imp_shift <- imp_shift[order(-abs(shift))]
print(imp_shift[, .(var_label, base = round(base, 2),
                    paired = round(paired, 2),
                    shift  = round(shift, 2))])


# ============================================================
# S5: Partial dependence plots — socioeconomic variables
# ============================================================
# Sweeps each variable across its 2nd–98th percentile range
# while holding all others at observed values (marginal PDP).
# Background sample: 3000 nonregrowth rows with complete socio vars.

message("\n--- S5: Partial dependence plots ---")

# ---- Define which variables to PDP ---------------------------
# Update this list as new variables are activated

pdp_socio_vars <- c(
  "highres_cultivated_grass_density_2km",
  "highres_dist_cultivated_grass",
  "highres_dist_cropland",
  "highres_ghm",                  
  "highres_travel",
  "highres_worldpop_density_2km"
)

pdp_socio_labels <- c(
  highres_cultivated_grass_density_2km = "Pasture density (2km)",
  highres_dist_cultivated_grass        = "Distance to pasture (m)",
  highres_dist_cropland                = "Distance to cropland (m)",
  highres_ghm                          = "Human modification index",
  highres_travel                       = "Travel time to city (min)",
  highres_worldpop_density_2km         = "Population density (2km)"
)

# ---- Background sample ---------------------------------------
set.seed(789)
pdp_bg_socio <- dt_socio[sample_set == "nonregrowth"][sample(.N, min(.N, 3000))]
pdp_bg_socio <- pdp_bg_socio[complete.cases(pdp_bg_socio[, ..socio_vars])]
message(sprintf("PDP background sample: %d rows", nrow(pdp_bg_socio)))

# ---- Compute -------------------------------------------------

compute_pdp_socio <- function(rf, data, vars_to_sweep, use_vars,
                              n_grid = 35, n_sample = 3000) {
  set.seed(123)
  data_s <- data[sample(min(.N, n_sample))]
  
  vars_to_cap <- c(
    "highres_dist_cropland",
    "highres_travel",
    "highres_dist_cultivated_grass",
    "highres_worldpop_density_2km"
  )
  
  rbindlist(lapply(vars_to_sweep, function(var) {
    
    if (var %in% vars_to_cap) {
      x_min <- min(data_s[[var]], na.rm = TRUE)
      x_max <- quantile(data_s[[var]], 0.90, na.rm = TRUE)
    } else {
      x_min <- quantile(data_s[[var]], 0.02, na.rm = TRUE)
      x_max <- quantile(data_s[[var]], 0.98, na.rm = TRUE)
    }
    
    grid_vals <- seq(x_min, x_max, length.out = n_grid)
    
    pdp_vals <- sapply(grid_vals, \(val) {
      data_tmp <- copy(data_s)
      data_tmp[, (var) := val]
      mean(predict(rf, data = data_tmp[, ..use_vars])$predictions[, "1"])
    })
    
    data.table(variable = var, x = grid_vals, prob = pdp_vals)
  }))
}

message("  PDP: loading rf_socio_50_50...")
rf_tmp             <- readRDS(file.path(output_dir, "rf_socio_50_50.rds"))
pdp_socio_baseline <- compute_pdp_socio(rf_tmp, pdp_bg_socio,
                                        pdp_socio_vars, socio_vars)
pdp_socio_baseline[, model := "50/50 baseline"]
rm(rf_tmp); gc()

message("  PDP: loading rf_socio_paired...")
rf_tmp            <- readRDS(file.path(output_dir, "rf_socio_paired.rds"))
pdp_socio_paired_dt <- compute_pdp_socio(rf_tmp, pdp_bg_socio,
                                         pdp_socio_vars, socio_vars)
pdp_socio_paired_dt[, model := "Paired"]
rm(rf_tmp, pdp_bg_socio); gc()

pdp_socio_all <- rbind(pdp_socio_baseline, pdp_socio_paired_dt)
pdp_socio_all[, model := factor(model, levels = c("50/50 baseline", "Paired"))]

# ---- Plot ----------------------------------------------------

make_pdp_socio_panel <- function(var, label, ylim_max = NA) {
  ggplot(pdp_socio_all[variable == var],
         aes(x = x, y = prob * 100, colour = model, linetype = model)) +
    geom_line(linewidth = 1.0) +
    { if (!is.na(ylim_max)) coord_cartesian(ylim = c(0, ylim_max)) } +
    scale_colour_manual(values = c("50/50 baseline" = col_baseline,
                                   "Paired"         = col_ce)) +
    scale_linetype_manual(values = c("50/50 baseline" = "solid",
                                     "Paired"         = "dashed")) +
    labs(title = label, x = NULL, y = NULL, colour = NULL, linetype = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position  = "none",
          plot.title       = element_text(size = 8, face = "bold"),
          panel.grid.minor = element_blank())
}

# Build panels dynamically from pdp_socio_vars
panels_socio <- lapply(pdp_socio_vars, function(var) {
  label <- pdp_socio_labels[var]
  if (is.na(label)) label <- var
  make_pdp_socio_panel(var, label)
})

# Layout: auto ncol based on count
n_pdp    <- length(panels_socio)
pdp_ncol <- min(n_pdp, 3)
pdp_nrow <- ceiling(n_pdp / pdp_ncol)

# Shared legend
legend_socio <- ggplot(pdp_socio_all[variable == pdp_socio_vars[1]],
                       aes(x = x, y = prob * 100, colour = model, linetype = model)) +
  geom_line() +
  scale_colour_manual(values = c("50/50 baseline" = col_baseline, "Paired" = col_ce)) +
  scale_linetype_manual(values = c("50/50 baseline" = "solid", "Paired" = "dashed")) +
  labs(colour = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

shared_legend_socio <- cowplot::get_legend(legend_socio)

fig_pdp_socio <- wrap_plots(panels_socio, ncol = pdp_ncol) /
  wrap_elements(shared_legend_socio) +
  plot_layout(heights = c(20, 1)) +
  plot_annotation(
    caption = "x: variable value    y: mean predicted P(regrowth) (%)"
  )

pdp_w <- min(180, pdp_ncol * 65)
pdp_h <- pdp_nrow * 70 + 15

ggsave(file.path(output_dir, "fig_supp_pdp_socio.pdf"),
       fig_pdp_socio, width = pdp_w, height = pdp_h, units = "mm",
       device = cairo_pdf)
ggsave(file.path(output_dir, "fig_supp_pdp_socio.png"),
       fig_pdp_socio, width = pdp_w, height = pdp_h, units = "mm", dpi = 300)
message(sprintf("Saved: fig_supp_pdp_socio.pdf/.png (%dx%d mm)", pdp_w, pdp_h))


# ============================================================
# S6: Sensitivity / specificity
# ============================================================
# Predicts on all regrowth pixels and the pred_sample (nonregrowth)
# to compute whole-Brazil sensitivity and specificity at threshold 0.5.
# Comparison models loaded from prediction caches (02_figures.R).

message("\n--- S6: Sensitivity / specificity ---")

# ---- Build prediction surfaces for socio models --------------
reg_pixels_socio <- dt[sample_set == "regrowth"]
reg_pixels_socio <- reg_pixels_socio[complete.cases(reg_pixels_socio[, ..socio_vars])]
message(sprintf("  Regrowth pixels (socio vars): %d", nrow(reg_pixels_socio)))

# Re-hydrate pred_sample with current columns
pred_sample_raw   <- readRDS(file.path(output_dir, "pred_sample.rds"))
pred_sample_ids   <- pred_sample_raw[, .(point_id)]
pred_sample_socio <- dt[pred_sample_ids, on = "point_id", nomatch = NULL]
pred_sample_socio <- pred_sample_socio[complete.cases(pred_sample_socio[, ..socio_vars])]
message(sprintf("  Pred sample (socio vars): %d / %d re-hydrated",
                nrow(pred_sample_socio), nrow(pred_sample_ids)))
rm(pred_sample_raw, pred_sample_ids); gc()

# ---- Print comparison table ----------------------------------

cat(sprintf("\n%-30s  %11s  %11s  %9s  %11s\n",
            "Model", "Sensitivity", "Specificity", "n_reg", "n_nonreg"))
cat(strrep("-", 78), "\n")

# Reference models from caches (02_figures.R)
ref_cache_specs <- list(
  list(cache = "preds_rf_prev_50_50.rds",   lbl = "Williams 50/50"),
  list(cache = "preds_rf_ce.rds",           lbl = "Williams paired"),
  list(cache = "preds_rf_ext_50_50.rds",    lbl = "Extended 50/50"),
  list(cache = "preds_rf_ext_paired.rds",   lbl = "Extended paired")
)

for (ref in ref_cache_specs) {
  cache_path <- file.path(output_dir, ref$cache)
  if (file.exists(cache_path)) {
    cache <- readRDS(cache_path)
    sens     <- mean(cache$reg    > 0.5)
    spec_val <- mean(cache$nonreg <= 0.5)
    cat(sprintf("%-30s  %11.4f  %11.4f  %9d  %11d\n",
                ref$lbl, sens, spec_val,
                length(cache$reg), length(cache$nonreg)))
  }
}

cat(strrep("-", 78), "\n")

# New socio models (predict directly)
socio_model_specs <- list(
  list(rds = "rf_socio_50_50.rds",  lbl = "Socio 50/50 baseline"),
  list(rds = "rf_socio_paired.rds", lbl = "Socio paired (50/40/10)")
)

for (spec in socio_model_specs) {
  rf <- readRDS(file.path(output_dir, spec$rds))
  
  p_reg    <- predict(rf, data = reg_pixels_socio[,    ..socio_vars])$predictions[, "1"]
  p_nonreg <- predict(rf, data = pred_sample_socio[, ..socio_vars])$predictions[, "1"]
  rm(rf); gc()
  
  sens     <- mean(p_reg    > 0.5)
  spec_val <- mean(p_nonreg <= 0.5)
  
  cat(sprintf("%-30s  %11.4f  %11.4f  %9d  %11d\n",
              spec$lbl, sens, spec_val,
              length(p_reg), length(p_nonreg)))
}

rm(reg_pixels_socio, pred_sample_socio); gc()
# ============================================================
# S7: OOB Metrics & Area Estimates (Baseline & Paired)
# ============================================================
# Computes OOB metrics for both socioeconomic models and recalibrates
# predictions to a 5% landscape prior using the Dal Pozzolo (2015) formula.
message("\n--- S7: OOB Metrics & Area Estimates ---")

# ---- Parameters and Load Data ----
pi_train <- 0.50
pi_new   <- 0.05
beta     <- (pi_new / (1 - pi_new)) / (pi_train / (1 - pi_train)) # ~ 1/19

total_area_ha           <- readRDS(file.path(output_dir, "total_domain_area_ha.rds"))
regrowth_domain_area_ha <- 933250.8113584173   # from GEE mask (same as 02_figures.R)

pred_sample_raw   <- readRDS(file.path(output_dir, "pred_sample.rds"))
pred_sample_socio <- dt[pred_sample_raw[, .(point_id)], on = "point_id", nomatch = NULL]
pred_sample_socio <- pred_sample_socio[complete.cases(pred_sample_socio[, ..socio_vars])]

# Regrowth pixels for total-domain area estimates
reg_pixels_socio_s7 <- dt_socio[sample_set == "regrowth"]

auc_roc <- function(y_true_01, prob_pos) {
  pos <- prob_pos[y_true_01 == 1L]
  neg <- prob_pos[y_true_01 == 0L]
  n1  <- as.numeric(length(pos)); n0 <- as.numeric(length(neg))
  r   <- rank(c(pos, neg))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

recalibrate_dalpozzolo <- function(p, beta) {
  (beta * p) / ((beta - 1) * p + 1)
}

# ---- Model 1: Socio 50/50 Baseline ----
message("  Processing Socio 50/50 Baseline...")
rf_50 <- readRDS(file.path(output_dir, "rf_socio_50_50.rds"))

set.seed(42)
n_pos_s <- min(round(n_total * 0.50), dt_socio[sample_set == "regrowth",    .N])
n_neg_s <- min(n_total - round(n_total * 0.50), dt_socio[sample_set == "nonregrowth", .N])

train_pos_s <- dt_socio[sample_set == "regrowth"]    |> _[sample(.N, n_pos_s)]
train_neg_s <- dt_socio[sample_set == "nonregrowth"] |> _[sample(.N, n_neg_s)]
train_s     <- rbind(train_pos_s, train_neg_s)
y_train_50  <- as.integer(train_s$y == "1")

# OOB Metrics
oob_probs_50 <- rf_50$predictions[, "1"]
acc_50       <- mean((oob_probs_50 > 0.5) == y_train_50)
auc_50       <- auc_roc(y_train_50, oob_probs_50)
brier_50     <- rf_50$prediction.error

# Full Sample Metrics & Area
raw_probs_50 <- predict(rf_50, data = pred_sample_socio[, ..socio_vars])$predictions[, "1"]
adj_probs_50 <- recalibrate_dalpozzolo(raw_probs_50, beta)

reg_probs_50 <- predict(rf_50, data = dt_socio[sample_set == "regrowth", ..socio_vars])$predictions[, "1"]
sens_50      <- mean(reg_probs_50 > 0.5)
spec_50      <- mean(raw_probs_50 <= 0.5)

area_raw_50  <- mean(raw_probs_50) * total_area_ha / 1e6 +
  mean(reg_probs_50) * regrowth_domain_area_ha / 1e6
adj_reg_50   <- recalibrate_dalpozzolo(reg_probs_50, beta)
area_adj_50  <- mean(adj_probs_50) * total_area_ha / 1e6 +
  mean(adj_reg_50) * regrowth_domain_area_ha / 1e6

rm(rf_50, train_pos_s, train_neg_s, train_s); gc()


# ---- Model 2: Socio Paired (50/40/10) ----
message("  Processing Socio Paired (50/40/10)...")
rf_paired <- readRDS(file.path(output_dir, "rf_socio_paired.rds"))

negatives_paired_raw <- readRDS(file.path(output_dir, "negatives_paired_p1.rds"))
matched_ids          <- negatives_paired_raw[, .(point_id)]
negatives_paired     <- dt[matched_ids, on = "point_id", nomatch = NULL]
negatives_paired_socio <- negatives_paired[complete.cases(negatives_paired[, ..socio_vars])]

set.seed(42)
n_pos_sp  <- min(round(n_total * 0.50), dt_socio[sample_set == "regrowth",    .N])
n_pair_sp <- min(round(n_total * 0.40), nrow(negatives_paired_socio))
n_wide_sp <- min(n_total - round(n_total * 0.50) - round(n_total * 0.40),
                 dt_socio[sample_set == "nonregrowth", .N])

train_pos_sp  <- dt_socio[sample_set == "regrowth"]    |> _[sample(.N, n_pos_sp)]
train_pair_sp <- negatives_paired_socio[sample(.N, n_pair_sp)]
train_wide_sp <- dt_socio[sample_set == "nonregrowth"] |> _[sample(.N, n_wide_sp)]
train_sp      <- rbind(train_pos_sp, train_pair_sp, train_wide_sp)
y_train_paired <- as.integer(train_sp$y == "1")

# OOB Metrics
oob_probs_paired <- rf_paired$predictions[, "1"]
acc_paired       <- mean((oob_probs_paired > 0.5) == y_train_paired)
auc_paired       <- auc_roc(y_train_paired, oob_probs_paired)
brier_paired     <- rf_paired$prediction.error

# Full Sample Metrics & Area
raw_probs_paired <- predict(rf_paired, data = pred_sample_socio[, ..socio_vars])$predictions[, "1"]
adj_probs_paired <- recalibrate_dalpozzolo(raw_probs_paired, beta)

reg_probs_paired <- predict(rf_paired, data = dt_socio[sample_set == "regrowth", ..socio_vars])$predictions[, "1"]
sens_paired      <- mean(reg_probs_paired > 0.5)
spec_paired      <- mean(raw_probs_paired <= 0.5)

area_raw_paired  <- mean(raw_probs_paired) * total_area_ha / 1e6 +
  mean(reg_probs_paired) * regrowth_domain_area_ha / 1e6
adj_reg_paired   <- recalibrate_dalpozzolo(reg_probs_paired, beta)
area_adj_paired  <- mean(adj_probs_paired) * total_area_ha / 1e6 +
  mean(adj_reg_paired) * regrowth_domain_area_ha / 1e6

rm(rf_paired, negatives_paired_raw, matched_ids, negatives_paired, negatives_paired_socio, train_pos_sp, train_pair_sp, train_wide_sp, train_sp); gc()


# ---- Print to Console ----
cat("\n=== SOCIOECONOMIC MODELS: Performance & Area ===\n")
cat(sprintf("%-28s  %4s  %4s  %5s  %4s  %4s  %6s  %6s\n", 
            "Model", "Acc", "AUC", "Brier", "Sens", "Spec", "Raw(M)", "Adj(M)"))
cat(strrep("-", 76), "\n")
cat(sprintf("%-28s  %.3f  %.3f  %.3f  %.3f  %.3f  %6.1f  %6.1f\n", 
            "Socio 50/50 Baseline", acc_50, auc_50, brier_50, sens_50, spec_50, area_raw_50, area_adj_50))
cat(sprintf("%-28s  %.3f  %.3f  %.3f  %.3f  %.3f  %6.1f  %6.1f\n", 
            "Socio Paired (50/40/10)", acc_paired, auc_paired, brier_paired, sens_paired, spec_paired, area_raw_paired, area_adj_paired))
cat(strrep("-", 76), "\n")


# ---- Build & Save LaTeX Table ----
tex_socio_table <- sprintf(
  "\\begin{tabular}{@{}lccccccc@{}}
\\toprule
& \\multicolumn{3}{c}{OOB (training set)} & \\multicolumn{2}{c}{Full sample} & \\multicolumn{2}{c}{Area (Mha)} \\\\
\\cmidrule(lr){2-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}
Model & Accuracy & AUC & Brier & Sensitivity & Specificity & Raw & Adj (5\\%%) \\\\
\\midrule
Socio 50/50 baseline          & %.3f & %.3f & %.3f & %.3f & %.3f & %.1f & %.1f \\\\
Socio paired (50/40/10)       & %.3f & %.3f & %.3f & %.3f & %.3f & %.1f & %.1f \\\\
\\bottomrule
\\end{tabular}",
  acc_50, auc_50, brier_50, sens_50, spec_50, area_raw_50, area_adj_50,
  acc_paired, auc_paired, brier_paired, sens_paired, spec_paired, area_raw_paired, area_adj_paired
)

writeLines(tex_socio_table, file.path(output_dir, "table_performance_socio.tex"))
message("Saved: table_performance_socio.tex")


# Save the adjusted probabilities for future mapping/comparisons
saveRDS(list(
  baseline_raw = raw_probs_50, baseline_adj = adj_probs_50,
  paired_raw = raw_probs_paired, paired_adj = adj_probs_paired
), file.path(output_dir, "socio_recalibrated_probs.rds"))

# ============================================================
# S8: Socioeconomic Counterfactual (Minimal Suppression)
# ============================================================
message("\n--- S8: Socioeconomic Counterfactual Analysis ---")

# 1. Isolate purely continuous socioeconomic variables 
# (Excludes categorical landcover, weights, and biophysical baselines)
pure_socio_vars <- setdiff(socio_vars, c(extended_vars, "cropland_weight", "bio_landcover_class"))

# 2. Build background sample for continuous sweeps
set.seed(789)
cf_bg <- dt[complete.cases(dt[, ..socio_vars]) & sample_set == "nonregrowth"] |> 
  _[sample(.N, min(.N, 3000))]

rf_paired <- readRDS(file.path(output_dir, "rf_socio_paired.rds"))

# 3. Sweep continuous socio variables for optimal (max prob) values (BATCHED)
message("  Sweeping continuous socio variables (Batched for speed)...")

best_values <- rbindlist(lapply(pure_socio_vars, function(v) {
  
  # Define 35 grid points from 2nd to 98th percentile
  grid_vals <- seq(quantile(cf_bg[[v]], 0.02, na.rm = TRUE), 
                   quantile(cf_bg[[v]], 0.98, na.rm = TRUE), length.out = 35)
  
  # Create one large batched data.table (105,000 rows)
  dt_batch <- cf_bg[rep(1:.N, each = length(grid_vals))]
  dt_batch[, (v) := rep(grid_vals, times = nrow(cf_bg))]
  
  # Predict ONCE for the whole batch (maxes out CPU cores)
  dt_batch[, pred := predict(rf_paired, data = dt_batch[, ..socio_vars])$predictions[, "1"]]
  
  # Average predictions for each grid value to find the peak
  mean_probs <- dt_batch[, .(prob = mean(pred)), by = v]
  data.table(variable = v, optimal_val = mean_probs[[v]][which.max(mean_probs$prob)])
}))

# 4. Create counterfactual dataset
dt_counter <- copy(pred_sample_socio)

# 4a. Apply optimal continuous values (Removes infrastructure/population friction)
for (v in best_values$variable) {
  dt_counter[, (v) := best_values[variable == v, optimal_val]]
}

# 4b. Apply land-use transition (Data-driven abandonment)
message("  Simulating optimal natural land abandonment (Pixel-by-pixel)...")

if ("bio_landcover_class" %in% names(dt_counter)) {
  # Identify indices of active agricultural pixels (Pure Crop = 1, Mosaic = 3)
  ag_idx <- dt_counter[bio_landcover_class %in% c("1", "3"), which = TRUE]
  
  if (length(ag_idx) > 0) {
    message(sprintf("    Testing optimal natural state for %d agricultural pixels...", length(ag_idx)))
    
    # Isolate the ag pixels
    dt_ag <- dt_counter[ag_idx]
    
    # We will test Shrubland (8) and Grassland (9)
    candidate_classes <- c("8", "9")
    
    # Create a batch: duplicate the ag pixels for each candidate class
    dt_ag_batch <- dt_ag[rep(1:.N, each = length(candidate_classes))]
    
    # Assign the candidate classes and fix factor levels
    dt_ag_batch[, bio_landcover_class := factor(rep(candidate_classes, times = nrow(dt_ag)), 
                                                levels = levels(dt$bio_landcover_class))]
    
    # Zero out the physical footprint weight for these abandoned states
    if ("cropland_weight" %in% names(dt_ag_batch)) {
      dt_ag_batch[, cropland_weight := 0]
    }
    
    # Predict regrowth probability for all candidates simultaneously
    dt_ag_batch[, pred_prob := predict(rf_paired, data = dt_ag_batch[, ..socio_vars])$predictions[, "1"]]
    
    # Group by the original pixel row, and keep the class that gave the HIGHEST probability
    dt_ag_batch[, pixel_id := rep(1:nrow(dt_ag), each = length(candidate_classes))]
    best_natural_states <- dt_ag_batch[, .SD[which.max(pred_prob)], by = pixel_id]
    
    # Push the optimal classes back into the main counterfactual dataset
    dt_counter[ag_idx, bio_landcover_class := best_natural_states$bio_landcover_class]
    if ("cropland_weight" %in% names(dt_counter)) {
      dt_counter[ag_idx, cropland_weight := 0]
    }
    
    # Brief summary of what the model preferred
    n_shrub <- sum(best_natural_states$bio_landcover_class == "8")
    n_grass <- sum(best_natural_states$bio_landcover_class == "9")
    message(sprintf("    Model preference: %d Shrubland, %d Grassland", n_shrub, n_grass))
  }
}

# Build counterfactual for regrowth pixels too
dt_counter_reg <- copy(reg_pixels_socio_s7)

# Apply same optimal continuous values
for (v in best_values$variable) {
  dt_counter_reg[, (v) := best_values[variable == v, optimal_val]]
}

# Apply same land-cover transition for any ag pixels in regrowth domain
if ("bio_landcover_class" %in% names(dt_counter_reg)) {
  ag_idx_reg <- dt_counter_reg[bio_landcover_class %in% c("1", "3"), which = TRUE]
  if (length(ag_idx_reg) > 0) {
    dt_ag_reg <- dt_counter_reg[ag_idx_reg]
    dt_ag_batch_reg <- dt_ag_reg[rep(1:.N, each = 2)]
    dt_ag_batch_reg[, bio_landcover_class := factor(rep(c("8", "9"), times = nrow(dt_ag_reg)),
                                                    levels = levels(dt$bio_landcover_class))]
    if ("cropland_weight" %in% names(dt_ag_batch_reg)) {
      dt_ag_batch_reg[, cropland_weight := 0]
    }
    dt_ag_batch_reg[, pred_prob := predict(rf_paired, data = dt_ag_batch_reg[, ..socio_vars])$predictions[, "1"]]
    dt_ag_batch_reg[, pixel_id := rep(1:nrow(dt_ag_reg), each = 2)]
    best_reg <- dt_ag_batch_reg[, .SD[which.max(pred_prob)], by = pixel_id]
    dt_counter_reg[ag_idx_reg, bio_landcover_class := best_reg$bio_landcover_class]
    if ("cropland_weight" %in% names(dt_counter_reg)) {
      dt_counter_reg[ag_idx_reg, cropland_weight := 0]
    }
  }
}

# Predict counterfactual on regrowth domain
raw_probs_cf_reg <- predict(rf_paired, data = dt_counter_reg[, ..socio_vars])$predictions[, "1"]
adj_probs_cf_reg <- recalibrate_dalpozzolo(raw_probs_cf_reg, beta)

# Total counterfactual area (full domain)
area_adj_cf <- mean(adj_probs_cf) * total_area_ha / 1e6 +
  mean(adj_probs_cf_reg) * regrowth_domain_area_ha / 1e6
uplift_mha  <- area_adj_cf - area_adj_paired


# 5. Predict & Recalibrate
message("  Predicting counterfactual landscape...")
raw_probs_cf <- predict(rf_paired, data = dt_counter[, ..socio_vars])$predictions[, "1"]
adj_probs_cf <- recalibrate_dalpozzolo(raw_probs_cf, beta)

# Calculate areas using total domain (nonreg + reg)
# Counterfactual only modifies non-regrowth pixels; regrowth component is unchanged
area_adj_cf <- mean(adj_probs_cf) * total_area_ha / 1e6 +
  mean(adj_reg_paired) * regrowth_domain_area_ha / 1e6
uplift_mha  <- area_adj_cf - area_adj_paired

# 6. Output Results
cat("\n=== COUNTERFACTUAL RESULTS ===\n")
cat(sprintf("Observed Socio Area (Adj):  %6.1f Mha\n", area_adj_paired))
cat(sprintf("Potential Socio Area (Adj): %6.1f Mha\n", area_adj_cf))
cat(sprintf("Estimated Socio Friction:   %6.1f Mha\n", uplift_mha))

saveRDS(list(cf_raw = raw_probs_cf, cf_adj = adj_probs_cf), 
        file.path(output_dir, "socio_counterfactual_probs.rds"))

# Clean up memory
rm(cf_bg, dt_counter, rf_paired); gc()
message("S8 complete.")


# Parameters from your environment
pi_train <- 0.50
pi_new   <- 0.05
beta     <- (pi_new / (1 - pi_new)) / (pi_train / (1 - pi_train))
total_area_ha           <- readRDS(file.path(output_dir, "total_domain_area_ha.rds"))
regrowth_domain_area_ha <- 933250.8113584173   # from GEE mask

recalibrate_dalpozzolo <- function(p, beta) {
  (beta * p) / ((beta - 1) * p + 1)
}

# Load the extended paired predictions (generated in 02_figures.R)
preds_ext <- readRDS(file.path(output_dir, "preds_rf_ext_paired.rds"))

# Recalibrate and calculate area (total domain = nonreg + reg)
adj_probs_ext     <- recalibrate_dalpozzolo(preds_ext$nonreg, beta)
adj_probs_ext_reg <- recalibrate_dalpozzolo(preds_ext$reg, beta)
area_ext_adj      <- mean(adj_probs_ext) * total_area_ha / 1e6 +
  mean(adj_probs_ext_reg) * regrowth_domain_area_ha / 1e6

cat(sprintf("Extended Biophysical Area (Adj): %6.1f Mha\n", area_ext_adj))


message("\n04_socioeconomic_models.R complete.")