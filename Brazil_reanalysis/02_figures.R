# ============================================================
# 02_figures.R
# Brazil Natural Regeneration — Williams et al. Replication
#
# Loads saved RDS files from 01_train_models.R and produces:
#   fig_main.pdf/.png  — three-panel figure for Matters Arising
#     Panel B (top-left):  Variable importance shift (baseline vs CE)
#     Panel A (top-right): Area estimate vs prevalence assumption
#     Panel C (bottom):    Partial dependence on bioclim PC1
#
# LaTeX caption for fig_main:
#   (a) Estimated regeneration area as mean predicted probability
#   over 797,076 random non-regrowth pixels multiplied by the total
#   non-regrowth domain (169.9 Mha). Training sample size fixed at
#   398,538 observations across all prevalence ratios; only the
#   positive:negative ratio varies.
#   (b) Variable permutation importance expressed as percentage of
#   total importance under the 50/50 baseline model (Williams et al.
#   design) and a climate-envelope-controlled (CE) model trained with
#   locally matched negatives. Red: between-climate classifiers
#   (Bioclim PC1, capturing the broad thermal gradient; Biome,
#   discrete climate zone; Soil pH, biome-correlated soil chemistry).
#   Blue: local within-climate variables (Forest density, Distance to
#   forest, Soil organic carbon).
#   (c) Partial dependence of predicted P(regrowth) on Bioclim PC1
#   (the dominant thermal gradient separating tropical from subtropical
#   Brazil), averaged over 5,000 random non-regrowth pixels. CE model
#   trained with greedy nearest-neighbour matched negatives (k=1000;
#   83.7% of regrowth points matched; median match distance 3,500 m).
#
# LaTeX caption for fig_supp_pdp:
#   Partial dependence plots for four key predictors under the 50/50
#   baseline model (red, solid) and the CE model (blue, dashed),
#   averaged over 2,000 random non-regrowth pixels. Distance to
#   existing forest and surrounding forest density show near-identical
#   functional relationships under both training designs, confirming
#   these local landscape signals are robust to the negative sampling
#   strategy. Bioclim PC2 captures dry-season severity within Brazil
#   and shows a modest positive relationship with P(regrowth) under
#   both models.
# ============================================================

library(data.table)
library(ranger)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(kableExtra)


setwd("/home/cannon/Dropbox/Regrowth/pnv_sampling/brazil_predict")

output_dir <- "rf_outputs"


# ---- Shared supporting objects (no ranger objects yet) ------

williams_vars        <- readRDS(file.path(output_dir, "williams_vars.rds"))
extended_vars        <- readRDS(file.path(output_dir, "extended_vars.rds"))
pred_sample          <- readRDS(file.path(output_dir, "pred_sample.rds"))
total_domain_area_ha <- readRDS(file.path(output_dir, "total_domain_area_ha.rds"))
regrowth_domain_area_ha <- 933250.8113584173   # from GEE mask

# ---- Build prediction surfaces (one read of the CSV) --------
# Load the full training data once to build:
#   (a) reg_pixels     — regrowth rows, complete for Williams vars
#   (b) reg_pixels_ext — regrowth rows, complete for extended vars
#   (c) pred_bg        — 10k nonregrowth sample for PDPs
# All ranger objects are loaded later, one at a time.

message("Loading training data and projecting PCA...")
dt_full <- fread("gee_outputs/brazil_training_data.csv")
dt_full[, bio_landcover_class := as.factor(bio_landcover_class)]
dt_full[, bio_biome_id        := as.factor(bio_biome_id)]

pca_fit  <- readRDS(file.path(output_dir, "pca_fit.rds"))
bio_cols <- paste0("bio", sprintf("%02d", 1:19))
pcs      <- predict(pca_fit, dt_full[, ..bio_cols])
dt_full[, bioclim_pc1 := pcs[, 1]]
dt_full[, bioclim_pc2 := pcs[, 2]]
dt_full[, bioclim_pc3 := pcs[, 3]]
dt_full[, bioclim_pc4 := pcs[, 4]]
dt_full[, bioclim_pc5 := pcs[, 5]]
rm(pcs); gc()

# Regrowth pixel subsets
reg_pixels     <- dt_full[sample_set == "regrowth"]
reg_pixels     <- reg_pixels[complete.cases(reg_pixels[, ..williams_vars])]
reg_pixels_ext <- dt_full[sample_set == "regrowth"]
reg_pixels_ext <- reg_pixels_ext[complete.cases(reg_pixels_ext[, ..extended_vars])]

# PDP background sample — nonregrowth, Williams vars, built here so we
# don't need to reload the CSV again later
set.seed(456)
pred_bg <- dt_full[sample_set == "nonregrowth"][sample(.N, 10000)]
pred_bg <- pred_bg[complete.cases(pred_bg[, ..williams_vars])]

rm(dt_full); gc()
message(sprintf("Regrowth pixels (Williams vars): %d  |  (extended vars): %d",
                nrow(reg_pixels), nrow(reg_pixels_ext)))

# ============================================================
# PREDICTION LOOP
# Each model is described by: display label, RDS filename, which
# variable set to use, and a per-model cache filename.
#
# For each model the loop:
#   1. Checks whether the cache file already exists — skips if so.
#   2. Loads the ranger object.
#   3. Predicts on pred_sample (non-regrowth) AND the appropriate
#      regrowth pixel table.
#   4. Saves both vectors to the per-model cache and drops the
#      ranger object immediately.
#
# On re-runs only missing caches trigger a ranger load, so
# partially-completed runs resume cheaply.
# ============================================================

model_specs <- list(
  list(label = "50/50",
       rds   = "rf_prev_50_50.rds",
       vars  = "williams",
       cache = "preds_rf_prev_50_50.rds"),
  list(label = "20/80",
       rds   = "rf_prev_20_80.rds",
       vars  = "williams",
       cache = "preds_rf_prev_20_80.rds"),
  list(label = "5/95",
       rds   = "rf_prev_5_95.rds",
       vars  = "williams",
       cache = "preds_rf_prev_5_95.rds"),
  list(label = "Paired",
       rds   = "rf_ce_50_40_10_nn.rds",
       vars  = "williams",
       cache = "preds_rf_ce.rds"),
  list(label = "ext_50/50",
       rds   = "rf_ext_50_50.rds",
       vars  = "extended",
       cache = "preds_rf_ext_50_50.rds"),
  list(label = "ext_Paired",
       rds   = "rf_ext_paired.rds",
       vars  = "extended",
       cache = "preds_rf_ext_paired.rds"),
  # Paired prevalence sensitivity (R7 models)
  list(label = "paired_20/80",
       rds   = "rf_paired_20_80.rds",
       vars  = "williams",
       cache = "preds_rf_paired_20_80.rds"),
  list(label = "paired_5/95",
       rds   = "rf_paired_5_95.rds",
       vars  = "williams",
       cache = "preds_rf_paired_5_95.rds"),
  list(label = "ext_paired_20/80",
       rds   = "rf_ext_paired_20_80.rds",
       vars  = "extended",
       cache = "preds_rf_ext_paired_20_80.rds"),
  list(label = "ext_paired_5/95",
       rds   = "rf_ext_paired_5_95.rds",
       vars  = "extended",
       cache = "preds_rf_ext_paired_5_95.rds")
)

pred_results <- list()

for (spec in model_specs) {
  
  cache_path <- file.path(output_dir, spec$cache)
  
  if (file.exists(cache_path)) {
    message(sprintf("  [cache] %s", spec$label))
    pred_results[[spec$label]] <- readRDS(cache_path)
    next
  }
  
  message(sprintf("  [run]   %s — loading %s", spec$label, spec$rds))
  rf_tmp   <- readRDS(file.path(output_dir, spec$rds))
  use_vars <- if (spec$vars == "williams") williams_vars else extended_vars
  reg_src  <- if (spec$vars == "williams") reg_pixels    else reg_pixels_ext
  
  p_nonreg <- predict(rf_tmp, data = pred_sample[, ..use_vars])$predictions[, "1"]
  p_reg    <- predict(rf_tmp, data = reg_src[,    ..use_vars])$predictions[, "1"]
  
  rm(rf_tmp); gc()
  
  result <- list(nonreg = p_nonreg, reg = p_reg)
  saveRDS(result, cache_path)
  pred_results[[spec$label]] <- result
  message(sprintf("            cached → %s", spec$cache))
}

# Keep slim coordinate tables for the spatial map panel — just lat/lon
# and point_id so we can join reg predictions back to space later.
# Only need Williams-var regrowth pixels (same set used for p_reg in
# the 50/50 and Paired predictions that feed into pD).
reg_coords <- reg_pixels[, .(point_id, latitude, longitude)]


# Unpack into named vectors for all downstream code
preds_50_50               <- pred_results[["50/50"]]$nonreg
preds_20_80               <- pred_results[["20/80"]]$nonreg
preds_5_95                <- pred_results[["5/95"]]$nonreg
preds_ce                  <- pred_results[["Paired"]]$nonreg
preds_ext_50_50           <- pred_results[["ext_50/50"]]$nonreg
preds_ext_paired          <- pred_results[["ext_Paired"]]$nonreg
preds_paired_20_80        <- pred_results[["paired_20/80"]]$nonreg
preds_paired_5_95         <- pred_results[["paired_5/95"]]$nonreg
preds_ext_paired_20_80    <- pred_results[["ext_paired_20/80"]]$nonreg
preds_ext_paired_5_95     <- pred_results[["ext_paired_5/95"]]$nonreg
# Non-regrowth extended (used in diagnostics)
nonreg_preds_ext_50_50    <- pred_results[["ext_50/50"]]$nonreg
nonreg_preds_ext_paired   <- pred_results[["ext_Paired"]]$nonreg
# Regrowth predictions
reg_preds_50_50               <- pred_results[["50/50"]]$reg
reg_preds_20_80               <- pred_results[["20/80"]]$reg
reg_preds_5_95                <- pred_results[["5/95"]]$reg
reg_preds_ce                  <- pred_results[["Paired"]]$reg
reg_preds_ext_50_50           <- pred_results[["ext_50/50"]]$reg
reg_preds_ext_paired          <- pred_results[["ext_Paired"]]$reg
reg_preds_paired_20_80        <- pred_results[["paired_20/80"]]$reg
reg_preds_paired_5_95         <- pred_results[["paired_5/95"]]$reg
reg_preds_ext_paired_20_80    <- pred_results[["ext_paired_20/80"]]$reg
reg_preds_ext_paired_5_95     <- pred_results[["ext_paired_5/95"]]$reg

rm(pred_results); gc()

# Attach to pred_sample for the spatial map panel
pred_sample[, p_baseline := preds_50_50]
pred_sample[, p_paired   := preds_ce]

# ---- Area estimates ----------------------------------------

area_dt <- data.table(
  model    = c("50/50", "20/80", "5/95",
               "Paired", "Paired_20/80", "Paired_5/95",
               "50/50_ext", "Paired_ext",
               "Paired_ext_20/80", "Paired_ext_5/95"),
  design   = c("Baseline random",  "Baseline random",  "Baseline random",
               "Williams paired",  "Williams paired",  "Williams paired",
               "Extended baseline","Extended paired",
               "Extended paired",  "Extended paired"),
  prev_pos = c(0.50, 0.20, 0.05,
               0.50, 0.20, 0.05,
               0.50, 0.50,
               0.20, 0.05),
  area_nonreg_Mha = c(
    mean(preds_50_50)            * total_domain_area_ha / 1e6,
    mean(preds_20_80)            * total_domain_area_ha / 1e6,
    mean(preds_5_95)             * total_domain_area_ha / 1e6,
    mean(preds_ce)               * total_domain_area_ha / 1e6,
    mean(preds_paired_20_80)     * total_domain_area_ha / 1e6,
    mean(preds_paired_5_95)      * total_domain_area_ha / 1e6,
    mean(preds_ext_50_50)        * total_domain_area_ha / 1e6,
    mean(preds_ext_paired)       * total_domain_area_ha / 1e6,
    mean(preds_ext_paired_20_80) * total_domain_area_ha / 1e6,
    mean(preds_ext_paired_5_95)  * total_domain_area_ha / 1e6
  ),
  area_reg_Mha = c(
    mean(reg_preds_50_50)            * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_20_80)            * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_5_95)             * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_ce)               * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_paired_20_80)     * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_paired_5_95)      * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_ext_50_50)        * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_ext_paired)       * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_ext_paired_20_80) * regrowth_domain_area_ha / 1e6,
    mean(reg_preds_ext_paired_5_95)  * regrowth_domain_area_ha / 1e6
  )
)
area_dt[, area_total_Mha := area_nonreg_Mha + area_reg_Mha]

cat("\nArea estimates (Mha):\n")
print(area_dt)

# ---- Variable importance (relative share) ------------------
# Load each ranger once, extract importance table, drop immediately.

compute_rel_imp <- function(rf, model_label) {
  vi    <- rf$variable.importance
  total <- sum(vi)
  data.table(
    model    = model_label,
    variable = names(vi),
    rel_imp  = 100 * vi / total
  )
}

message("Extracting variable importance...")
rf_tmp      <- readRDS(file.path(output_dir, "rf_prev_50_50.rds"))
vi_baseline <- compute_rel_imp(rf_tmp, "50/50 baseline")
rm(rf_tmp); gc()

rf_tmp   <- readRDS(file.path(output_dir, "rf_ce_50_40_10_nn.rds"))
vi_ce_dt <- compute_rel_imp(rf_tmp, "Paired")
rm(rf_tmp); gc()

# Short display names
clean_names <- c(
  bio_forest_density_1km2  = "Forest density",
  bioclim_pc1              = "Bioclim PC1",
  bio_dist_forest_2000     = "Dist. to forest",
  bio_soil_ph              = "Soil pH",
  bioclim_pc2              = "Bioclim PC2",
  bioclim_pc3              = "Bioclim PC3",
  bio_landcover_class      = "Landcover",
  bioclim_pc4              = "Bioclim PC4",
  bio_soil_organic_carbon  = "Soil org. C",
  bio_biome_id             = "Biome"
)

vi_both <- rbind(vi_baseline, vi_ce_dt)
vi_both[, var_label := clean_names[variable]]
vi_both[, model     := factor(model, levels = c("50/50 baseline", "Paired"))]

var_order <- vi_baseline[order(-rel_imp), clean_names[variable]]
vi_both[, var_label := factor(var_label, levels = rev(var_order))]

# Direction flags for three-colour scheme:
#   red  (down) = between-climate classifiers: PC1, Biome, Soil pH
#   blue (up)   = local within-climate variables: Forest density,
#                 Dist. to forest, Soil org. C
#   grey        = other (PC2, PC3, PC4, Landcover)
vi_both[, direction := fifelse(
  var_label %in% c("Bioclim PC1", "Biome", "Soil pH"), "down",
  fifelse(var_label %in% c("Forest density", "Dist. to forest", "Soil org. C"),
          "up", "neutral")
)]

cat("\nPC1 relative share:\n")
cat(sprintf("  50/50 baseline: %.1f%%\n",
            vi_baseline[variable == "bioclim_pc1", rel_imp]))
cat(sprintf("  Paired:         %.1f%%\n",
            vi_ce_dt[variable == "bioclim_pc1", rel_imp]))

# ---- Partial dependence (PC1 only) -------------------------
# pred_bg was built above during the single CSV load and is still in memory.
# Load each ranger, compute PDP, drop immediately.

message("Computing partial dependence for bioclim_pc1...")

compute_pdp <- function(rf, data, variable, use_vars, n_grid = 50, n_sample = 10000) {
  set.seed(123)
  data_s    <- data[sample(min(.N, n_sample))]
  grid_vals <- seq(quantile(data_s[[variable]], 0.02),
                   quantile(data_s[[variable]], 0.98),
                   length.out = n_grid)
  pdp_vals  <- sapply(grid_vals, \(val) {
    data_tmp <- copy(data_s)
    data_tmp[, (variable) := val]
    mean(predict(rf, data = data_tmp[, ..use_vars])$predictions[, "1"])
  })
  data.table(x = grid_vals, prob = pdp_vals)
}

rf_tmp    <- readRDS(file.path(output_dir, "rf_prev_50_50.rds"))
pdp_50_50 <- compute_pdp(rf_tmp, pred_bg, "bioclim_pc1", williams_vars)
pdp_50_50[, model := "50/50 baseline"]
rm(rf_tmp); gc()

rf_tmp     <- readRDS(file.path(output_dir, "rf_ce_50_40_10_nn.rds"))
pdp_ce_pc1 <- compute_pdp(rf_tmp, pred_bg, "bioclim_pc1", williams_vars)
pdp_ce_pc1[, model := "Paired"]
rm(rf_tmp); gc()

pdp_pc1 <- rbind(pdp_50_50, pdp_ce_pc1)
rm(pred_bg); gc()

# ============================================================
# FIGURE ASSEMBLY
# Layout: (pB | pA) / pC
#   pB top-left  — variable importance slope (wider, needs room for labels)
#   pA top-right — area vs prevalence (compact)
#   pC bottom    — PC1 partial dependence (full width)
# ============================================================

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

# ---- Panel A: Area vs prevalence ---------------------------
# Three series: total (non-reg + reg), non-regrowth only, regrowth only.
# Regrowth domain is ~0.93 Mha so its absolute contribution is small
# but shown explicitly to flag the omission in Williams et al.

area_prev_dt <- area_dt[!is.na(prev_pos)]
area_prev_dt <- area_dt[design == "Baseline random"]
area_prev_dt[, prev_label := c("50/50\n(Williams)", "20/80", "5/95\n(landscape)")]
area_prev_dt[, prev_label := factor(prev_label,
                                    levels = c("50/50\n(Williams)",
                                               "20/80",
                                               "5/95\n(landscape)"))]

area_long <- melt(area_prev_dt,
                  id.vars       = c("model", "prev_pos", "prev_label"),
                  measure.vars  = c("area_total_Mha", "area_nonreg_Mha", "area_reg_Mha"),
                  variable.name = "series",
                  value.name    = "area_Mha")

area_long[, series_label := factor(series,
  levels = c("area_total_Mha", "area_nonreg_Mha", "area_reg_Mha"),
  labels = c("Total (non-reg + reg)", "Non-regrowth only", "Regrowth only")
)]

series_colours <- c(
  "Total (non-reg + reg)" = "black",
  "Non-regrowth only"     = col_baseline,
  "Regrowth only"         = col_ce
)
series_ltys <- c(
  "Total (non-reg + reg)" = "solid",
  "Non-regrowth only"     = "dashed",
  "Regrowth only"         = "dotted"
)

pA <- ggplot(area_long, aes(x = prev_label, y = area_Mha,
                            colour = series_label, linetype = series_label,
                            group = series_label)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  annotate("text", x = 2.5, y = max(area_prev_dt$area_total_Mha) * 0.95,
           label = sprintf("%.1fx range",
                           max(area_prev_dt$area_total_Mha) /
                             min(area_prev_dt$area_total_Mha)),
           size = 3, hjust = 1, colour = "grey30") +
  scale_colour_manual(values = series_colours, name = NULL) +
  scale_linetype_manual(values = series_ltys,  name = NULL) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "a  Estimated area by training prevalence",
    x     = "Prevalence (regrowth / non-regrowth)",
    y     = "Estimated area (Mha)"
  ) +
  nature_theme +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

# ---- Panel B: Variable importance slope --------------------

vi_labels <- vi_both[model == "Paired"]

pB <- ggplot(vi_both, aes(x = model, y = rel_imp, group = var_label)) +
  geom_line(aes(colour = direction), linewidth = 0.7, alpha = 0.85) +
  geom_point(aes(colour = direction), size = 2) +
  scale_colour_manual(
    values = c("down"    = col_baseline,
               "up"      = col_ce,
               "neutral" = "grey70"),
    labels = c("down"    = "Climate classifiers (↓)",
               "up"      = "Local variables (↑)",
               "neutral" = "Other"),
    name   = NULL
  ) +
  geom_text_repel(
    data               = vi_labels,
    aes(label = var_label, colour = direction),
    hjust              = 0,
    direction          = "y",
    nudge_x            = 0.08,
    segment.size       = 0.3,
    segment.alpha      = 0.5,
    size               = 2.4,
    show.legend        = FALSE,
    force              = 3,
    min.segment.length = 0
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.45))) +
  labs(
    title = "c  Variable importance: baseline vs paired",
    x     = NULL,
    y     = "Relative importance (% of total)"
  ) +
  nature_theme +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

# ---- Panel C: PC1 partial dependence ----------------------

pC <- ggplot(pdp_pc1, aes(x = x, y = prob * 100,
                          colour = model, linetype = model)) +
  geom_line(linewidth = 1.0) +
  scale_colour_manual(
    values = c("50/50 baseline" = col_baseline, "Paired" = col_ce),
    name   = NULL
  ) +
  scale_linetype_manual(
    values = c("50/50 baseline" = "solid", "Paired" = "dashed"),
    name   = NULL
  ) +
  labs(
    title = "b  Marginal effect of Bioclim PC1 on P(regrowth)",
    x     = "Bioclim PC1 (cool/subtropical \u2190  \u2192 hot/tropical)",
    y     = "Mean predicted P(regrowth) (%)"
  ) +
  scale_y_continuous(limits = c(0, 35), # or c(0, 40) to be safe
                     expand = expansion(mult = c(0, 0))) +
  nature_theme +
  theme(
    legend.position = "bottom"
  ) 


# ---- Panel D: Spatial redistribution map -------------------
# Combines non-regrowth pred_sample with regrowth pixels, weighting
# each pixel by the area it represents (domain_area / n_pixels) so
# that the heavily oversampled regrowth domain doesn't dominate
# within-bin averages.

ha_per_nonreg <- total_domain_area_ha    / nrow(pred_sample)
ha_per_reg    <- regrowth_domain_area_ha / length(reg_preds_50_50)

nonreg_spatial <- data.table(
  latitude   = pred_sample$latitude,
  longitude  = pred_sample$longitude,
  p_baseline = preds_50_50,
  p_paired   = preds_ce,
  w          = ha_per_nonreg
)

reg_spatial <- data.table(
  latitude   = reg_coords$latitude,
  longitude  = reg_coords$longitude,
  p_baseline = reg_preds_50_50,
  p_paired   = reg_preds_ce,
  w          = ha_per_reg
)

spatial_all <- rbind(nonreg_spatial, reg_spatial)
rm(nonreg_spatial, reg_spatial); gc()

spatial_all[, lat_bin := round(latitude  * 15) / 15]
spatial_all[, lon_bin := round(longitude * 15) / 15]

grid_summary <- spatial_all[, .(
  delta = weighted.mean(p_paired - p_baseline, w),
  n     = .N
), by = .(lat_bin, lon_bin)]

# can increase number to drop bins with fewer pixels
grid_summary <- grid_summary[n >= 1]

pD <- ggplot(grid_summary, aes(x = lon_bin, y = lat_bin, fill = delta)) +
  geom_tile() +
  scale_fill_gradient2(
    low      = col_baseline,
    mid      = "#fffae6",
    high     = col_ce,
    midpoint = 0,
    name     = "Delta\n(paired \u2212 baseline)"
  ) +
  coord_fixed() +
  labs(
    title = "d  Spatial redistribution: paired - baseline",
    x     = "Longitude",
    y     = "Latitude"
  ) +
  nature_theme +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(1.0, "cm"))


# ---- Combine with patchwork --------------------------------
# Layout: pB full height left | right column (pA / pC / pD)

# Note the order: Row 1 (pA, pC), then Row 2 (pB, pD)
fig_main <- (pA + pC + pB + pD) + 
  plot_layout(
    ncol = 2, 
    # Col 1 (A, B) is narrow; Col 2 (C, D) is wide
    widths = c(1, 1.5),   
    # Row 1 (A, C) is short; Row 2 (B, D) is tall
    heights = c(0.7, 1.1)   
  )

ggsave(file.path(output_dir, "fig_main.pdf"),
       fig_main, device = cairo_pdf, width = 200, height = 180, units = "mm")

# ggsave(file.path(output_dir, "fig_main.png"),
#        fig_main, width = 180, height = 200, units = "mm", dpi = 300)
# message("Saved: fig_main.png")

print(fig_main)

# ---- Supplementary: full PDP four-panel --------------------
# ---- Build pdp_supp_full for supplementary PDP figure ------
# Load each ranger, compute PDPs for all 8 variables, drop ranger.
# pred_bg was dropped after the main PC1 PDP — reload a fresh sample here.

message("Building supplementary PDPs...")

pdp_vars_supp <- c(
  "bioclim_pc1", "bioclim_pc2", "bioclim_pc3", "bioclim_pc4",
  "bio_dist_forest_2000", "bio_soil_organic_carbon",
  "bio_soil_ph", "bio_forest_density_1km2"
)

# Fresh background sample (nonregrowth, 10k rows, Williams vars)
dt_tmp <- fread("gee_outputs/brazil_training_data.csv")
dt_tmp[, bio_landcover_class := as.factor(bio_landcover_class)]
dt_tmp[, bio_biome_id        := as.factor(bio_biome_id)]
pcs_tmp <- predict(pca_fit, dt_tmp[, ..bio_cols])
dt_tmp[, bioclim_pc1 := pcs_tmp[, 1]]
dt_tmp[, bioclim_pc2 := pcs_tmp[, 2]]
dt_tmp[, bioclim_pc3 := pcs_tmp[, 3]]
dt_tmp[, bioclim_pc4 := pcs_tmp[, 4]]
rm(pcs_tmp); gc()
set.seed(789)
pred_bg_supp <- dt_tmp[sample_set == "nonregrowth"][sample(.N, 2000)]
pred_bg_supp <- pred_bg_supp[complete.cases(pred_bg_supp[, ..williams_vars])]
rm(dt_tmp); gc()

compute_pdp_supp <- function(rf, data, use_vars, n_grid = 35, n_sample = 2000) {
  set.seed(123)
  data_s <- data[sample(min(.N, n_sample))]
  rbindlist(lapply(pdp_vars_supp, function(var) {
    grid_vals <- seq(quantile(data_s[[var]], 0.02),
                     quantile(data_s[[var]], 0.98),
                     length.out = n_grid)
    pdp_vals  <- sapply(grid_vals, \(val) {
      data_tmp <- copy(data_s)
      data_tmp[, (var) := val]
      mean(predict(rf, data = data_tmp[, ..use_vars])$predictions[, "1"])
    })
    data.table(variable = var, x = grid_vals, prob = pdp_vals)
  }))
}

message("  Supp PDP: loading rf_prev_50_50...")
rf_tmp         <- readRDS(file.path(output_dir, "rf_prev_50_50.rds"))
pdp_s_baseline <- compute_pdp_supp(rf_tmp, pred_bg_supp, williams_vars)
pdp_s_baseline[, model := "50/50 baseline"]
rm(rf_tmp); gc()

message("  Supp PDP: loading rf_ce...")
rf_tmp    <- readRDS(file.path(output_dir, "rf_ce_50_40_10_nn.rds"))
pdp_s_ce  <- compute_pdp_supp(rf_tmp, pred_bg_supp, williams_vars)
pdp_s_ce[, model := "Paired"]
rm(rf_tmp, pred_bg_supp); gc()

pdp_supp_full <- rbind(pdp_s_baseline, pdp_s_ce)
pdp_supp_full[, model := factor(model, levels = c("50/50 baseline", "Paired"))]


make_pdp_panel <- function(var, label, ylim_max) {
  ggplot(pdp_supp_full[variable == var],
         aes(x = x, y = prob * 100, colour = model, linetype = model)) +
    geom_line(linewidth = 1.0) +
    coord_cartesian(ylim = c(0, ylim_max)) +
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

panels <- list(
  make_pdp_panel("bioclim_pc1",             "Bioclim PC1 (thermal gradient)",    35),
  make_pdp_panel("bioclim_pc2",             "Bioclim PC2 (dry season severity)", 35),
  make_pdp_panel("bioclim_pc3",             "Bioclim PC3 (temp. stability)",     35),
  make_pdp_panel("bioclim_pc4",             "Bioclim PC4 (hot-wet season)",      35),
  make_pdp_panel("bio_dist_forest_2000",    "Distance to forest 2000 (m)",       NA),
  make_pdp_panel("bio_soil_organic_carbon", "Soil organic carbon",               35),
  make_pdp_panel("bio_soil_ph",             "Soil pH",                           35),
  make_pdp_panel("bio_forest_density_1km2", "Forest density 1km\u00b2 (%)",      55)
)

# Shared legend
legend_plot <- ggplot(pdp_supp_full[variable == "bioclim_pc1"],
                      aes(x = x, y = prob * 100, colour = model, linetype = model)) +
  geom_line() +
  scale_colour_manual(values = c("50/50 baseline" = col_baseline, "Paired" = col_ce)) +
  scale_linetype_manual(values = c("50/50 baseline" = "solid", "Paired" = "dashed")) +
  labs(colour = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

shared_legend <- cowplot::get_legend(legend_plot)

p_supp_full <- wrap_plots(panels, ncol = 3) /
  wrap_elements(shared_legend) +
  plot_layout(heights = c(20, 1)) +
  plot_annotation(
    caption = "x: variable value    y: mean predicted P(regrowth) (%)"
  )

ggsave(file.path(output_dir, "fig_supp_pdp.pdf"),
       p_supp_full, width = 180, height = 220, units = "mm")
ggsave(file.path(output_dir, "fig_supp_pdp.png"),
       p_supp_full, width = 180, height = 220, units = "mm", dpi = 300)
message("Saved: fig_supp_pdp.pdf / .png")


# ---- Panel B extended: variable importance (extended model) ----------------

# extended_vars already in memory from prediction loop
rf_tmp_ext    <- readRDS(file.path(output_dir, "rf_ext_50_50.rds"))
vi_ext_base   <- compute_rel_imp(rf_tmp_ext, "50/50 baseline")
rm(rf_tmp_ext); gc()

rf_tmp_ext    <- readRDS(file.path(output_dir, "rf_ext_paired.rds"))
vi_ext_paired <- compute_rel_imp(rf_tmp_ext, "Paired")
rm(rf_tmp_ext); gc()

clean_names_ext <- c(
  clean_names,  # reuse existing names from main figure
  bio_npp_mean               = "NPP",
  bio_fire_freq              = "Fire frequency",
  bio_slope                  = "Slope",
  bio_elevation              = "Elevation",
  bio_soil_sand              = "Soil sand",
  bio_soil_clay              = "Soil clay",
  bio_soil_bulk_density_fine = "Soil bulk density",
  bio_soil_water_33kpa       = "Soil water capacity"
)

vi_ext_both <- rbind(vi_ext_base, vi_ext_paired)
vi_ext_both[, var_label := clean_names_ext[variable]]
vi_ext_both[, model     := factor(model, levels = c("50/50 baseline", "Paired"))]

var_order_ext <- vi_ext_base[order(-rel_imp), clean_names_ext[variable]]
vi_ext_both[, var_label := factor(var_label, levels = rev(var_order_ext))]

# Direction flags — same logic as main figure plus new variables
vi_ext_both[, direction := fifelse(
  var_label %in% c("Bioclim PC1", "Biome", "Soil pH", "Elevation"), "down",
  fifelse(var_label %in% c("Forest density", "Dist. to forest",
                           "Soil org. C", "NPP", "Soil bulk density",
                           "Soil sand", "Fire frequency", "Slope"),
          "up", "neutral")
)]

vi_ext_labels <- vi_ext_both[model == "Paired"]

pB_ext <- ggplot(vi_ext_both, aes(x = model, y = rel_imp, group = var_label)) +
  geom_line(aes(colour = direction), linewidth = 0.7, alpha = 0.85) +
  geom_point(aes(colour = direction), size = 2) +
  scale_colour_manual(
    values = c("down"    = col_baseline,
               "up"      = col_ce,
               "neutral" = "grey70"),
    labels = c("down"    = "Climate classifiers (↓)",
               "up"      = "Local variables (↑)",
               "neutral" = "Other"),
    name   = NULL
  ) +
  geom_text_repel(
    data               = vi_ext_labels,
    aes(label = var_label, colour = direction),
    hjust              = 0,
    direction          = "y",
    nudge_x            = 0.08,
    segment.size       = 0.3,
    segment.alpha      = 0.5,
    size               = 2.4,
    show.legend        = FALSE,
    force              = 3,
    min.segment.length = 0
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.45))) +
  labs(
    title = "Variable importance: extended model (baseline vs paired)",
    x     = NULL,
    y     = "Relative importance (% of total)"
  ) +
  nature_theme +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

ggsave(file.path(output_dir, "fig_supp_importance_ext.pdf"),
       pB_ext, width = 120, height = 140, units = "mm",
       device = cairo_pdf)
ggsave(file.path(output_dir, "fig_supp_importance_ext.png"),
       pB_ext, width = 120, height = 140, units = "mm", dpi = 300)
message("Saved: fig_supp_importance_ext.pdf/.png")

# % LaTeX caption for fig_supp_importance_ext:
#   % Relative permutation importance under the extended biophysical model
# % (Williams et al. variables plus slope, elevation, NPP, fire frequency,
#    % and four additional soil properties: sand, clay, bulk density, water
#    % capacity at 33kPa). Colour coding as in main Figure 1b. The Williams
# % et al. variable selection procedure ranked predictors by importance
# % under a 50/50 balanced design; variables that functioned as
# % between-climate classifiers (notably Bioclim PC1, -7.6 percentage
#                                % points) therefore ranked highly and genuinely local variables ranked
# % poorly and were excluded. Under the paired design, soil bulk density
# % (+3.0 pp) and NPP (+1.4 pp) emerge as meaningful local predictors,
# % consistent with compaction from cattle grazing and within-envelope
# % productivity variation as barriers to regeneration independent of
# % climate zone.

# ---- Supplementary: predicted probability distributions ---
# Shows that the paired model shifts probability mass from near-zero
# toward intermediate values (0.1-0.4) relative to the baseline.
# The baseline confidently excludes climatically marginal pixels
# (large spike at ~0); the paired model cannot lean on the climate
# boundary signal and assigns these pixels intermediate probability
# instead. Both models produce similar aggregate area estimates
# (mean prediction × domain area) but through different spatial
# allocations — the distribution plot makes this redistribution
# visible in a way the aggregate cannot.

df_hist <- data.table(
  prob  = c(preds_50_50, preds_ce),
  model = rep(c("50/50 baseline", "Paired"),
              times = c(length(preds_50_50), length(preds_ce)))
)

p_hist <- ggplot(df_hist, aes(x = prob, fill = model)) +
  geom_histogram(
    position  = "identity",
    alpha     = 0.6,
    bins      = 30,
    colour    = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c("50/50 baseline" = col_baseline, "Paired" = col_ce),
    name   = NULL
  ) +
  labs(
    title = "Distribution of predicted P(regrowth): baseline vs paired",
    x     = "Predicted P(regrowth)",
    y     = "Frequency"
  ) +
  nature_theme +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "fig_supp_hist.pdf"),
       p_hist, width = 120, height = 100, units = "mm",
       device = cairo_pdf)
ggsave(file.path(output_dir, "fig_supp_hist.png"),
       p_hist, width = 120, height = 100, units = "mm", dpi = 300)
message("Saved: fig_supp_hist.pdf/.png")


# ---- Supplementary: regrowth pixel predicted probabilities ---------
# Two panels:
#   (a) Mean predicted P(regrowth) on regrowth pixels across models
#   (b) % of regrowth pixels classified as regrowth at 0.5 threshold
# Covers both Williams-var models (50/50, 20/80, 5/95, Paired) and
# extended-var models (50/50 and Paired only).

reg_diag_williams <- data.table(
  model     = c("50/50\n(Williams)", "20/80", "5/95\n(landscape)", "Paired"),
  var_set   = "Williams vars (10)",
  mean_prob = c(mean(reg_preds_50_50), mean(reg_preds_20_80),
                mean(reg_preds_5_95),  mean(reg_preds_ce)),
  pct_above = c(mean(reg_preds_50_50 > 0.5), mean(reg_preds_20_80 > 0.5),
                mean(reg_preds_5_95  > 0.5), mean(reg_preds_ce    > 0.5)) * 100
)

reg_diag_ext <- data.table(
  model     = c("50/50\n(Williams)", "Paired"),
  var_set   = "Extended vars",
  mean_prob = c(mean(reg_preds_ext_50_50), mean(reg_preds_ext_paired)),
  pct_above = c(mean(reg_preds_ext_50_50  > 0.5),
                mean(reg_preds_ext_paired  > 0.5)) * 100
)

reg_diag <- rbind(reg_diag_williams, reg_diag_ext)

model_order <- c("50/50\n(Williams)", "20/80", "5/95\n(landscape)", "Paired")
reg_diag[, model   := factor(model,   levels = model_order)]
reg_diag[, var_set := factor(var_set, levels = c("Williams vars (10)", "Extended vars"))]

cat("\nRegrowth pixel diagnostics:\n")
print(reg_diag)

p_reg_prob <- ggplot(reg_diag,
                     aes(x = model, y = mean_prob * 100,
                         colour = var_set, shape = var_set, group = var_set)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("Williams vars (10)" = col_baseline,
                                 "Extended vars"      = col_ce),
                      name = NULL) +
  scale_shape_manual(values  = c("Williams vars (10)" = 16,
                                 "Extended vars"      = 17),
                     name = NULL) +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "a  Mean predicted P(regrowth) on regrowth pixels",
    x     = NULL,
    y     = "Mean predicted P(regrowth) (%)"
  ) +
  nature_theme +
  theme(legend.position = "bottom")

p_reg_pct <- ggplot(reg_diag,
                    aes(x = model, y = pct_above,
                        colour = var_set, shape = var_set, group = var_set)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  annotate("text", x = 0.6, y = 52, label = "50%", size = 2.5, colour = "grey50") +
  scale_colour_manual(values = c("Williams vars (10)" = col_baseline,
                                 "Extended vars"      = col_ce),
                      name = NULL) +
  scale_shape_manual(values  = c("Williams vars (10)" = 16,
                                 "Extended vars"      = 17),
                     name = NULL) +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "b  % of regrowth pixels classified as regrowth (threshold = 0.5)",
    x     = "Model",
    y     = "% classified as regrowth"
  ) +
  nature_theme +
  theme(legend.position = "bottom")

fig_supp_regrowth <- p_reg_prob / p_reg_pct +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_supp_regrowth.pdf"),
       fig_supp_regrowth, device = cairo_pdf, width = 120, height = 160, units = "mm")
ggsave(file.path(output_dir, "fig_supp_regrowth.png"),
       fig_supp_regrowth, width = 120, height = 160, units = "mm", dpi = 300)
message("Saved: fig_supp_regrowth.pdf / .png")

# LaTeX caption for fig_supp_regrowth:
#   Predicted P(regrowth) evaluated on the regrowth pixel set across
#   model variants. (a) Mean predicted probability. (b) Percentage of
#   regrowth pixels classified as regrowth at the 0.5 decision
#   threshold. Red circles: Williams et al. 10-variable specification;
#   blue triangles: extended biophysical specification. Models shown
#   are the three prevalence variants (50/50, 20/80, 5/95) and the
#   paired (climate-envelope) model; the extended specification is
#   fitted only at 50/50 and paired. Higher prevalence of positives in
#   training directly inflates predicted probabilities and the
#   resulting classification rate, independently of the underlying
#   spatial signal.

# LaTeX caption for fig_supp_hist:
#   Distribution of predicted P(regrowth) across 797,076 random
#   non-regrowth pixels under the 50/50 baseline model (red) and
#   the paired model (blue). The baseline model concentrates
#   probability mass near zero, reflecting confident exclusion of
#   climatically marginal pixels; the paired model redistributes
#   this mass toward intermediate probabilities (0.1-0.4), consistent
#   with reduced reliance on the broad thermal gradient (Bioclim PC1)
#   as a discriminating signal. Both models produce similar aggregate
#   area estimates but through different spatial allocations, as shown
#   in Figure Xd.





# ---- Supplementary: regrowth + non-regrowth pixel diagnostics ------
# Four panels:
#   (a) Mean predicted P(regrowth) — regrowth pixels
#   (b) % classified as regrowth at 0.5 threshold — regrowth pixels
#   (c) Mean predicted P(regrowth) — non-regrowth pixels
#   (d) % classified as regrowth at 0.5 threshold — non-regrowth pixels
# Each panel shows Williams-var models (50/50, 20/80, 5/95, Paired)
# and extended-var models (50/50 and Paired only) as two series.

model_order <- c("50/50\n(Williams)", "20/80", "5/95\n(landscape)", "Paired")

# ---- Regrowth pixel diagnostics ----
reg_diag_williams <- data.table(
  model     = c("50/50\n(Williams)", "20/80", "5/95\n(landscape)", "Paired"),
  var_set   = "Williams vars (10)",
  pixel_set = "Regrowth pixels",
  mean_prob = c(mean(reg_preds_50_50), mean(reg_preds_20_80),
                mean(reg_preds_5_95),  mean(reg_preds_ce)),
  pct_above = c(mean(reg_preds_50_50 > 0.5), mean(reg_preds_20_80 > 0.5),
                mean(reg_preds_5_95  > 0.5), mean(reg_preds_ce    > 0.5)) * 100
)
reg_diag_ext <- data.table(
  model     = c("50/50\n(Williams)", "Paired"),
  var_set   = "Extended vars",
  pixel_set = "Regrowth pixels",
  mean_prob = c(mean(reg_preds_ext_50_50), mean(reg_preds_ext_paired)),
  pct_above = c(mean(reg_preds_ext_50_50  > 0.5),
                mean(reg_preds_ext_paired  > 0.5)) * 100
)

# ---- Non-regrowth pixel diagnostics ----
nonreg_diag_williams <- data.table(
  model     = c("50/50\n(Williams)", "20/80", "5/95\n(landscape)", "Paired"),
  var_set   = "Williams vars (10)",
  pixel_set = "Non-regrowth pixels",
  mean_prob = c(mean(preds_50_50), mean(preds_20_80),
                mean(preds_5_95),  mean(preds_ce)),
  pct_above = c(mean(preds_50_50 > 0.5), mean(preds_20_80 > 0.5),
                mean(preds_5_95  > 0.5), mean(preds_ce    > 0.5)) * 100
)
nonreg_diag_ext <- data.table(
  model     = c("50/50\n(Williams)", "Paired"),
  var_set   = "Extended vars",
  pixel_set = "Non-regrowth pixels",
  mean_prob = c(mean(nonreg_preds_ext_50_50), mean(nonreg_preds_ext_paired)),
  pct_above = c(mean(nonreg_preds_ext_50_50  > 0.5),
                mean(nonreg_preds_ext_paired  > 0.5)) * 100
)

pixel_diag <- rbind(reg_diag_williams, reg_diag_ext,
                    nonreg_diag_williams, nonreg_diag_ext)
pixel_diag[, model     := factor(model,     levels = model_order)]
pixel_diag[, var_set   := factor(var_set,   levels = c("Williams vars (10)", "Extended vars"))]
pixel_diag[, pixel_set := factor(pixel_set, levels = c("Regrowth pixels", "Non-regrowth pixels"))]

cat("\nPixel diagnostics (regrowth + non-regrowth):\n")
print(pixel_diag)

# Shared scale definitions
diag_colours <- c("Williams vars (10)" = col_baseline, "Extended vars" = col_ce)
diag_shapes  <- c("Williams vars (10)" = 16,           "Extended vars" = 17)

make_diag_panel <- function(data, metric, ylab, title_prefix, add_ref = FALSE) {
  p <- ggplot(data,
              aes(x = model, y = .data[[metric]],
                  colour = var_set, shape = var_set, group = var_set)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = diag_colours, name = NULL) +
    scale_shape_manual(values  = diag_shapes,  name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    facet_wrap(~ pixel_set) +
    labs(title = title_prefix, x = NULL, y = ylab) +
    nature_theme +
    theme(legend.position  = "bottom",
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text       = element_text(size = 8, face = "bold"))
  if (add_ref)
    p <- p + geom_hline(yintercept = 50, linetype = "dashed",
                        colour = "grey60", linewidth = 0.4)
  p
}

p_diag_prob <- make_diag_panel(pixel_diag, "mean_prob",
                               "Mean predicted P(regrowth)",
                               "a  Mean predicted P(regrowth)")
p_diag_pct  <- make_diag_panel(pixel_diag, "pct_above",
                               "% classified as regrowth",
                               "b  % classified as regrowth (threshold = 0.5)",
                               add_ref = TRUE)

fig_supp_regrowth <- p_diag_prob / p_diag_pct +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_supp_regrowth.pdf"),
       fig_supp_regrowth, device = cairo_pdf, width = 180, height = 160, units = "mm")
ggsave(file.path(output_dir, "fig_supp_regrowth.png"),
       fig_supp_regrowth, width = 180, height = 160, units = "mm", dpi = 300)
message("Saved: fig_supp_regrowth.pdf / .png")

# LaTeX caption for fig_supp_regrowth:
#   Predicted P(regrowth) on regrowth (left) and non-regrowth (right)
#   pixel sets across model variants. (a) Mean predicted probability.
#   (b) Percentage of pixels classified as regrowth at the 0.5
#   decision threshold. Red circles: Williams et al. 10-variable
#   specification; blue triangles: extended biophysical specification.
#   Models shown are the three prevalence variants (50/50, 20/80,
#   5/95) and the paired (climate-envelope) model; the extended
#   specification is fitted only at 50/50 and paired. For regrowth
#   pixels, higher training prevalence inflates predicted probability
#   and classification rate directly. For non-regrowth pixels, the
#   mirror pattern (false positive rate) reveals how probability mass
#   is redistributed across the full prediction surface as the
#   prevalence assumption changes.

# LaTeX caption for fig_supp_hist:
#   Distribution of predicted P(regrowth) across 797,076 random
#   non-regrowth pixels under the 50/50 baseline model (red) and
#   the paired model (blue). The baseline model concentrates
#   probability mass near zero, reflecting confident exclusion of
#   climatically marginal pixels; the paired model redistributes
#   this mass toward intermediate probabilities (0.1-0.4), consistent
#   with reduced reliance on the broad thermal gradient (Bioclim PC1)
#   as a discriminating signal. Both models produce similar aggregate
#   area estimates but through different spatial allocations, as shown
#   in Figure Xd.


# ---- Supplementary: spatial redistribution by pixel type ---
# Two maps side by side: non-regrowth domain (left) and regrowth
# domain (right), both showing delta = paired - baseline.
# Shared colour scale so the two panels are directly comparable.

# Non-regrowth grid
pred_sample[, lat_bin := round(latitude  * 5) / 5]
pred_sample[, lon_bin := round(longitude * 5) / 5]

grid_nonreg <- pred_sample[, .(
  delta = mean(p_paired) - mean(p_baseline),
  n     = .N
), by = .(lat_bin, lon_bin)]
grid_nonreg <- grid_nonreg[n >= 1]

# Regrowth grid
reg_spatial <- data.table(
  lat_bin    = round(reg_coords$latitude  * 5) / 5,
  lon_bin    = round(reg_coords$longitude * 5) / 5,
  p_baseline = reg_preds_50_50,
  p_paired   = reg_preds_ce
)

grid_reg <- reg_spatial[, .(
  delta = mean(p_paired) - mean(p_baseline),
  n     = .N
), by = .(lat_bin, lon_bin)]
grid_reg <- grid_reg[n >= 1]

# Shared colour scale limits across both panels
delta_lim <- max(abs(c(grid_nonreg$delta, grid_reg$delta)))

make_spatial_panel <- function(grid, title) {
  ggplot(grid, aes(x = lon_bin, y = lat_bin, fill = delta)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = col_baseline,
      mid      = "#fffae6",
      high     = col_ce,
      midpoint = 0,
      limits   = c(-delta_lim, delta_lim),
      name     = "Delta\n(paired \u2212 baseline)"
    ) +
    coord_fixed() +
    labs(title = title, x = "Longitude", y = "Latitude") +
    nature_theme +
    theme(legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.key.width = unit(1.0, "cm"))
}

p_spatial_nonreg <- make_spatial_panel(
  grid_nonreg,
  "a  Non-regrowth domain: paired \u2212 baseline"
)

p_spatial_reg <- make_spatial_panel(
  grid_reg,
  "b  Regrowth domain: paired \u2212 baseline"
)

fig_supp_spatial <- p_spatial_nonreg + p_spatial_reg +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_supp_spatial.pdf"),
       fig_supp_spatial, device = cairo_pdf,
       width = 200, height = 110, units = "mm")
ggsave(file.path(output_dir, "fig_supp_spatial.png"),
       fig_supp_spatial, width = 200, height = 110, units = "mm", dpi = 300)
message("Saved: fig_supp_spatial.pdf / .png")

# LaTeX caption for fig_supp_spatial:
#   Spatial distribution of the shift in predicted P(regrowth) between
#   the paired and 50/50 baseline models (\(\Delta = \hat{p}_\text{paired}
#   - \hat{p}_\text{baseline}\)) across (a) 797,076 random non-regrowth
#   pixels and (b) all regrowth pixels, aggregated to
#   \(\approx 7\,\text{km}\) grid cells. Colour scale is shared across
#   both panels (red: baseline assigns higher probability; blue: paired
#   assigns higher probability; midpoint = 0). The non-regrowth panel
#   shows how the paired model redistributes probability mass spatially
#   relative to the baseline — notably reducing predictions in
#   climatically marginal areas where the baseline relies on the thermal
#   gradient signal. The regrowth panel shows that the paired model
#   assigns broadly higher probability to known regrowth locations,
#   consistent with reduced suppression of within-climate signals.
#   Note: the two domains differ greatly in size (169.9 Mha vs 0.93 Mha)
#   and sampling density; panels are not area-comparable.



# ---- Supplementary: 5/95 baseline vs 5/95 extended paired --
# Two maps: delta = ext_paired_5/95 minus baseline_5/95
# Left: non-regrowth pixels (pred_sample)
# Right: regrowth pixels (reg_coords)
# Shared colour scale across both panels.

# Non-regrowth grid
grid_nonreg_595 <- data.table(
  lat_bin   = round(pred_sample$latitude  * 5) / 5,
  lon_bin   = round(pred_sample$longitude * 5) / 5,
  p_base    = preds_5_95,
  p_ext_pai = preds_ext_paired_5_95
)[, .(
  delta = mean(p_ext_pai) - mean(p_base),
  n     = .N
), by = .(lat_bin, lon_bin)][n >= 1]

# Regrowth grid
grid_reg_595 <- data.table(
  lat_bin   = round(reg_coords$latitude  * 5) / 5,
  lon_bin   = round(reg_coords$longitude * 5) / 5,
  p_base    = reg_preds_5_95,
  p_ext_pai = reg_preds_ext_paired_5_95
)[, .(
  delta = mean(p_ext_pai) - mean(p_base),
  n     = .N
), by = .(lat_bin, lon_bin)][n >= 1]

# Shared colour scale
delta_lim_595 <- max(abs(c(grid_nonreg_595$delta, grid_reg_595$delta)))

make_spatial_panel_595 <- function(grid, title) {
  ggplot(grid, aes(x = lon_bin, y = lat_bin, fill = delta)) +
    geom_tile() +
    scale_fill_gradient2(
      low      = col_baseline,
      mid      = "#fffae6",
      high     = col_ce,
      midpoint = 0,
      limits   = c(-delta_lim_595, delta_lim_595),
      name     = "Delta\n(ext. paired \u2212 baseline)"
    ) +
    coord_fixed() +
    labs(title = title, x = "Longitude", y = "Latitude") +
    nature_theme +
    theme(legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.key.width = unit(1.0, "cm"))
}

p_595_nonreg <- make_spatial_panel_595(
  grid_nonreg_595,
  "a  Non-regrowth: ext. paired 5/95 \u2212 baseline 5/95"
)

p_595_reg <- make_spatial_panel_595(
  grid_reg_595,
  "b  Regrowth: ext. paired 5/95 \u2212 baseline 5/95"
)

fig_supp_spatial_595 <- p_595_nonreg + p_595_reg +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_supp_spatial_595.pdf"),
       fig_supp_spatial_595, device = cairo_pdf,
       width = 200, height = 110, units = "mm")
ggsave(file.path(output_dir, "fig_supp_spatial_595.png"),
       fig_supp_spatial_595, width = 200, height = 110, units = "mm", dpi = 300)
message("Saved: fig_supp_spatial_595.pdf / .png")

# LaTeX caption for fig_supp_spatial_595:
#   Spatial distribution of the difference in predicted P(regrowth)
#   between the extended-variable paired model and the baseline random
#   model, both trained at 5/95 prevalence
#   (\(\Delta = \hat{p}_\text{ext. paired} - \hat{p}_\text{baseline}\)),
#   across (a) 797,076 random non-regrowth pixels and (b) all regrowth
#   pixels, aggregated to \(\approx 7\,\text{km}\) grid cells.
#   Blue indicates regions where the extended paired model assigns
#   higher predicted probability; red indicates regions where the
#   baseline assigns higher probability. At 5/95 prevalence both
#   models are strongly suppressed relative to 50/50, but the spatial
#   pattern of residual disagreement reveals where the additional
#   biophysical variables and climate-envelope control jointly
#   redistribute probability mass.


# ---- Supplementary: area sensitivity — three design lines --
# Three lines across prevalence points (50/50, 20/80, 5/95):
#   1. Baseline random sampling (Williams vars)
#   2. Williams-var paired (multi-pass NN matched negatives)
#   3. Extended-var paired
# x-axis is prevalence (pos fraction); y is total area (Mha).
# Panel A is not changed — this is a standalone supplementary figure.

prev_levels  <- c("50/50", "20/80", "5/95")
prev_numeric <- c(0.50,    0.20,    0.05)

area_supp_dt <- data.table(
  prev_label = rep(prev_levels, 3),
  prev_pos   = rep(prev_numeric, 3),
  design     = c(rep("Baseline random",   3),
                 rep("Williams paired",   3),
                 rep("Extended paired",   3)),
  area_Mha   = c(
    # Baseline random
    area_dt[model == "50/50",    area_total_Mha],
    area_dt[model == "20/80",    area_total_Mha],
    area_dt[model == "5/95",     area_total_Mha],
    # Williams paired
    area_dt[model == "Paired",       area_total_Mha],
    area_dt[model == "Paired_20/80", area_total_Mha],
    area_dt[model == "Paired_5/95",  area_total_Mha],
    # Extended paired
    area_dt[model == "Paired_ext",       area_total_Mha],
    area_dt[model == "Paired_ext_20/80", area_total_Mha],
    area_dt[model == "Paired_ext_5/95",  area_total_Mha]
  )
)

area_supp_dt[, prev_label := factor(prev_label, levels = prev_levels)]
area_supp_dt[, design     := factor(design,
  levels = c("Baseline random", "Williams paired", "Extended paired"))]

design_colours <- c(
  "Baseline random"  = col_baseline,
  "Williams paired"  = col_ce,
  "Extended paired"  = "#4DAF4A"
)
design_shapes <- c(
  "Baseline random"  = 16,
  "Williams paired"  = 17,
  "Extended paired"  = 15
)
design_ltys <- c(
  "Baseline random"  = "solid",
  "Williams paired"  = "dashed",
  "Extended paired"  = "dotted"
)

p_area_supp <- ggplot(area_supp_dt,
                      aes(x = prev_label, y = area_Mha,
                          colour = design, shape = design,
                          linetype = design, group = design)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = design_colours, name = NULL) +
  scale_shape_manual(values  = design_shapes,  name = NULL) +
  scale_linetype_manual(values = design_ltys,  name = NULL) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Estimated total area by prevalence and training design",
    x     = "Prevalence (regrowth : non-regrowth)",
    y     = "Estimated area (Mha)"
  ) +
  nature_theme +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

ggsave(file.path(output_dir, "fig_supp_area.pdf"),
       p_area_supp, device = cairo_pdf, width = 120, height = 100, units = "mm")
ggsave(file.path(output_dir, "fig_supp_area.png"),
       p_area_supp, width = 120, height = 100, units = "mm", dpi = 300)
message("Saved: fig_supp_area.pdf / .png")

# LaTeX caption for fig_supp_area:
#   Estimated total regeneration area (non-regrowth domain +
#   regrowth domain) as a function of training prevalence under
#   three modelling designs: baseline random sampling (red, circles),
#   Williams-variable paired design with multi-pass nearest-neighbour
#   matched negatives (blue, triangles), and extended biophysical
#   paired design (green, squares). Points at 50/50 for the paired
#   designs correspond to the 50\% positive models trained in R5b/R6b;
#   points at 20/80 and 5/95 use the R7 prevalence-sensitivity models.


# ---- Supplementary: z-score normalized DiD spatial map -----
# Normalizes each paired-minus-baseline map to z-scores before
# differencing, so the DiD reflects differences in spatial ordering
# rather than differences in scale between prevalence levels.
#
# DiD = z(delta_50/50) - z(delta_5/95)
# where z(x) = (x - mean(x)) / sd(x) computed over all pixels
# separately for non-regrowth and regrowth domains.

# ---- Non-regrowth ----
nonreg_did_dt <- data.table(
  lat_bin     = round(pred_sample$latitude  * 5) / 5,
  lon_bin     = round(pred_sample$longitude * 5) / 5,
  delta_50_50 = preds_ce          - preds_50_50,
  delta_5_95  = preds_paired_5_95 - preds_5_95
)

# Bin first, then z-score the binned means so the normalization
# reflects the spatial distribution not individual pixel noise
grid_did_nonreg <- nonreg_did_dt[, .(
  d_50_50 = mean(delta_50_50),
  d_5_95  = mean(delta_5_95),
  n       = .N
), by = .(lat_bin, lon_bin)][n >= 1]

grid_did_nonreg[, z_50_50 := (d_50_50 - mean(d_50_50)) / sd(d_50_50)]
grid_did_nonreg[, z_5_95  := (d_5_95  - mean(d_5_95))  / sd(d_5_95)]
grid_did_nonreg[, did     := z_50_50 - z_5_95]

# ---- Regrowth ----
reg_did_dt <- data.table(
  lat_bin     = round(reg_coords$latitude  * 5) / 5,
  lon_bin     = round(reg_coords$longitude * 5) / 5,
  delta_50_50 = reg_preds_ce          - reg_preds_50_50,
  delta_5_95  = reg_preds_paired_5_95 - reg_preds_5_95
)

grid_did_reg <- reg_did_dt[, .(
  d_50_50 = mean(delta_50_50),
  d_5_95  = mean(delta_5_95),
  n       = .N
), by = .(lat_bin, lon_bin)][n >= 1]

grid_did_reg[, z_50_50 := (d_50_50 - mean(d_50_50)) / sd(d_50_50)]
grid_did_reg[, z_5_95  := (d_5_95  - mean(d_5_95))  / sd(d_5_95)]
grid_did_reg[, did     := z_50_50 - z_5_95]

library(scales)

# 1. Lock limits to 5 and set saturation point at 2.5
did_limit <- 5 
sat_lim   <- 2.5 

# 2. Rescale custom breakpoints to a 0-1 range for the gradient mapping
color_vals <- scales::rescale(c(-did_limit, -sat_lim, 0, sat_lim, did_limit))

# 3. Apply scale_fill_gradientn with explicit breaks
make_did_panel <- function(grid, title) {
  ggplot(grid, aes(x = lon_bin, y = lat_bin, fill = did)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c(col_baseline, col_baseline, "#fffae6", col_ce, col_ce),
      values = color_vals,
      limits = c(-did_limit, did_limit),
      breaks = c(-5, -2.5, 0, 2.5, 5), # Hardcoded to show exactly these numbers
      oob    = scales::squish,
      name   = "DiD\n(z-score\nunits)"
    ) +
    coord_fixed() +
    labs(title = title, x = "Longitude", y = "Latitude") +
    nature_theme +
    theme(legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.key.width = unit(1.0, "cm"))
}

p_did_nonreg <- make_did_panel(
  grid_did_nonreg,
  "a  Non-regrowth: z(\u0394 50/50) \u2212 z(\u0394 5/95)"
)

p_did_reg <- make_did_panel(
  grid_did_reg,
  "b  Regrowth: z(\u0394 50/50) \u2212 z(\u0394 5/95)"
)

fig_supp_did <- p_did_nonreg + p_did_reg +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_supp_did.pdf"),
       fig_supp_did, device = cairo_pdf,
       width = 200, height = 110, units = "mm")
ggsave(file.path(output_dir, "fig_supp_did.png"),
       fig_supp_did, width = 200, height = 110, units = "mm", dpi = 300)
message("Saved: fig_supp_did.pdf / .png")


# LaTeX caption:
#
# Difference-in-differences of z-score normalised spatial
# redistribution maps for (a) non-regrowth and (b) regrowth pixels.
# Each delta map (paired minus baseline) is standardised to z-scores
# within its prevalence level before differencing:
# DiD = z(Delta_50/50) - z(Delta_5/95). Normalisation removes scale
# differences between prevalence levels so that only differences in
# spatial ordering are shown. Custard regions indicate pixels whose
# relative rank in the pairing effect is stable across prevalence
# levels; blue regions rank higher at 50/50 than at 5/95; red regions
# rank higher at 5/95. Within-zone Spearman correlations between
# z(Delta_50/50) and z(Delta_5/95) are reported in fig_supp_split.

# Interpretation note (not for publication):
#
# The non-uniform spatial structure confirms that prevalence and
# pairing cannot be treated as separable sensitivity axes.
# The Amazon interior is relatively stable across both design choices
# because the concentration of observed regrowth there during the
# study period is strong enough that all models learn biome membership
# as a proxy for potential — reinforced by the fact that non-regrowth
# pixels in the Amazon are spatially confined to narrow road and river
# corridors, leaving little within-biome variation to disagree over.
# The SE shows the most structure because it has more non-regrowth
# area across a wider range of conditions, so the choice of prevalence
# and negative sampling design jointly determines which pixels get
# elevated predictions.

# ---- Supplementary: Climate Zone Split & Correlations ------

grid_did_nonreg[, zone := fcase(
  lon_bin < (-38 + (lat_bin * 0.83) - (0.00008 * lat_bin^4) + (0.008 * lat_bin^2) + (0.0002 * lat_bin^3) ), "NW",
  default = "SE"
)]
grid_did_reg[, zone := fcase(
  lon_bin < (-38 + (lat_bin * 0.83) - (0.00008 * lat_bin^4) + (0.008 * lat_bin^2) + (0.0002 * lat_bin^3) ), "NW",
  default = "SE"
)]


get_top_k_overlap <- function(dt, k = 0.10) {
  calc <- function(d) {
    q_50 <- quantile(d$z_50_50, 1 - k)
    q_95 <- quantile(d$z_5_95,  1 - k)
    top_50 <- d$z_50_50 >= q_50
    top_95 <- d$z_5_95  >= q_95
    sum(top_50 & top_95) / sum(top_50)
  }
  rbind(
    data.table(zone = "Overall", overlap = calc(dt)),
    data.table(zone = "NW",      overlap = calc(dt[zone == "NW"])),
    data.table(zone = "SE",      overlap = calc(dt[zone == "SE"]))
  )
}


# Compute overlap statistics
overlap_nonreg <- get_top_k_overlap(grid_did_nonreg, k = 0.10)
overlap_reg    <- get_top_k_overlap(grid_did_reg,    k = 0.10)

get_annotations_v2 <- function(dt, overlap_dt) {
  r_all <- cor(dt$z_50_50, dt$z_5_95, method = "spearman")
  r_nw  <- dt[zone == "NW", cor(z_50_50, z_5_95, method = "spearman")]
  r_se  <- dt[zone == "SE", cor(z_50_50, z_5_95, method = "spearman")]
  
  ov_all <- overlap_dt[zone == "Overall", overlap]
  ov_nw  <- overlap_dt[zone == "NW",      overlap]
  ov_se  <- overlap_dt[zone == "SE",      overlap]
  
  cor_label <- sprintf(
    "Spearman r\nOverall: %.2f\nNW:      %.2f\nSE:      %.2f",
    r_all, r_nw, r_se
  )
  overlap_label <- sprintf(
    "Top-10%% pixel overlap\nOverall: %.2f\nNW:      %.2f\nSE:      %.2f",
    ov_all, ov_nw, ov_se
  )
  list(cor = cor_label, overlap = overlap_label)
}

annots_nonreg <- get_annotations_v2(grid_did_nonreg, overlap_nonreg)
annots_reg    <- get_annotations_v2(grid_did_reg,    overlap_reg)

make_split_panel <- function(grid, title_text, annots) {
  
  x_min <- min(grid$lon_bin, na.rm = TRUE)
  x_max <- max(grid$lon_bin, na.rm = TRUE)
  y_min <- min(grid$lat_bin, na.rm = TRUE)
  y_max <- max(grid$lat_bin, na.rm = TRUE)
  
  ggplot(grid, aes(x = lon_bin, y = lat_bin, fill = zone)) +
    geom_tile() +
    scale_fill_manual(
      values = c("NW" = col_ce, "SE" = col_baseline),
      name   = "Climate zone"
    ) +
    geom_abline(
      intercept = -38,
      slope     = 0.8,
      linetype  = "dashed",
      colour    = "white",
      linewidth = 0.5
    ) +
    # Correlations — top right
    annotate("label",
             x = x_max, y = y_max,
             label         = annots$cor,
             hjust         = 1, vjust = 1,
             size          = 2.6,
             family        = "mono",
             fill          = "white",
             label.size    = 0.3,
             label.padding = unit(0.3, "lines")) +
    # Top-k overlap — bottom left
    annotate("label",
             x = x_min, y = y_min,
             label         = annots$overlap,
             hjust         = 0, vjust = 0,
             size          = 2.6,
             family        = "mono",
             fill          = "white",
             label.size    = 0.3,
             label.padding = unit(0.3, "lines")) +
    coord_fixed() +
    labs(title = title_text, x = "Longitude", y = "Latitude") +
    nature_theme +
    theme(legend.position  = "bottom",
          legend.key.size  = unit(0.35, "cm"),
          legend.key.width = unit(0.5, "cm"))
}

p_split_nonreg <- make_split_panel(
  grid_did_nonreg,
  "a  Non-regrowth: climate zone split",
  annots_nonreg
)

p_split_reg <- make_split_panel(
  grid_did_reg,
  "b  Regrowth: climate zone split",
  annots_reg
)

fig_supp_split <- p_split_nonreg + p_split_reg +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_supp_split.pdf"),
       fig_supp_split, device = cairo_pdf,
       width = 200, height = 110, units = "mm")
ggsave(file.path(output_dir, "fig_supp_split.png"),
       fig_supp_split, width = 200, height = 110, units = "mm", dpi = 300)
message("Saved: fig_supp_split.pdf / .png")

# LaTeX caption:
#
# \caption{\textbf{Supplementary Fig.~\thefigure{} $|$ Geographic stability of 
#   the paired-vs-baseline spatial redistribution effect across prevalence levels.} 
#   Each panel shows whether the spatial pattern of the pairing effect 
#   (paired minus baseline predictions) is consistent between the 50/50 and 5/95 
#   prevalence levels, after z-score normalising each delta map within its 
#   prevalence level. Tile colour indicates climate zone assignment (NW: Amazon and 
#                                                                    Cerrado transition; SE: Atlantic Forest and Caatinga). Annotations report the 
#   Spearman correlation between $z(\Delta_{50/50})$ and $z(\Delta_{5/95})$---the 
#   degree to which pixels that are strongly affected by paired sampling at one 
#   prevalence level are also strongly affected at the other---and the top-10\% 
#   overlap, defined as the fraction of pixels in the top decile of the pairing 
#   effect under 50/50 prevalence that also fall in the top decile under 5/95 
#   prevalence. A value of 1.0 would indicate that the same pixels are most 
#   affected by the pairing design regardless of prevalence; the random baseline 
#   is 0.10.}

# Interpretation note (not for publication):
#
# The NW zone shows high stability across all comparisons (0.73--0.80
# for regrowth), reflecting the dominance of the Amazon distributional
# signal rather than model quality. The SE tells the opposite story:
# prevalence alone drops regrowth overlap to 0.29, and the joint
# change in prevalence and extended paired design reduces it to 0.13
# — fewer than one in eight priority pixels are stable. This
# asymmetry between zones is the clearest evidence that a single
# model trained across Brazil learns a compromise between two
# fundamentally different regeneration regimes, and that the
# headline area estimates are driven primarily by how well the
# model captures the Amazon rather than Brazil as a whole.

# Global domain areas (in hectares)

area_nr_total <- total_domain_area_ha - regrowth_domain_area_ha
area_r_total  <- regrowth_domain_area_ha

# Helper: Compute Area-Weighted Quantile
get_wtd_quantile <- function(vals, wts, prob = 0.95) {
  ok <- !is.na(vals) & !is.na(wts)
  vals <- vals[ok]
  wts <- wts[ok]
  
  ord <- order(vals)
  v_ord <- vals[ord]
  w_ord <- wts[ord]
  
  cum_w <- cumsum(w_ord) / sum(w_ord)
  return(v_ord[which.max(cum_w >= prob)])
}

# Helper: Compute Area-Weighted Percentiles (0 to 100)
get_wtd_percentile <- function(vals, wts) {
  ord <- order(vals)
  w_ord <- wts[ord]
  
  # Calculate cumulative area (using mid-point of the weight to avoid ties pushing to 100)
  cum_area <- cumsum(w_ord) - (w_ord / 2)
  
  # Map back to original order and scale 0-100
  pctiles <- numeric(length(vals))
  pctiles[ord] <- (cum_area / sum(wts)) * 100
  return(pctiles)
}

# Pixel-level zone assignment (using data.table's fast fifelse)
pred_sample[, zone := fifelse(
  longitude < (-38 + (latitude * 0.83) - (0.00008 * latitude^4) + (0.008 * latitude^2) + (0.0002 * latitude^3)), "NW", "SE"
)]
reg_zone_vec <- fifelse(
  reg_coords$longitude < (-38 + (reg_coords$latitude * 0.83) - (0.00008 * reg_coords$latitude^4) + (0.008 * reg_coords$latitude^2) + (0.0002 * reg_coords$latitude^3)), "NW", "SE"
)
combined_zone_vec <- c(pred_sample$zone, reg_zone_vec)

comparisons <- list(
  list(label = "Prevalence only (base 50/50 vs base 5/95)",
       A_nr = preds_50_50,       A_r = reg_preds_50_50,
       B_nr = preds_5_95,        B_r = reg_preds_5_95),
  list(label = "Design only @ 50/50 (base vs ext-paired)",
       A_nr = preds_50_50,       A_r = reg_preds_50_50,
       B_nr = preds_ext_paired,  B_r = reg_preds_ext_paired),
  list(label = "Design only @ 5/95 (base vs ext-paired)",
       A_nr = preds_5_95,        A_r = reg_preds_5_95,
       B_nr = preds_ext_paired_5_95, B_r = reg_preds_ext_paired_5_95),
  list(label = "Joint (base 50/50 vs ext-paired 5/95)",
       A_nr = preds_50_50,       A_r = reg_preds_50_50,
       B_nr = preds_ext_paired_5_95, B_r = reg_preds_ext_paired_5_95),
  list(label = "Full design sensitivity (ext-paired 50/50 vs ext-paired 5/95)",
       A_nr = preds_ext_paired,  A_r = reg_preds_ext_paired,
       B_nr = preds_ext_paired_5_95, B_r = reg_preds_ext_paired_5_95)
)

# Weights are constant across comparisons — compute once
wt_nr <- area_nr_total / length(comparisons[[1]]$A_nr)  
wt_r  <- area_r_total  / length(comparisons[[1]]$A_r)   

all_stats <- rbindlist(lapply(comparisons, function(comp) {
  
  comb_A <- c(comp$A_nr, comp$A_r)
  comb_B <- c(comp$B_nr, comp$B_r)
  comb_wts <- c(rep(wt_nr, length(comp$A_nr)), rep(wt_r, length(comp$A_r)))
  
  get_zone_stats <- function(zone_name) {
    if (zone_name == "Overall") {
      cA <- comb_A; cB <- comb_B; cW <- comb_wts
      rA <- comp$A_r; rB <- comp$B_r
    } else {
      idx_c <- combined_zone_vec == zone_name
      idx_r <- reg_zone_vec == zone_name
      
      cA <- comb_A[idx_c]; cB <- comb_B[idx_c]; cW <- comb_wts[idx_c]
      rA <- comp$A_r[idx_r]; rB <- comp$B_r[idx_r]
    }
    
    # Top-5% threshold from the area-weighted combined landscape for this zone
    q_A <- get_wtd_quantile(cA, cW, prob = 0.95)
    q_B <- get_wtd_quantile(cB, cW, prob = 0.95)
    
    # Combined: fraction of model-A top-5% AREA also top-5% in model-B
    # Combined: fraction of model-A top-5% AREA also top-5% in model-B
    calc_combined <- function(a, b, w) {
      if (length(a) == 0) return(list(overlap = NA, pctile_shift = NA, spearman_r = NA))
      
      # Overlap (already area-weighted in your previous fix)
      overlap <- (sum(w[a >= q_A & b >= q_B], na.rm = TRUE) / sum(w[a >= q_A], na.rm = TRUE)) * 100
      
      # True Landscape Percentiles (0-100)
      pctile_a <- get_wtd_percentile(a, w)
      pctile_b <- get_wtd_percentile(b, w)
      
      # Area-weighted Mean Absolute Percentile Shift
      pctile_shift <- weighted.mean(abs(pctile_a - pctile_b), w = w, na.rm = TRUE)
      
      # Area-weighted Spearman (Pearson correlation of the weighted percentiles)
      spearman_r <- cov.wt(cbind(pctile_a, pctile_b), wt = w, cor = TRUE)$cor[1, 2]
      
      list(overlap = overlap, pctile_shift = pctile_shift, spearman_r = spearman_r)
    }
    
    # Regrowth recall: Note we use unweighted length(a) for the denominator 
    # here because we are explicitly counting regrowth *pixels*, and every 
    # regrowth pixel inherently carries the exact same area weight (wt_r).
    calc_regrowth_recall <- function(a, b) {
      if (length(a) == 0) return(list(overlap = NA, pctile_shift = NA, spearman_r = NA))
      
      overlap      <- (sum(a >= q_A & b >= q_B, na.rm = TRUE) / length(a)) * 100
      pctile_shift <- (mean(abs(rank(a) - rank(b))) / length(a)) * 100
      spearman_r   <- cor(a, b, method = "spearman")
      
      list(overlap = overlap, pctile_shift = pctile_shift, spearman_r = spearman_r)
    }
    
    m_comb <- calc_combined(cA, cB, cW)
    m_reg  <- calc_regrowth_recall(rA, rB)
    
    rbind(
      data.table(sample_type = "Combined", zone = zone_name, comparison = comp$label,
                 overlap = m_comb$overlap, pctile_shift = m_comb$pctile_shift,
                 spearman_r = m_comb$spearman_r),
      data.table(sample_type = "Regrowth", zone = zone_name, comparison = comp$label,
                 overlap = m_reg$overlap,  pctile_shift = m_reg$pctile_shift,
                 spearman_r = m_reg$spearman_r)
    )
  }
  
  rbind(get_zone_stats("Overall"), get_zone_stats("NW"), get_zone_stats("SE"))
}))

# Pivot to wide format
table_wide <- dcast(
  all_stats,
  sample_type + comparison ~ factor(zone, levels = c("Overall", "NW", "SE")),
  value.var = c("overlap", "pctile_shift", "spearman_r")
)

print(table_wide)

# Set row order and clean labels
row_order <- c(
  "Prevalence only (base 50/50 vs base 5/95)",
  "Design only @ 50/50 (base vs ext-paired)",
  "Design only @ 5/95 (base vs ext-paired)",
  "Joint (base 50/50 vs ext-paired 5/95)",
  "Full design sensitivity (ext-paired 50/50 vs ext-paired 5/95)"
)
table_wide[, comparison  := factor(comparison,  levels = row_order)]
table_wide[, sample_type := factor(sample_type, levels = c("Combined", "Regrowth"))]
setorder(table_wide, sample_type, comparison)

table_wide[, comparison := rep(c(
  "Prevalence only",
  "Design only (50/50)",
  "Design only (5/95)",
  "Joint",
  "Prevalence only (paired)"
), 2)]

# Helper: generate LaTeX rows (SE column bolded per metric group)
make_tex_rows <- function(dt_subset) {
  paste(
    apply(dt_subset, 1, function(r) {
      sprintf(
        "%s & %.1f & %.1f & \\textbf{%.1f} & %.1f & %.1f & \\textbf{%.1f} & %.3f & %.3f & \\textbf{%.3f} \\\\",
        r[2],
        as.numeric(r[3]), as.numeric(r[4]), as.numeric(r[5]),
        as.numeric(r[6]), as.numeric(r[7]), as.numeric(r[8]),
        as.numeric(r[9]), as.numeric(r[10]), as.numeric(r[11])
      )
    }),
    collapse = "\n"
  )
}

# Build LaTeX table (no siunitx — plain r columns with @{} outer padding removed)
tex_table <- sprintf(
  '\\begin{tabular}{@{}p{4.0cm}rrrrrrrrr@{}}
\\toprule
 & \\multicolumn{3}{c}{Overlap (\\%%)}
 & \\multicolumn{3}{c}{Percentile Shift}
 & \\multicolumn{3}{c}{Spearman $r$} \\\\
\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}
\\textbf{Comparison}
 & \\textbf{All} & \\textbf{NW} & \\textbf{SE}
 & \\textbf{All} & \\textbf{NW} & \\textbf{SE}
 & \\textbf{All} & \\textbf{NW} & \\textbf{SE} \\\\
\\midrule
\\multicolumn{10}{@{}l}{\\textit{Combined landscape} --- top-5\\%% overlap} \\\\[2pt]
%s
\\addlinespace[6pt]
\\multicolumn{10}{@{}l}{\\textit{Candidate regrowth} --- \\%% of zone regrowth pixels in both models\\textquotesingle{} top-5\\%%} \\\\[2pt]
%s
\\bottomrule
\\end{tabular}',
  make_tex_rows(table_wide[sample_type == "Combined"]),
  make_tex_rows(table_wide[sample_type == "Regrowth"])
)

writeLines(tex_table, file.path(output_dir, "table_stability.tex"))
message("Saved: table_stability.tex")

# ============================================================
# LaTeX caption (Included in the \caption{} block above):
#
# Stability of predicted spatial hierarchy across model specifications. 
# Overlap and rank metrics for the combined landscape (Top-10% threshold) 
# and candidate regrowth sites (Top-5% threshold), separated by climate 
# zone (Overall/NW/SE). Overlap is the percentage of top-tier pixels 
# under the baseline model retained under the alternative specification. 
# Percentile shift measures the mean absolute displacement in rank, 
# scaled to a 100-point index. Comparisons isolate the effect of training 
# prevalence alone, absence sampling design alone, and both changes jointly. 
# NW: Amazon and Cerrado transition; SE: Atlantic Forest and Caatinga.
# ============================================================

# Interpretation note (Not for publication):
#
# All Combined metrics are area-weighted to correct for the ~180x size
# difference between the non-regrowth (169 Mha) and regrowth (0.9 Mha)
# domains. Without weighting, the regrowth domain would be massively
# overrepresented in any unweighted top-5% threshold. Regrowth recall
# remains unweighted because every regrowth pixel carries the same
# area weight (wt_r), so weighting would cancel.
#
# The NW zone shows high stability across all metrics because the dense
# Amazon training data acts as an anchor. The models generally agree on
# priorities here (Regrowth overlap ~[UPDATE]) with relatively minor
# rank shuffling (~[UPDATE] percentiles).
#
# The SE tells the opposite story. Under the joint change (Row 4 in
# Regrowth), the SE regrowth overlap drops to [UPDATE]%, meaning over
# [UPDATE]% of the top restoration sites under the baseline are rejected
# under the paired 5/95 model. The average priority pixel in the SE is
# displaced by [UPDATE] percentiles.
#
# This asymmetry confirms that in the fragmented Atlantic Forest and
# Caatinga biomes, the priority map is an artifact of global model
# design rather than a stable ecological signal. The bolded SE columns
# provide the numerical evidence for the spatial redistribution maps
# shown in the main text.
# ============================================================


# ---- Supplementary: National Priority Allocation and Additive Survival ----
# Evaluates how national top-5% priorities are spatially distributed across 
# climate zones (NW/SE), and decomposes the overall survival rate into 
# additive contributions from each zone.

message("\n--- Calculating Additive National Allocation by Zone ---")

# Safe label mapping to prevent alphabetical sorting bugs
comp_labels <- c(
  "Prevalence only (base 50/50 vs base 5/95)" = "Prevalence only",
  "Design only @ 50/50 (base vs ext-paired)"  = "Design only (50/50)",
  "Design only @ 5/95 (base vs ext-paired)"   = "Design only (5/95)",
  "Joint (base 50/50 vs ext-paired 5/95)"     = "Joint",
  "Full design sensitivity (ext-paired 50/50 vs ext-paired 5/95)" = "Prevalence only (paired)"
)

allocation_stats <- rbindlist(lapply(comparisons, function(comp) {
  
  comb_A <- c(comp$A_nr, comp$A_r)
  comb_B <- c(comp$B_nr, comp$B_r)
  comb_wts <- c(rep(wt_nr, length(comp$A_nr)), rep(wt_r, length(comp$A_r)))
  is_regrowth <- c(rep(FALSE, length(comp$A_nr)), rep(TRUE, length(comp$A_r)))
  
  # 1. Establish NATIONAL Top 5% Thresholds
  q_A_nat <- get_wtd_quantile(comb_A, comb_wts, prob = 0.95)
  q_B_nat <- get_wtd_quantile(comb_B, comb_wts, prob = 0.95)
  
  top_A <- comb_A >= q_A_nat
  top_B <- comb_B >= q_B_nat
  
  comp_short <- comp_labels[[comp$label]]
  
  get_alloc_stats <- function(z) {
    if (z == "Overall") {
      z_mask <- rep(TRUE, length(comb_A))
    } else {
      z_mask <- combined_zone_vec == z
    }
    
    # --- COMBINED LANDSCAPE (Area-Weighted) ---
    area_top_A_nat <- sum(comb_wts[top_A], na.rm = TRUE)
    area_top_B_nat <- sum(comb_wts[top_B], na.rm = TRUE)
    
    area_top_A_z <- sum(comb_wts[top_A & z_mask], na.rm = TRUE)
    area_top_B_z <- sum(comb_wts[top_B & z_mask], na.rm = TRUE)
    area_surv_z  <- sum(comb_wts[top_A & top_B & z_mask], na.rm = TRUE)
    
    share_A_comb <- (area_top_A_z / area_top_A_nat) * 100
    share_B_comb <- (area_top_B_z / area_top_B_nat) * 100
    # Additive Survival: Denominator is always the NATIONAL total for Model A
    surv_comb    <- ifelse(area_top_A_nat > 0, (area_surv_z / area_top_A_nat) * 100, NA)
    
    # --- REGROWTH PIXELS (Unweighted Count) ---
    reg_top_A_nat <- sum(top_A & is_regrowth, na.rm = TRUE)
    reg_top_B_nat <- sum(top_B & is_regrowth, na.rm = TRUE)
    
    reg_top_A_z <- sum(top_A & z_mask & is_regrowth, na.rm = TRUE)
    reg_top_B_z <- sum(top_B & z_mask & is_regrowth, na.rm = TRUE)
    reg_surv_z  <- sum(top_A & top_B & z_mask & is_regrowth, na.rm = TRUE)
    
    share_A_reg <- ifelse(reg_top_A_nat > 0, (reg_top_A_z / reg_top_A_nat) * 100, 0)
    share_B_reg <- ifelse(reg_top_B_nat > 0, (reg_top_B_z / reg_top_B_nat) * 100, 0)
    # Additive Survival: Denominator is always the NATIONAL total for Model A
    surv_reg    <- ifelse(reg_top_A_nat > 0, (reg_surv_z / reg_top_A_nat) * 100, NA)
    
    rbind(
      data.table(sample_type = "Combined", zone = z, comparison = comp_short,
                 share_A = share_A_comb, share_B = share_B_comb, survival = surv_comb),
      data.table(sample_type = "Regrowth", zone = z, comparison = comp_short,
                 share_A = share_A_reg, share_B = share_B_reg, survival = surv_reg)
    )
  }
  
  rbind(get_alloc_stats("Overall"), get_alloc_stats("NW"), get_alloc_stats("SE"))
}))

# Pivot to wide format
alloc_wide <- dcast(
  allocation_stats,
  sample_type + comparison ~ factor(zone, levels = c("Overall", "NW", "SE")),
  value.var = c("share_A", "share_B", "survival")
)

# Define exact row order
row_order <- c(
  "Prevalence only",
  "Design only (50/50)",
  "Design only (5/95)",
  "Joint",
  "Prevalence only (paired)"
)
alloc_wide[, comparison := factor(comparison, levels = row_order)]
alloc_wide[, sample_type := factor(sample_type, levels = c("Combined", "Regrowth"))]
setorder(alloc_wide, sample_type, comparison)

# Helper: generate LaTeX rows (Dropping the redundant 100% "All" share columns)
make_alloc_tex_rows <- function(dt_subset) {
  paste(
    apply(dt_subset, 1, function(r) {
      # r layout: sample_type, comparison, share_A_All, A_NW, A_SE, share_B_All, B_NW, B_SE, surv_All, surv_NW, surv_SE
      # We drop indices 3 (share_A_All) and 6 (share_B_All) to save space
      sprintf(
        "%s & %.1f & \\textbf{%.1f} & %.1f & \\textbf{%.1f} & %.1f & %.1f & \\textbf{%.1f} \\\\",
        r[2],
        as.numeric(r[4]), as.numeric(r[5]),   # Share A (NW, SE)
        as.numeric(r[7]), as.numeric(r[8]),   # Share B (NW, SE)
        as.numeric(r[9]), as.numeric(r[10]), as.numeric(r[11])  # Survival (Overall, NW, SE)
      )
    }),
    collapse = "\n"
  )
}

# Build LaTeX table (Tabular only, no wrappers, for \input{})
tex_alloc_table <- sprintf(
  "\\begin{tabular}{@{}lrrrrrrr@{}}
\\toprule
& \\multicolumn{2}{c}{Share Model A (\\%%)} & \\multicolumn{2}{c}{Share Model B (\\%%)} & \\multicolumn{3}{c}{Retained National Share (\\%%)} \\\\
\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-8}
\\textbf{Comparison} & \\textbf{NW} & \\textbf{SE} & \\textbf{NW} & \\textbf{SE} & \\textbf{Overall} & \\textbf{NW} & \\textbf{SE} \\\\
\\midrule
\\multicolumn{8}{@{}l}{\\textit{Combined landscape} --- Area-weighted composition of national top-5\\%%} \\\\[2pt]
%s
\\addlinespace[6pt]
\\multicolumn{8}{@{}l}{\\textit{Candidate regrowth} --- Composition of national top-5\\%% regrowth pixels} \\\\[2pt]
%s
\\bottomrule
\\end{tabular}",
  make_alloc_tex_rows(alloc_wide[sample_type == "Combined"]),
  make_alloc_tex_rows(alloc_wide[sample_type == "Regrowth"])
)

writeLines(tex_alloc_table, file.path(output_dir, "table_national_allocation.tex"))
message("Saved: table_national_allocation.tex")



# ---- Supplementary: regrowth + non-regrowth diagnostics v2 --
# Updated version of fig_supp_regrowth with three lines instead of two:
#   1. Baseline random (Williams vars) — 50/50, 20/80, 5/95
#   2. Williams-var paired             — 50/50, 20/80, 5/95
#   3. Extended-var paired             — 50/50, 20/80, 5/95
# x-axis: prevalence point (50/50 → 20/80 → 5/95)
# Two panels: (a) regrowth pixels  (b) non-regrowth pixels
# Each panel shows mean P(regrowth) and % classified > 0.5.

prev_labels3 <- c("50/50", "20/80", "5/95")

make_diag3 <- function(mean_vals, pct_vals, design_label, pixel_label) {
  data.table(
    prev_label = prev_labels3,
    design     = design_label,
    pixel_set  = pixel_label,
    mean_prob  = mean_vals,
    pct_above  = pct_vals
  )
}

diag3 <- rbindlist(list(
  # ---- Regrowth pixels ----
  make_diag3(
    c(mean(reg_preds_50_50),        mean(reg_preds_20_80),        mean(reg_preds_5_95)),
    c(mean(reg_preds_50_50 > 0.5),  mean(reg_preds_20_80 > 0.5),  mean(reg_preds_5_95  > 0.5)) * 100,
    "Baseline random", "Regrowth pixels"
  ),
  make_diag3(
    c(mean(reg_preds_ce),            mean(reg_preds_paired_20_80),        mean(reg_preds_paired_5_95)),
    c(mean(reg_preds_ce > 0.5),      mean(reg_preds_paired_20_80 > 0.5),  mean(reg_preds_paired_5_95  > 0.5)) * 100,
    "Williams paired", "Regrowth pixels"
  ),
  make_diag3(
    c(mean(reg_preds_ext_paired),        mean(reg_preds_ext_paired_20_80),        mean(reg_preds_ext_paired_5_95)),
    c(mean(reg_preds_ext_paired > 0.5),  mean(reg_preds_ext_paired_20_80 > 0.5),  mean(reg_preds_ext_paired_5_95  > 0.5)) * 100,
    "Extended paired", "Regrowth pixels"
  ),
  # ---- Non-regrowth pixels ----
  make_diag3(
    c(mean(preds_50_50),        mean(preds_20_80),        mean(preds_5_95)),
    c(mean(preds_50_50 > 0.5),  mean(preds_20_80 > 0.5),  mean(preds_5_95  > 0.5)) * 100,
    "Baseline random", "Non-regrowth pixels"
  ),
  make_diag3(
    c(mean(preds_ce),            mean(preds_paired_20_80),        mean(preds_paired_5_95)),
    c(mean(preds_ce > 0.5),      mean(preds_paired_20_80 > 0.5),  mean(preds_paired_5_95  > 0.5)) * 100,
    "Williams paired", "Non-regrowth pixels"
  ),
  make_diag3(
    c(mean(preds_ext_paired),        mean(preds_ext_paired_20_80),        mean(preds_ext_paired_5_95)),
    c(mean(preds_ext_paired > 0.5),  mean(preds_ext_paired_20_80 > 0.5),  mean(preds_ext_paired_5_95  > 0.5)) * 100,
    "Extended paired", "Non-regrowth pixels"
  )
))

diag3[, prev_label := factor(prev_label, levels = prev_labels3)]
diag3[, design     := factor(design,
  levels = c("Baseline random", "Williams paired", "Extended paired"))]
diag3[, pixel_set  := factor(pixel_set,
  levels = c("Regrowth pixels", "Non-regrowth pixels"))]

cat("\nDiagnostics v2 (3 designs x 2 pixel sets):\n")
print(diag3)

make_diag3_panel <- function(metric, ylab, title_str, add_ref = FALSE) {
  p <- ggplot(diag3,
              aes(x = prev_label, y = .data[[metric]],
                  colour = design, shape = design,
                  linetype = design, group = design)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = design_colours, name = NULL) +
    scale_shape_manual(values  = design_shapes,  name = NULL) +
    scale_linetype_manual(values = design_ltys,  name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    facet_wrap(~ pixel_set) +
    labs(title = title_str, x = "Prevalence", y = ylab) +
    nature_theme +
    theme(legend.position  = "bottom",
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text       = element_text(size = 8, face = "bold"))
  if (add_ref)
    p <- p + geom_hline(yintercept = 50, linetype = "dashed",
                        colour = "grey60", linewidth = 0.4)
  p
}

p_diag3_prob <- make_diag3_panel("mean_prob",
                                 "Mean predicted P(regrowth)",
                                 "a  Mean predicted P(regrowth)")
p_diag3_pct  <- make_diag3_panel("pct_above",
                                 "% classified as regrowth",
                                 "b  % classified as regrowth (threshold = 0.5)",
                                 add_ref = TRUE)

fig_supp_diag3 <- p_diag3_prob / p_diag3_pct +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_supp_diag3.pdf"),
       fig_supp_diag3, device = cairo_pdf, width = 180, height = 160, units = "mm")
ggsave(file.path(output_dir, "fig_supp_diag3.png"),
       fig_supp_diag3, width = 180, height = 160, units = "mm", dpi = 300)
message("Saved: fig_supp_diag3.pdf / .png")

# LaTeX caption for fig_supp_diag3:
#   Predicted P(regrowth) on regrowth (left) and non-regrowth (right)
#   pixels as a function of training prevalence under three designs:
#   baseline random sampling (red), Williams-variable paired (blue),
#   and extended-variable paired (green). (a) Mean predicted
#   probability. (b) Percentage classified as regrowth at threshold
#   0.5. For regrowth pixels, higher prevalence inflates predicted
#   probability regardless of design, but the paired designs show
#   attenuated sensitivity and higher absolute recovery rates at
#   reduced prevalence, consistent with the climate-envelope control
#   reducing between-biome suppression. For non-regrowth pixels, the
#   mirror pattern (false positive rate) is similarly attenuated under
#   paired designs.
