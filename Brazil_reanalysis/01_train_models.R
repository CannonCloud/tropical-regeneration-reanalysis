# ============================================================
# 01_train_models.R
# Brazil Natural Regeneration — Williams et al. Replication
#
# Trains four RF models and saves RDS files:
#   rf_prev_50_50.rds
#   rf_prev_20_80.rds
#   rf_prev_5_95.rds
#   rf_ce_50_40_10_nn.rds
#   pca_fit.rds
#   pred_sample.rds
#
# Run time: ~2-3 hours at full sample
# RAM: ~16GB peak during spatial matching (R6)
# ============================================================

library(data.table)
library(ranger)
options(ranger.num.threads = parallel::detectCores() - 3)
library(sf)
library(nngeo)

setwd("/home/cannon/Dropbox/Regrowth/pnv_sampling/brazil_predict")

output_dir <- "rf_outputs"
dir.create(output_dir, showWarnings = FALSE)

# ---- R1: Load & clean --------------------------------------

message("Loading data...")
dt <- fread("gee_outputs/brazil_training_data.csv")
message(sprintf("Loaded: %d rows x %d cols", nrow(dt), ncol(dt)))

dt[, bio_forest_density_2018 := NULL]
dt[, bio_landcover_class     := as.factor(bio_landcover_class)]
dt[, bio_biome_id            := as.factor(bio_biome_id)]
dt[, y                       := as.factor(y)]

print(dt[, .N, by = sample_set])

# ---- R2: Bioclim PCA ---------------------------------------

bio_cols <- paste0("bio", sprintf("%02d", 1:19))

message("Fitting PCA on regrowth + nonregrowth...")
pca_data <- dt[sample_set %in% c("regrowth", "nonregrowth"), ..bio_cols]
pca_fit  <- prcomp(pca_data, center = TRUE, scale. = TRUE)

cat("\nVariance explained by PC1-PC5:\n")
print(summary(pca_fit)$importance[, 1:5])

pcs <- predict(pca_fit, dt[, ..bio_cols])
dt[, bioclim_pc1 := pcs[, 1]]
dt[, bioclim_pc2 := pcs[, 2]]
dt[, bioclim_pc3 := pcs[, 3]]
dt[, bioclim_pc4 := pcs[, 4]]
dt[, bioclim_pc5 := pcs[, 5]]

rm(pca_data, pcs); gc()

williams_vars <- c(
  "bio_forest_density_1km2",
  "bio_dist_forest_2000",
  "bio_landcover_class",
  "bio_soil_organic_carbon",
  "bio_soil_ph",
  "bio_biome_id",
  "bioclim_pc1", "bioclim_pc2", "bioclim_pc3", "bioclim_pc4"
)

dt <- dt[complete.cases(dt[, ..williams_vars])]
message(sprintf("Rows after NA drop: %d", nrow(dt)))

# Save PCA fit for use in 02_figures.R
saveRDS(pca_fit,       file.path(output_dir, "pca_fit.rds"))
saveRDS(williams_vars, file.path(output_dir, "williams_vars.rds"))
message("Saved: pca_fit.rds, williams_vars.rds")

# ---- R3: Prevalence sensitivity loop -----------------------
# n_total fixed across prevalence ratios so only the ratio changes.
# 50/50 = Williams' unjustified choice; 5/95 ~ true landscape prevalence.

set.seed(42)
prevalences    <- c(0.50, 0.20, 0.05)
n_regrowth_all <- dt[sample_set == "regrowth", .N]
n_total        <- 2 * n_regrowth_all   # ~398k; set to 20000 for pilot

message(sprintf("\nn_total = %d (n_regrowth = %d)", n_total, n_regrowth_all))

for (prev in prevalences) {

  label <- sprintf("prev_%.0f_%.0f", prev * 100, (1 - prev) * 100)
  message(sprintf("\n--- Fitting RF: %s ---", label))

  n_pos <- min(round(n_total * prev),       dt[sample_set == "regrowth",    .N])
  n_neg <- min(n_total - round(n_total * prev), dt[sample_set == "nonregrowth", .N])

  train_pos <- dt[sample_set == "regrowth"]    |> _[sample(.N, n_pos)]
  train_neg <- dt[sample_set == "nonregrowth"] |> _[sample(.N, n_neg)]
  train_dt  <- rbind(train_pos, train_neg)
  rm(train_pos, train_neg); gc()

  message(sprintf("  Training: %d pos, %d neg (%.1f%% positive)",
                  n_pos, n_neg, 100 * n_pos / (n_pos + n_neg)))

  t0 <- proc.time()
  rf <- ranger(
    formula     = y ~ .,
    data        = train_dt[, c("y", williams_vars), with = FALSE],
    num.trees   = 500,
    mtry        = floor(sqrt(length(williams_vars))),
    probability = TRUE,
    importance  = "permutation",
    seed        = 42
  )
  message(sprintf("  Fitted in %.1f seconds", (proc.time() - t0)["elapsed"]))

  saveRDS(rf, file.path(output_dir, sprintf("rf_%s.rds", label)))
  message(sprintf("  Saved: rf_%s.rds", label))

  rm(rf, train_dt); gc()
}

# ---- R4: Build pred_sample (fixed prediction surface) ------
# Used in 02_figures.R for area estimates.
# Random sample of nonregrowth domain — large enough that mean pred is stable.

total_domain_area_ha <- 169944473.57   # from GEE pixelArea at 60m
n_pred_sample        <- 797076

set.seed(99)
pred_sample <- dt[sample_set == "nonregrowth"][sample(.N, n_pred_sample)]
saveRDS(pred_sample,          file.path(output_dir, "pred_sample.rds"))
saveRDS(total_domain_area_ha, file.path(output_dir, "total_domain_area_ha.rds"))
message("Saved: pred_sample.rds, total_domain_area_ha.rds")
# ---- R5: Spatial matching — multi-pass greedy NN -----------
# Runs st_nn once with k = K_NN candidates per regrowth point,
# then extracts up to 3 greedy passes from the result.
# Each pass assigns each regrowth point its next-closest unassigned
# negative, giving three distinct matched negative pools that are
# used across the paired prevalence models in R5b and R7.
#
# Saved RDS files:
#   negatives_paired_p1.rds  — 1st-nearest match per regrowth point
#   negatives_paired_p2.rds  — 2nd-nearest match per regrowth point
#   negatives_paired_p3.rds  — 3rd-nearest match per regrowth point
#
# Tuning: reduce K_NN if st_nn crashes due to RAM.
# K_NN = 3000 is the first attempt; try 2000, 1500, 1000 if needed.
# With k=1000 only pass 1 is reliable; passes 2-3 will have more
# unmatched points but the code handles this gracefully.

K_NN <- 2000

message(sprintf("\n--- Spatial matching: greedy NN, k=%d ---", K_NN))

candidate_pool_dt <- dt[sample_set %in% c("annulus", "nonregrowth")]
candidate_pool_sf <- st_as_sf(candidate_pool_dt,
                               coords = c("longitude", "latitude"), crs = 4326) |>
                     st_transform(5641)

regrowth_sf <- st_as_sf(dt[sample_set == "regrowth"],
                        coords = c("longitude", "latitude"), crs = 4326)
regrowth_m  <- st_transform(regrowth_sf, 5641)

message(sprintf("  Running st_nn k=%d (this is the slow step)...", K_NN))
t0        <- proc.time()
nn_result <- st_nn(regrowth_m, candidate_pool_sf,
                   k = K_NN, returnDist = TRUE, progress = TRUE)
message(sprintf("  Done in %.1f seconds", (proc.time() - t0)["elapsed"]))

rm(candidate_pool_sf, regrowth_sf, regrowth_m); gc()

# Save nn_result so R7 can be re-run without repeating st_nn
saveRDS(nn_result, file.path(output_dir, "nn_result.rds"))
message("  Saved: nn_result.rds")

# ---- Greedy pass function -----------------------------------
# Takes nn_result and a logical vector of already-used negative
# indices, assigns each regrowth point its closest available
# negative, returns (assigned indices, updated used_negs).

greedy_pass <- function(nn_result, used_negs, pass_label) {
  n_regrowth <- length(nn_result$nn)
  assigned   <- integer(n_regrowth)
  unmatched  <- 0L

  for (i in seq_len(n_regrowth)) {
    matched <- FALSE
    for (neg_idx in nn_result$nn[[i]]) {
      if (!used_negs[neg_idx]) {
        assigned[i]        <- neg_idx
        used_negs[neg_idx] <- TRUE
        matched            <- TRUE
        break
      }
    }
    if (!matched) {
      assigned[i] <- NA_integer_
      unmatched   <- unmatched + 1L
    }
  }

  cat(sprintf("  Pass %s — assigned: %d / %d (%.1f%%)  unmatched: %d\n",
              pass_label,
              sum(!is.na(assigned)), n_regrowth,
              100 * mean(!is.na(assigned)),
              unmatched))

  list(assigned = assigned, used_negs = used_negs)
}

# Three passes — each uses negatives not yet claimed by earlier passes
used_negs <- logical(nrow(candidate_pool_dt))

res1      <- greedy_pass(nn_result, used_negs, "1")
used_negs <- res1$used_negs

res2      <- greedy_pass(nn_result, used_negs, "2")
used_negs <- res2$used_negs

res3      <- greedy_pass(nn_result, used_negs, "3")
used_negs <- res3$used_negs

res4      <- greedy_pass(nn_result, used_negs, "4")
used_negs <- res4$used_negs

res5      <- greedy_pass(nn_result, used_negs, "5")
used_negs <- res5$used_negs

res6      <- greedy_pass(nn_result, used_negs, "6")
used_negs <- res6$used_negs

res7      <- greedy_pass(nn_result, used_negs, "7")

rm(nn_result, used_negs); gc()

# Extract negative rows for each pass
negatives_paired_p1 <- candidate_pool_dt[res1$assigned[!is.na(res1$assigned)]]
negatives_paired_p2 <- candidate_pool_dt[res2$assigned[!is.na(res2$assigned)]]
negatives_paired_p3 <- candidate_pool_dt[res3$assigned[!is.na(res3$assigned)]]
negatives_paired_p4 <- candidate_pool_dt[res4$assigned[!is.na(res4$assigned)]]
negatives_paired_p5 <- candidate_pool_dt[res5$assigned[!is.na(res5$assigned)]]
negatives_paired_p6 <- candidate_pool_dt[res6$assigned[!is.na(res6$assigned)]]
negatives_paired_p7 <- candidate_pool_dt[res7$assigned[!is.na(res7$assigned)]]
rm(res1, res2, res3, res4, res5, res6, res7, candidate_pool_dt); gc()

# Report source composition for each pass
for (nm in c("negatives_paired_p1", "negatives_paired_p2", "negatives_paired_p3", "negatives_paired_p4", "negatives_paired_p5", "negatives_paired_p6", "negatives_paired_p7")) {
  neg <- get(nm)
  cat(sprintf("  %s: n=%d  annulus=%.1f%%  nonregrowth=%.1f%%\n",
              nm, nrow(neg),
              100 * mean(neg$sample_set == "annulus"),
              100 * mean(neg$sample_set == "nonregrowth")))
}

saveRDS(negatives_paired_p1, file.path(output_dir, "negatives_paired_p1.rds"))
saveRDS(negatives_paired_p2, file.path(output_dir, "negatives_paired_p2.rds"))
saveRDS(negatives_paired_p3, file.path(output_dir, "negatives_paired_p3.rds"))
saveRDS(negatives_paired_p4, file.path(output_dir, "negatives_paired_p4.rds"))
saveRDS(negatives_paired_p5, file.path(output_dir, "negatives_paired_p5.rds"))
saveRDS(negatives_paired_p6, file.path(output_dir, "negatives_paired_p6.rds"))
saveRDS(negatives_paired_p7, file.path(output_dir, "negatives_paired_p7.rds"))
message("  Saved: negatives_paired_p1-p7.rds")

# Keep p1 in memory as negatives_paired for R5b and R6 compatibility
negatives_paired <- negatives_paired_p1

# ---- R5b: Fit CE model (50/40/10, pass-1 matched) ----------

message("\n--- Fitting RF: CE paired (50/40/10, greedy NN pass 1) ---")

set.seed(42)
n_ce      <- n_total
n_pos_ce  <- min(round(n_ce * 0.50), dt[sample_set == "regrowth",    .N])
n_pair_ce <- min(round(n_ce * 0.40), nrow(negatives_paired))
n_wide_ce <- min(n_ce - round(n_ce * 0.50) - round(n_ce * 0.40),
                 dt[sample_set == "nonregrowth", .N])

message(sprintf("  Available matched negatives: %d (need %d for 40%%)",
                nrow(negatives_paired), round(n_ce * 0.40)))

train_pos_ce  <- dt[sample_set == "regrowth"]    |> _[sample(.N, n_pos_ce)]
train_pair_ce <- negatives_paired[sample(.N, n_pair_ce)]
train_wide_ce <- dt[sample_set == "nonregrowth"] |> _[sample(.N, n_wide_ce)]
train_ce      <- rbind(train_pos_ce, train_pair_ce, train_wide_ce)
rm(train_pos_ce, train_pair_ce, train_wide_ce); gc()

message(sprintf("  Training: %d regrowth, %d matched, %d wide",
                n_pos_ce, n_pair_ce, n_wide_ce))
message(sprintf("  Total: %d rows (%.1f%% positive)",
                nrow(train_ce), 100 * n_pos_ce / nrow(train_ce)))

t0 <- proc.time()
rf_ce <- ranger(
  formula     = y ~ .,
  data        = train_ce[, c("y", williams_vars), with = FALSE],
  num.trees   = 500,
  mtry        = floor(sqrt(length(williams_vars))),
  probability = TRUE,
  importance  = "permutation",
  seed        = 42
)
message(sprintf("  Fitted in %.1f seconds", (proc.time() - t0)["elapsed"]))

saveRDS(rf_ce, file.path(output_dir, "rf_ce_50_40_10_nn.rds"))
message("  Saved: rf_ce_50_40_10_nn.rds")
rm(train_ce, rf_ce); gc()

message("\n01_train_models.R complete. All RDS files saved to rf_outputs/")


# ---- R6: Extended biophysical RF (baseline and paired) -----
# Williams et al. use only 10 predictors. Here we add the remaining
# biophysical variables available in the dataset to test whether
# within-climate variables (NPP, fire, slope, elevation, dist to water)
# gain relative importance under the paired design.
# Models saved as:
#   rf_ext_50_50.rds  — extended baseline (50/50)
#   rf_ext_paired.rds — extended paired (50/40/10, greedy NN pass 1)

extended_vars <- c(
  williams_vars,
  "bio_npp_mean",
  "bio_fire_freq",
  "bio_slope",
  "bio_elevation",
  #"paper_dist_water", # never downloaded it
  "bioclim_pc5",
  "bio_soil_sand",
  "bio_soil_clay",
  "bio_soil_bulk_density_fine",
  "bio_soil_water_33kpa"
)

# Verify all vars exist in dt
missing_ext <- extended_vars[!extended_vars %in% names(dt)]
if (length(missing_ext) > 0) {
  warning("Missing extended vars: ", paste(missing_ext, collapse = ", "))
  extended_vars <- extended_vars[extended_vars %in% names(dt)]
}

# Drop NAs for extended predictor set
dt_ext <- dt[complete.cases(dt[, ..extended_vars])]
message(sprintf("Rows after NA drop for extended vars: %d", nrow(dt_ext)))

# ---- R6a: Extended baseline (50/50) ------------------------

message("\n--- Fitting extended RF: 50/50 baseline ---")

set.seed(42)
n_pos_ext <- min(round(n_total * 0.50), dt_ext[sample_set == "regrowth",    .N])
n_neg_ext <- min(n_total - round(n_total * 0.50), dt_ext[sample_set == "nonregrowth", .N])

train_pos_ext <- dt_ext[sample_set == "regrowth"]    |> _[sample(.N, n_pos_ext)]
train_neg_ext <- dt_ext[sample_set == "nonregrowth"] |> _[sample(.N, n_neg_ext)]
train_ext     <- rbind(train_pos_ext, train_neg_ext)
rm(train_pos_ext, train_neg_ext); gc()

t0 <- proc.time()
rf_ext_50_50 <- ranger(
  formula     = y ~ .,
  data        = train_ext[, c("y", extended_vars), with = FALSE],
  num.trees   = 500,
  mtry        = floor(sqrt(length(extended_vars))),
  probability = TRUE,
  importance  = "permutation",
  seed        = 42
)
message(sprintf("  Fitted in %.1f seconds", (proc.time() - t0)["elapsed"]))
rm(train_ext); gc()

saveRDS(rf_ext_50_50, file.path(output_dir, "rf_ext_50_50.rds"))
saveRDS(extended_vars, file.path(output_dir, "extended_vars.rds"))
message("  Saved: rf_ext_50_50.rds")
rm(rf_ext_50_50); gc()

# ---- R6b: Extended paired (50/40/10, pass 1) ---------------

message("\n--- Fitting extended RF: paired (50/40/10) ---")

negatives_paired_ext <- negatives_paired[
  complete.cases(negatives_paired[, ..extended_vars])
]

set.seed(42)
n_pos_ep  <- min(round(n_total * 0.50), dt_ext[sample_set == "regrowth",    .N])
n_pair_ep <- min(round(n_total * 0.40), nrow(negatives_paired_ext))
n_wide_ep <- min(n_total - round(n_total * 0.50) - round(n_total * 0.40),
                 dt_ext[sample_set == "nonregrowth", .N])

train_pos_ep  <- dt_ext[sample_set == "regrowth"]    |> _[sample(.N, n_pos_ep)]
train_pair_ep <- negatives_paired_ext[sample(.N, n_pair_ep)]
train_wide_ep <- dt_ext[sample_set == "nonregrowth"] |> _[sample(.N, n_wide_ep)]
train_ep      <- rbind(train_pos_ep, train_pair_ep, train_wide_ep)
rm(train_pos_ep, train_pair_ep, train_wide_ep); gc()

message(sprintf("  Training: %d regrowth, %d matched, %d wide",
                n_pos_ep, n_pair_ep, n_wide_ep))

t0 <- proc.time()
rf_ext_paired <- ranger(
  formula     = y ~ .,
  data        = train_ep[, c("y", extended_vars), with = FALSE],
  num.trees   = 500,
  mtry        = floor(sqrt(length(extended_vars))),
  probability = TRUE,
  importance  = "permutation",
  seed        = 42
)
message(sprintf("  Fitted in %.1f seconds", (proc.time() - t0)["elapsed"]))
rm(train_ep); gc()

saveRDS(rf_ext_paired, file.path(output_dir, "rf_ext_paired.rds"))
message("  Saved: rf_ext_paired.rds")
rm(rf_ext_paired); gc()

message("\nR6 complete. Extended models saved.")

# ---- R7: Paired prevalence sensitivity (20/80 and 5/95) ----
# Uses multi-pass matched negatives from R5 to train paired models
# at reduced positive prevalence. Negative pools:
#
#   20/80: 20% regrowth + 40% p1 + 30% p2 + 10% wide
#   5/95:  5% regrowth  + 40% p1 + 22% p2 + 14% p3 + 10% p4
#                       +  8% p5 +  4% p6 +  2% p7 + 0% wide
#          (wide=0 because paired negatives cover the full 95%)
#
# The positive subsample is always drawn randomly from the full
# ~199k regrowth set, so spatial coverage stays consistent.
#
# Standalone execution requires:
#   - dt and dt_ext in memory (or re-source R1-R2 + R6 setup)
#   - negatives_paired_p1 through p7 in rf_outputs/
#   - extended_vars.rds and williams_vars.rds in rf_outputs/
#
# Models saved as:
#   rf_paired_20_80.rds
#   rf_paired_5_95.rds
#   rf_ext_paired_20_80.rds
#   rf_ext_paired_5_95.rds

message("\n\n============================================================")
message("R7: Paired prevalence sensitivity (20/80 and 5/95)")
message("============================================================")

# Guards for standalone execution
if (!exists("dt")) {
  stop("dt not in memory. Source R1-R2 first or load from CSV.")
}

if (!exists("negatives_paired_p1")) {
  message("Loading matched negative pools from disk...")
  for (i in 1:7) {
    assign(sprintf("negatives_paired_p%d", i),
           readRDS(file.path(output_dir, sprintf("negatives_paired_p%d.rds", i))))
  }
  message("  Loaded p1-p7.")
}

if (!exists("extended_vars")) {
  extended_vars <- readRDS(file.path(output_dir, "extended_vars.rds"))
}

if (!exists("dt_ext")) {
  message("Rebuilding dt_ext...")
  dt_ext <- dt[complete.cases(dt[, ..extended_vars])]
  message(sprintf("  Rows: %d", nrow(dt_ext)))
}

if (!exists("n_total")) {
  n_total <- 2 * dt[sample_set == "regrowth", .N]
  message(sprintf("  n_total set to %d", n_total))
}

if (!exists("williams_vars")) {
  williams_vars <- readRDS(file.path(output_dir, "williams_vars.rds"))
}

# ---- Design table ------------------------------------------

designs <- list(
  list(label    = "paired_20_80",     rds = "rf_paired_20_80.rds",
       vars     = "williams",         pos_frac = 0.20),
  list(label    = "paired_5_95",      rds = "rf_paired_5_95.rds",
       vars     = "williams",         pos_frac = 0.05),
  list(label    = "ext_paired_20_80", rds = "rf_ext_paired_20_80.rds",
       vars     = "extended",         pos_frac = 0.20),
  list(label    = "ext_paired_5_95",  rds = "rf_ext_paired_5_95.rds",
       vars     = "extended",         pos_frac = 0.05)
)

for (d in designs) {
  
  use_vars <- if (d$vars == "williams") williams_vars else extended_vars
  use_dt   <- if (d$vars == "williams") dt            else dt_ext
  
  # Pre-filter all pools once for this var set
  pools <- lapply(1:7, function(i) {
    p <- get(sprintf("negatives_paired_p%d", i))
    p[complete.cases(p[, ..use_vars])]
  })
  
  # Fixed components
  n_pos  <- min(round(n_total * d$pos_frac), use_dt[sample_set == "regrowth",    .N])
  n_wide <- min(round(n_total * 0.10),       use_dt[sample_set == "nonregrowth", .N])
  
  # Greedy fill of remaining budget from p1 → p7 in order
  budget    <- n_total - n_pos - n_wide
  n_pools   <- integer(7)
  remaining <- budget
  
  for (i in 1:7) {
    if (remaining <= 0) break
    n_pools[i] <- min(remaining, nrow(pools[[i]]))
    remaining  <- remaining - n_pools[i]
  }
  
  message(sprintf("\n--- Fitting RF: %s ---", d$label))
  message(sprintf("  pos=%d  %s  wide=%d  total=%d  budget_remaining=%d",
                  n_pos,
                  paste(sprintf("p%d=%d", 1:7, n_pools), collapse = "  "),
                  n_wide,
                  n_pos + sum(n_pools) + n_wide,
                  remaining))
  message(sprintf("  %.1f%% positive",
                  100 * n_pos / (n_pos + sum(n_pools) + n_wide)))
  
  # Build training set
  set.seed(42)
  train_parts <- list(
    use_dt[sample_set == "regrowth"]    |> _[sample(.N, n_pos)],
    use_dt[sample_set == "nonregrowth"] |> _[sample(.N, n_wide)]
  )
  for (i in 1:7) {
    if (n_pools[i] > 0)
      train_parts[[length(train_parts) + 1]] <- pools[[i]][sample(.N, n_pools[i])]
  }
  rm(pools); gc()
  
  train_dt <- rbindlist(train_parts)
  rm(train_parts); gc()
  
  t0 <- proc.time()
  rf_tmp <- ranger(
    formula     = y ~ .,
    data        = train_dt[, c("y", use_vars), with = FALSE],
    num.trees   = 500,
    mtry        = floor(sqrt(length(use_vars))),
    probability = TRUE,
    importance  = "permutation",
    seed        = 42
  )
  message(sprintf("  Fitted in %.1f seconds", (proc.time() - t0)["elapsed"]))
  rm(train_dt); gc()
  
  saveRDS(rf_tmp, file.path(output_dir, d$rds))
  message(sprintf("  Saved: %s", d$rds))
  rm(rf_tmp); gc()
}

message("\nR7 complete.")
message("Saved: rf_paired_20_80.rds, rf_paired_5_95.rds,")
message("       rf_ext_paired_20_80.rds, rf_ext_paired_5_95.rds")