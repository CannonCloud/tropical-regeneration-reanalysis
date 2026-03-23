# ============================================================
# oob_metrics_brazil.R
#
# Extracts OOB accuracy + AUC for Brazil baseline and paired models.
# Run from:  # set your working directory here
#
# Requires:
#   gee_outputs/brazil_training_data.csv
#   rf_outputs/rf_prev_50_50.rds
#   rf_outputs/rf_prev_20_80.rds
#   rf_outputs/rf_prev_5_95.rds
#   rf_outputs/rf_ce_50_40_10_nn.rds
#   rf_outputs/negatives_paired_p1.rds
#   rf_outputs/williams_vars.rds
#
# All models were trained with probability = TRUE, so rf$prediction.error
# is the OOB Brier score. OOB accuracy and AUC are computed here by
# reconstructing the training y vector with the same seeds used in
# 01_train_models.R and pairing it with rf$predictions.
# ============================================================

library(data.table)

output_dir <- "rf_outputs"

# ---- 1. Rebuild dt (mirrors R1 in 01_train_models.R) --------

message("Loading training data and projecting PCA...")
dt <- fread("gee_outputs/brazil_training_data.csv")
dt[, bio_forest_density_2018 := NULL]
dt[, bio_landcover_class     := as.factor(bio_landcover_class)]
dt[, bio_biome_id            := as.factor(bio_biome_id)]
dt[, y                       := as.factor(y)]

pca_fit  <- readRDS(file.path(output_dir, "pca_fit.rds"))
bio_cols <- paste0("bio", sprintf("%02d", 1:19))
pcs      <- predict(pca_fit, dt[, ..bio_cols])
dt[, bioclim_pc1 := pcs[, 1]]
dt[, bioclim_pc2 := pcs[, 2]]
dt[, bioclim_pc3 := pcs[, 3]]
dt[, bioclim_pc4 := pcs[, 4]]
dt[, bioclim_pc5 := pcs[, 5]]
rm(pcs); gc()

williams_vars <- readRDS(file.path(output_dir, "williams_vars.rds"))
dt <- dt[complete.cases(dt[, ..williams_vars])]
message(sprintf("dt: %d rows after NA drop", nrow(dt)))

safe_rbind <- function(...) {
  parts <- list(...)
  common <- Reduce(intersect, lapply(parts, names))
  rbindlist(lapply(parts, function(x) x[, ..common]))
}

n_regrowth_all <- dt[sample_set == "regrowth", .N]
n_total        <- 2 * n_regrowth_all

# ---- 2. Rank-based AUC (no extra packages) ------------------

auc_roc <- function(y_true_01, prob_pos) {
  pos <- prob_pos[y_true_01 == 1L]
  neg <- prob_pos[y_true_01 == 0L]
  n1  <- as.numeric(length(pos)); n0 <- as.numeric(length(neg))
  r   <- rank(c(pos, neg))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- 3. Metric extractor ------------------------------------

get_oob_metrics <- function(rds_name, y_train, label) {
  rf       <- readRDS(file.path(output_dir, rds_name))
  # rf$predictions: n_train x 2 matrix; columns named by factor levels ("0","1")
  oob_prob <- rf$predictions[, "1"]          # P(regrowth) per training obs
  y_01     <- as.integer(y_train == "1")
  
  acc   <- mean((oob_prob > 0.5) == y_01)
  auc   <- auc_roc(y_01, oob_prob)
  brier <- rf$prediction.error               # stored directly in ranger object
  
  cat(sprintf("%-28s  Acc=%.4f  AUC=%.4f  Brier=%.5f  n=%d\n",
              label, acc, auc, brier, length(y_01)))
  
  invisible(list(label = label, accuracy = acc, auc = auc, brier = brier,
                 n = length(y_01)))
}

# ---- 4. Prevalence models -----------------------------------
#
# IMPORTANT: 01_train_models.R calls set.seed(42) ONCE before the
# prevalence loop (line 83). The three sample() draws share the same
# sequential random state. We replicate that exactly here.

message("\n--- Prevalence models (R3) ---")

results <- list()
set.seed(42)   # single seed before the loop, same as 01_train_models.R

for (prev in c(0.50, 0.20, 0.05)) {
  n_pos <- min(round(n_total * prev),           dt[sample_set == "regrowth",    .N])
  n_neg <- min(n_total - round(n_total * prev), dt[sample_set == "nonregrowth", .N])
  
  train_pos <- dt[sample_set == "regrowth"]    [sample(.N, n_pos)]
  train_neg <- dt[sample_set == "nonregrowth"] [sample(.N, n_neg)]
  train_dt  <- rbind(train_pos, train_neg)
  rm(train_pos, train_neg)
  
  lbl <- sprintf("prev_%.0f_%.0f", prev * 100, (1 - prev) * 100)
  rds <- sprintf("rf_%s.rds", lbl)
  
  results[[lbl]] <- get_oob_metrics(rds, train_dt$y, lbl)
  rm(train_dt)
}

# ---- 5. CE paired model (50/40/10, pass-1 matched) ----------
#
# R5b in 01_train_models.R uses its own set.seed(42) at line 273,
# independent of the prevalence loop above.

message("\n--- CE paired model (R5b) ---")

negatives_paired <- readRDS(file.path(output_dir, "negatives_paired_p1.rds"))

set.seed(42)   # mirrors line 273 in 01_train_models.R
n_pos_ce  <- min(round(n_total * 0.50), dt[sample_set == "regrowth",    .N])
n_pair_ce <- min(round(n_total * 0.40), nrow(negatives_paired))
n_wide_ce <- min(n_total - round(n_total * 0.50) - round(n_total * 0.40),
                 dt[sample_set == "nonregrowth", .N])

train_pos_ce  <- dt[sample_set == "regrowth"]    [sample(.N, n_pos_ce)]
train_pair_ce <- negatives_paired[sample(.N, n_pair_ce)]
train_wide_ce <- dt[sample_set == "nonregrowth"] [sample(.N, n_wide_ce)]
train_ce      <- safe_rbind(train_pos_ce, train_pair_ce, train_wide_ce)
rm(train_pos_ce, train_pair_ce, train_wide_ce)

results[["ce_paired_50_40_10"]] <- get_oob_metrics(
  "rf_ce_50_40_10_nn.rds", train_ce$y, "CE paired (50/40/10)"
)
rm(train_ce, negatives_paired)

# ---- 6. Paired prevalence sensitivity (20/80 and 5/95) ------
#
# R7 in 01_train_models.R uses set.seed(42) INSIDE each loop iteration,
# so every design gets the same seed independently.
# Negative pools are filled greedily p1 → p7.

message("\n--- Paired prevalence sensitivity models (R7) ---")

# Load all matched negative pools
pools <- lapply(1:7, function(i)
  readRDS(file.path(output_dir, sprintf("negatives_paired_p%d.rds", i))))

for (cfg in list(
  list(pos_frac = 0.20, rds = "rf_paired_20_80.rds", lbl = "paired_20_80"),
  list(pos_frac = 0.05, rds = "rf_paired_5_95.rds",  lbl = "paired_5_95")
)) {
  
  n_pos  <- min(round(n_total * cfg$pos_frac), dt[sample_set == "regrowth",    .N])
  n_wide <- min(round(n_total * 0.10),         dt[sample_set == "nonregrowth", .N])
  
  # Greedy fill remaining budget from p1 → p7
  budget    <- n_total - n_pos - n_wide
  n_pools   <- integer(7)
  remaining <- budget
  for (i in 1:7) {
    if (remaining <= 0) break
    n_pools[i] <- min(remaining, nrow(pools[[i]]))
    remaining  <- remaining - n_pools[i]
  }
  
  set.seed(42)   # matches set.seed(42) inside R7 loop
  train_parts <- list(
    dt[sample_set == "regrowth"]    [sample(.N, n_pos)],
    dt[sample_set == "nonregrowth"] [sample(.N, n_wide)]
  )
  for (i in 1:7) {
    if (n_pools[i] > 0)
      train_parts[[length(train_parts) + 1]] <- pools[[i]][sample(.N, n_pools[i])]
  }
  train_dt <- do.call(safe_rbind, train_parts)
  rm(train_parts)
  
  message(sprintf("  %s: pos=%d  wide=%d  %s  total=%d",
                  cfg$lbl, n_pos, n_wide,
                  paste(sprintf("p%d=%d", 1:7, n_pools), collapse = " "),
                  nrow(train_dt)))
  
  results[[cfg$lbl]] <- get_oob_metrics(cfg$rds, train_dt$y, cfg$lbl)
  rm(train_dt)
}
rm(pools)

# ---- 7. Sensitivity / specificity from prediction caches ----
#
# 02_figures.R already saved per-model prediction caches:
#   cache$reg    — P(regrowth) for every regrowth pixel
#   cache$nonreg — P(regrowth) for 797k pred_sample nonregrowth pixels
#
# Sensitivity = mean(reg    > 0.5)   how often regrowth    is called regrowth
# Specificity = mean(nonreg <= 0.5)  how often non-regrowth is called non-regrowth
#
# No ranger objects loaded; no predict() calls needed.

message("\n--- Sensitivity / specificity from prediction caches ---")

cache_specs <- list(
  list(cache = "preds_rf_prev_50_50.rds",    lbl = "prev_50_50"),
  list(cache = "preds_rf_prev_20_80.rds",    lbl = "prev_20_80"),
  list(cache = "preds_rf_prev_5_95.rds",     lbl = "prev_5_95"),
  list(cache = "preds_rf_ce.rds",            lbl = "CE paired (50/40/10)"),
  list(cache = "preds_rf_paired_20_80.rds",  lbl = "paired_20_80"),
  list(cache = "preds_rf_paired_5_95.rds",   lbl = "paired_5_95")
)

full_results <- list()

for (m in cache_specs) {
  cache <- readRDS(file.path(output_dir, m$cache))
  
  sensitivity <- mean(cache$reg    > 0.5)
  specificity <- mean(cache$nonreg <= 0.5)
  
  cat(sprintf("%-28s  Sens=%.4f  Spec=%.4f  n_reg=%d  n_nonreg=%d\n",
              m$lbl, sensitivity, specificity,
              length(cache$reg), length(cache$nonreg)))
  
  full_results[[m$lbl]] <- list(
    label       = m$lbl,
    sensitivity = sensitivity,
    specificity = specificity
  )
}
# ---- 8. Summary tables --------------------------------------

cat("\n\n=== OOB Performance Summary (training-set) — Brazil Models ===\n")
cat(sprintf("%-28s  %8s  %8s  %8s  %8s\n",
            "Model", "Accuracy", "AUC", "Brier", "n_train"))
cat(strrep("-", 72), "\n")
for (r in results) {
  cat(sprintf("%-28s  %8.4f  %8.4f  %8.5f  %8d\n",
              r$label, r$accuracy, r$auc, r$brier, r$n))
}

cat("\n\n=== Full-sample Accuracy — Brazil Models ===\n")
cat(sprintf("  Sensitivity = P(pred regrowth     | true regrowth)    [from $reg cache]\n"))
cat(sprintf("  Specificity = P(pred non-regrowth | true non-regrowth) [from $nonreg cache, 797k sample]\n\n"))
cat(sprintf("%-28s  %11s  %11s\n", "Model", "Sensitivity", "Specificity"))
cat(strrep("-", 54), "\n")
for (r in full_results) {
  cat(sprintf("%-28s  %11.4f  %11.4f\n", r$label, r$sensitivity, r$specificity))
}

cat("\nNotes:\n")
cat("  OOB metrics use held-out (out-of-bag) training rows only.\n")
cat("  Sensitivity/specificity use pre-computed prediction caches from 02_figures.R:\n")
cat("    $reg    = all regrowth pixels (complete Williams vars)\n")
cat("    $nonreg = 797k pred_sample nonregrowth pixels\n")
cat("  Both are in-sample for the training rows but no ranger objects are reloaded.\n")
cat("  AUC is threshold-free and is the preferred discrimination metric.\n")

# ---- 9. LaTeX tabular output --------------------------------
# Writes just the \tabular{} block (no \begin{table}, no \caption)
# so it can be \input{} directly into the paper.
# ---- Recalibration function ----------------------------------
recalibrate <- function(p, pi_train, pi_new = 0.05) {
  beta <- (pi_new / (1 - pi_new)) / (pi_train / (1 - pi_train))
  (beta * p) / (beta * p + (1 - p))
}

full_domain_ha <- total_domain_area_ha + regrowth_domain_area_ha

# ---- Area computation per model ------------------------------
# Map each row_order key to its prediction vectors and training prevalence
area_specs <- list(
  prev_50_50         = list(nr = preds_50_50,            rg = reg_preds_50_50,            pi = 0.50),
  prev_20_80         = list(nr = preds_20_80,            rg = reg_preds_20_80,            pi = 0.20),
  prev_5_95          = list(nr = preds_5_95,             rg = reg_preds_5_95,             pi = 0.05),
  ce_paired_50_40_10 = list(nr = preds_ce,               rg = reg_preds_ce,               pi = 0.50),
  paired_20_80       = list(nr = preds_paired_20_80,     rg = reg_preds_paired_20_80,     pi = 0.20),
  paired_5_95        = list(nr = preds_paired_5_95,      rg = reg_preds_paired_5_95,      pi = 0.05)
)

compute_total_area <- function(preds_nr, preds_rg) {
  (mean(preds_nr) * total_domain_area_ha + mean(preds_rg) * regrowth_domain_area_ha) / 1e6
}

compute_adj_area <- function(preds_nr, preds_rg, pi_train) {
  (mean(recalibrate(preds_nr, pi_train)) * total_domain_area_ha + 
     mean(recalibrate(preds_rg, pi_train)) * regrowth_domain_area_ha) / 1e6
}

# ---- LaTeX table with area columns ---------------------------
format_spec <- function(x) {
  if (x > 0.9999) return("$\\sim$1.000")
  sprintf("%.3f", x)
}

row_order <- list(
  list(oob_key = "prev_50_50",         full_key = "prev_50_50",           label = "Baseline 50/50",      spacer = FALSE),
  list(oob_key = "prev_20_80",         full_key = "prev_20_80",           label = "Baseline 20/80",      spacer = FALSE),
  list(oob_key = "prev_5_95",          full_key = "prev_5_95",            label = "Baseline 5/95",       spacer = TRUE),
  list(oob_key = "ce_paired_50_40_10", full_key = "CE paired (50/40/10)", label = "CE paired 50/40/10",  spacer = FALSE),
  list(oob_key = "paired_20_80",       full_key = "paired_20_80",         label = "Paired 20/80",        spacer = FALSE),
  list(oob_key = "paired_5_95",        full_key = "paired_5_95",          label = "Paired 5/95",         spacer = FALSE)
)

tex_rows <- character(0)
for (rw in row_order) {
  o  <- results[[rw$oob_key]]
  f  <- full_results[[rw$full_key]]
  a  <- area_specs[[rw$oob_key]]
  
  raw_area <- compute_total_area(a$nr, a$rg)
  adj_area <- compute_adj_area(a$nr, a$rg, a$pi)
  
  tex_rows <- c(tex_rows, sprintf(
    "%-22s & %.3f & %.3f & %.3f & %.3f & %s & %.1f & %.1f \\\\",
    rw$label,
    o$accuracy, o$auc, o$brier,
    f$sensitivity, format_spec(f$specificity),
    raw_area, adj_area
  ))
  if (rw$spacer) tex_rows <- c(tex_rows, "\\addlinespace")
}

tex_tabular <- c(
  "\\begin{tabular}{@{}lccccccc@{}}",
  "\\toprule",
  "& \\multicolumn{3}{c}{OOB (training set)} & \\multicolumn{2}{c}{Full sample} & \\multicolumn{2}{c}{Area (Mha)} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}",
  "Model & Accuracy & AUC & Brier & Sensitivity & Specificity & Raw & Adj (5\\%) \\\\",
  "\\midrule",
  tex_rows,
  "\\bottomrule",
  "\\end{tabular}"
)

tex_path <- file.path(output_dir, "tab_performance_tabular.tex")
writeLines(tex_tabular, tex_path)
message(sprintf("Saved: %s", tex_path))

# ---- 10. Match quality diagnostics --------------------------
# Loads nn_result.rds and re-runs pass-1 greedy assignment, recording:
#   (a) nearest-neighbour distance (rank-1, before competition) — what
#       was probably cited as "1.2 km" in the manuscript
#   (b) assigned-match distance (the rank actually used after greedy
#       competition) — larger in dense areas
#
# Also reverse-engineers whether negatives_paired_p1.rds on disk was
# built from this nn_result (current k) or an earlier run with different k.

message("\n--- Match quality diagnostics ---")
message("  Loading nn_result.rds (this may take a moment)...")
nn_result  <- readRDS(file.path(output_dir, "nn_result.rds"))

n_regrowth <- length(nn_result$nn)
K_NN_used  <- length(nn_result$nn[[1]])
message(sprintf("  nn_result: n_regrowth=%d  k=%d", n_regrowth, K_NN_used))

# Nearest-neighbour distance (rank 1 for every point, ignoring competition)
nn1_dist_km <- sapply(nn_result$dist, function(d) d[1]) / 1000
cat(sprintf("\n  Nearest-neighbour dist (rank-1, no competition):\n"))
cat(sprintf("    median=%.2f km  IQR=[%.2f, %.2f]  p5=%.2f  p95=%.2f\n",
            median(nn1_dist_km),
            quantile(nn1_dist_km, 0.25), quantile(nn1_dist_km, 0.75),
            quantile(nn1_dist_km, 0.05), quantile(nn1_dist_km, 0.95)))

# Greedy pass 1 — assigned match distance (accounts for competition)
assigned_idx  <- integer(n_regrowth)
match_dist_m  <- numeric(n_regrowth)
assigned_rank <- integer(n_regrowth)   # which rank in the list was used
used_negs     <- logical(max(unlist(nn_result$nn)))

message("  Running greedy pass 1 with distance + rank recording...")
for (i in seq_len(n_regrowth)) {
  matched <- FALSE
  for (j in seq_along(nn_result$nn[[i]])) {
    neg_idx <- nn_result$nn[[i]][j]
    if (!used_negs[neg_idx]) {
      assigned_idx[i]    <- neg_idx
      match_dist_m[i]    <- nn_result$dist[[i]][j]
      assigned_rank[i]   <- j
      used_negs[neg_idx] <- TRUE
      matched            <- TRUE
      break
    }
  }
  if (!matched) {
    assigned_idx[i]  <- NA_integer_
    match_dist_m[i]  <- NA_real_
    assigned_rank[i] <- NA_integer_
  }
}
rm(nn_result, used_negs); gc()

matched_mask  <- !is.na(assigned_idx)
match_dist_km <- match_dist_m[matched_mask] / 1000
n_matched     <- sum(matched_mask)
n_unmatched   <- sum(!matched_mask)
match_rate    <- n_matched / n_regrowth

cat(sprintf("\n  Match rate (pass 1): %d / %d = %.1f%% matched  |  %.1f%% unmatched\n",
            n_matched, n_regrowth, 100 * match_rate, 100 * (1 - match_rate)))
cat(sprintf("  Assigned-match dist:  median=%.2f km  IQR=[%.2f, %.2f]  p5=%.2f  p95=%.2f\n",
            median(match_dist_km),
            quantile(match_dist_km, 0.25), quantile(match_dist_km, 0.75),
            quantile(match_dist_km, 0.05), quantile(match_dist_km, 0.95)))
cat(sprintf("  Assigned rank:        median=%.0f  IQR=[%.0f, %.0f]  p95=%.0f\n",
            median(assigned_rank[matched_mask]),
            quantile(assigned_rank[matched_mask], 0.25),
            quantile(assigned_rank[matched_mask], 0.75),
            quantile(assigned_rank[matched_mask], 0.95)))

# ---- Reverse-engineer k used to build negatives_paired_p1.rds ----
# If nrow(p1) == n_matched from this run, p1 came from this nn_result.
# If nrow(p1) is different, p1 was built from a run with a different k.
neg_p1_nrow <- nrow(readRDS(file.path(output_dir, "negatives_paired_p1.rds")))
cat(sprintf("\n  Reverse-engineering check:\n"))
cat(sprintf("    nrow(negatives_paired_p1.rds) on disk : %d\n", neg_p1_nrow))
cat(sprintf("    n_matched from current nn_result (k=%d): %d\n", K_NN_used, n_matched))
if (neg_p1_nrow == n_matched) {
  cat(sprintf("    --> MATCH: p1 was built from k=%d nn_result\n", K_NN_used))
} else {
  cat(sprintf("    --> MISMATCH: p1 was built from a DIFFERENT run (probably k=%d or k=%d)\n",
              ifelse(K_NN_used == 2000, 1000, 2000),
              ifelse(K_NN_used == 2000, 1500, 1000)))
  cat(sprintf("    The CE model used matched negatives from the earlier run.\n"))
  cat(sprintf("    This affects match distances but likely not covariate balance materially.\n"))
}

# How many matched candidates are available vs used in CE training
n_total_ce   <- n_regrowth * 2
n_pair_drawn <- min(round(n_total_ce * 0.40), n_matched)
cat(sprintf("\n  CE training usage of matched negatives:\n"))
cat(sprintf("    Available from pass 1 : %d (%.1f%% of regrowth)\n", n_matched, 100*match_rate))
cat(sprintf("    Drawn for 40%% budget  : %d (%.1f%% of available)\n",
            n_pair_drawn, 100 * n_pair_drawn / n_matched))

saveRDS(list(match_dist_km = match_dist_km,
             nn1_dist_km   = nn1_dist_km,
             assigned_rank = assigned_rank[matched_mask]),
        file.path(output_dir, "match_dist_km_pass1.rds"))
message("  Saved: match_dist_km_pass1.rds")

# Distance distribution plot — both nearest and assigned on same plot
png(file.path(output_dir, "fig_supp_match_dist.png"),
    width = 160, height = 90, units = "mm", res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
hist(nn1_dist_km, breaks = 100, col = "steelblue", border = "white",
     main = "Nearest-neighbour distance\n(rank 1, no competition)",
     xlab = "Distance (km)", ylab = "Count")
abline(v = median(nn1_dist_km), col = "firebrick", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Median=%.2f km", median(nn1_dist_km)),
       col = "firebrick", lty = 2, lwd = 2, bty = "n", cex = 0.75)
hist(match_dist_km, breaks = 100, col = "darkorange", border = "white",
     main = "Assigned-match distance\n(after greedy competition)",
     xlab = "Distance (km)", ylab = "Count")
abline(v = median(match_dist_km), col = "firebrick", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Median=%.2f km", median(match_dist_km)),
       col = "firebrick", lty = 2, lwd = 2, bty = "n", cex = 0.75)
dev.off()
message("  Saved: fig_supp_match_dist.png")

# ---- 11. Covariate balance ----------------------------------

message("\n--- Covariate balance: regrowth vs pass-1 matched negatives ---")

neg_p1     <- readRDS(file.path(output_dir, "negatives_paired_p1.rds"))
reg_pixels <- dt[sample_set == "regrowth"]

balance_vars <- c("bioclim_pc1", "bioclim_pc2",
                  "bio_soil_ph", "bio_soil_organic_carbon",
                  "bio_forest_density_1km2", "bio_dist_forest_2000")

cat(sprintf("\n  %-30s  %8s  %8s  %8s\n", "Variable", "Regrowth", "Matched", "Std.Diff"))
cat(strrep("-", 62), "\n")

balance_rows <- list()
for (v in balance_vars) {
  if (!v %in% names(reg_pixels) || !v %in% names(neg_p1)) {
    cat(sprintf("  %-30s  [not found in data]\n", v)); next
  }
  r_vals   <- reg_pixels[[v]][is.finite(reg_pixels[[v]])]
  m_vals   <- neg_p1[[v]][is.finite(neg_p1[[v]])]
  std_diff <- (mean(r_vals) - mean(m_vals)) / sqrt((var(r_vals) + var(m_vals)) / 2)
  cat(sprintf("  %-30s  %8.3f  %8.3f  %8.3f\n", v, mean(r_vals), mean(m_vals), std_diff))
  balance_rows[[v]] <- list(variable = v, mean_reg = mean(r_vals),
                            mean_neg = mean(m_vals), std_diff = std_diff)
}
cat("\n  |Std.Diff|: <0.1 good balance, 0.1-0.2 moderate, >0.2 poor\n")
cat("  Note: forest density/dist std.diff expected to be large — that IS the\n")
cat("  treatment contrast of interest, not a balance failure.\n")

# ---- 12. Corrected methods text -----------------------------
# Fill-in values computed above. The key corrections vs original text:
#   - k is K_NN_used (2000), not 1000
#   - 86.9% means that fraction had a match AVAILABLE in pass 1;
#     only min(40% of n_total, n_matched) were actually drawn for training
#   - "nearest" distance and "assigned" distance are distinguished

cat(sprintf("
=== CORRECTED METHODS TEXT ===

Spatial matching used a greedy nearest-neighbour algorithm implemented
via st_nn (nngeo package) in projected coordinates (SIRGAS 2000 /
Brazil Polyconic, EPSG:5641). For each of the %d regrowth points,
the %d nearest candidates from a combined pool of annulus and
Brazil-wide non-regrowth pixels were identified. A single greedy pass
then assigned each regrowth point its closest spatially available
candidate not already claimed by another regrowth point.

%.1f%%%% of regrowth points (n = %d) received a unique match within their
%d-candidate pool; the remaining %.1f%%%% (n = %d) were unmatched in
pass 1 because all %d of their nearest candidates had already been
claimed by spatially proximate regrowth points. The median distance
from each regrowth point to its nearest candidate (before competition)
was %.2f km (IQR: %.2f--%.2f km); the median distance to the
actually-assigned match (after competition) was %.2f km
(IQR: %.2f--%.2f km), reflecting displacement in densely clustered
areas such as the Amazon.

Of the %d matched negatives available from pass 1, a random subsample
of %d (%.1f%%%% of total training) was drawn to fill the 40%%%% negative
budget for the main CE model. The remaining 10%%%% of negatives were
drawn Brazil-wide to retain representation of the broader landscape.
Total training n = %d (50%%%% regrowth / 40%%%% matched / 10%%%% wide).

Seven sequential greedy passes over the same k = %d candidate list
produced seven distinct matched-negative pools (p1--p7), used to fill
negative budgets in the paired prevalence-sensitivity models (20/80,
5/95).
==============================
",
            n_regrowth, K_NN_used,
            100 * match_rate, n_matched, K_NN_used,
            100 * (1 - match_rate), n_unmatched, K_NN_used,
            median(nn1_dist_km), quantile(nn1_dist_km, 0.25), quantile(nn1_dist_km, 0.75),
            median(match_dist_km), quantile(match_dist_km, 0.25), quantile(match_dist_km, 0.75),
            n_matched, n_pair_drawn, 100 * n_pair_drawn / n_matched,
            n_regrowth * 2,
            K_NN_used
))
# ---- 13. Prevalence recalibration: 50/50 model → known prior ----------------
#
# The 50/50 model was trained on a balanced sample (pi_train = 0.50) but the
# true landscape prevalence of regrowth is much lower. Raw predicted
# probabilities from this model are therefore upward-biased and must be
# corrected before converting to area estimates.
#
# Two algebraically equivalent corrections are applied:
#
#   (a) Saerens et al. (2002) / Dal Pozzolo et al. (2015) general formula:
#         p_adj = [p * (pi_new / pi_train)] /
#                 [p * (pi_new / pi_train) + (1 - p) * ((1 - pi_new) / (1 - pi_train))]
#
#   (b) Dal Pozzolo et al. (2015) compact form (equivalent when pi_train = 0.5):
#         p_adj = (beta * p) / [(beta - 1) * p + 1]
#       where beta = (pi_new / (1 - pi_new)) / (pi_train / (1 - pi_train))
#                  = odds_new / odds_train
#
# With pi_train = 0.50 and pi_new = 0.05:
#   beta = (0.05/0.95) / (0.50/0.50) = 1/19
#
# The two formulas produce identical results (verified below).
# The 5/95 model (trained at natural prevalence) is reported alongside
# as an independent cross-check.

message("\n--- Prevalence recalibration (50/50 model → pi = 0.05) ---")

# ---- Prior assumptions ----
pi_train <- 0.50   # training prevalence (balanced 50/50 sample)
pi_new   <- 0.05   # estimated true landscape prevalence of regrowth

# ---- (a) Saerens / Dal Pozzolo general formula ----
recalibrate_saerens <- function(p, pi_train, pi_new) {
  num <- p * (pi_new / pi_train)
  den <- num + (1 - p) * ((1 - pi_new) / (1 - pi_train))
  num / den
}

# ---- (b) Dal Pozzolo compact form ----
beta <- (pi_new / (1 - pi_new)) / (pi_train / (1 - pi_train))   # = 1/19

recalibrate_dalpozzolo <- function(p, beta) {
  (beta * p) / ((beta - 1) * p + 1)
}

# ---- Load inputs ----
rf_50        <- readRDS(file.path(output_dir, "rf_prev_50_50.rds"))
pred_sample  <- readRDS(file.path(output_dir, "pred_sample.rds"))
total_area_ha <- readRDS(file.path(output_dir, "total_domain_area_ha.rds"))

# Raw predicted P(regrowth) from the 50/50 model over the non-regrowth sample
raw_probs <- predict(rf_50, pred_sample)$predictions[, "1"]

# ---- Apply corrections ----
adj_probs_saerens <- recalibrate_saerens(raw_probs, pi_train, pi_new)
adj_probs_dp      <- recalibrate_dalpozzolo(raw_probs, beta)

# Confirm the two methods agree (should be < 1e-10)
max_diff <- max(abs(adj_probs_saerens - adj_probs_dp))
message(sprintf("  Max abs diff between Saerens and Dal Pozzolo formulas: %.2e  %s",
                max_diff, ifelse(max_diff < 1e-8, "(equivalent — expected)", "WARNING: mismatch")))

# ---- Area estimates ----
raw_area_ha       <- mean(raw_probs)        * total_area_ha
adj_area_saerens  <- mean(adj_probs_saerens) * total_area_ha
adj_area_dp       <- mean(adj_probs_dp)      * total_area_ha

# 5/95 model as independent cross-check (no recalibration needed)
preds_5_95     <- readRDS(file.path(output_dir, "preds_rf_prev_5_95.rds"))
area_5_95_ha   <- mean(preds_5_95$nonreg) * total_area_ha

# ---- Report ----
cat(sprintf("\n  Prior assumptions: pi_train = %.2f  pi_new = %.2f  beta = %.6f (= 1/%.1f)\n",
            pi_train, pi_new, beta, 1 / beta))
cat(sprintf("\n  %-38s  %12s  %12s\n",  "Method",               "Mean prob", "Area (ha)"))
cat(strrep("-", 66), "\n")
cat(sprintf("  %-38s  %12.4f  %12s\n",
            "Raw 50/50 (uncorrected)",
            mean(raw_probs), format(round(raw_area_ha),      big.mark = ",")))
cat(sprintf("  %-38s  %12.4f  %12s\n",
            "Saerens recalibration",
            mean(adj_probs_saerens), format(round(adj_area_saerens), big.mark = ",")))
cat(sprintf("  %-38s  %12.4f  %12s\n",
            "Dal Pozzolo recalibration (beta = 1/19)",
            mean(adj_probs_dp), format(round(adj_area_dp),   big.mark = ",")))
cat(sprintf("  %-38s  %12.4f  %12s\n",
            "5/95 model (direct, no recalibration)",
            mean(preds_5_95$nonreg), format(round(area_5_95_ha),     big.mark = ",")))
cat(strrep("-", 66), "\n")
cat(sprintf("\n  Note: Saerens and Dal Pozzolo formulas are algebraically identical\n"))
cat(sprintf("  when pi_train = 0.5; any difference is floating-point rounding only.\n"))
cat(sprintf("  The 5/95 model provides an independent (unrecalibrated) benchmark.\n"))
