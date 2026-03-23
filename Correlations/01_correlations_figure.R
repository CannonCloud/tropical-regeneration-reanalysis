# ============================================================
# 01_correlations_figure.R
# Brazil Natural Regeneration — Williams et al. Replication
#
# Produces:
#   fig_corr_main.pdf/.png   — focused 6x11 heatmap for main text
#   fig_corr_global.pdf/.png — full matrix, global (appendix)
#   fig_corr_3deg.pdf/.png   — full matrix, within 3-degree FE (appendix)
#   fig_corr_1deg.pdf/.png   — full matrix, within 1-degree FE (appendix)
#
# Row structure in appendix figures:
#   OUTCOME            — PNV probability (Williams et al. prediction target)
#   SOCIOECONOMIC NEW  — variables not included in Williams et al.
#   SOCIOECONOMIC OLD  — variables included in Williams et al.
#
# LaTeX caption for fig_corr_main:
#   Pearson correlations between six socioeconomic variables (rows) and
#   eleven biophysical predictors used by Williams et al. (columns),
#   computed at N=200,000 randomly sampled tropical forest pixels.
#   Colour intensity reflects absolute correlation magnitude; sign shown
#   as cell label. Both positive and negative correlations contribute to
#   the omitted-variable problem regardless of direction. Key finding:
#   Dist. to forest, one of the most important predictors in the
#   Williams et al. model, correlates 0.55 with distance to plantation,
#   meaning the model partly measures proximity to sterile monoculture
#   rather than natural forest seed sources. Human modification index
#   correlates -0.42 with net primary productivity, embedding economic
#   land-use intensity into an ostensibly biophysical predictor. Full
#   correlation matrices at global and within-grid scales shown in
#   Supplementary Figures X-Z.
#
# LaTeX caption for appendix figures:
#   Pearson correlations between socioeconomic variables (rows) and
#   biophysical predictors (columns). Colour intensity = |r|, sign
#   shown as cell label. Top panel (OUTCOME): correlations between the
#   Williams et al. prediction target (PNV probability) and each
#   biophysical predictor. Middle panel (SOCIOECONOMIC NEW): variables
#   not included in Williams et al. Bottom panel (SOCIOECONOMIC OLD):
#   variables included in Williams et al. (a) Global. (b) Within
#   3-degree grid fixed effects (N>=10). (c) Within 1-degree grid
#   fixed effects (N>=10). Correlations persisting under fixed effects
#   reflect within-landscape structure rather than broad geographic
#   confounding.
# ============================================================

library(data.table)
library(ggplot2)
library(stringr)

setwd("# replace with your WD")

output_dir <- "corr_outputs"
dir.create(output_dir, showWarnings = FALSE)

# ============================================================
# 1. DATA PREP
# ============================================================

df <- fread("gee_outputs/pnv_final_dataset.csv")

df[, place_id_3deg := paste0(round(longitude / 3) * 3, "_",
                              round(latitude  / 3) * 3)]
df[, place_id_1deg := paste0(round(longitude), "_", round(latitude))]

df[, `:=`(
  biome_rainforest = as.integer(bio_biome_id == 1),
  biome_dry_forest = as.integer(bio_biome_id == 2),
  biome_coniferous = as.integer(bio_biome_id == 3)
)]

# ============================================================
# 2. VARIABLE SETS
# ============================================================

# Biophysical — PNV excluded here, handled separately as outcome row
bio_global_no_pnv <- c(
  "biome_rainforest", "biome_dry_forest", "biome_coniferous",
  "bio_min_temp_coldest", "bio_precip_annual", "bio_temp_mean",
  "bio_elevation", "bio_slope", "bio_npp_mean", "bio_fire_freq",
  "bio_soil_organic_carbon", "bio_soil_ph", "bio_soil_clay",
  "bio_dist_forest_2018", "bio_forest_density_1km2",
  "paper_dist_water", "bio_canopy_height_1km_mean"
)

bio_local_no_pnv <- c(
  "bio_min_temp_coldest", "bio_precip_annual", "bio_temp_mean",
  "bio_elevation", "bio_slope", "bio_npp_mean", "bio_fire_freq",
  "bio_soil_organic_carbon", "bio_soil_ph", "bio_soil_clay",
  "bio_dist_forest_2018", "bio_forest_density_1km2",
  "paper_dist_water", "bio_canopy_height_1km_mean"
)

# Main figure biophysical columns (11 selected, no PNV as column)
bio_main <- c(
  "bio_min_temp_coldest", "bio_precip_annual", "bio_temp_mean",
  "bio_elevation", "bio_slope", "bio_npp_mean",
  "bio_soil_organic_carbon", "bio_soil_ph",
  "bio_dist_forest_2018", "bio_forest_density_1km2",
  "bio_canopy_height_1km_mean"
)

socio_old <- c(
  "paper_ghs_pop", "paper_gdp_2015", "paper_roads_density",
  "paper_protected_binary", "paper_dist_urban", "bio_cropland_density"
)

socio_new <- c(
  "highres_travel", "highres_nightlight_density_2km",
  "highres_worldpop_density_5km", "highres_ghm",
  "econ_opp_cost_usd_ha", "highres_dist_settlement",
  "highres_dist_agriculture", "highres_agriculture_density_1km",
  "highres_dist_cultivated_grass", "highres_cultivated_grass_density_2km",
  "highres_dist_cropland", "highres_cropland_density_2km",
  "bio_dist_plantation", "bio_plantation_density_1km"
)

socio_main <- c(
  "highres_travel",
  "highres_ghm",
  "highres_cropland_density_2km",
  "highres_cultivated_grass_density_2km",
  "bio_cropland_density",
  "bio_dist_plantation"
)

# ============================================================
# 3. CLEAN DISPLAY NAMES
# ============================================================

se_labels <- c(
  pnv_probability                      = "PNV\nprobability",
  highres_travel                        = "Travel time\nto market",
  highres_nightlight_density_2km        = "Nightlight density\n(2km)",
  highres_worldpop_density_5km          = "Population density\n(5km)",
  highres_ghm                           = "Human\nmodification\nindex",
  econ_opp_cost_usd_ha                  = "Opportunity cost\n(USD/ha)",
  highres_dist_settlement               = "Dist. to\nsettlement",
  highres_dist_agriculture              = "Dist. to\nagriculture",
  highres_agriculture_density_1km       = "Agriculture density\n(1km)",
  highres_dist_cultivated_grass         = "Dist. to\npasture",
  highres_cultivated_grass_density_2km  = "Pasture density\n(2km)",
  highres_dist_cropland                 = "Dist. to\ncropland",
  highres_cropland_density_2km          = "Cropland density\n(2km)",
  bio_dist_plantation                   = "Dist. to\nplantation",
  bio_plantation_density_1km            = "Plantation density\n(1km)",
  paper_ghs_pop                         = "Population\n(Williams)",
  paper_gdp_2015                        = "GDP 2015\n(Williams)",
  paper_roads_density                   = "Road density\n(Williams)",
  paper_protected_binary                = "Protected area\n(Williams)",
  paper_dist_urban                      = "Dist. to urban\n(Williams)",
  bio_cropland_density                  = "Cropland density\n(Williams, 5km)"
)

bio_labels <- c(
  biome_rainforest           = "Biome: rainforest",
  biome_dry_forest           = "Biome: dry forest",
  biome_coniferous           = "Biome: coniferous",
  bio_min_temp_coldest       = "Min temp coldest",
  bio_precip_annual          = "Annual precip",
  bio_temp_mean              = "Mean temp",
  bio_elevation              = "Elevation",
  bio_slope                  = "Slope",
  bio_npp_mean               = "NPP",
  bio_fire_freq              = "Fire frequency",
  bio_soil_organic_carbon    = "Soil org. C",
  bio_soil_ph                = "Soil pH",
  bio_soil_clay              = "Soil clay",
  bio_dist_forest_2018       = "Dist. to forest",
  bio_forest_density_1km2    = "Forest density",
  paper_dist_water           = "Dist. to water",
  bio_canopy_height_1km_mean = "Canopy height"
)

# ============================================================
# 4. ANALYSIS FUNCTION
# ============================================================

analyze_scale <- function(data, group_col = NULL, bio_vars,
                           se_vars_old, se_vars_new) {

  required_vars <- unique(c("pnv_probability", bio_vars,
                             se_vars_old, se_vars_new))
  dt <- data[, ..required_vars]

  if (!is.null(group_col)) {
    dt[, group_id := data[[group_col]]]
    valid_groups <- dt[, .N, by = group_id][N >= 10, group_id]
    dt <- dt[group_id %in% valid_groups]
    dt[, (required_vars) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)),
       by = group_id, .SDcols = required_vars]
    dt[, group_id := NULL]
  }

  # SE vs bio correlations
  m_old <- suppressWarnings(
    cor(dt[, ..se_vars_old], dt[, ..bio_vars], use = "pairwise.complete.obs")
  )
  m_new <- suppressWarnings(
    cor(dt[, ..se_vars_new], dt[, ..bio_vars], use = "pairwise.complete.obs")
  )

  # PNV outcome row vs bio
  m_pnv <- suppressWarnings(
    cor(dt[["pnv_probability"]],
        as.matrix(dt[, ..bio_vars]),
        use = "pairwise.complete.obs")
  )
  rownames(m_pnv) <- "pnv_probability"

  # Combine with category labels
  m_combined <- rbind(m_pnv, m_old, m_new)
  rownames(m_combined) <- c(
    "RESULT_pnv_probability",
    paste0("OLD_", rownames(m_old)),
    paste0("NEW_", rownames(m_new))
  )

  as.data.table(m_combined, keep.rownames = "Socioeconomic") |>
    melt(id.vars      = "Socioeconomic",
         variable.name = "Biophysical",
         value.name    = "Correlation") |>
    transform(Category = ifelse(str_detect(Socioeconomic, "^RESULT_"), "RESULT",
                         ifelse(str_detect(Socioeconomic, "^OLD_"), "OLD", "NEW")))
}

# ============================================================
# 5. SHARED PLOT FUNCTION
# ============================================================

plot_corr <- function(data, se_order, bio_order,
                      facet = TRUE, title_text = NULL) {
  
  dt <- as.data.table(data)
  
  dt[, se_raw    := str_remove(Socioeconomic, "^(RESULT|OLD|NEW)_")]
  dt[, se_clean  := se_labels[se_raw]]
  dt[, bio_clean := bio_labels[as.character(Biophysical)]]
  
  dt[, Category := factor(Category,
                          levels = c("RESULT", "NEW", "OLD"),
                          labels = c("RESULT",
                                     "SOCIOECONOMIC — NEW",
                                     "SOCIOECONOMIC — OLD"))]
  
  se_clean_order  <- rev(se_labels[se_order])
  bio_clean_order <- bio_labels[bio_order]
  
  dt[, se_clean  := factor(se_clean,  levels = se_clean_order)]
  dt[, bio_clean := factor(bio_clean, levels = bio_clean_order)]
  dt[, abs_corr  := abs(Correlation)]
  dt[, label     := sprintf("%.2f", round(Correlation, 2))]
  
  p <- ggplot(dt, aes(x = bio_clean, y = se_clean, fill = abs_corr)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 2.4, colour = "grey20") +
    scale_fill_gradient(
      low    = "white",
      high   = "#2166AC",
      limits = c(0, 1),
      name   = "|r|",
      # Fix for legend label cutting off:
      guide  = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y       = element_text(size = 7, lineheight = 0.9),
      panel.grid        = element_blank(),
      legend.position   = "right",
      legend.key.height = unit(1.0, "cm"),
      plot.title        = element_text(face = "bold", size = 9),
      plot.subtitle     = element_text(size = 8, colour = "grey40", margin = margin(b=10))
    ) +
    labs(title    = title_text,
         subtitle = "Biophysical predictors \u2192",
         x = NULL, 
         y = "\u2190 Socioeconomic variables") # Added side label
  
  if (facet) {
    p <- p +
      facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
      theme(
        strip.text   = element_text(face = "bold", size = 7.5),
        strip.background = element_rect(fill = "grey92", color = "grey70"),
        panel.border = element_rect(color = "black", fill = NA)
      )
  }
  
  p
}

# ============================================================
# 6. COMPUTE CORRELATIONS
# ============================================================

message("Computing global correlations...")
data_global <- analyze_scale(df, group_col = NULL,
                              bio_vars    = bio_global_no_pnv,
                              se_vars_old = socio_old,
                              se_vars_new = socio_new)

message("Computing 3-degree fixed effects...")
data_3deg <- analyze_scale(df, group_col  = "place_id_3deg",
                            bio_vars    = bio_local_no_pnv,
                            se_vars_old = socio_old,
                            se_vars_new = socio_new)

message("Computing 1-degree fixed effects...")
data_1deg <- analyze_scale(df, group_col  = "place_id_1deg",
                            bio_vars    = bio_local_no_pnv,
                            se_vars_old = socio_old,
                            se_vars_new = socio_new)

# ============================================================
# 7. MAIN FIGURE
# — 6 SE rows, 11 bio columns, no facet, no PNV outcome row
# ============================================================

message("Building main figure...")

data_main <- data_global[
  str_remove(Socioeconomic, "^(RESULT|OLD|NEW)_") %in% socio_main &
    Biophysical %in% bio_main
]

p_main <- plot_corr(
  data       = data_main,
  se_order   = socio_main,
  bio_order  = bio_main,
  facet      = FALSE,
  title_text = NULL
)

ggsave(file.path(output_dir, "fig_corr_main.pdf"),
       p_main, width = 160, height = 90, units = "mm",
       device = cairo_pdf)
ggsave(file.path(output_dir, "fig_corr_main.png"),
       p_main, width = 160, height = 90, units = "mm", dpi = 300)
message("Saved: fig_corr_main.pdf/.png")

# ============================================================
# 8. APPENDIX FIGURES — full matrices with 3-level facet
# ============================================================

message("Building appendix figures...")

p_global <- plot_corr(
  data       = data_global,
  se_order   = c("pnv_probability", socio_new, socio_old),
  bio_order  = bio_global_no_pnv,
  facet      = TRUE,
  title_text = "Global correlations"
)

p_3deg <- plot_corr(
  data       = data_3deg,
  se_order   = c("pnv_probability", socio_new, socio_old),
  bio_order  = bio_local_no_pnv,
  facet      = TRUE,
  title_text = "Within 3-degree grid (fixed effects, N\u226510)"
)

p_1deg <- plot_corr(
  data       = data_1deg,
  se_order   = c("pnv_probability", socio_new, socio_old),
  bio_order  = bio_local_no_pnv,
  facet      = TRUE,
  title_text = "Within 1-degree grid (fixed effects, N\u226510)"
)

ggsave(file.path(output_dir, "fig_corr_global.pdf"),
       p_global, width = 280, height = 230, units = "mm",
       device = cairo_pdf)
ggsave(file.path(output_dir, "fig_corr_global.png"),
       p_global, width = 280, height = 230, units = "mm", dpi = 200)
message("Saved: fig_corr_global.pdf/.png")

ggsave(file.path(output_dir, "fig_corr_3deg.pdf"),
       p_3deg, width = 260, height = 230, units = "mm",
       device = cairo_pdf)
ggsave(file.path(output_dir, "fig_corr_3deg.png"),
       p_3deg, width = 260, height = 230, units = "mm", dpi = 200)
message("Saved: fig_corr_3deg.pdf/.png")

ggsave(file.path(output_dir, "fig_corr_1deg.pdf"),
       p_1deg, width = 260, height = 230, units = "mm",
       device = cairo_pdf)
ggsave(file.path(output_dir, "fig_corr_1deg.png"),
       p_1deg, width = 260, height = 230, units = "mm", dpi = 200)
message("Saved: fig_corr_1deg.pdf/.png")

message("\n01_correlations_figure.R complete.")
