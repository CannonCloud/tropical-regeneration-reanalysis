library(data.table)
library(ggplot2)
library(stringr)

# ============================================================================
# 1. DATA PREP
# ============================================================================
df <- fread("pnv_final_dataset.csv")

# Create IDs
df[, place_id_3deg := paste0(round(longitude / 3) * 3, "_", round(latitude / 3) * 3)]
df[, place_id_1deg := paste0(round(longitude), "_", round(latitude))]

# Create Dummies (Only Biomes needed now, since we are dropping LC)
df[, `:=`(
  biome_rainforest   = as.integer(bio_biome_id == 1),
  biome_dry_forest   = as.integer(bio_biome_id == 2),
  biome_coniferous   = as.integer(bio_biome_id == 3)
)]

# ============================================================================
# 2. DEFINE VARIABLE SETS (REMOVING lc_ VARS)
# ============================================================================

# Base Biophysical list (No Landcover 'lc_' variables)
bio_vars_base <- c(
  "pnv_probability", 
  "bio_min_temp_coldest", "bio_precip_annual", "bio_temp_mean", "bio_elevation", 
  "bio_slope", "bio_npp_mean", "bio_fire_freq", "bio_soil_organic_carbon", 
  "bio_soil_ph", "bio_soil_clay", "bio_dist_forest_2018", "bio_forest_density_1km2",
  "paper_dist_water", "bio_canopy_height_1km_mean"
)

# Global gets Biomes + Base
bio_global <- c("biome_rainforest", "biome_dry_forest", "biome_coniferous", bio_vars_base)

# Local gets just Base (No Biomes, No Landcover)
bio_local <- bio_vars_base

socio_old <- c("paper_ghs_pop", "paper_gdp_2015", "paper_roads_density", "paper_protected_binary", "paper_dist_urban", "bio_cropland_density")
socio_new <- c(
  "highres_travel", "highres_nightlight_density_2km", "highres_worldpop_density_5km",
  "highres_ghm", "econ_opp_cost_usd_ha", "highres_dist_settlement", 
  "highres_dist_agriculture", "highres_agriculture_density_1km",
  "highres_dist_cultivated_grass", "highres_cultivated_grass_density_2km",
  "highres_dist_cropland", "highres_cropland_density_2km",
  "bio_dist_plantation", "bio_plantation_density_1km"
)

# ============================================================================
# 3. ANALYSIS FUNCTION
# ============================================================================

analyze_scale <- function(data, group_col = NULL, bio_vars) {
  
  # 1. Filter Variables
  required_vars <- c(bio_vars, socio_old, socio_new)
  dt <- data[, ..required_vars]
  
  # 2. Handle Grouping (Global vs Local)
  if (!is.null(group_col)) {
    dt[, group_id := data[[group_col]]]
    
    # FILTER: Keep only groups with N >= 10
    valid_groups <- dt[, .N, by = group_id][N >= 10, group_id]
    dt <- dt[group_id %in% valid_groups]
    
    # DEMEAN (Fixed Effects)
    dt[, (required_vars) := lapply(.SD, function(x) x - mean(x, na.rm=TRUE)), 
       by = group_id, .SDcols = required_vars]
    
    dt[, group_id := NULL] 
  }
  
  # 3. Calculate Correlations
  m_old <- suppressWarnings(cor(dt[, ..socio_old], dt[, ..bio_vars], use = "pairwise.complete.obs"))
  m_new <- suppressWarnings(cor(dt[, ..socio_new], dt[, ..bio_vars], use = "pairwise.complete.obs"))
  
  # 4. Format for Plotting (FIXED HERE)
  # Combine matrices
  m_combined <- rbind(m_old, m_new)
  
  # Explicitly rename rows so the regex works later
  rownames(m_combined) <- c(
    paste0("OLD_", rownames(m_old)), 
    paste0("NEW_", rownames(m_new))
  )
  
  # Convert to long format
  as.data.table(m_combined, keep.rownames = "Socioeconomic") |>
    melt(id.vars = "Socioeconomic", variable.name = "Biophysical", value.name = "Correlation") |>
    transform(Category = ifelse(str_detect(Socioeconomic, "^OLD_"), "OLD", "NEW"))
}

# ============================================================================
# 4. PLOTTING FUNCTION
# ============================================================================

plot_corr <- function(data, title_text) {
  ggplot(data, aes(x = Biophysical, y = Socioeconomic, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 2.5, na.rm = TRUE) + # <--- ADDED BACK
    scale_fill_gradient2(low = "#D73027", mid = "white", high = "#1A9850", 
                         midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.border = element_rect(color = "black", fill = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(title = title_text, x = NULL, y = NULL)
}

# ============================================================================
# 5. GENERATE THE 3 FIGURES
# ============================================================================

# Figure 1: Global
data_global <- analyze_scale(df, group_col = NULL, bio_vars = bio_global)
p1 <- plot_corr(data_global, "1. Global Correlations (Standard)")

# Figure 2: 3-Degree (Local)
data_3deg <- analyze_scale(df, group_col = "place_id_3deg", bio_vars = bio_local)
p2 <- plot_corr(data_3deg, "2. Within 3-Degree Grid (Fixed Effects, N>=10)")

# Figure 3: 1-Degree (Local)
data_1deg <- analyze_scale(df, group_col = "place_id_1deg", bio_vars = bio_local)
p3 <- plot_corr(data_1deg, "3. Within 1-Degree Grid (Fixed Effects, N>=10)")

# Display them
print(p1)
print(p2)
print(p3)