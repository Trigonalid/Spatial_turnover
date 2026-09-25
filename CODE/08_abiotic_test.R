# ============================================================
# 08_abiotic_gradient_annual.R
# Study: Spatial turnover amplifies with trophic level in
#        hyperdiverse food webs (Libra et al.)
#
# Abiotic (climate) gradient, ANNUAL resolution:
#   - annual aggregation of TerraClimate variables (ppt summed, temps averaged)
#   - Fig. S1: maps + violins
#   - Fig. S2: PCA
#   - Fig. S3: spatial regressions (lat + lon)
#   - Table S7: LMM spatial trends (Wald 95% CIs)
#
# Input:  DATA/monthly_data.csv   Output: output/fig/
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ggplot2)
library(rstatix)
library(patchwork)
library(ggrepel)
library(viridis)
library(maps)
library(lme4)
library(broom.mixed)
library(gt)
library(svglite)     # SVG export without cairo/XQuartz

monthly_csv <- "DATA/monthly_data.csv"     # input climate file
dir.create("output/fig", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Colour palette & locality order
# ------------------------------------------------------------
locality_levels <- c(
  "Elem","Morox","Niksek","Ohu",
  "Utai","Wamangu","Wanang","Yapsiei"
)
palette_localities <- c(
  Elem    = "#440154",
  Morox   = "#482878",
  Niksek  = "#3E4989",
  Ohu     = "#31688E",
  Utai    = "#26828E",
  Wamangu = "#1F9E89",
  Wanang  = "#35B779",
  Yapsiei = "#6CCE59"
)

# ------------------------------------------------------------
# 2. Study localities
# ------------------------------------------------------------
point_data <- tibble(
  locality = locality_levels,
  lat = c(
    -4.8167,-4.0167,-4.7000,-5.2333,
    -3.3841,-3.7871,-5.2309,-4.6283
  ),
  lon = c(
    143.9167,144.1000,142.5333,145.6833,
    141.5859,143.6521,145.1818,141.0973
  )
)

# ------------------------------------------------------------
# 3. Load extracted climate data
# ------------------------------------------------------------
monthly_data <- read_csv(monthly_csv, show_col_types = FALSE) %>%
  mutate(
    date     = as.Date(date),
    locality = factor(locality, levels = locality_levels)
  )

# ------------------------------------------------------------
# 4. Annual aggregation (ppt summed, temperatures averaged)
# ------------------------------------------------------------
annual <- monthly_data %>%
  group_by(locality, variable, year) %>%
  summarise(
    annual_value = if (first(variable) == "ppt") {
      sum(value, na.rm = TRUE)
    } else {
      mean(value, na.rm = TRUE)
    },
    .groups = "drop"
  )

annual <- annual %>%
  pivot_wider(names_from = variable, values_from = annual_value) %>%
  mutate(tavg = (tmin + tmax) / 2) %>%
  pivot_longer(
    cols = c(tmin, tmax, tavg, ppt),
    names_to = "variable",
    values_to = "annual_value"
  )

# ------------------------------------------------------------
# 5. Dataset for plotting (exclude 1998)
# ------------------------------------------------------------
annual_plot <- annual %>%
  filter(year != 1998) %>%
  mutate(
    variable_code = variable,
    variable = factor(
      variable,
      levels = c("tmin","tmax","tavg","ppt"),
      labels = c(
        "Minimum temperature (°C)",
        "Maximum temperature (°C)",
        "Average temperature (°C)",
        "Average precipitation (mm yr⁻¹)"
      )
    ),
    locality = factor(locality, levels = locality_levels)
  )

# ------------------------------------------------------------
# 6. Violin plot function
# ------------------------------------------------------------
plot_violin <- function(varname) {
  
  df <- filter(annual_plot, variable == varname)
  
  subtitle_txt <- case_when(
    varname == "Minimum temperature (°C)" ~ "Annual minimum temperature",
    varname == "Maximum temperature (°C)" ~ "Annual maximum temperature",
    varname == "Average temperature (°C)" ~ "Annual mean temperature",
    varname == "Average precipitation (mm yr⁻¹)" ~ "Annual precipitation"
  )
  
  y_lab <- if (varname == "Average precipitation (mm yr⁻¹)") "Precipitation (mm)" else "Temperature (°C)"
  
  ggplot(df, aes(locality, annual_value, fill = locality)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 1) +
    scale_fill_manual(values = palette_localities, limits = locality_levels, drop = FALSE) +
    theme_bw(base_size = 18) +
    theme(
      axis.text.x   = element_text(color = "black", angle = 45, hjust = 1),
      axis.text.y   = element_text(color = "black"),
      axis.title.y  = element_text(color = "black"),
      plot.subtitle = element_text(color = "black"),
      legend.position = "none"
    ) +
    labs(x = "", y = y_lab, subtitle = subtitle_txt)
}

# ------------------------------------------------------------
# 7. Maps
# ------------------------------------------------------------
world_map <- map_data("world")

plot_map <- function(varname) {
  
  df_map <- annual_plot %>%
    filter(variable == varname) %>%
    group_by(locality) %>%
    summarise(mean_value = mean(annual_value, na.rm = TRUE), .groups = "drop") %>%
    left_join(point_data, by = "locality")
  
  ggplot() +
    geom_polygon(data = world_map, aes(long, lat, group = group),
                 fill = "grey95", color = "grey70", linewidth = 0.3) +
    geom_point(data = df_map, aes(lon, lat, color = mean_value), size = 4) +
    geom_text_repel(data = df_map, aes(lon, lat, label = locality),
                    size = 4, fontface = "bold", color = "black") +
    scale_color_viridis_c(option = "D", name = "Mean annual value") +
    coord_cartesian(xlim = c(140,147), ylim = c(-6,-2)) +
    labs(title = varname, x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 18) +
    theme(
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title   = element_text(color = "black", face = "bold")
    )
}

# ------------------------------------------------------------
# 8. Fig. S1 - combined MAP + VIOLIN figure
# ------------------------------------------------------------
figure_maps_violins <-
  (plot_map("Average temperature (°C)") | plot_map("Average precipitation (mm yr⁻¹)")) /
  (plot_violin("Minimum temperature (°C)") | plot_violin("Maximum temperature (°C)")) /
  (plot_violin("Average temperature (°C)") | plot_violin("Average precipitation (mm yr⁻¹)"))

figure_maps_violins <- figure_maps_violins +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 20, face = "bold")))

ggsave("output/fig/Figure_S1_climate_maps_violins.svg", figure_maps_violins,
       width = 21, height = 18, device = svglite::svglite)
ggsave("output/fig/Figure_S1_climate_maps_violins.pdf", figure_maps_violins,
       width = 21, height = 18)
ggsave("output/fig/Figure_S1_climate_maps_violins.png", figure_maps_violins,
       width = 21, height = 18, dpi = 300)

# ------------------------------------------------------------
# 9. Fig. S2 - PCA (all variables, all years excl. 1998)
# ------------------------------------------------------------
run_pca_plot <- function(data, variables, title) {
  
  pca_data <- data %>%
    summarise(annual_value = mean(annual_value, na.rm = TRUE),
              .by = c(locality, year, variable_code)) %>%
    pivot_wider(names_from = variable_code, values_from = annual_value) %>%
    drop_na()
  
  pca_res <- prcomp(pca_data %>% select(all_of(variables)), center = TRUE, scale. = TRUE)
  
  scores <- as_tibble(pca_res$x) %>%
    bind_cols(pca_data %>% select(locality, year)) %>%
    mutate(locality = factor(locality, levels = locality_levels))
  
  hulls <- scores %>% group_by(locality) %>% slice(chull(PC1, PC2))
  
  centroids <- scores %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .by = locality) %>%
    mutate(locality = factor(locality, levels = locality_levels))
  
  var_exp <- summary(pca_res)$importance["Proportion of Variance", 1:2] * 100
  
  ggplot(scores, aes(PC1, PC2, color = locality)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_polygon(data = hulls, aes(fill = locality), alpha = 0.15, color = NA) +
    geom_point(data = centroids, size = 5, shape = 4, stroke = 2) +
    scale_color_manual(values = palette_localities, limits = locality_levels, drop = FALSE) +
    scale_fill_manual(values = palette_localities, limits = locality_levels, drop = FALSE) +
    theme_bw(base_size = 20) +
    labs(title = title,
         x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
         color = "Locality", fill = "Locality")
}

pca_all <- run_pca_plot(annual_plot, c("tmin","tmax","tavg","ppt"), title = "")

ggsave("output/fig/Figure_S2_climate_PCA.svg", pca_all,
       width = 10, height = 8, device = svglite::svglite)
ggsave("output/fig/Figure_S2_climate_PCA.pdf", pca_all, width = 10, height = 8)
ggsave("output/fig/Figure_S2_climate_PCA.png", pca_all, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------
# 10. LMM table: spatial trends (Wald 95% CIs)
# ------------------------------------------------------------
reg_data <- annual_plot %>%
  left_join(point_data, by = "locality") %>%
  mutate(lat_z = as.numeric(scale(lat)),
         lon_z = as.numeric(scale(lon)))

models <- reg_data %>%
  group_by(variable) %>%
  group_map(~ lmer(annual_value ~ lat_z + lon_z + (1 | locality),
                   data = .x, REML = FALSE))
names(models) <- levels(reg_data$variable)

lm_table <- bind_rows(
  lapply(models, tidy, effects = "fixed", conf.int = TRUE),
  .id = "variable"
)

gt_table <- lm_table %>%
  filter(term %in% c("lat_z", "lon_z")) %>%
  mutate(term = recode(term,
                       lat_z = "Latitude (standardised)",
                       lon_z = "Longitude (standardised)")) %>%
  select(variable, term, estimate, conf.low, conf.high) %>%
  gt(groupname_col = "variable") %>%
  fmt_number(columns = c(estimate, conf.low, conf.high), decimals = 3) %>%
  cols_label(estimate = "Estimate", conf.low = "CI low", conf.high = "CI high") %>%
  tab_header(title = "Table S7: Spatial trends in annual climate variables",
             subtitle = "Linear mixed-effects models (Wald 95% confidence intervals)")

print(gt_table)
write.csv(lm_table, "output/rds/Table_S7_LMM_spatial_trends.csv", row.names = FALSE)
gt::gtsave(gt_table, "output/Table_S7_LMM_spatial_trends.docx")

# ------------------------------------------------------------
# 11. Fig. S3 - spatial regressions (lat + lon)
# ------------------------------------------------------------
plot_regression <- function(varname, axis = c("lat", "lon")) {
  
  axis  <- match.arg(axis)
  df    <- reg_data %>% filter(variable == varname)
  x_var <- if (axis == "lat") "lat" else "lon"
  x_lab <- if (axis == "lat") "Latitude" else "Longitude"
  y_lab <- if (varname == "Average precipitation (mm yr⁻¹)") "Precipitation (mm)" else "Temperature (°C)"
  
  subtitle_txt <- dplyr::case_when(
    varname == "Minimum temperature (°C)" ~ "Annual minimum temperature",
    varname == "Maximum temperature (°C)" ~ "Annual maximum temperature",
    varname == "Average temperature (°C)" ~ "Annual mean temperature",
    varname == "Average precipitation (mm yr⁻¹)" ~ "Annual precipitation"
  )
  
  fml    <- as.formula(paste("annual_value ~", x_var))
  lm_fit <- lm(fml, data = df)
  r2_lab <- paste0("R² = ", round(summary(lm_fit)$r.squared, 2))
  
  ggplot(df, aes_string(x = x_var, y = "annual_value", color = "locality")) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(aes(group = locality), method = "lm", se = FALSE, linetype = "dashed") +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, linetype = "dashed", color = "black") +
    annotate("text", x = Inf, y = Inf, label = r2_lab,
             hjust = 1.1, vjust = 1.3, size = 5, color = "black") +
    scale_color_manual(values = palette_localities, limits = locality_levels, drop = FALSE) +
    theme_bw(base_size = 18) +
    theme(
      axis.text     = element_text(color = "black"),
      axis.title    = element_text(color = "black"),
      plot.subtitle = element_text(color = "black"),
      legend.position = "right"
    ) +
    labs(x = x_lab, y = y_lab, subtitle = subtitle_txt, color = "Locality")
}

figure_regressions <-
  (plot_regression("Minimum temperature (°C)", "lat") |
     plot_regression("Maximum temperature (°C)", "lat") |
     plot_regression("Average temperature (°C)", "lat") |
     plot_regression("Average precipitation (mm yr⁻¹)", "lat")) /
  (plot_regression("Minimum temperature (°C)", "lon") |
     plot_regression("Maximum temperature (°C)", "lon") |
     plot_regression("Average temperature (°C)", "lon") |
     plot_regression("Average precipitation (mm yr⁻¹)", "lon"))

figure_regressions <- figure_regressions +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 18, face = "bold", color = "black"))) +
  plot_layout(guides = "collect")

ggsave("output/fig/Figure_S3_climate_regressions.svg", figure_regressions,
       width = 24, height = 12, device = svglite::svglite)
ggsave("output/fig/Figure_S3_climate_regressions.pdf", figure_regressions,
       width = 24, height = 12)
ggsave("output/fig/Figure_S3_climate_regressions.png", figure_regressions,
       width = 24, height = 12, dpi = 300)