# ============================================================
# COMPLETE CLIMATE ANALYSIS SCRIPT
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------

library(QBMS)
library(terra)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)
library(rstatix)
library(patchwork)
library(ggrepel)
library(viridis)
library(maps)

theme_set(
  theme_bw(base_size = 18) +
    theme(
      text            = element_text(color = "black"),
      axis.text       = element_text(color = "black"),
      axis.title      = element_text(color = "black"),
      plot.title      = element_text(color = "black", face = "bold"),
      plot.subtitle   = element_text(color = "black"),
      legend.text     = element_text(color = "black"),
      legend.title    = element_text(color = "black"),
      strip.text      = element_text(color = "black"),
      plot.tag        = element_text(color = "black", face = "bold"),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
)

# ------------------------------------------------------------
# CONTROL FLAGS (IMPORTANT)
# ------------------------------------------------------------

RUN_DOWNLOAD  <- FALSE  # step 5
RUN_EXTRACTION <- TRUE # step 6
RUN_ANNUAL    <- TRUE  # step 7



# ------------------------------------------------------------
# 1. Color palette & locality order (GLOBAL)
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
# 2. Directories
# ------------------------------------------------------------

base_dir <- getwd()
data_dir <- file.path(base_dir, "terraclimate_download")
dir.create(data_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Study localities
# ------------------------------------------------------------

point_data <- tibble(
  locality = locality_levels,
  lat = c(-4.8167,-4.0167,-4.7000,-5.2333,
          -3.3841,-3.7871,-5.2309,-4.6283),
  lon = c(143.9167,144.1000,142.5333,145.6833,
          141.5859,143.6521,145.1818,141.0973)
)

pts <- terra::vect(point_data[, c("lon","lat")], crs = "EPSG:4326")

# ------------------------------------------------------------
# 4. Variables and years
# ------------------------------------------------------------

vars <- c("tmin","tmax","ppt")
start_year <- 1970
end_year   <- 2006


# ------------------------------------------------------------
# 5. Download TerraClimate data (RUN ONCE)
# ------------------------------------------------------------

if (RUN_DOWNLOAD) {
  
  for (v in vars) {
    for (yr in start_year:end_year) {
      
      fname <- paste0("TerraClimate_", v, "_", yr, ".nc")
      fpath <- file.path(data_dir, fname)
      
      if (!file.exists(fpath)) {
        
        message("Downloading ", fname)
        
        QBMS::ini_terraclimate(
          from = paste0(yr, "-01-01"),
          to   = paste0(yr, "-12-31"),
          clim_vars = v,
          data_path = data_dir
        )
      }
    }
  }
}


# ------------------------------------------------------------
# 6. Extract monthly data
# ------------------------------------------------------------

# ------------------------------------------------------------
# 6. Extract monthly data
# ------------------------------------------------------------

if (RUN_EXTRACTION) {
  
  extract_var <- function(varname) {
    
    files <- list.files(
      data_dir,
      pattern = paste0("^TerraClimate_", varname, "_\\d{4}\\.nc$"),
      full.names = TRUE
    )
    
    out <- list()
    
    for (f in files) {
      
      r <- terra::rast(f)
      tt <- as.Date(terra::time(r))
      vals <- terra::extract(r, pts)[, -1]
      
      df <- expand.grid(
        locality = point_data$locality,
        date = tt
      )
      
      df$value <- as.vector(t(vals))
      df$year <- year(df$date)
      df$month <- month(df$date)
      df$variable <- varname
      
      out[[length(out) + 1]] <- as_tibble(df)
    }
    
    bind_rows(out)
  }
  
  monthly_data <- bind_rows(lapply(vars, extract_var))
  saveRDS(monthly_data, "monthly_data.rds")
  
} else {
  
  monthly_data <- readRDS("monthly_data.rds")
}


# ------------------------------------------------------------
# 7. Correct ANNUAL aggregation + tavg
# ------------------------------------------------------------

if (RUN_ANNUAL) {
  
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
    filter(variable %in% c("tmin","tmax","ppt")) %>%
    pivot_wider(names_from = variable, values_from = annual_value) %>%
    mutate(tavg = (tmin + tmax) / 2) %>%
    pivot_longer(
      cols = c(tmin, tmax, tavg, ppt),
      names_to = "variable",
      values_to = "annual_value"
    )
  
  saveRDS(annual, "annual_data.rds")
  
} else {
  
  annual <- readRDS("annual_data.rds")
}


# ------------------------------------------------------------
# 8. Dataset for plotting & PCA
# ------------------------------------------------------------

annual_plot <- annual %>%
  filter(year != 1998) %>%   # ⬅⬅⬅ TADY
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
# 9. Violin plot function
# ------------------------------------------------------------

plot_violin <- function(varname) {
  
  df <- filter(annual_plot, variable == varname)
  
  subtitle_txt <- case_when(
    varname == "Minimum temperature (°C)" ~ "Annual minimum temperature",
    varname == "Maximum temperature (°C)" ~ "Annual maximum temperature",
    varname == "Average temperature (°C)" ~ "Annual mean temperature",
    varname == "Average precipitation (mm yr⁻¹)" ~ "Annual precipitation"
  )
  
  y_lab <- if (varname == "Average precipitation (mm yr⁻¹)") {
    "Precipitation (mm)"
  } else {
    "Temperature (°C)"
  }
  
  ggplot(df, aes(locality, annual_value, fill = locality)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 1) +
    scale_fill_manual(
      values = palette_localities,
      limits = locality_levels,
      drop = FALSE
    ) +
    theme_bw(base_size = 18) +
    theme(
      axis.text.x  = element_text(color = "black", angle = 45, hjust = 1),
      axis.text.y  = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.subtitle = element_text(color = "black"),
      legend.position = "none"
    ) +
    labs(
      x = "",
      y = y_lab,
      subtitle = subtitle_txt
    )
}


# ------------------------------------------------------------
# 10. Maps
# ------------------------------------------------------------

world_map <- map_data("world")

map_data <- annual_plot %>%
  filter(variable %in% c(
    "Average temperature (°C)",
    "Average precipitation (mm yr⁻¹)"
  )) %>%
  group_by(locality, variable) %>%
  summarise(mean_value = mean(annual_value), .groups = "drop") %>%
  left_join(point_data, by = "locality")

plot_map <- function(varname) {
  
  df_map <- annual_plot %>%
    filter(variable == varname) %>%
    group_by(locality) %>%
    summarise(
      mean_value = mean(annual_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(point_data, by = "locality")
  
  ggplot() +
    geom_polygon(
      data = world_map,
      aes(long, lat, group = group),
      fill = "grey95",
      color = "grey70",
      linewidth = 0.3
    ) +
    geom_point(
      data = df_map,
      aes(lon, lat, color = mean_value),
      size = 4
    ) +
    geom_text_repel(
      data = df_map,
      aes(lon, lat, label = locality),
      size = 4,
      fontface = "bold",
      color = "black"
    ) +
    scale_color_viridis_c(
      option = "D",
      name = "Mean annual value"
    ) +
    coord_cartesian(xlim = c(140,147), ylim = c(-6,-2)) +
    labs(
      title = varname,
      x = "Longitude",
      y = "Latitude"
    ) +
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
# 11. Combined MAP + VIOLIN figure
# ------------------------------------------------------------

figure_maps_violins <-
  (plot_map("Average temperature (°C)") |
     plot_map("Average precipitation (mm yr⁻¹)")) /
  (plot_violin("Minimum temperature (°C)") |
     plot_violin("Maximum temperature (°C)")) /
  (plot_violin("Average temperature (°C)") |
     plot_violin("Average precipitation (mm yr⁻¹)"))

figure_maps_violins <-figure_maps_violins +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold")
    )
  )
figure_maps_violins


ggsave(
  "Figure_maps_violins.svg",
  figure_maps_violins,
  width = 21,
  height = 18, device = "pdf"
)
ggsave(
  "Figure_maps_violins.png",
  figure_maps_violins,
  width = 21,
  height = 18
)
# ------------------------------------------------------------
# 12. PCA (all variables, all years excl. 1998)
# ------------------------------------------------------------

run_pca_plot <- function(data, variables, title) {
  
  pca_data <- data %>%
    summarise(
      annual_value = mean(annual_value, na.rm = TRUE),
      .by = c(locality, year, variable_code)
    ) %>%
    pivot_wider(names_from = variable_code, values_from = annual_value) %>%
    drop_na()
  
  pca_res <- prcomp(
    pca_data %>% select(all_of(variables)),
    center = TRUE,
    scale. = TRUE
  )
  
  scores <- as_tibble(pca_res$x) %>%
    bind_cols(pca_data %>% select(locality, year)) %>%
    mutate(locality = factor(locality, levels = locality_levels))
  
  hulls <- scores %>%
    group_by(locality) %>%
    slice(chull(PC1, PC2))
  
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
    labs(
      title = title,
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
      color = "Locality",
      fill  = "Locality"
    )
  
}

pca_all <- run_pca_plot(
  annual_plot,
  c("tmin","tmax","tavg","ppt"),
  title = ""
)

pca_all



ggsave(
  "Figure_PCA_all_variables.svg",
  pca_all,
  width = 10,
  height = 8,
  device = "pdf"
)

ggsave(
  "Figure_PCA_all_variables.png",
  pca_all,
  width = 10,
  height = 8
)

library(lme4)
library(broom.mixed)
library(gt)

reg_data <- annual_plot %>%
  left_join(point_data, by = "locality") %>%
  mutate(
    lat_z = as.numeric(scale(lat)),
    lon_z = as.numeric(scale(lon))
  )

# ---- Mixed-effects models: spatial trends
models <- reg_data %>%
  group_by(variable) %>%
  group_map(
    ~ lmer(
      annual_value ~ lat_z + lon_z + (1 | locality),
      data = .x,
      REML = FALSE
    )
  )

names(models) <- levels(reg_data$variable)

# ---- Tidy fixed effects
lm_table <- bind_rows(
  lapply(models, tidy, effects = "fixed", conf.int = TRUE),
  .id = "variable"
)

gt_table <- lm_table %>%
  filter(term %in% c("lat_z", "lon_z")) %>%
  mutate(
    term = recode(
      term,
      lat_z = "Latitude (standardised)",
      lon_z = "Longitude (standardised)"
    )
  ) %>%
  select(variable, term, estimate, conf.low, conf.high) %>%
  gt(groupname_col = "variable") %>%
  fmt_number(
    columns = c(estimate, conf.low, conf.high),
    decimals = 3
  ) %>%
  cols_label(
    estimate = "Estimate",
    conf.low = "CI low",
    conf.high = "CI high"
  ) %>%
  tab_header(
    title = "Spatial trends in annual climate variables",
    subtitle = "Linear mixed-effects models (95% confidence intervals)"
  )

gt_table
# ------------------------------------------------------------
# 13. Spatial regressions: combined LAT + LON figure
# ------------------------------------------------------------

plot_regression <- function(varname, axis = c("lat", "lon")) {
  
  axis <- match.arg(axis)
  
  df <- reg_data %>%
    filter(variable == varname)
  
  x_var <- if (axis == "lat") "lat" else "lon"
  x_lab <- if (axis == "lat") "Latitude" else "Longitude"
  
  # ---- simplified Y-axis labels
  y_lab <- if (varname == "Average precipitation (mm yr⁻¹)") {
    "Precipitation (mm)"
  } else {
    "Temperature (°C)"
  }
  
  subtitle_txt <- dplyr::case_when(
    varname == "Minimum temperature (°C)" ~ "Annual minimum temperature",
    varname == "Maximum temperature (°C)" ~ "Annual maximum temperature",
    varname == "Average temperature (°C)" ~ "Annual mean temperature",
    varname == "Average precipitation (mm yr⁻¹)" ~ "Annual precipitation"
  )
  
  # ---- global regression for R²
  fml <- as.formula(paste("annual_value ~", x_var))
  lm_fit <- lm(fml, data = df)
  r2_lab <- paste0("R² = ", round(summary(lm_fit)$r.squared, 2))
  
  ggplot(df, aes_string(x = x_var, y = "annual_value", color = "locality")) +
    geom_point(alpha = 0.7, size = 2.5) +   # ⬅ bigger points
    
    # locality-specific regressions
    geom_smooth(
      aes(group = locality),
      method = "lm",
      se = FALSE,
      linetype = "dashed"
    ) +
    
    # global regression
    geom_smooth(
      aes(group = 1),
      method = "lm",
      se = TRUE,
      linetype = "dashed",
      color = "black"
    ) +
    
    annotate(
      "text",
      x = Inf, y = Inf,
      label = r2_lab,
      hjust = 1.1, vjust = 1.3,
      size = 5,
      color = "black"
    ) +
    
    scale_color_manual(
      values = palette_localities,
      limits = locality_levels,
      drop = FALSE
    ) +
    
    theme_bw(base_size = 18) +
    theme(
      axis.text  = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.subtitle = element_text(color = "black"),
      legend.position = "right"
    ) +
    
    labs(
      x = x_lab,
      y = y_lab,
      subtitle = subtitle_txt,
      color = "Locality"
    )
}

# ------------------------------------------------------------
# Combined regression figure: LAT (top) + LON (bottom)
# ------------------------------------------------------------

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
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 18, face = "bold", color = "black")
    )
  ) +
  plot_layout(guides = "collect")
figure_regressions
# ------------------------------------------------------------
# Export
# ------------------------------------------------------------

ggsave(
  "Figure_spatial_regressions_lat_lon.svg",
  figure_regressions,
  width = 24,
  height = 12,
  device = "pdf"
  
)
