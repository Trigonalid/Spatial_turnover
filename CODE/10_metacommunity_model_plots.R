# ============================================================
# 10_metacommunity_model_plots.R   >>> PLOTS ONLY <<<
# Heatmaps (spatial + temporal Sorensen turnover across the colonisation
# grid) and spatial-vs-temporal scatter plots for the metacommunity model.
# Reads output/rds/model_beta_spatial.rds and model_beta_temporal.rds
# from 10_metacommunity_model.R -> run that first.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(ggplot2)
library(cowplot)     # get_legend(), plot_grid(), ggdraw(), draw_plot()

# ------------------------------------------------------------
# 1. Colour ramps  >>> EXACT COLOURS - DO NOT CHANGE <<<
# ------------------------------------------------------------
ramp_hosts <- c("#FDEBD0", "#F0A97A", "#E9711A", "#B85414", "#7D380D")  # resources
ramp_paras <- c("#FEF9E7", "#F9E4A0", "#EDB439", "#B8880A", "#7D5C07")  # consumers
ramp_diff  <- c("white",   "#92C5DE", "#4393C3", "#2166AC", "#053061")  # difference
pt_hosts   <- "#E9711A"   # resources scatter points
pt_paras   <- "#EDB439"   # consumers scatter points

# ------------------------------------------------------------
# 2. Load analysis output
# ------------------------------------------------------------
results      <- readRDS("output/rds/model_beta_spatial.rds")
results_temp <- readRDS("output/rds/model_beta_temporal.rds")

dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 3. Heatmap builder
# NOTE: fill scale limits are c(0, 1); for the difference panels any
# negative value (resource turnover > consumer) falls outside and is
# drawn as na.value grey - same grey as truly empty cells. Check whether
# negative differences should be shown distinctly.
# ------------------------------------------------------------
heatmap_panel <- function(df, fill_var, ramp, legend_name, title, base = 22) {
  ggplot(df, aes(x = para_col_rate, y = host_col_rate, fill = .data[[fill_var]])) +
    geom_tile(colour = "white", size = 0.1) +
    scale_fill_stepsn(
      colours = ramp, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1),
      name = legend_name, na.value = "grey80",
      guide = guide_coloursteps(even.steps = TRUE, show.limits = TRUE)) +
    scale_x_continuous(breaks = seq(0.2, 2.5, by = 0.4)) +
    scale_y_continuous(breaks = seq(0.2, 2.5, by = 0.4)) +
    labs(title = title, x = "Consumer colonisation rate", y = "Resource colonisation rate") +
    theme_bw(base_size = base) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
}

# ------------------------------------------------------------
# 4. Spatial heatmaps — Fig. 3B (resources), 3C (consumers), 3D (difference)
# ------------------------------------------------------------
p_hosts <- heatmap_panel(results, "host_SOR", ramp_hosts, "Spatial Turnover",
                         "Spatial species turnover\n- Resources", base = 24)
p_paras <- heatmap_panel(results, "para_SOR", ramp_paras, "Spatial Turnover",
                         "Spatial species turnover\n- Consumers")
p_diff  <- heatmap_panel(results, "diff_SOR", ramp_diff,
                         "delta Spatial Turnover\n(consumer - resource)",
                         "Spatial species turnover\n- Consumer vs Resource")

ggsave("output/fig/Figure_3B_spatial_turnover_resources.png",  p_hosts, width = 7, height = 6, dpi = 300)
ggsave("output/fig/Figure_3C_spatial_turnover_consumers.png",  p_paras, width = 7, height = 6, dpi = 300)
ggsave("output/fig/Figure_3D_spatial_turnover_difference.png", p_diff,  width = 7, height = 6, dpi = 300)
ggsave("output/fig/Figure_3B_spatial_turnover_resources.svg",  p_hosts, width = 7, height = 6, dpi = 300)
ggsave("output/fig/Figure_3C_spatial_turnover_consumers.svg",  p_paras, width = 7, height = 6, dpi = 300)
ggsave("output/fig/Figure_3D_spatial_turnover_difference.svg", p_diff,  width = 7, height = 6, dpi = 300)

# ------------------------------------------------------------
# 6. Spatial vs temporal scatter — Fig. 3F (resources), 3G (consumers)
# ------------------------------------------------------------
combined <- merge(results, results_temp, by = c("file", "host_col_rate", "para_col_rate"))

scatter_panel <- function(x_var, y_var, point_col, title, base = 22) {
  ggplot(combined, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(colour = point_col, alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", colour = "black", se = TRUE) +
    labs(title = title, x = "Spatial Turnover (t=100)",
         y = "Temporal Turnover (t=50 vs t=100)") +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw(base_size = base) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 22, hjust = 0.5),
          axis.title = element_text(size = 22), axis.text = element_text(size = 20))
}
p_scatter_hosts <- scatter_panel("host_SOR", "host_SOR_temp", pt_hosts,
                                 "Spatial vs Temporal turnover - Resources")
p_scatter_paras <- scatter_panel("para_SOR", "para_SOR_temp", pt_paras,
                                 "Spatial vs Temporal turnover - Consumers", base = 24)

ggsave("output/fig/Figure_3F_spatial_vs_temporal_resources.png", p_scatter_hosts, width = 6, height = 6, dpi = 300)
ggsave("output/fig/Figure_3G_spatial_vs_temporal_consumers.png", p_scatter_paras, width = 6, height = 6, dpi = 300)
ggsave("output/fig/Figure_3F_spatial_vs_temporal_resources.svg", p_scatter_hosts, width = 6, height = 6, dpi = 300)
ggsave("output/fig/Figure_3G_spatial_vs_temporal_consumers.svg", p_scatter_paras, width = 6, height = 6, dpi = 300)

# ------------------------------------------------------------
# 6. Spatial vs temporal scatter (Resources, Consumers)
# ------------------------------------------------------------
combined <- merge(results, results_temp, by = c("file", "host_col_rate", "para_col_rate"))

scatter_panel <- function(x_var, y_var, point_col, title, base = 22) {
  ggplot(combined, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(colour = point_col, alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", colour = "black", se = TRUE) +
    labs(title = title, x = "Spatial Turnover (t=100)",
         y = "Temporal Turnover (t=50 vs t=100)") +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw(base_size = base) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 22, hjust = 0.5),
          axis.title = element_text(size = 22), axis.text = element_text(size = 20))
}
p_scatter_hosts <- scatter_panel("host_SOR", "host_SOR_temp", pt_hosts,
                                 "Spatial vs Temporal turnover - Resources")
p_scatter_paras <- scatter_panel("para_SOR", "para_SOR_temp", pt_paras,
                                 "Spatial vs Temporal turnover - Consumers", base = 24)

ggsave("output/fig/scatter_hosts_spatial_temporal.png", p_scatter_hosts, width = 6, height = 6, dpi = 300)
ggsave("output/fig/scatter_paras_spatial_temporal.png", p_scatter_paras, width = 6, height = 6, dpi = 300)
ggsave("output/fig/scatter_hosts_spatial_temporal.svg", p_scatter_hosts, width = 6, height = 6, dpi = 300)
ggsave("output/fig/scatter_paras_spatial_temporal.svg", p_scatter_paras, width = 6, height = 6, dpi = 300)

# ------------------------------------------------------------
# 7. Combined spatial heatmaps (each with its own legend below)
# ------------------------------------------------------------
attach_legend <- function(p) {
  plot_grid(p + theme(legend.position = "none"), get_legend(p),
            ncol = 1, rel_heights = c(1, 0.12))
}
combined_spatial <- plot_grid(
  attach_legend(p_hosts), attach_legend(p_paras), attach_legend(p_diff), ncol = 3
)
ggsave("output/fig/heatmap_spatial_combined.png", combined_spatial, width = 18, height = 7, dpi = 300)
ggsave("output/fig/heatmap_spatial_combined.svg", combined_spatial, width = 18, height = 7, dpi = 300)

combined_spatial_padded <- ggdraw() +
  draw_plot(combined_spatial, x = 0.1, y = 0.15, width = 0.8, height = 0.7)
ggsave("output/fig/heatmap_spatial_combined_padded.svg",
       combined_spatial_padded, width = 24, height = 10, dpi = 300, bg = "white")

