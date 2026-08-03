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
# 4. Spatial heatmaps
# ------------------------------------------------------------
p_hosts <- heatmap_panel(results, "host_SOR", ramp_hosts, "Spatial Turnover",
                         "Spatial species turnover\n- Resources", base = 24)
p_paras <- heatmap_panel(results, "para_SOR", ramp_paras, "Spatial Turnover",
                         "Spatial species turnover\n- Consumers")
p_diff  <- heatmap_panel(results, "diff_SOR", ramp_diff,
                         "delta Spatial Turnover\n(consumer - resource)",
                         "Spatial species turnover\n- Consumer vs Resource")

ggsave("output/fig/heatmap_hosts_spatial.png",       p_hosts, width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_parasitoids_spatial.png", p_paras, width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_diff_spatial.png",        p_diff,  width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_hosts_spatial.svg",       p_hosts, width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_parasitoids_spatial.svg", p_paras, width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_diff_spatial.svg",        p_diff,  width = 7, height = 6, dpi = 300)

# ------------------------------------------------------------
# 5. Temporal heatmaps
# ------------------------------------------------------------
p_hosts_temp <- heatmap_panel(results_temp, "host_SOR_temp", ramp_hosts, "Temporal Turnover",
                              "Temporal species turnover\n- Resources", base = 16)
p_paras_temp <- heatmap_panel(results_temp, "para_SOR_temp", ramp_paras, "Temporal Turnover",
                              "Temporal species turnover\n- Consumers", base = 16)
p_diff_temp  <- heatmap_panel(results_temp, "diff_SOR_temp", ramp_diff,
                              "delta Temporal Turnover\n(consumer - resource)",
                              "Temporal species turnover\n- Consumer vs Resource")

ggsave("output/fig/heatmap_hosts_temporal.png",       p_hosts_temp, width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_parasitoids_temporal.png", p_paras_temp, width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_diff_temporal.png",        p_diff_temp,  width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_hosts_temporal.svg",       p_hosts_temp, width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_parasitoids_temporal.svg", p_paras_temp, width = 7, height = 6, dpi = 300)
ggsave("output/fig/heatmap_diff_temporal.svg",        p_diff_temp,  width = 7, height = 6, dpi = 300)

combined_temporal <- plot_grid(p_hosts_temp, p_paras_temp, p_diff_temp, ncol = 3)
ggsave("output/fig/heatmap_temporal_combined.png", combined_temporal, width = 18, height = 7, dpi = 300)
ggsave("output/fig/heatmap_temporal_combined.svg", combined_temporal, width = 18, height = 7, dpi = 300)

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

# ------------------------------------------------------------
# 8. Final 2x3 grid: spatial heatmaps (top) + scatters (bottom)
# (consolidated from earlier duplicated blocks; legends removed)
#
# !!! PLACEHOLDER: the bottom-left panel currently duplicates the
#     Resources scatter. Replace 'p_scatter_placeholder' with the
#     intended third panel before using this figure in the paper.
# ------------------------------------------------------------
unify_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 22, colour = "black"),
  axis.title = element_text(size = 20, colour = "black"),
  axis.text  = element_text(size = 18, colour = "black"),
  legend.position = "none"
)

p_scatter_placeholder <- p_scatter_hosts   # <-- TODO: replace with real panel

top_row    <- plot_grid(p_hosts + unify_theme, p_paras + unify_theme,
                        p_diff + unify_theme, ncol = 3)
bottom_row <- plot_grid(p_scatter_placeholder + unify_theme,
                        p_scatter_hosts + unify_theme,
                        p_scatter_paras + unify_theme, ncol = 3)
combined_final <- plot_grid(top_row, bottom_row, nrow = 2)

ggsave("output/fig/combined_figure_spatial_scatter.png", combined_final, width = 18, height = 12, dpi = 300)
ggsave("output/fig/combined_figure_spatial_scatter.svg", combined_final, width = 18, height = 12, dpi = 300)
