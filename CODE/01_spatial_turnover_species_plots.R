# ============================================================
# 01_spatial_turnover_species_plots.R   >>> PLOTS ONLY <<<
# Builds Fig. S4 (three indices) and Fig. 1C (Bray-Curtis).
# Reads output/rds/turnover_long.rds produced by
# 01_spatial_turnover_species.R -> run that first.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)      # ggarrange(), stat_compare_means()
library(scales)      # alpha()

# ------------------------------------------------------------
# 1. Colour palette, shapes, theme  >>> EXACT COLOURS - DO NOT CHANGE <<<
# ------------------------------------------------------------
col_parasitoids      <- "#EDB439"
col_caterpillars     <- "#F25E0D"
col_caterpillars_sub <- "#F018DA"   # subsampled caterpillars

fill_cols <- c(
  "Parasitoids"             = col_parasitoids,
  "Caterpillars"            = col_caterpillars,
  "Caterpillars-subsampled" = col_caterpillars_sub
)
shape_vals <- c(
  "Caterpillars"            = 21,   # filled circle with outline
  "Caterpillars-subsampled" = 22,   # filled square with outline
  "Parasitoids"             = 24    # filled triangle with outline
)
x_labels <- c(
  "Caterpillars",
  "Caterpillars\nsubsampled",
  "Parasitoids"
)
desired_order <- c("Caterpillars", "Caterpillars-subsampled", "Parasitoids")

theme_pub_black <- theme_classic(base_size = 20) +
  theme(
    text                  = element_text(color = "black"),
    axis.line             = element_line(color = "black", linewidth = 0.8),
    axis.ticks            = element_line(color = "black", linewidth = 0.8),
    axis.text             = element_text(color = "black"),
    axis.title            = element_text(color = "black"),
    legend.text           = element_text(color = "black"),
    legend.title          = element_text(color = "black"),
    legend.box.background = element_rect(color = NA, fill = NA),
    plot.tag              = element_text(size = 28, face = "bold", color = "black")
  )

# ------------------------------------------------------------
# 2. Load analysis output
# ------------------------------------------------------------
turnover_long <- readRDS("output/rds/turnover_long.rds")
dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 3. FIGURE S4 — violin plots of all three indices
# ============================================================
figure_data <- turnover_long %>%
  mutate(guild = factor(guild, levels = desired_order))

my_comparisons <- list(
  c("Caterpillars", "Parasitoids"),
  c("Caterpillars", "Caterpillars-subsampled"),
  c("Parasitoids",  "Caterpillars-subsampled")
)

plot_beta_panel <- function(index_name, y_label, tag_letter, ylim_upper = 1.01) {
  
  df <- figure_data %>% filter(index == index_name)
  
  ggplot(df, aes(x = guild, y = dissimilarity)) +
    geom_violin(aes(fill = guild), trim = TRUE, color = "black") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
    geom_jitter(
      aes(shape = guild),
      position    = position_jitter(0.2),
      size        = 2,
      color       = "black",
      fill        = alpha("black", 0.4),
      stroke      = 0.6,
      show.legend = FALSE
    ) +
    stat_compare_means(
      comparisons = my_comparisons,
      method      = "wilcox.test",
      label       = "p.signif",
      color       = "black",
      tip.length  = 0.01,
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols   = c("***", "**", "*", "ns")
      )
    ) +
    scale_fill_manual(values = fill_cols) +
    scale_shape_manual(values = shape_vals) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = c(0, ylim_upper)) +
    labs(x = "", y = y_label, tag = tag_letter) +
    guides(
      fill  = guide_legend(title = "Community"),
      shape = guide_legend(title = "Community")
    ) +
    theme_pub_black
}

BC_FIG  <- plot_beta_panel("Bray-Curtis",   "Bray\u2013Curtis dissimilarity",   "a", 1.01)
CS_FIG  <- plot_beta_panel("Chao-Sorensen", "Chao\u2013S\u00f8rensen dissimilarity", "b", 1.00)
SOR_FIG <- plot_beta_panel("Sorensen",      "S\u00f8rensen dissimilarity",      "c", 1.00)

Fig_S4 <- ggpubr::ggarrange(
  BC_FIG, CS_FIG, SOR_FIG,
  common.legend = TRUE, legend = "bottom", nrow = 1
); Fig_S4

ggsave("output/fig/Figure_S4_beta_diversity.pdf", Fig_S4, width = 16, height = 10, units = "in")
ggsave("output/fig/Figure_S4_beta_diversity.png", Fig_S4, width = 16, height = 10, units = "in", dpi = 300)
ggsave("output/fig/Figure_S4_beta_diversity.svg", Fig_S4, width = 16, height = 10, units = "in")

# ============================================================
# 4. FIGURE 1C — Bray-Curtis only, caterpillars vs. parasitoids
# ============================================================
comparison_1c <- list(c("Caterpillars", "Parasitoids"))

fig1c_data <- turnover_long %>%
  filter(index == "Bray-Curtis",
         guild %in% c("Caterpillars", "Parasitoids")) %>%
  mutate(guild = factor(guild, levels = c("Caterpillars", "Parasitoids")))

Fig_1C <- ggplot(fig1c_data, aes(x = guild, y = dissimilarity, fill = guild)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(
    
    position    = position_jitter(0.2),
    size        = 5,
    color       = "black",
    fill        = alpha("grey", 0.4),
    stroke      = 0.3,
    show.legend = FALSE
  ) +
  stat_compare_means(
    comparisons = comparison_1c,
    method      = "wilcox.test",
    label       = "p.signif",
    color       = "black",
    tip.length  = 0.01,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols   = c("***", "**", "*", "ns")
    )
  ) +
  scale_fill_manual(values = c(
    "Caterpillars" = col_caterpillars,
    "Parasitoids"  = col_parasitoids
  )) +
  scale_shape_manual(values = c(
    "Caterpillars" = 21,
    "Parasitoids"  = 24
  )) +
  coord_cartesian(ylim = c(0, 1.01)) +
  labs(x = "", y = "Spatial turnover\nof species\n(Bray-Curtis dissimilarity)") +
  theme_pub_black +
  theme(legend.position = "none") ;Fig_1C

ggsave("output/fig/Figure_1C_BrayCurtis.pdf", Fig_1C, width = 6, height = 8, units = "in")
ggsave("output/fig/Figure_1C_BrayCurtis.png", Fig_1C, width = 6, height = 8, units = "in", dpi = 300)
ggsave("output/fig/Figure_1C_BrayCurtis.svg", Fig_1C, width = 6, height = 8, units = "in")
