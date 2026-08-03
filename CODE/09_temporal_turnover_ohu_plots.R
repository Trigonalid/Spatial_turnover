# ============================================================
# 09_temporal_turnover_ohu_plots.R   >>> PLOTS ONLY <<<
# Fig. S6A — Bray-Curtis temporal dissimilarity (caterpillars, parasitoids)
# Fig. S6B — WN temporal dissimilarity (par-cat, cat-plant, cat-plant subsampled)
# Reads output/rds/temporal_bray.rds and output/rds/temporal_wn.rds
# from 09_temporal_turnover_ohu.R -> run that first.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(grid)        # unit() for theme sizing
library(svglite)     # clean SVG export (Affinity compatible)
library(patchwork)   # combining panels
library(scales)      # alpha()

# ------------------------------------------------------------
# 1. Colours + theme  >>> EXACT COLOURS - DO NOT CHANGE <<<
# ------------------------------------------------------------
col_parasitoids  <- "#EDB439"
col_caterpillars <- "#F25E0D"
col_par_cat      <- "#4682B4"
col_cat_plant    <- "#B4464B"

theme_pub_black <- theme_classic(base_size = 18) +
  theme(
    text              = element_text(color = "black"),
    axis.line         = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    axis.ticks        = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(4, "pt"),
    axis.text         = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    legend.text       = element_text(color = "black"),
    legend.title      = element_text(color = "black"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.tag          = element_text(size = 28, face = "bold", color = "black")
  )

# ------------------------------------------------------------
# 2. Load analysis output
# ------------------------------------------------------------
bray_data_all <- readRDS("output/rds/temporal_bray.rds")
plot_data_wn  <- readRDS("output/rds/temporal_wn.rds")

dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 3. FIGURE S6A — Bray-Curtis temporal dissimilarity
# ------------------------------------------------------------
Fig_S6A <- ggplot(bray_data_all, aes(x = x_group, y = value)) +
  geom_violin(
    data  = filter(bray_data_all, type == "Subsampled"),
    aes(fill = guild), trim = FALSE, color = "black", width = 0.7
  ) +
  geom_boxplot(
    data = filter(bray_data_all, type == "Subsampled"),
    width = 0.12, outlier.shape = NA, fill = "white", color = "black"
  ) +
  geom_jitter(
    data = filter(bray_data_all, type == "Subsampled"),
    shape = 21, fill = alpha("black", 0.15),
    width = 0.1, size = 1.5, color = alpha("black", 0.3), stroke = 0.3
  ) +
  geom_point(
    data = filter(bray_data_all, type == "Observed"),
    aes(fill = guild), size = 5, shape = 21, color = "black", stroke = 1.5
  ) +
  scale_fill_manual(values = c(
    "Caterpillars" = col_caterpillars,
    "Parasitoids"  = col_parasitoids
  )) +
  scale_x_discrete(limits = c(
    "Caterpillars\nObserved", "Caterpillars\nSubsampled",
    "Parasitoids\nObserved",  "Parasitoids\nSubsampled"
  )) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Bray\u2013Curtis dissimilarity", tag = "A") +
  guides(fill = "none") +
  theme_pub_black +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------------------------------------
# 4. FIGURE S6B — WN temporal dissimilarity
# ------------------------------------------------------------
Fig_S6B <- ggplot(plot_data_wn, aes(x = dataset, y = WN)) +
  geom_violin(
    data = filter(plot_data_wn, dataset == "Caterpillar\u2013Plant (subsampled)"),
    fill = col_cat_plant, color = "black", trim = FALSE, width = 0.6
  ) +
  geom_boxplot(
    data = filter(plot_data_wn, dataset == "Caterpillar\u2013Plant (subsampled)"),
    width = 0.12, outlier.shape = NA, fill = "white", color = "black"
  ) +
  geom_jitter(
    data = filter(plot_data_wn, dataset == "Caterpillar\u2013Plant (subsampled)"),
    shape = 21, fill = alpha("black", 0.15),
    width = 0.12, size = 1.5, color = alpha("black", 0.3), stroke = 0.3
  ) +
  geom_point(
    data = filter(plot_data_wn, dataset != "Caterpillar\u2013Plant (subsampled)"),
    aes(fill = dataset), size = 5, shape = 21, color = "black", stroke = 1.5
  ) +
  scale_fill_manual(values = c(
    "Caterpillar\u2013Plant"      = col_cat_plant,
    "Parasitoid\u2013Caterpillar" = col_par_cat
  )) +
  scale_x_discrete(
    limits = c(
      "Caterpillar\u2013Plant",
      "Caterpillar\u2013Plant (subsampled)",
      "Parasitoid\u2013Caterpillar"
    ),
    labels = c(
      "Caterpillar\u2013Plant",
      "Caterpillar\u2013Plant\n(subsampled)",
      "Parasitoid\u2013Caterpillar"
    )
  ) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(x = NULL, y = "Overall dissimilarity (WN)", tag = "B") +
  guides(fill = "none") +
  theme_pub_black +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------------------------------------
# 5. Combined Figure S6
# ------------------------------------------------------------
Fig_S6 <- Fig_S6A + Fig_S6B + plot_layout(ncol = 2)

ggsave("output/fig/Figure_S6_temporal_Ohu.svg", Fig_S6,
       width = 25, height = 15, units = "cm", device = svglite)
ggsave("output/fig/Figure_S6_temporal_Ohu.pdf", Fig_S6,
       width = 25, height = 15, units = "cm")
ggsave("output/fig/Figure_S6_temporal_Ohu.png", Fig_S6,
       width = 25, height = 15, units = "cm", dpi = 300)
