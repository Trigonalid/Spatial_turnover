# ============================================================
# 02_rarity_parasitoids_plots.R   >>> PLOTS ONLY <<<
# Fig. S5 — dissimilarity by parasitoid rarity group,
# faceted by index (columns) and dataset (rows).
# Reads output/rds/rarity_long.rds from 02_rarity_parasitoids.R -> run that first.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggsignif)    # geom_signif() for significance brackets in facets
library(scales)      # alpha()

# ------------------------------------------------------------
# 1. Colours + theme  >>> EXACT COLOURS - DO NOT CHANGE <<<
# ------------------------------------------------------------
col_all    <- "#EDB439"   # all parasitoids
col_common <- "#772CA0"   # common parasitoids - purple
col_rare   <- "#4682B4"   # rare parasitoids - light blue

fill_cols_rarity <- c("All" = col_all, "Common" = col_common, "Rare" = col_rare)
shape_vals_rarity <- c("All" = 21, "Common" = 22, "Rare" = 24)

theme_pub_black <- theme_classic(base_size = 20) +
  theme(
    text              = element_text(color = "black"),
    axis.line         = element_line(color = "black", linewidth = 0.8),
    axis.ticks        = element_line(color = "black", linewidth = 0.8),
    axis.text         = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    legend.text       = element_text(color = "black"),
    legend.title      = element_text(color = "black"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.tag          = element_text(size = 28, face = "bold", color = "black")
  )

# ------------------------------------------------------------
# 2. Load analysis output + set factor levels / labels
# ------------------------------------------------------------
plot_data <- readRDS("output/rds/rarity_long.rds") %>%
  mutate(
    guild_clean = factor(guild_clean, levels = c("All", "Common", "Rare")),
    dataset     = factor(dataset, levels = c("All data",
                                             "Common Plants Only",
                                             "Common Caterpillars Only")),
    index_type  = factor(index_type,
                         levels = c("bray", "chao-sorensen", "sorensen"),
                         labels = c("Bray\u2013Curtis", "Chao\u2013S\u00f8rensen", "S\u00f8rensen"))
  )

dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 3. Figure S5
# ------------------------------------------------------------
Fig_S5 <- ggplot(plot_data,
                 aes(x = guild_clean, y = dissimilarity, fill = guild_clean)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(
    aes(shape = guild_clean),
    position    = position_jitter(0.2),
    size        = 2,
    color       = "black",
    fill        = alpha("black", 0.4),   # unified with 01 (outside aes)
    stroke      = 0.6,
    show.legend = FALSE
  ) +
  geom_signif(comparisons = list(c("All", "Common")),
              test = "wilcox.test", map_signif_level = TRUE,
              textsize = 5, color = "black", tip.length = 0.01, y_position = 0.92) +
  geom_signif(comparisons = list(c("Common", "Rare")),
              test = "wilcox.test", map_signif_level = TRUE,
              textsize = 5, color = "black", tip.length = 0.01, y_position = 1.00) +
  geom_signif(comparisons = list(c("All", "Rare")),
              test = "wilcox.test", map_signif_level = TRUE,
              textsize = 5, color = "black", tip.length = 0.01, y_position = 1.08) +
  scale_fill_manual(values = fill_cols_rarity) +
  scale_shape_manual(values = shape_vals_rarity) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.2)) +
  facet_grid(dataset ~ index_type) +
  labs(x = "Parasitoid group", y = "Dissimilarity", fill = "Parasitoid group") +
  guides(fill = guide_legend(title = "Parasitoid group")) +
  theme_pub_black +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    strip.text      = element_text(size = 14),
    legend.position = "bottom"
  )

ggsave("output/fig/Figure_S5_rarity.svg", Fig_S5, width = 14, height = 14, units = "in")
ggsave("output/fig/Figure_S5_rarity.pdf", Fig_S5, width = 14, height = 14, units = "in")
ggsave("output/fig/Figure_S5_rarity.png", Fig_S5, width = 14, height = 14, units = "in", dpi = 300)
