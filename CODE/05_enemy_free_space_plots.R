# ============================================================
# 05_enemy_free_space_plots.R   >>> PLOTS ONLY <<<
# Fig. S8A — proportion of sites free from ALL parasitoids (per caterpillar)
# Fig. S8B — proportion of sites free from a SPECIFIC parasitoid (per pair)
# Reads output/rds/efs_per_cat.rds and output/rds/efs_per_pair.rds
# from 05_enemy_free_space.R -> run that first.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)      # ggarrange()
library(scales)      # alpha()

# ------------------------------------------------------------
# 1. Colours + theme  >>> EXACT COLOURS - DO NOT CHANGE <<<
# NOTE: these differ from the guild palette used in Fig. 1/S4/S6
#       (there parasitoids = #EDB439, caterpillars = #F25E0D).
#       Confirm this separate scheme for Fig. S8 is intended.
# ------------------------------------------------------------
col_parasitoids  <- "#009E73"
col_caterpillars <- "#0E2E71"

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
# 2. Load analysis output
# ------------------------------------------------------------
summary_per_cat  <- readRDS("output/rds/efs_per_cat.rds")
summary_per_pair <- readRDS("output/rds/efs_per_pair.rds")

dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 3. Figure S8
# ------------------------------------------------------------
set.seed(1234)   # reproducible jitter

plot_efs_panel <- function(data_source, fill_color, y_label) {
  ggplot(data_source, aes(x = "group", y = localities_prop)) +
    geom_violin(trim = TRUE, fill = fill_color, color = "black") +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
    geom_jitter(
      shape = 21, fill = alpha("black", 0.4),
      width = 0.15, height = 0, size = 2, color = "black", stroke = 0.6
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_discrete(labels = c("")) +
    labs(x = "", y = y_label) +
    theme_pub_black +
    theme(legend.position = "none")
}

Fig_S8A <- plot_efs_panel(summary_per_cat,  col_caterpillars,
                          "Proportion of sites without any parasitoid")
Fig_S8B <- plot_efs_panel(summary_per_pair, col_parasitoids,
                          "Proportion of sites without\na particular parasitoid")

Fig_S8 <- ggpubr::ggarrange(Fig_S8A, Fig_S8B, ncol = 2,
                            labels = c("A", "B"), hjust = -0.2)

ggsave("output/fig/Figure_S8_enemy_free_space.svg", Fig_S8, width = 10, height = 7, units = "in")
ggsave("output/fig/Figure_S8_enemy_free_space.pdf", Fig_S8, width = 10, height = 7, units = "in")
ggsave("output/fig/Figure_S8_enemy_free_space.png", Fig_S8, width = 10, height = 7, units = "in", dpi = 300)
