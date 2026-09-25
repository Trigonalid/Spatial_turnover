# ============================================================
# 06_glmm_parasitism_plots.R   >>> PLOTS ONLY <<<
# Parasitism rate vs number of attacking parasitoid species.
# Points coloured by locality (shows the relationship holds across all
# 8 sites, not driven by one), size = N_total, black GLMM-predicted curve.
# Reads output/rds/glmm_*.rds from 06_glmm_parasitism.R -> run that first.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Colours + theme  >>> EXACT COLOURS - DO NOT CHANGE <<<
# ------------------------------------------------------------
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999")   # locality palette

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
plot_df <- readRDS("output/rds/glmm_plot_df.rds")
pred_df <- readRDS("output/rds/glmm_pred_df.rds")

dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

locality_levels <- sort(unique(plot_df$locality))
locality_cols   <- setNames(okabe_ito[seq_along(locality_levels)], locality_levels)
plot_df <- plot_df %>% mutate(locality = factor(locality, levels = locality_levels))

# ------------------------------------------------------------
# 3. Figure: parasitism rate vs N_par_species (points by locality)
# ------------------------------------------------------------
fig_pd <- ggplot(plot_df, aes(x = N_par_species, y = parasitism_rate)) +
  geom_jitter(aes(size = N_total, fill = locality),
              width = 0.15, height = 0, alpha = 0.80,
              colour = "black", shape = 21, stroke = 0.3) +
  geom_line(data = pred_df, aes(x = N_par_species, y = pred),
            inherit.aes = FALSE, colour = "black", linewidth = 1.1) +
  scale_size_continuous(name = expression(italic(N)[total]), range = c(1.5, 7),
                        breaks = c(10, 25, 50, 100, 200)) +
  scale_fill_manual(name = "Locality", values = locality_cols) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 5)),
         size = guide_legend(order = 2)) +
  scale_x_continuous(breaks = seq(0, max(plot_df$N_par_species), by = 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of parasitoid species attacking the caterpillar",
       y = "Parasitism rate") +
  theme_pub_black + theme(legend.position = "right")

ggsave("output/fig/Figure_S10_parasitism_diversity.svg", fig_pd, width = 11, height = 8, units = "in")
ggsave("output/fig/Figure_S10_parasitism_diversity.pdf", fig_pd, width = 11, height = 8, units = "in")
ggsave("output/fig/Figure_S10_parasitism_diversity.png", fig_pd, width = 11, height = 8, units = "in", dpi = 300)