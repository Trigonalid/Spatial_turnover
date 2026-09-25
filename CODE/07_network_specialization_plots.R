# ============================================================
# 07_network_specialization_plots.R   >>> PLOTS ONLY <
# Fig. S7 — network-level specialization H2' per locality
#   (paired Wilcoxon, par-cat vs cat-plant).
# Single panel. Main = min5_cons; robustness = full, min5_both.
# Reads output/rds/spec_H2_local.rds from 07_network_specialization.R
# -> run that first.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)    # geom_signif()
library(scales)      # alpha()

# ------------------------------------------------------------
# 1. Colours + theme  >>> EXACT COLOURS - DO NOT CHANGE <<<
# ------------------------------------------------------------
col_par_cat   <- "#4682B4"
col_cat_plant <- "#B4464B"
fill_cols_fw <- c(
  "Parasitoid-Caterpillar" = col_par_cat,
  "Caterpillar-Plant"      = col_cat_plant
)

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
    plot.background   = element_rect(fill = "white", color = NA)
  )

# ------------------------------------------------------------
# 2. Load analysis output
# ------------------------------------------------------------
H2_local_all <- readRDS("output/rds/spec_H2_local.rds")
dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 3. H2' panel (per locality); paired Wilcoxon significance
# ------------------------------------------------------------
# Paired stars (par-cat vs cat-plant), matched by locality, to match Table S6
paired_stars <- function(df) {
  wide <- df %>% pivot_wider(names_from = network, values_from = H2)
  p <- tryCatch(
    wilcox.test(wide$`Parasitoid-Caterpillar`, wide$`Caterpillar-Plant`,
                paired = TRUE, exact = FALSE)$p.value,
    error = function(e) NA_real_)
  if      (is.na(p))  "NA"
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else                "NS"
}

make_panel_H2 <- function(df_H2) {
  df_plot <- df_H2 %>%
    mutate(network = factor(network, levels = c("Caterpillar-Plant", "Parasitoid-Caterpillar")))
  
  ggplot(df_plot, aes(x = network, y = H2, fill = network)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
    geom_jitter(aes(shape = network),
                position = position_jitter(0.15), size = 2,
                color = "black", fill = alpha("black", 0.4), stroke = 0.6, show.legend = FALSE) +
    geom_signif(annotations = paired_stars(df_plot),
                xmin = 1, xmax = 2, y_position = 0.98,
                textsize = 6, color = "black", tip.length = 0.01) +
    scale_fill_manual(values = fill_cols_fw) +
    scale_shape_manual(values = c("Caterpillar-Plant" = 21, "Parasitoid-Caterpillar" = 24)) +
    scale_x_discrete(labels = c("Caterpillar-Plant" = "Caterpillar-\nPlant",
                                "Parasitoid-Caterpillar" = "Parasitoid-\nCaterpillar")) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(x = "", y = expression("Network specialization " ~ italic(H[2]*"'"))) +
    theme_pub_black +
    theme(legend.position = "none")
}

# ------------------------------------------------------------
# 4. Save H2' figure per variant
#    main = >=5 reared specimens (Fig. S7); all species = robustness
# ------------------------------------------------------------
save_H2 <- function(variant_name, file_name) {
  p <- make_panel_H2(H2_local_all %>% filter(variant == variant_name))
  ggsave(paste0("output/fig/", file_name, ".svg"), p, width = 6, height = 8, units = "in")
  ggsave(paste0("output/fig/", file_name, ".pdf"), p, width = 6, height = 8, units = "in")
  ggsave(paste0("output/fig/", file_name, ".png"), p, width = 6, height = 8, units = "in", dpi = 300)
  p
}

save_H2("\u22655 reared specimens (parasitoids & caterpillars)", "Figure_S7_specialization")            # MAIN (Fig. S7)
save_H2("All species (no threshold)",                            "Figure_S7_specialization_nothreshold") # robustness