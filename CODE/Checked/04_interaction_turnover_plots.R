# ============================================================
# 04_interaction_turnover_plots.R   >>> PLOTS ONLY <<<
# Fig. 1D  — WN dissimilarity: Par-Cat vs Cat-Plant (full data)
# Fig. 2   — WN / ST / OS: three groups incl. subsampled
# Fig. S6 — WN by dataset type, dodged violins
# Reads output/rds/interaction_turnover_main.rds and
#       output/rds/interaction_turnover_wn_by_type.rds
# from 04_interaction_turnover.R -> run that first.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)      # ggarrange()
library(ggsignif)    # geom_signif()
library(scales)      # alpha()

# ------------------------------------------------------------
# 1. Colours + theme  >>> EXACT COLOURS - DO NOT CHANGE <<<
# ------------------------------------------------------------
col_par_cat       <- "#4682B4"   # parasitoid-caterpillar interactions
col_cat_plant     <- "#B4464B"   # caterpillar-plant interactions
col_cat_plant_sub <- "#E5A7AB"   # caterpillar-plant subsampled (lighter)

fill_cols_fw <- c(
  "Parasitoid-Caterpillar"       = col_par_cat,
  "Caterpillar-Plant"            = col_cat_plant,
  "Caterpillar-Plant subsampled" = col_cat_plant_sub
)

# Fig. S7 palette: two networks x three dataset types
fill_cols_s7 <- c(
  "Caterpillar-Plant | full"                          = col_cat_plant,
  "Caterpillar-Plant | common plants only"            = "#e49d99",
  "Caterpillar-Plant | common caterpillars only"      = "#ffd6d3",
  "Parasitoid-Caterpillar | full"                     = col_par_cat,
  "Parasitoid-Caterpillar | common plants only"       = "#a3c8e9",
  "Parasitoid-Caterpillar | common caterpillars only" = "#a7f5ff"
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
    plot.background   = element_rect(fill = "white", color = NA),
    plot.tag          = element_text(size = 28, face = "bold", color = "black")
  )

# ------------------------------------------------------------
# 2. Load analysis output
# ------------------------------------------------------------
plot_data_main <- readRDS("output/rds/interaction_turnover_main.rds") %>%
  mutate(dataset = factor(dataset, levels = c(
    "Parasitoid-Caterpillar", "Caterpillar-Plant", "Caterpillar-Plant subsampled"
  )))

wn_by_type <- readRDS("output/rds/interaction_turnover_wn_by_type.rds")

dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Paired Wilcoxon signed-rank significance, so figure brackets
# match the paired tests in Table S5.
# ------------------------------------------------------------
wilcox_stars <- function(df, d1, d2) {
  a <- df %>% filter(dataset == d1) %>% dplyr::select(i, j, va = value)
  b <- df %>% filter(dataset == d2) %>% dplyr::select(i, j, vb = value)
  m <- dplyr::inner_join(a, b, by = c("i", "j"))
  if (nrow(m) < 2) return("NA")
  p <- wilcox.test(m$va, m$vb, paired = TRUE, exact = FALSE)$p.value
  if      (is.na(p))  "NA"
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else                "NS"
}

# ============================================================
# 3. FIGURE 1D — WN dissimilarity: Par-Cat vs Cat-Plant (full only)
# ============================================================
fig1d_data <- plot_data_main %>%
  filter(metric == "WN", type == "full") %>%
  mutate(dataset = factor(dataset, levels = c("Caterpillar-Plant", "Parasitoid-Caterpillar")))

Fig_1D <- ggplot(fig1d_data, aes(x = dataset, y = value, fill = dataset)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(
    aes(shape = dataset),
    position = position_jitter(0.2), size = 2,
    color = "black", fill = alpha("black", 0.4), stroke = 0.6, show.legend = FALSE
  ) +
  geom_signif(
    annotations = wilcox_stars(fig1d_data, "Caterpillar-Plant", "Parasitoid-Caterpillar"),
    xmin = 1, xmax = 2, y_position = 1.0,
    textsize = 6, color = "black", tip.length = 0.01
  ) +
  scale_fill_manual(values = c(
    "Parasitoid-Caterpillar" = col_par_cat,
    "Caterpillar-Plant"      = col_cat_plant
  )) +
  scale_shape_manual(values = c(
    "Caterpillar-Plant"      = 21,
    "Parasitoid-Caterpillar" = 24
  )) +
  coord_cartesian(ylim = c(0, 1.1)) +
  labs(x = "", y = "Overall interaction dissimilarity (WN)") +
  theme_pub_black +
  theme(legend.position = "none")

ggsave("output/fig/Figure_1D_WN_dissimilarity.svg", Fig_1D, width = 6, height = 8, units = "in")
ggsave("output/fig/Figure_1D_WN_dissimilarity.pdf", Fig_1D, width = 6, height = 8, units = "in")
ggsave("output/fig/Figure_1D_WN_dissimilarity.png", Fig_1D, width = 6, height = 8, units = "in", dpi = 300)

# ============================================================
# 4. FIGURE 2 — WN, ST, OS: all three groups incl. subsampled
# ============================================================
plot_fw_panel <- function(metric_name, y_label, tag_letter) {

  df <- plot_data_main %>%
    filter(metric == metric_name) %>%
    mutate(dataset = factor(dataset, levels = c(
      "Caterpillar-Plant", "Caterpillar-Plant subsampled", "Parasitoid-Caterpillar"
    )))

  ggplot(df, aes(x = dataset, y = value, fill = dataset)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
    geom_jitter(
      aes(shape = dataset),
      position = position_jitter(0.2), size = 2,
      color = "black", fill = alpha("black", 0.4), stroke = 0.6, show.legend = FALSE
    ) +
    geom_signif(annotations = wilcox_stars(df, "Caterpillar-Plant", "Caterpillar-Plant subsampled"),
                xmin = 1, xmax = 2, y_position = 0.88,
                textsize = 5, color = "black", tip.length = 0.01) +
    geom_signif(annotations = wilcox_stars(df, "Caterpillar-Plant subsampled", "Parasitoid-Caterpillar"),
                xmin = 2, xmax = 3, y_position = 0.96,
                textsize = 5, color = "black", tip.length = 0.01) +
    geom_signif(annotations = wilcox_stars(df, "Caterpillar-Plant", "Parasitoid-Caterpillar"),
                xmin = 1, xmax = 3, y_position = 1.06,
                textsize = 5, color = "black", tip.length = 0.01) +
    scale_fill_manual(values = fill_cols_fw) +
    scale_shape_manual(values = c(
      "Parasitoid-Caterpillar"       = 24,
      "Caterpillar-Plant"            = 21,
      "Caterpillar-Plant subsampled" = 22
    )) +
    scale_x_discrete(labels = c(
      "Caterpillar-Plant"            = "Caterpillar-\nPlant",
      "Caterpillar-Plant subsampled" = "Caterpillar-Plant\nsubsampled",
      "Parasitoid-Caterpillar"       = "Parasitoid-\nCaterpillar"
    )) +
    coord_cartesian(ylim = c(0, 1.15)) +
    labs(x = "", y = y_label, tag = tag_letter) +
    guides(fill = guide_legend(title = "Network")) +
    theme_pub_black +
    theme(legend.position = "bottom")
}

panel_WN <- plot_fw_panel("WN", "Overall dissimilarity in interactions",       "a")
panel_ST <- plot_fw_panel("ST", "Dissimilarity explained by species turnover", "b")
panel_OS <- plot_fw_panel("OS", "Dissimilarity explained by rewiring",         "c")

Fig_2 <- ggpubr::ggarrange(panel_WN, panel_ST, panel_OS,
                           common.legend = TRUE, legend = "bottom", nrow = 1)

ggsave("output/fig/Figure_2_network_metrics.svg", Fig_2, width = 16, height = 10, units = "in")
ggsave("output/fig/Figure_2_network_metrics.pdf", Fig_2, width = 16, height = 10, units = "in")
ggsave("output/fig/Figure_2_network_metrics.png", Fig_2, width = 16, height = 10, units = "in", dpi = 300)

# ============================================================
# 5. FIGURE S6 — WN by dataset type, dodged violins
# ============================================================
plot_wn_combined <- wn_by_type %>%
  mutate(
    network  = factor(network, levels = c("Caterpillar-Plant", "Parasitoid-Caterpillar")),
    type     = factor(type,    levels = c("full", "common plants only", "common caterpillars only")),
    fill_key = paste(network, type, sep = " | ")
  )

pd_s6 <- position_dodge(width = 0.8)

# Unpaired Wilcoxon among the three dataset types WITHIN each network
# (matches Fig. S6 caption + Table S9 unpaired variant). geom_signif does
# not support dodged groups, so brackets are built manually per dodge side.
x_ctr      <- c("full" = 1, "common plants only" = 2, "common caterpillars only" = 3)
net_offset <- c("Caterpillar-Plant" = -0.2, "Parasitoid-Caterpillar" = 0.2)
cmp_s6     <- list(c("full", "common plants only"),
                   c("common plants only", "common caterpillars only"),
                   c("full", "common caterpillars only"))
y_base     <- c("Caterpillar-Plant" = 0.92, "Parasitoid-Caterpillar" = 1.14)
y_step     <- 0.08

brackets_s6 <- purrr::imap_dfr(net_offset, function(off, net) {
  dfn <- plot_wn_combined %>% filter(network == net)
  purrr::imap_dfr(cmp_s6, function(pr, k) {
    xv <- dfn %>% filter(type == pr[1]) %>% pull(value)
    yv <- dfn %>% filter(type == pr[2]) %>% pull(value)
    wt <- wilcox.test(xv, yv, exact = FALSE)          # UNPAIRED (Mann-Whitney), matches Table S9
    tibble(
      network = net,
      x1 = unname(x_ctr[pr[1]] + off),
      x2 = unname(x_ctr[pr[2]] + off),
      y  = unname(y_base[net] + (k - 1) * y_step),
      label = case_when(wt$p.value < 0.001 ~ "***",
                        wt$p.value < 0.01  ~ "**",
                        wt$p.value < 0.05  ~ "*",
                        TRUE               ~ "ns")
    )
  })
})

Fig_S6 <- ggplot(plot_wn_combined,
                 aes(x = type, y = value, fill = fill_key,
                     group = interaction(type, network))) +
  geom_violin(position = pd_s6, trim = FALSE, width = 0.75, color = "black") +
  geom_boxplot(position = pd_s6, width = 0.12, outlier.shape = NA,
               fill = "white", color = "black") +
  geom_jitter(
    aes(shape = network),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2, color = "black", fill = alpha("black", 0.4), stroke = 0.6, show.legend = FALSE
  ) +
  geom_segment(data = brackets_s6, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_segment(data = brackets_s6, aes(x = x1, xend = x1, y = y - 0.015, yend = y),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_segment(data = brackets_s6, aes(x = x2, xend = x2, y = y - 0.015, yend = y),
               inherit.aes = FALSE, color = "black", linewidth = 0.5) +
  geom_text(data = brackets_s6, aes(x = (x1 + x2) / 2, y = y + 0.02, label = label),
            inherit.aes = FALSE, color = "black", size = 5) +
  scale_fill_manual(values = fill_cols_s7, name = "Network | Dataset", labels = names(fill_cols_s7)) +
  scale_shape_manual(values = c("Caterpillar-Plant" = 21, "Parasitoid-Caterpillar" = 24),
                     name = "Interaction group") +
  scale_x_discrete(labels = c(
    "full" = "Full", "common plants only" = "Common\nplants",
    "common caterpillars only" = "Common\ncaterpillars"
  )) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  coord_cartesian(ylim = c(0, 1.4)) +
  labs(x = "", y = "Dissimilarity between pairs of networks") +
  guides(shape = guide_legend(override.aes = list(size = 4)), fill = "none") +
  theme_pub_black +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("output/fig/Figure_S6_WN_by_dataset.svg", Fig_S6, width = 12, height = 8, units = "in")
ggsave("output/fig/Figure_S6_WN_by_dataset.pdf", Fig_S6, width = 12, height = 8, units = "in")
ggsave("output/fig/Figure_S6_WN_by_dataset.png", Fig_S6, width = 12, height = 8, units = "in", dpi = 300)
