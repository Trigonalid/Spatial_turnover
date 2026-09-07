# ============================================================
# SCRIPT: 04_enemy_free_space.R
#
# STUDY:
#   Approaching the limit of diversity: spatial turnover
#   amplifies with trophic level
#   Libra et al.
#
# ANALYSES:
#   1. Enemy-free space from all parasitoids combined:
#      proportion of occupied sites where caterpillar species
#      had no parasitoid recorded (Fig. S8A)
#   2. Enemy-free space from individual parasitoid species:
#      proportion of occupied sites where caterpillar species
#      was free from a specific parasitoid (Fig. S8B)
#   3. Summary statistics (mean, median) for both metrics
#   4. Export of raw interaction summary table
#
# KEY THRESHOLDS:
#   - Only caterpillar species parasitised at least once
#     across the entire study area are included (n = 127)
#   - Ohu2 excluded (temporal replicate)
#
# AUTHORS: Libra, M., Mottl, O.
# DATE:    [date]
# ============================================================


# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------

library(dplyr)       # data wrangling
library(tidyr)       # pivot functions
library(readr)       # read_delim
library(readxl)      # read_excel
library(here)        # here() for file paths
library(ggplot2)     # plotting
library(ggpubr)      # ggarrange()
library(writexl)     # write_xlsx()


# ------------------------------------------------------------
# 1. Global colour palette and theme
# ------------------------------------------------------------

col_parasitoids  <- "#009E73"
col_caterpillars <- "#0E2E71"

theme_pub_black <- theme_classic(base_size = 20) +
  theme(
    text                    = element_text(color = "black"),
    axis.line               = element_line(color = "black", linewidth = 0.8),
    axis.ticks              = element_line(color = "black", linewidth = 0.8),
    axis.text               = element_text(color = "black"),
    axis.title              = element_text(color = "black"),
    legend.text             = element_text(color = "black"),
    legend.title            = element_text(color = "black"),
    legend.background       = element_rect(fill = "white", color = NA),
    plot.background         = element_rect(fill = "white", color = NA),
    plot.tag                = element_text(size = 28, face = "bold", color = "black")
  )


# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

# Main dataset; Ohu2 excluded (temporal replicate, analysed separately)
MASTER <- read_excel(here("DATA/MASTER.xlsx")) %>%
  as_tibble() %>%
  dplyr::filter(locality != "Ohu2")

# Keep only caterpillar species that were parasitised at least once
# anywhere across the study area
vec_cat_with_par <- MASTER %>%
  group_by(CAT_sp) %>%
  summarise(
    has_parasitoid = any(!is.na(PAR_sp) & PAR_sp != ""),
    .groups = "drop"
  ) %>%
  filter(has_parasitoid == TRUE) %>%
  pull(CAT_sp)

# Working dataset: only caterpillars with at least one parasitoid record
data_work <- MASTER %>%
  filter(CAT_sp %in% vec_cat_with_par)


# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

# Fig. S8A: proportion of sites where caterpillar has NO parasitoid at all
# For each caterpillar species: (sites with cat but no par) / (total sites with cat)
get_summary_per_cat <- function(data_source) {
  data_source %>%
    group_by(CAT_sp) %>%
    summarise(
      .groups              = "drop",
      localities_total     = n_distinct(locality),
      localities_with_par  = n_distinct(locality[!is.na(PAR_sp)]),
      localities_no_par    = localities_total - localities_with_par,
      localities_prop      = localities_no_par / localities_total
    ) %>%
    select(CAT_sp, localities_total, localities_prop)
}

# Fig. S8B: proportion of sites where caterpillar is free from a SPECIFIC parasitoid
# For each caterpillar-parasitoid pair:
# (sites with cat but without that specific par) / (total sites with cat)
get_summary_per_pair <- function(data_source) {
  
  # Total localities per caterpillar species
  cat_localities <- data_source %>%
    group_by(CAT_sp) %>%
    summarise(localities_total_cat = n_distinct(locality), .groups = "drop")
  
  # Localities where each specific cat-par pair co-occurs
  pair_localities <- data_source %>%
    filter(!is.na(PAR_sp)) %>%
    group_by(CAT_sp, PAR_sp) %>%
    summarise(localities_with_pair = n_distinct(locality), .groups = "drop")
  
  # Join and compute proportion free from specific parasitoid
  pair_localities %>%
    left_join(cat_localities, by = "CAT_sp") %>%
    mutate(
      localities_prop = (localities_total_cat - localities_with_pair) /
        localities_total_cat
    ) %>%
    select(CAT_sp, PAR_sp, localities_prop)
}


# ------------------------------------------------------------
# 4. Run analyses
# ------------------------------------------------------------

summary_per_cat  <- get_summary_per_cat(data_work)
summary_per_pair <- get_summary_per_pair(data_work)


# ------------------------------------------------------------
# 5. Summary statistics
# ------------------------------------------------------------

cat("\n--- Enemy-free space: all parasitoids (Fig. S8A) ---\n")
summary_per_cat %>%
  summarise(
    mean   = round(mean(localities_prop,   na.rm = TRUE), 3),
    median = round(median(localities_prop, na.rm = TRUE), 3),
    n      = n()
  ) %>%
  print()

cat("\n--- Enemy-free space: per parasitoid species (Fig. S8B) ---\n")
summary_per_pair %>%
  summarise(
    mean   = round(mean(localities_prop,   na.rm = TRUE), 3),
    median = round(median(localities_prop, na.rm = TRUE), 3),
    n      = n()
  ) %>%
  print()


# ------------------------------------------------------------
# 6. Export raw interaction summary table
# ------------------------------------------------------------

# Full cat-par interaction table with locality counts
tab_efs_raw <- MASTER %>%
  filter(!is.na(PAR_sp)) %>%
  group_by(CAT_sp, PAR_sp) %>%
  summarise(localities_total = n_distinct(locality), .groups = "drop")

write_xlsx(tab_efs_raw,      here("output/Table_S_EFS_raw.xlsx"))
write_xlsx(summary_per_cat,  here("output/Table_S_EFS_per_cat.xlsx"))
write_xlsx(summary_per_pair, here("output/Table_S_EFS_per_pair.xlsx"))


# ------------------------------------------------------------
# 7. FIGURE S8 — Enemy-free space violin plots
# A: proportion of sites free from all parasitoids
# B: proportion of sites free from specific parasitoid
# Points: black outline, 40% transparent black fill.
# ------------------------------------------------------------

set.seed(1234)

plot_efs_panel <- function(data_source, fill_color, y_label) {
  ggplot(data_source, aes(x = "group", y = localities_prop)) +
    geom_violin(trim = TRUE, fill = fill_color, color = "black") +
    geom_boxplot(width = 0.15, fill = "white", color = "black",
                 outlier.shape = NA) +
    geom_jitter(
      aes(fill = I(alpha("black", 0.4))),
      shape    = 21,
      width    = 0.15,
      height   = 0,
      size     = 2,
      color    = "black",
      stroke   = 0.6
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_discrete(labels = c("")) +
    labs(x = "", y = y_label) +
    theme_pub_black +
    theme(legend.position = "none")
}

Fig_S8A <- plot_efs_panel(
  summary_per_cat,
  fill_color = col_caterpillars,
  y_label    = "Proportion of sites without any parasitoid"
)

Fig_S8B <- plot_efs_panel(
  summary_per_pair,
  fill_color = col_parasitoids,
  y_label    = "Proportion of sites without\na particular parasitoid"
)

Fig_S8 <- ggpubr::ggarrange(
  Fig_S8A, Fig_S8B,
  ncol   = 2,
  labels = c("a", "b"),
  hjust  = -0.2
)

ggsave("output/fig/Figure_S8_enemy_free_space.svg",
       Fig_S8, width = 10, height = 7, units = "in")
ggsave("output/fig/Figure_S8_enemy_free_space.pdf",
       Fig_S8, width = 10, height = 7, units = "in")
ggsave("output/fig/Figure_S8_enemy_free_space.png",
       Fig_S8, width = 10, height = 7, units = "in", dpi = 300)