# ============================================================
# 06b_parasitism_heatmap.R   >>> STANDALONE (analysis + figure) <<<
# Study: Spatial turnover amplifies with trophic level in
#        hyperdiverse food webs (Libra et al.)
#
# Fig. S9: heatmap of parasitism rate for common caterpillar species
# (>= 50 individuals across the study area), across all localities.
# Species ordered by decreasing total abundance (most abundant on top);
# parasitism rate binned in 10% intervals (1-9, 10-19, ..., 90-99, 100).
#
# NOTE: Ohu2 is KEPT here (all localities shown).
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(readxl)
library(here)
library(grid)       # unit()
library(svglite)    # clean SVG export

# ------------------------------------------------------------
# 1. Load data (robust load; Ohu2 kept)
# ------------------------------------------------------------
MASTER <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(
    CAT_sp = CAT_scientific_name,
    PAR_sp = PAR_species_code
  ) %>%
  mutate(PAR_sp = na_if(as.character(PAR_sp), ""))

all_localities <- sort(unique(MASTER$locality))   # includes Ohu2

# ------------------------------------------------------------
# 2. Common caterpillars (>= 50 individuals) + parasitism flag
# ------------------------------------------------------------
common_cats <- MASTER %>%
  filter(!is.na(CAT_sp)) %>%
  count(CAT_sp, name = "total_abundance") %>%
  filter(total_abundance >= 50)

MASTER_common_cats <- MASTER %>%
  filter(CAT_sp %in% common_cats$CAT_sp) %>%
  mutate(parasitized = ifelse(!is.na(PAR_sp) & PAR_sp != "", 1L, 0L))

# ------------------------------------------------------------
# 3. Parasitism rate per species x locality
# ------------------------------------------------------------
species_locality_raw <- MASTER_common_cats %>%
  group_by(CAT_sp, locality) %>%
  summarise(
    n_total         = n(),
    n_parasitized   = sum(parasitized, na.rm = TRUE),
    parasitism_rate = ifelse(n_total > 0, 100 * n_parasitized / n_total, NA_real_),
    .groups = "drop"
  )

# Complete grid: explicit absences (species not at a locality) -> NA rate
species_locality_complete <- species_locality_raw %>%
  group_by(CAT_sp) %>%
  complete(locality = all_localities) %>%
  ungroup() %>%
  mutate(
    n_total         = replace_na(n_total, 0L),
    n_parasitized   = replace_na(n_parasitized, 0L),
    parasitism_rate = ifelse(n_total == 0, NA_real_, parasitism_rate)
  )

# ------------------------------------------------------------
# 4. Heatmap data + ordering (most abundant on top)
# ------------------------------------------------------------
heatmap_df <- species_locality_complete %>%
  left_join(common_cats, by = "CAT_sp") %>%
  mutate(
    CAT_sp   = fct_reorder(CAT_sp, total_abundance),
    locality = factor(locality, levels = all_localities)
  )

# ------------------------------------------------------------
# 5. Binning (10% intervals: 1-9, 10-19, ..., 90-99, 100) + 0%
# ------------------------------------------------------------
breaks_vec <- c(0, 9, seq(19, 99, by = 10), 100)
labels_vec <- c(
  "1-9%",
  paste0(seq(10, 90, by = 10), "-", seq(19, 99, by = 10), "%"),
  "100%"
)
possible_bins <- c("0%", labels_vec)

heatmap_df <- heatmap_df %>%
  mutate(
    parasitism_bin = case_when(
      is.na(parasitism_rate) ~ NA_character_,
      parasitism_rate == 0   ~ "0%",
      TRUE ~ as.character(cut(parasitism_rate, breaks = breaks_vec,
                              include.lowest = TRUE, right = TRUE, labels = labels_vec))
    ),
    parasitism_bin = factor(parasitism_bin, levels = possible_bins),
    label_text = case_when(
      n_total == 0         ~ "NA",
      parasitism_rate == 0 ~ "0%",
      TRUE                 ~ paste0(round(parasitism_rate), "%")
    )
  )

# ------------------------------------------------------------
# 6. Colourblind-safe palette  >>> EXACT COLOURS - DO NOT CHANGE <<<
# ------------------------------------------------------------
palette_cb_safe <- c(
  "0%"     = "#FFFFFF",
  "1-9%"   = "#E0F3F8",
  "10-19%" = "#ABD9E9",
  "20-29%" = "#74ADD1",
  "30-39%" = "#4575B4",
  "40-49%" = "#313695",
  "50-59%" = "#FEE090",
  "60-69%" = "#FDAE61",
  "70-79%" = "#F46D43",
  "80-89%" = "#D73027",
  "90-99%" = "#A50026",
  "100%"   = "#67001F"
)

# ------------------------------------------------------------
# 7. Figure S9
# ------------------------------------------------------------
# Invisible dummy layer: one tile per bin so EVERY colour appears in the legend,
# even bins with no data. alpha = 0 keeps them invisible in the plot.
legend_dummy <- data.frame(
  locality       = factor(all_localities[1], levels = all_localities),
  CAT_sp         = factor(levels(heatmap_df$CAT_sp)[1], levels = levels(heatmap_df$CAT_sp)),
  parasitism_bin = factor(possible_bins, levels = possible_bins)
)

p_heatmap <- ggplot(heatmap_df,
                    aes(x = locality, y = CAT_sp, fill = parasitism_bin)) +
  geom_tile(data = legend_dummy, alpha = 0) +          # invisible; forces full legend
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = label_text), size = 3.0, color = "black", na.rm = TRUE) +
  scale_fill_manual(values = palette_cb_safe,
                    breaks = names(palette_cb_safe),   # force ALL bins into the legend
                    limits = names(palette_cb_safe),   # full palette as the scale domain
                    drop = FALSE, na.value = "grey90", na.translate = FALSE,
                    labels = function(x) gsub("%", "", x)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, colour = "grey40", linewidth = 0.3))) +
  coord_fixed(ratio = 0.4) +
  labs(x = "Locality", y = "Caterpillar species", fill = "Parasitism rate (%)") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x       = element_text(angle = 70, hjust = 1, vjust = 1, size = 11),
    axis.text.y       = element_text(size = 13, color = "black"),
    axis.title.x      = element_text(size = 13, color = "black"),
    axis.title.y      = element_text(size = 13, color = "black"),
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 10),
    panel.grid        = element_blank(),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width  = unit(1.3, "cm"),
    plot.background   = element_rect(fill = "white", color = NA)
  )

# ------------------------------------------------------------
# 8. Export
# ------------------------------------------------------------
dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)
ggsave("output/fig/Figure_S9_parasitism_heatmap.svg", p_heatmap,
       width = 220, height = 280, units = "mm", device = svglite)
ggsave("output/fig/Figure_S9_parasitism_heatmap.pdf", p_heatmap,
       width = 220, height = 280, units = "mm")
ggsave("output/fig/Figure_S9_parasitism_heatmap.png", p_heatmap,
       width = 220, height = 280, units = "mm", dpi = 300)