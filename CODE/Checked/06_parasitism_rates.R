# =====================================================================
# 06_parasitism_rates.R
# Heatmap of parasitism rates of common caterpillar species
# Ordered by total abundance (most abundant on top)
# Colourblind-safe
# Fig. S9
# =====================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(readxl)
library(here)

# ---------------------------------------------------------------------
# 0) Load data (same load as the rest of the pipeline; Ohu2 stays in —
#    kept as a separate locality here, not merged with Ohu1)
# ---------------------------------------------------------------------
MASTER <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(
    CAT_sp = CAT_scientific_name,
    PAR_sp = PAR_species_code
  ) %>%
  mutate(PAR_sp = na_if(as.character(PAR_sp), ""))

all_localities <- sort(unique(MASTER$locality))

# ---------------------------------------------------------------------
# 1) Select common caterpillars (>= 50 individuals region-wide) and
#    define parasitism. (Replaces the old `fifty` column, which existed
#    only in MASTER.csv and is not present in MASTER.xlsx — the threshold
#    is now computed directly from the data.)
# ---------------------------------------------------------------------
cat_counts <- MASTER %>%
  filter(!is.na(CAT_sp)) %>%
  group_by(CAT_sp) %>%
  summarise(n_individuals = n(), .groups = "drop") %>%
  mutate(common_caterpillar = n_individuals >= 50)

MASTER_common_cats <- MASTER %>%
  left_join(cat_counts %>% select(CAT_sp, common_caterpillar), by = "CAT_sp") %>%
  filter(common_caterpillar == TRUE) %>%
  mutate(
    parasitized = ifelse(!is.na(PAR_sp) & PAR_sp != "", 1L, 0L)
  )

# ---------------------------------------------------------------------
# 2) Aggregate parasitism by species × locality
# ---------------------------------------------------------------------
species_locality_raw <- MASTER_common_cats %>%
  group_by(CAT_sp, locality) %>%
  summarise(
    n_total = n(),
    n_parasitized = sum(parasitized, na.rm = TRUE),
    parasitism_rate = ifelse(
      n_total > 0,
      100 * n_parasitized / n_total,
      NA_real_
    ),
    .groups = "drop"
  )

# ---------------------------------------------------------------------
# 3) Complete species × locality matrix (explicit absences)
# ---------------------------------------------------------------------
species_locality_complete <- species_locality_raw %>%
  group_by(CAT_sp) %>%
  complete(locality = all_localities) %>%
  ungroup() %>%
  mutate(
    n_total = replace_na(n_total, 0L),
    n_parasitized = replace_na(n_parasitized, 0L),
    parasitism_rate = ifelse(n_total == 0, NA_real_, parasitism_rate)
  )

# ---------------------------------------------------------------------
# 4) Total abundance per species (ordering variable)
# ---------------------------------------------------------------------
species_abundance <- MASTER_common_cats %>%
  group_by(CAT_sp) %>%
  summarise(
    total_abundance = n(),
    .groups = "drop"
  )

# ---------------------------------------------------------------------
# 5) Heatmap dataframe + ORDERING (most abundant on top)
# ---------------------------------------------------------------------
heatmap_df <- species_locality_complete %>%
  left_join(species_abundance, by = "CAT_sp") %>%
  mutate(
    CAT_sp   = fct_reorder(CAT_sp, total_abundance),
    locality = factor(locality, levels = all_localities)
  )

# ---------------------------------------------------------------------
# 6) Binning of parasitism rates (10% intervals)
# ---------------------------------------------------------------------
breaks_vec <- c(0, 9, seq(19, 99, by = 10), 100)

labels_vec <- c(
  "1-9%",
  paste0(seq(10, 90, by = 10), "-", seq(19, 99, by = 10), "%"),
  "100%"
)

heatmap_df <- heatmap_df %>%
  mutate(
    parasitism_bin = case_when(
      is.na(parasitism_rate) ~ NA_character_,
      parasitism_rate == 0 ~ "0%",
      TRUE ~ as.character(
        cut(
          parasitism_rate,
          breaks = breaks_vec,
          include.lowest = TRUE,
          right = TRUE,
          labels = labels_vec
        )
      )
    )
  )

# ---------------------------------------------------------------------
# 7) Colourblind-safe palette
# ---------------------------------------------------------------------
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

possible_bins <- c("0%", labels_vec)

heatmap_df$parasitism_bin <- factor(
  heatmap_df$parasitism_bin,
  levels = possible_bins
)

# ---------------------------------------------------------------------
# 8) Text inside cells
# ---------------------------------------------------------------------
heatmap_df <- heatmap_df %>%
  mutate(
    label_text = case_when(
      n_total == 0 ~ "NA",
      parasitism_rate == 0 ~ "0%",
      TRUE ~ paste0(round(parasitism_rate), "%")
    )
  )

# ---------------------------------------------------------------------
# 9) Plot heatmap
# ---------------------------------------------------------------------
p_heatmap <- ggplot(
  heatmap_df,
  aes(x = locality, y = CAT_sp, fill = parasitism_bin)
) +
  geom_tile(color = "white", size = 0.25) +
  geom_text(
    aes(label = label_text),
    size = 3.0,
    color = "black",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = palette_cb_safe,
    na.value = "grey90",
    drop = FALSE,
    labels = function(x) gsub("%", "", x)
  ) +
  coord_fixed(ratio = 0.4) +
  labs(
    x = "Locality",
    y = "Caterpillar species",
    fill = "Parasitism rate (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 70, hjust = 1, vjust = 1, size = 11
    ),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 13, color = "black"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank(),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width  = unit(1.3, "cm")
  )

print(p_heatmap)

# ---------------------------------------------------------------------
# 10) Export — Fig. S9
# ---------------------------------------------------------------------
dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)
ggsave("output/fig/Figure_S9.svg", p_heatmap, width = 220, height = 280, units = "mm")
ggsave("output/fig/Figure_S9.pdf", p_heatmap, width = 220, height = 280, units = "mm")
ggsave("output/fig/Figure_S9.png", p_heatmap, width = 220, height = 280, units = "mm", dpi = 300)