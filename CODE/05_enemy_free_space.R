# ============================================================
# 05_enemy_free_space.R   >>> ANALYSIS ONLY <<<
# Study: Approaching the limit of diversity: spatial turnover
#        amplifies with trophic level (Libra et al.)
#
# Analysis 5: enemy-free space.
#   S8A: proportion of a caterpillar's occupied sites with NO parasitoid.
#   S8B: proportion of a caterpillar's occupied sites free from a SPECIFIC parasitoid.
# Only caterpillars parasitised at least once across the study area are included.
#
# Outputs: EFS tables (csv) + summary-stats table (Word + csv) +
#          output/rds/ for the plotting script.
# Figures are built in 05_enemy_free_space_plots.R (Fig. S8).
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)       # data wrangling
library(tidyr)       # pivot functions
library(readxl)      # read_excel
library(here)        # here() for file paths
library(flextable)   # flextable() for Word export
library(officer)     # read_docx()

# ------------------------------------------------------------
# 1. Load data (robust load; Ohu2 excluded, temporal replicate)
# ------------------------------------------------------------
MASTER <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(
    PLANT_sp = PLANT_species_code,
    CAT_sp   = CAT_scientific_name,
    PAR_sp   = PAR_species_code
  ) %>%
  mutate(PAR_sp = na_if(as.character(PAR_sp), "")) %>%
  filter(locality != "Ohu2")

dir.create("output",     showWarnings = FALSE, recursive = TRUE)
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)

# Keep only caterpillar species parasitised at least once anywhere
vec_cat_with_par <- MASTER %>%
  group_by(CAT_sp) %>%
  summarise(has_parasitoid = any(!is.na(PAR_sp)), .groups = "drop") %>%
  filter(has_parasitoid) %>%
  pull(CAT_sp)

data_work <- MASTER %>% filter(CAT_sp %in% vec_cat_with_par)

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------
# S8A: per caterpillar, proportion of its sites without any parasitoid
get_summary_per_cat <- function(data_source) {
  data_source %>%
    group_by(CAT_sp) %>%
    summarise(
      .groups             = "drop",
      localities_total    = n_distinct(locality),
      localities_with_par = n_distinct(locality[!is.na(PAR_sp)]),
      localities_no_par   = localities_total - localities_with_par,
      localities_prop     = localities_no_par / localities_total
    ) %>%
    select(CAT_sp, localities_total, localities_prop)
}

# S8B: per cat-par pair, proportion of the caterpillar's sites free from that parasitoid
get_summary_per_pair <- function(data_source) {
  cat_localities <- data_source %>%
    group_by(CAT_sp) %>%
    summarise(localities_total_cat = n_distinct(locality), .groups = "drop")

  pair_localities <- data_source %>%
    filter(!is.na(PAR_sp)) %>%
    group_by(CAT_sp, PAR_sp) %>%
    summarise(localities_with_pair = n_distinct(locality), .groups = "drop")

  pair_localities %>%
    left_join(cat_localities, by = "CAT_sp") %>%
    mutate(localities_prop = (localities_total_cat - localities_with_pair) / localities_total_cat) %>%
    select(CAT_sp, PAR_sp, localities_prop)
}

# ------------------------------------------------------------
# 3. Run analyses
# ------------------------------------------------------------
summary_per_cat  <- get_summary_per_cat(data_work)
summary_per_pair <- get_summary_per_pair(data_work)

saveRDS(summary_per_cat,  "output/rds/efs_per_cat.rds")
saveRDS(summary_per_pair, "output/rds/efs_per_pair.rds")

# ------------------------------------------------------------
# 4. Summary statistics
# ------------------------------------------------------------
efs_summary_stats <- bind_rows(
  summary_per_cat %>%
    summarise(panel = "S8A: without any parasitoid",
              mean = round(mean(localities_prop, na.rm = TRUE), 3),
              median = round(median(localities_prop, na.rm = TRUE), 3),
              n = n()),
  summary_per_pair %>%
    summarise(panel = "S8B: without a particular parasitoid",
              mean = round(mean(localities_prop, na.rm = TRUE), 3),
              median = round(median(localities_prop, na.rm = TRUE), 3),
              n = n())
)
print(efs_summary_stats)

# ------------------------------------------------------------
# 5. Export tables (csv + Word)
# ------------------------------------------------------------
# Raw cat-par interaction table with locality counts
tab_efs_raw <- MASTER %>%
  filter(!is.na(PAR_sp)) %>%
  group_by(CAT_sp, PAR_sp) %>%
  summarise(localities_total = n_distinct(locality), .groups = "drop")

write.csv(tab_efs_raw,       "output/rds/Table_S_EFS_raw.csv",      row.names = FALSE)
write.csv(summary_per_cat,   "output/rds/Table_S_EFS_per_cat.csv",  row.names = FALSE)
write.csv(summary_per_pair,  "output/rds/Table_S_EFS_per_pair.csv", row.names = FALSE)
write.csv(efs_summary_stats, "output/rds/Table_S_EFS_summary.csv",  row.names = FALSE)

doc_efs <- officer::read_docx() %>%
  officer::body_add_par("Table S8: Enemy-free space summary statistics",
                        style = "heading 1") %>%
  flextable::body_add_flextable(flextable(efs_summary_stats))
print(doc_efs, target = "output/Table_S8_enemy_free_space_summary.docx")
