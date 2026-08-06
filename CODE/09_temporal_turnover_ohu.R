# ============================================================
# 09_temporal_turnover_ohu.R   >>> ANALYSIS ONLY <<<
# Study: Approaching the limit of diversity: spatial turnover
#        amplifies with trophic level (Libra et al.)
#
# Analysis 9: temporal turnover at Ohu (two time points).
#   Ohu1 = first time point (larger dataset)
#   Ohu2 = second time point (smaller dataset)
#
#   Panel A: Bray-Curtis community dissimilarity (caterpillars, parasitoids),
#            observed + subsampled (larger time point -> smaller, 999x).
#   Panel B: WN interaction dissimilarity (par-cat; cat-plant observed;
#            cat-plant subsampled to par-cat size per time point, 999x).
#
# Outputs: Table S6 (Word + csv) + output/rds/ for the plotting script.
# Figures are built in 09_temporal_turnover_ohu_plots.R (Fig. S6).
# ------------------------------------------------------------
# METHODS NOTES:
#   - MASTER is NOT filtered here: this analysis needs both Ohu1 and Ohu2.
#   - Subsampling: 999 iterations, set.seed(1234), with replacement.
#   - Two subsampling schemes: Bray = larger time point down to smaller;
#     WN = caterpillar-plant down to parasitoid-caterpillar size.

# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)       # data wrangling
library(tidyr)       # pivot functions
library(purrr)       # map functions
library(tibble)      # tibble, column_to_rownames
library(readxl)      # read_excel
library(here)        # here() for file paths
library(vegan)       # vegdist()
library(bipartite)   # betalinkr_multi()
library(flextable)   # flextable() for Word export
library(officer)     # read_docx()

set.seed(1234)
n_rand <- 999

# ------------------------------------------------------------
# 1. Load data (robust load; Ohu2 kept for the temporal analysis)
# ------------------------------------------------------------
MASTER <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(
    PLANT_sp = PLANT_species_code,
    CAT_sp   = CAT_scientific_name,
    PAR_sp   = PAR_species_code
  ) %>%
  mutate(PAR_sp = na_if(as.character(PAR_sp), ""))

dir.create("output",     showWarnings = FALSE, recursive = TRUE)
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)

# Ohu temporal replicates only
data_ohu <- MASTER %>% filter(locality %in% c("Ohu1", "Ohu2"))

# Caterpillar-Plant interactions
data_cat_plant <- data_ohu %>%
  filter(!is.na(CAT_sp), !is.na(PLANT_sp)) %>%
  select(locality, caterpillar = CAT_sp, plant = PLANT_sp)

# Parasitoid-Caterpillar interactions
data_par_cat <- data_ohu %>%
  filter(!is.na(CAT_sp), !is.na(PAR_sp)) %>%
  select(locality, caterpillar = CAT_sp, parasitoid = PAR_sp)

# ------------------------------------------------------------
# 2. Bray-Curtis helpers (Panel A)
# ------------------------------------------------------------
get_species_counts <- function(data) {
  data %>%
    count(locality, sp) %>%
    pivot_wider(names_from = locality, values_from = n, values_fill = 0) %>%
    as.data.frame()
}
get_species_matrix <- function(df) {
  df %>% column_to_rownames("sp") %>% t() %>% as.matrix()
}
get_bray <- function(mat) {
  as.matrix(vegdist(mat, method = "bray"))[1, 2]
}
# Downsample the larger time point to the smaller (with replacement)
resample_bray <- function(data) {
  counts    <- data %>% count(locality)
  small_n   <- min(counts$n)
  big_loc   <- counts %>% filter(n == max(n)) %>% pull(locality)
  small_loc <- counts %>% filter(n == min(n)) %>% pull(locality)

  bind_rows(
    data %>% filter(locality == big_loc)   %>% sample_n(small_n, replace = TRUE),
    data %>% filter(locality == small_loc)
  ) %>%
    get_species_counts() %>% get_species_matrix() %>% get_bray()
}
make_bray_data <- function(data, guild_name) {
  observed   <- data %>% get_species_counts() %>% get_species_matrix() %>% get_bray()
  subsampled <- map_dbl(seq_len(n_rand), ~ resample_bray(data))
  bind_rows(
    tibble(type = "Observed",   value = observed),
    tibble(type = "Subsampled", value = subsampled)
  ) %>%
    mutate(guild = guild_name)
}

# ------------------------------------------------------------
# 3. Bray-Curtis data (caterpillars, parasitoids)
# ------------------------------------------------------------
data_cat_bc <- data_ohu %>% filter(!is.na(CAT_sp)) %>% select(locality, sp = CAT_sp)
data_par_bc <- data_ohu %>% filter(!is.na(PAR_sp)) %>% select(locality, sp = PAR_sp)

bray_data_all <- bind_rows(
  make_bray_data(data_cat_bc, "Caterpillars"),
  make_bray_data(data_par_bc, "Parasitoids")
) %>%
  mutate(
    x_group = case_when(
      guild == "Caterpillars" & type == "Observed"   ~ "Caterpillars\nObserved",
      guild == "Caterpillars" & type == "Subsampled" ~ "Caterpillars\nSubsampled",
      guild == "Parasitoids"  & type == "Observed"   ~ "Parasitoids\nObserved",
      guild == "Parasitoids"  & type == "Subsampled" ~ "Parasitoids\nSubsampled"
    ),
    x_group = factor(x_group, levels = c(
      "Caterpillars\nObserved", "Caterpillars\nSubsampled",
      "Parasitoids\nObserved",  "Parasitoids\nSubsampled"
    ))
  )

saveRDS(bray_data_all, "output/rds/temporal_bray.rds")

# ------------------------------------------------------------
# 4. WN helpers (Panel B)
# ------------------------------------------------------------
make_matrix_full <- function(df, lower, upper, all_l, all_u) {
  mat_local <- df %>%
    count(.data[[lower]], .data[[upper]]) %>%
    pivot_wider(names_from = .data[[upper]], values_from = n, values_fill = 0) %>%
    column_to_rownames(lower) %>%
    as.matrix()

  mat_full <- matrix(0, nrow = length(all_l), ncol = length(all_u),
                     dimnames = list(all_l, all_u))
  mat_full[rownames(mat_local), colnames(mat_local)] <- mat_local
  mat_full
}
# WN dissimilarity between two matrices (same betalinkr settings as script 04)
compute_wn <- function(m1, m2) {
  webs      <- array(dim = c(dim(m1), 2))
  webs[, , 1] <- m1
  webs[, , 2] <- m2
  bipartite::betalinkr_multi(webs, partitioning = "commondenom", partition.st = TRUE)$WN
}

# ------------------------------------------------------------
# 5. WN calculations
# ------------------------------------------------------------
# --- Parasitoid-Caterpillar (observed) ---
pc_1 <- data_par_cat %>% filter(locality == "Ohu1")
pc_2 <- data_par_cat %>% filter(locality == "Ohu2")
all_cat_pc <- union(pc_1$caterpillar, pc_2$caterpillar)
all_par    <- union(pc_1$parasitoid,  pc_2$parasitoid)
wn_pc_raw <- compute_wn(
  make_matrix_full(pc_1, "caterpillar", "parasitoid", all_cat_pc, all_par),
  make_matrix_full(pc_2, "caterpillar", "parasitoid", all_cat_pc, all_par)
)

# --- Caterpillar-Plant (observed) ---
cpl_1 <- data_cat_plant %>% filter(locality == "Ohu1")
cpl_2 <- data_cat_plant %>% filter(locality == "Ohu2")
all_cat_cp <- union(cpl_1$caterpillar, cpl_2$caterpillar)
all_plants <- union(cpl_1$plant,       cpl_2$plant)
wn_cpl_raw <- compute_wn(
  make_matrix_full(cpl_1, "plant", "caterpillar", all_plants, all_cat_cp),
  make_matrix_full(cpl_2, "plant", "caterpillar", all_plants, all_cat_cp)
)

# --- Caterpillar-Plant subsampled to par-cat size per time point ---
n_pc_ohu1 <- nrow(pc_1)
n_pc_ohu2 <- nrow(pc_2)
wn_cpl_sub <- map_dbl(seq_len(n_rand), function(i) {
  sub_1 <- cpl_1 %>% sample_n(n_pc_ohu1, replace = TRUE)
  sub_2 <- cpl_2 %>% sample_n(n_pc_ohu2, replace = TRUE)
  compute_wn(
    make_matrix_full(sub_1, "plant", "caterpillar", all_plants, all_cat_cp),
    make_matrix_full(sub_2, "plant", "caterpillar", all_plants, all_cat_cp)
  )
})

# ------------------------------------------------------------
# 6. WN plot data
# ------------------------------------------------------------
plot_data_wn <- bind_rows(
  tibble(dataset = "Caterpillar\u2013Plant",                WN = wn_cpl_raw),
  tibble(dataset = "Caterpillar\u2013Plant (subsampled)",   WN = wn_cpl_sub),
  tibble(dataset = "Parasitoid\u2013Caterpillar",           WN = wn_pc_raw)
) %>%
  mutate(dataset = factor(dataset, levels = c(
    "Caterpillar\u2013Plant",
    "Caterpillar\u2013Plant (subsampled)",
    "Parasitoid\u2013Caterpillar"
  )))

saveRDS(plot_data_wn, "output/rds/temporal_wn.rds")

# ------------------------------------------------------------
# 7. Summary tables (Table S6) — Word + csv
# ------------------------------------------------------------
summary_bc <- bray_data_all %>%
  group_by(guild, type) %>%
  summarise(mean = round(mean(value, na.rm = TRUE), 3),
            sd   = round(sd(value,   na.rm = TRUE), 3),
            n    = n(), .groups = "drop") %>%
  rename(community = guild)

summary_wn <- plot_data_wn %>%
  group_by(dataset) %>%
  summarise(mean = round(mean(WN), 3),
            sd   = ifelse(n() > 1, round(sd(WN), 3), NA_real_),
            n    = n(), .groups = "drop")

write.csv(summary_bc, "output/rds/Table_S6_bray_communities.csv", row.names = FALSE)
write.csv(summary_wn, "output/rds/Table_S6_wn_interactions.csv",  row.names = FALSE)

doc_s6 <- officer::read_docx() %>%
  officer::body_add_par("Table S6a: Temporal Bray-Curtis dissimilarity (Ohu)",
                        style = "heading 1") %>%
  flextable::body_add_flextable(flextable(summary_bc)) %>%
  officer::body_add_par("Table S6b: Temporal WN interaction dissimilarity (Ohu)",
                        style = "heading 1") %>%
  flextable::body_add_flextable(flextable(summary_wn))
print(doc_s6, target = "output/Table_S6_temporal_summary.docx")
