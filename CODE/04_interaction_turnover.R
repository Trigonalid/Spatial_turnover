# ============================================================
# 04_interaction_turnover.R   >>> ANALYSIS ONLY <<<
# Study: Approaching the limit of diversity: spatial turnover
#        amplifies with trophic level (Libra et al.)
#
# Analysis 4: interaction (network) turnover, decomposed into
# species turnover (ST) and rewiring (OS), whole-network total (WN),
# via {bipartite} betalinkr (partitioning = "commondenom").
#
# Networks: parasitoid-caterpillar and caterpillar-plant.
# Datasets: all data / common plants only / common caterpillars only.
# Robustness: caterpillar-plant subsampled to parasitoid-caterpillar size.
#
# Outputs: Table S4 (summary) + Table S4b (Wilcoxon) as Word + csv,
#          and output/rds/ files for the plotting script.
# Figures are built in 04_interaction_turnover_plots.R
#   (Fig. 1D, Fig. 2, Fig. S7).
# ------------------------------------------------------------
# METHODS NOTES (consistent with 01/02):
#   - Mantel / distance-decay NOT reported (removed). betalinkr already
#     returns each site-pair once, so no matrix reconstruction is needed.
#   - Subsampling: 999 iterations, set.seed(1234), with replacement.
#   - Decomposition uses the common-denominator variant of betalinkr
#     (ST + OS = WN), building on the Poisot et al. turnover/rewiring
#     framework -> check the Methods wording matches this.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)       # data wrangling
library(tidyr)       # pivot / nest
library(purrr)       # map functions
library(tibble)      # tibble, column_to_rownames
library(readxl)      # read_excel
library(here)        # here() for file paths
library(rlang)       # set_names(), chuck()
library(bipartite)   # webs2array(), betalinkr_multi()
library(flextable)   # flextable() for Word export
library(officer)     # read_docx()

# ------------------------------------------------------------
# 1. Load data (same robust load as 01/02)
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

# ------------------------------------------------------------
# 2. Rarity categories (consistent with 02)
# ------------------------------------------------------------
cat_counts <- MASTER %>%
  filter(!is.na(CAT_sp)) %>%
  group_by(CAT_sp) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(common_caterpillar = n >= 50)               # >= 50 individuals

plant_counts <- MASTER %>%
  filter(!is.na(PLANT_sp)) %>%
  group_by(PLANT_sp) %>%
  summarise(n_localities = n_distinct(locality), .groups = "drop") %>%
  mutate(common_plant = n_localities >= 5)           # present at >= 5 localities

MASTER <- MASTER %>%
  left_join(cat_counts   %>% select(CAT_sp,   common_caterpillar), by = "CAT_sp") %>%
  left_join(plant_counts %>% select(PLANT_sp, common_plant),       by = "PLANT_sp")

# ------------------------------------------------------------
# 3. Filtered datasets
# ------------------------------------------------------------
data_full         <- MASTER %>% select(locality, PAR_sp, CAT_sp, PLANT_sp)
data_common_cat   <- MASTER %>% filter(common_caterpillar == TRUE) %>%
  select(locality, PAR_sp, CAT_sp, PLANT_sp)
data_common_plant <- MASTER %>% filter(common_plant == TRUE) %>%
  select(locality, PAR_sp, CAT_sp, PLANT_sp)

# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------
# Build one contingency table (species_A x species_B) per locality
get_contingency_table <- function(data_source, guild_a, guild_b) {
  data_source %>%
    dplyr::group_by(locality) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      contingency_table = purrr::map(
        .x = data,
        .f = ~ table(purrr::chuck(.x, guild_a), purrr::chuck(.x, guild_b)) %>% as.matrix()
      )
    )
}

# Network dissimilarities across all locality pairs (each pair once)
get_network_dissimilarities <- function(data_source) {
  data_source %>%
    purrr::chuck("contingency_table") %>%
    rlang::set_names(data_source$locality) %>%
    bipartite::webs2array() %>%
    bipartite::betalinkr_multi(partitioning = "commondenom", partition.st = TRUE)
}

metrics_keep <- c("WN", "ST", "OS")

# betalinkr output -> long (WN/ST/OS), each site-pair once
to_long <- function(nd, dataset_label, type_label) {
  nd %>%
    as_tibble() %>%
    dplyr::select(i, j, dplyr::all_of(metrics_keep)) %>%
    tidyr::pivot_longer(dplyr::all_of(metrics_keep),
                        names_to = "metric", values_to = "value") %>%
    dplyr::mutate(dataset = dataset_label, type = type_label)
}

# Mean +/- sd summary for one network x dataset
to_summary <- function(nd, dataset_label, type_label) {
  to_long(nd, dataset_label, type_label) %>%
    dplyr::group_by(dataset, type, metric) %>%
    dplyr::summarise(mean = round(mean(value), 3),
                     sd   = round(sd(value),   3), .groups = "drop")
}

# ------------------------------------------------------------
# 5. Run network dissimilarity analyses (2 networks x 3 datasets)
# ------------------------------------------------------------
nd_cat_par      <- get_network_dissimilarities(get_contingency_table(data_full,         "CAT_sp", "PAR_sp"))
nd_cat_plant    <- get_network_dissimilarities(get_contingency_table(data_full,         "CAT_sp", "PLANT_sp"))
nd_cc_cat_par   <- get_network_dissimilarities(get_contingency_table(data_common_cat,   "CAT_sp", "PAR_sp"))
nd_cc_cat_plant <- get_network_dissimilarities(get_contingency_table(data_common_cat,   "CAT_sp", "PLANT_sp"))
nd_cp_cat_par   <- get_network_dissimilarities(get_contingency_table(data_common_plant, "CAT_sp", "PAR_sp"))
nd_cp_cat_plant <- get_network_dissimilarities(get_contingency_table(data_common_plant, "CAT_sp", "PLANT_sp"))

# ------------------------------------------------------------
# 6. Subsampling — caterpillar-plant to parasitoid-caterpillar size
# At each locality, CAT-PLANT interactions are resampled WITH replacement
# to match the number of PAR-CAT interactions. 999 iterations; per-pair
# mean reported. set.seed(1234). Same betalinkr settings as above.
# ------------------------------------------------------------
set.seed(1234)
n_rand         <- 999
localities_all <- sort(unique(data_full$locality))

# Target size = number of PAR-CAT interactions per locality
n_par_cat_per_locality <- data_full %>%
  filter(!is.na(PAR_sp), !is.na(CAT_sp)) %>%
  count(locality, name = "N")

# CAT-PLANT interactions to resample from
data_cat_plant_full <- data_full %>% filter(!is.na(CAT_sp), !is.na(PLANT_sp))

# Build one locality matrix with fixed row/col species set
create_matrix_local <- function(data_loc, vec_rows, vec_cols, col_row, col_col) {
  mat_local <- data_loc %>%
    dplyr::count(.data[[col_row]], .data[[col_col]]) %>%
    tidyr::pivot_wider(names_from = col_col, values_from = n, values_fill = 0) %>%
    tibble::column_to_rownames(col_row) %>%
    as.matrix()
  
  mat_full <- matrix(0, nrow = length(vec_rows), ncol = length(vec_cols),
                     dimnames = list(vec_rows, vec_cols))
  common_rows <- intersect(rownames(mat_local), vec_rows)
  common_cols <- intersect(colnames(mat_local), vec_cols)
  mat_full[common_rows, common_cols] <- mat_local[common_rows, common_cols]
  mat_full
}

res_sub <- vector("list", n_rand)
for (i in seq_len(n_rand)) {
  
  resampled <- purrr::map_df(localities_all, function(loc) {
    N_target <- n_par_cat_per_locality %>% filter(locality == loc) %>% pull(N)
    if (length(N_target) == 0 || is.na(N_target) || N_target == 0) return(tibble())
    data_cat_plant_full %>% filter(locality == loc) %>% sample_n(size = N_target, replace = TRUE)
  })
  
  vec_plants       <- sort(unique(resampled$PLANT_sp))
  vec_caterpillars <- sort(unique(resampled$CAT_sp))
  loc_names        <- sort(unique(resampled$locality))
  
  loc_matrices <- rlang::set_names(loc_names) %>%
    purrr::map(~ resampled %>% filter(locality == .x) %>%
                 create_matrix_local(vec_plants, vec_caterpillars, "PLANT_sp", "CAT_sp"))
  
  # One betalinkr_multi call per iteration (same settings as main analysis)
  webs <- bipartite::webs2array(loc_matrices)
  res_sub[[i]] <- bipartite::betalinkr_multi(webs, partitioning = "commondenom",
                                             partition.st = TRUE) %>%
    as_tibble() %>% mutate(iteration = i)
}

# Per-pair mean across iterations -> long (WN/ST/OS)
sub_mean <- bind_rows(res_sub) %>%
  group_by(i, j) %>%
  summarise(across(all_of(metrics_keep), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(all_of(metrics_keep), names_to = "metric", values_to = "value") %>%
  mutate(dataset = "Caterpillar-Plant subsampled", type = "subsampled")

# ------------------------------------------------------------
# 7. Assemble long-format data for plotting + save
# ------------------------------------------------------------
# Fig. 1D / Fig. 2: Par-Cat vs Cat-Plant vs Cat-Plant subsampled (full data)
plot_data_main <- bind_rows(
  to_long(nd_cat_par,   "Parasitoid-Caterpillar", "full"),
  to_long(nd_cat_plant, "Caterpillar-Plant",      "full"),
  sub_mean
)
saveRDS(plot_data_main, "output/rds/interaction_turnover_main.rds")

# Fig. S7: WN by dataset type, both networks
wn_by_type <- bind_rows(
  to_long(nd_cat_par,      "Parasitoid-Caterpillar", "full"),
  to_long(nd_cc_cat_par,   "Parasitoid-Caterpillar", "common caterpillars only"),
  to_long(nd_cp_cat_par,   "Parasitoid-Caterpillar", "common plants only"),
  to_long(nd_cat_plant,    "Caterpillar-Plant",      "full"),
  to_long(nd_cc_cat_plant, "Caterpillar-Plant",      "common caterpillars only"),
  to_long(nd_cp_cat_plant, "Caterpillar-Plant",      "common plants only")
) %>%
  filter(metric == "WN") %>%
  rename(network = dataset)
saveRDS(wn_by_type, "output/rds/interaction_turnover_wn_by_type.rds")

# ------------------------------------------------------------
# 8. Summary table (Table S4) — Word + csv
# ------------------------------------------------------------
summary_sub <- sub_mean %>%
  group_by(dataset, type, metric) %>%
  summarise(mean = round(mean(value), 3), sd = round(sd(value), 3), .groups = "drop") %>%
  mutate(type = "Subsampled")

results_merged <- bind_rows(
  to_summary(nd_cat_par,      "Parasitoid-Caterpillar", "All data"),
  to_summary(nd_cat_plant,    "Caterpillar-Plant",      "All data"),
  to_summary(nd_cc_cat_par,   "Parasitoid-Caterpillar", "Common caterpillars only"),
  to_summary(nd_cc_cat_plant, "Caterpillar-Plant",      "Common caterpillars only"),
  to_summary(nd_cp_cat_par,   "Parasitoid-Caterpillar", "Common plants only"),
  to_summary(nd_cp_cat_plant, "Caterpillar-Plant",      "Common plants only"),
  summary_sub
) %>%
  select(dataset, type, metric, mean, sd)

write.csv(results_merged, "output/rds/Table_S4_network_summary.csv", row.names = FALSE)
doc_network <- officer::read_docx() %>%
  officer::body_add_par("Table S4: Summary of network dissimilarity indices",
                        style = "heading 1") %>%
  flextable::body_add_flextable(flextable(results_merged))
print(doc_network, target = "output/Table_S4_network_summary.docx")

# ------------------------------------------------------------
# 9. Wilcoxon tests (Table S4b) — Word + csv
# PAIRED Wilcoxon signed-rank tests: values are matched by site-pair (i, j),
# since both networks are compared over the same pairs of localities.
# Statistic is V (signed-rank).
# ------------------------------------------------------------
wilcox_pairs <- list(
  c("Parasitoid-Caterpillar", "Caterpillar-Plant"),
  c("Parasitoid-Caterpillar", "Caterpillar-Plant subsampled"),
  c("Caterpillar-Plant",      "Caterpillar-Plant subsampled")
)

fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

wilcox_fw <- purrr::map_df(metrics_keep, function(m) {
  df <- plot_data_main %>% filter(metric == m)
  purrr::map_df(wilcox_pairs, function(pair) {
    a <- df %>% filter(dataset == pair[1]) %>% dplyr::select(i, j, value_a = value)
    b <- df %>% filter(dataset == pair[2]) %>% dplyr::select(i, j, value_b = value)
    paired_df <- dplyr::inner_join(a, b, by = c("i", "j"))   # align by site-pair
    if (nrow(paired_df) < 2) {
      return(tibble(metric = m, group_1 = pair[1], group_2 = pair[2], n_pairs = nrow(paired_df),
                    V_statistic = NA_real_, p_value = NA_character_, significance = NA_character_))
    }
    wt <- wilcox.test(paired_df$value_a, paired_df$value_b, paired = TRUE, exact = FALSE)
    tibble(
      metric       = m,
      group_1      = pair[1],
      group_2      = pair[2],
      n_pairs      = nrow(paired_df),
      V_statistic  = round(wt$statistic, 1),
      p_value      = fmt_p(wt$p.value),
      significance = case_when(
        wt$p.value < 0.001 ~ "***",
        wt$p.value < 0.01  ~ "**",
        wt$p.value < 0.05  ~ "*",
        TRUE               ~ "ns"
      )
    )
  })
})

print(wilcox_fw)
write.csv(wilcox_fw, "output/rds/Table_S4b_network_wilcoxon.csv", row.names = FALSE)
doc_wilcox_fw <- officer::read_docx() %>%
  officer::body_add_par("Table S4b: Paired Wilcoxon signed-rank tests - network type comparisons",
                        style = "heading 1") %>%
  flextable::body_add_flextable(flextable(wilcox_fw))
print(doc_wilcox_fw, target = "output/Table_S4b_network_wilcoxon.docx")
# ------------------------------------------------------------
# 10. Wilcoxon tests for Fig. S6 (rare-removal) — Word + csv
# WN interaction dissimilarity, PAIRED signed-rank (matched by site-pair i,j):
#   (A) BETWEEN networks (Cat-Plant vs Par-Cat) at each dataset level
#   (B) WITHIN each network, across dataset reductions
# ------------------------------------------------------------
type_lab <- c("full" = "All",
              "common plants only" = "Common plants",
              "common caterpillars only" = "Common caterpillars")

# general paired WN test between two (network, type) groups, aligned by (i, j)
wn_paired <- function(net_a, type_a, net_b, type_b, comparison_label) {
  a <- wn_by_type %>% dplyr::filter(network == net_a, type == type_a) %>%
    dplyr::select(i, j, va = value)
  b <- wn_by_type %>% dplyr::filter(network == net_b, type == type_b) %>%
    dplyr::select(i, j, vb = value)
  m <- dplyr::inner_join(a, b, by = c("i", "j"))
  g1 <- paste0(net_a, " (", type_lab[type_a], ")")
  g2 <- paste0(net_b, " (", type_lab[type_b], ")")
  if (nrow(m) < 2) {
    return(tibble(comparison = comparison_label, group_1 = g1, group_2 = g2,
                  n_pairs = nrow(m), V_statistic = NA_real_,
                  p_value = NA_character_, significance = NA_character_))
  }
  wt <- wilcox.test(m$va, m$vb, paired = TRUE, exact = FALSE)
  tibble(
    comparison   = comparison_label,
    group_1      = g1,
    group_2      = g2,
    n_pairs      = nrow(m),
    V_statistic  = round(wt$statistic, 1),
    p_value      = fmt_p(wt$p.value),
    significance = case_when(wt$p.value < 0.001 ~ "***",
                             wt$p.value < 0.01  ~ "**",
                             wt$p.value < 0.05  ~ "*",
                             TRUE               ~ "ns")
  )
}

types      <- c("full", "common plants only", "common caterpillars only")
type_pairs <- list(c("full", "common plants only"),
                   c("full", "common caterpillars only"),
                   c("common plants only", "common caterpillars only"))

# (A) between networks, one per dataset level
between_net <- purrr::map_df(types, function(t)
  wn_paired("Parasitoid-Caterpillar", t, "Caterpillar-Plant", t,
            "Between networks"))

# (B) within each network, across dataset reductions
within_net <- purrr::map_df(c("Parasitoid-Caterpillar", "Caterpillar-Plant"), function(net)
  purrr::map_df(type_pairs, function(pr)
    wn_paired(net, pr[1], net, pr[2], paste0("Within ", net))))

wilcox_s6 <- dplyr::bind_rows(between_net, within_net)
print(wilcox_s6)

write.csv(wilcox_s6, "output/rds/Table_S10_rareremoval_wilcoxon.csv", row.names = FALSE)
doc_wilcox_s6 <- officer::read_docx() %>%
  officer::body_add_par(
    paste("Table S10: Paired Wilcoxon signed-rank tests - effect of removing rare plant and",
          "caterpillar species on interaction dissimilarity (Fig. S6). Whole-network (WN)",
          "dissimilarity compared (A) between the parasitoid-caterpillar and caterpillar-plant",
          "networks at each dataset level, and (B) among dataset reductions (All / Common plants /",
          "Common caterpillars) within each network. All tests matched by pair of localities.",
          "V, signed-rank statistic; n_pairs, number of locality pairs; P, two-sided P value."),
    style = "heading 1") %>%
  flextable::body_add_flextable(autofit(flextable(wilcox_s6)))
print(doc_wilcox_s6, target = "output/Table_S10_rareremoval_wilcoxon.docx")