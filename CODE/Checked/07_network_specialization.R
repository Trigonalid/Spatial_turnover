# ============================================================
# 07_network_specialization.R   >>> ANALYSIS ONLY <<<
# Study: Spatial turnover amplifies with trophic level in
#        hyperdiverse food webs (Libra et al.)
#
# Analysis 7: network-level specialization H2' (Bluthgen, sample-size
# corrected) per locality, for parasitoid-caterpillar and caterpillar-plant
# food webs. Wilcoxon rank-sum across localities. Fig. S7, Table S6.
#
# Species-inclusion variants (as in Methods):
#   main        : >= 5 reared specimens -> parasitoids AND caterpillars >= 5
#                 (plants unrestricted). This is the MAIN analysis / Fig. S10.
#   no_threshold: all species (no filtering) -> robustness
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readxl)
library(here)
library(bipartite)   # networklevel()
library(flextable)
library(officer)

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

# ------------------------------------------------------------
# 2. Species with >= 5 specimens (regional)
# ------------------------------------------------------------
par_counts <- MASTER %>%
  filter(!is.na(PAR_sp)) %>% group_by(PAR_sp) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(min5_par = n >= 5)
cat_counts_all <- MASTER %>%
  filter(!is.na(CAT_sp)) %>% group_by(CAT_sp) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(min5_cat = n >= 5)
plant_counts <- MASTER %>%
  filter(!is.na(PLANT_sp)) %>% group_by(PLANT_sp) %>%
  summarise(n = n(), .groups = "drop") %>% mutate(min5_plant = n >= 5)

MASTER <- MASTER %>%
  left_join(par_counts     %>% select(PAR_sp,   min5_par),   by = "PAR_sp") %>%
  left_join(cat_counts_all %>% select(CAT_sp,   min5_cat),   by = "CAT_sp") %>%
  left_join(plant_counts   %>% select(PLANT_sp, min5_plant), by = "PLANT_sp")

cat("\n--- Species counts at >= 5 specimens threshold ---\n")
cat(sprintf("  Parasitoids >= 5 specimens:  %d of %d species\n",
            sum(par_counts$min5_par),     nrow(par_counts)))
cat(sprintf("  Caterpillars >= 5 specimens: %d of %d species\n",
            sum(cat_counts_all$min5_cat), nrow(cat_counts_all)))
cat(sprintf("  Plants >= 5 specimens:       %d of %d species\n\n",
            sum(plant_counts$min5_plant), nrow(plant_counts)))

# ------------------------------------------------------------
# 3. Datasets for the three variants
# ------------------------------------------------------------
# main = >= 5 reared specimens (parasitoids & caterpillars); plants unrestricted
data_par_cat_main   <- MASTER %>% filter(!is.na(PAR_sp), !is.na(CAT_sp), min5_par, min5_cat) %>% select(locality, PAR_sp, CAT_sp)
data_cat_plant_main <- MASTER %>% filter(!is.na(CAT_sp), !is.na(PLANT_sp), min5_cat) %>% select(locality, CAT_sp, PLANT_sp)

# no_threshold = all species
data_par_cat_full   <- MASTER %>% filter(!is.na(PAR_sp), !is.na(CAT_sp)) %>% select(locality, PAR_sp, CAT_sp)
data_cat_plant_full <- MASTER %>% filter(!is.na(CAT_sp), !is.na(PLANT_sp)) %>% select(locality, CAT_sp, PLANT_sp)

# ------------------------------------------------------------
# 4. Helpers: build bipartite web + H2' per locality
# ------------------------------------------------------------
build_web <- function(df, row_var, col_var) {
  if (nrow(df) == 0) return(matrix(numeric(0), 0, 0))
  df %>%
    filter(!is.na(.data[[row_var]]), !is.na(.data[[col_var]])) %>%
    count(.data[[row_var]], .data[[col_var]]) %>%
    pivot_wider(names_from = all_of(col_var), values_from = n, values_fill = 0) %>%
    column_to_rownames(row_var) %>%
    as.matrix()
}
safe_H2 <- function(web) {
  if (is.null(web) || nrow(web) < 2 || ncol(web) < 2) return(NA_real_)
  tryCatch(unname(bipartite::networklevel(web, index = "H2")),
           error = function(e) NA_real_)
}

# H2' per locality for one dataset variant
run_variant <- function(df, row_var, col_var, variant_name, network_name) {
  df %>%
    filter(!is.na(.data[[row_var]]), !is.na(.data[[col_var]])) %>%
    group_by(locality) %>% nest() %>%
    mutate(web = purrr::map(data, ~ build_web(.x, row_var, col_var)),
           H2  = purrr::map_dbl(web, safe_H2),
           variant = variant_name, network = network_name) %>%
    ungroup() %>%
    select(locality, variant, network, H2)
}

# ------------------------------------------------------------
# 5. H2' per locality (2 networks x 3 variants)
# ------------------------------------------------------------
H2_local_all <- bind_rows(
  run_variant(data_par_cat_main, "CAT_sp",   "PAR_sp", "\u22655 reared specimens (parasitoids & caterpillars)", "Parasitoid-Caterpillar"),
  run_variant(data_par_cat_full, "CAT_sp",   "PAR_sp", "All species (no threshold)",                          "Parasitoid-Caterpillar"),
  run_variant(data_cat_plant_main, "PLANT_sp", "CAT_sp", "\u22655 reared specimens (parasitoids & caterpillars)", "Caterpillar-Plant"),
  run_variant(data_cat_plant_full, "PLANT_sp", "CAT_sp", "All species (no threshold)",                          "Caterpillar-Plant")
)

saveRDS(H2_local_all, "output/rds/spec_H2_local.rds")

# ------------------------------------------------------------
# 6. Summary: mean +/- SD H2' per variant x network
# ------------------------------------------------------------
variant_levels <- c("\u22655 reared specimens (parasitoids & caterpillars)",
                    "All species (no threshold)")
network_levels <- c("Caterpillar-Plant", "Parasitoid-Caterpillar")

H2_local_summary <- H2_local_all %>%
  group_by(variant, network) %>%
  summarise(n_localities = sum(!is.na(H2)),
            H2_mean = round(mean(H2, na.rm = TRUE), 3),
            H2_sd   = round(sd(H2,   na.rm = TRUE), 3), .groups = "drop") %>%
  mutate(variant = factor(variant, levels = variant_levels),
         network = factor(network, levels = network_levels)) %>%
  arrange(variant, network)

cat("--- H2' per locality (mean +/- SD) ---\n"); print(H2_local_summary); cat("\n")

# ------------------------------------------------------------
# 7. Wilcoxon rank-sum (Mann-Whitney) tests (per variant, across localities)
# ------------------------------------------------------------
fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
sig_label <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")

wilcox_H2 <- H2_local_all %>%
  group_by(variant) %>%
  group_map(function(df, keys) {
    wide <- df %>% pivot_wider(names_from = network, values_from = H2)
    wt <- tryCatch(
      wilcox.test(wide$`Parasitoid-Caterpillar`, wide$`Caterpillar-Plant`, paired = FALSE, exact = FALSE),
      error = function(e) NULL)
    if (is.null(wt)) {
      return(tibble(variant = keys$variant, comparison = "H2' per locality (unpaired)",
                    W_statistic = NA_real_, p_value = NA_character_, significance = NA_character_))
    }
    tibble(variant = keys$variant, comparison = "H2' per locality (unpaired)",
           W_statistic = round(unname(wt$statistic), 1),
           p_value = fmt_p(wt$p.value),
           significance = sig_label(wt$p.value))
  }) %>% bind_rows()

cat("--- Wilcoxon rank-sum (H2', par-cat vs cat-plant) ---\n"); print(wilcox_H2); cat("\n")

# ------------------------------------------------------------
# 8. Export Table S10 (Word + csv)
# ------------------------------------------------------------

write.csv(H2_local_summary, "output/rds/Table_S6_H2_summary.csv",  row.names = FALSE)
write.csv(wilcox_H2,        "output/rds/Table_S6_H2_wilcoxon.csv", row.names = FALSE)
write.csv(H2_local_all,     "output/rds/Table_S6_H2_values.csv",   row.names = FALSE)

doc_s6 <- officer::read_docx() %>%
  officer::body_add_par("Table S6: Network-level specialization (H2') per locality", style = "heading 1") %>%
  officer::body_add_par("a) Mean H2' (+/- SD) per network and variant", style = "heading 2") %>%
  flextable::body_add_flextable(flextable(H2_local_summary)) %>%
  officer::body_add_par("b) Wilcoxon rank-sum (Mann-Whitney) tests (per variant)", style = "heading 2") %>%
  flextable::body_add_flextable(flextable(wilcox_H2))
print(doc_s6, target = "output/Table_S6_specialization.docx")
