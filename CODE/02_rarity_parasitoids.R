# ============================================================
# 02_rarity_parasitoids.R   >>> ANALYSIS ONLY <<<
# Study: Approaching the limit of diversity: spatial turnover
#        amplifies with trophic level (Libra et al.)
#
# Analysis 2: robustness to rare vs. common species.
# Recomputes spatial turnover (Bray-Curtis, Chao-Sorensen, Sorensen)
# for parasitoid rarity groups (All / Common / Rare) across three
# datasets (all data / common plants / common caterpillars).
#
# Rarity thresholds (across study area):
#   Rare parasitoids    : <= 9 reared specimens
#   Common parasitoids  : >= 10 reared specimens
#   Common caterpillars : >= 50 individuals
#   Common plants       : present at >= 5 localities (13 of 25 spp.)
#
# Outputs: Table S3 (summary) + Table S3b (Wilcoxon) as Word + csv,
#          and output/rds/rarity_long.rds for the plotting script.
# Figures are built separately in 02_rarity_parasitoids_plots.R (Fig. S5).
# ------------------------------------------------------------
# METHODS NOTES (consistent with 01):
#   - Bray-Curtis via {vegan}; Mantel / distance-decay NOT reported.
#   - Each site-pair counted once (upper triangle of the matrix).
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)       # data wrangling
library(purrr)       # map functions
library(tidyr)       # pivot functions
library(tibble)      # tibble
library(readxl)      # read_excel
library(here)        # here() for file paths
library(reshape2)    # melt()
library(vegan)       # vegdist()
library(betapart)    # beta.pair() for Sorensen
library(CommEcol)    # dis.chao() for Chao-Sorensen
library(flextable)   # flextable() for Word export
library(officer)     # read_docx()

# ------------------------------------------------------------
# 1. Load data (same robust load as 01)
# ------------------------------------------------------------
# guess_max scans all rows so the sparse PAR_species_code column is not
# mis-typed as logical. Ohu2 excluded (temporal replicate, analysed separately).
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

# Canonical locality order (used to align every subset matrix)
site_levels <- sort(unique(MASTER$locality))

# ------------------------------------------------------------
# 2. Classify species into rarity categories
# Computed directly from MASTER — not pre-assigned.
# ------------------------------------------------------------
# --- 2.1 Parasitoid rarity (total rearings per species) ---
par_counts <- MASTER %>%
  filter(!is.na(PAR_sp)) %>%
  group_by(PAR_sp) %>%
  summarise(n_rearings = n(), .groups = "drop") %>%
  mutate(par_rarity = if_else(n_rearings <= 9, "Rare", "Common"))

# --- 2.2 Caterpillar rarity (total individuals per species) ---
cat_counts <- MASTER %>%
  filter(!is.na(CAT_sp)) %>%
  group_by(CAT_sp) %>%
  summarise(n_individuals = n(), .groups = "drop") %>%
  mutate(common_caterpillar = n_individuals >= 50)

# --- 2.3 Plant rarity (number of localities per species) ---
plant_counts <- MASTER %>%
  filter(!is.na(PLANT_sp)) %>%
  group_by(PLANT_sp) %>%
  summarise(n_localities = n_distinct(locality), .groups = "drop") %>%
  mutate(common_plant = n_localities >= 5)

# --- 2.4 Join categories back to MASTER ---
MASTER <- MASTER %>%
  left_join(par_counts   %>% select(PAR_sp,   par_rarity),         by = "PAR_sp") %>%
  left_join(cat_counts   %>% select(CAT_sp,   common_caterpillar), by = "CAT_sp") %>%
  left_join(plant_counts %>% select(PLANT_sp, common_plant),       by = "PLANT_sp")

cat("\n--- Parasitoid rarity summary ---\n");        print(par_counts   %>% count(par_rarity))
cat("\n--- Common caterpillars (>= 50 ind.) ---\n"); print(cat_counts   %>% count(common_caterpillar))
cat("\n--- Common plants (>= 5 localities) ---\n");  print(plant_counts %>% count(common_plant))

# ------------------------------------------------------------
# 3. Filtered datasets
# ------------------------------------------------------------
MASTER_common_plants <- MASTER %>% filter(common_plant       == TRUE)
MASTER_common_cats   <- MASTER %>% filter(common_caterpillar == TRUE)

# ------------------------------------------------------------
# 4. Helpers: community matrix + dissimilarity
# ------------------------------------------------------------
# Rows = localities (aligned to site_levels), columns = parasitoid species.
build_para_matrix <- function(data) {
  mat <- as.matrix(table(data$locality, data$PAR_sp))
  
  missing <- setdiff(site_levels, rownames(mat))
  if (length(missing) > 0) {
    empty <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                    dimnames = list(missing, colnames(mat)))
    mat <- rbind(mat, empty)
  }
  mat[site_levels, , drop = FALSE]
}

compute_dissimilarity <- function(mat, method) {
  if (method == "bray") {
    vegan::vegdist(mat, method = "bray") %>% as.matrix()
  } else if (method == "chao-sorensen") {
    CommEcol::dis.chao(mat, index = "sorensen", version = "rare") %>% as.matrix()
  } else if (method == "sorensen") {
    mat_pa <- (mat > 0) * 1L
    betapart::beta.pair(mat_pa)$beta.sor %>% as.matrix()
  }
}

# ------------------------------------------------------------
# 5. Run one group x dataset x index combination
# Returns a long-format data frame (upper triangle) + a summary row.
# ------------------------------------------------------------
run_analysis <- function(data, group_label, dataset_label, method) {
  
  par_data <- data %>% filter(!is.na(PAR_sp))
  if (group_label == "Common") {
    par_data <- par_data %>% filter(par_rarity == "Common")
  } else if (group_label == "Rare") {
    par_data <- par_data %>% filter(par_rarity == "Rare")
  }
  # "All" uses unfiltered par_data
  
  mat  <- build_para_matrix(par_data)
  dmat <- compute_dissimilarity(mat, method)
  
  # Long format: keep upper triangle only -> each site-pair once
  dmat[lower.tri(dmat, diag = TRUE)] <- NA
  long <- reshape2::melt(dmat,
                         varnames   = c("Locality_A", "Locality_B"),
                         value.name = "dissimilarity") %>%
    as_tibble() %>%
    filter(!is.na(dissimilarity)) %>%    # drops self-pairs, duplicates and NaN
    mutate(
      Locality_A  = as.character(Locality_A),
      Locality_B  = as.character(Locality_B),
      guild_clean = group_label,
      dataset     = dataset_label,
      index_type  = method
    )
  
  summary <- tibble(
    index_type         = method,
    dataset            = dataset_label,
    group              = group_label,
    mean_dissimilarity = round(mean(long$dissimilarity), 3),
    sd_dissimilarity   = round(sd(long$dissimilarity),   3),
    n_species          = ncol(mat)
  )
  
  list(long = long, summary = summary)
}

# ------------------------------------------------------------
# 6. Run all combinations (3 indices x 3 datasets x 3 groups = 27)
# ------------------------------------------------------------
methods  <- c("bray", "chao-sorensen", "sorensen")
datasets <- list(
  "All data"                 = MASTER,
  "Common Plants Only"       = MASTER_common_plants,
  "Common Caterpillars Only" = MASTER_common_cats
)
groups <- c("All", "Common", "Rare")

all_long    <- list()
all_summary <- list()

for (method in methods) {
  for (ds_name in names(datasets)) {
    for (grp in groups) {
      res <- run_analysis(datasets[[ds_name]], grp, ds_name, method)
      all_long[[length(all_long) + 1]]       <- res$long
      all_summary[[length(all_summary) + 1]] <- res$summary
    }
  }
}

rarity_long   <- bind_rows(all_long)      # character columns; factors set in plots script
summary_table <- bind_rows(all_summary)

saveRDS(rarity_long, "output/rds/rarity_long.rds")

# ------------------------------------------------------------
# 7. Wilcoxon tests: pairwise rarity-group comparisons
# per index x dataset combination. PAIRED signed-rank tests: dissimilarity
# values are matched by site-pair (Locality_A, Locality_B), since the rarity
# groups are compared over the same pairs of localities. Statistic is V.
# (Only site-pairs present in both groups are used -> n_pairs.)
# ------------------------------------------------------------
wilcox_pairs <- list(c("All", "Common"), c("All", "Rare"), c("Common", "Rare"))

fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

wilcox_results <- rarity_long %>%
  group_by(index_type, dataset) %>%
  group_map(function(df, keys) {
    purrr::map_df(wilcox_pairs, function(pair) {
      a <- df %>% filter(guild_clean == pair[1]) %>%
        dplyr::select(Locality_A, Locality_B, value_a = dissimilarity)
      b <- df %>% filter(guild_clean == pair[2]) %>%
        dplyr::select(Locality_A, Locality_B, value_b = dissimilarity)
      paired_df <- dplyr::inner_join(a, b, by = c("Locality_A", "Locality_B"))
      if (nrow(paired_df) < 2) {
        return(tibble(index_type = keys$index_type, dataset = keys$dataset,
                      group_1 = pair[1], group_2 = pair[2], n_pairs = nrow(paired_df),
                      V_statistic = NA_real_, p_value = NA_character_,
                      significance = NA_character_))
      }
      wt <- wilcox.test(paired_df$value_a, paired_df$value_b, paired = TRUE, exact = FALSE)
      tibble(
        index_type   = keys$index_type,
        dataset      = keys$dataset,
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
  }) %>%
  bind_rows()

print(wilcox_results)

# ------------------------------------------------------------
# 8. Export tables (Word + csv)
# ------------------------------------------------------------
# Table S3 — summary of dissimilarity indices by rarity group
write.csv(summary_table, "output/rds/Table_S3_rarity_summary.csv", row.names = FALSE)
doc_summary <- officer::read_docx() %>%
  officer::body_add_par("Table S3: Summary of dissimilarity indices by rarity group",
                        style = "heading 1") %>%
  flextable::body_add_flextable(flextable(summary_table))
print(doc_summary, target = "output/Table_S3_rarity_summary.docx")

# Table S3b — Wilcoxon rarity-group comparisons
write.csv(wilcox_results, "output/rds/Table_S3b_rarity_wilcoxon.csv", row.names = FALSE)
doc_wilcox <- officer::read_docx() %>%
  officer::body_add_par("Table S3b: Paired Wilcoxon signed-rank tests - rarity group comparisons",
                        style = "heading 1") %>%
  flextable::body_add_flextable(flextable(wilcox_results))
print(doc_wilcox, target = "output/Table_S3b_rarity_wilcoxon.docx")