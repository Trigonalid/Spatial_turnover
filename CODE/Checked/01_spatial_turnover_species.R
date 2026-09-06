# ============================================================
# 01_spatial_turnover_species.R   >>> ANALYSIS ONLY <<<
# Analysis 1: spatial turnover of species (Bray-Curtis, Chao-Sorensen, Sorensen)
# Analysis 3: subsampling of caterpillars to parasitoid sample size
# Outputs: Table S2 (Word) + output/rds/turnover_long.rds (for plotting)
# Figures are built separately in 01_spatial_turnover_species_plots.R
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(dplyr)       # data wrangling
library(tidyr)       # pivot_wider / pivot_longer
library(purrr)       # map functions for downsampling loop
library(tibble)      # tibble, column_to_rownames
library(readxl)      # read_excel
library(here)        # here() for file paths
library(reshape2)    # melt() for matrices
library(vegan)       # vegdist()
library(betapart)    # beta.pair() for Sorensen
library(CommEcol)    # dis.chao() for Chao-Sorensen
library(flextable)   # flextable() for Word table export
library(officer)     # read_docx() for Word export

# ------------------------------------------------------------
# 1. Load data + create guild (loaded once, reused throughout)
# ------------------------------------------------------------
# MASTER: one row = one rearing event (caterpillar + plant + optional parasitoid).
# Excel columns are renamed to the short codes used below.
# guild is already present in MASTER (values "PAR" / "CAT") - not recreated here.
# Ohu2 excluded (temporal replicate, handled in a separate script).
MASTER <- read_excel(
    here("DATA/MASTER.xlsx"),
    guess_max = 1048576    # scan ALL rows when guessing column types.
  ) %>%                    # PAR_species_code is sparse; with the default guess it
                           # was typed as logical and its text codes coerced to NA.
  as_tibble() %>%
  rename(
    PLANT_sp = PLANT_species_code,
    CAT_sp   = CAT_scientific_name,
    PAR_sp   = PAR_species_code
  ) %>%
  mutate(
    PAR_sp = na_if(as.character(PAR_sp), "")    # ensure character; empty -> NA
  ) %>%                                          # NOTE: guild already exists in MASTER
  filter(locality != "Ohu2")

# Create output directories if missing
dir.create("output",     showWarnings = FALSE, recursive = TRUE)
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)

# NOTE: geographic distances (Distance.csv) were only needed for Mantel,
# which is no longer reported. Loading kept commented out in case a
# distance-decay analysis is added later.
# Distance <- read.csv2(here("DATA/Distance.csv"), row.names = 1)

# ------------------------------------------------------------
# 2. Community matrices (site x species)
# ------------------------------------------------------------
# --- Parasitoids --- (table() ignores NA, so only rows with PAR_sp are counted)
para_matrix    <- as.matrix(table(MASTER$locality, MASTER$PAR_sp))
para_matrix_PA <- (para_matrix > 0) * 1L   # presence/absence for Sorensen

# --- Caterpillars --- (every row has a CAT_sp)
cat_matrix    <- as.matrix(table(MASTER$locality, MASTER$CAT_sp))
cat_matrix_PA <- (cat_matrix > 0) * 1L     # presence/absence for Sorensen

# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------
# Melt a dissimilarity matrix to long format, keeping the UPPER triangle
# only (no diagonal) -> each site-pair appears exactly once.
matrix_to_long <- function(mat, guild_name, index_name) {

  mat[lower.tri(mat, diag = TRUE)] <- NA   # drops self-pairs and duplicate B-A

  reshape2::melt(mat, varnames = c("Locality_A", "Locality_B"),
                 value.name = "dissimilarity") %>%
    as_tibble() %>%
    filter(!is.na(dissimilarity)) %>%
    mutate(
      Locality_A = as.character(Locality_A),
      Locality_B = as.character(Locality_B),
      guild      = guild_name,
      index      = index_name
    )
}

# Summary row (mean +/- sd) for the summary table
summary_row <- function(index_name, guild_name, values_vec, dataset = "All data") {
  data.frame(
    index     = index_name,
    Community = guild_name,
    mean      = round(mean(values_vec, na.rm = TRUE), 3),
    sd        = round(sd(values_vec,   na.rm = TRUE), 3),
    Dataset   = dataset,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 4. BRAY-CURTIS DISSIMILARITY
# Abundance-based. {vegan} vegdist(method = "bray").
# ============================================================
bc_para <- vegan::vegdist(para_matrix, method = "bray") %>% as.matrix()
bc_cat  <- vegan::vegdist(cat_matrix,  method = "bray") %>% as.matrix()

BC_par <- matrix_to_long(bc_para, "Parasitoids",  "Bray-Curtis")
BC_cat <- matrix_to_long(bc_cat,  "Caterpillars", "Bray-Curtis")

result_bc_para <- summary_row("Bray-Curtis", "Parasitoids",  BC_par$dissimilarity)
result_bc_cat  <- summary_row("Bray-Curtis", "Caterpillars", BC_cat$dissimilarity)

wilcox_bc <- wilcox.test(BC_par$dissimilarity, BC_cat$dissimilarity); print(wilcox_bc)

# ============================================================
# 5. CHAO-SORENSEN DISSIMILARITY
# Corrects for unsampled rare species. {CommEcol} dis.chao()
# ============================================================
cs_para <- CommEcol::dis.chao(para_matrix, index = "sorensen", version = "rare") %>% as.matrix()
cs_cat  <- CommEcol::dis.chao(cat_matrix,  index = "sorensen", version = "rare") %>% as.matrix()

CS_par <- matrix_to_long(cs_para, "Parasitoids",  "Chao-Sorensen")
CS_cat <- matrix_to_long(cs_cat,  "Caterpillars", "Chao-Sorensen")

result_cs_para <- summary_row("Chao-Sorensen", "Parasitoids",  CS_par$dissimilarity)
result_cs_cat  <- summary_row("Chao-Sorensen", "Caterpillars", CS_cat$dissimilarity)

wilcox_cs <- wilcox.test(CS_par$dissimilarity, CS_cat$dissimilarity); print(wilcox_cs)

# ============================================================
# 6. SORENSEN DISSIMILARITY
# Presence/absence based. {betapart} beta.pair().
# Baseline; less robust with many rare species.
# ============================================================
sor_para <- betapart::beta.pair(para_matrix_PA)$beta.sor %>% as.matrix()
sor_cat  <- betapart::beta.pair(cat_matrix_PA)$beta.sor  %>% as.matrix()

SOR_par <- matrix_to_long(sor_para, "Parasitoids",  "Sorensen")
SOR_cat <- matrix_to_long(sor_cat,  "Caterpillars", "Sorensen")

result_sor_para <- summary_row("Sorensen", "Parasitoids",  SOR_par$dissimilarity)
result_sor_cat  <- summary_row("Sorensen", "Caterpillars", SOR_cat$dissimilarity)

wilcox_sor <- wilcox.test(SOR_par$dissimilarity, SOR_cat$dissimilarity); print(wilcox_sor)

# ============================================================
# 7. DOWNSAMPLING — caterpillars to parasitoid sample size
# At each site, caterpillar individuals (= rearing events) are randomly
# resampled WITH replacement to match the number of parasitoid rearings
# at that site. Repeated 999x; mean dissimilarity reported. set.seed(1234).
# ============================================================
set.seed(1234)
n_rand      <- 999
site_levels <- sort(unique(MASTER$locality))   # canonical order (= row order of table())

# Target sample size = number of parasitoid rearings per site
n_par_per_site <- MASTER %>%
  filter(guild == "PAR") %>%
  group_by(locality) %>%
  summarise(N = n(), .groups = "drop")

res_bc  <- vector("list", n_rand)
res_cs  <- vector("list", n_rand)
res_sor <- vector("list", n_rand)

for (i in seq_len(n_rand)) {

  # Resample caterpillar rows at each site to the target N (with replacement)
  resampled <- purrr::map_df(site_levels, function(loc) {
    N_target <- n_par_per_site %>% filter(locality == loc) %>% pull(N)
    if (length(N_target) == 0 || is.na(N_target) || N_target == 0) return(tibble())
    MASTER %>% filter(locality == loc) %>% sample_n(size = N_target, replace = TRUE)
  })

  # Wide abundance matrix; rows forced to canonical order
  df_wide <- resampled %>%
    count(locality, CAT_sp) %>%
    pivot_wider(names_from = CAT_sp, values_from = n, values_fill = 0) %>%
    column_to_rownames("locality")
  df_wide <- df_wide[site_levels, , drop = FALSE]

  df_pa <- (df_wide > 0) * 1L

  res_bc[[i]]  <- vegan::vegdist(df_wide, method = "bray") %>% as.matrix()
  res_cs[[i]]  <- CommEcol::dis.chao(df_wide, index = "sorensen", version = "rare") %>% as.matrix()
  res_sor[[i]] <- betapart::beta.pair(df_pa)$beta.sor %>% as.matrix()
}

# Element-wise mean across iterations (3D array)
mean_matrix <- function(res_list) {
  arr <- array(do.call(cbind, res_list),
               dim = c(dim(res_list[[1]]), length(res_list)))
  apply(arr, c(1, 2), mean, na.rm = TRUE)
}

bc_sub_mean  <- mean_matrix(res_bc)
cs_sub_mean  <- mean_matrix(res_cs)
sor_sub_mean <- mean_matrix(res_sor)

# Assign locality names to the averaged matrices
loc_names <- site_levels
dimnames(bc_sub_mean)  <- list(loc_names, loc_names)
dimnames(cs_sub_mean)  <- list(loc_names, loc_names)
dimnames(sor_sub_mean) <- list(loc_names, loc_names)

BC_sub  <- matrix_to_long(bc_sub_mean,  "Caterpillars-subsampled", "Bray-Curtis")
CS_sub  <- matrix_to_long(cs_sub_mean,  "Caterpillars-subsampled", "Chao-Sorensen")
SOR_sub <- matrix_to_long(sor_sub_mean, "Caterpillars-subsampled", "Sorensen")

result_bc_sub  <- summary_row("Bray-Curtis",   "Caterpillars", BC_sub$dissimilarity,  "Subsampled")
result_cs_sub  <- summary_row("Chao-Sorensen", "Caterpillars", CS_sub$dissimilarity,  "Subsampled")
result_sor_sub <- summary_row("Sorensen",      "Caterpillars", SOR_sub$dissimilarity, "Subsampled")

# ============================================================
# 8. Save long-format results for plotting (01_..._plots.R)
# ============================================================
turnover_long <- bind_rows(
  BC_par,  BC_cat,  BC_sub,
  CS_par,  CS_cat,  CS_sub,
  SOR_par, SOR_cat, SOR_sub
)
saveRDS(turnover_long, "output/rds/turnover_long.rds")

# ============================================================
# 9. SUMMARY TABLE — dissimilarity indices (Table S2)
# ============================================================
all_results <- bind_rows(
  result_bc_para,  result_bc_cat,  result_bc_sub,
  result_cs_para,  result_cs_cat,  result_cs_sub,
  result_sor_para, result_sor_cat, result_sor_sub
)
print(all_results)
write.csv(all_results, "output/rds/summary_dissimilarity.csv", row.names = FALSE)

ft  <- flextable(all_results)
doc <- officer::read_docx() %>%
  officer::body_add_par("Table S2: Summary of dissimilarity indices",
                        style = "heading 1") %>%
  flextable::body_add_flextable(ft)
print(doc, target = "output/Table_S2_dissimilarity_summary_new.docx")

# ============================================================
# 10. WILCOXON TABLE — species spatial turnover (unpaired)
#     Three comparisons per index (parallels rarity S3b / network S4b):
#       Parasitoids vs Caterpillars
#       Parasitoids vs Caterpillars-subsampled
#       Caterpillars vs Caterpillars-subsampled
# ============================================================
fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

wilcox_pairs <- list(
  c("Parasitoids",  "Caterpillars"),
  c("Parasitoids",  "Caterpillars-subsampled"),
  c("Caterpillars", "Caterpillars-subsampled")
)
indices_vec <- c("Bray-Curtis", "Chao-Sorensen", "Sorensen")

wilcox_species <- dplyr::bind_rows(lapply(indices_vec, function(idx) {
  dplyr::bind_rows(lapply(wilcox_pairs, function(pr) {
    a <- turnover_long %>% dplyr::filter(index == idx, guild == pr[1]) %>%
      dplyr::select(Locality_A, Locality_B, va = dissimilarity)
    b <- turnover_long %>% dplyr::filter(index == idx, guild == pr[2]) %>%
      dplyr::select(Locality_A, Locality_B, vb = dissimilarity)
    m <- dplyr::inner_join(a, b, by = c("Locality_A", "Locality_B"))
    if (nrow(m) < 2) {
      return(data.frame(index = idx, group_1 = pr[1], group_2 = pr[2],
                        n_pairs = nrow(m), W_statistic = NA_real_,
                        p_value = NA_character_, stringsAsFactors = FALSE))
    }
    wt <- wilcox.test(m$va, m$vb, paired = FALSE, exact = FALSE)
    data.frame(index = idx, group_1 = pr[1], group_2 = pr[2],
               n_pairs = nrow(m),
               W_statistic = round(unname(wt$statistic), 1),
               p_value = fmt_p(wt$p.value), stringsAsFactors = FALSE)
  }))
}))

print(wilcox_species)
write.csv(wilcox_species, "output/rds/Table_S2b_species_wilcoxon.csv", row.names = FALSE)

doc_w <- officer::read_docx() %>%
  officer::body_add_par("Table S2b: Wilcoxon rank-sum (Mann-Whitney) tests - species spatial turnover",
                        style = "heading 1") %>%
  flextable::body_add_flextable(autofit(flextable(wilcox_species)))
print(doc_w, target = "output/Table_S2b_species_wilcoxon_new.docx")
