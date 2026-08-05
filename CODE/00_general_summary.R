# ============================================================
# 00_table_s1_dataset_summary.R   >>> TABLE S1 <<<
# Study: Spatial turnover amplifies with trophic level in
#        hyperdiverse food webs (Libra et al.)
#
# Table S1: overview of the reared dataset -- abundance, species richness,
# subfamily richness, and family richness for caterpillars and parasitoids.
#   Caterpillar abundance = number of reared specimens (rows).
#   Parasitoid abundance  = sum(PAR_abundance) (total parasitoid individuals).
#
# Reference (main text): caterpillars ~27,537 / 422 sp.; parasitoids 304 sp.
# ============================================================

library(dplyr)
library(tibble)
library(readxl)
library(here)
library(flextable)
library(officer)

# Table S1 describes the FULL reared dataset -> do not drop Ohu2 by default.
FILTER_OHU2 <- FALSE

# ------------------------------------------------------------
# 1. Load data (robust load; keep taxonomy columns)
# ------------------------------------------------------------
MASTER <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(
    CAT_sp = CAT_scientific_name,
    PAR_sp = PAR_species_code
  ) %>%
  mutate(PAR_sp = na_if(as.character(PAR_sp), ""))

if (FILTER_OHU2) MASTER <- MASTER %>% filter(locality != "Ohu2")

# ------------------------------------------------------------
# 2. Counts
# ------------------------------------------------------------
n_distinct_nona <- function(x) n_distinct(x[!is.na(x)])

cat_rows <- MASTER %>% filter(!is.na(CAT_sp))
par_rows <- MASTER %>% filter(guild == "PAR")   # parasitoids defined by guild

cat("Guild counts:\n"); print(count(MASTER, guild)); cat("\n")

# Parasitoid abundance from PAR_abundance over guild == "PAR" rows
if ("PAR_abundance" %in% names(MASTER)) {
  par_individuals <- sum(par_rows$PAR_abundance, na.rm = TRUE)
} else {
  message("PAR_abundance column not found; using number of guild=='PAR' rows instead.")
  par_individuals <- nrow(par_rows)
}

table_s1 <- tibble(
  Group          = c("Caterpillars (Lepidoptera)", "Parasitoids"),
  Individuals    = c(nrow(cat_rows), par_individuals),
  Rearing_events = c(nrow(cat_rows), nrow(par_rows)),
  Families       = c(n_distinct_nona(cat_rows$CAT_family),    n_distinct_nona(par_rows$PAR_family)),
  Subfamilies    = c(n_distinct_nona(cat_rows$CAT_subfamily), n_distinct_nona(par_rows$PAR_Subfamily)),
  Species        = c(n_distinct_nona(cat_rows$CAT_sp),        n_distinct_nona(par_rows$PAR_sp))
)

# ------------------------------------------------------------
# 3. Report + sanity check
# ------------------------------------------------------------
cat("--- Table S1 ---\n")
print(table_s1)
cat("\nExpected (main text): caterpillars 27,537 individuals / 422 species; parasitoids 304 species.\n")
cat("(Parasitoid 'Individuals' = sum of PAR_abundance; 'Rearing_events' = rows with a parasitoid,\n")
cat(" i.e. parasitized caterpillars, expected ~2,875.)\n")

# ------------------------------------------------------------
# 4. Export (Word + csv)
# ------------------------------------------------------------
dir.create("output",     showWarnings = FALSE, recursive = TRUE)
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)

write.csv(table_s1, "output/rds/Table_S1_dataset_summary.csv", row.names = FALSE)

doc_s1 <- officer::read_docx() %>%
  officer::body_add_par(
    paste("Table S1. Overview of the reared dataset: abundance, species richness,",
          "subfamily richness, and family richness of caterpillars and parasitoids."),
    style = "heading 2") %>%
  flextable::body_add_flextable(flextable(table_s1))
print(doc_s1, target = "output/Table_S1_dataset_summary.docx")