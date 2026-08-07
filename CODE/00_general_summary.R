# ============================================================
# 00b_table_s1_community_composition.R   >>> TABLE S1 (per locality) <<<
# Study: Spatial turnover amplifies with trophic level in
#        hyperdiverse food webs (Libra et al.)
#
# Community composition per locality: species richness + abundance of
# plants, caterpillars, and parasitoids, plus BIN counts for caterpillars
# and parasitoids. Ohu1 + Ohu2 are combined into "Ohu".
# Total row = overall distinct richness / summed abundance across ALL data
# (NOT the sum of per-locality richness, since species are shared).
# ============================================================

library(dplyr)
library(tibble)
library(readxl)
library(here)
library(flextable)
library(officer)

# ------------------------------------------------------------
# 1. Load data (robust load; Ohu kept and merged)
# ------------------------------------------------------------
MASTER <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(
    PLANT_sp = PLANT_species_code,
    CAT_sp   = CAT_scientific_name,
    PAR_sp   = PAR_species_code
  ) %>%
  mutate(
    PAR_sp = na_if(as.character(PAR_sp), ""),
    locality_grp = if_else(locality %in% c("Ohu1", "Ohu2"), "Ohu", locality)
  )

# ------------------------------------------------------------
# 2. Reusable summary (richness / abundance / BINs)
#    plant abundance = number of rearing records on a host plant
#    parasitoid abundance = sum(PAR_abundance) (total individuals)
# ------------------------------------------------------------
summarise_composition <- function(df) {
  df %>% summarise(
    plant_S  = n_distinct(PLANT_sp[!is.na(PLANT_sp)]),
    plant_N  = sum(!is.na(PLANT_sp)),
    cat_S    = n_distinct(CAT_sp[!is.na(CAT_sp)]),
    cat_N    = sum(!is.na(CAT_sp)),
    cat_BIN  = n_distinct(CAT_BIN[!is.na(CAT_BIN)]),
    par_S    = n_distinct(PAR_sp[!is.na(PAR_sp)]),
    par_N    = sum(PAR_abundance[!is.na(PAR_sp)], na.rm = TRUE),
    par_BIN  = n_distinct(PAR_BIN[!is.na(PAR_BIN)]),
    .groups  = "drop"
  )
}

per_loc <- MASTER %>%
  group_by(Locality = locality_grp) %>%
  summarise_composition() %>%
  arrange(Locality)

total_row <- MASTER %>%
  summarise_composition() %>%
  mutate(Locality = "Total") %>%
  select(Locality, everything())

table_s1 <- bind_rows(per_loc, total_row)

cat("--- Table S1: community composition per locality ---\n")
print(table_s1)

# ------------------------------------------------------------
# 3. Export: csv + Word (grouped header: Plants / Caterpillars / Parasitoids)
# ------------------------------------------------------------
dir.create("output",     showWarnings = FALSE, recursive = TRUE)
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)
write.csv(table_s1, "output/rds/Table_S1_community_composition.csv", row.names = FALSE)

caption_txt <- paste(
  "Tab. S1: Community composition across localities of Papua New Guinea.",
  "Values represent species richness and abundance of plants, caterpillars,",
  "and parasitoids, and the number of Barcode Index Numbers (BINs) for",
  "caterpillars and parasitoids. Total values indicate overall richness and",
  "abundance across all localities."
)

ft <- flextable(table_s1) %>%
  set_header_labels(
    Locality = "Locality",
    plant_S = "Richness", plant_N = "Abundance",
    cat_S = "Richness", cat_N = "Abundance", cat_BIN = "BINs",
    par_S = "Richness", par_N = "Abundance", par_BIN = "BINs"
  ) %>%
  add_header_row(values = c("", "Plants", "Caterpillars", "Parasitoids"),
                 colwidths = c(1, 2, 3, 3)) %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  bold(part = "header") %>%
  bold(i = ~ Locality == "Total") %>%
  autofit()

doc_s1 <- officer::read_docx() %>%
  officer::body_add_par(caption_txt, style = "Normal") %>%
  officer::body_add_par("") %>%
  flextable::body_add_flextable(ft)
print(doc_s1, target = "output/Table_S1_community_composition.docx")