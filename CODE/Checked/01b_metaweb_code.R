# ============================================================
# 01b_meta_food_web_plot.R
# Fig. 1B: combined tritrophic metaweb (plant -> caterpillar -> parasitoid),
# pooled across all eight localities, via bipartite::plotweb2.
#
# Colours matched to Fig. 1C/1D:
#   caterpillars    = #F25E0D  (col.low / col.pred2)
#   parasitoids     = #EDB439  (col.pred2 upper level)
#   plant-cat links = #B4464B  (col.interaction, web 1)
#   cat-par links   = #4682B4  (col.interaction2, web 2)
#   plants          = green    (col.low, web 1)
# ============================================================

library(tidyverse)
library(bipartite)
library(here)

# ------------------------------------------------------------
# 1. Load data (same robust load as 01/02/04)
# ------------------------------------------------------------
raw_data <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(PLANT_sp = PLANT_species_code,
         CAT_sp   = CAT_scientific_name,
         PAR_sp   = PAR_species_code) %>%
  mutate(PAR_sp = na_if(as.character(PAR_sp), "")) %>%
  filter(locality != "Ohu2")

dir.create("output/fig", showWarnings = FALSE, recursive = TRUE)

use_freq <- TRUE   # TRUE = interaction counts (abundance), FALSE = presence/absence

# ------------------------------------------------------------
# 2. Build pooled matrices (all localities combined)
# ------------------------------------------------------------
plant_caterpillar <- raw_data %>%
  filter(!is.na(PLANT_sp) & !is.na(CAT_sp)) %>%
  { if (use_freq) group_by(., PLANT_sp, CAT_sp) %>% summarise(freq = n(), .groups = "drop")
    else distinct(., PLANT_sp, CAT_sp) %>% mutate(freq = 1) } %>%
  pivot_wider(names_from = CAT_sp, values_from = freq, values_fill = 0) %>%
  column_to_rownames("PLANT_sp") %>%
  as.matrix()

caterpillar_parasitoid <- raw_data %>%
  filter(!is.na(CAT_sp) & !is.na(PAR_sp)) %>%
  { if (use_freq) group_by(., CAT_sp, PAR_sp) %>% summarise(freq = n(), .groups = "drop")
    else distinct(., CAT_sp, PAR_sp) %>% mutate(freq = 1) } %>%
  pivot_wider(names_from = PAR_sp, values_from = freq, values_fill = 0) %>%
  column_to_rownames("CAT_sp") %>%
  as.matrix()

if (is.null(dim(plant_caterpillar))     || any(dim(plant_caterpillar) == 0) ||
    is.null(dim(caterpillar_parasitoid)) || any(dim(caterpillar_parasitoid) == 0)) {
  stop("One of the matrices is empty — check raw_data (PLANT_sp, CAT_sp, PAR_sp).")
}

rownames(plant_caterpillar)      <- make.unique(rownames(plant_caterpillar))
colnames(plant_caterpillar)      <- make.unique(colnames(plant_caterpillar))
rownames(caterpillar_parasitoid) <- make.unique(rownames(caterpillar_parasitoid))
colnames(caterpillar_parasitoid) <- make.unique(colnames(caterpillar_parasitoid))

# ------------------------------------------------------------
# 3. Unify the middle level (caterpillars) across both matrices
#    -> every caterpillar species appears as a column in web 1
#       AND as a row in web 2, even if only recorded in one of them
# ------------------------------------------------------------
all_caterpillars <- sort(unique(c(colnames(plant_caterpillar),
                                  rownames(caterpillar_parasitoid))))

plant_species <- rownames(plant_caterpillar)
unified_plant_cat <- matrix(0, nrow = length(plant_species), ncol = length(all_caterpillars),
                            dimnames = list(plant_species, all_caterpillars))
unified_plant_cat[plant_species, colnames(plant_caterpillar)] <- plant_caterpillar
plant_caterpillar <- unified_plant_cat

parasitoid_species <- colnames(caterpillar_parasitoid)
unified_cat_par <- matrix(0, nrow = length(all_caterpillars), ncol = length(parasitoid_species),
                          dimnames = list(all_caterpillars, parasitoid_species))
unified_cat_par[rownames(caterpillar_parasitoid), parasitoid_species] <- caterpillar_parasitoid
caterpillar_parasitoid <- unified_cat_par

# ------------------------------------------------------------
# 4. Export pooled interaction matrices (for reference / SI)
# ------------------------------------------------------------
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)
write.csv(as.data.frame(plant_caterpillar),
          "output/rds/Figure_1B_plant_caterpillar_matrix.csv", row.names = TRUE)
write.csv(as.data.frame(caterpillar_parasitoid),
          "output/rds/Figure_1B_caterpillar_parasitoid_matrix.csv", row.names = TRUE)

# ------------------------------------------------------------
# 5. Hide outer-level labels (plants, parasitoids); keep caterpillar
#    names so plotweb2 can link web 1 and web 2 correctly
# ------------------------------------------------------------
rownames(plant_caterpillar)         <- rep(" ", nrow(plant_caterpillar))
colnames(caterpillar_parasitoid)    <- rep(" ", ncol(caterpillar_parasitoid))

# ------------------------------------------------------------
# 6. Colours matched to Fig. 1C/1D  >>> DO NOT CHANGE <
# ------------------------------------------------------------
col_link_plant_cat <- "#B4464B"   # plant -> caterpillar links (Fig. 1D)
col_link_cat_par   <- "#4682B4"   # caterpillar -> parasitoid links (Fig. 1D)
col_plants         <- "forestgreen"
col_caterpillars   <- "#F25E0D"   # Fig. 1C
col_parasitoids    <- "#EDB439"   # Fig. 1C

# ------------------------------------------------------------
# 7. Plot combined tritrophic metaweb (Fig. 1B)
# ------------------------------------------------------------
plot_metaweb <- function() {
  plotweb2(
    web  = plant_caterpillar,
    web2 = caterpillar_parasitoid,
    method  = "cca",
    method2 = "cca",
    empty  = FALSE,
    empty2 = FALSE,
    labsize = 0.01,
    ybig = 1,
    y_width  = 0.05,
    spacing  = 0.05,
    spacing2 = 0.05,
    arrow  = "down",
    arrow2 = "up",
    col.interaction  = col_link_plant_cat,
    col.interaction2 = col_link_cat_par,
    low.abun.col   = col_plants,
    high.abun.col  = col_caterpillars,
    low.abun.col2  = col_caterpillars,
    high.abun.col2 = col_parasitoids,
    col.prey  = col_plants,
    col.pred  = col_caterpillars,
    col.prey2 = col_caterpillars,
    col.pred2 = col_parasitoids,
    lab.space = 1
  )
}

svg("output/fig/Figure_1B_metaweb.svg", width = 12, height = 8)
plot_metaweb(); dev.off()

png("output/fig/Figure_1B_metaweb.png", width = 2400, height = 1600, res = 200)
plot_metaweb(); dev.off()

pdf("output/fig/Figure_1B_metaweb.pdf", width = 12, height = 8)
plot_metaweb(); dev.off()

message("Fig. 1B saved to output/fig/Figure_1B_metaweb.{svg,pdf,png}")