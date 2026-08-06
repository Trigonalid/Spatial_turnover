# ============================================================
# food_web_metaweb.R  -- combined tritrophic metaweb (all localities)
# plant -> caterpillar -> parasitoid, via bipartite::plotweb2.
#
# Colours matched to the manuscript figures (Fig. 1C/1D):
#   plants        = green   (bottom level)
#   caterpillars  = #F25E0D (Fig. 1C caterpillars)
#   parasitoids   = #EDB439 (Fig. 1C parasitoids)
#   plant-cat links = #B4464B (Fig. 1D)
#   cat-par links   = #4682B4 (Fig. 1D)
# No black outlines: interaction borders match the link colour.
#
# ============================================================

library(tidyverse)
library(bipartite)
raw_data <- read_excel(here("DATA/MASTER.xlsx"), guess_max = 1048576) %>%
  as_tibble() %>%
  rename(PLANT_sp = PLANT_species_code,
         CAT_sp   = CAT_scientific_name,
         PAR_sp   = PAR_species_code) %>%
  mutate(PAR_sp = na_if(as.character(PAR_sp), "")) %>%
  filter(locality != "Ohu2")
out_dir <- "FW_outputs_metaweb"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

use_freq <- TRUE   # TRUE = interaction counts, FALSE = presence/absence
library(tidyverse) # Or dplyr, tidyr, and tibble if you prefer non-meta packages
library(bipartite)

##########
# --- Output folder ---
out_dir <- "FW_outputs_bipartite"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper to make safe filenames from locality values
safe_name <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^[:alnum:]_\\-]", "", x)
  if (nchar(x) == 0) x <- "unknown"
  x
}

# --- ASSUMPTION ---
# This code assumes 'raw_data' is already loaded and is a data frame with at least
# the columns: 'locality', 'PLANT_sp', 'CAT_sp', and 'PAR_sp'.
# If 'raw_data' is not loaded, the script will fail.

# Example placeholder for 'raw_data' structure if needed:
# raw_data <- data.frame(
#   locality = rep(c("Elem", "Forest"), each = 10),
#   PLANT_sp = sample(c("Oak", "Pine", NA), 20, replace = TRUE),
#   CAT_sp = sample(c("MothA", "MothB", NA), 20, replace = TRUE),
#   PAR_sp = sample(c("WaspX", "FlyY", NA), 20, replace = TRUE)
# )

localities <- sort(unique(raw_data$locality))

for (loc in localities) {
  message("Processing locality: ", loc)
  dat_loc <- raw_data %>% filter(locality == loc)
  
  # Build matrices for this locality
  plant_caterpillar_matrix <- dat_loc %>%
    filter(!is.na(PLANT_sp) & !is.na(CAT_sp)) %>%
    group_by(PLANT_sp, CAT_sp) %>%
    summarise(freq = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = CAT_sp, values_from = freq, values_fill = 0) %>%
    tibble::column_to_rownames(var = "PLANT_sp") %>%
    as.matrix()
  
  caterpillar_parasitoid_matrix <- dat_loc %>%
    filter(!is.na(CAT_sp) & !is.na(PAR_sp)) %>%
    group_by(CAT_sp, PAR_sp) %>%
    summarise(freq = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = PAR_sp, values_from = freq, values_fill = 0) %>%
    tibble::column_to_rownames(var = "CAT_sp") %>%
    as.matrix()
  
  # Skip if any matrix is empty
  if (is.null(dim(plant_caterpillar_matrix)) ||
      any(dim(plant_caterpillar_matrix) == 0) ||
      is.null(dim(caterpillar_parasitoid_matrix)) ||
      any(dim(caterpillar_parasitoid_matrix) == 0)) {
    message("  Skipped (no interactions for this locality).")
    next
  }
  
  # ------------------------------
  # 1) Ensure dimnames exist and are unique (internal names)
  # ------------------------------
  # row/col names can never be NULL for plotting internals
  if (is.null(rownames(plant_caterpillar_matrix))) {
    rownames(plant_caterpillar_matrix) <- paste0("PL_row_", seq_len(nrow(plant_caterpillar_matrix)))
  }
  if (is.null(colnames(plant_caterpillar_matrix))) {
    colnames(plant_caterpillar_matrix) <- paste0("PL_col_", seq_len(ncol(plant_caterpillar_matrix)))
  }
  if (is.null(rownames(caterpillar_parasitoid_matrix))) {
    rownames(caterpillar_parasitoid_matrix) <- paste0("CA_row_", seq_len(nrow(caterpillar_parasitoid_matrix)))
  }
  if (is.null(colnames(caterpillar_parasitoid_matrix))) {
    colnames(caterpillar_parasitoid_matrix) <- paste0("CA_col_", seq_len(ncol(caterpillar_parasitoid_matrix)))
  }
  
  # make unique to avoid duplicate-name pitfalls
  # IMPORTANT: These unique names must be preserved for the caterpillar level
  rownames(plant_caterpillar_matrix)      <- make.unique(rownames(plant_caterpillar_matrix))
  colnames(plant_caterpillar_matrix)      <- make.unique(colnames(plant_caterpillar_matrix))
  rownames(caterpillar_parasitoid_matrix) <- make.unique(rownames(caterpillar_parasitoid_matrix))
  colnames(caterpillar_parasitoid_matrix) <- make.unique(colnames(caterpillar_parasitoid_matrix))
  
  # -----------------------------------------------------------
  # 2) Replace visible labels with a single space " "
  #    -> ONLY FOR THE OUTER TROPHIC LEVELS (PLANTS & PARASITOIDS)
  #    -> The caterpillar names (middle level) MUST be kept unique for plotweb2
  # -----------------------------------------------------------
  
  # Plants (web 1, rows) -> Set to " "
  rownames(plant_caterpillar_matrix) <- rep(" ", nrow(plant_caterpillar_matrix))
  
  # Caterpillars (web 1, columns AND web 2, rows) -> Names kept (NO CHANGE HERE)
  # These unique names link the two webs together.
  
  # Parasitoids (web 2, columns) -> Set to " "
  colnames(caterpillar_parasitoid_matrix) <- rep(" ", ncol(caterpillar_parasitoid_matrix))
  
  # File paths
  base <- file.path(out_dir, paste0("tritrophic_", safe_name(loc)))
  png_file <- paste0(base, ".png")
  svg_file <- paste0(base, ".svg")
  
  # Optional: diagnostic dump... (commented out as in original)
  # saveRDS(list(loc=loc,
  #              dim_plant = dim(plant_caterpillar_matrix),
  #              dim_catpar = dim(caterpillar_parasitoid_matrix),
  #              rownames_plant = rownames(plant_caterpillar_matrix),
  #              colnames_plant = colnames(plant_caterpillar_matrix)),
  #         file = paste0(base, "_diag.rds"))
  
  # ---- Save SVG ----
  tryCatch({
    svg(svg_file, width = 12, height = 8)  # inches
    plotweb2(
      web  = plant_caterpillar_matrix,
      web2 = caterpillar_parasitoid_matrix,
      method = "cca",
      method2 = "cca",
      empty = FALSE,
      empty2 = FALSE,
      labsize = 0.01,      # tiny but labels are just spaces anyway
      ybig = 1,
      y_width = 0.05,
      spacing = 0.05,
      spacing2 = 0.05,
      arrow = "down",
      arrow2 = "up",
      col.interaction  = "grey80",
      col.interaction2 = "grey80",
      low.abun.col = "orange",   
      high.abun.col = "darkgreen",
      low.abun.col2  = "orange",
      high.abun.col2 = "darkgreen",
      col.prey  = "darkgreen",
      col.pred  = "orange",
      col.prey2 = "orange",
      col.pred2 = "lightblue",
      lab.space = 1
    )
    dev.off()
  }, error = function(e) {
    # Check if a device is open before calling dev.off()
    if (!is.null(dev.list())) dev.off()
    message("  ERROR while saving SVG for ", loc, " : ", conditionMessage(e))
    # continue to PNG attempt or next locality
  })
  
  # ---- Save PNG ----
  tryCatch({
    png(png_file, width = 2400, height = 1600, res = 200)  # pixels @ 200 dpi
    bipartite::plotweb2(
      web  = plant_caterpillar_matrix,
      web2 = caterpillar_parasitoid_matrix,
      method = "cca",
      method2 = "cca",
      empty = FALSE,
      empty2 = FALSE,
      labsize = 0.01,      # tiny but labels are just spaces anyway
      ybig = 1,
      y_width = 0.05,
      spacing = 0.05,
      spacing2 = 0.05,
      arrow = "down",
      arrow2 = "up",
      col.interaction  = "grey80",
      col.interaction2 = "grey80",
      low.abun.col = "orange",   
      high.abun.col = "darkgreen",
      low.abun.col2  = "orange",
      high.abun.col2 = "darkgreen",
      col.prey  = "darkgreen",
      col.pred  = "orange",
      col.prey2 = "orange",
      col.pred2 = "lightblue",
      lab.space = 1
    )
    dev.off()
  }, error = function(e) {
    # Check if a device is open before calling dev.off()
    if (!is.null(dev.list())) dev.off()
    message("  ERROR while saving PNG for ", loc, " : ", conditionMessage(e))
  })
  
  message("  Saved (or attempted): ", basename(svg_file), " & ", basename(png_file))
}



####################
library(tidyverse)
library(bipartite)

##########
# --- Output folder ---
out_dir <- "FW_outputs_bipartite_2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper to make safe filenames from locality values
safe_name <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^[:alnum:]_\\-]", "", x)
  if (nchar(x) == 0) x <- "unknown"
  x
}

# --- ASSUMPTION ---
# This code assumes 'raw_data' is already loaded and is a data frame with at least
# the columns: 'locality', 'PLANT_sp', 'CAT_sp', and 'PAR_sp'.

localities <- sort(unique(raw_data$locality))

for (loc in localities) {
  message("Processing locality: ", loc)
  dat_loc <- raw_data %>% filter(locality == loc)
  
  # Build matrices for this locality
  plant_caterpillar_matrix <- dat_loc %>%
    filter(!is.na(PLANT_sp) & !is.na(CAT_sp)) %>%
    group_by(PLANT_sp, CAT_sp) %>%
    summarise(freq = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = CAT_sp, values_from = freq, values_fill = 0) %>%
    tibble::column_to_rownames(var = "PLANT_sp") %>%
    as.matrix()
  
  caterpillar_parasitoid_matrix <- dat_loc %>%
    filter(!is.na(CAT_sp) & !is.na(PAR_sp)) %>%
    group_by(CAT_sp, PAR_sp) %>%
    summarise(freq = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = PAR_sp, values_from = freq, values_fill = 0) %>%
    tibble::column_to_rownames(var = "CAT_sp") %>%
    as.matrix()
  
  # Skip if any matrix is empty
  if (is.null(dim(plant_caterpillar_matrix)) ||
      any(dim(plant_caterpillar_matrix) == 0) ||
      is.null(dim(caterpillar_parasitoid_matrix)) ||
      any(dim(caterpillar_parasitoid_matrix) == 0)) {
    message("  Skipped (no interactions for this locality).")
    next
  }
  
  # ------------------------------
  # 1) Ensure dimnames exist and are unique (internal names)
  # ------------------------------
  if (is.null(rownames(plant_caterpillar_matrix))) {
    rownames(plant_caterpillar_matrix) <- paste0("PL_row_", seq_len(nrow(plant_caterpillar_matrix)))
  }
  if (is.null(colnames(plant_caterpillar_matrix))) {
    colnames(plant_caterpillar_matrix) <- paste0("PL_col_", seq_len(ncol(plant_caterpillar_matrix)))
  }
  if (is.null(rownames(caterpillar_parasitoid_matrix))) {
    rownames(caterpillar_parasitoid_matrix) <- paste0("CA_row_", seq_len(nrow(caterpillar_parasitoid_matrix)))
  }
  if (is.null(colnames(caterpillar_parasitoid_matrix))) {
    colnames(caterpillar_parasitoid_matrix) <- paste0("CA_col_", seq_len(ncol(caterpillar_parasitoid_matrix)))
  }
  
  # make unique to avoid duplicate-name pitfalls
  rownames(plant_caterpillar_matrix)      <- make.unique(rownames(plant_caterpillar_matrix))
  colnames(plant_caterpillar_matrix)      <- make.unique(colnames(plant_caterpillar_matrix))
  rownames(caterpillar_parasitoid_matrix) <- make.unique(rownames(caterpillar_parasitoid_matrix))
  colnames(caterpillar_parasitoid_matrix) <- make.unique(colnames(caterpillar_parasitoid_matrix))
  
  # ----------------------------------------------------------------------
  # 3) Sjednocení housenek (pro kompletní a spojitou střední úroveň) 🐛
  # ----------------------------------------------------------------------
  
  # 3.1) Identifikace všech housenek pro sjednocení dimenzí
  all_caterpillars <- sort(unique(c(
    colnames(plant_caterpillar_matrix),
    rownames(caterpillar_parasitoid_matrix)
  )))
  
  # 3.2) Vytvoření sjednocené Plant-Caterpillar matice (web)
  plant_species <- rownames(plant_caterpillar_matrix)
  unified_plant_cat_matrix <- matrix(0, 
                                     nrow = length(plant_species), 
                                     ncol = length(all_caterpillars),
                                     dimnames = list(plant_species, all_caterpillars))
  unified_plant_cat_matrix[plant_species, colnames(plant_caterpillar_matrix)] <- plant_caterpillar_matrix
  plant_caterpillar_matrix <- unified_plant_cat_matrix
  
  # 3.3) Vytvoření sjednocené Caterpillar-Parasitoid matice (web2)
  parasitoid_species <- colnames(caterpillar_parasitoid_matrix)
  unified_cat_par_matrix <- matrix(0, 
                                   nrow = length(all_caterpillars), 
                                   ncol = length(parasitoid_species),
                                   dimnames = list(all_caterpillars, parasitoid_species))
  unified_cat_par_matrix[rownames(caterpillar_parasitoid_matrix), parasitoid_species] <- caterpillar_parasitoid_matrix
  caterpillar_parasitoid_matrix <- unified_cat_par_matrix
  
  # -----------------------------------------------------------
  # 2) Replace visible labels with a single space " " 
  # -----------------------------------------------------------
  # Plants (web 1, rows) -> Set to " "
  rownames(plant_caterpillar_matrix) <- rep(" ", nrow(plant_caterpillar_matrix))
  
  # Caterpillars (web 1, columns AND web 2, rows) -> Names are KEPT (NO CHANGE)
  
  # Parasitoids (web 2, columns) -> Set to " "
  colnames(caterpillar_parasitoid_matrix) <- rep(" ", ncol(caterpillar_parasitoid_matrix))
  
  # File paths
  base <- file.path(out_dir, paste0("tritrophic_", safe_name(loc)))
  png_file <- paste0(base, ".png")
  svg_file <- paste0(base, ".svg")
  
  # ---- Save SVG ----
  tryCatch({
    svg(svg_file, width = 12, height = 8)  # inches
    plotweb2(
      web  = plant_caterpillar_matrix,
      web2 = caterpillar_parasitoid_matrix,
      method = NULL,      # ZMĚNA: Vypne řazení první matice (Planty)
      method2 = NULL,     # ZMĚNA: Vypne řazení druhé matice (Parazitoidi)
      empty = FALSE,
      empty2 = FALSE,
      labsize = 0.01,
      ybig = 1,
      y_width = 0.05,
      spacing = 0.05,
      spacing2 = 0.05,
      arrow = "down",
      arrow2 = "up",
      col.interaction  = "grey80",
      col.interaction2 = "grey80",
      low.abun.col = "orange",
      high.abun.col = "darkgreen",
      low.abun.col2  = "orange",
      high.abun.col2 = "darkgreen",
      col.prey  = "darkgreen",
      col.pred  = "orange",
      col.prey2 = "orange",
      col.pred2 = "lightblue",
      lab.space = 1
    )
    dev.off()
  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()
    message("  ERROR while saving SVG for ", loc, " : ", conditionMessage(e))
  })
  
  # ---- Save PNG ----
  tryCatch({
    png(png_file, width = 2400, height = 1600, res = 200)  # pixels @ 200 dpi
    bipartite::plotweb2(
      web  = plant_caterpillar_matrix,
      web2 = caterpillar_parasitoid_matrix,
      method = NULL,      # ZMĚNA: Vypne řazení první matice (Planty)
      method2 = NULL,     # ZMĚNA: Vypne řazení druhé matice (Parazitoidi)
      empty = FALSE,
      empty2 = FALSE,
      labsize = 0.01,
      ybig = 1,
      y_width = 0.05,
      spacing = 0.05,
      spacing2 = 0.05,
      arrow = "down",
      arrow2 = "up",
      col.interaction  = "grey80",
      col.interaction2 = "grey80",
      low.abun.col = "orange",
      high.abun.col = "darkgreen",
      low.abun.col2  = "orange",
      high.abun.col2 = "darkgreen",
      col.prey  = "darkgreen",
      col.pred  = "orange",
      col.prey2 = "orange",
      col.pred2 = "lightblue",
      lab.space = 1
    )
    dev.off()
  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()
    message("  ERROR while saving PNG for ", loc, " : ", conditionMessage(e))
  })
  
  message("  Saved (or attempted): ", basename(svg_file), " & ", basename(png_file))
}


#####################################


## --- KOMBINOVANÝ TRITROFICKÝ WEB + EXPORT CSV + COLORBLIND FRIENDLY BARVY ---
library(tidyverse)
library(bipartite)

out_dir_comb <- "FW_outputs_bipartite_combined"
dir.create(out_dir_comb, recursive = TRUE, showWarnings = FALSE)

# Zda použít četnosti (TRUE) nebo pouze presence/absence (FALSE)
use_freq <- TRUE

# --- 1) Sestavení celkových matic ---
plant_caterpillar <- raw_data %>%
  filter(!is.na(PLANT_sp) & !is.na(CAT_sp)) %>%
  { if (use_freq) group_by(., PLANT_sp, CAT_sp) %>% summarise(freq = n(), .groups = "drop")
    else distinct(., PLANT_sp, CAT_sp) %>% mutate(freq = 1) } %>%
  pivot_wider(names_from = CAT_sp, values_from = freq, values_fill = 0) %>%
  tibble::column_to_rownames(var = "PLANT_sp") %>%
  as.matrix()

caterpillar_parasitoid <- raw_data %>%
  filter(!is.na(CAT_sp) & !is.na(PAR_sp)) %>%
  { if (use_freq) group_by(., CAT_sp, PAR_sp) %>% summarise(freq = n(), .groups = "drop")
    else distinct(., CAT_sp, PAR_sp) %>% mutate(freq = 1) } %>%
  pivot_wider(names_from = PAR_sp, values_from = freq, values_fill = 0) %>%
  tibble::column_to_rownames(var = "CAT_sp") %>%
  as.matrix()

# Kontrola
if (is.null(dim(plant_caterpillar)) || any(dim(plant_caterpillar) == 0) ||
    is.null(dim(caterpillar_parasitoid)) || any(dim(caterpillar_parasitoid) == 0)) {
  stop("Jedna z matic je prázdná. Zkontroluj raw_data (PLANT_sp, CAT_sp, PAR_sp).")
}

# --- 2) Zajistit unikátní interní jména (nutné pro plotweb2) ---
if (is.null(rownames(plant_caterpillar))) rownames(plant_caterpillar) <- paste0("PL_row_", seq_len(nrow(plant_caterpillar)))
if (is.null(colnames(plant_caterpillar))) colnames(plant_caterpillar) <- paste0("PL_col_", seq_len(ncol(plant_caterpillar)))
if (is.null(rownames(caterpillar_parasitoid))) rownames(caterpillar_parasitoid) <- paste0("CA_row_", seq_len(nrow(caterpillar_parasitoid)))
if (is.null(colnames(caterpillar_parasitoid))) colnames(caterpillar_parasitoid) <- paste0("CA_col_", seq_len(ncol(caterpillar_parasitoid)))

rownames(plant_caterpillar)      <- make.unique(rownames(plant_caterpillar))
colnames(plant_caterpillar)      <- make.unique(colnames(plant_caterpillar))
rownames(caterpillar_parasitoid) <- make.unique(rownames(caterpillar_parasitoid))
colnames(caterpillar_parasitoid) <- make.unique(colnames(caterpillar_parasitoid))

# --- 3) Sjednocení střední úrovně (všechny housenky) ---
all_caterpillars <- sort(unique(c(colnames(plant_caterpillar), rownames(caterpillar_parasitoid))))

# Rozšíření plant->cat (sloupce = all_caterpillars)
plant_species <- rownames(plant_caterpillar)
unified_plant_cat <- matrix(0,
                            nrow = length(plant_species),
                            ncol = length(all_caterpillars),
                            dimnames = list(plant_species, all_caterpillars))
if (ncol(plant_caterpillar) > 0) {
  unified_plant_cat[plant_species, colnames(plant_caterpillar)] <- plant_caterpillar
}
plant_caterpillar <- unified_plant_cat

# Rozšíření cat->par (řádky = all_caterpillars)
parasitoid_species <- colnames(caterpillar_parasitoid)
unified_cat_par <- matrix(0,
                          nrow = length(all_caterpillars),
                          ncol = length(parasitoid_species),
                          dimnames = list(all_caterpillars, parasitoid_species))
if (nrow(caterpillar_parasitoid) > 0) {
  unified_cat_par[rownames(caterpillar_parasitoid), parasitoid_species] <- caterpillar_parasitoid
}
caterpillar_parasitoid <- unified_cat_par

# --- 4) Export CSV (pro další analýzy) ---
write.csv(as.data.frame(plant_caterpillar), file = file.path(out_dir_comb, "plant_caterpillar_combined_matrix.csv"), row.names = TRUE)
write.csv(as.data.frame(caterpillar_parasitoid), file = file.path(out_dir_comb, "caterpillar_parasitoid_combined_matrix.csv"), row.names = TRUE)

# --- 5) Estetika: skrytí jmen u rostlin a parazitoidů (volitelné) ---
# POZOR: NEPŘEPISUJ jména střední úrovně (housenky)
rownames(plant_caterpillar) <- rep(" ", nrow(plant_caterpillar))
colnames(caterpillar_parasitoid) <- rep(" ", ncol(caterpillar_parasitoid))

# --- 6) COLORBLIND-FRIENDLY barvy pro interakce a abundance ---
# Okabe–Ito / obecně používané colorblind friendly odstíny:
# zde volíme: plant–caterpillar = oranžovo-červená; caterpillar–parasitoid = tmavě modrá
col_plant_cat_link <- "#D55E00"  # červenooranžová (plant -> caterpillar)
col_cat_par_link   <- "#0072B2"  # tmavě modrá (caterpillar -> parasitoid)

# Low/high barvy (jemný gradient laděný k hlavní barvě)
low_abun_1  <- "#FFDAB3"  # světlá varianta pro plant-cat
high_abun_1 <- "#D55E00"  # tmavší pro plant-cat
low_abun_2  <- "#CFEAF7"  # světlá varianta pro cat-par
high_abun_2 <- "#0072B2"  # tmavší pro cat-par

# --- 7) Vykreslení (SVG + PNG) ---
base_comb <- file.path(out_dir_comb, "tritrophic_combined")
svg_file <- paste0(base_comb, ".svg")
png_file <- paste0(base_comb, ".png")

# Pokud chceš řadit uzly, nastav method/ method2 na "cca" / "null" podle potřeby
method1 <- "cca"
method2 <- "cca"

# SVG
tryCatch({
  svg(svg_file, width = 12, height = 8)
  plotweb2(
    web  = plant_caterpillar,
    web2 = caterpillar_parasitoid,
    method = method1,
    method2 = method2,
    empty = FALSE,
    empty2 = FALSE,
    labsize = 0.01,
    ybig = 1,
    y_width = 0.05,
    spacing = 0.05,
    spacing2 = 0.05,
    arrow = "down",
    arrow2 = "up",
    # DŮLEŽITÉ: col.interaction = barva linků pro web (web = plant->cat)
    #            col.interaction2= barva linků pro web2 (cat->par)
    # barvy pro abundance uzlů (laděné k link barvám)
    # barvy uzlů (prey/pred) - ladíme decentně
    low.abun.col = "orange",
    high.abun.col = "darkgreen",
    low.abun.col2  = "orange",
    high.abun.col2 = "darkgreen",
    col.prey  = "darkgreen",
    col.pred  = "orange",
    col.prey2 = "orange",
    col.pred2 = "lightblue",
    col.interaction = "red",
    bor.col.interaction = "red",
    bor.col.low = "purple",
    lab.space = 1
  )
  dev.off()
  message("SVG uložen: ", svg_file)
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  message("Chyba při ukládání SVG: ", conditionMessage(e))
})

# PNG
tryCatch({
  png(png_file, width = 2400, height = 1600, res = 200)
  bipartite::plotweb2(
    web  = plant_caterpillar,
    web2 = caterpillar_parasitoid,
    method = method1,
    method2 = method2,
    empty = FALSE,
    empty2 = FALSE,
    labsize = 0.01,
    ybig = 1,
    y_width = 0.05,
    spacing = 0.05,
    spacing2 = 0.05,
    arrow = "down",
    arrow2 = "up",
    col.interaction  = col_plant_cat_link,
    col.interaction2 = col_cat_par_link,
    low.abun.col = "orange",
    high.abun.col = "darkgreen",
    low.abun.col2  = "orange",
    high.abun.col2 = "darkgreen",
    col.prey  = "darkgreen",
    col.pred  = "orange",
    col.prey2 = "orange",
    col.pred2 = "lightblue",
    lab.space = 1
  )
  dev.off()
  message("PNG uložen: ", png_file)
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  message("Chyba při ukládání PNG: ", conditionMessage(e))
})

############################
# --- BLOC: EXPORT CSV + SAMOSTATNÉ PLOTY (plotweb) + KOMBINOVANÝ TRITROFICKÝ (plotweb2) ---
library(tidyverse)
library(bipartite)

out_dir_comb <- "FW_outputs_bipartite_combined"
dir.create(out_dir_comb, recursive = TRUE, showWarnings = FALSE)

use_freq <- TRUE  # TRUE = counts, FALSE = presence/absence

# --- 1) Sestavení matice plant->caterpillar a cat->parasitoid (celkově) ---
plant_caterpillar_matrix <- raw_data %>%
  filter(!is.na(PLANT_sp) & !is.na(CAT_sp)) %>%
  { if (use_freq) group_by(., PLANT_sp, CAT_sp) %>% summarise(freq = n(), .groups = "drop")
    else distinct(., PLANT_sp, CAT_sp) %>% mutate(freq = 1) } %>%
  pivot_wider(names_from = CAT_sp, values_from = freq, values_fill = 0) %>%
  tibble::column_to_rownames(var = "PLANT_sp") %>%
  as.matrix()

caterpillar_parasitoid_matrix <- raw_data %>%
  filter(!is.na(CAT_sp) & !is.na(PAR_sp)) %>%
  { if (use_freq) group_by(., CAT_sp, PAR_sp) %>% summarise(freq = n(), .groups = "drop")
    else distinct(., CAT_sp, PAR_sp) %>% mutate(freq = 1) } %>%
  pivot_wider(names_from = PAR_sp, values_from = freq, values_fill = 0) %>%
  tibble::column_to_rownames(var = "CAT_sp") %>%
  as.matrix()

# Kontrola
if (is.null(dim(plant_caterpillar_matrix)) || any(dim(plant_caterpillar_matrix) == 0) ||
    is.null(dim(caterpillar_parasitoid_matrix)) || any(dim(caterpillar_parasitoid_matrix) == 0)) {
  stop("Jedna z matic je prázdná. Zkontroluj raw_data (PLANT_sp, CAT_sp, PAR_sp).")
}

# --- 2) Ujistit se, že dimnames jsou přítomné a unikátní ---
if (is.null(rownames(plant_caterpillar_matrix))) rownames(plant_caterpillar_matrix) <- paste0("PL_row_", seq_len(nrow(plant_caterpillar_matrix)))
if (is.null(colnames(plant_caterpillar_matrix))) colnames(plant_caterpillar_matrix) <- paste0("PL_col_", seq_len(ncol(plant_caterpillar_matrix)))
if (is.null(rownames(caterpillar_parasitoid_matrix))) rownames(caterpillar_parasitoid_matrix) <- paste0("CA_row_", seq_len(nrow(caterpillar_parasitoid_matrix)))
if (is.null(colnames(caterpillar_parasitoid_matrix))) colnames(caterpillar_parasitoid_matrix) <- paste0("CA_col_", seq_len(ncol(caterpillar_parasitoid_matrix)))

rownames(plant_caterpillar_matrix)      <- make.unique(rownames(plant_caterpillar_matrix))
colnames(plant_caterpillar_matrix)      <- make.unique(colnames(plant_caterpillar_matrix))
rownames(caterpillar_parasitoid_matrix) <- make.unique(rownames(caterpillar_parasitoid_matrix))
colnames(caterpillar_parasitoid_matrix) <- make.unique(colnames(caterpillar_parasitoid_matrix))

# --- 3) Sjednotit střední úroveň (všechny housenky musí být v obou maticích) ---
all_caterpillars <- sort(unique(c(colnames(plant_caterpillar_matrix), rownames(caterpillar_parasitoid_matrix))))

# Rozšířit plant->cat (sloupce = all_caterpillars)
plant_species <- rownames(plant_caterpillar_matrix)
unified_plant_cat <- matrix(0,
                            nrow = length(plant_species),
                            ncol = length(all_caterpillars),
                            dimnames = list(plant_species, all_caterpillars))
if (ncol(plant_caterpillar_matrix) > 0) {
  unified_plant_cat[plant_species, colnames(plant_caterpillar_matrix)] <- plant_caterpillar_matrix
}
plant_caterpillar_matrix <- unified_plant_cat

# Rozšířit cat->par (řádky = all_caterpillars)
parasitoid_species <- colnames(caterpillar_parasitoid_matrix)
unified_cat_par <- matrix(0,
                          nrow = length(all_caterpillars),
                          ncol = length(parasitoid_species),
                          dimnames = list(all_caterpillars, parasitoid_species))
if (nrow(caterpillar_parasitoid_matrix) > 0) {
  unified_cat_par[rownames(caterpillar_parasitoid_matrix), parasitoid_species] <- caterpillar_parasitoid_matrix
}
caterpillar_parasitoid_matrix <- unified_cat_par

# --- 4) Export CSV ---
write.csv(as.data.frame(plant_caterpillar_matrix),
          file = file.path(out_dir_comb, "plant_caterpillar_combined_matrix.csv"),
          row.names = TRUE)
write.csv(as.data.frame(caterpillar_parasitoid_matrix),
          file = file.path(out_dir_comb, "caterpillar_parasitoid_combined_matrix.csv"),
          row.names = TRUE)

# --- 5) PLOTTY: použité barvy podle tvého příkladu (hex hodnoty) ---
# Caterpillar - Parasitoid (podle tvého příkladu)
cp_cols <- list(
  col.interaction = "#4682B4",   # link color
  bor.col.interaction = "#4682B4",
  col.high = "#1B9E77",          # parasitoids (higher trophic level)
  col.low  = "#D95F02",          # caterpillars (lower trophic level)
  bor.col.high = "#1B9E77",
  bor.col.low  = "#D95F02"
)

# Plant - Caterpillar (podle tvého příkladu)
pc_cols <- list(
  col.interaction = "#B4464B",   # link color (plant->caterpillar)
  bor.col.interaction = "#B4464B",
  col.high = "#D95F02",          # caterpillars (higher trophic level)
  col.low  = "#7570B3",          # plants (lower trophic level)
  bor.col.high = "#D95F02",
  bor.col.low  = "#7570B3"
)

# Skryjeme názvy rostlin/parazitoidů, pokud chceš (volitelné)
# POZOR: nikdy nepřepisuj jména housenek, ty musí zůstat unikátní v objektech
rownames(plant_caterpillar_matrix) <- rep(" ", nrow(plant_caterpillar_matrix))
colnames(caterpillar_parasitoid_matrix) <- rep(" ", ncol(caterpillar_parasitoid_matrix))

# --- 6) Uložení samostatných plotů (Caterpillar-Parasitoid) ---
svg_cp <- file.path(out_dir_comb, "plot_cp_caterpillar_parasitoid.svg")
png_cp <- file.path(out_dir_comb, "plot_cp_caterpillar_parasitoid.png")

tryCatch({
  svg(svg_cp, width = 10, height = 6)
  plotweb(
    caterpillar_parasitoid_matrix,
    method = "normal",
    empty = FALSE,
    col.interaction = cp_cols$col.interaction,
    bor.col.interaction = cp_cols$bor.col.interaction,
    col.high = cp_cols$col.high,
    col.low  = cp_cols$col.low,
    bor.col.high = cp_cols$bor.col.high,
    bor.col.low  = cp_cols$bor.col.low,
    labsize = 0.01,
    high.lab.dis = 0.05,
    low.lab.dis = 0.05
  )
  title("Caterpillar - Parasitoid Interactions", line = 2.5)
  dev.off()
  message("SVG uložen: ", svg_cp)
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  message("Chyba při ukládání CP SVG: ", conditionMessage(e))
})

tryCatch({
  png(png_cp, width = 2000, height = 1400, res = 200)
  plotweb(
    caterpillar_parasitoid_matrix,
    method = "normal",
    empty = FALSE,
    col.interaction = cp_cols$col.interaction,
    bor.col.interaction = cp_cols$bor.col.interaction,
    col.high = cp_cols$col.high,
    col.low  = cp_cols$col.low,
    bor.col.high = cp_cols$bor.col.high,
    bor.col.low  = cp_cols$bor.col.low,
    labsize = 0.01,
    high.lab.dis = 0.05,
    low.lab.dis = 0.05
  )
  title("Caterpillar - Parasitoid Interactions", line = 2.5)
  dev.off()
  message("PNG uložen: ", png_cp)
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  message("Chyba při ukládání CP PNG: ", conditionMessage(e))
})

# --- 7) Uložení samostatných plotů (Plant-Caterpillar) ---
svg_pc <- file.path(out_dir_comb, "plot_pc_plant_caterpillar.svg")
png_pc <- file.path(out_dir_comb, "plot_pc_plant_caterpillar.png")

tryCatch({
  svg(svg_pc, width = 10, height = 6)
  plotweb(
    plant_caterpillar_matrix,
    method = "normal",
    empty = FALSE,
    col.interaction = pc_cols$col.interaction,
    bor.col.interaction = pc_cols$bor.col.interaction,
    col.high = pc_cols$col.high,
    col.low  = pc_cols$col.low,
    bor.col.high = pc_cols$bor.col.high,
    bor.col.low  = pc_cols$bor.col.low,
    labsize = 0.01,
    high.lab.dis = 0.05,
    low.lab.dis = 0.05
  )
  title("Plant - Caterpillar Interactions", line = 2.5)
  dev.off()
  message("SVG uložen: ", svg_pc)
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  message("Chyba při ukládání PC SVG: ", conditionMessage(e))
})

tryCatch({
  png(png_pc, width = 2000, height = 1400, res = 200)
  plotweb(
    plant_caterpillar_matrix,
    method = "normal",
    empty = FALSE,
    col.interaction = pc_cols$col.interaction,
    bor.col.interaction = pc_cols$bor.col.interaction,
    col.high = pc_cols$col.high,
    col.low  = pc_cols$col.low,
    bor.col.high = pc_cols$bor.col.high,
    bor.col.low  = pc_cols$bor.col.low,
    labsize = 0.01,
    high.lab.dis = 0.05,
    low.lab.dis = 0.05
  )
  title("Plant - Caterpillar Interactions", line = 2.5)
  dev.off()
  message("PNG uložen: ", png_pc)
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  message("Chyba při ukládání PC PNG: ", conditionMessage(e))
})

# --- 8) KOMBINOVANÝ TRITROFICKÝ (plotweb2) - barvy sladěné s individuálními ploty ---
# Mapuji linky plant->cat na barvu pc_cols$col.interaction,
# mapuji linky cat->par na barvu cp_cols$col.interaction.
# Pro barvy uzlů použiji podobné odstíny, aby to bylo konzistentní a colorblind-friendly.

col_link_plant_cat <- pc_cols$col.interaction
col_link_cat_par   <- cp_cols$col.interaction

# Uzly (prey/pred) barvy laděné k tvému výběru
# plants (low of plant-cat) -> použiju pc_cols$col.low ("#7570B3")
# caterpillars -> použiju pc_cols$col.high / cp_cols$col.low ("#D95F02")
# parasitoids -> použiju cp_cols$col.high ("#1B9E77")
col_prey_1  <- pc_cols$col.low      # plants (low level, left bar)
col_pred_1  <- pc_cols$col.high     # caterpillars (right bar of web1)
col_prey_2  <- cp_cols$col.low      # caterpillars (left bar of web2)
col_pred_2  <- cp_cols$col.high     # parasitoids (right bar of web2)

svg_trit <- file.path(out_dir_comb, "tritrophic_combined_plotweb2.svg")
png_trit <- file.path(out_dir_comb, "tritrophic_combined_plotweb2.png")

tryCatch({
  svg(svg_trit, width = 12, height = 8)
  plotweb2(
    web  = plant_caterpillar_matrix,
    web2 = caterpillar_parasitoid_matrix,
    method = "cca",      # nebo NULL pokud nechceš řazení
    method2 = "cca",
    empty = FALSE,
    empty2 = FALSE,
    labsize = 0.01,
    ybig = 1,
    y_width = 0.05,
    spacing = 0.05,
    spacing2 = 0.05,
    arrow = "down",
    arrow2 = "up",
    col.interaction  = col_link_plant_cat,
    col.interaction2 = col_link_cat_par,
    low.abun.col = "#FFE6E6",   # světlá varianta pro plant->cat abundance (laděná k #B4464B)
    high.abun.col = pc_cols$col.high,
    low.abun.col2 = "#E6F3FF",  # světlá varianta pro cat->par abundance (laděná k #4682B4)
    high.abun.col2 = cp_cols$col.high,
    col.prey  = col_prey_1,
    col.pred  = col_pred_1,
    col.prey2 = col_prey_2,
    col.pred2 = col_pred_2,
    lab.space = 1
  )
  dev.off()
  message("SVG uložen: ", svg_trit)
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  message("Chyba při ukládání tritrophic SVG: ", conditionMessage(e))
})

tryCatch({
  png(png_trit, width = 2400, height = 1600, res = 200)
  bipartite::plotweb2(
    web  = plant_caterpillar_matrix,
    web2 = caterpillar_parasitoid_matrix,
    method = "cca",
    method2 = "cca",
    empty = FALSE,
    empty2 = FALSE,
    labsize = 0.01,
    ybig = 1,
    y_width = 0.05,
    spacing = 0.05,
    spacing2 = 0.05,
    arrow = "down",
    arrow2 = "up",
    col.interaction  = col_link_plant_cat,
    col.interaction2 = col_link_cat_par,
    low.abun.col = "#FFE6E6",
    high.abun.col = pc_cols$col.high,
    low.abun.col2 = "#E6F3FF",
    high.abun.col2 = cp_cols$col.high,
    col.prey  = col_prey_1,
    col.pred  = col_pred_1,
    col.prey2 = col_prey_2,
    col.pred2 = col_pred_2,
    lab.space = 1
  )
  dev.off()
  message("PNG uložen: ", png_trit)
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  message("Chyba při ukládání tritrophic PNG: ", conditionMessage(e))
})

# --- KONEC bloku ---
message("Hotovo — CSV + jednotlivé ploty (PC, CP) + kombinovaný tritrofický plot uloženy do: ", normalizePath(out_dir_comb))

