# ============================================================
# 10_metacommunity_model.R   >>> ANALYSIS ONLY <<<
# Study: Approaching the limit of diversity: spatial turnover
#        amplifies with trophic level (Libra et al.)
#
# Analysis 10: post-processing of the spatially explicit metacommunity
# simulation (C++ output). Computes Sorensen beta diversity across the
# colonisation-rate parameter grid, for resources ("hosts") and
# consumers ("paras"), both spatially and temporally.
#
#   SPATIAL : beta.pair() Sorensen across all 100 sites at t = 100.
#             Empty species columns removed; NaN pairs excluded (na.rm).
#   TEMPORAL: beta.pair() Sorensen between t = 50 and t = 100, per site
#             (100 pairs per file), then averaged.
#
# Each file = one colonisation-rate combination (parsed from filename).
# Outputs: results tables (csv) + output/rds/ for the plotting script.
# Figures are built in 10_metacommunity_model_plots.R.
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(betapart)    # beta.pair()
library(here)        # here() for file paths

# ------------------------------------------------------------
# 1. Paths + load simulation files
# ------------------------------------------------------------
data_dir <- here("Data_model")   # adjust if the model output lives elsewhere
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)

files <- list.files(data_dir, pattern = "LocalData.*N0\\.txt$", full.names = TRUE)
cat("Files found:", length(files), "\n")

# ------------------------------------------------------------
# 2. Parse colonisation rates from filename
# ------------------------------------------------------------
parse_colrates <- function(fname) {
  code <- sub(".*LocalData([0-9]+)N0\\.txt$", "\\1", fname)
  code <- formatC(as.integer(code), width = 4, flag = "0")
  x <- as.integer(substr(code, 1, 2))
  y <- as.integer(substr(code, 3, 4))
  list(host_col = 0.2 + (x - 1) * 0.05,
       para_col = 0.2 + (y - 1) * 0.05)
}

# ------------------------------------------------------------
# 3. Spatial beta diversity (t = 100)
# ------------------------------------------------------------
results <- data.frame(
  file = character(), host_col_rate = numeric(), para_col_rate = numeric(),
  host_SOR = numeric(), para_SOR = numeric(), stringsAsFactors = FALSE
)
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (i in seq_along(files)) {
  f <- files[i]
  tryCatch({
    dat <- read.table(f, header = FALSE, sep = "", fill = TRUE)
    dat_t100 <- dat[dat[, 1] == 100, ]
    if (nrow(dat_t100) != 100) { cat("Warning:", basename(f), "\n"); next }
    
    rates <- parse_colrates(basename(f))
    
    hosts <- dat_t100[, 4:103]
    hosts <- hosts[, colSums(hosts) > 0]
    paras <- dat_t100[, 104:203]
    paras <- paras[, colSums(paras) > 0]
    
    host_sor <- tryCatch(
      if (ncol(hosts) >= 1)
        mean(as.dist(beta.pair(hosts, index.family = "sorensen")$beta.sor), na.rm = TRUE) else NA,
      error = function(e) NA)
    para_sor <- tryCatch(
      if (ncol(paras) >= 1)
        mean(as.dist(beta.pair(paras, index.family = "sorensen")$beta.sor), na.rm = TRUE) else NA,
      error = function(e) NA)
    
    results <- rbind(results, data.frame(
      file = basename(f), host_col_rate = rates$host_col, para_col_rate = rates$para_col,
      host_SOR = host_sor, para_SOR = para_sor, stringsAsFactors = FALSE))
  }, error = function(e) cat("Error:", basename(f), "-", conditionMessage(e), "\n"))
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\nSPATIAL done. Rows:", nrow(results), "\n")
cat("NA in host_SOR:", sum(is.na(results$host_SOR)), "\n")
cat("NA in para_SOR:", sum(is.na(results$para_SOR)), "\n")

# ------------------------------------------------------------
# 4. Temporal beta diversity (t = 50 vs t = 100)
# ------------------------------------------------------------
results_temp <- data.frame(
  file = character(), host_col_rate = numeric(), para_col_rate = numeric(),
  host_SOR_temp = numeric(), para_SOR_temp = numeric(), stringsAsFactors = FALSE
)
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (i in seq_along(files)) {
  f <- files[i]
  tryCatch({
    dat <- read.table(f, header = FALSE, sep = "", fill = TRUE)
    dat_t50  <- dat[dat[, 1] == 50, ]
    dat_t100 <- dat[dat[, 1] == 100, ]
    if (nrow(dat_t50) != 100 || nrow(dat_t100) != 100) { cat("Warning:", basename(f), "\n"); next }
    
    rates <- parse_colrates(basename(f))
    
    hosts_t50  <- dat_t50[, 4:103];  hosts_t100 <- dat_t100[, 4:103]
    paras_t50  <- dat_t50[, 104:203]; paras_t100 <- dat_t100[, 104:203]
    
    host_betas <- sapply(1:100, function(j) tryCatch({
      pm <- rbind(hosts_t50[j, ], hosts_t100[j, ]); pm <- pm[, colSums(pm) > 0]
      if (ncol(pm) < 1) return(NA)
      as.numeric(beta.pair(pm, index.family = "sorensen")$beta.sor)
    }, error = function(e) NA))
    
    para_betas <- sapply(1:100, function(j) tryCatch({
      pm <- rbind(paras_t50[j, ], paras_t100[j, ]); pm <- pm[, colSums(pm) > 0]
      if (ncol(pm) < 1) return(NA)
      as.numeric(beta.pair(pm, index.family = "sorensen")$beta.sor)
    }, error = function(e) NA))
    
    results_temp <- rbind(results_temp, data.frame(
      file = basename(f), host_col_rate = rates$host_col, para_col_rate = rates$para_col,
      host_SOR_temp = mean(host_betas, na.rm = TRUE),
      para_SOR_temp = mean(para_betas, na.rm = TRUE), stringsAsFactors = FALSE))
  }, error = function(e) cat("Error:", basename(f), "-", conditionMessage(e), "\n"))
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\nTEMPORAL done. Rows:", nrow(results_temp), "\n")
cat("NA in host_SOR_temp:", sum(is.na(results_temp$host_SOR_temp)), "\n")
cat("NA in para_SOR_temp:", sum(is.na(results_temp$para_SOR_temp)), "\n")

# ------------------------------------------------------------
# 5. Difference (consumer - resource) + save
# ------------------------------------------------------------
results$diff_SOR           <- results$para_SOR      - results$host_SOR
results_temp$diff_SOR_temp <- results_temp$para_SOR_temp - results_temp$host_SOR_temp

saveRDS(results,      "output/rds/model_beta_spatial.rds")
saveRDS(results_temp, "output/rds/model_beta_temporal.rds")
write.csv(results,      "output/rds/model_beta_spatial_t100.csv",      row.names = FALSE)
write.csv(results_temp, "output/rds/model_beta_temporal_t50_t100.csv", row.names = FALSE)
