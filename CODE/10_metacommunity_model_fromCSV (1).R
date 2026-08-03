# ============================================================
# 10_metacommunity_model_fromCSV.R   >>> TEST / ALTERNATIVE INPUT <<<
# Same analysis as 10_metacommunity_model.R, but reads the consolidated
# LocalData_all.csv.gz (from 10a_consolidate_localdata.R) instead of the
# many raw .txt files.
#
# Kept as a SEPARATE script for testing: it writes to *_fromCSV output
# names (so it never overwrites the originals) and, if the original
# results exist, checks that the two paths agree.
# ============================================================

library(data.table)
library(betapart)

# ------------------------------------------------------------
# 1. Load consolidated CSV
# ------------------------------------------------------------
csv_file <- "LocalData_all.csv.gz"   # <- output of 10a (adjust path if needed)
dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)

big <- fread(csv_file)

host_cols <- sprintf("host_%03d", 1:100)
para_cols <- sprintf("para_%03d", 1:100)
if (!all(c("file", "time", "site", host_cols, para_cols) %in% names(big))) {
  stop("CSV is missing expected columns (file/time/site/host_*/para_*). ",
       "Was it written by 10a with the 203-column layout?")
}
setkey(big, file, time, site)   # keyed subsets are fast and site-ordered

# ------------------------------------------------------------
# 2. Parse colonisation rates from the file name (same as original)
# ------------------------------------------------------------
parse_colrates <- function(fname) {
  code <- sub(".*LocalData([0-9]+)N0\\.txt$", "\\1", fname)
  code <- formatC(as.integer(code), width = 4, flag = "0")
  x <- as.integer(substr(code, 1, 2)); y <- as.integer(substr(code, 3, 4))
  list(host_col = 0.2 + (x - 1) * 0.05, para_col = 0.2 + (y - 1) * 0.05)
}

files <- unique(big$file)
cat("Files in CSV:", length(files), "\n")

sor_mean <- function(mat) {
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (ncol(mat) < 1) return(NA_real_)
  tryCatch(mean(as.dist(beta.pair(mat, index.family = "sorensen")$beta.sor), na.rm = TRUE),
           error = function(e) NA_real_)
}

# ------------------------------------------------------------
# 3. Spatial beta diversity (t = 100)
# ------------------------------------------------------------
results <- data.frame(file = character(), host_col_rate = numeric(), para_col_rate = numeric(),
                      host_SOR = numeric(), para_SOR = numeric(), stringsAsFactors = FALSE)
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (i in seq_along(files)) {
  f <- files[i]
  d100 <- big[.(f, 100)]
  if (nrow(d100) != 100) { cat("Warning:", f, "\n"); setTxtProgressBar(pb, i); next }
  rates <- parse_colrates(f)
  results <- rbind(results, data.frame(
    file = f, host_col_rate = rates$host_col, para_col_rate = rates$para_col,
    host_SOR = sor_mean(as.matrix(d100[, ..host_cols])),
    para_SOR = sor_mean(as.matrix(d100[, ..para_cols])), stringsAsFactors = FALSE))
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\nSPATIAL done. Rows:", nrow(results), "\n")

# ------------------------------------------------------------
# 4. Temporal beta diversity (t = 50 vs t = 100)
# ------------------------------------------------------------
pair_sor <- function(m50, m100) {
  sapply(seq_len(nrow(m50)), function(j) {
    pm <- rbind(m50[j, ], m100[j, ]); pm <- pm[, colSums(pm) > 0, drop = FALSE]
    if (ncol(pm) < 2) return(NA_real_)   # match original: single-species pair -> NA
    tryCatch(as.numeric(beta.pair(pm, index.family = "sorensen")$beta.sor),
             error = function(e) NA_real_)
  })
}

results_temp <- data.frame(file = character(), host_col_rate = numeric(), para_col_rate = numeric(),
                           host_SOR_temp = numeric(), para_SOR_temp = numeric(), stringsAsFactors = FALSE)
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (i in seq_along(files)) {
  f <- files[i]
  d50  <- big[.(f, 50)]
  d100 <- big[.(f, 100)]
  if (nrow(d50) != 100 || nrow(d100) != 100) { cat("Warning:", f, "\n"); setTxtProgressBar(pb, i); next }
  rates <- parse_colrates(f)
  results_temp <- rbind(results_temp, data.frame(
    file = f, host_col_rate = rates$host_col, para_col_rate = rates$para_col,
    host_SOR_temp = mean(pair_sor(as.matrix(d50[, ..host_cols]),  as.matrix(d100[, ..host_cols])),  na.rm = TRUE),
    para_SOR_temp = mean(pair_sor(as.matrix(d50[, ..para_cols]),  as.matrix(d100[, ..para_cols])),  na.rm = TRUE),
    stringsAsFactors = FALSE))
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\nTEMPORAL done. Rows:", nrow(results_temp), "\n")

# ------------------------------------------------------------
# 5. Difference + save (to *_fromCSV names, so originals are untouched)
# ------------------------------------------------------------
results$diff_SOR           <- results$para_SOR      - results$host_SOR
results_temp$diff_SOR_temp <- results_temp$para_SOR_temp - results_temp$host_SOR_temp

saveRDS(results,      "output/rds/model_beta_spatial_fromCSV.rds")
saveRDS(results_temp, "output/rds/model_beta_temporal_fromCSV.rds")
write.csv(results,      "output/rds/model_beta_spatial_t100_fromCSV.csv",      row.names = FALSE)
write.csv(results_temp, "output/rds/model_beta_temporal_t50_t100_fromCSV.csv", row.names = FALSE)

# ------------------------------------------------------------
# 6. Validation: compare against the original file-based results (if present)
# ------------------------------------------------------------
compare_to_original <- function(new, orig_path, cols) {
  if (!file.exists(orig_path)) { cat("  (no original at", orig_path, "- skipping)\n"); return(invisible()) }
  orig <- readRDS(orig_path)
  m <- merge(new, orig, by = "file", suffixes = c(".csv", ".orig"))
  for (cc in cols) {
    d <- max(abs(m[[paste0(cc, ".csv")]] - m[[paste0(cc, ".orig")]]), na.rm = TRUE)
    cat(sprintf("  max |diff| in %-14s = %.3g\n", cc, d))
  }
}
cat("\n--- Validation vs original file-based results ---\n")
cat("Spatial:\n");  compare_to_original(results,      "output/rds/model_beta_spatial.rds",  c("host_SOR", "para_SOR"))
cat("Temporal:\n"); compare_to_original(results_temp, "output/rds/model_beta_temporal.rds", c("host_SOR_temp", "para_SOR_temp"))
cat("\nIf the max |diff| values are ~0 (e.g. < 1e-10), the CSV path reproduces the originals.\n")
