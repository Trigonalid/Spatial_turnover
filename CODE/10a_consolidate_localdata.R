# ============================================================
# 10a_consolidate_localdata.R
# Consolidates the many raw metacommunity-model output files
# (Data_model/LocalData<set>N0.txt) into a single compressed CSV,
# DATA/LocalData_all.csv.gz, for fast downstream use by
# 10_metacommunity_model_fromCSV.R.
#
# Only the N0 replicate is used (matching 10_metacommunity_model.R's
# file pattern), and only the two recorded time points (50 and 100),
# matching the raw files as produced by the model (see
# https://github.com/AceRNorth/HostParasitoidModel).
#
# Raw file format (one row per patch per recorded time point, no header):
#   col 1       = T (recorded time point; 50 or 100)
#   col 2, 3    = x, y patch coordinates (not carried over)
#   col 4:103   = Hpres_1..Hpres_100 (host presence/absence, 0/1)
#   col 104:203 = Ppres_1..Ppres_100 (parasitoid presence/absence, 0/1)
#
# Output columns: file, time, site, host_001..host_100, para_001..para_100
# ============================================================

# ------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------
library(data.table)
library(here)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
data_dir <- here("Data_model")
out_file <- here("DATA/LocalData_all.csv.gz")
dir.create(here("DATA"), showWarnings = FALSE, recursive = TRUE)

files <- list.files(data_dir, pattern = "LocalData.*N0\\.txt$", full.names = TRUE)
cat("Files found:", length(files), "\n")
if (length(files) == 0) {
  stop("No LocalData*N0.txt files found in ", data_dir,
       " — check that the raw model output has been placed there.")
}

# ------------------------------------------------------------
# 2. Column names for the 203-column raw layout
# ------------------------------------------------------------
host_cols <- sprintf("host_%03d", 1:100)
para_cols <- sprintf("para_%03d", 1:100)
raw_names <- c("time", "x", "y", host_cols, para_cols)

# ------------------------------------------------------------
# 3. Read, tag with file + site, keep only t = 50 and t = 100
# ------------------------------------------------------------
read_one <- function(f) {
  dat <- fread(f, header = FALSE, sep = "", fill = TRUE)

  if (ncol(dat) != length(raw_names)) {
    warning("Skipping ", basename(f), ": expected ", length(raw_names),
            " columns, found ", ncol(dat))
    return(NULL)
  }
  setnames(dat, raw_names)

  dat <- dat[time %in% c(50, 100)]
  if (nrow(dat) == 0) {
    warning("Skipping ", basename(f), ": no rows at time = 50 or 100")
    return(NULL)
  }

  # site index: 1..100 within each time block, in the row order as written
  dat[, site := seq_len(.N), by = time]
  dat[, file := basename(f)]
  dat[, c("x", "y") := NULL]   # coordinates not needed downstream

  setcolorder(dat, c("file", "time", "site", host_cols, para_cols))
  dat
}

pb  <- txtProgressBar(min = 0, max = length(files), style = 3)
out <- vector("list", length(files))
for (i in seq_along(files)) {
  out[[i]] <- read_one(files[i])
  setTxtProgressBar(pb, i)
}
close(pb)

big <- rbindlist(out, use.names = TRUE)
cat("\nConsolidated rows:", nrow(big), " (expected ", length(files) * 200, " if every file has",
    " 100 sites x 2 time points and none were skipped)\n", sep = "")

# ------------------------------------------------------------
# 4. Write compressed CSV
# ------------------------------------------------------------
fwrite(big, out_file, compress = "gzip")
cat("Written to", out_file, "\n")
