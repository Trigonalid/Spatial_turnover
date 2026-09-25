# Spatial turnover amplifies with trophic level in hyperdiverse food webs

Analysis code for **Libra et al.** — the R pipeline reproducing all analyses, tables, and figures in the manuscript and supplement. The C++ metacommunity simulation used to generate the Fig. 3 model data lives in its own repository, linked below, with its own README covering its system requirements, installation, and usage.

> **Metacommunity model repository:** <https://github.com/AceRNorth/HostParasitoidModel>

> **Citation / DOI:** *add on acceptance* **Code repository:** <https://github.com/Trigonalid/Spatial_turnover> **License:** MIT — see "License" section below

------------------------------------------------------------------------

## Overview

The study quantifies how **spatial (and temporal) species turnover changes across trophic levels** (plants → caterpillars → parasitoids), using a plant–caterpillar–parasitoid rearing dataset from eight lowland rainforest localities in Papua New Guinea, plus climate data and a spatially explicit metacommunity simulation.

> **Note on the data in this repository:** `DATA/MASTER.xlsx`, as shipped both here on [github.com/Trigonalid/Spatial_turnover](https://github.com/Trigonalid/Spatial_turnover) and used by every R script in this repository, is a **reduced, demo-sized subset (\~2,700 rearing records)** of the full study dataset (\~27,536 records, as reported in the manuscript, Table S1). Running the pipeline as shipped reproduces the pipeline's logic on the subsampled data, not the manuscript's exact reported values. **The full, unreduced dataset — and the raw metacommunity model output (`Data_model/*.txt`, `DATA/LocalData_all.csv.gz`) — are deposited on figshare: <https://doi.org/10.6084/m9.figshare.33449206>.** Substitute the full `MASTER.xlsx` in place of the demo file, keeping the same file name and column format, to reproduce the manuscript's results exactly.

Each analysis is split into an **analysis script** (`NN_name.R`, writes results to `output/rds/`) and a **plotting script** (`NN_name_plots.R`, writes figures to `output/fig/`). Run the analysis script before its plotting script.

## Repository structure

```         
DATA/                     # inputs: MASTER.xlsx (demo-sized), monthly_data.csv (pre-downloaded), LocalData_all.csv.gz
Data_model/               # raw model output (LocalData*N0.txt), consolidated by 10a — see below
output/rds/                # intermediate results (.rds) + result tables (.csv)
output/fig/                # figures (.pdf/.png/.svg)
NN_*.R                     # analysis + plotting scripts
```

`output/` is created automatically. The C++ metacommunity simulation itself (source code, `hosts.csv`/`parasitoids.csv`, compilation) is **not** part of this repository — see <https://github.com/AceRNorth/HostParasitoidModel>. This repository only consumes its output files (`Data_model/LocalData<set>N<rep>.txt`, or the consolidated `DATA/LocalData_all.csv.gz`) in `10a_consolidate_localdata.R` → `10_metacommunity_model_fromCSV.R` → `10_metacommunity_model_plots.R`.

------------------------------------------------------------------------

## 1. System requirements

**Operating systems tested:** - macOS Tahoe 26 (Apple silicon / ARM64)

**Software dependencies:** - **R** ≥ 4.1 (developed and tested on **R 4.3.1**) - R packages (see full list under "Installation guide" below); consider running `sessionInfo()` after installation and saving the output to `sessionInfo.txt` for the record

**Non-standard hardware:** - The R analysis pipeline (all `NN_*.R` scripts) runs on a standard desktop/laptop; no special hardware is required. - The C++ metacommunity simulation used to generate the Fig. 3 model data is a **separate repository** with its own system requirements (it needs more than a standard desktop for the full parameter grid) — see <https://github.com/AceRNorth/HostParasitoidModel> for details. This repository only needs its output files, not the simulation itself, to reproduce Fig. 3.

------------------------------------------------------------------------

## 2. Installation guide

**Instructions:**

1.  Install R (≥ 4.1) from [CRAN](https://cran.r-project.org/).
2.  Install required packages:

``` r
install.packages(c(
  "dplyr","tidyr","purrr","tibble","readxl","here","data.table","R.utils","forcats",
  "vegan","betapart","CommEcol","bipartite","lme4","broom.mixed","terra",
  "ggplot2","ggpubr","ggsignif","patchwork","ggrepel","viridis","maps","scales","cowplot","svglite",
  "flextable","officer","writexl","gt","rstatix"
))
```

Optional: `QBMS` (downloading raw TerraClimate rasters), `iNEXT.beta3D` (exploratory `01b` only)

**Typical install time on a normal desktop computer:** approximately 10–20 minutes (mostly R package compilation/download time, depending on internet speed).

The C++ metacommunity simulation is not compiled or run from this repository — see <https://github.com/AceRNorth/HostParasitoidModel> for its own installation guide. The simulation run for long period, better to use cluster computer.

------------------------------------------------------------------------

## 3. Instructions for use

`DATA/MASTER.xlsx` as shipped is the demo-sized subset described above. To confirm the pipeline runs correctly, run e.g.:

``` r
source("01_spatial_turnover_species.R")
source("01_spatial_turnover_species_plots.R")
```

This produces `output/rds/Table_S2_dissimilarity_summary_new.docx`/`.csv` (spatial dissimilarity by guild) and `output/fig/Figure_1C_BrayCurtis.*`/`Figure_S4_beta_diversity.*` (showing higher spatial turnover in parasitoids than caterpillars), confirming the analysis and plotting code execute correctly end to end — though, on the demo-sized data, these will not match the manuscript's exact reported values (see note above). Most of the runtime is the 999-iteration subsampling loop in step 7 of `01_spatial_turnover_species.R`; time it yourself with `system.time(source("01_spatial_turnover_species.R"))` — *add the measured elapsed time here*.

**To reproduce the manuscript's results:** replace the demo-sized `DATA/MASTER.xlsx` with the full, unreduced dataset from figshare (<https://doi.org/10.6084/m9.figshare.33449206>), keeping `monthly_data.csv` and `LocalData_all.csv.gz` in `DATA/` as described under "Data" below, then run each analysis script followed by its plotting script, in the order given in the table below.

**To run the pipeline on your own data:** replace `DATA/MASTER.xlsx` with a dataset in the same column format (one row per reared caterpillar, with `locality`, `PLANT_species_code`, `CAT_scientific_name`, and an optional `PAR_species_code`), and re-run the scripts in order. No code changes are required as long as the column names and the `locality` values match.

### (Optional) Reproduction instructions — mapping to manuscript figures/tables

| Script (analysis → plots) | Manuscript output |
|------------------------------------|------------------------------------|
| `00_general_summary.R` | Table S1 |
| `01a_study_area_map_plot.R` | Fig. 1a |
| `01b_metaweb_code.R` | Fig. 1b |
| `01_spatial_turnover_species.R` → `_plots.R` | Fig. 1c, Fig. S4, Table S2 |
| `02_rarity_parasitoids.R` → `_plots.R` | Fig. S5, Table S3, Table S10 |
| `04_interaction_turnover.R` → `_plots.R` | Fig. 1d, Fig. 2, Fig. S6, Table S4, Table S5, Table S9 |
| `06_parasitism_rates.R` | Fig. S9 |
| `06_glmm_parasitism.R` → `_plots.R` | Fig. S10, Table S11 |
| `07_network_specialization.R` → `_plots.R` | Fig. S7, Table S6 |
| `08_abiotic_test.R` (or `run_08.R`) | Fig. S1, Fig. S2, Fig. S3, Table S7 |
| `09_temporal_turnover_ohu.R` → `_plots.R` | Fig. S11 |
| `11_enemy_free_space.R` | Fig. S8 |
| `10a_consolidate_localdata.R` → `10_metacommunity_model_fromCSV.R` → `10_metacommunity_model_plots.R` | Fig. 3b, c, d, f, g, Table S8 |

Set the working directory to the repository root, ensure `DATA/` holds the inputs, then run each analysis script followed by its plotting script.

## Data

| File | Description |
|-----------------------|-------------------------------------------------|
| `DATA/MASTER.xlsx` | One row per reared caterpillar: locality, plant, caterpillar, optional parasitoid, guild. **This repository ships a demo-sized subset (\~2,700 records); the full dataset (\~27,536 records, as reported in the manuscript) is deposited at <https://doi.org/10.6084/m9.figshare.33449206>.** Input for analyses 1–9 and 11. |
| `DATA/monthly_data.csv` | TerraClimate `tmin`/`tmax`/`ppt` at the 8 localities (monthly, 1970–2006). **Provided pre-downloaded and pre-processed into this CSV** — downloading and extracting the raw TerraClimate `.nc` rasters (via `QBMS`, in `08a`) takes a very long time, so this repository ships the already-extracted result rather than requiring reviewers to re-download it. |
| `DATA/LocalData_all.csv.gz` | All metacommunity-model output consolidated into one compressed CSV (times 50 and 100); produced by `10a`. |

`MASTER.xlsx` is loaded with `read_excel(..., guess_max = 1048576)` and the parasitoid column coerced with `as.character()` — the column is sparse and otherwise gets silently dropped. **The raw `Data_model/*.txt` files (one per colonisation-rate combination, \~2,209 files) are the main metacommunity-model output and are provided in the figshare data deposit** alongside the full `MASTER.xlsx` (<https://doi.org/10.6084/m9.figshare.33449206>) — too numerous for GitHub

## Conventions

- Analysis and plotting are separate scripts (`08_abiotic_gradient_combined.R` is the one exception).
- Subsampling (caterpillars in `01`, interactions in `04`) uses **999 iterations**, `set.seed(1234)`.
- Tables exported as Word (`.docx`) + CSV; all code and object names in English.

## Notes on methods

- Bray–Curtis dissimilarity via **`vegan`**; Sørensen/Chao–Sørensen via `betapart`/`CommEcol`.
- Mixed-model estimates use **Wald 95% confidence intervals** (`broom.mixed`); climate PCA uses base `prcomp()`.
- Climate variables at **full monthly resolution** (month × year, 1998 excluded).
- Interaction turnover via `bipartite::betalinkr_multi()` with `partitioning = "commondenom"`.

------------------------------------------------------------------------

## Metacommunity Simulation — Host–Parasitoid Food Web

The stochastic patch-occupancy simulation (Gillespie algorithm) used to generate the Fig. 3 model data is maintained as its own repository, with its own README covering the model description, system requirements, compilation, parameters, and usage:

> [**https://github.com/AceRNorth/HostParasitoidModel**](https://github.com/AceRNorth/HostParasitoidModel){.uri}

This repository does not run the simulation — it only reads its output files, for post-processing and figure generation (turnover calculations, heatmaps, and scatter plots for Fig. 3), in `10a_consolidate_localdata.R` → `10_metacommunity_model_fromCSV.R` → `10_metacommunity_model_plots.R`. The linked model repository covers the simulation itself but not this downstream analysis or plotting, which is why those steps live here instead. For everything about the simulation (compiling the model, choosing `delta`/`gamma`/`pcol`/`numvil`, running the parameter grid), see the linked repository.

**Why the raw output files are provided rather than re-run:** the full parameter-grid simulation (2,209 colonisation-rate combinations) takes substantially more than a normal desktop can complete in a practical time frame (see "Non-standard hardware" above). To let reviewers reproduce Fig. 3 without re-running the simulation themselves, **the raw per-combination model output (`Data_model/LocalData<set>N0.txt`, one file per combination) is the main data artifact provided on figshare** (<https://doi.org/10.6084/m9.figshare.33449206>; see also "Data" below) — `10a_consolidate_localdata.R` only needs to read and consolidate these already-computed files, which takes seconds to minutes, not days.

### Output file format expected by the R scripts

One log file per replicate, named `LocalData<set>N<replicate>.txt`, with one row per patch per recorded time point:

```         
T   x   y   Hpres_0   Hpres_1   ...   Ppres_0   Ppres_1   ...
```

`T` is the recorded time point, `x`/`y` are the patch's spatial coordinates (not used by the R scripts), and `Hpres_i`/`Ppres_i` are 0/1 presence flags for host/parasitoid species `i` (100 of each, matching `nh = np = 100` in Table S8). `10_metacommunity_model.R` reads columns 4–103 as host presence and 104–203 as parasitoid presence, and expects two time points per file (50 and 100, matching the pre-processed data in `DATA/LocalData_all.csv.gz`).

## License

This code is released under the **MIT License** — permissive, allowing free use, modification, and redistribution (including commercial use), provided the copyright notice and license text are kept. See `LICENSE` in the repository root for the full text.

**Citation:** if you use this code (in whole or in part), please cite the associated manuscript: *add full citation on acceptance*.

*(Note: an open-source license cannot legally enforce citation of a paper — that is a scholarly norm, not a license term, and OSI-approved licenses do not carry citation-mandatory clauses. MIT + an explicit citation request, as above, is the standard way academic software combines free reuse with a citation expectation.)*

## Code availability statement (for the manuscript)

For the Nature Code and Software Submission Checklist item "Your manuscript should include a complete, detailed description of the code's functionality" — check **"Methods section"**. All standard statistical methods (Bray–Curtis/Sørensen/Chao–Sørensen dissimilarity, GLMM, network specialization H2', PCA) are described there by name, package, and parameters, which is sufficient since these are established methods from cited R packages rather than novel algorithms. The one piece of custom software — the C++ metacommunity simulation — is described at pseudocode/equation level in the Methods subsections "Metacommunity model", "Resources", "Consumers", and "Simulations" (colonisation and extinction rates for resource and consumer species, and the simulated parameter range). Its own repository (<https://github.com/AceRNorth/HostParasitoidModel>) additionally documents the implementation.

## Contact

[martin.libra.cz\@gmail.dom](mailto:martin.libra.cz@gmail.dom){.email}
