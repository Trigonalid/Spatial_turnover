# Spatial turnover amplifies with trophic level in hyperdiverse food webs

Analysis code for **Libra et al.** — the R pipeline reproducing all analyses, tables, and figures
in the manuscript and supplement, plus the C++ metacommunity simulation that generates the model data.

> **Citation / DOI:** _add on acceptance_

---

## Overview

The study quantifies how **spatial (and temporal) species turnover changes across trophic levels**
(plants → caterpillars → parasitoids), using a plant–caterpillar–parasitoid rearing dataset from
eight lowland rainforest localities in Papua New Guinea, plus climate data and a spatially explicit
metacommunity simulation.

Each analysis is split into an **analysis script** (`NN_name.R`, writes results to `output/rds/`)
and a **plotting script** (`NN_name_plots.R`, writes figures to `output/fig/`). Run the analysis
script before its plotting script.

## Repository structure

```
DATA/                     # inputs: MASTER.xlsx, monthly_data.csv, LocalData_all.csv.gz
Data_model/               # raw model output (LocalData*N0.txt) — optional once consolidated
model/                    # C++ metacommunity simulation (main.cpp, hosts.csv, parasitoids.csv, params.txt)
output/rds/               # intermediate results (.rds) + result tables (.csv)
output/fig/               # figures (.pdf/.png/.svg)
NN_*.R                    # analysis + plotting scripts
```

`output/` is created automatically.

## Data

| File | Description |
|------|-------------|
| `DATA/MASTER.xlsx` | One row per reared caterpillar: locality, plant, caterpillar, optional parasitoid, guild. Input for analyses 1–7 and 9. |
| `DATA/monthly_data.csv` | TerraClimate `tmin`/`tmax`/`ppt` at the 8 sites (monthly, 1970–2006); produced once by `08a`. |
| `DATA/LocalData_all.csv.gz` | All metacommunity-model output consolidated into one compressed CSV (times 50 and 100); produced by `10a`. |

`MASTER.xlsx` is loaded with `read_excel(..., guess_max = 1048576)` and the parasitoid column coerced
with `as.character()` — the column is sparse and otherwise gets silently dropped. The raw
`Data_model/*.txt` files are too numerous for GitHub; consolidate them once with `10a`.

## Requirements

R ≥ 4.1. Install packages once:

```r
install.packages(c(
  "dplyr","tidyr","purrr","tibble","readxl","here","data.table","R.utils","forcats",
  "vegan","betapart","CommEcol","bipartite","lme4","broom.mixed","terra",
  "ggplot2","ggpubr","ggsignif","patchwork","ggrepel","viridis","maps","scales","cowplot","svglite",
  "flextable","officer"
))
```

Optional: `QBMS` (downloading raw TerraClimate rasters), `iNEXT.beta3D` (exploratory `01b` only).

## How to run

Set the working directory to the repo root, ensure `DATA/` holds the inputs, then run each analysis
then its plotting script.

| # | Analysis | Scripts |
|---|----------|---------|
| 1, 3 | Spatial turnover of species + caterpillar subsampling | `01_spatial_turnover_species.R` → `_plots.R` |
| 2 | Rarity robustness | `02_rarity_parasitoids.R` → `_plots.R` |
| 4 | Interaction turnover (+ subsampling) | `04_interaction_turnover.R` → `_plots.R` |
| 5 | Enemy-free space | `05_enemy_free_space.R` → `_plots.R` |
| 6 | Parasitism-rate heatmap + descriptive stats | `06_parasitism_rates.R` → `_plots.R` |
| 6 | GLMM: parasitism ~ parasitoid diversity | `06_glmm_parasitism.R` → `_plots.R` |
| 7 | Network specialization (H₂′, d′) | `07_network_specialization.R` → `_plots.R` |
| 8 | Abiotic (climate) gradient | `08a_…` → `08_abiotic_gradient.R` → `_plots.R` (or `run_08.R`) |
| 9 | Temporal turnover at Ohu | `09_temporal_turnover_ohu.R` → `_plots.R` |
| 10 | Metacommunity model | `10a_…` → `10_metacommunity_model_fromCSV.R` → `10_…_plots.R` |

For analysis 8, `run_08.R` chains download → analysis → plots; `08_abiotic_gradient_combined.R` is a
single-file version. For analysis 10, `10a_consolidate_localdata.R` builds the CSV first; the raw
`LocalData*N0.txt` files are produced by the C++ simulation below.

## Conventions

- Analysis and plotting are separate scripts (`08_abiotic_gradient_combined.R` is the one exception).
- Subsampling (caterpillars in `01`, interactions in `04`) uses **999 iterations**, `set.seed(1234)`.
- Tables exported as Word (`.docx`) + CSV; all code and object names in English.

## Notes on methods

- Bray–Curtis dissimilarity via **`vegan`**; Sørensen/Chao–Sørensen via `betapart`/`CommEcol`.
- Mixed-model CIs are **Wald 95%** (`broom.mixed`); climate PCA uses base `prcomp()`.
- Climate variables at **full monthly resolution** (month × year, 1998 excluded).
- Interaction turnover via `bipartite::betalinkr_multi()` with `partitioning = "commondenom"`.

---

## Metacommunity Simulation — Host–Parasitoid Food Web

A stochastic patch-occupancy simulation of host–parasitoid metacommunities using the Gillespie
algorithm. Species colonize and go extinct across a spatially explicit landscape, with parasitoids
depending on hosts for persistence.

### Model Overview

The simulation tracks presence/absence of host and parasitoid species across a set of patches
arranged in a 2D space with **periodic (toroidal) boundary conditions**. Dynamics follow four event
types:

| Event | Description |
|---|---|
| Host colonization | A host colonizes an empty patch from a connected occupied patch |
| Host extinction | A host goes locally extinct (rate increases with parasitoid load) |
| Parasitoid colonization | A parasitoid colonizes a patch with a suitable host present |
| Parasitoid extinction | A parasitoid goes locally extinct (rate decreases with more host species) |

Timing between events is drawn from an exponential distribution (Gillespie algorithm). Connectivity
between patches is determined by dispersal distance `L` — patches within distance `L` of each other
are linked.

### Dependencies

- C++11 or later
- Standard library only (no external dependencies)

Compile with:

```bash
g++ -O2 -std=c++11 main.cpp -o metacom
```

### Input Files

#### 1. Standard input (simulation parameters)

Parameters are read from `stdin` in this exact order:

| Parameter | Type | Description |
|---|---|---|
| `interval` | int | Time step interval |
| `totalruntime` | int | Wall-clock time limit (seconds) |
| `maxT` | int | Maximum simulation time steps |
| `rec` | int | Recording interval (time steps between snapshots) |
| `N` | int | Number of independent replicates |
| `set` | int | Simulation set ID (used in output filenames) |
| `delta` | double | Effect of each parasitoid on host extinction rate |
| `gamma` | double | Effect of hosts on parasitoid survival (reduces extinction) |
| `U` | double | Side length of spatial domain |
| `numvil` | int | Number of patches |
| `pcol` | int | Parasitoid colonization mode: `0` = fixed rate, `1` = density-dependent |
| `hostfile` | string | Path to host species CSV file |
| `parafile` | string | Path to parasitoid species CSV file |

#### 2. Host species file (CSV)

One row per host species:

```
extinction_rate, colonization_rate, dispersal_distance, parasitoid_index_1, parasitoid_index_2, ...
```

Example (`hosts.csv`):

```
0.1, 0.3, 0.4, 0, 1
0.1, 0.3, 0.4, 1, 2
0.1, 0.3, 0.4, 0, 2
```

#### 3. Parasitoid species file (CSV)

One row per parasitoid species:

```
extinction_rate, colonization_rate, dispersal_distance, host_index_1, host_index_2, ...
```

Example (`parasitoids.csv`):

```
0.2, 0.5, 0.2, 0, 1
0.2, 0.5, 0.2, 1, 2
0.2, 0.5, 0.2, 0, 2
```

> Indices in host and parasitoid files are **zero-based** and must be mutually consistent (host row
> `i` lists its parasitoids; parasitoid row `j` lists its hosts).

### Output Files

#### Console output (`stdout`)

One line per time step per replicate:

```
T   HTot_0   HTot_1   ...   PTot_0   PTot_1   ...   Hcol   Hext   Pcol   Pext
```

Where `HTot_i` and `PTot_i` are the number of patches occupied by host/parasitoid species `i`, and
`Hcol`, `Hext`, `Pcol`, `Pext` are total colonization/extinction rates.

#### Log files (`LocalData<set>N<replicate>.txt`)

Spatial snapshot of all patches at each recording interval:

```
T   x   y   Hpres_0   Hpres_1   ...   Ppres_0   Ppres_1   ...
```

One file is created per replicate, e.g. `LocalData1N0.txt`, `LocalData1N1.txt`, etc.

### Usage Example

#### Parameter file (`params.txt`)

```
1        # interval
3600     # totalruntime (seconds)
10000    # maxT
100      # rec
5        # N (replicates)
1        # set
0.5      # delta
0.3      # gamma
1.0      # U (spatial domain size)
50       # numvil (patches)
0        # pcol (0 = fixed colonization rate)
hosts.csv
parasitoids.csv
```

#### Run

```bash
./metacom < params.txt > output.txt
```

#### Redirect replicate log files

Log files are written automatically to the working directory. To keep things tidy:

```bash
mkdir results
cd results
../metacom < ../params.txt > output.txt
```

### Parameter Guide

#### `delta` — parasitoid pressure on hosts

Higher values increase host extinction when parasitoids are present. Set to `0` to decouple trophic
levels.

#### `gamma` — host support for parasitoids

Increases parasitoid persistence when more host species are present in a patch. The per-parasitoid
extinction rate is `xP + gamma / n_hosts`. Higher `gamma` = stronger host-dependence.

#### `pcol` — parasitoid colonization mode

- `0`: Each parasitoid colonizes at a fixed rate `mP` regardless of source patch abundance.
- `1`: Colonization rate scales with the number of occupied source patches (`Pcol_count × mP`), i.e.
  propagule rain.

#### `U` and `numvil`

Patches are placed randomly in a `U × U` square with toroidal boundaries. Increasing `U` while
keeping `numvil` constant reduces connectivity. Dispersal distances `LH` and `LP` in the species
files determine which patch pairs are connected.


### Notes

- The simulation exits early if wall-clock time exceeds `0.9 × totalruntime`.
- Patches are re-initialized at the start of each replicate; spatial coordinates are re-randomized.
- A parasitoid is immediately extinguished if it colonizes a patch with no suitable hosts.

---

## Contact

_add corresponding author / contact_

## License

_add license (e.g. MIT, CC-BY)_
