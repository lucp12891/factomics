# factomics <img src="man/figures/logo.png" align="right" height="139" alt="" />

> A unified R package for \*\*multi-omics latent factor analysis\*\* supporting
> FABIA, MOFA2, GFA and MFA with a single consistent interface, method
> comparison, and publication-ready plots.

\---

## Overview

`factomics` wraps four established factor-analysis methods behind one
function call:

|Method|R package|Description|
|-|-|-|
|**FABIA**|`fabia`|Factor Analysis for Bicluster Acquisition|
|**MOFA2**|`MOFA2`|Multi-Omics Factor Analysis v2|
|**GFA**|`GFA`|Group Factor Analysis|
|**MFA**|`FactoMineR`|Multiple Factor Analysis|

A single input (a named list of feature × sample matrices) is accepted by
**all** methods.  Results share a common structure (`scores`, `loadings`,
`meta`) so that downstream comparison, visualisation, and integration with
sample metadata are trivial.

\---

## Installation

### From source

```r
# install.packages("devtools")
devtools::install\_github("yourname/factomics")
```

### Required dependencies (all methods)

```r
install.packages(c("fabia", "GFA", "FactoMineR", "clue",
                   "ggplot2", "reshape2", "patchwork",
                   "scales", "tibble", "dplyr", "tidyr"))
BiocManager::install("MOFA2")          # MOFA2 also needs mofapy2 (Python)
```

### MOFA2 Python backend

```bash
pip install mofapy2
```

\---

## Quick start

```r
library(factomics)

# ── 1. Prepare data ──────────────────────────────────────────────
#  Named list of numeric matrices, each features × samples
data\_list <- list(
  microbiome = as.matrix(t(dat\_16s)),    # features × samples
  metabolome = as.matrix(t(dat\_met)),
  lipidome   = as.matrix(t(dat\_lip))
)

# ── 2. Run methods ────────────────────────────────────────────────
results <- run\_factor\_pipeline(
  X               = data\_list,
  methods         = c("FABIA", "MOFA", "GFA", "MFA"),   # or a subset
  k               = 3,
  scale\_features  = TRUE,
  seed            = 21
)

# Access scores/loadings for a specific method
results$MOFA$scores     # samples × k matrix
results$MOFA$loadings   # (view::feature) × k matrix

# ── 3. Attach metadata ────────────────────────────────────────────
results <- attach\_scores\_metadata(results, metadata,
                                   sample\_col = "Mouse")

# ── 4. Plot factor scores coloured by treatment group ─────────────
plot\_scores\_scatter\_grid(
  bench      = results,
  factor\_idx = 1,
  group\_col  = "Treatment",
  methods    = c("MOFA", "FABIA", "MFA", "GFA"),
  align\_to   = "MOFA"
)

# ── 5. Benchmark / compare methods ───────────────────────────────
bench <- benchmark\_factor\_methods(results,
                                   methods = c("FABIA", "MOFA", "GFA", "MFA"))
print(bench$summary\_table\_scores)
print(bench$all\_scores\_heatmap)
```

\---

## API reference

### `run\_factor\_pipeline()`

|Argument|Default|Description|
|-|-|-|
|`X`|—|Matrix or named list of matrices (features × samples)|
|`methods`|`c("FABIA","MOFA","GFA","MFA")`|Methods to run (case-insensitive)|
|`k`|`10`|Number of latent factors|
|`scale\_features`|`TRUE`|Z-score features before fitting|
|`fabia\_alpha`|`0.01`|FABIA sparseness|
|`mofa\_outfile`|`"mofa\_model.hdf5"`|MOFA2 HDF5 output path|
|`mofa\_convergence`|`"slow"`|MOFA2 convergence mode|
|`gfa\_iter\_burnin`|`500`|GFA burn-in iterations|
|`gfa\_iter\_max`|`1000`|GFA max iterations|
|`seed`|`21`|Random seed|

**Returns:** Named list — one entry per method with `$scores`, `$loadings`, `$meta`.

\---

### `benchmark\_factor\_methods()`

Compares all fitted methods pairwise:

* Hungarian-matched factor correlations
* Per-omic loading correlations
* Heatmaps and scatter plots
* Optional AUC against ground-truth labels

\---

### `attach\_scores\_metadata()`

Joins sample metadata into `results\[\[method]]$scores\_with\_meta`.

\---

### `plot\_scores\_scatter\_grid()` / `plot\_scores\_lollipop\_grid()`

2 × 2 grids comparing factor scores across 4 methods, coloured by any
metadata variable.

\---

## Docker

### Build

```bash
docker build -t omicsfactor:latest .
```

### Run interactive R

```bash
docker run --rm -it \\
  -v $(pwd)/data:/work/data \\
  -v $(pwd)/output:/work/output \\
  omicsfactor:latest R
```

### Docker Compose (recommended)

```bash
# Interactive R shell
docker compose up r-session

# RStudio Server at http://localhost:8787  (user: rstudio  pwd: omicspassword)
docker compose --profile rstudio up rstudio
```

### Mount your data

Place your data files in `./data/` on the host; they appear inside the
container at `/work/data/`.

\---

## Output structure

```
results/
├── FABIA/
│   ├── scores    # samples × k matrix
│   ├── loadings  # (view::feature) × k matrix
│   └── meta      # alpha, cyc, n\_factors
├── MOFA/
│   ├── scores
│   ├── loadings
│   ├── meta      # variance\_explained, hdf5 path
│   └── plots     # variance\_by\_factor, factors\_by\_group
├── GFA/
│   ├── scores
│   ├── loadings
│   └── meta
└── MFA/
    ├── scores
    ├── loadings
    └── meta
```

\---

## Citation

If you use this package please cite the underlying method packages:

* **FABIA**: Hochreiter et al. (2010) *Bioinformatics*
* **MOFA2**: Argelaguet et al. (2020) *Genome Biology*
* **GFA**: Virtanen et al. (2012) *JMLR*
* **MFA**: Escofier \& Pagès (1994) *Computational Statistics \& Data Analysis*

\---

## License

MIT © 2024

