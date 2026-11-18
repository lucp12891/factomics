# Run multi-omics factor analysis pipeline

Runs one or more latent factor methods on a (list of) numeric
matrix/matrices and returns scores and loadings in a single, unified
format, so that any downstream function (benchmarking, embeddings,
clustering, plots) works the same way whatever method produced the
factors.

## Usage

``` r
run_factor_pipeline(
  X,
  methods = c("FABIA", "MOFA", "GFA", "MFA"),
  k = 10,
  scale_features = TRUE,
  mfa_group_sizes = NULL,
  fabia_alpha = 0.01,
  fabia_cyc = 1000,
  fabia_ens_control = list(),
  mofa_outfile = file.path(tempdir(), "mofa_model.hdf5"),
  mofa_convergence = "slow",
  mofa_maxiter = 2000,
  mofa_use_basilisk = FALSE,
  mofa_python = NULL,
  gfa_iter_burnin = 500,
  gfa_iter_max = 1000,
  pca_block_weight = TRUE,
  sgcca_keep = NULL,
  jive_rank_individual = 1,
  rabfa_graph = NULL,
  rabfa_iter = 2000,
  rabfa_burnin = 1000,
  rabfa_control = list(),
  seed = 21,
  verbose = TRUE,
  stop_on_error = FALSE
)
```

## Arguments

- X:

  Either a single numeric matrix (**features x samples**) or a named
  list of such matrices. All matrices must share the same sample names
  (column names).

- methods:

  Character vector of method names to run, case-insensitive. Any subset
  of
  `c("PCA", "MFA", "FABIA", "FABIA_ENS", "GFA", "MOFA", "SGCCA", "JIVE", "RABFA")`;
  see
  [`list_factor_methods()`](https://lucp12891.github.io/factomics/reference/list_factor_methods.md).
  Defaults to `c("FABIA", "MOFA", "GFA", "MFA")`.

- k:

  Integer. Target number of latent factors (default `10`).

- scale_features:

  Logical. If `TRUE` (default), each feature is z-scored before fitting.
  Features with zero variance become 0.

- mfa_group_sizes:

  Integer vector. Column group sizes passed to
  [`FactoMineR::MFA()`](https://rdrr.io/pkg/FactoMineR/man/MFA.html). If
  `NULL` (default), derived from view sizes.

- fabia_alpha:

  Numeric. Sparseness parameter for FABIA (default `0.01`).

- fabia_cyc:

  Integer. Number of FABIA cycles (default `1000`). Also used by
  `"FABIA_ENS"`, together with `fabia_alpha`.

- fabia_ens_control:

  List of settings for ensemble FABIA (`"FABIA_ENS"`): number of runs,
  subsampling, similarity and reliability thresholds; see
  [`fabia_ensemble_control()`](https://lucp12891.github.io/factomics/reference/fabia_ensemble_control.md).

- mofa_outfile:

  Character. Path for the MOFA2 HDF5 output file. Defaults to a file in
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- mofa_convergence:

  Character. Convergence mode for MOFA2: `"fast"`, `"medium"`, or
  `"slow"` (default `"slow"`).

- mofa_maxiter:

  Integer. Maximum MOFA2 iterations (default `2000`).

- mofa_use_basilisk:

  Logical. Let MOFA2 install its own Python environment through basilisk
  (default `FALSE`, i.e. use the Python found by reticulate, which must
  have `mofapy2` installed).

- mofa_python:

  Optional path to the Python binary that has `mofapy2` installed;
  passed to
  [`reticulate::use_python()`](https://rstudio.github.io/reticulate/reference/use_python.html)
  before MOFA2 runs.

- gfa_iter_burnin:

  Integer. GFA burn-in iterations (default `500`).

- gfa_iter_max:

  Integer. GFA maximum iterations (default `1000`).

- pca_block_weight:

  Logical. For `"PCA"`, divide each view by the square root of its
  number of features so that large views do not dominate (default
  `TRUE`).

- sgcca_keep:

  Optional integer vector (one value per view) giving the number of
  features kept per component in SGCCA. `NULL` (default) keeps all
  features (no sparsity).

- jive_rank_individual:

  Integer. Rank of each view's individual structure for JIVE (default
  `1`).

- rabfa_graph:

  Optional feature graph for RABFA (edge data frame, adjacency matrix,
  or per-view list); see
  [`rabfa()`](https://lucp12891.github.io/factomics/reference/rabfa.md)
  and
  [`build_feature_graph()`](https://lucp12891.github.io/factomics/reference/build_feature_graph.md).
  `NULL` (default) fits RABFA without a graph.

- rabfa_iter, rabfa_burnin:

  Total MCMC iterations and burn-in for RABFA (defaults `2000` and
  `1000`).

- rabfa_control:

  List of RABFA hyperparameters; see
  [`rabfa_control()`](https://lucp12891.github.io/factomics/reference/rabfa_control.md).

- seed:

  Integer. Random seed (default `21`).

- verbose:

  Logical. Print progress messages (default `TRUE`).

- stop_on_error:

  Logical. If `FALSE` (default) a method that fails is skipped with a
  warning and the other methods are still returned.

## Value

A named list of class `factomics_results`, with one entry per fitted
method. Each entry contains:

- `scores`:

  Matrix of sample factor scores (**samples x k**).

- `loadings`:

  Matrix of feature loadings (**features x k**). Row names use
  `"view::feature"` notation.

- `meta`:

  Method-specific metadata (number of factors, settings, etc.).

- `plots`:

  (MOFA only) List of diagnostic ggplot objects.

The preprocessed views are kept in `attr(x, "views")` and are used by
[`variance_explained()`](https://lucp12891.github.io/factomics/reference/variance_explained.md)
and
[`embed_results()`](https://lucp12891.github.io/factomics/reference/embed_results.md).

## Examples

``` r
sim <- simulate_multiomics(n_samples = 40, seed = 1)
res <- run_factor_pipeline(sim$data, methods = c("PCA", "MFA"), k = 3,
                           verbose = FALSE)
res
#> <factomics_results>
#>   40 samples, 3 view(s): transcriptome (200), metabolome (80), microbiome (50)
#>   PCA       3 factors
#>   MFA       3 factors
head(res$MFA$scores)
#>              F1         F2         F3
#> S001 -1.6155192  0.7839443  0.4786979
#> S002 -0.1738509 -1.9503153  0.2820183
#> S003  2.6580348  0.3970507 -0.3917322
#> S004 -3.4817133  1.2120201  1.0907031
#> S005  0.2965065 -1.4153488 -0.9479438
#> S006  2.7317508  1.1183909  1.1629432
```
