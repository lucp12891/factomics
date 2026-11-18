# Robust and adaptive-shrinkage Bayesian factor analysis (RABFA)

Fits the graph-guided Bayesian factor model of Zhang et al. (2024) to
one or more omics views. Each loading has its own Laplace shrinkage
parameter; the log shrinkage parameters of features that are connected
in a feature graph (for example genes in the same pathway) are linked
through a sparse precision matrix that is itself learnt, which makes the
model robust to noisy or wrong edges. Without a graph it reduces to
adaptive-shrinkage Bayesian factor analysis.

## Usage

``` r
rabfa(
  X,
  k = 5,
  graph = NULL,
  iter = 2000,
  burnin = 1000,
  control = list(),
  scale_features = TRUE,
  seed = 21,
  verbose = TRUE
)
```

## Arguments

- X:

  Named list of **features x samples** matrices (or one matrix).

- k:

  Number of factors.

- graph:

  Optional feature graph: an edge data frame with columns `from` and
  `to` (see
  [`build_feature_graph()`](https://lucp12891.github.io/factomics/reference/build_feature_graph.md)),
  a feature x feature adjacency matrix, or a named list with one
  adjacency matrix per view. Feature names can be plain (`"TP53"`) or
  view-prefixed (`"rna::TP53"`). Edges between different views are
  ignored.

- iter:

  Total MCMC iterations (default `2000`).

- burnin:

  Burn-in iterations discarded (default `1000`).

- control:

  List of hyperparameters; see
  [`rabfa_control()`](https://lucp12891.github.io/factomics/reference/rabfa_control.md).

- scale_features:

  Z-score each feature first (default `TRUE`).

- seed:

  Random seed.

- verbose:

  Print progress.

## Value

A `factomics_results` object with one entry, `RABFA`, holding `scores`,
`loadings` (posterior means) and `meta`, which includes the posterior
mean shrinkage parameters (`lambda`), the loadings whose 95% credible
interval excludes zero (`selected`), a BIC-type criterion (`bic`) for
choosing `k` or hyperparameters, and the log-likelihood trace
(`loglik_trace`) for checking convergence.

## Details

This is the same model as `run_factor_pipeline(X, methods = "RABFA")`;
`rabfa()` is a convenience entry point. Inference is by MCMC: use more
iterations (the paper used 5,000-10,000) for final analyses.

## References

Zhang Q, Chang C, Shen L, Long Q (2024). Incorporating graph information
in Bayesian factor analysis with robust and adaptive shrinkage priors.
*Biometrics* 80(1), ujad014.
[doi:10.1093/biomtc/ujad014](https://doi.org/10.1093/biomtc/ujad014)

## Examples

``` r
sim <- simulate_multiomics(n_samples = 30,
                           views = c(rna = 30, met = 15), seed = 1)
g   <- build_feature_graph(sim$data, cutoff = 0.6)
fit <- rabfa(sim$data, k = 3, graph = g, iter = 300, burnin = 150,
             verbose = FALSE)
fit$RABFA$meta$n_selected
#> [1] 26
```
