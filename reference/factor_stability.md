# Stability of factor methods across random seeds

Fits the chosen methods once per seed and measures how similar the
results are: for every pair of seeds the factors are matched one-to-one
(Hungarian algorithm on the absolute correlation of loadings) and the
matched correlations are averaged. A method whose result does not depend
on the seed has stability 1. If the true factor scores are known (for
example from
[`simulate_multiomics()`](https://lucp12891.github.io/factomics/reference/simulate_multiomics.md)),
the recovery of the true factors is reported for every seed as well.

## Usage

``` r
factor_stability(
  X,
  methods = c("FABIA", "FABIA_ENS"),
  k = 5,
  seeds = 1:5,
  truth = NULL,
  verbose = FALSE,
  ...
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

- seeds:

  Integer vector of seeds (at least two).

- truth:

  Optional matrix of true factor scores (samples x factors).

- verbose:

  Logical. Print progress messages (default `TRUE`).

- ...:

  Further arguments passed to
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

## Value

A list of class `factomics_stability` with `pairs` (one row per method
and pair of seeds), `recovery` (one row per method and seed, if `truth`
is given), `summary` (per method) and `fits` (the results per seed).

## Details

This is how plain FABIA and ensemble FABIA (`"FABIA_ENS"`) can be
compared: a single FABIA run depends on its starting values, the
ensemble much less.

## See also

[`plot_factor_stability()`](https://lucp12891.github.io/factomics/reference/plot_factor_stability.md)

## Examples

``` r
# \donttest{
sim <- simulate_multiomics(n_samples = 40, seed = 1)
st <- factor_stability(sim$data, methods = c("PCA", "MFA"), k = 3,
                       seeds = 1:3, truth = sim$truth$Z, verbose = FALSE)
st$summary
#>   method mean_stability min_stability mean_recovery  sd_recovery
#> 1    PCA              1             1     0.8059550 0.000000e+00
#> 2    MFA              1             1     0.8046031 5.744449e-13
# }
```
