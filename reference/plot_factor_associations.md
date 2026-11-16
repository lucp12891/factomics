# Plot factor-covariate associations

Plot factor-covariate associations

## Usage

``` r
plot_factor_associations(
  results,
  metadata,
  covariates = NULL,
  methods = names(results),
  sample_col = "sample"
)
```

## Arguments

- results:

  Object returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- metadata:

  Sample annotations.

- covariates:

  Metadata columns to test (default: all except the sample column).

- methods:

  Methods to include (default: all).

- sample_col:

  Sample-ID column of `metadata`.

## Value

A ggplot heatmap of -log10 adjusted p-values.
