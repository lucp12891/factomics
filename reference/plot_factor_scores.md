# Scatter plot of two factors for one or more methods

Scatter plot of two factors for one or more methods

## Usage

``` r
plot_factor_scores(
  results,
  methods = names(results),
  factors = c(1, 2),
  metadata = NULL,
  color_by = NULL,
  sample_col = "sample"
)
```

## Arguments

- results:

  Object returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- methods:

  Methods to plot (default: all).

- factors:

  Integer vector of length 2: which factors to plot.

- metadata:

  Optional sample annotations.

- color_by:

  Column of `metadata` used for colour.

- sample_col:

  Sample-ID column of `metadata`.

## Value

A ggplot object.
