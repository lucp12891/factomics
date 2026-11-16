# Plot variance explained per factor and view

Plot variance explained per factor and view

## Usage

``` r
plot_variance_explained(
  results,
  methods = names(results),
  type = c("heatmap", "total"),
  data = NULL
)
```

## Arguments

- results:

  Object returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- methods:

  Methods to include (default: all).

- type:

  `"heatmap"` (factor x view R^2 per method) or `"total"` (total R^2 per
  view as bars).

- data:

  Optional input data, if `results` no longer carries it.

## Value

A ggplot object.
