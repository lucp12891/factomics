# Heatmap of factor scores, ordered by a sample annotation

Heatmap of factor scores, ordered by a sample annotation

## Usage

``` r
plot_factor_heatmap(
  results,
  method,
  metadata = NULL,
  order_by = NULL,
  sample_col = "sample"
)
```

## Arguments

- results:

  Object returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- method:

  Method name.

- metadata:

  Optional sample annotations.

- order_by:

  Optional metadata column used to order (and label) samples.

- sample_col:

  Sample-ID column of `metadata`.

## Value

A ggplot object.
