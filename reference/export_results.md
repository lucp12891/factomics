# Write factor scores, loadings and summaries to CSV files

Write factor scores, loadings and summaries to CSV files

## Usage

``` r
export_results(x, dir = "factomics_output")
```

## Arguments

- x:

  Output of
  [`factomics()`](https://lucp12891.github.io/factomics/reference/factomics.md)
  or
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- dir:

  Output directory (created if needed).

## Value

Invisibly, the paths of the files written.
