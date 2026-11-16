# Plot the top-weighted features of a factor

Plot the top-weighted features of a factor

## Usage

``` r
plot_top_loadings(results, method, factor = 1, n = 20, view = NULL)
```

## Arguments

- results:

  Object returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- method:

  Method name.

- factor:

  Factor index (default `1`).

- n:

  Number of features to show (default `20`).

- view:

  Optional view name to restrict to.

## Value

A ggplot object (bars coloured by view).
