# Variance explained by each factor in each view

Computes, for every method, the fraction of variance (R^2) of each
(centred) view explained by each factor, by regressing the view on the
factor scores. This works for any method, so methods can be compared on
the same scale.

## Usage

``` r
variance_explained(results, methods = names(results), data = NULL)
```

## Arguments

- results:

  Object returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- methods:

  Methods to include (default: all).

- data:

  Optional input data, if `results` no longer carries it.

## Value

A data frame with columns `method`, `view`, `factor` (`"F1"`, ..., or
`"Total"` for all factors jointly) and `R2`.
