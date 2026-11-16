# Build a loading comparison object aligned across methods

Takes the `results` list returned by
[`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md)
and returns a named list of loading matrices with a common (union)
feature set, padded with zeros where a method does not have a particular
feature.

## Usage

``` r
build_loading_comparison(results)
```

## Arguments

- results:

  Named list as returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

## Value

Named list (one entry per method) of **features x k** matrices sharing
the same row-name order.
