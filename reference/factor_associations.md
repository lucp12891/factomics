# Test association between factors and sample covariates

For every method, factor and covariate, tests whether the factor scores
differ by the covariate: Kruskal-Wallis for categorical covariates,
Spearman correlation for numeric ones.

## Usage

``` r
factor_associations(
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

A data frame with `method`, `factor`, `covariate`, `test`, `statistic`,
`p_value` and BH-adjusted `p_adj`.
