# Extract scores as a labelled data frame for one method

Pulls `scores_with_meta` from a method result entry, identifies the
factor column, and returns a tidy tibble suitable for plotting.

## Usage

``` r
scores_df_for_method(bench, method, factor_idx = 1, group_col = NULL)
```

## Arguments

- bench:

  Named list returned by
  [`attach_scores_metadata()`](https://lucp12891.github.io/factomics/reference/attach_scores_metadata.md).

- method:

  Character. Method name.

- factor_idx:

  Integer. Factor index (default `1`).

- group_col:

  Character or `NULL`. Metadata grouping column.

## Value

A [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with columns `sample`, `.sample_idx`, `.score`, `.group`, `method`.
