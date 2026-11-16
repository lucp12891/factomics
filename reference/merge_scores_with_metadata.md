# Merge factor scores with sample metadata

Performs a join between the scores matrices in a factor-analysis result
object and a metadata data frame, returning one wide data frame per
method as well as a stacked long table.

## Usage

``` r
merge_scores_with_metadata(
  results,
  metadata,
  sample_col = "sample",
  methods = names(results),
  join = c("left", "inner")
)
```

## Arguments

- results:

  Named list returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- metadata:

  Data frame (or coercible object) with at least one column identifying
  samples. Can also use `rownames` as sample IDs.

- sample_col:

  Character. Column name in `metadata` holding sample IDs (default
  `"sample"`). If not present, `rownames(metadata)` is used.

- methods:

  Character vector of method names to process (default: all methods
  found in `results`).

- join:

  `"left"` (default) or `"inner"`. Join type between scores and
  metadata.

## Value

A list with:

- `by_method`:

  Named list of data frames (one per method).

- `all_long`:

  Single row-stacked data frame with a `method` column.
