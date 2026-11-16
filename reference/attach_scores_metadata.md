# Attach sample metadata to a factor-analysis results object (in-place)

For each requested method, merges the score matrix with a metadata data
frame and stores the result in `results[[method]]$scores_with_meta`.
Optionally replaces the original score matrix.

## Usage

``` r
attach_scores_metadata(
  results,
  metadata,
  sample_col = "sample",
  methods = names(results),
  store_name = "scores_with_meta",
  replace_scores = FALSE
)
```

## Arguments

- results:

  Named list returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- metadata:

  Metadata data frame (see
  [`merge_scores_with_metadata()`](https://lucp12891.github.io/factomics/reference/merge_scores_with_metadata.md)).

- sample_col:

  Character. Column holding sample IDs (default `"sample"`).

- methods:

  Character vector of methods to process.

- store_name:

  Character. Name under which the merged data frame is stored in each
  method's sub-list (default `"scores_with_meta"`).

- replace_scores:

  Logical. If `TRUE`, replace `scores` with the merged data frame
  (original matrix is saved as `scores_matrix`).

## Value

The modified `results` list.
