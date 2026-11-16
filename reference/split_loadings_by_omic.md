# Split a loading data frame by feature prefix

Given a data frame with a `feature` column whose values follow the
`"view::feature_id"` convention, return a named list of subsets—one per
pattern match.

## Usage

``` r
split_loadings_by_omic(df, column = "feature", patterns = NULL)
```

## Arguments

- df:

  Data frame containing a feature identifier column.

- column:

  Character. Name of the column to match against (default `"feature"`).

- patterns:

  Named character vector of regex patterns. Names become the list
  element names.

## Value

Named list of data frames.

## Examples

``` r
if (FALSE) { # \dontrun{
split_loadings_by_omic(
  loading_df,
  patterns = c(microbiome = "^microbiome",
               metabolome = "^metabolome",
               lipidome   = "^lipidome")
)
} # }
```
