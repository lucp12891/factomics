# Compare clusterings across methods (and against known labels)

Compare clusterings across methods (and against known labels)

## Usage

``` r
compare_clusterings(clusters, truth = NULL)
```

## Arguments

- clusters:

  Output of
  [`cluster_results()`](https://lucp12891.github.io/factomics/reference/cluster_results.md).

- truth:

  Optional named vector of true labels (names = sample IDs).

## Value

A list with `agreement`, a methods x methods matrix of adjusted Rand
indices, and (if `truth` is given) `vs_truth`, a data frame of ARI per
method against the known labels.
