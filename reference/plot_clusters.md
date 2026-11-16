# Plot clustering quality across methods

Plot clustering quality across methods

## Usage

``` r
plot_clusters(clusters, type = c("silhouette", "agreement"), truth = NULL)
```

## Arguments

- clusters:

  Output of
  [`cluster_results()`](https://lucp12891.github.io/factomics/reference/cluster_results.md).

- type:

  `"silhouette"` (mean silhouette width against `k` for each method) or
  `"agreement"` (heatmap of pairwise adjusted Rand indices).

- truth:

  Optional named vector of true labels; adds a "Truth" row and column to
  the agreement heatmap.

## Value

A ggplot object.
