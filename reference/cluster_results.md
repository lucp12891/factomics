# Cluster the samples on each method's factor scores

Cluster the samples on each method's factor scores

## Usage

``` r
cluster_results(
  results,
  methods = names(results),
  cluster_method = c("kmeans", "hclust", "pam"),
  k = 2:6,
  include_data = TRUE,
  ...
)
```

## Arguments

- results:

  Object returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- methods:

  Methods to use (default: all in `results`).

- cluster_method:

  Clustering algorithm: `"kmeans"`, `"hclust"` or `"pam"` (see
  [`cluster_samples()`](https://lucp12891.github.io/factomics/reference/cluster_samples.md)).

- k:

  Candidate numbers of clusters (default `2:6`).

- include_data:

  Logical. Also cluster the concatenated input data as a baseline
  (default `TRUE`).

- ...:

  Passed to
  [`cluster_samples()`](https://lucp12891.github.io/factomics/reference/cluster_samples.md)
  (e.g. `scale`, `linkage`, `seed`).

## Value

A list of class `factomics_clusters` with `assignments` (long data
frame: `sample`, `method`, `cluster`), `summary` (chosen `k` and
silhouette per method), `by_k` (silhouette for every candidate `k`) and
`fits` (the individual
[`cluster_samples()`](https://lucp12891.github.io/factomics/reference/cluster_samples.md)
results).

## Examples

``` r
sim <- simulate_multiomics(n_samples = 45, seed = 1)
res <- run_factor_pipeline(sim$data, methods = c("PCA", "MFA"), k = 3,
                           verbose = FALSE)
cl  <- cluster_results(res, k = 2:5)
cl$summary
#>              method algorithm k silhouette
#> 1               PCA    kmeans 4  0.4677575
#> 2               MFA    kmeans 4  0.4661464
#> 3 Concatenated data    kmeans 3  0.1333660
compare_clusterings(cl, truth = setNames(sim$metadata$group,
                                         sim$metadata$sample))
#> $agreement
#>                         PCA       MFA Concatenated data
#> PCA               1.0000000 1.0000000         0.8181818
#> MFA               1.0000000 1.0000000         0.8181818
#> Concatenated data 0.8181818 0.8181818         1.0000000
#> 
#> $vs_truth
#>              method       ARI
#> 1               PCA 0.8181818
#> 2               MFA 0.8181818
#> 3 Concatenated data 1.0000000
#> 
```
