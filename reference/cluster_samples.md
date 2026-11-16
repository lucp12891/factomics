# Cluster samples, choosing the number of clusters by silhouette

Cluster samples, choosing the number of clusters by silhouette

## Usage

``` r
cluster_samples(
  x,
  method = c("kmeans", "hclust", "pam"),
  k = 2:6,
  scale = FALSE,
  linkage = "ward.D2",
  nstart = 25,
  seed = 21
)
```

## Arguments

- x:

  Numeric matrix, **samples x features** (typically factor scores).

- method:

  `"kmeans"`, `"hclust"` (Ward linkage by default) or `"pam"`.

- k:

  Integer or integer vector of candidate cluster numbers. When several
  values are given, the one with the highest mean silhouette width is
  kept.

- scale:

  Logical. Standardise columns first (default `FALSE`, which keeps each
  factor's own variance so weak or noise factors weigh less).

- linkage:

  Linkage for `"hclust"` (default `"ward.D2"`).

- nstart:

  Number of random starts for k-means (default `25`).

- seed:

  Random seed (default `21`).

## Value

A list of class `factomics_clustering` with `cluster` (named integer
vector), `k`, `silhouette` (mean width of the chosen solution), `by_k`
(a table of mean silhouette width for every candidate `k`) and `method`.

## Examples

``` r
sim <- simulate_multiomics(n_samples = 45, seed = 1)
cl <- cluster_samples(t(sim$data$transcriptome), k = 2:5)
cl$k
#> [1] 3
table(cl$cluster, sim$metadata$group)
#>    
#>     G1 G2 G3
#>   1 15  0  0
#>   2  0  0 15
#>   3  0 15  0
```
