# Run a complete multi-omics factor analysis workflow in one call

Fits the requested factor methods with
[`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md),
then (when asked) attaches metadata, embeds the factor scores (PCA /
UMAP / t-SNE), clusters the samples on each method's factors, and
benchmarks the methods against each other.

## Usage

``` r
factomics(
  X,
  metadata = NULL,
  sample_col = "sample",
  methods = c("PCA", "MFA", "FABIA", "GFA"),
  k = 5,
  embeddings = c("umap", "tsne"),
  cluster_method = "kmeans",
  n_clusters = 2:6,
  benchmark = TRUE,
  seed = 21,
  verbose = TRUE,
  ...
)
```

## Arguments

- X:

  Named list of **features x samples** matrices (or one matrix).

- metadata:

  Optional data frame of sample annotations.

- sample_col:

  Sample-ID column of `metadata` (default `"sample"`).

- methods:

  Factor methods to run; see
  [`list_factor_methods()`](https://lucp12891.github.io/factomics/reference/list_factor_methods.md).

- k:

  Number of factors (default `5`).

- embeddings:

  Any of `"umap"`, `"tsne"`, `"pca"`; `NULL` to skip.

- cluster_method:

  `"kmeans"`, `"hclust"` or `"pam"`; `NULL` to skip.

- n_clusters:

  Candidate numbers of clusters (default `2:6`).

- benchmark:

  Logical. Run
  [`benchmark_factor_methods()`](https://lucp12891.github.io/factomics/reference/benchmark_factor_methods.md)
  when at least two methods were fitted (default `TRUE`).

- seed:

  Random seed (default `21`).

- verbose:

  Logical. Print progress (default `TRUE`).

- ...:

  Further arguments for
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md)
  (e.g. `fabia_alpha`, `gfa_iter_max`, `mofa_python`).

## Value

An object of class `factomics` (a list) with elements `results`,
`metadata`, `embeddings`, `clusters`, `cluster_agreement` and
`benchmark`.

## Examples

``` r
sim <- simulate_multiomics(n_samples = 45, seed = 1)
fx  <- factomics(sim$data, metadata = sim$metadata,
                 methods = c("PCA", "MFA"), k = 3,
                 embeddings = "pca", benchmark = FALSE, verbose = FALSE)
fx
#> <factomics workflow>
#> <factomics_results>
#>   45 samples, 3 view(s): transcriptome (200), metabolome (80), microbiome (50)
#>   PCA       3 factors
#>   MFA       3 factors
#>   Embeddings: PCA 
#>   Clustering:
#>             method algorithm k silhouette
#>                PCA    kmeans 4      0.468
#>                MFA    kmeans 4      0.466
#>  Concatenated data    kmeans 3      0.133
```
