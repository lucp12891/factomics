# Embed the factor scores of every method (and the raw data)

Runs
[`compute_embedding()`](https://lucp12891.github.io/factomics/reference/compute_embedding.md)
on the factor scores of each method in a
[`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md)
result, and optionally on the concatenated input data, so that methods
can be compared visually in one figure.

## Usage

``` r
embed_results(
  results,
  type = c("umap", "tsne", "pca"),
  methods = names(results),
  include_data = TRUE,
  ...
)
```

## Arguments

- results:

  Object returned by
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- type:

  Character vector, any of `"umap"`, `"tsne"`, `"pca"`.

- methods:

  Methods to embed (default: all in `results`).

- include_data:

  Logical. Also embed the concatenated preprocessed views (labelled
  `"Concatenated data"`), as a baseline (default `TRUE`).

- ...:

  Passed to
  [`compute_embedding()`](https://lucp12891.github.io/factomics/reference/compute_embedding.md).

## Value

A long data frame with columns `sample`, `Dim1`, `Dim2`, `method` and
`type`, of class `factomics_embedding`.

## Examples

``` r
sim <- simulate_multiomics(n_samples = 40, seed = 1)
res <- run_factor_pipeline(sim$data, methods = c("PCA", "MFA"), k = 4,
                           verbose = FALSE)
emb <- embed_results(res, type = "pca")
plot_embedding(emb, metadata = sim$metadata, color_by = "group")
```
