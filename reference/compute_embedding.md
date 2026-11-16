# Compute a 2-D (or n-D) embedding of a samples x features matrix

A thin, consistent wrapper around PCA
([`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)), UMAP
([`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html))
and t-SNE ([`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html))
with sensible defaults for the small sample sizes typical of multi-omics
studies.

## Usage

``` r
compute_embedding(
  x,
  type = c("umap", "tsne", "pca"),
  n_components = 2,
  scale = FALSE,
  n_neighbors = NULL,
  min_dist = 0.1,
  perplexity = NULL,
  seed = 21,
  ...
)
```

## Arguments

- x:

  Numeric matrix, **samples x features** (e.g. factor scores).

- type:

  One of `"umap"`, `"tsne"` or `"pca"`.

- n_components:

  Number of output dimensions (default `2`).

- scale:

  Logical. Standardise columns before embedding (default `FALSE`, which
  keeps each factor's own variance).

- n_neighbors:

  UMAP neighbourhood size. Default `min(15, n_samples - 1)`.

- min_dist:

  UMAP `min_dist` (default `0.1`).

- perplexity:

  t-SNE perplexity. Default chosen from the sample size (at most 30, and
  below `(n_samples - 1) / 3`).

- seed:

  Random seed (default `21`).

- ...:

  Further arguments passed to
  [`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html)
  or [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html).

## Value

A data frame with a `sample` column and `Dim1`, `Dim2`, ...

## Examples

``` r
sim <- simulate_multiomics(n_samples = 40, seed = 1)
x <- t(do.call(rbind, sim$data))
head(compute_embedding(x, "pca"))
#>   sample       Dim1       Dim2
#> 1   S001 -18.492770  10.975741
#> 2   S002   1.748707 -25.579452
#> 3   S003  34.324237   7.862819
#> 4   S004 -47.844100  13.332657
#> 5   S005   2.307836 -17.974998
#> 6   S006  42.956048  19.336523
```
