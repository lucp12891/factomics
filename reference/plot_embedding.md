# Plot embeddings coloured by metadata or clusters

Plot embeddings coloured by metadata or clusters

## Usage

``` r
plot_embedding(
  emb,
  metadata = NULL,
  color_by = NULL,
  shape_by = NULL,
  sample_col = "sample",
  clusters = NULL,
  point_size = 2.2,
  ncol = NULL
)
```

## Arguments

- emb:

  Data frame from
  [`embed_results()`](https://lucp12891.github.io/factomics/reference/embed_results.md)
  or
  [`compute_embedding()`](https://lucp12891.github.io/factomics/reference/compute_embedding.md).

- metadata:

  Optional data frame of sample annotations.

- color_by:

  Column of `metadata` (or of `emb`) used for colour.

- shape_by:

  Optional column used for point shape.

- sample_col:

  Name of the sample-ID column in `metadata` (default `"sample"`).

- clusters:

  Optional output of
  [`cluster_results()`](https://lucp12891.github.io/factomics/reference/cluster_results.md);
  when given, points are coloured by each method's own cluster
  assignment.

- point_size:

  Point size (default `2.2`).

- ncol:

  Number of facet columns (default: automatic).

## Value

A ggplot object.
