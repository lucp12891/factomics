# Build a feature graph for RABFA from the data

RABFA uses a graph of features (genes, proteins, metabolites) that are
expected to act together, usually from a pathway database such as KEGG.
When no such knowledge is available, this helper builds a data-driven
graph by connecting features within the same view whose absolute
correlation is at least `cutoff`, keeping at most `max_neighbors`
strongest edges per feature. Prefer a knowledge-based graph when you
have one: a graph learnt from the same data partly re-uses the signal.

## Usage

``` r
build_feature_graph(X, cutoff = 0.5, max_neighbors = 5)
```

## Arguments

- X:

  Named list of **features x samples** matrices (or one matrix).

- cutoff:

  Minimum absolute Pearson correlation for an edge (default `0.5`).

- max_neighbors:

  Maximum number of edges kept per feature (default `5`).

## Value

A data frame of edges with columns `from` and `to`, using
`"view::feature"` names, ready for `rabfa_graph`.

## Examples

``` r
sim <- simulate_multiomics(n_samples = 30, seed = 1)
g <- build_feature_graph(sim$data, cutoff = 0.6)
head(g)
#>                      from                      to
#> 1 transcriptome::gene_006 transcriptome::gene_015
#> 2 transcriptome::gene_010 transcriptome::gene_015
#> 3 transcriptome::gene_009 transcriptome::gene_017
#> 4 transcriptome::gene_001 transcriptome::gene_025
#> 5 transcriptome::gene_016 transcriptome::gene_031
#> 6 transcriptome::gene_009 transcriptome::gene_040
```
