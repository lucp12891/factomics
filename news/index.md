# Changelog

## factomics 0.1.0

First public release.

### Graph-guided factor analysis

- New method **RABFA**: the graph-guided Bayesian factor analysis with
  robust and adaptive shrinkage priors of Zhang, Chang, Shen and Long
  (2024, *Biometrics*). Available as `methods = "RABFA"` in
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md)
  and
  [`factomics()`](https://lucp12891.github.io/factomics/reference/factomics.md),
  or directly with
  [`rabfa()`](https://lucp12891.github.io/factomics/reference/rabfa.md).
  Pass a feature graph (for example pathway edges) with `rabfa_graph`;
  [`build_feature_graph()`](https://lucp12891.github.io/factomics/reference/build_feature_graph.md)
  builds a data-driven graph when no pathway knowledge is at hand.
  Hyperparameters are set with
  [`rabfa_control()`](https://lucp12891.github.io/factomics/reference/rabfa_control.md).
  This is an independent implementation of the published model (Gaussian
  data).
- `factomics-cli` gains `--graph FILE` (edge CSV with `from`, `to`) for
  RABFA.

### New methods

- [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md)
  now also runs block-weighted **PCA**, **SGCCA** (mixOmics) and
  **JIVE** (r.jive), alongside FABIA, MOFA2, GFA and MFA.
- [`list_factor_methods()`](https://lucp12891.github.io/factomics/reference/list_factor_methods.md)
  lists the methods and whether each is installed.
- A failing method is skipped with a warning (`stop_on_error = FALSE`)
  instead of stopping the whole pipeline.
- New `mofa_python` and `mofa_use_basilisk` arguments for MOFA2’s Python
  backend.

### New downstream tools

- Embeddings:
  [`compute_embedding()`](https://lucp12891.github.io/factomics/reference/compute_embedding.md),
  [`embed_results()`](https://lucp12891.github.io/factomics/reference/embed_results.md)
  and
  [`plot_embedding()`](https://lucp12891.github.io/factomics/reference/plot_embedding.md)
  (PCA, UMAP, t-SNE).
- Clustering:
  [`cluster_samples()`](https://lucp12891.github.io/factomics/reference/cluster_samples.md),
  [`cluster_results()`](https://lucp12891.github.io/factomics/reference/cluster_results.md),
  [`compare_clusterings()`](https://lucp12891.github.io/factomics/reference/compare_clusterings.md),
  [`adjusted_rand_index()`](https://lucp12891.github.io/factomics/reference/adjusted_rand_index.md)
  and
  [`plot_clusters()`](https://lucp12891.github.io/factomics/reference/plot_clusters.md)
  (k-means, hierarchical, PAM, with k chosen by silhouette).
- Interpretation:
  [`variance_explained()`](https://lucp12891.github.io/factomics/reference/variance_explained.md),
  [`plot_variance_explained()`](https://lucp12891.github.io/factomics/reference/plot_variance_explained.md),
  [`factor_associations()`](https://lucp12891.github.io/factomics/reference/factor_associations.md),
  [`plot_factor_associations()`](https://lucp12891.github.io/factomics/reference/plot_factor_associations.md),
  [`plot_top_loadings()`](https://lucp12891.github.io/factomics/reference/plot_top_loadings.md),
  [`plot_factor_heatmap()`](https://lucp12891.github.io/factomics/reference/plot_factor_heatmap.md),
  [`plot_factor_scores()`](https://lucp12891.github.io/factomics/reference/plot_factor_scores.md).
- [`factomics()`](https://lucp12891.github.io/factomics/reference/factomics.md)
  runs the whole workflow in one call;
  [`export_results()`](https://lucp12891.github.io/factomics/reference/export_results.md)
  writes CSV tables.
- [`simulate_multiomics()`](https://lucp12891.github.io/factomics/reference/simulate_multiomics.md)
  generates multi-omics data with known factors.
- Docker image with every method pre-installed (including mofapy2) and a
  command-line tool, `factomics-cli`.

### Fixes

- The package can be installed again: analysis code that had been left
  at the top level of `R/pipeline.R` ran (and failed) at build time.
- The package now sits at the repository root, so
  `remotes::install_github("lucp12891/factomics")` works.
- [`attach_scores_metadata()`](https://lucp12891.github.io/factomics/reference/attach_scores_metadata.md),
  [`merge_scores_with_metadata()`](https://lucp12891.github.io/factomics/reference/merge_scores_with_metadata.md)
  and
  [`benchmark_factor_methods()`](https://lucp12891.github.io/factomics/reference/benchmark_factor_methods.md)
  use every method in the results by default, not only the original
  four, and keep the original sample order.
- Factor matching in the grid plots matches samples by ID, not by row
  position.
- [`plot_scores_scatter_grid()`](https://lucp12891.github.io/factomics/reference/plot_scores_scatter_grid.md)
  and
  [`plot_scores_lollipop_grid()`](https://lucp12891.github.io/factomics/reference/plot_scores_lollipop_grid.md)
  accept any number of methods; `align_to = "auto"` picks MOFA when
  present.
- MOFA2 scores are returned in the input sample order, and its HDF5 file
  is written to [`tempdir()`](https://rdrr.io/r/base/tempfile.html) by
  default instead of the working directory.

### Earlier development

- FABIA, MOFA2, GFA and MFA behind
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md),
  with benchmarking and score plots.
