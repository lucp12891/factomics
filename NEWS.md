# factomics 0.2.0

* New method **RABFA**: the graph-guided Bayesian factor analysis with robust
  and adaptive shrinkage priors of Zhang, Chang, Shen and Long (2024,
  *Biometrics*). Available as `methods = "RABFA"` in `run_factor_pipeline()`
  and `factomics()`, or directly with `rabfa()`. Pass a feature graph (for
  example pathway edges) with `rabfa_graph`; `build_feature_graph()` builds a
  data-driven graph when no pathway knowledge is at hand. Hyperparameters are
  set with `rabfa_control()`. This is an independent implementation of the
  published model (Gaussian data).
* `factomics-cli` gains `--graph FILE` (edge CSV with `from`, `to`) for RABFA.

# factomics 0.1.0

First public release.

## New methods
* `run_factor_pipeline()` now also runs block-weighted **PCA**, **SGCCA**
  (mixOmics) and **JIVE** (r.jive), alongside FABIA, MOFA2, GFA and MFA.
* `list_factor_methods()` lists the methods and whether each is installed.
* A failing method is skipped with a warning (`stop_on_error = FALSE`) instead of
  stopping the whole pipeline.
* New `mofa_python` and `mofa_use_basilisk` arguments for MOFA2's Python backend.

## New downstream tools
* Embeddings: `compute_embedding()`, `embed_results()` and `plot_embedding()`
  (PCA, UMAP, t-SNE).
* Clustering: `cluster_samples()`, `cluster_results()`, `compare_clusterings()`,
  `adjusted_rand_index()` and `plot_clusters()` (k-means, hierarchical, PAM,
  with k chosen by silhouette).
* Interpretation: `variance_explained()`, `plot_variance_explained()`,
  `factor_associations()`, `plot_factor_associations()`, `plot_top_loadings()`,
  `plot_factor_heatmap()`, `plot_factor_scores()`.
* `factomics()` runs the whole workflow in one call; `export_results()` writes
  CSV tables.
* `simulate_multiomics()` generates multi-omics data with known factors.
* Docker image with every method pre-installed (including mofapy2) and a
  command-line tool, `factomics-cli`.

## Fixes
* The package can be installed again: analysis code that had been left at the
  top level of `R/pipeline.R` ran (and failed) at build time.
* The package now sits at the repository root, so
  `remotes::install_github("lucp12891/factomics")` works.
* `attach_scores_metadata()`, `merge_scores_with_metadata()` and
  `benchmark_factor_methods()` use every method in the results by default, not
  only the original four, and keep the original sample order.
* Factor matching in the grid plots matches samples by ID, not by row position.
* `plot_scores_scatter_grid()` and `plot_scores_lollipop_grid()` accept any number
  of methods; `align_to = "auto"` picks MOFA when present.
* MOFA2 scores are returned in the input sample order, and its HDF5 file is
  written to `tempdir()` by default instead of the working directory.

## Earlier development

* FABIA, MOFA2, GFA and MFA behind `run_factor_pipeline()`, with
  benchmarking and score plots.
