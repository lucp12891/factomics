# Package index

## One-call workflow

- [`factomics()`](https://lucp12891.github.io/factomics/reference/factomics.md)
  : Run a complete multi-omics factor analysis workflow in one call
- [`export_results()`](https://lucp12891.github.io/factomics/reference/export_results.md)
  : Write factor scores, loadings and summaries to CSV files

## Fitting factor models

- [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md)
  : Run multi-omics factor analysis pipeline
- [`list_factor_methods()`](https://lucp12891.github.io/factomics/reference/list_factor_methods.md)
  : List the factor-analysis methods supported by factomics
- [`fabia_ensemble()`](https://lucp12891.github.io/factomics/reference/fabia_ensemble.md)
  : Ensemble FABIA: robust biclusters and factors from many FABIA runs
- [`fabia_ensemble_control()`](https://lucp12891.github.io/factomics/reference/fabia_ensemble_control.md)
  : Settings for ensemble FABIA
- [`plot_fabia_ensemble()`](https://lucp12891.github.io/factomics/reference/plot_fabia_ensemble.md)
  : Plot the reliability of ensemble FABIA biclusters
- [`rabfa()`](https://lucp12891.github.io/factomics/reference/rabfa.md)
  : Robust and adaptive-shrinkage Bayesian factor analysis (RABFA)
- [`rabfa_control()`](https://lucp12891.github.io/factomics/reference/rabfa_control.md)
  : Default control settings for RABFA
- [`build_feature_graph()`](https://lucp12891.github.io/factomics/reference/build_feature_graph.md)
  : Build a feature graph for RABFA from the data
- [`simulate_multiomics()`](https://lucp12891.github.io/factomics/reference/simulate_multiomics.md)
  : Simulate a multi-omics data set with known latent factors

## Embeddings

- [`compute_embedding()`](https://lucp12891.github.io/factomics/reference/compute_embedding.md)
  : Compute a 2-D (or n-D) embedding of a samples x features matrix
- [`embed_results()`](https://lucp12891.github.io/factomics/reference/embed_results.md)
  : Embed the factor scores of every method (and the raw data)
- [`plot_embedding()`](https://lucp12891.github.io/factomics/reference/plot_embedding.md)
  : Plot embeddings coloured by metadata or clusters

## Clustering

- [`cluster_samples()`](https://lucp12891.github.io/factomics/reference/cluster_samples.md)
  : Cluster samples, choosing the number of clusters by silhouette
- [`cluster_results()`](https://lucp12891.github.io/factomics/reference/cluster_results.md)
  : Cluster the samples on each method's factor scores
- [`compare_clusterings()`](https://lucp12891.github.io/factomics/reference/compare_clusterings.md)
  : Compare clusterings across methods (and against known labels)
- [`adjusted_rand_index()`](https://lucp12891.github.io/factomics/reference/adjusted_rand_index.md)
  : Adjusted Rand index between two partitions
- [`plot_clusters()`](https://lucp12891.github.io/factomics/reference/plot_clusters.md)
  : Plot clustering quality across methods

## Interpretation

- [`variance_explained()`](https://lucp12891.github.io/factomics/reference/variance_explained.md)
  : Variance explained by each factor in each view
- [`plot_variance_explained()`](https://lucp12891.github.io/factomics/reference/plot_variance_explained.md)
  : Plot variance explained per factor and view
- [`factor_associations()`](https://lucp12891.github.io/factomics/reference/factor_associations.md)
  : Test association between factors and sample covariates
- [`plot_factor_associations()`](https://lucp12891.github.io/factomics/reference/plot_factor_associations.md)
  : Plot factor-covariate associations
- [`plot_top_loadings()`](https://lucp12891.github.io/factomics/reference/plot_top_loadings.md)
  : Plot the top-weighted features of a factor
- [`plot_factor_heatmap()`](https://lucp12891.github.io/factomics/reference/plot_factor_heatmap.md)
  : Heatmap of factor scores, ordered by a sample annotation
- [`plot_factor_scores()`](https://lucp12891.github.io/factomics/reference/plot_factor_scores.md)
  : Scatter plot of two factors for one or more methods

## Comparing methods

- [`benchmark_factor_methods()`](https://lucp12891.github.io/factomics/reference/benchmark_factor_methods.md)
  : Compare and benchmark multiple factor analysis results
- [`build_loading_comparison()`](https://lucp12891.github.io/factomics/reference/build_loading_comparison.md)
  : Build a loading comparison object aligned across methods
- [`plot_scores_scatter_grid()`](https://lucp12891.github.io/factomics/reference/plot_scores_scatter_grid.md)
  : Plot factor scores for multiple methods (scatter grid)
- [`plot_scores_lollipop_grid()`](https://lucp12891.github.io/factomics/reference/plot_scores_lollipop_grid.md)
  : Plot factor scores as lollipop charts for multiple methods
- [`factor_stability()`](https://lucp12891.github.io/factomics/reference/factor_stability.md)
  : Stability of factor methods across random seeds
- [`plot_factor_stability()`](https://lucp12891.github.io/factomics/reference/plot_factor_stability.md)
  : Plot the stability of factor methods across seeds

## Metadata and helpers

- [`attach_scores_metadata()`](https://lucp12891.github.io/factomics/reference/attach_scores_metadata.md)
  : Attach sample metadata to a factor-analysis results object
  (in-place)
- [`merge_scores_with_metadata()`](https://lucp12891.github.io/factomics/reference/merge_scores_with_metadata.md)
  : Merge factor scores with sample metadata
- [`scores_df_for_method()`](https://lucp12891.github.io/factomics/reference/scores_df_for_method.md)
  : Extract scores as a labelled data frame for one method
- [`split_loadings_by_omic()`](https://lucp12891.github.io/factomics/reference/split_loadings_by_omic.md)
  : Split a loading data frame by feature prefix
- [`theme_factomics()`](https://lucp12891.github.io/factomics/reference/theme_factomics.md)
  [`factomics_palette()`](https://lucp12891.github.io/factomics/reference/theme_factomics.md)
  : The factomics ggplot2 theme and colour palette
- [`factomics-package`](https://lucp12891.github.io/factomics/reference/factomics-package.md)
  : factomics: one toolbox for multi-omics factor analysis
