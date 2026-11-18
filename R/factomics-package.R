#' factomics: one toolbox for multi-omics factor analysis
#'
#' factomics wraps several multi-omics latent factor methods (PCA, MFA, FABIA,
#' GFA, MOFA2, SGCCA, JIVE and the graph-guided RABFA of Zhang et al. 2024)
#' behind a single function,
#' [run_factor_pipeline()], and returns their results in one common format.
#' On top of that it provides embeddings (PCA, UMAP, t-SNE), clustering,
#' method benchmarking, factor interpretation and plotting, and a one-call
#' workflow, [factomics()].
#'
#' @keywords internal
#' @importFrom ggplot2 .data
"_PACKAGE"
