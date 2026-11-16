# ============================================================
#  factomics — one-call workflow and export
# ============================================================

#' Run a complete multi-omics factor analysis workflow in one call
#'
#' Fits the requested factor methods with [run_factor_pipeline()], then (when
#' asked) attaches metadata, embeds the factor scores (PCA / UMAP / t-SNE),
#' clusters the samples on each method's factors, and benchmarks the methods
#' against each other.
#'
#' @param X Named list of **features x samples** matrices (or one matrix).
#' @param metadata Optional data frame of sample annotations.
#' @param sample_col Sample-ID column of `metadata` (default `"sample"`).
#' @param methods Factor methods to run; see [list_factor_methods()].
#' @param k Number of factors (default `5`).
#' @param embeddings Any of `"umap"`, `"tsne"`, `"pca"`; `NULL` to skip.
#' @param cluster_method `"kmeans"`, `"hclust"` or `"pam"`; `NULL` to skip.
#' @param n_clusters Candidate numbers of clusters (default `2:6`).
#' @param benchmark Logical. Run [benchmark_factor_methods()] when at least two
#'   methods were fitted (default `TRUE`).
#' @param seed Random seed (default `21`).
#' @param verbose Logical. Print progress (default `TRUE`).
#' @param ... Further arguments for [run_factor_pipeline()] (e.g.
#'   `fabia_alpha`, `gfa_iter_max`, `mofa_python`).
#'
#' @return An object of class `factomics` (a list) with elements `results`,
#'   `metadata`, `embeddings`, `clusters`, `cluster_agreement` and `benchmark`.
#' @examples
#' sim <- simulate_multiomics(n_samples = 45, seed = 1)
#' fx  <- factomics(sim$data, metadata = sim$metadata,
#'                  methods = c("PCA", "MFA"), k = 3,
#'                  embeddings = "pca", benchmark = FALSE, verbose = FALSE)
#' fx
#' @export
factomics <- function(X, metadata = NULL, sample_col = "sample",
                      methods = c("PCA", "MFA", "FABIA", "GFA"), k = 5,
                      embeddings = c("umap", "tsne"),
                      cluster_method = "kmeans", n_clusters = 2:6,
                      benchmark = TRUE, seed = 21, verbose = TRUE, ...) {
  res <- run_factor_pipeline(X, methods = methods, k = k, seed = seed,
                             verbose = verbose, ...)
  if (!length(res)) stop("No method could be fitted.", call. = FALSE)

  if (!is.null(metadata))
    res <- attach_scores_metadata(res, metadata, sample_col = sample_col)

  emb <- NULL
  if (length(embeddings)) {
    if (verbose) message("[factomics] Embedding factor scores: ",
                         paste(embeddings, collapse = ", "))
    emb <- embed_results(res, type = embeddings, seed = seed)
  }

  cl <- agree <- NULL
  if (!is.null(cluster_method)) {
    if (verbose) message("[factomics] Clustering samples (", cluster_method,
                         ")")
    cl <- cluster_results(res, cluster_method = cluster_method, k = n_clusters,
                          seed = seed)
    agree <- compare_clusterings(cl)$agreement
  }

  bench <- NULL
  if (benchmark && length(res) >= 2) {
    if (verbose) message("[factomics] Benchmarking methods")
    bench <- tryCatch(benchmark_factor_methods(res, methods = names(res)),
                      error = function(e) {
                        warning("Benchmark failed: ", conditionMessage(e),
                                call. = FALSE)
                        NULL
                      })
  }

  structure(list(results = res,
                 metadata = if (!is.null(metadata))
                   .prep_metadata(metadata, sample_col) else NULL,
                 embeddings = emb, clusters = cl,
                 cluster_agreement = agree, benchmark = bench),
            class = "factomics")
}

#' @export
print.factomics <- function(x, ...) {
  cat("<factomics workflow>\n")
  print(x$results)
  if (!is.null(x$embeddings))
    cat("  Embeddings:", paste(levels(x$embeddings$type), collapse = ", "),
        "\n")
  if (!is.null(x$clusters)) {
    cat("  Clustering:\n")
    print(x$clusters$summary, row.names = FALSE, digits = 3)
  }
  if (!is.null(x$benchmark)) {
    cat("  Benchmark (mean |r| of matched factor scores):\n")
    print(x$benchmark$summary_table_scores[, c("pair",
                                               "mean_abs_corr_scores")],
          row.names = FALSE, digits = 3)
  }
  invisible(x)
}


#' Write factor scores, loadings and summaries to CSV files
#'
#' @param x Output of [factomics()] or [run_factor_pipeline()].
#' @param dir Output directory (created if needed).
#' @return Invisibly, the paths of the files written.
#' @export
export_results <- function(x, dir = "factomics_output") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  res <- if (inherits(x, "factomics")) x$results else x
  files <- character(0)
  wr <- function(df, name, row_names = TRUE) {
    f <- file.path(dir, name)
    utils::write.csv(df, f, row.names = row_names)
    files <<- c(files, f)
  }
  for (m in names(res)) {
    wr(res[[m]]$scores,   paste0(m, "_scores.csv"))
    wr(res[[m]]$loadings, paste0(m, "_loadings.csv"))
  }
  if (!is.null(attr(res, "views")))
    wr(variance_explained(res), "variance_explained.csv", FALSE)
  if (inherits(x, "factomics")) {
    if (!is.null(x$clusters)) {
      wr(x$clusters$assignments, "clusters.csv", FALSE)
      wr(x$clusters$summary, "cluster_summary.csv", FALSE)
    }
    if (!is.null(x$embeddings))
      wr(as.data.frame(x$embeddings), "embeddings.csv", FALSE)
    if (!is.null(x$benchmark))
      wr(x$benchmark$summary_table_scores, "benchmark_scores.csv", FALSE)
  }
  invisible(files)
}
