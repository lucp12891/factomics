# ============================================================
#  factomics — low-dimensional embeddings (PCA, UMAP, t-SNE)
# ============================================================

#' Compute a 2-D (or n-D) embedding of a samples x features matrix
#'
#' A thin, consistent wrapper around PCA ([stats::prcomp()]), UMAP
#' ([uwot::umap()]) and t-SNE ([Rtsne::Rtsne()]) with sensible defaults for
#' the small sample sizes typical of multi-omics studies.
#'
#' @param x Numeric matrix, **samples x features** (e.g. factor scores).
#' @param type One of `"umap"`, `"tsne"` or `"pca"`.
#' @param n_components Number of output dimensions (default `2`).
#' @param scale Logical. Standardise columns before embedding (default
#'   `FALSE`, which keeps each factor's own variance).
#' @param n_neighbors UMAP neighbourhood size. Default
#'   `min(15, n_samples - 1)`.
#' @param min_dist UMAP `min_dist` (default `0.1`).
#' @param perplexity t-SNE perplexity. Default chosen from the sample size
#'   (at most 30, and below `(n_samples - 1) / 3`).
#' @param seed Random seed (default `21`).
#' @param ... Further arguments passed to [uwot::umap()] or [Rtsne::Rtsne()].
#'
#' @return A data frame with a `sample` column and `Dim1`, `Dim2`, ...
#' @examples
#' sim <- simulate_multiomics(n_samples = 40, seed = 1)
#' x <- t(do.call(rbind, sim$data))
#' head(compute_embedding(x, "pca"))
#' @export
compute_embedding <- function(x, type = c("umap", "tsne", "pca"),
                              n_components = 2, scale = FALSE,
                              n_neighbors = NULL, min_dist = 0.1,
                              perplexity = NULL, seed = 21, ...) {
  type <- match.arg(tolower(type), c("umap", "tsne", "pca"))
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(rownames(x))) rownames(x) <- paste0("S", seq_len(nrow(x)))
  keep <- apply(x, 2, function(v) all(is.finite(v)) && stats::var(v) > 0)
  x <- x[, keep, drop = FALSE]
  if (ncol(x) < 1) stop("No informative columns to embed.", call. = FALSE)
  if (scale) x <- scale(x)
  n <- nrow(x)

  set.seed(seed)
  emb <- switch(
    type,
    pca = {
      nc <- min(n_components, ncol(x), n - 1)
      stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = nc)$x
    },
    umap = {
      .stop_pkg("uwot")
      nn <- if (is.null(n_neighbors)) min(15, n - 1) else n_neighbors
      uwot::umap(x, n_neighbors = nn, n_components = n_components,
                 min_dist = min_dist, seed = seed, ...)
    },
    tsne = {
      .stop_pkg("Rtsne")
      pp <- if (is.null(perplexity))
        max(1, min(30, floor((n - 1) / 3) - 1)) else perplexity
      Rtsne::Rtsne(x, dims = n_components, perplexity = pp,
                   check_duplicates = FALSE, pca = ncol(x) > 50,
                   ...)$Y
    }
  )
  emb <- as.matrix(emb)
  colnames(emb) <- paste0("Dim", seq_len(ncol(emb)))
  data.frame(sample = rownames(x), emb, row.names = NULL,
             stringsAsFactors = FALSE)
}


#' Embed the factor scores of every method (and the raw data)
#'
#' Runs [compute_embedding()] on the factor scores of each method in a
#' [run_factor_pipeline()] result, and optionally on the concatenated input
#' data, so that methods can be compared visually in one figure.
#'
#' @param results Object returned by [run_factor_pipeline()].
#' @param type Character vector, any of `"umap"`, `"tsne"`, `"pca"`.
#' @param methods Methods to embed (default: all in `results`).
#' @param include_data Logical. Also embed the concatenated preprocessed views
#'   (labelled `"Concatenated data"`), as a baseline (default `TRUE`).
#' @param ... Passed to [compute_embedding()].
#'
#' @return A long data frame with columns `sample`, `Dim1`, `Dim2`, `method`
#'   and `type`, of class `factomics_embedding`.
#' @examples
#' sim <- simulate_multiomics(n_samples = 40, seed = 1)
#' res <- run_factor_pipeline(sim$data, methods = c("PCA", "MFA"), k = 4,
#'                            verbose = FALSE)
#' emb <- embed_results(res, type = "pca")
#' plot_embedding(emb, metadata = sim$metadata, color_by = "group")
#' @export
embed_results <- function(results, type = c("umap", "tsne", "pca"),
                          methods = names(results), include_data = TRUE,
                          ...) {
  type <- tolower(type)
  .check_methods(results, methods)
  srcs <- lapply(stats::setNames(methods, methods),
                 function(m) results[[m]]$scores)
  if (include_data) {
    V <- attr(results, "views")
    if (!is.null(V)) srcs[["Concatenated data"]] <- t(.stack_views(V))
  }
  out <- list()
  for (ty in type) for (s in names(srcs)) {
    e <- compute_embedding(srcs[[s]], type = ty, ...)
    e$method <- s
    e$type   <- switch(ty, umap = "UMAP", tsne = "t-SNE", pca = "PCA")
    out[[length(out) + 1]] <- e
  }
  out <- do.call(rbind, out)
  out$method <- factor(out$method, levels = names(srcs))
  out$type   <- factor(out$type, levels = unique(out$type))
  class(out) <- c("factomics_embedding", "data.frame")
  out
}


#' Plot embeddings coloured by metadata or clusters
#'
#' @param emb Data frame from [embed_results()] or [compute_embedding()].
#' @param metadata Optional data frame of sample annotations.
#' @param color_by Column of `metadata` (or of `emb`) used for colour.
#' @param shape_by Optional column used for point shape.
#' @param sample_col Name of the sample-ID column in `metadata`
#'   (default `"sample"`).
#' @param clusters Optional output of [cluster_results()]; when given, points
#'   are coloured by each method's own cluster assignment.
#' @param point_size Point size (default `2.2`).
#' @param ncol Number of facet columns (default: automatic).
#'
#' @return A ggplot object.
#' @export
plot_embedding <- function(emb, metadata = NULL, color_by = NULL,
                           shape_by = NULL, sample_col = "sample",
                           clusters = NULL, point_size = 2.2, ncol = NULL) {
  df <- as.data.frame(emb)
  if (!is.null(metadata)) {
    meta <- .prep_metadata(metadata, sample_col)
    extra <- setdiff(names(meta), names(df))
    df <- cbind(df, meta[match(df$sample, meta$sample), extra, drop = FALSE])
  }
  if (!is.null(clusters)) {
    a <- clusters$assignments
    key <- paste(a$method, a$sample)
    df$cluster <- factor(a$cluster[match(paste(df$method, df$sample), key)])
    if (is.null(color_by)) color_by <- "cluster"
  }
  for (col in c(color_by, shape_by))
    if (!col %in% names(df))
      stop("Column '", col, "' not found in embedding or metadata.",
           call. = FALSE)

  aes_args <- list(x = quote(.data[["Dim1"]]), y = quote(.data[["Dim2"]]))
  if (!is.null(color_by)) aes_args$colour <- substitute(.data[[cc]],
                                                       list(cc = color_by))
  if (!is.null(shape_by)) aes_args$shape <- substitute(.data[[ss]],
                                                      list(ss = shape_by))
  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = point_size, alpha = 0.85) +
    ggplot2::labs(x = "Dimension 1", y = "Dimension 2") +
    theme_factomics()
  if (!is.null(color_by) && !is.numeric(df[[color_by]]))
    p <- p + ggplot2::scale_colour_manual(
      values = factomics_palette(length(unique(df[[color_by]]))))
  has_m <- "method" %in% names(df) && length(unique(df$method)) > 1
  has_t <- "type" %in% names(df) && length(unique(df$type)) > 1
  if (has_m && has_t) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data[["type"]]),
                                 cols = ggplot2::vars(.data[["method"]]),
                                 scales = "free")
  } else if (has_m) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[["method"]]),
                                 scales = "free", ncol = ncol)
  } else if (has_t) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[["type"]]),
                                 scales = "free", ncol = ncol)
  }
  p
}
