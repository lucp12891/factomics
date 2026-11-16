# ============================================================
#  factomics — clustering on factor scores
# ============================================================

#' Cluster samples, choosing the number of clusters by silhouette
#'
#' @param x Numeric matrix, **samples x features** (typically factor scores).
#' @param method `"kmeans"`, `"hclust"` (Ward linkage by default) or `"pam"`.
#' @param k Integer or integer vector of candidate cluster numbers. When
#'   several values are given, the one with the highest mean silhouette width
#'   is kept.
#' @param scale Logical. Standardise columns first (default `FALSE`, which
#'   keeps each factor's own variance so weak or noise factors weigh less).
#' @param linkage Linkage for `"hclust"` (default `"ward.D2"`).
#' @param nstart Number of random starts for k-means (default `25`).
#' @param seed Random seed (default `21`).
#'
#' @return A list of class `factomics_clustering` with `cluster` (named
#'   integer vector), `k`, `silhouette` (mean width of the chosen solution),
#'   `by_k` (a table of mean silhouette width for every candidate `k`) and
#'   `method`.
#' @examples
#' sim <- simulate_multiomics(n_samples = 45, seed = 1)
#' cl <- cluster_samples(t(sim$data$transcriptome), k = 2:5)
#' cl$k
#' table(cl$cluster, sim$metadata$group)
#' @export
cluster_samples <- function(x, method = c("kmeans", "hclust", "pam"),
                            k = 2:6, scale = FALSE, linkage = "ward.D2",
                            nstart = 25, seed = 21) {
  method <- match.arg(method)
  x <- as.matrix(x)
  if (is.null(rownames(x))) rownames(x) <- paste0("S", seq_len(nrow(x)))
  keep <- apply(x, 2, function(v) all(is.finite(v)) && stats::var(v) > 0)
  x <- x[, keep, drop = FALSE]
  if (scale) x <- scale(x)
  k <- sort(unique(as.integer(k)))
  k <- k[k >= 2 & k < nrow(x)]
  if (!length(k)) stop("No valid k: need 2 <= k < n_samples.", call. = FALSE)

  d  <- stats::dist(x)
  hc <- if (method == "hclust") stats::hclust(d, method = linkage) else NULL
  fit_k <- function(kk) {
    set.seed(seed)
    cl <- switch(method,
                 kmeans = stats::kmeans(x, centers = kk, nstart = nstart)$cluster,
                 hclust = stats::cutree(hc, k = kk),
                 pam    = cluster::pam(d, k = kk, diss = TRUE)$clustering)
    sil <- cluster::silhouette(cl, d)
    list(cluster = cl, sil = mean(sil[, "sil_width"]))
  }
  fits  <- lapply(k, fit_k)
  sils  <- vapply(fits, `[[`, numeric(1), "sil")
  best  <- which.max(sils)
  cl    <- fits[[best]]$cluster
  names(cl) <- rownames(x)
  structure(list(cluster = cl, k = k[best], silhouette = sils[best],
                 by_k = data.frame(k = k, silhouette = sils),
                 method = method, hclust = hc),
            class = "factomics_clustering")
}

#' @export
print.factomics_clustering <- function(x, ...) {
  cat(sprintf("<factomics_clustering> %s, k = %d, mean silhouette = %.3f\n",
              x$method, x$k, x$silhouette))
  print(table(cluster = x$cluster))
  invisible(x)
}


#' Cluster the samples on each method's factor scores
#'
#' @param results Object returned by [run_factor_pipeline()].
#' @param methods Methods to use (default: all in `results`).
#' @param cluster_method Clustering algorithm: `"kmeans"`, `"hclust"` or
#'   `"pam"` (see [cluster_samples()]).
#' @param k Candidate numbers of clusters (default `2:6`).
#' @param include_data Logical. Also cluster the concatenated input data as a
#'   baseline (default `TRUE`).
#' @param ... Passed to [cluster_samples()] (e.g. `scale`, `linkage`,
#'   `seed`).
#'
#' @return A list of class `factomics_clusters` with `assignments` (long data
#'   frame: `sample`, `method`, `cluster`), `summary` (chosen `k` and
#'   silhouette per method), `by_k` (silhouette for every candidate `k`) and
#'   `fits` (the individual [cluster_samples()] results).
#' @examples
#' sim <- simulate_multiomics(n_samples = 45, seed = 1)
#' res <- run_factor_pipeline(sim$data, methods = c("PCA", "MFA"), k = 3,
#'                            verbose = FALSE)
#' cl  <- cluster_results(res, k = 2:5)
#' cl$summary
#' compare_clusterings(cl, truth = setNames(sim$metadata$group,
#'                                          sim$metadata$sample))
#' @export
cluster_results <- function(results, methods = names(results),
                            cluster_method = c("kmeans", "hclust", "pam"),
                            k = 2:6, include_data = TRUE, ...) {
  cluster_method <- match.arg(cluster_method)
  .check_methods(results, methods)
  srcs <- lapply(stats::setNames(methods, methods),
                 function(m) results[[m]]$scores)
  if (include_data) {
    V <- attr(results, "views")
    if (!is.null(V)) srcs[["Concatenated data"]] <- t(.stack_views(V))
  }
  fits <- lapply(srcs, cluster_samples, method = cluster_method, k = k, ...)
  assignments <- do.call(rbind, lapply(names(fits), function(m)
    data.frame(sample = names(fits[[m]]$cluster), method = m,
               cluster = unname(fits[[m]]$cluster),
               stringsAsFactors = FALSE)))
  summary <- data.frame(
    method     = names(fits),
    algorithm  = vapply(fits, `[[`, character(1), "method"),
    k          = vapply(fits, `[[`, integer(1), "k"),
    silhouette = vapply(fits, `[[`, numeric(1), "silhouette"),
    row.names = NULL, stringsAsFactors = FALSE)
  by_k <- do.call(rbind, lapply(names(fits), function(m)
    cbind(method = m, fits[[m]]$by_k)))
  structure(list(assignments = assignments, summary = summary, by_k = by_k,
                 fits = fits),
            class = "factomics_clusters")
}

#' @export
print.factomics_clusters <- function(x, ...) {
  cat("<factomics_clusters>\n")
  print(x$summary, row.names = FALSE, digits = 3)
  invisible(x)
}


#' Adjusted Rand index between two partitions
#'
#' @param a,b Vectors of cluster labels of equal length.
#' @return A single number; 1 for identical partitions, about 0 for random.
#' @examples
#' adjusted_rand_index(c(1, 1, 2, 2), c("a", "a", "b", "b"))
#' @export
adjusted_rand_index <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  tab <- table(a[ok], b[ok])
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  c2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(c2(tab))
  sum_a  <- sum(c2(rowSums(tab)))
  sum_b  <- sum(c2(colSums(tab)))
  expected <- sum_a * sum_b / c2(n)
  maxidx   <- (sum_a + sum_b) / 2
  if (maxidx == expected) return(if (sum_ij == expected) 1 else 0)
  (sum_ij - expected) / (maxidx - expected)
}


#' Compare clusterings across methods (and against known labels)
#'
#' @param clusters Output of [cluster_results()].
#' @param truth Optional named vector of true labels (names = sample IDs).
#'
#' @return A list with `agreement`, a methods x methods matrix of adjusted
#'   Rand indices, and (if `truth` is given) `vs_truth`, a data frame of ARI
#'   per method against the known labels.
#' @export
compare_clusterings <- function(clusters, truth = NULL) {
  a <- clusters$assignments
  ms <- unique(a$method)
  lab <- lapply(stats::setNames(ms, ms), function(m) {
    s <- a[a$method == m, ]
    stats::setNames(s$cluster, s$sample)
  })
  samples <- names(lab[[1]])
  M <- outer(seq_along(ms), seq_along(ms), Vectorize(function(i, j)
    adjusted_rand_index(lab[[i]][samples], lab[[j]][samples])))
  dimnames(M) <- list(ms, ms)
  out <- list(agreement = M)
  if (!is.null(truth)) {
    if (is.null(names(truth)))
      stop("'truth' must be a named vector (names = sample IDs).",
           call. = FALSE)
    out$vs_truth <- data.frame(
      method = ms,
      ARI = vapply(ms, function(m)
        adjusted_rand_index(lab[[m]][samples], truth[samples]), numeric(1)),
      row.names = NULL)
  }
  out
}


#' Plot clustering quality across methods
#'
#' @param clusters Output of [cluster_results()].
#' @param type `"silhouette"` (mean silhouette width against `k` for each
#'   method) or `"agreement"` (heatmap of pairwise adjusted Rand indices).
#' @param truth Optional named vector of true labels; adds a "Truth" row and
#'   column to the agreement heatmap.
#' @return A ggplot object.
#' @export
plot_clusters <- function(clusters, type = c("silhouette", "agreement"),
                          truth = NULL) {
  type <- match.arg(type)
  if (type == "silhouette") {
    df <- clusters$by_k
    best <- clusters$summary
    return(
      ggplot2::ggplot(df, ggplot2::aes(.data[["k"]], .data[["silhouette"]],
                                       colour = .data[["method"]])) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_point(data = best, size = 4.5, shape = 21,
                            stroke = 1.2, fill = NA) +
        ggplot2::scale_colour_manual(
          values = factomics_palette(length(unique(df$method)))) +
        ggplot2::scale_x_continuous(breaks = unique(df$k)) +
        ggplot2::labs(x = "Number of clusters (k)",
                      y = "Mean silhouette width", colour = "Method",
                      title = "Choosing k (circled: selected)") +
        theme_factomics()
    )
  }
  if (!is.null(truth)) {
    clusters$assignments <- rbind(
      clusters$assignments,
      data.frame(sample = names(truth), method = "Truth",
                 cluster = as.integer(factor(truth))))
  }
  M <- compare_clusterings(clusters)$agreement
  df <- data.frame(m1 = factor(rep(rownames(M), ncol(M)), levels = rownames(M)),
                   m2 = factor(rep(colnames(M), each = nrow(M)),
                               levels = rev(colnames(M))),
                   ARI = as.vector(M))
  ggplot2::ggplot(df, ggplot2::aes(.data[["m1"]], .data[["m2"]],
                                   fill = .data[["ARI"]])) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data[["ARI"]]),
                                    colour = .data[["ARI"]] > 0.55),
                       size = 3.4, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white",
                                            `FALSE` = "#1A2A4F"),
                                 guide = "none") +
    ggplot2::scale_fill_gradient(low = "#F4EDDE", high = "#1A2A4F",
                                 limits = c(min(0, min(df$ARI)), 1)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, fill = "ARI",
                  title = "Cluster agreement (adjusted Rand index)") +
    theme_factomics() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   panel.grid = ggplot2::element_blank())
}
