# ============================================================
#  factomics — factor visualisation and interpretation
# ============================================================

#' The factomics ggplot2 theme and colour palette
#'
#' @param base_size Base font size.
#' @param n Number of colours.
#' @return `theme_factomics()` returns a ggplot2 theme;
#'   `factomics_palette()` returns a character vector of colours.
#' @export
theme_factomics <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", colour = "#1A2A4F"),
      strip.text       = ggplot2::element_text(face = "bold", colour = "#1A2A4F"),
      strip.background = ggplot2::element_rect(fill = "#F4EDDE", colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(face = "bold")
    )
}

#' @rdname theme_factomics
#' @export
factomics_palette <- function(n) {
  base <- c("#1A2A4F", "#E8366F", "#12B886", "#F59F00", "#4C6EF5",
            "#AE3EC9", "#868E96", "#2B8A3E", "#C92A2A", "#1098AD")
  if (n <= length(base)) base[seq_len(n)] else
    grDevices::colorRampPalette(base)(n)
}

.prep_metadata <- function(metadata, sample_col = "sample") {
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!(sample_col %in% names(meta))) {
    if (!.has_default_rownames(meta)) meta[[sample_col]] <- rownames(meta)
    else stop("`metadata` must have column '", sample_col,
              "' or non-default rownames as sample IDs.", call. = FALSE)
  }
  meta$sample <- as.character(meta[[sample_col]])
  meta[!duplicated(meta$sample), , drop = FALSE]
}

.get_views <- function(results, data = NULL) {
  if (!is.null(data)) return(.prepare_views(data, scale_features = TRUE))
  V <- attr(results, "views")
  if (is.null(V))
    stop("The input data are not stored in `results`; pass them with `data =`.",
         call. = FALSE)
  V
}


#' Variance explained by each factor in each view
#'
#' Computes, for every method, the fraction of variance (R^2) of each
#' (centred) view explained by each factor, by regressing the view on the
#' factor scores. This works for any method, so methods can be compared on
#' the same scale.
#'
#' @param results Object returned by [run_factor_pipeline()].
#' @param methods Methods to include (default: all).
#' @param data Optional input data, if `results` no longer carries it.
#'
#' @return A data frame with columns `method`, `view`, `factor` (`"F1"`, ...,
#'   or `"Total"` for all factors jointly) and `R2`.
#' @export
variance_explained <- function(results, methods = names(results),
                               data = NULL) {
  V <- .get_views(results, data)
  out <- list()
  for (m in methods) {
    Z <- as.matrix(results[[m]]$scores)
    Z <- scale(Z[colnames(V[[1]]), , drop = FALSE], scale = FALSE)
    for (v in names(V)) {
      Xv  <- scale(t(V[[v]]), scale = FALSE)
      tss <- sum(Xv^2)
      r2f <- vapply(seq_len(ncol(Z)), function(j) {
        z <- Z[, j]
        if (sum(z^2) == 0) return(0)
        b <- crossprod(z, Xv) / sum(z^2)
        1 - sum((Xv - z %*% b)^2) / tss
      }, numeric(1))
      fit <- stats::lm.fit(Z, Xv)
      r2t <- 1 - sum(fit$residuals^2) / tss
      out[[length(out) + 1]] <- data.frame(
        method = m, view = v,
        factor = c(colnames(Z), "Total"),
        R2 = c(r2f, r2t), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}


#' Plot variance explained per factor and view
#'
#' @inheritParams variance_explained
#' @param type `"heatmap"` (factor x view R^2 per method) or `"total"`
#'   (total R^2 per view as bars).
#' @return A ggplot object.
#' @export
plot_variance_explained <- function(results, methods = names(results),
                                    type = c("heatmap", "total"),
                                    data = NULL) {
  type <- match.arg(type)
  ve <- variance_explained(results, methods, data)
  ve$method <- factor(ve$method, levels = methods)
  if (type == "total") {
    df <- ve[ve$factor == "Total", ]
    return(
      ggplot2::ggplot(df, ggplot2::aes(.data[["view"]], 100 * .data[["R2"]],
                                       fill = .data[["method"]])) +
        ggplot2::geom_col(position = ggplot2::position_dodge(0.8),
                          width = 0.75) +
        ggplot2::scale_fill_manual(values = factomics_palette(length(methods))) +
        ggplot2::labs(x = NULL, y = "Variance explained (%)", fill = "Method",
                      title = "Total variance explained per view") +
        theme_factomics()
    )
  }
  df <- ve[ve$factor != "Total", ]
  df$factor <- factor(df$factor, levels = rev(unique(df$factor)))
  ggplot2::ggplot(df, ggplot2::aes(.data[["view"]], .data[["factor"]],
                                   fill = 100 * .data[["R2"]])) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f", 100 * .data[["R2"]])),
                       size = 3) +
    ggplot2::scale_fill_gradient(low = "#F4EDDE", high = "#C92A64",
                                 name = "R\u00b2 (%)") +
    ggplot2::facet_wrap(ggplot2::vars(.data[["method"]]), nrow = 1) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Variance explained by each factor") +
    theme_factomics() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   panel.grid = ggplot2::element_blank())
}


#' Scatter plot of two factors for one or more methods
#'
#' @param results Object returned by [run_factor_pipeline()].
#' @param methods Methods to plot (default: all).
#' @param factors Integer vector of length 2: which factors to plot.
#' @param metadata Optional sample annotations.
#' @param color_by Column of `metadata` used for colour.
#' @param sample_col Sample-ID column of `metadata`.
#' @return A ggplot object.
#' @export
plot_factor_scores <- function(results, methods = names(results),
                               factors = c(1, 2), metadata = NULL,
                               color_by = NULL, sample_col = "sample") {
  df <- do.call(rbind, lapply(methods, function(m) {
    Z <- as.matrix(results[[m]]$scores)
    if (ncol(Z) < max(factors)) return(NULL)
    data.frame(sample = rownames(Z), method = m,
               Dim1 = Z[, factors[1]], Dim2 = Z[, factors[2]],
               stringsAsFactors = FALSE)
  }))
  if (is.null(df)) stop("No method has enough factors.", call. = FALSE)
  df$method <- factor(df$method, levels = methods)
  plot_embedding(df, metadata = metadata, color_by = color_by,
                 sample_col = sample_col) +
    ggplot2::labs(x = paste0("Factor ", factors[1]),
                  y = paste0("Factor ", factors[2]))
}


#' Plot the top-weighted features of a factor
#'
#' @param results Object returned by [run_factor_pipeline()].
#' @param method Method name.
#' @param factor Factor index (default `1`).
#' @param n Number of features to show (default `20`).
#' @param view Optional view name to restrict to.
#' @return A ggplot object (bars coloured by view).
#' @export
plot_top_loadings <- function(results, method, factor = 1, n = 20,
                              view = NULL) {
  L <- as.matrix(results[[method]]$loadings)
  if (factor > ncol(L)) stop("Factor ", factor, " not available for ", method,
                             call. = FALSE)
  df <- data.frame(id = rownames(L), weight = L[, factor],
                   view = sub("::.*$", "", rownames(L)),
                   feature = sub("^.*::", "", rownames(L)),
                   stringsAsFactors = FALSE)
  if (!is.null(view)) df <- df[df$view %in% view, ]
  df <- utils::head(df[order(-abs(df$weight)), ], n)
  df$feature <- factor(df$id, levels = rev(df$id),
                       labels = rev(df$feature))
  ggplot2::ggplot(df, ggplot2::aes(.data[["weight"]], .data[["feature"]],
                                   fill = .data[["view"]])) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::scale_fill_manual(values = factomics_palette(
      length(unique(df$view)))) +
    ggplot2::labs(x = "Loading", y = NULL, fill = "View",
                  title = sprintf("%s - top features of F%d", method, factor)) +
    theme_factomics()
}


#' Heatmap of factor scores, ordered by a sample annotation
#'
#' @param results Object returned by [run_factor_pipeline()].
#' @param method Method name.
#' @param metadata Optional sample annotations.
#' @param order_by Optional metadata column used to order (and label) samples.
#' @param sample_col Sample-ID column of `metadata`.
#' @return A ggplot object.
#' @export
plot_factor_heatmap <- function(results, method, metadata = NULL,
                                order_by = NULL, sample_col = "sample") {
  Z <- scale(as.matrix(results[[method]]$scores))
  ord <- rownames(Z)
  grp <- NULL
  if (!is.null(metadata) && !is.null(order_by)) {
    meta <- .prep_metadata(metadata, sample_col)
    grp  <- meta[[order_by]][match(rownames(Z), meta$sample)]
    ord  <- rownames(Z)[order(grp)]
  }
  df <- data.frame(sample = factor(rep(rownames(Z), ncol(Z)), levels = ord),
                   factor = factor(rep(colnames(Z), each = nrow(Z)),
                                   levels = rev(colnames(Z))),
                   z = as.vector(Z))
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[["sample"]], .data[["factor"]],
                                        fill = .data[["z"]])) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#1A2A4F", mid = "white",
                                  high = "#E8366F", name = "z-score") +
    ggplot2::labs(x = "Samples", y = NULL,
                  title = paste0(method, " - factor scores")) +
    theme_factomics() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  if (!is.null(grp)) {
    bar <- data.frame(sample = factor(rownames(Z), levels = ord),
                      group = as.character(grp))
    p <- p + ggplot2::geom_tile(
      data = bar,
      ggplot2::aes(x = .data[["sample"]], y = ncol(Z) + 0.9,
                   colour = .data[["group"]]),
      inherit.aes = FALSE, height = 0.5, fill = NA, linewidth = 2) +
      ggplot2::scale_colour_manual(
        values = factomics_palette(length(unique(bar$group))),
        name = order_by)
  }
  p
}


#' Test association between factors and sample covariates
#'
#' For every method, factor and covariate, tests whether the factor scores
#' differ by the covariate: Kruskal-Wallis for categorical covariates,
#' Spearman correlation for numeric ones.
#'
#' @param results Object returned by [run_factor_pipeline()].
#' @param metadata Sample annotations.
#' @param covariates Metadata columns to test (default: all except the sample
#'   column).
#' @param methods Methods to include (default: all).
#' @param sample_col Sample-ID column of `metadata`.
#' @return A data frame with `method`, `factor`, `covariate`, `test`,
#'   `statistic`, `p_value` and BH-adjusted `p_adj`.
#' @export
factor_associations <- function(results, metadata, covariates = NULL,
                                methods = names(results),
                                sample_col = "sample") {
  meta <- .prep_metadata(metadata, sample_col)
  if (is.null(covariates))
    covariates <- setdiff(names(meta), c("sample", sample_col))
  out <- list()
  for (m in methods) {
    Z <- as.matrix(results[[m]]$scores)
    mm <- meta[match(rownames(Z), meta$sample), , drop = FALSE]
    for (cv in covariates) for (j in seq_len(ncol(Z))) {
      y <- mm[[cv]]; z <- Z[, j]
      ok <- !is.na(y)
      if (sum(ok) < 3) next
      if (is.numeric(y)) {
        tt <- suppressWarnings(stats::cor.test(z[ok], y[ok],
                                               method = "spearman"))
        test <- "spearman"
      } else {
        if (length(unique(y[ok])) < 2) next
        tt <- stats::kruskal.test(z[ok], factor(y[ok]))
        test <- "kruskal"
      }
      out[[length(out) + 1]] <- data.frame(
        method = m, factor = colnames(Z)[j], covariate = cv, test = test,
        statistic = unname(tt$statistic), p_value = tt$p.value,
        stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out)
  res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
  res
}


#' Plot factor-covariate associations
#'
#' @inheritParams factor_associations
#' @return A ggplot heatmap of -log10 adjusted p-values.
#' @export
plot_factor_associations <- function(results, metadata, covariates = NULL,
                                     methods = names(results),
                                     sample_col = "sample") {
  fa <- factor_associations(results, metadata, covariates, methods, sample_col)
  fa$method <- factor(fa$method, levels = methods)
  fa$score  <- -log10(pmax(fa$p_adj, 1e-300))
  fa$factor <- factor(fa$factor, levels = rev(unique(fa$factor)))
  ggplot2::ggplot(fa, ggplot2::aes(.data[["covariate"]], .data[["factor"]],
                                   fill = .data[["score"]])) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data[["p_adj"]] < 0.05,
                                                   "*", "")), size = 5) +
    ggplot2::scale_fill_gradient(low = "#F4EDDE", high = "#1A2A4F",
                                 name = "-log10(p adj)") +
    ggplot2::facet_wrap(ggplot2::vars(.data[["method"]]), nrow = 1) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Factor-covariate association (* p adj < 0.05)") +
    theme_factomics() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   panel.grid = ggplot2::element_blank())
}
