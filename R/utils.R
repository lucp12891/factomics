# ============================================================
#  factomics — internal utilities
#  These functions are not exported (no @export tag).
# ============================================================

# ---- package dependency check ----
.stop_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is required but not installed. ",
         "Install it with: install.packages('", p, "')", call. = FALSE)
}

# ---- view helpers ----
.stack_views <- function(lst) {
  do.call(rbind, lapply(names(lst), function(v) {
    M  <- lst[[v]]
    rn <- rownames(M)
    if (is.null(rn)) rn <- paste0("f", seq_len(nrow(M)))
    rownames(M) <- paste0(v, "::", rn)
    M
  }))
}

.as_list_views <- function(X) {
  if (is.list(X)) unclass(X) else list(omic1 = X)
}

# ---- sign alignment ----
.align_sign_vec <- function(ref, v) {
  r <- suppressWarnings(stats::cor(ref, v, use = "pairwise.complete.obs"))
  if (is.finite(r) && r < 0) -v else v
}

# ---- pairs plot — pure ggplot2, no base graphics device risk ----
#
# Replaces the old base-R pairs() approach entirely.
# Returns a patchwork object (or NULL if < 2 non-degenerate columns).

.gg_pairs_plot <- function(mat, main = "") {

  # ---- keep only columns with enough finite variance ----
  ok_col <- vapply(as.data.frame(mat), function(v) {
    fin <- v[is.finite(v)]
    length(fin) >= 3L && stats::var(fin) > .Machine$double.eps
  }, logical(1L))

  mat <- mat[, ok_col, drop = FALSE]
  if (ncol(mat) < 2L) {
    warning("[factomics] pairs plot skipped: fewer than 2 non-degenerate columns.",
            call. = FALSE)
    return(NULL)
  }

  nms  <- colnames(mat)
  n    <- length(nms)
  df   <- as.data.frame(mat, check.names = FALSE)

  # ---- build one ggplot per (row, col) cell ----
  plots <- vector("list", n * n)

  for (row in seq_len(n)) {
    for (col in seq_len(n)) {
      idx <- (row - 1L) * n + col
      nm_x <- nms[col]; nm_y <- nms[row]
      x <- df[[nm_x]]; y <- df[[nm_y]]
      ok <- is.finite(x) & is.finite(y)

      if (row == col) {
        # ---- diagonal: density histogram ----
        p <- ggplot2::ggplot(df[ok, , drop = FALSE],
                              ggplot2::aes(x = .data[[nm_x]])) +
          ggplot2::geom_histogram(bins = 15, fill = "grey75",
                                   colour = "white", linewidth = 0.2) +
          ggplot2::annotate("text", x = Inf, y = Inf,
                            label = nm_x, hjust = 1.1, vjust = 1.3,
                            size = 3.2, fontface = "bold") +
          ggplot2::theme_classic(base_size = 9) +
          ggplot2::theme(axis.title = ggplot2::element_blank(),
                         axis.text  = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank())

      } else if (row < col) {
        # ---- upper triangle: Pearson r ----
        r <- suppressWarnings(
          stats::cor(x[ok], y[ok], use = "pairwise.complete.obs"))
        lbl <- if (is.finite(r)) sprintf("%.2f", abs(r)) else "NA"
        sz  <- if (is.finite(r)) 3 + 4 * abs(r) else 3

        p <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5, label = lbl,
                            size = sz, fontface = "bold") +
          ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
          ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
          ggplot2::theme_void()

      } else {
        # ---- lower triangle: scatter + optional regression line ----
        dfc <- data.frame(x = x[ok], y = y[ok])
        p   <- ggplot2::ggplot(dfc, ggplot2::aes(x = .data[["x"]],
                                                   y = .data[["y"]])) +
          ggplot2::geom_point(size = 0.7, alpha = 0.6, shape = 8) +
          ggplot2::theme_classic(base_size = 9) +
          ggplot2::theme(axis.title = ggplot2::element_blank(),
                         axis.text  = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank())

        # Regression line: fit on clean data, check coefs before adding
        if (nrow(dfc) >= 3L &&
            stats::var(dfc$x) > .Machine$double.eps &&
            stats::var(dfc$y) > .Machine$double.eps) {
          fit <- tryCatch(stats::lm(y ~ x, data = dfc), error = function(e) NULL)
          if (!is.null(fit)) {
            cf <- stats::coef(fit)
            if (length(cf) == 2L && all(is.finite(cf)))
              p <- p + ggplot2::geom_abline(intercept = cf[1L], slope = cf[2L],
                                             colour = "red", linewidth = 0.5,
                                             linetype = "dashed")
          }
        }
      }
      plots[[idx]] <- p
    }
  }

  # ---- assemble with patchwork ----
  pw <- patchwork::wrap_plots(plots, nrow = n, ncol = n) +
    patchwork::plot_annotation(title = main,
                               theme = ggplot2::theme(
                                 plot.title = ggplot2::element_text(
                                   face = "bold", size = 11, hjust = 0.5)))
  pw
}

# Keep old name as alias so all call-sites stay unchanged
.make_pairs_plot <- .gg_pairs_plot

# ---- Hungarian factor alignment (requires clue) ----
.match_factors <- function(C) {
  if (is.null(C) || any(dim(C) == 0)) return(NULL)
  C[!is.finite(C)] <- NA_real_
  absC <- abs(C); absC[is.na(absC)] <- 0
  ka   <- nrow(absC); kb <- ncol(absC)
  maxv <- if (all(absC == 0)) 0 else max(absC, na.rm = TRUE)
  cost <- maxv - absC
  pad_val <- suppressWarnings(max(cost, na.rm = TRUE))
  pad_val <- if (is.finite(pad_val)) pad_val + 1 else 1
  if (ka > kb) cost <- cbind(cost, matrix(pad_val, nrow = ka, ncol = ka - kb))
  if (kb > ka) cost <- rbind(cost, matrix(pad_val, nrow = kb - ka, ncol = kb))
  storage.mode(cost) <- "double"
  ass       <- clue::solve_LSAP(cost)
  real_rows <- seq_len(ka)
  j_all     <- as.integer(ass[real_rows])
  keep      <- which(j_all <= kb)
  if (!length(keep)) return(NULL)
  i_keep <- real_rows[keep]; j_keep <- j_all[keep]
  data.frame(
    A    = i_keep,
    B    = j_keep,
    corr = mapply(function(i, j) C[i, j], i_keep, j_keep),
    stringsAsFactors = FALSE
  )
}

# ---- reference-based alignment helpers ----
.align_to_reference_scores <- function(scores_list, ref, match_factors_fun) {
  ref_mat <- scores_list[[ref]]
  ref_k   <- ncol(ref_mat)
  maps    <- list()
  for (m in names(scores_list)) {
    if (m == ref) next
    common <- intersect(rownames(ref_mat), rownames(scores_list[[m]]))
    A <- ref_mat[common, , drop = FALSE]
    B <- scores_list[[m]][common, , drop = FALSE]
    C <- stats::cor(A, B, use = "pairwise.complete.obs")
    ms  <- match_factors_fun(C)
    map <- rep(NA_integer_, ref_k)
    if (!is.null(ms) && nrow(ms) > 0) map[ms$A] <- ms$B
    maps[[m]] <- map
  }
  maps
}

.align_to_reference_loadings <- function(loadings_list, ref,
                                          get_prefix, get_suffix,
                                          match_factors_fun) {
  ref_mat  <- loadings_list[[ref]]
  ref_k    <- ncol(ref_mat)
  pref_ref <- get_prefix(rownames(ref_mat))
  ds       <- unique(pref_ref)
  out      <- stats::setNames(vector("list", length(ds)), ds)

  for (om in ds) {
    idx_ref <- which(pref_ref == om)
    ids_ref <- get_suffix(rownames(ref_mat)[idx_ref])
    maps_om <- list()
    for (m in names(loadings_list)) {
      if (m == ref) next
      pref_m <- get_prefix(rownames(loadings_list[[m]]))
      idx_m  <- which(pref_m == om)
      if (!length(idx_m)) { maps_om[[m]] <- rep(NA_integer_, ref_k); next }
      ids_m  <- get_suffix(rownames(loadings_list[[m]])[idx_m])
      common <- intersect(ids_ref, ids_m)
      if (!length(common)) { maps_om[[m]] <- rep(NA_integer_, ref_k); next }
      A <- ref_mat[idx_ref[match(common, ids_ref)], , drop = FALSE]
      B <- loadings_list[[m]][idx_m[match(common, ids_m)], , drop = FALSE]
      C <- stats::cor(A, B, use = "pairwise.complete.obs")
      ms  <- match_factors_fun(C)
      map <- rep(NA_integer_, ref_k)
      if (!is.null(ms) && nrow(ms) > 0) map[ms$A] <- ms$B
      maps_om[[m]] <- map
    }
    out[[om]] <- maps_om
  }
  out
}

# ---- palette builder ----
.build_palette <- function(levels, pal = NULL, na_color = "#999999") {
  levels <- unique(as.character(levels))
  if (!"NA" %in% levels) levels <- c(levels, "NA")
  non_na <- setdiff(levels, "NA")
  auto   <- if (length(non_na) > 0)
    stats::setNames(scales::hue_pal()(length(non_na)), non_na)
  else character(0)
  out <- c(auto, "NA" = na_color)
  if (!is.null(pal)) {
    if (is.list(pal)) pal <- unlist(pal, use.names = TRUE)
    nms <- names(pal)          # save names BEFORE as.character() (safe for all R versions)
    pal <- as.character(pal)
    names(pal) <- nms
    if (is.null(names(pal)) || all(nzchar(names(pal)) == FALSE))
      stop("Custom palette must be a *named* vector (names = group levels).", call. = FALSE)
    out[names(pal)] <- pal
  }
  out[levels]
}

# ---- coerce to numeric matrix (from data.frame with possible id cols) ----
.to_numeric_matrix <- function(obj, method_name,
                                id_pref = c("sample", "feature")) {
  df <- as.data.frame(obj, check.names = FALSE, stringsAsFactors = FALSE)

  # Prefer an explicit id column over rownames (default integer rownames
  # "1","2",... are not informative and should not be used as sample IDs).
  hit <- intersect(id_pref, names(df))
  if (length(hit)) {
    rn <- as.character(df[[hit[1]]])
  } else {
    rn <- rownames(df)
    if (is.null(rn) || anyNA(rn) || any(rn == "")) rn <- NULL
  }

  # prefer columns already labelled "F1 (METHOD)"
  f_like <- grepl("^F[0-9]+\\s*\\(", names(df))
  keep   <- f_like & vapply(df, is.numeric, logical(1))
  if (!any(keep)) {
    drop_names <- c("sample", "feature",
                    grep("^signal_", names(df), value = TRUE))
    keep <- vapply(df, is.numeric, logical(1)) & !names(df) %in% drop_names
  }
  if (!any(keep))
    stop("No numeric columns found for method ", method_name,
         ". Columns: ", paste(names(df), collapse = ", "), call. = FALSE)
  mat <- as.matrix(df[, keep, drop = FALSE])
  storage.mode(mat) <- "double"
  if (!is.null(rn)) rownames(mat) <- rn
  mat
}

# ---- check requested methods exist in a results object ----
.check_methods <- function(results, methods) {
  miss <- setdiff(methods, names(results))
  if (length(miss))
    stop("Method(s) not found in results: ", paste(miss, collapse = ", "),
         ". Available: ", paste(names(results), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}
