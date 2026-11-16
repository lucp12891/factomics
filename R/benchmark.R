# ============================================================
#  factomics — benchmarking / comparison
# ============================================================

#' Compare and benchmark multiple factor analysis results
#'
#' Given the output of [run_factor_pipeline()] (or the metadata-merged variant
#' from [attach_scores_metadata()]), this function:
#' \enumerate{
#'   \item collects scores and loadings for all available methods,
#'   \item optionally z-scores the factors,
#'   \item computes pairwise factor correlations (Hungarian matching),
#'   \item builds correlation heatmaps and scatter plots for scores and
#'         per-omic loadings,
#'   \item builds an all-vs-all score correlation heatmap,
#'   \item (optionally) computes AUC against ground-truth feature/sample labels.
#' }
#'
#' @param results Named list returned by [run_factor_pipeline()] or
#'   [attach_scores_metadata()].
#' @param methods Character vector of method names to compare.  Defaults to
#'   every method in `results`.
#' @param rename_cols Logical.  If `TRUE` (default), factor columns are
#'   labelled `"F1 (METHOD)"`.
#' @param center_scale Logical.  If `TRUE` (default), z-score each factor
#'   column before comparison.
#' @param feature_id_split Character.  Separator used to split
#'   `"view::feature"` row names when aligning loadings (default `"::"`).
#' @param ground_truth_feature Optional named numeric/integer 0/1 vector of
#'   feature labels for AUC calculation.
#' @param ground_truth_sample Optional named numeric/integer 0/1 vector of
#'   sample labels for AUC calculation.
#' @param plot_theme_base_size Numeric.  Base font size for ggplot themes
#'   (default `11`).
#'
#' @return A named list containing:
#' \describe{
#'   \item{`available_methods`}{Character vector of methods actually found.}
#'   \item{`scores`}{Named list of scaled score matrices.}
#'   \item{`loadings`}{Named list of scaled loading matrices.}
#'   \item{`pairwise_scores`}{Pairwise summary stats and correlation matrices
#'         (scores).}
#'   \item{`pairwise_loadings`}{Same for loadings, per dataset.}
#'   \item{`pair_heatmaps_scores`}{Named list of ggplot heatmaps (scores).}
#'   \item{`pair_heatmaps_loadings`}{Named list of ggplot heatmaps (loadings).}
#'   \item{`pair_scatter_scores`}{Scatter plots for matched score pairs.}
#'   \item{`pair_scatter_loadings`}{Scatter plots for matched loading pairs.}
#'   \item{`all_scores_heatmap`}{Global all-method score correlation heatmap.}
#'   \item{`summary_table_scores`}{Combined summary data frame (scores).}
#'   \item{`summary_table_loadings`}{Combined summary data frame (loadings).}
#'   \item{`auc_results`}{AUC table (if ground truth provided).}
#'   \item{`pairs_scores_by_factor`}{Pairs plots per factor (scores).}
#'   \item{`pairs_loadings_by_dataset_and_factor`}{Pairs plots per dataset and
#'         factor (loadings).}
#' }
#'
#' @examples
#' sim   <- simulate_multiomics(n_samples = 40, seed = 1)
#' res   <- run_factor_pipeline(sim$data, methods = c("PCA", "MFA"), k = 3,
#'                              verbose = FALSE)
#' bench <- benchmark_factor_methods(res)
#' bench$summary_table_scores
#' bench$all_scores_heatmap
#'
#' @export
benchmark_factor_methods <- function(
    results,
    methods               = names(results),
    rename_cols           = TRUE,
    center_scale          = TRUE,
    feature_id_split      = "::",
    ground_truth_feature  = NULL,
    ground_truth_sample   = NULL,
    plot_theme_base_size  = 11
) {
  # ---- deps ----
  for (p in c("ggplot2", "reshape2", "clue"))
    .stop_pkg(p)

  # ---- safe helpers ----
  safe_scale <- function(X) {
    if (!center_scale) return(as.matrix(X))
    X  <- as.matrix(X)
    if (ncol(X) == 0L) return(X)
    ok <- colSums(!is.na(X)) > 0
    if (!any(ok)) return(X)
    tmp         <- scale(X[, ok, drop = FALSE])
    X[, ok]     <- as.matrix(tmp)
    X
  }

  rename_factor_cols <- function(X, mname) {
    if (rename_cols) colnames(X) <- paste0("F", seq_len(ncol(X)), " (", mname, ")")
    else if (is.null(colnames(X))) colnames(X) <- paste0("F", seq_len(ncol(X)))
    X
  }

  # ---- identify available methods ----
  available <- intersect(methods, names(results))
  if (length(available) < 2)
    stop("Need at least two methods in 'results'. Found: ",
         paste(available, collapse = ", "), call. = FALSE)

  # ---- collect matrices ----
  get_scores <- function(m) {
    X <- .to_numeric_matrix(results[[m]]$scores, m, id_pref = "sample")
    X <- rename_factor_cols(X, m)
    safe_scale(X)
  }
  get_loadings <- function(m) {
    X <- .to_numeric_matrix(results[[m]]$loadings, m, id_pref = "feature")
    X <- rename_factor_cols(X, m)
    safe_scale(X)
  }
  scores_list   <- lapply(available, get_scores);   names(scores_list)   <- available
  loadings_list <- lapply(available, get_loadings); names(loadings_list) <- available

  # ---- prefix/suffix helpers ----
  get_prefix    <- function(rn) sub("::.*$", "", rn)
  get_suffix    <- function(rn) sub("^.*::", "", rn)
  normalize_ids <- function(xnames) {
    if (is.null(xnames) || is.null(feature_id_split)) return(xnames)
    parts <- strsplit(xnames, split = feature_id_split, fixed = TRUE)
    vapply(parts, function(p) if (length(p)) utils::tail(p, 1) else "", character(1))
  }

  # ---- per-omic loading correlation helper ----
  cor_mat_peromic <- function(A, B) {
    pA <- get_prefix(rownames(A)); pB <- get_prefix(rownames(B))
    omics <- intersect(unique(pA), unique(pB))
    out <- vector("list", length(omics)); names(out) <- omics
    for (om in omics) {
      idxA <- which(pA == om); idxB <- which(pB == om)
      if (!length(idxA) || !length(idxB)) next
      common <- intersect(get_suffix(rownames(A)[idxA]),
                          get_suffix(rownames(B)[idxB]))
      if (!length(common)) next
      AA <- A[idxA[match(common, get_suffix(rownames(A)[idxA]))], , drop = FALSE]
      BB <- B[idxB[match(common, get_suffix(rownames(B)[idxB]))], , drop = FALSE]
      out[[om]] <- stats::cor(AA, BB, use = "pairwise.complete.obs")
    }
    Filter(Negate(is.null), out)
  }

  get_aligned_loading_vecs <- function(A, B, dataset, i, j) {
    pA <- get_prefix(rownames(A)); pB <- get_prefix(rownames(B))
    idxA <- which(pA == dataset); idxB <- which(pB == dataset)
    if (!length(idxA) || !length(idxB)) return(NULL)
    common <- intersect(get_suffix(rownames(A)[idxA]),
                        get_suffix(rownames(B)[idxB]))
    if (!length(common)) return(NULL)
    x <- A[idxA[match(common, get_suffix(rownames(A)[idxA]))], i]
    y <- B[idxB[match(common, get_suffix(rownames(B)[idxB]))], j]
    list(x = x, y = y)
  }

  # ---- scatter helper ----
  make_scatter <- function(x, y, title, xlab, ylab) {
    df  <- data.frame(x = as.numeric(x), y = as.numeric(y))
    ok  <- is.finite(df$x) & is.finite(df$y)
    r   <- suppressWarnings(
      stats::cor(df$x[ok], df$y[ok], use = "pairwise.complete.obs"))
    r_label <- if (is.finite(r)) sprintf("r = %.3f", r) else "r = NA"

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
      ggplot2::geom_point(alpha = 0.6, size = 1, na.rm = TRUE) +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = r_label,
                        hjust = 1.05, vjust = 1.5, size = 4) +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_minimal(base_size = plot_theme_base_size) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

    # Regression line only when clean data is sufficient and coefficients
    # are finite — avoids 'a and b must be finite' crashes with degenerate
    # or NA-containing vectors
    if (sum(ok) >= 3 &&
        stats::var(df$x[ok]) > .Machine$double.eps &&
        stats::var(df$y[ok]) > .Machine$double.eps) {
      fit <- tryCatch(
        stats::lm(y ~ x, data = df[ok, , drop = FALSE]),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        coefs <- stats::coef(fit)
        if (length(coefs) == 2L && all(is.finite(coefs)))
          p <- p + ggplot2::geom_abline(intercept = coefs[1L],
                                         slope     = coefs[2L],
                                         linetype  = "dashed",
                                         linewidth = 0.4)
      }
    }
    p
  }

  # ---- pairwise comparisons ----
  method_pairs             <- t(utils::combn(available, 2))
  pair_summaries_scores    <- list()
  pair_summaries_loadings  <- list()
  pair_heatmaps_scores     <- list()
  pair_heatmaps_loadings   <- list()
  pair_scatter_scores      <- list()
  pair_scatter_loadings    <- list()

  for (r in seq_len(nrow(method_pairs))) {
    m1 <- method_pairs[r, 1]; m2 <- method_pairs[r, 2]
    pair_key <- paste(m1, m2, sep = "|")

    # --- scores ---
    common_samp <- intersect(rownames(scores_list[[m1]]),
                             rownames(scores_list[[m2]]))
    A_s <- scores_list[[m1]][common_samp, , drop = FALSE]
    B_s <- scores_list[[m2]][common_samp, , drop = FALSE]
    Cs  <- stats::cor(A_s, B_s, use = "pairwise.complete.obs")
    ms  <- .match_factors(Cs)

    # score scatters — sign-flip m2 factor when raw correlation is negative
    # (factor signs are arbitrary; |r| is the true match quality measure)
    ss_plots <- list()
    if (!is.null(ms) && nrow(ms) > 0) {
      for (ki in seq_len(nrow(ms))) {
        i    <- ms$A[ki]; j <- ms$B[ki]
        x_sc <- A_s[, i]
        y_sc <- B_s[, j]
        if (is.finite(ms$corr[ki]) && ms$corr[ki] < 0) y_sc <- -y_sc
        ss_plots[[ki]] <- make_scatter(
          x_sc, y_sc,
          title = sprintf("Scores: %s[%s] vs %s[%s]",
                          m1, colnames(A_s)[i], m2, colnames(B_s)[j]),
          xlab  = colnames(A_s)[i], ylab = colnames(B_s)[j]
        )
      }
      # abs-ify stored correlations now that sign-flip has been applied
      ms$corr <- abs(ms$corr)
    }
    pair_scatter_scores[[pair_key]] <- ss_plots

    # summary built AFTER the abs-ification so scores_match$corr is always >= 0
    pair_summaries_scores[[pair_key]] <- list(
      summary = data.frame(
        pair                   = paste(m1, "vs", m2),
        k_m1_scores            = ncol(A_s),
        k_m2_scores            = ncol(B_s),
        k_matched_scores       = if (!is.null(ms)) nrow(ms) else 0L,
        mean_abs_corr_scores   = if (!is.null(ms)) mean(ms$corr,            na.rm = TRUE) else NA_real_,
        median_abs_corr_scores = if (!is.null(ms)) stats::median(ms$corr,   na.rm = TRUE) else NA_real_,
        max_abs_corr_scores    = if (!is.null(ms)) max(ms$corr,              na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      ),
      scores_cor   = Cs,   # raw full cross-corr matrix kept for reference
      scores_match = ms    # corr column is now always >= 0
    )

    dfh  <- reshape2::melt(Cs,
                           varnames = c(paste0("Factor (", m1, ")"),
                                        paste0("Factor (", m2, ")")),
                           value.name = "Correlation")
    xcol <- names(dfh)[2]; ycol <- names(dfh)[1]
    pair_heatmaps_scores[[pair_key]] <-
      ggplot2::ggplot(dfh, ggplot2::aes(x = .data[[xcol]],
                                         y = .data[[ycol]],
                                         fill = .data[["Correlation"]])) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "Corr") +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = paste0("Scores correlation: ", m1, " vs ", m2),
                    x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = plot_theme_base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid  = ggplot2::element_blank(),
                     plot.title  = ggplot2::element_text(face = "bold"))

    # --- loadings (per dataset) ---
    Cl_list <- cor_mat_peromic(loadings_list[[m1]], loadings_list[[m2]])
    ml_list <- lapply(Cl_list, .match_factors)

    rows       <- list()
    scat_by_ds <- list()
    for (om in names(Cl_list)) {
      ml <- ml_list[[om]]

      # loading scatters — sign-flip the m2 factor when raw correlation is negative
      ds_plots <- list()
      if (!is.null(ml) && nrow(ml) > 0) {
        for (ki in seq_len(nrow(ml))) {
          i <- ml$A[ki]; j <- ml$B[ki]
          vecs <- get_aligned_loading_vecs(loadings_list[[m1]],
                                           loadings_list[[m2]], om, i, j)
          if (is.null(vecs)) next
          if (is.finite(ml$corr[ki]) && ml$corr[ki] < 0) vecs$y <- -vecs$y
          ds_plots[[ki]] <- make_scatter(
            vecs$x, vecs$y,
            title = sprintf("Loadings (%s): %s[%s] vs %s[%s]",
                            om, m1, colnames(loadings_list[[m1]])[i],
                            m2, colnames(loadings_list[[m2]])[j]),
            xlab  = colnames(loadings_list[[m1]])[i],
            ylab  = colnames(loadings_list[[m2]])[j]
          )
        }
        # abs-ify now that sign-flip has been applied
        ml_list[[om]]$corr <- abs(ml_list[[om]]$corr)
        ml <- ml_list[[om]]   # refresh local reference
      }
      scat_by_ds[[om]] <- ds_plots

      # summary row built AFTER abs-ification so corr is always >= 0
      rows[[length(rows) + 1]] <- data.frame(
        pair                     = paste(m1, "vs", m2),
        dataset                  = om,
        k_m1_loadings            = ncol(loadings_list[[m1]]),
        k_m2_loadings            = ncol(loadings_list[[m2]]),
        k_matched_loadings       = if (!is.null(ml)) nrow(ml) else 0L,
        mean_abs_corr_loadings   = if (!is.null(ml)) mean(ml$corr,            na.rm = TRUE) else NA_real_,
        median_abs_corr_loadings = if (!is.null(ml)) stats::median(ml$corr,   na.rm = TRUE) else NA_real_,
        max_abs_corr_loadings    = if (!is.null(ml)) max(ml$corr,              na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )

      dfhL <- reshape2::melt(Cl_list[[om]],
                             varnames = c(paste0("Factor (", m1, ")"),
                                          paste0("Factor (", m2, ")")),
                             value.name = "Correlation")
      xL <- names(dfhL)[2]; yL <- names(dfhL)[1]
      pair_heatmaps_loadings[[paste(m1, m2, om, sep = "|")]] <-
        ggplot2::ggplot(dfhL,
                        ggplot2::aes(x = .data[[xL]], y = .data[[yL]],
                                     fill = .data[["Correlation"]])) +
        ggplot2::geom_tile(color = "white", linewidth = 0.2) +
        ggplot2::scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "Corr") +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = paste0("Loadings (", om, "): ", m1, " vs ", m2),
                      x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = plot_theme_base_size) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       panel.grid  = ggplot2::element_blank(),
                       plot.title  = ggplot2::element_text(face = "bold"))
    }

    pair_summaries_loadings[[pair_key]] <- list(
      cor     = Cl_list,
      match   = ml_list,
      summary = if (length(rows)) do.call(rbind, rows) else NULL
    )
    pair_scatter_loadings[[pair_key]] <- scat_by_ds
  }

  # helper: remove columns with zero (or near-zero) variance
  .drop_degenerate_cols <- function(df) {
    keep <- vapply(df, function(v) {
      ok <- is.finite(v)
      sum(ok) >= 3 && stats::var(v[ok]) > .Machine$double.eps
    }, logical(1))
    df[, keep, drop = FALSE]
  }

  # ---- pairs-style multi-method plots ----
  ref_method <- if ("FABIA" %in% available) "FABIA" else available[1]

  maps_scores        <- .align_to_reference_scores(scores_list, ref_method,
                                                    .match_factors)
  common_samples_all <- Reduce(intersect, lapply(scores_list, rownames))
  ref_k              <- ncol(scores_list[[ref_method]])

  pairs_scores_by_factor <- tryCatch({
    out <- list()
    for (i in seq_len(ref_k)) {
      cols    <- list()
      ref_vec <- scores_list[[ref_method]][common_samples_all, i]
      for (m in available) {
        if (m == ref_method) {
          cols[[paste0("Score (", m, ")")]] <- ref_vec
        } else {
          j <- maps_scores[[m]][i]
          if (is.na(j)) next
          v <- scores_list[[m]][common_samples_all, j]
          v <- .align_sign_vec(ref_vec, v)
          cols[[paste0("Score (", m, ")")]] <- v
        }
      }
      if (length(cols) >= 2) {
        mat <- .drop_degenerate_cols(as.data.frame(cols, check.names = FALSE))
        if (ncol(mat) >= 2)
          out[[paste0("F", i)]] <-
            .make_pairs_plot(mat, main = "Factor scores (Pearson cor.)")
      }
    }
    out
  }, error = function(e) {
    message("[factomics] pairs_scores_by_factor skipped: ", conditionMessage(e))
    list()
  })

  maps_loadings    <- .align_to_reference_loadings(
    loadings_list, ref_method, get_prefix, get_suffix, .match_factors)
  datasets_ref     <- unique(get_prefix(rownames(loadings_list[[ref_method]])))

  pairs_loadings_by_dataset_and_factor <- tryCatch({
    out <- list()
    for (om in datasets_ref) {
      pref_ref   <- get_prefix(rownames(loadings_list[[ref_method]]))
      idx_ref    <- which(pref_ref == om)
      ids_ref    <- get_suffix(rownames(loadings_list[[ref_method]])[idx_ref])
      common_ids <- ids_ref
      for (m in available) {
        if (m == ref_method) next
        pref_m     <- get_prefix(rownames(loadings_list[[m]]))
        ids_m      <- get_suffix(rownames(loadings_list[[m]])[pref_m == om])
        common_ids <- intersect(common_ids, ids_m)
      }
      if (!length(common_ids)) next

      row_idx <- list()
      row_idx[[ref_method]] <- idx_ref[match(common_ids, ids_ref)]
      for (m in available) {
        if (m == ref_method) next
        pref_m        <- get_prefix(rownames(loadings_list[[m]]))
        idx_m         <- which(pref_m == om)
        ids_m         <- get_suffix(rownames(loadings_list[[m]])[idx_m])
        row_idx[[m]]  <- idx_m[match(common_ids, ids_m)]
      }

      for (i in seq_len(ref_k)) {
        cols    <- list()
        ref_vec <- loadings_list[[ref_method]][row_idx[[ref_method]], i]
        for (m in available) {
          if (m == ref_method) {
            cols[[m]] <- ref_vec
          } else {
            j <- maps_loadings[[om]][[m]][i]
            if (is.na(j)) next
            v <- loadings_list[[m]][row_idx[[m]], j]
            v <- .align_sign_vec(ref_vec, v)
            cols[[m]] <- v
          }
        }
        if (length(cols) >= 2) {
          mat <- .drop_degenerate_cols(as.data.frame(cols, check.names = FALSE))
          if (ncol(mat) >= 2)
            out[[om]][[paste0("F", i)]] <-
              .make_pairs_plot(mat,
                               main = paste0("Loadings (", om, ") - Pearson cor."))
        }
      }
    }
    out
  }, error = function(e) {
    message("[factomics] pairs_loadings_by_dataset_and_factor skipped: ",
            conditionMessage(e))
    list()
  })

  # ---- global all-in-one heatmap ----
  merge_scores <- function(lst) {
    X <- lst[[1]]; colnames(X) <- make.unique(colnames(X))
    for (i in seq(2, length(lst))) {
      Y       <- lst[[i]]; colnames(Y) <- make.unique(colnames(Y))
      common  <- intersect(rownames(X), rownames(Y))
      X       <- cbind(X[common, , drop = FALSE], Y[common, , drop = FALSE])
    }
    X
  }
  all_scores_mat <- merge_scores(scores_list)
  all_cor        <- stats::cor(all_scores_mat, use = "pairwise.complete.obs")

  get_method <- function(nm) sub(".*\\(([^)]+)\\)\\s*$", "\\1", nm)
  get_fnum   <- function(nm) suppressWarnings(
    as.integer(sub("^\\s*F(\\d+).*", "\\1", nm)))
  mth  <- vapply(colnames(all_cor), get_method, character(1))
  fno  <- vapply(colnames(all_cor), get_fnum, numeric(1))
  lev  <- unique(mth)
  ord  <- order(match(mth, lev), fno, colnames(all_cor))
  all_cor <- all_cor[ord, ord, drop = FALSE]

  dfall         <- reshape2::melt(all_cor,
                                  varnames = c("Factor1", "Factor2"),
                                  value.name = "Correlation")
  dfall$Factor1 <- factor(dfall$Factor1, levels = rownames(all_cor))
  dfall$Factor2 <- factor(dfall$Factor2, levels = colnames(all_cor))

  tab    <- as.data.frame(table(mth[ord]), stringsAsFactors = FALSE)
  tab$cum <- cumsum(tab$Freq)
  bounds <- tab$cum + 0.5

  p_all <- ggplot2::ggplot(dfall,
                            ggplot2::aes(x = .data[["Factor2"]],
                                         y = .data[["Factor1"]],
                                         fill = .data[["Correlation"]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "Corr") +
    ggplot2::geom_vline(xintercept = bounds, linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = bounds, linewidth = 0.3) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "All methods - factor correlation (scores)",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = plot_theme_base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid  = ggplot2::element_blank(),
                   plot.title  = ggplot2::element_text(face = "bold"))

  # ---- optional AUC ----
  auc_results <- NULL
  if (!is.null(ground_truth_feature) || !is.null(ground_truth_sample)) {
    .stop_pkg("pROC")
    calc_auc <- function(scores, truth) {
      truth <- as.numeric(truth)
      if (!length(scores) || !length(truth) ||
          length(scores) != length(truth)) return(NA_real_)
      roc <- pROC::roc(response = truth, predictor = scores,
                       quiet = TRUE, direction = "<")
      as.numeric(pROC::auc(roc))
    }
    rows_auc <- list()
    if (!is.null(ground_truth_feature)) {
      for (m in available) {
        L   <- loadings_list[[m]]
        s   <- apply(abs(L), 1, max, na.rm = TRUE)
        idsL <- normalize_ids(names(s))
        idsT <- names(ground_truth_feature)
        if (!is.null(idsT)) idsT <- normalize_ids(idsT)
        common <- intersect(idsL, idsT)
        if (length(common))
          rows_auc[[length(rows_auc) + 1]] <- data.frame(
            target = "feature", method = m,
            AUC    = calc_auc(s[match(common, idsL)],
                              ground_truth_feature[match(common, idsT)]),
            stringsAsFactors = FALSE
          )
      }
    }
    if (!is.null(ground_truth_sample)) {
      for (m in available) {
        S      <- scores_list[[m]]
        ss     <- apply(abs(S), 1, max, na.rm = TRUE)
        common <- intersect(names(ss), names(ground_truth_sample))
        if (length(common))
          rows_auc[[length(rows_auc) + 1]] <- data.frame(
            target = "sample", method = m,
            AUC    = calc_auc(ss[common], ground_truth_sample[common]),
            stringsAsFactors = FALSE
          )
      }
    }
    auc_results <- if (length(rows_auc)) do.call(rbind, rows_auc) else NULL
  }

  # ---- assemble summaries ----
  summary_table_scores <- do.call(
    rbind,
    lapply(pair_summaries_scores, `[[`, "summary")
  )
  rownames(summary_table_scores) <- NULL

  summary_table_loadings <- do.call(
    rbind,
    lapply(pair_summaries_loadings, function(x)
      if (!is.null(x$summary)) x$summary else NULL)
  )
  rownames(summary_table_loadings) <- NULL

  list(
    available_methods                     = available,
    scores                                = scores_list,
    loadings                              = loadings_list,
    pairwise_scores                       = pair_summaries_scores,
    pairwise_loadings                     = pair_summaries_loadings,
    pair_heatmaps_scores                  = pair_heatmaps_scores,
    pair_heatmaps_loadings                = pair_heatmaps_loadings,
    pair_scatter_scores                   = pair_scatter_scores,
    pair_scatter_loadings                 = pair_scatter_loadings,
    all_scores_heatmap                    = p_all,
    summary_table_scores                  = summary_table_scores,
    summary_table_loadings                = summary_table_loadings,
    auc_results                           = auc_results,
    pairs_scores_by_factor                = pairs_scores_by_factor,
    pairs_loadings_by_dataset_and_factor  = pairs_loadings_by_dataset_and_factor
  )
}
