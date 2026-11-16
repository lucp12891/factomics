# ============================================================
#  factomics — visualisations
# ============================================================

# ---- internal: resolve factor column name in a data frame ----
.resolve_factor_col <- function(df, method, factor_idx = 1) {
  cand  <- paste0("F", factor_idx, " (", method, ")")
  short <- paste0("F", factor_idx)
  if (cand  %in% names(df)) return(cand)
  if (short %in% names(df)) return(short)
  hit <- grep(paste0("^F", factor_idx, "(\\b|\\s*\\()"), names(df), value = TRUE)
  if (length(hit)) hit[1] else stop("Could not find F", factor_idx,
                                     " for ", method, call. = FALSE)
}

# ---- internal: find which factor of target_method best matches
#      factor factor_idx of ref_method (by max |Pearson r|) ----
#
# Returns the matched integer factor index AND the sign (+1/-1) needed
# so that the matched factor is positively correlated with the reference.
.find_matched_factor <- function(bench, ref_method, factor_idx, target_method) {
  ref_df  <- bench[[ref_method]][["scores_with_meta"]]
  tgt_df  <- bench[[target_method]][["scores_with_meta"]]
  if (is.null(ref_df) || is.null(tgt_df))
    return(list(idx = factor_idx, sign = 1L))

  ref_col <- .resolve_factor_col(ref_df, ref_method, factor_idx)
  ref_vec <- as.numeric(ref_df[[ref_col]])

  # all factor columns in the target (e.g. "F1", "F2", ... or "F1 (GFA)")
  all_cols   <- names(tgt_df)
  fcols      <- grep("^F[0-9]+(\\s*\\(|$)", all_cols, value = TRUE)
  if (!length(fcols)) return(list(idx = factor_idx, sign = 1L))

  # correlation of each target factor with the reference factor
  cors <- vapply(fcols, function(col) {
    v      <- as.numeric(tgt_df[[col]])
    common <- stats::complete.cases(ref_vec, v)
    if (sum(common) < 3L) return(0)
    stats::cor(ref_vec[common], v[common])
  }, numeric(1))

  best_pos <- which.max(abs(cors))
  best_col <- fcols[best_pos]

  # extract numeric index from the column name
  best_idx  <- suppressWarnings(
    as.integer(sub("^F([0-9]+).*", "\\1", best_col)))
  if (is.na(best_idx)) best_idx <- factor_idx

  best_sign <- if (is.finite(cors[best_pos]) && cors[best_pos] < 0) -1L else 1L
  list(idx = best_idx, sign = best_sign)
}

# ---- internal: build tidy data frame for one method ----
.scores_df_for_method <- function(bench, method, factor_idx = 1,
                                   group_col = NULL) {
  x <- bench[[method]][["scores_with_meta"]]
  if (is.null(x))
    stop("'scores_with_meta' missing for '", method,
         "'. Run attach_scores_metadata() first.", call. = FALSE)
  score_col <- .resolve_factor_col(x, method, factor_idx)
  grp_raw   <- if (!is.null(group_col) && group_col %in% names(x))
    x[[group_col]] else NA
  grp_chr   <- ifelse(is.na(grp_raw), "NA", as.character(grp_raw))
  tibble::tibble(
    sample      = if ("sample" %in% names(x)) as.character(x$sample)
                  else as.character(seq_len(nrow(x))),
    .sample_idx = seq_len(nrow(x)),
    .score      = as.numeric(x[[score_col]]),
    .group      = grp_chr,
    method      = method
  )
}

# ---- internal: build matched + sign-aligned tidy df for one method ----
# ref_method / ref_factor_idx  — the reference for factor matching
# If method == ref_method the df is returned unchanged.
.scores_df_matched <- function(bench, method, ref_method, ref_factor_idx,
                                group_col = NULL) {
  if (method == ref_method) {
    df <- .scores_df_for_method(bench, method, ref_factor_idx, group_col)
    df$.matched_factor <- ref_factor_idx
    return(df)
  }
  m   <- .find_matched_factor(bench, ref_method, ref_factor_idx, method)
  df  <- .scores_df_for_method(bench, method, m$idx, group_col)
  df$.score          <- df$.score * m$sign   # sign-align to reference
  df$.matched_factor <- m$idx
  df
}

# ---- internal: sign alignment (legacy, used when align_to matches by sign only) ----
.align_factor_signs <- function(dfs, ref_method) {
  stopifnot(ref_method %in% names(dfs))
  ref <- dfs[[ref_method]][, c("sample", ".score")]
  names(ref)[2] <- "ref_score"
  out <- dfs
  for (m in names(dfs)) {
    if (m == ref_method) next
    joined <- merge(out[[m]][, c("sample", ".score")], ref, by = "sample")
    ok <- stats::complete.cases(joined$.score, joined$ref_score)
    if (sum(ok) < 3) next
    r <- stats::cor(joined$.score[ok], joined$ref_score[ok])
    if (is.finite(r) && r < 0) out[[m]]$.score <- -out[[m]]$.score
  }
  out
}

# ---- internal: single scatter panel ----
.plot_scores_panel <- function(df, y_lab, panel_tag = NULL,
                                tag_position = c("tl_in", "none"),
                                palette = NULL, legend_title = "group") {
  tag_position <- match.arg(tag_position)
  df$.group  <- as.character(df$.group)
  df$.group[is.na(df$.group)] <- "NA"
  levs       <- unique(df$.group)
  if (!"NA" %in% levs) levs <- c(levs, "NA")
  df$.group  <- factor(df$.group, levels = levs)
  pal        <- .build_palette(levs, pal = palette)

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = .data[[".sample_idx"]],
                                    y = .data[[".score"]],
                                    fill = .data[[".group"]],
                                    color = .data[[".group"]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_point(shape = 21, size = 2.7, stroke = 0.8, alpha = 0.95) +
    ggplot2::scale_fill_manual(values = pal,  drop = FALSE, name = legend_title) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE, guide = "none") +
    ggplot2::labs(x = "Samples", y = y_lab) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title.x       = ggplot2::element_text(size = 14, colour = "black"),
      axis.title.y       = ggplot2::element_text(size = 14, colour = "black"),
      plot.margin        = ggplot2::margin(6, 6, 6, 6)
    )
  if (tag_position == "tl_in" && !is.null(panel_tag) && nzchar(panel_tag))
    p <- p + ggplot2::annotate("label", x = -Inf, y = Inf, label = panel_tag,
                               hjust = -0.1, vjust = 1.1, size = 4.5,
                               fill = "grey90",
                               label.r = grid::unit(0.15, "lines"))
  p
}

# ---- internal: title grob ----
.title_grob <- function(title_text) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = title_text,
                      fontface = 2, size = 5.5, colour = "#6a1b9a") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 6, 0))
}


#' Plot factor scores for multiple methods (scatter grid)
#'
#' Produces a 2 × 2 grid of scatter plots showing factor scores per sample,
#' coloured by a metadata grouping variable.  When `align_to` is set, each
#' non-reference method automatically displays its **matched** factor — i.e.
#' the factor most correlated (by |r|) with factor `factor_idx` of the
#' reference method — and its sign is flipped if the raw correlation is
#' negative.  Use `align_to = NULL` to disable matching and plot each method's
#' factor `factor_idx` as-is.
#'
#' @param bench Named list returned by [attach_scores_metadata()].
#' @param factor_idx Integer.  Factor index to display for the **reference**
#'   method (`align_to`). Other methods will show their best-matching factor.
#' @param group_col Character.  Column from `scores_with_meta` for colouring.
#' @param methods Character vector of exactly 4 method names to plot.
#' @param panel_tags Optional character vector of length 4 for panel labels.
#' @param tag_position `"tl_in"` (default) or `"none"`.
#' @param title_text Character.  Figure title.
#' @param palette Optional named colour vector for groups.
#' @param legend_position Legend position string (default `"right"`).
#' @param align_to Character.  Reference method for factor matching and sign
#'   alignment.  Set `NULL` to disable (default `"MOFA"`).
#'
#' @return A patchwork ggplot object.
#'
#' @export
plot_scores_scatter_grid <- function(
    bench,
    factor_idx      = 1,
    group_col       = "Treatment",
    methods         = c("MOFA", "FABIA", "MFA", "GFA"),
    panel_tags      = NULL,
    tag_position    = c("tl_in", "none"),
    title_text      = "Factor scores",
    palette         = NULL,
    legend_position = "right",
    align_to        = "MOFA"
) {
  if (length(methods) != 4)
    stop("'methods' must have exactly 4 elements.", call. = FALSE)
  tag_position <- match.arg(tag_position)

  # -----------------------------------------------------------------
  # Build per-method data frames using MATCHED factors when align_to
  # is set; fall back to factor_idx for all methods when it is NULL.
  # -----------------------------------------------------------------
  ref <- if (!is.null(align_to)) {
    if (!align_to %in% methods)
      stop("'align_to' must be one of: ", paste(methods, collapse = ", "),
           call. = FALSE)
    align_to
  } else {
    methods[1]   # reference only used for labelling; no matching done
  }

  dfs <- lapply(methods, function(m) {
    if (!is.null(align_to)) {
      .scores_df_matched(bench, m, ref_method = ref,
                          ref_factor_idx = factor_idx, group_col = group_col)
    } else {
      .scores_df_for_method(bench, m, factor_idx, group_col)
    }
  })
  names(dfs) <- methods

  # y-axis label includes which factor is actually shown
  .y_label <- function(m, df) {
    fi <- if (".matched_factor" %in% names(df)) df$.matched_factor[1] else factor_idx
    paste0("Score (", m, " \u2014 F", fi, ")")
  }

  tags      <- if (is.null(panel_tags)) rep("", 4) else rep_len(panel_tags, 4)
  leg_theme <- ggplot2::theme(legend.position = legend_position,
                              legend.title    = ggplot2::element_text(face = "bold"))

  panels <- lapply(seq_along(methods), function(i) {
    m <- methods[i]
    p <- .plot_scores_panel(
      dfs[[m]],
      y_lab        = .y_label(m, dfs[[m]]),
      panel_tag    = if (tag_position == "tl_in") tags[i] else NULL,
      tag_position = tag_position,
      palette      = palette,
      legend_title = if (is.null(group_col)) "group" else group_col
    )
    p + leg_theme
  })

  panel_grid <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]])
  .title_grob(title_text) / panel_grid +
    patchwork::plot_layout(heights = c(0.15, 1), guides = "collect")
}


#' Plot factor scores as lollipop charts for multiple methods
#'
#' Produces a 2 × 2 grid of lollipop plots showing each sample's factor score
#' as a segment from a baseline.  When `align_to` is set, each non-reference
#' method displays its **best-matched** factor (by max |r| with the reference)
#' and the sign is aligned automatically.
#'
#' @inheritParams plot_scores_scatter_grid
#' @param baseline Baseline: `"method_mean"` (default), `"group_mean"`, or
#'   `"overall_median"`.
#'
#' @return A patchwork ggplot object.
#'
#' @export
plot_scores_lollipop_grid <- function(
    bench,
    factor_idx      = 1,
    group_col       = "Treatment",
    methods         = c("MOFA", "FABIA", "MFA", "GFA"),
    baseline        = c("method_mean", "group_mean", "overall_median"),
    palette         = NULL,
    legend_position = "right",
    align_to        = "MOFA",
    title_text      = "Factor scores"
) {
  if (length(methods) != 4)
    stop("'methods' must have exactly 4 elements.", call. = FALSE)
  baseline <- match.arg(baseline)

  if (!is.null(align_to) && !align_to %in% methods)
    stop("'align_to' must be one of: ", paste(methods, collapse = ", "),
         call. = FALSE)
  ref <- if (!is.null(align_to)) align_to else methods[1]

  # ------------------------------------------------------------------
  # Build matched, sign-aligned score data frames
  # ------------------------------------------------------------------
  dfs <- lapply(methods, function(m) {
    df <- if (!is.null(align_to)) {
      .scores_df_matched(bench, m, ref_method = ref,
                          ref_factor_idx = factor_idx, group_col = group_col)
    } else {
      .scores_df_for_method(bench, m, factor_idx, group_col)
    }
    df$.group <- ifelse(is.na(df$.group), "NA", as.character(df$.group))

    # Compute baseline AFTER sign alignment (sign is already baked in)
    if (baseline == "group_mean") {
      grp_means <- tapply(df$.score, df$.group, mean, na.rm = TRUE)
      df$base   <- grp_means[df$.group]
    } else if (baseline == "method_mean") {
      df$base <- mean(df$.score, na.rm = TRUE)
    } else {
      df$base <- stats::median(df$.score, na.rm = TRUE)
    }
    df
  })
  names(dfs) <- methods

  # y-axis label shows which factor is actually displayed
  .y_label <- function(m, df) {
    fi <- if (".matched_factor" %in% names(df)) df$.matched_factor[1] else factor_idx
    paste0("Score (", m, " \u2014 F", fi, ")")
  }

  .lolli_panel <- function(df, y_lab) {
    levs <- unique(c(as.character(df$.group), "NA"))
    pal  <- .build_palette(levs, pal = palette)
    df$.group <- factor(df$.group, levels = names(pal))
    ggplot2::ggplot(df,
                    ggplot2::aes(x     = .data[[".sample_idx"]],
                                 y     = .data[[".score"]],
                                 fill  = .data[[".group"]],
                                 color = .data[[".group"]])) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = .data[[".sample_idx"]],
                     y    = .data[["base"]],
                     yend = .data[[".score"]]),
        linewidth = 0.7, alpha = 0.9) +
      ggplot2::geom_point(shape = 21, size = 2.7, stroke = 0.8, alpha = 0.95) +
      ggplot2::geom_point(ggplot2::aes(y = .data[["base"]]),
                          shape = 95, size = 6, alpha = 0.7) +
      ggplot2::scale_fill_manual(values = pal, drop = FALSE,
                                  name = if (is.null(group_col)) "group" else group_col) +
      ggplot2::scale_color_manual(values = pal, drop = FALSE, guide = "none") +
      ggplot2::labs(x = "Samples", y = y_lab) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.title.x       = ggplot2::element_text(size = 14, colour = "black"),
        axis.title.y       = ggplot2::element_text(size = 14, colour = "black"),
        plot.margin        = ggplot2::margin(6, 6, 6, 6)
      )
  }

  leg_theme <- ggplot2::theme(legend.position = legend_position,
                              legend.title    = ggplot2::element_text(face = "bold"))
  panels    <- lapply(seq_along(methods), function(i) {
    m <- methods[i]
    .lolli_panel(dfs[[m]], .y_label(m, dfs[[m]])) + leg_theme
  })

  panel_grid <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]])
  .title_grob(title_text) / panel_grid +
    patchwork::plot_layout(heights = c(0.15, 1), guides = "collect")
}


#' Extract scores as a labelled data frame for one method
#'
#' Pulls `scores_with_meta` from a method result entry, identifies the factor
#' column, and returns a tidy tibble suitable for plotting.
#'
#' @param bench Named list returned by [attach_scores_metadata()].
#' @param method Character.  Method name.
#' @param factor_idx Integer.  Factor index (default `1`).
#' @param group_col Character or `NULL`.  Metadata grouping column.
#'
#' @return A [tibble::tibble] with columns `sample`, `.sample_idx`, `.score`,
#'   `.group`, `method`.
#'
#' @export
scores_df_for_method <- function(bench, method, factor_idx = 1,
                                  group_col = NULL) {
  .scores_df_for_method(bench, method, factor_idx, group_col)
}
