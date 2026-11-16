# ============================================================
#  factomics — metadata merging utilities
# ============================================================

# Internal helper: returns TRUE when rownames are the default "1","2","3",... sequence
.has_default_rownames <- function(df) {
  rn <- rownames(df)
  is.null(rn) || identical(rn, as.character(seq_len(nrow(df))))
}

#' Merge factor scores with sample metadata
#'
#' Performs a join between the scores matrices in a factor-analysis result
#' object and a metadata data frame, returning one wide data frame per method
#' as well as a stacked long table.
#'
#' @param results Named list returned by [run_factor_pipeline()].
#' @param metadata Data frame (or coercible object) with at least one column
#'   identifying samples.  Can also use `rownames` as sample IDs.
#' @param sample_col Character.  Column name in `metadata` holding sample IDs
#'   (default `"sample"`).  If not present, `rownames(metadata)` is used.
#' @param methods Character vector of method names to process (default: all
#'   methods found in `results`).
#' @param join `"left"` (default) or `"inner"`.  Join type between scores and
#'   metadata.
#'
#' @return A list with:
#' \describe{
#'   \item{`by_method`}{Named list of data frames (one per method).}
#'   \item{`all_long`}{Single row-stacked data frame with a `method` column.}
#' }
#'
#' @export
merge_scores_with_metadata <- function(
    results,
    metadata,
    sample_col = "sample",
    methods    = names(results),
    join       = c("left", "inner")
) {
  join <- match.arg(join)
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!(sample_col %in% names(meta))) {
    if (!.has_default_rownames(meta)) {
      meta[[sample_col]] <- rownames(meta)
    } else {
      stop("`metadata` must have column '", sample_col,
           "' or non-default rownames as sample IDs.", call. = FALSE)
    }
  }
  meta[[sample_col]] <- as.character(meta[[sample_col]])
  meta <- meta[!duplicated(meta[[sample_col]]), , drop = FALSE]

  by_method <- list()
  for (m in methods) {
    if (is.null(results[[m]]$scores)) next
    S    <- results[[m]]$scores
    Sdf  <- as.data.frame(S, check.names = FALSE, stringsAsFactors = FALSE)
    if (!("sample" %in% names(Sdf))) {
      rn <- rownames(Sdf)
      if (is.null(rn))
        stop("Scores for '", m, "' have no rownames.", call. = FALSE)
      Sdf$sample <- rn
    }
    Sdf$sample <- as.character(Sdf$sample)
    by_map     <- stats::setNames(sample_col, "sample")
    out <- merge(Sdf, meta, by.x = "sample", by.y = sample_col,
                 all.x = (join == "left"), sort = FALSE)
    out <- out[order(match(out$sample, Sdf$sample)), , drop = FALSE]
    rownames(out) <- NULL
    fac_cols <- grep("^F\\d+(\\s*\\(|$)", names(Sdf), value = TRUE)
    front    <- c("sample", fac_cols)
    out      <- out[, c(intersect(front, names(out)),
                        setdiff(names(out), front)), drop = FALSE]
    out$method    <- m
    by_method[[m]] <- out
  }

  all_long <- do.call(rbind, unname(by_method))
  list(by_method = by_method, all_long = all_long)
}


#' Attach sample metadata to a factor-analysis results object (in-place)
#'
#' For each requested method, merges the score matrix with a metadata data
#' frame and stores the result in `results[[method]]$scores_with_meta`.
#' Optionally replaces the original score matrix.
#'
#' @param results Named list returned by [run_factor_pipeline()].
#' @param metadata Metadata data frame (see [merge_scores_with_metadata()]).
#' @param sample_col Character.  Column holding sample IDs (default `"sample"`).
#' @param methods Character vector of methods to process.
#' @param store_name Character.  Name under which the merged data frame is
#'   stored in each method's sub-list (default `"scores_with_meta"`).
#' @param replace_scores Logical.  If `TRUE`, replace `scores` with the merged
#'   data frame (original matrix is saved as `scores_matrix`).
#'
#' @return The modified `results` list.
#'
#' @export
attach_scores_metadata <- function(
    results,
    metadata,
    sample_col     = "sample",
    methods        = names(results),
    store_name     = "scores_with_meta",
    replace_scores = FALSE
) {
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!(sample_col %in% names(meta))) {
    if (!.has_default_rownames(meta)) {
      meta[[sample_col]] <- rownames(meta)
    } else {
      stop("`metadata` must have column '", sample_col,
           "' or non-default rownames as sample IDs.", call. = FALSE)
    }
  }
  meta[[sample_col]] <- as.character(meta[[sample_col]])
  meta <- meta[!duplicated(meta[[sample_col]]), , drop = FALSE]

  for (m in methods) {
    if (is.null(results[[m]]$scores)) next
    Sdf        <- as.data.frame(results[[m]]$scores,
                                check.names = FALSE, stringsAsFactors = FALSE)
    if (!("sample" %in% names(Sdf))) {
      rn <- rownames(Sdf)
      if (is.null(rn))
        stop("Scores for '", m, "' have no rownames.", call. = FALSE)
      Sdf$sample <- rn
    }
    Sdf$sample <- as.character(Sdf$sample)
    merged     <- merge(Sdf, meta, by.x = "sample",
                        by.y = sample_col, all.x = TRUE, sort = FALSE)
    merged     <- merged[match(Sdf$sample, merged$sample), , drop = FALSE]
    rownames(merged) <- NULL
    merged$method <- m

    if (isTRUE(replace_scores)) {
      results[[m]]$scores_matrix <- results[[m]]$scores
      results[[m]]$scores        <- merged
    } else {
      results[[m]][[store_name]] <- merged
    }
  }
  results
}


#' Split a loading data frame by feature prefix
#'
#' Given a data frame with a `feature` column whose values follow the
#' `"view::feature_id"` convention, return a named list of subsets—one per
#' pattern match.
#'
#' @param df Data frame containing a feature identifier column.
#' @param column Character.  Name of the column to match against (default
#'   `"feature"`).
#' @param patterns Named character vector of regex patterns.  Names become the
#'   list element names.
#'
#' @return Named list of data frames.
#'
#' @examples
#' \dontrun{
#' split_loadings_by_omic(
#'   loading_df,
#'   patterns = c(microbiome = "^microbiome",
#'                metabolome = "^metabolome",
#'                lipidome   = "^lipidome")
#' )
#' }
#'
#' @export
split_loadings_by_omic <- function(df,
                                    column   = "feature",
                                    patterns = NULL) {
  if (!(column %in% names(df)))
    stop("Column '", column, "' not found in data frame.", call. = FALSE)
  if (is.null(names(patterns)))
    stop("'patterns' must be a *named* character vector.", call. = FALSE)
  lapply(patterns, function(pat) df[grepl(pat, df[[column]]), , drop = FALSE])
}
