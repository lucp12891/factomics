# tests/testthat/test-metadata.R

# ---- shared minimal result fixture ----
.make_mfa_res <- function() {
  skip_if_not_installed("FactoMineR")
  run_factor_pipeline(.data_list, methods = "MFA", k = 2, verbose = FALSE)
}

# ---- attach_scores_metadata ----

test_that("attach_scores_metadata adds scores_with_meta slot", {
  res    <- .make_mfa_res()
  merged <- attach_scores_metadata(res, .metadata, sample_col = "sample")
  expect_true(!is.null(merged$MFA$scores_with_meta))
  expect_s3_class(merged$MFA$scores_with_meta, "data.frame")
})

test_that("attach_scores_metadata preserves all samples", {
  res    <- .make_mfa_res()
  merged <- attach_scores_metadata(res, .metadata, sample_col = "sample")
  swd    <- merged$MFA$scores_with_meta
  expect_equal(nrow(swd), .n_samples)
})

test_that("attach_scores_metadata joins metadata columns correctly", {
  res    <- .make_mfa_res()
  merged <- attach_scores_metadata(res, .metadata, sample_col = "sample")
  swd    <- merged$MFA$scores_with_meta
  expect_true("Treatment" %in% names(swd))
  expect_true("Batch"     %in% names(swd))
})

test_that("attach_scores_metadata works with rownames as sample IDs", {
  meta2             <- .metadata[, -1, drop = FALSE]   # drop sample col
  rownames(meta2)   <- .metadata$sample
  res               <- .make_mfa_res()
  expect_no_error(
    attach_scores_metadata(res, meta2, sample_col = "sample")
  )
})

test_that("attach_scores_metadata replace_scores=TRUE stores matrix backup", {
  res    <- .make_mfa_res()
  merged <- attach_scores_metadata(res, .metadata, sample_col = "sample",
                                    replace_scores = TRUE)
  expect_true(!is.null(merged$MFA$scores_matrix))
  expect_true(is.matrix(merged$MFA$scores_matrix))
  expect_s3_class(merged$MFA$scores, "data.frame")
})

test_that("attach_scores_metadata errors on missing sample column", {
  res      <- .make_mfa_res()
  # data.frame with no sample col AND only default integer rownames → must error
  meta_bad <- data.frame(x = 1:5)
  expect_error(
    attach_scores_metadata(res, meta_bad, sample_col = "sample"),
    regexp = "sample|non-default rownames"
  )
})

# ---- merge_scores_with_metadata ----

test_that("merge_scores_with_metadata returns by_method and all_long", {
  res  <- .make_mfa_res()
  out  <- merge_scores_with_metadata(res, .metadata, sample_col = "sample")
  expect_named(out, c("by_method", "all_long"))
  expect_true("MFA" %in% names(out$by_method))
  expect_true("method" %in% names(out$all_long))
})

test_that("merge_scores_with_metadata inner join reduces rows on mismatch", {
  res        <- .make_mfa_res()
  meta_small <- .metadata[1:6, ]   # only first 6 samples
  out_left   <- merge_scores_with_metadata(res, meta_small,
                                            sample_col = "sample",
                                            join = "left")
  out_inner  <- merge_scores_with_metadata(res, meta_small,
                                            sample_col = "sample",
                                            join = "inner")
  expect_equal(nrow(out_left$by_method$MFA),  .n_samples)
  expect_equal(nrow(out_inner$by_method$MFA), 6L)
})

# ---- split_loadings_by_omic ----

test_that("split_loadings_by_omic returns named list", {
  res     <- .make_mfa_res()
  ldf     <- data.frame(res$MFA$loadings, check.names = FALSE)
  ldf$feature <- rownames(ldf)
  splits  <- split_loadings_by_omic(
    ldf,
    patterns = c(microbiome = "^microbiome::",
                 metabolome = "^metabolome::",
                 lipidome   = "^lipidome::")
  )
  expect_named(splits, c("microbiome", "metabolome", "lipidome"))
  expect_equal(nrow(splits$microbiome), .n_feat1)
  expect_equal(nrow(splits$metabolome), .n_feat2)
  expect_equal(nrow(splits$lipidome),   .n_feat3)
})

test_that("split_loadings_by_omic errors on missing column", {
  df <- data.frame(x = 1:3)
  expect_error(split_loadings_by_omic(df, column = "feature",
                                       patterns = c(a = "^a")),
               regexp = "feature")
})

test_that("split_loadings_by_omic errors on unnamed patterns", {
  res <- .make_mfa_res()
  ldf <- data.frame(res$MFA$loadings, check.names = FALSE)
  ldf$feature <- rownames(ldf)
  expect_error(
    split_loadings_by_omic(ldf, patterns = c("^microbiome::")),
    regexp = "named"
  )
})
