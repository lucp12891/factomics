# tests/testthat/test-pipeline.R

test_that("run_factor_pipeline rejects bad input", {
  # non-list, non-matrix inputs
  expect_error(run_factor_pipeline("not_a_matrix"),
               regexp = "matrix or a named list")
  expect_error(run_factor_pipeline(42L),
               regexp = "matrix or a named list")

  # mismatched sample names across views
  expect_error(
    run_factor_pipeline(list(
      a = matrix(1:4, 2, 2, dimnames = list(NULL, c("S1", "S2"))),
      b = matrix(1:6, 3, 2, dimnames = list(NULL, c("S1", "S3")))
    )),
    regexp = "same sample names"
  )

  # matrix without column names
  expect_error(
    run_factor_pipeline(matrix(1:4, 2, 2)),
    regexp = "column.*sample.*names|sample.*names"
  )
})

test_that("run_factor_pipeline returns correct structure for MFA", {
  skip_if_not_installed("FactoMineR")
  res <- run_factor_pipeline(.data_list, methods = "MFA", k = 2,
                              scale_features = TRUE, verbose = FALSE)
  expect_type(res, "list")
  expect_named(res, "MFA")
  expect_true(is.matrix(res$MFA$scores))
  expect_true(is.matrix(res$MFA$loadings))
  expect_equal(nrow(res$MFA$scores), .n_samples)
  expect_equal(ncol(res$MFA$scores), 2L)
  expect_equal(rownames(res$MFA$scores), .snames)
  # loadings: one row per feature across all views
  total_feats <- .n_feat1 + .n_feat2 + .n_feat3
  expect_equal(nrow(res$MFA$loadings), total_feats)
})

test_that("MFA loadings have view:: prefix in row names", {
  skip_if_not_installed("FactoMineR")
  res <- run_factor_pipeline(.data_list, methods = "MFA", k = 2,
                              verbose = FALSE)
  rn <- rownames(res$MFA$loadings)
  expect_true(all(grepl("::", rn)))
  prefixes <- unique(sub("::.*$", "", rn))
  expect_setequal(prefixes, names(.data_list))
})

test_that("run_factor_pipeline works with a single matrix (not a list)", {
  # MFA requires multiple groups, so use FABIA for single-matrix test
  skip_if_not_installed("fabia")
  res <- run_factor_pipeline(.single_mat, methods = "FABIA", k = 2,
                              verbose = FALSE)
  expect_named(res, "FABIA")
  expect_equal(nrow(res$FABIA$scores), .n_samples)
  expect_equal(rownames(res$FABIA$scores), .snames)
})

test_that("run_factor_pipeline with scale_features=FALSE still works", {
  skip_if_not_installed("FactoMineR")
  res <- run_factor_pipeline(.data_list, methods = "MFA", k = 2,
                              scale_features = FALSE, verbose = FALSE)
  expect_equal(nrow(res$MFA$scores), .n_samples)
})

test_that("run_factor_pipeline ignores unknown methods with a warning", {
  skip_if_not_installed("FactoMineR")
  expect_warning(
    run_factor_pipeline(.data_list, methods = c("MFA", "UNICORN"),
                         k = 2, verbose = FALSE),
    regexp = "Unknown method"
  )
})

test_that("run_factor_pipeline method names are case-insensitive", {
  skip_if_not_installed("FactoMineR")
  res_upper <- run_factor_pipeline(.data_list, methods = "MFA",  k = 2,
                                    verbose = FALSE)
  res_lower <- run_factor_pipeline(.data_list, methods = "mfa",  k = 2,
                                    verbose = FALSE)
  res_mixed <- run_factor_pipeline(.data_list, methods = "Mfa",  k = 2,
                                    verbose = FALSE)
  expect_named(res_upper, "MFA")
  expect_named(res_lower, "MFA")
  expect_named(res_mixed, "MFA")
})

test_that("run_factor_pipeline factor columns are named F1, F2, ...", {
  skip_if_not_installed("FactoMineR")
  res <- run_factor_pipeline(.data_list, methods = "MFA", k = 3,
                              verbose = FALSE)
  expect_equal(colnames(res$MFA$scores),   paste0("F", 1:3))
  expect_equal(colnames(res$MFA$loadings), paste0("F", 1:3))
})

test_that("run_factor_pipeline with FABIA returns expected structure", {
  skip_if_not_installed("fabia")
  res <- run_factor_pipeline(.data_list, methods = "FABIA", k = 2,
                              verbose = FALSE)
  expect_named(res, "FABIA")
  expect_true(is.matrix(res$FABIA$scores))
  expect_true(is.matrix(res$FABIA$loadings))
  expect_equal(nrow(res$FABIA$scores), .n_samples)
  expect_equal(rownames(res$FABIA$scores), .snames)
})

test_that("build_loading_comparison pads missing features with zeros", {
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("fabia")
  res <- run_factor_pipeline(.data_list,
                              methods = c("FABIA", "MFA"), k = 2,
                              verbose = FALSE)
  cmp <- build_loading_comparison(res)
  expect_named(cmp, c("FABIA", "MFA"))
  # both matrices must have the same row names (union)
  expect_equal(rownames(cmp$FABIA), rownames(cmp$MFA))
})
