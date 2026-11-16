# tests/testthat/test-benchmark.R

.make_two_method_res <- function() {
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("fabia")
  run_factor_pipeline(.data_list, methods = c("FABIA", "MFA"),
                      k = 2, verbose = FALSE)
}

test_that("benchmark_factor_methods errors with fewer than 2 methods", {
  skip_if_not_installed("FactoMineR")
  res <- run_factor_pipeline(.data_list, methods = "MFA", k = 2,
                             verbose = FALSE)
  expect_error(
    benchmark_factor_methods(res, methods = "MFA"),
    regexp = "at least two"
  )
})

test_that("benchmark_factor_methods returns expected top-level names", {
  res   <- .make_two_method_res()
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"))
  
  expected_names <- c(
    "available_methods", "scores", "loadings",
    "pairwise_scores", "pairwise_loadings",
    "pair_heatmaps_scores", "pair_heatmaps_loadings",
    "pair_scatter_scores", "pair_scatter_loadings",
    "all_scores_heatmap",
    "summary_table_scores", "summary_table_loadings",
    "auc_results",
    "pairs_scores_by_factor",
    "pairs_loadings_by_dataset_and_factor"
  )
  expect_named(bench, expected_names, ignore.order = TRUE)
})

test_that("benchmark_factor_methods available_methods matches input", {
  res   <- .make_two_method_res()
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"))
  expect_setequal(bench$available_methods, c("FABIA", "MFA"))
})

test_that("summary_table_scores is a data frame with expected columns", {
  res   <- .make_two_method_res()
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"))
  tbl   <- bench$summary_table_scores
  expect_s3_class(tbl, "data.frame")
  expect_true("pair"                   %in% names(tbl))
  expect_true("mean_abs_corr_scores"   %in% names(tbl))
  expect_true("median_abs_corr_scores" %in% names(tbl))
  expect_true("max_abs_corr_scores"    %in% names(tbl))
})

test_that("summary_table_loadings has a dataset column", {
  res   <- .make_two_method_res()
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"))
  tbl   <- bench$summary_table_loadings
  expect_s3_class(tbl, "data.frame")
  expect_true("dataset" %in% names(tbl))
})

test_that("pair_heatmaps_scores contains ggplot objects", {
  skip_if_not_installed("ggplot2")
  res   <- .make_two_method_res()
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"))
  hp    <- bench$pair_heatmaps_scores[["FABIA|MFA"]]
  expect_s3_class(hp, "gg")
})

test_that("all_scores_heatmap is a ggplot object", {
  skip_if_not_installed("ggplot2")
  res   <- .make_two_method_res()
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"))
  expect_s3_class(bench$all_scores_heatmap, "gg")
})

test_that("scores matrices in benchmark are scaled when center_scale=TRUE", {
  res   <- .make_two_method_res()
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"),
                                    center_scale = TRUE)
  for (m in bench$available_methods) {
    col_means <- colMeans(bench$scores[[m]], na.rm = TRUE)
    expect_true(all(abs(col_means) < 1e-10),
                label = paste0(m, " scores should be zero-mean after scaling"))
  }
})

test_that("center_scale=FALSE leaves matrices unscaled", {
  res         <- .make_two_method_res()
  bench_raw   <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"),
                                          center_scale = FALSE)
  bench_scaled <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"),
                                           center_scale = TRUE)
  # They should differ for at least one method
  differ <- any(vapply(bench_raw$available_methods, function(m) {
    !isTRUE(all.equal(bench_raw$scores[[m]], bench_scaled$scores[[m]]))
  }, logical(1)))
  expect_true(differ)
})

test_that("auc_results is NULL when no ground truth supplied", {
  res   <- .make_two_method_res()
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"))
  expect_null(bench$auc_results)
})

test_that("auc_results is computed when ground_truth_sample is supplied", {
  skip_if_not_installed("pROC")
  res <- .make_two_method_res()
  gt  <- stats::setNames(as.integer(seq_len(.n_samples) <= 6), .snames)
  bench <- benchmark_factor_methods(res, methods = c("FABIA", "MFA"),
                                    ground_truth_sample = gt)
  expect_false(is.null(bench$auc_results))
  expect_true("AUC" %in% names(bench$auc_results))
  expect_true(all(bench$auc_results$AUC >= 0 & bench$auc_results$AUC <= 1,
                  na.rm = TRUE))
})