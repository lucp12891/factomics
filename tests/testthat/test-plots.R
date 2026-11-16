# tests/testthat/test-plots.R

.make_bench_with_meta <- function() {
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("fabia")
  res <- run_factor_pipeline(.data_list,
                             methods = c("FABIA", "MFA"), k = 2,
                             verbose = FALSE)
  attach_scores_metadata(res, .metadata, sample_col = "sample")
}

# ---- scores_df_for_method ----

test_that("scores_df_for_method returns a tibble with expected columns", {
  bench <- .make_bench_with_meta()
  df    <- scores_df_for_method(bench, "MFA", factor_idx = 1,
                                group_col = "Treatment")
  expect_s3_class(df, "data.frame")
  expect_true(all(c("sample", ".sample_idx", ".score", ".group", "method")
                  %in% names(df)))
  expect_equal(nrow(df), .n_samples)
})

test_that("scores_df_for_method errors when scores_with_meta is missing", {
  skip_if_not_installed("FactoMineR")
  res <- run_factor_pipeline(.data_list, methods = "MFA", k = 2,
                             verbose = FALSE)
  expect_error(
    scores_df_for_method(res, "MFA", factor_idx = 1),
    regexp = "scores_with_meta"
  )
})

# ---- plot_scores_scatter_grid ----

test_that("plot_scores_scatter_grid works with any number of methods", {
  bench <- .make_bench_with_meta()
  p <- plot_scores_scatter_grid(bench, methods = c("FABIA", "MFA"),
                                group_col = "Treatment")
  expect_s3_class(p, "patchwork")
  p1 <- plot_scores_scatter_grid(bench, methods = "MFA",
                                 group_col = "Treatment")
  expect_s3_class(p1, "patchwork")
  expect_error(plot_scores_scatter_grid(bench, align_to = "MOFA"),
               regexp = "align_to")
})

test_that("plot_scores_scatter_grid returns a ggplot/patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("fabia")
  
  # Need 4 methods — use same method repeated with different names via 2 methods
  # Build a bench that has FABIA & MFA, then duplicate slots for testing
  bench  <- .make_bench_with_meta()
  bench$GFA  <- bench$FABIA   # duplicate as stand-in
  bench$MOFA <- bench$MFA
  
  p <- plot_scores_scatter_grid(
    bench, factor_idx = 1, group_col = "Treatment",
    methods = c("FABIA", "MFA", "GFA", "MOFA"),
    align_to = NULL
  )
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

# ---- plot_scores_lollipop_grid ----

test_that("plot_scores_lollipop_grid defaults to all methods with metadata", {
  bench <- .make_bench_with_meta()
  p <- plot_scores_lollipop_grid(bench, group_col = "Treatment")
  expect_s3_class(p, "patchwork")
})

test_that("plot_scores_lollipop_grid baseline options are accepted", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("fabia")
  
  bench       <- .make_bench_with_meta()
  bench$GFA   <- bench$FABIA
  bench$MOFA  <- bench$MFA
  
  for (bl in c("method_mean", "group_mean", "overall_median")) {
    expect_no_error(
      plot_scores_lollipop_grid(
        bench, factor_idx = 1, group_col = "Treatment",
        methods  = c("FABIA", "MFA", "GFA", "MOFA"),
        baseline = bl, align_to = NULL
      )
    )
  }
})