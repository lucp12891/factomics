# tests/testthat/test-methods.R — additional factor methods

.check_fit <- function(res, m, n, p, k) {
  expect_true(m %in% names(res))
  expect_equal(dim(res[[m]]$scores), c(n, k))
  expect_equal(nrow(res[[m]]$loadings), p)
  expect_equal(colnames(res[[m]]$scores), paste0("F", seq_len(k)))
  expect_true(all(grepl("::", rownames(res[[m]]$loadings))))
}

test_that("list_factor_methods describes every method", {
  lm <- list_factor_methods()
  expect_setequal(lm$method,
                  c("PCA", "MFA", "FABIA", "FABIA_ENS", "GFA", "MOFA",
                    "SGCCA", "JIVE", "RABFA"))
  expect_type(lm$installed, "logical")
})

test_that("PCA runs without optional packages and recovers group signal", {
  res <- run_factor_pipeline(.sim$data, methods = "PCA", k = 3,
                             verbose = FALSE)
  .check_fit(res, "PCA", 30, 60, 3)
  expect_s3_class(res, "factomics_results")
  expect_gt(max(abs(cor(res$PCA$scores, .sim$truth$Z[, 1]))), 0.8)
  expect_length(res$PCA$meta$variance_percent, 3)
})

test_that("SGCCA returns unified output", {
  skip_if_not_installed("mixOmics")
  res <- run_factor_pipeline(.sim$data, methods = "SGCCA", k = 2,
                             verbose = FALSE)
  .check_fit(res, "SGCCA", 30, 60, 2)
})

test_that("JIVE returns unified output", {
  skip_if_not_installed("r.jive")
  res <- run_factor_pipeline(.sim$data, methods = "JIVE", k = 2,
                             verbose = FALSE)
  .check_fit(res, "JIVE", 30, 60, 2)
})

test_that("GFA returns unified output", {
  skip_if_not_installed("GFA")
  res <- suppressWarnings(run_factor_pipeline(
    .sim$data, methods = "GFA", k = 2, verbose = FALSE,
    gfa_iter_burnin = 50, gfa_iter_max = 100))
  .check_fit(res, "GFA", 30, 60, 2)
})

test_that("a failing method is skipped with a warning", {
  skip_if_not_installed("FactoMineR")
  one_view <- list(a = .sim$data[[1]])
  expect_warning(
    res <- run_factor_pipeline(one_view, methods = c("PCA", "MFA"), k = 2,
                               verbose = FALSE),
    regexp = "MFA failed")
  expect_named(res, "PCA")
  expect_error(run_factor_pipeline(one_view, methods = "MFA", k = 2,
                                   verbose = FALSE, stop_on_error = TRUE),
               regexp = "two views")
})

test_that("subsetting results keeps class and stored views", {
  res <- run_factor_pipeline(.sim$data, methods = "PCA", k = 2,
                             verbose = FALSE)
  sub <- res["PCA"]
  expect_s3_class(sub, "factomics_results")
  expect_false(is.null(attr(sub, "views")))
})

test_that("simulate_multiomics returns consistent structure", {
  sim <- simulate_multiomics(n_samples = 20, views = c(a = 10, b = 5),
                             seed = 3)
  expect_named(sim$data, c("a", "b"))
  expect_equal(dim(sim$data$a), c(10, 20))
  expect_equal(sim$metadata$sample, colnames(sim$data$a))
  expect_length(sim$truth$feature_labels, 15)
  expect_true(all(sim$truth$feature_labels %in% 0:1))
})
