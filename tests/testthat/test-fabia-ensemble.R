# tests/testthat/test-fabia-ensemble.R — ensemble FABIA and seed stability

test_that("FABIA_ENS is registered and settings are validated", {
  expect_true("FABIA_ENS" %in% list_factor_methods()$method)
  ctl <- fabia_ensemble_control(runs = 5, sample_fraction = 0.8)
  expect_equal(ctl$runs, 5L)
  expect_error(fabia_ensemble_control(runs = 1))
  expect_error(fabia_ensemble_control(sample_fraction = 1.5))
})

test_that("bicluster Jaccard index counts row-column combinations", {
  R <- cbind(c(1, 1, 0, 0), c(1, 1, 0, 0), c(0, 0, 1, 1))
  C <- cbind(c(1, 1, 0), c(1, 0, 0), c(0, 0, 1))
  J <- factomics:::.bicluster_jaccard(R, C)
  expect_equal(J[1, 2], 2 / 4)     # 4 cells vs 2 cells, all shared
  expect_equal(J[1, 3], 0)
  expect_equal(diag(J), rep(1, 3))
})

test_that("FABIA_ENS returns robust factors in the unified format", {
  skip_if_not_installed("fabia")
  res <- run_factor_pipeline(.sim$data, methods = "FABIA_ENS", k = 3,
                             fabia_cyc = 300,
                             fabia_ens_control = list(runs = 8),
                             verbose = FALSE)
  f <- res$FABIA_ENS
  expect_s3_class(res, "factomics_results")
  expect_equal(rownames(f$scores), colnames(.sim$data[[1]]))
  expect_equal(nrow(f$loadings), 60)
  expect_true(ncol(f$scores) >= 1 && ncol(f$scores) <= 3)
  expect_false(anyNA(f$scores))
  m <- f$meta
  expect_equal(m$runs, 8L)
  expect_true(all(m$reliability >= 0.5 & m$reliability <= 1))
  expect_equal(names(m$biclusters), colnames(f$scores))
  expect_true(all(diff(m$reliability) <= 0))
  expect_s3_class(m$tree, "hclust")
  expect_gt(max(abs(cor(f$scores, .sim$truth$Z[, 1]))), 0.8)
  expect_s3_class(plot_fabia_ensemble(res), "ggplot")
})

test_that("subsampling still gives scores for every sample", {
  skip_if_not_installed("fabia")
  fit <- fabia_ensemble(.sim$data, k = 3, cyc = 300, verbose = FALSE,
                        control = fabia_ensemble_control(
                          runs = 8, sample_fraction = 0.7,
                          min_reliability = 0.3))
  expect_equal(nrow(fit$scores), 30)
  expect_false(anyNA(fit$scores))
})

test_that("an impossible reliability threshold gives a clear error", {
  skip_if_not_installed("fabia")
  expect_error(
    fabia_ensemble(.sim$data, k = 3, cyc = 100, verbose = FALSE,
                   control = fabia_ensemble_control(runs = 3, similarity = 1,
                                                    min_reliability = 1,
                                                    membership = 1,
                                                    thresZ = 50)),
    "bicluster")
})

test_that("factor_stability compares methods across seeds", {
  st <- factor_stability(.sim$data, methods = c("PCA", "MFA"), k = 2,
                         seeds = 1:3, truth = .sim$truth$Z)
  expect_s3_class(st, "factomics_stability")
  expect_equal(nrow(st$pairs), 6)          # 3 pairs x 2 methods
  expect_equal(nrow(st$recovery), 6)       # 3 seeds x 2 methods
  # PCA does not depend on the seed
  expect_equal(st$summary$mean_stability[st$summary$method == "PCA"], 1,
               tolerance = 1e-6)
  expect_s3_class(plot_factor_stability(st), "ggplot")
  expect_s3_class(plot_factor_stability(st, "recovery"), "ggplot")
  expect_error(factor_stability(.sim$data, methods = "PCA", seeds = 1))
})
