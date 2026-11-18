# tests/testthat/test-rabfa.R — graph-guided Bayesian factor analysis

.true_graph <- function(sim) {
  W <- sim$truth$W
  do.call(rbind, lapply(names(W), function(v) do.call(rbind, lapply(
    seq_len(ncol(W[[v]])), function(f) {
      on <- rownames(W[[v]])[W[[v]][, f] != 0]
      data.frame(from = paste0(v, "::", on[1]), to = paste0(v, "::", on[-1]))
    }))))
}

test_that("RABFA is registered", {
  expect_true("RABFA" %in% list_factor_methods()$method)
  expect_true(list_factor_methods()$installed[
    list_factor_methods()$method == "RABFA"])
})

test_that("RABFA without a graph returns unified output and finds factors", {
  res <- run_factor_pipeline(.sim$data, methods = "RABFA", k = 3,
                             rabfa_iter = 300, rabfa_burnin = 150,
                             verbose = FALSE)
  expect_equal(dim(res$RABFA$scores), c(30, 3))
  expect_equal(nrow(res$RABFA$loadings), 60)
  expect_equal(rownames(res$RABFA$scores), colnames(.sim$data[[1]]))
  expect_true(all(grepl("::", rownames(res$RABFA$loadings))))
  expect_gt(max(abs(cor(res$RABFA$scores, .sim$truth$Z[, 1]))), 0.8)
  m <- res$RABFA$meta
  expect_equal(dim(m$selected), c(60, 3))
  expect_true(all(m$lambda > 0))
  expect_length(m$loglik_trace, 300)
  expect_true(is.finite(m$bic))
  expect_equal(m$n_edges, 0)
})

test_that("RABFA accepts edge lists, matrices and per-view lists", {
  g <- .true_graph(.sim)
  r1 <- rabfa(.sim$data, k = 2, graph = g, iter = 120, burnin = 60,
              verbose = FALSE)
  expect_equal(r1$RABFA$meta$n_edges, nrow(unique(g)))

  feats <- rownames(factomics:::.stack_views(.sim$data))
  A <- matrix(0, length(feats), length(feats), dimnames = list(feats, feats))
  A[cbind(g$from, g$to)] <- 1
  r2 <- rabfa(.sim$data, k = 2, graph = A + t(A), iter = 120, burnin = 60,
              verbose = FALSE)
  expect_equal(r2$RABFA$meta$n_edges, r1$RABFA$meta$n_edges)

  # per-view list with plain feature names
  v <- names(.sim$data)[1]
  fv <- rownames(.sim$data[[v]])
  Av <- matrix(0, length(fv), length(fv), dimnames = list(fv, fv))
  Av[1, 2:4] <- Av[2:4, 1] <- 1
  r3 <- rabfa(.sim$data, k = 2, graph = setNames(list(Av), v), iter = 120,
              burnin = 60, verbose = FALSE)
  expect_equal(r3$RABFA$meta$n_edges, 3)
})

test_that("cross-view edges are dropped with a message", {
  f1 <- paste0(names(.sim$data)[1], "::", rownames(.sim$data[[1]])[1])
  f2 <- paste0(names(.sim$data)[2], "::", rownames(.sim$data[[2]])[1])
  expect_message(
    r <- rabfa(.sim$data, k = 2, graph = data.frame(from = f1, to = f2),
               iter = 60, burnin = 30, verbose = FALSE),
    "between different views")
  expect_equal(r$RABFA$meta$n_edges, 0)
})

test_that("RABFA works alongside other methods and downstream tools", {
  res <- run_factor_pipeline(.sim$data, methods = c("PCA", "RABFA"), k = 3,
                             rabfa_iter = 200, rabfa_burnin = 100,
                             verbose = FALSE)
  expect_named(res, c("PCA", "RABFA"))
  expect_s3_class(cluster_results(res, k = 2:4), "factomics_clusters")
  ve <- variance_explained(res, "RABFA")
  expect_true(all(ve$R2 > -1e-8 & ve$R2 < 1 + 1e-8))
})

test_that("build_feature_graph returns within-view edges", {
  g <- build_feature_graph(.sim$data, cutoff = 0.5, max_neighbors = 3)
  expect_named(g, c("from", "to"))
  expect_true(all(sub("::.*", "", g$from) == sub("::.*", "", g$to)))
  expect_true(nrow(g) > 0)
})

test_that("bad graph input and settings are rejected", {
  expect_error(rabfa(.sim$data, k = 2, graph = "kegg", iter = 20,
                     burnin = 10, verbose = FALSE), "rabfa_graph")
  expect_error(rabfa(.sim$data, k = 2, iter = 20, burnin = 20,
                     verbose = FALSE), "burnin")
})
