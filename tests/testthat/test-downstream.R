# tests/testthat/test-downstream.R — embeddings, clustering, plots, workflow

.res <- function() run_factor_pipeline(.sim$data, methods = "PCA", k = 3,
                                       verbose = FALSE)
.truth <- stats::setNames(.sim$metadata$group, .sim$metadata$sample)

test_that("compute_embedding returns sample + Dim columns", {
  x <- .res()$PCA$scores
  e <- compute_embedding(x, "pca")
  expect_named(e, c("sample", "Dim1", "Dim2"))
  expect_equal(e$sample, rownames(x))
  skip_if_not_installed("uwot")
  expect_equal(nrow(compute_embedding(x, "umap")), 30)
  skip_if_not_installed("Rtsne")
  expect_equal(nrow(compute_embedding(x, "tsne")), 30)
})

test_that("embed_results covers methods and the concatenated data", {
  emb <- embed_results(.res(), type = "pca")
  expect_s3_class(emb, "factomics_embedding")
  expect_setequal(levels(emb$method), c("PCA", "Concatenated data"))
  p <- plot_embedding(emb, metadata = .sim$metadata, color_by = "group")
  expect_s3_class(p, "ggplot")
  expect_error(plot_embedding(emb, metadata = .sim$metadata,
                              color_by = "nope"), "not found")
  expect_error(embed_results(.res(), methods = "MOFA"), "not found")
})

test_that("cluster_samples picks k by silhouette and recovers groups", {
  cl <- cluster_samples(.res()$PCA$scores, k = 2:5)
  expect_s3_class(cl, "factomics_clustering")
  expect_true(cl$k %in% 3:4)
  expect_equal(cl$k, cl$by_k$k[which.max(cl$by_k$silhouette)])
  expect_gt(adjusted_rand_index(cl$cluster, .truth[names(cl$cluster)]), 0.8)
  for (m in c("hclust", "pam"))
    expect_length(cluster_samples(.res()$PCA$scores, method = m, k = 3)$cluster,
                  30)
})

test_that("cluster_results / compare_clusterings / plot_clusters work", {
  cl <- cluster_results(.res(), cluster_method = "hclust", k = 2:4)
  expect_s3_class(cl, "factomics_clusters")
  expect_equal(nrow(cl$summary), 2)
  cmp <- compare_clusterings(cl, truth = .truth)
  expect_equal(dim(cmp$agreement), c(2, 2))
  expect_equal(unname(diag(cmp$agreement)), c(1, 1))
  expect_named(cmp$vs_truth, c("method", "ARI"))
  expect_s3_class(plot_clusters(cl), "ggplot")
  expect_s3_class(plot_clusters(cl, "agreement", truth = .truth), "ggplot")
})

test_that("adjusted_rand_index behaves", {
  expect_equal(adjusted_rand_index(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
  expect_lt(abs(adjusted_rand_index(rep(1:2, 50), rep(1:2, each = 50))), 0.05)
})

test_that("variance_explained is bounded and has a Total row", {
  ve <- variance_explained(.res())
  expect_true(all(ve$R2 >= -1e-8 & ve$R2 <= 1 + 1e-8))
  tot <- ve[ve$factor == "Total", ]
  expect_equal(nrow(tot), 2)
  expect_s3_class(plot_variance_explained(.res()), "ggplot")
  expect_s3_class(plot_variance_explained(.res(), type = "total"), "ggplot")
})

test_that("factor plots and associations work", {
  r <- .res()
  expect_s3_class(plot_factor_scores(r, metadata = .sim$metadata,
                                     color_by = "group"), "ggplot")
  expect_s3_class(plot_top_loadings(r, "PCA", 1, n = 10), "ggplot")
  expect_s3_class(plot_factor_heatmap(r, "PCA", .sim$metadata,
                                      order_by = "group"), "ggplot")
  fa <- factor_associations(r, .sim$metadata, covariates = c("group", "age"))
  expect_true(all(c("p_value", "p_adj") %in% names(fa)))
  expect_lt(min(fa$p_adj[fa$covariate == "group"]), 0.01)
  expect_s3_class(plot_factor_associations(r, .sim$metadata), "ggplot")
})

test_that("factomics() runs the full workflow and exports results", {
  skip_if_not_installed("FactoMineR")
  fx <- factomics(.sim$data, metadata = .sim$metadata,
                  methods = c("PCA", "MFA"), k = 3, embeddings = "pca",
                  n_clusters = 2:4, verbose = FALSE)
  expect_s3_class(fx, "factomics")
  expect_false(is.null(fx$benchmark))
  expect_false(is.null(fx$results$PCA$scores_with_meta))
  expect_output(print(fx), "factomics workflow")
  files <- export_results(fx, tempfile())
  expect_true(all(file.exists(files)))
  expect_true(any(grepl("clusters.csv", files)))
})
