# tests/testthat/test-utils.R
# Tests for internal helpers — accessed via factomics:::.fn_name

test_that(".stack_views prefixes rownames with view name", {
  lst <- list(
    viewA = matrix(1:6, nrow = 2,
                   dimnames = list(c("f1", "f2"), c("S1", "S2", "S3"))),
    viewB = matrix(7:9, nrow = 1,
                   dimnames = list("g1", c("S1", "S2", "S3")))
  )
  stacked <- factomics:::.stack_views(lst)
  expect_equal(rownames(stacked), c("viewA::f1", "viewA::f2", "viewB::g1"))
  expect_equal(nrow(stacked), 3L)
})

test_that(".as_list_views wraps a single matrix", {
  m   <- matrix(1:4, 2, 2, dimnames = list(NULL, c("S1","S2")))
  out <- factomics:::.as_list_views(m)
  expect_type(out, "list")
  expect_named(out, "omic1")
  expect_identical(out$omic1, m)
})

test_that(".as_list_views passes a list through unchanged", {
  lst <- list(a = matrix(1), b = matrix(2))
  expect_identical(factomics:::.as_list_views(lst), lst)
})

test_that(".align_sign_vec flips sign when correlation is negative", {
  ref <- c(1, 2, 3, 4, 5)
  v   <- -ref
  out <- factomics:::.align_sign_vec(ref, v)
  expect_equal(out, ref)
})

test_that(".align_sign_vec keeps sign when correlation is positive", {
  ref <- c(1, 2, 3, 4, 5)
  v   <- ref * 0.9
  out <- factomics:::.align_sign_vec(ref, v)
  expect_equal(out, v)
})

test_that(".match_factors returns NULL for empty matrix", {
  expect_null(factomics:::.match_factors(matrix(nrow = 0, ncol = 0)))
  expect_null(factomics:::.match_factors(NULL))
})

test_that(".match_factors produces a data frame with columns A, B, corr", {
  skip_if_not_installed("clue")
  set.seed(1)
  A   <- matrix(rnorm(30), 10, 3)
  B   <- A + matrix(rnorm(30, 0, 0.1), 10, 3)   # nearly identical
  C   <- cor(A, B)
  res <- factomics:::.match_factors(C)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("A", "B", "corr") %in% names(res)))
  # Perfect factors should match to themselves
  expect_equal(res$A, res$B)
})

test_that(".build_palette always includes NA colour", {
  pal <- factomics:::.build_palette(c("CRC", "CRC+RT"))
  expect_true("NA" %in% names(pal))
})

test_that(".build_palette respects custom overrides", {
  pal <- factomics:::.build_palette(
    c("CRC", "CRC+RT"),
    pal = c("CRC" = "#FF0000")
  )
  expect_equal(pal[["CRC"]], "#FF0000")
})

test_that(".build_palette errors on unnamed custom palette", {
  expect_error(
    factomics:::.build_palette(c("A", "B"), pal = c("#FF0000", "#00FF00")),
    regexp = "named"
  )
})

test_that(".to_numeric_matrix extracts numeric columns", {
  df <- data.frame(
    sample = paste0("S", 1:4),
    `F1 (MFA)` = c(1.1, 2.2, 3.3, 4.4),
    `F2 (MFA)` = c(5.5, 6.6, 7.7, 8.8),
    check.names = FALSE
  )
  m <- factomics:::.to_numeric_matrix(df, "MFA", id_pref = "sample")
  expect_true(is.matrix(m))
  expect_equal(nrow(m), 4L)
  expect_equal(ncol(m), 2L)
  expect_equal(rownames(m), paste0("S", 1:4))
})

test_that(".to_numeric_matrix errors when no numeric columns found", {
  df <- data.frame(sample = paste0("S", 1:3),
                   label  = c("a", "b", "c"),
                   stringsAsFactors = FALSE)
  expect_error(
    factomics:::.to_numeric_matrix(df, "X"),
    regexp = "No numeric columns"
  )
})