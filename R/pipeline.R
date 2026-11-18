# ============================================================
#  factomics — factor pipeline
# ============================================================

# Registry of supported methods ------------------------------------------------
.method_registry <- function() {
  data.frame(
    method      = c("PCA", "MFA", "FABIA", "GFA", "MOFA", "SGCCA", "JIVE",
                    "RABFA"),
    package     = c("stats", "FactoMineR", "fabia", "GFA", "MOFA2",
                    "mixOmics", "r.jive", "factomics"),
    family      = c("Matrix factorisation", "Matrix factorisation",
                    "Bayesian / biclustering", "Bayesian",
                    "Bayesian", "Correlation-based", "Matrix factorisation",
                    "Bayesian / graph-guided"),
    description = c(
      "Block-weighted PCA on concatenated views",
      "Multiple Factor Analysis (Escofier & Pages)",
      "Factor Analysis for Bicluster Acquisition",
      "Group Factor Analysis (Bayesian, view-specific sparsity)",
      "Multi-Omics Factor Analysis v2 (needs Python mofapy2)",
      "Sparse Generalised Canonical Correlation Analysis",
      "Joint and Individual Variation Explained",
      "Graph-guided Bayesian FA with robust adaptive shrinkage (Zhang 2024)"
    ),
    stringsAsFactors = FALSE
  )
}

#' List the factor-analysis methods supported by factomics
#'
#' @return A data frame with one row per method: its name (as passed to
#'   `methods` in [run_factor_pipeline()]), the R package that implements it,
#'   its family, a short description and whether that package is installed.
#' @examples
#' list_factor_methods()
#' @export
list_factor_methods <- function() {
  reg <- .method_registry()
  reg$installed <- vapply(reg$package, requireNamespace, logical(1),
                          quietly = TRUE)
  reg
}


#' Run multi-omics factor analysis pipeline
#'
#' Runs one or more latent factor methods on a (list of) numeric
#' matrix/matrices and returns scores and loadings in a single, unified
#' format, so that any downstream function (benchmarking, embeddings,
#' clustering, plots) works the same way whatever method produced the factors.
#'
#' @param X Either a single numeric matrix (**features x samples**) or a named
#'   list of such matrices.  All matrices must share the same sample names
#'   (column names).
#' @param methods Character vector of method names to run, case-insensitive.
#'   Any subset of
#'   `c("PCA", "MFA", "FABIA", "GFA", "MOFA", "SGCCA", "JIVE", "RABFA")`;
#'   see [list_factor_methods()]. Defaults to `c("FABIA", "MOFA", "GFA", "MFA")`.
#' @param k Integer.  Target number of latent factors (default `10`).
#' @param scale_features Logical.  If `TRUE` (default), each feature is
#'   z-scored before fitting. Features with zero variance become 0.
#' @param mfa_group_sizes Integer vector.  Column group sizes passed to
#'   [FactoMineR::MFA()].  If `NULL` (default), derived from view sizes.
#' @param fabia_alpha Numeric.  Sparseness parameter for FABIA (default `0.01`).
#' @param fabia_cyc Integer.  Number of FABIA cycles (default `1000`).
#' @param mofa_outfile Character.  Path for the MOFA2 HDF5 output file.
#'   Defaults to a file in [tempdir()].
#' @param mofa_convergence Character. Convergence mode for MOFA2: `"fast"`,
#'   `"medium"`, or `"slow"` (default `"slow"`).
#' @param mofa_maxiter Integer.  Maximum MOFA2 iterations (default `2000`).
#' @param mofa_use_basilisk Logical. Let MOFA2 install its own Python
#'   environment through basilisk (default `FALSE`, i.e. use the Python found
#'   by reticulate, which must have `mofapy2` installed).
#' @param mofa_python Optional path to the Python binary that has `mofapy2`
#'   installed; passed to [reticulate::use_python()] before MOFA2 runs.
#' @param gfa_iter_burnin Integer.  GFA burn-in iterations (default `500`).
#' @param gfa_iter_max Integer.  GFA maximum iterations (default `1000`).
#' @param pca_block_weight Logical. For `"PCA"`, divide each view by the square
#'   root of its number of features so that large views do not dominate
#'   (default `TRUE`).
#' @param sgcca_keep Optional integer vector (one value per view) giving the
#'   number of features kept per component in SGCCA. `NULL` (default) keeps
#'   all features (no sparsity).
#' @param jive_rank_individual Integer. Rank of each view's individual
#'   structure for JIVE (default `1`).
#' @param rabfa_graph Optional feature graph for RABFA (edge data frame,
#'   adjacency matrix, or per-view list); see [rabfa()] and
#'   [build_feature_graph()]. `NULL` (default) fits RABFA without a graph.
#' @param rabfa_iter,rabfa_burnin Total MCMC iterations and burn-in for RABFA
#'   (defaults `2000` and `1000`).
#' @param rabfa_control List of RABFA hyperparameters; see [rabfa_control()].
#' @param seed Integer.  Random seed (default `21`).
#' @param verbose Logical.  Print progress messages (default `TRUE`).
#' @param stop_on_error Logical. If `FALSE` (default) a method that fails is
#'   skipped with a warning and the other methods are still returned.
#'
#' @return A named list of class `factomics_results`, with one entry per fitted
#'   method.  Each entry contains:
#' \describe{
#'   \item{`scores`}{Matrix of sample factor scores (**samples x k**).}
#'   \item{`loadings`}{Matrix of feature loadings (**features x k**).
#'     Row names use `"view::feature"` notation.}
#'   \item{`meta`}{Method-specific metadata (number of factors, settings, etc.).}
#'   \item{`plots`}{(MOFA only) List of diagnostic ggplot objects.}
#' }
#' The preprocessed views are kept in `attr(x, "views")` and are used by
#' [variance_explained()] and [embed_results()].
#'
#' @examples
#' sim <- simulate_multiomics(n_samples = 40, seed = 1)
#' res <- run_factor_pipeline(sim$data, methods = c("PCA", "MFA"), k = 3,
#'                            verbose = FALSE)
#' res
#' head(res$MFA$scores)
#'
#' @export
run_factor_pipeline <- function(
    X,
    methods              = c("FABIA", "MOFA", "GFA", "MFA"),
    k                    = 10,
    scale_features       = TRUE,
    mfa_group_sizes      = NULL,
    fabia_alpha          = 0.01,
    fabia_cyc            = 1000,
    mofa_outfile         = file.path(tempdir(), "mofa_model.hdf5"),
    mofa_convergence     = "slow",
    mofa_maxiter         = 2000,
    mofa_use_basilisk    = FALSE,
    mofa_python          = NULL,
    gfa_iter_burnin      = 500,
    gfa_iter_max         = 1000,
    pca_block_weight     = TRUE,
    sgcca_keep           = NULL,
    jive_rank_individual = 1,
    rabfa_graph          = NULL,
    rabfa_iter           = 2000,
    rabfa_burnin         = 1000,
    rabfa_control        = list(),
    seed                 = 21,
    verbose              = TRUE,
    stop_on_error        = FALSE
) {
  methods <- unique(toupper(trimws(methods)))
  V <- .prepare_views(X, scale_features)
  sample_names <- colnames(V[[1]])

  known   <- .method_registry()$method
  unknown <- setdiff(methods, known)
  if (length(unknown))
    warning("Unknown method(s) ignored: ", paste(unknown, collapse = ", "),
            ". See list_factor_methods().", call. = FALSE)
  methods <- intersect(methods, known)

  fitters <- list(
    FABIA = function() .fit_fabia(V, k, fabia_alpha, fabia_cyc, seed, verbose),
    MOFA  = function() .fit_mofa(V, k, sample_names, mofa_outfile,
                                 mofa_convergence, mofa_maxiter,
                                 mofa_use_basilisk, mofa_python, seed),
    GFA   = function() .fit_gfa(V, k, gfa_iter_burnin, gfa_iter_max, seed,
                                verbose),
    MFA   = function() .fit_mfa(V, k, mfa_group_sizes, seed),
    PCA   = function() .fit_pca(V, k, pca_block_weight),
    SGCCA = function() .fit_sgcca(V, k, sgcca_keep),
    JIVE  = function() .fit_jive(V, k, jive_rank_individual, seed),
    RABFA = function() .fit_rabfa(V, k, rabfa_graph, rabfa_iter, rabfa_burnin,
                                  rabfa_control, seed, verbose)
  )

  results <- list()
  for (m in methods) {
    pkg <- .method_registry()$package[.method_registry()$method == m]
    if (pkg != "factomics") .stop_pkg(pkg)
    if (verbose) message("[factomics] Running ", m, " ...")
    fit <- if (stop_on_error) fitters[[m]]() else
      tryCatch(fitters[[m]](), error = function(e) {
        warning("[factomics] ", m, " failed and was skipped: ",
                conditionMessage(e), call. = FALSE)
        NULL
      })
    if (!is.null(fit)) results[[m]] <- .finalise_fit(fit, sample_names)
  }

  if (verbose)
    message("[factomics] Done. Methods fit: ",
            if (length(results)) paste(names(results), collapse = ", ") else "none")
  structure(results, class = c("factomics_results", "list"),
            views = V, sample_names = sample_names)
}


#' @export
print.factomics_results <- function(x, ...) {
  V <- attr(x, "views")
  cat("<factomics_results>\n")
  if (!is.null(V)) {
    cat(sprintf("  %d samples, %d view(s): %s\n", ncol(V[[1]]), length(V),
                paste(sprintf("%s (%d)", names(V), vapply(V, nrow, 1L)),
                      collapse = ", ")))
  }
  for (m in names(x)) {
    sc <- x[[m]]$scores
    k  <- if (is.matrix(sc)) ncol(sc) else NA
    cat(sprintf("  %-6s %s factors\n", m, k))
  }
  invisible(x)
}

#' @export
`[.factomics_results` <- function(x, i) {
  out <- unclass(x)[i]
  attributes(out) <- c(attributes(out)[c("names")],
                       list(class = class(x), views = attr(x, "views"),
                            sample_names = attr(x, "sample_names")))
  out
}


# ---- preprocessing ---------------------------------------------------------
.prepare_views <- function(X, scale_features = TRUE) {
  if (is.data.frame(X))
    stop("X must be a matrix or a named list of matrices, not a data.frame. ",
         "Convert with as.matrix() first.", call. = FALSE)
  if (!(is.matrix(X) || is.list(X)))
    stop("X must be a matrix or a named list of matrices.", call. = FALSE)
  V <- .as_list_views(X)
  if (!length(V)) stop("X is an empty list.", call. = FALSE)
  if (is.null(names(V)) || any(!nzchar(names(V))))
    names(V) <- paste0("omic", seq_along(V))
  if (anyDuplicated(names(V)))
    stop("View names must be unique.", call. = FALSE)

  V <- lapply(V, function(m) {
    if (is.data.frame(m)) m <- as.matrix(m)
    if (!is.matrix(m) || !is.numeric(m))
      stop("Every view must be a numeric matrix (features x samples).",
           call. = FALSE)
    m
  })

  sample_names <- colnames(V[[1]])
  if (is.null(sample_names))
    stop("All matrices must have column (sample) names.", call. = FALSE)
  for (nm in names(V)) {
    if (is.null(colnames(V[[nm]])))
      stop("Matrix '", nm, "' has no column names.", call. = FALSE)
    if (!setequal(colnames(V[[nm]]), sample_names))
      stop("All views must share the same sample names.", call. = FALSE)
    if (is.null(rownames(V[[nm]])))
      rownames(V[[nm]]) <- paste0(nm, "_f", seq_len(nrow(V[[nm]])))
  }
  if (length(sample_names) < 3)
    stop("At least 3 samples are required.", call. = FALSE)

  V <- lapply(V, function(m) {
    m <- m[, sample_names, drop = FALSE]
    storage.mode(m) <- "double"
    m
  })

  if (scale_features) {
    V <- lapply(V, function(m) {
      m[] <- t(scale(t(m)))
      m[!is.finite(m)] <- 0
      m
    })
  }
  V
}

.quietly <- function(expr, verbose) {
  if (verbose) return(expr)
  res <- NULL
  utils::capture.output(res <- suppressMessages(expr))
  res
}

.fcols <- function(M) {
  if (NCOL(M) > 0) colnames(M) <- paste0("F", seq_len(NCOL(M)))
  M
}

.finalise_fit <- function(fit, sample_names) {
  Z <- as.matrix(fit$scores)
  if (!is.null(rownames(Z)) && all(sample_names %in% rownames(Z)))
    Z <- Z[sample_names, , drop = FALSE]
  else rownames(Z) <- sample_names
  fit$scores   <- .fcols(Z)
  if (!is.null(fit$loadings)) fit$loadings <- .fcols(as.matrix(fit$loadings))
  fit$meta$n_factors <- ncol(fit$scores)
  fit
}


# ---- FABIA -------------------------------------------------------------------
.fit_fabia <- function(V, k, alpha, cyc, seed, verbose) {
  M <- .stack_views(V)
  set.seed(seed)
  fit <- .quietly(fabia::fabia(X = M, p = min(nrow(M), k), alpha = alpha,
                               cyc = cyc, spz = 0.5, spl = 0.1, lap = 1.0,
                               nL = 2), verbose)
  L <- as.matrix(fit@L)
  Z <- as.matrix(fit@Z)
  if (nrow(L) != nrow(M) && ncol(L) == nrow(M)) L <- t(L)
  if (nrow(Z) != ncol(M) && ncol(Z) == ncol(M)) Z <- t(Z)
  rownames(L) <- rownames(M)
  rownames(Z) <- colnames(M)
  list(scores = Z, loadings = L,
       meta = list(alpha = alpha, cyc = cyc))
}


# ---- MOFA2 -------------------------------------------------------------------
.fit_mofa <- function(V, k, sample_names, outfile, convergence, maxiter,
                      use_basilisk, python, seed) {
  if (!is.null(python)) {
    .stop_pkg("reticulate")
    reticulate::use_python(python, required = TRUE)
  }
  mofa_obj <- MOFA2::create_mofa(V)

  data_opts <- MOFA2::get_default_data_options(mofa_obj)
  data_opts$scale_views <- TRUE

  model_opts <- MOFA2::get_default_model_options(mofa_obj)
  model_opts$num_factors <- k

  train_opts <- MOFA2::get_default_training_options(mofa_obj)
  train_opts$seed             <- seed
  train_opts$convergence_mode <- convergence
  train_opts$maxiter          <- maxiter

  mofa_obj <- MOFA2::prepare_mofa(mofa_obj, data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts)
  fit <- MOFA2::run_mofa(mofa_obj, outfile, use_basilisk = use_basilisk)

  Z <- as.matrix(MOFA2::get_factors(fit, factors = "all",
                                    as.data.frame = FALSE)[[1]])
  W <- MOFA2::get_weights(fit, views = "all", factors = "all",
                          as.data.frame = FALSE)
  L <- .stack_views(W)

  varexp <- tryCatch(MOFA2::calculate_variance_explained(fit),
                     error = function(e) NULL)
  p_var  <- tryCatch(MOFA2::plot_variance_explained(fit, max_r2 = 15),
                     error = function(e) NULL)

  list(scores = Z, loadings = L,
       meta  = list(trained_hdf5 = outfile, variance_explained = varexp),
       plots = list(variance_by_factor = p_var))
}


# ---- GFA ---------------------------------------------------------------------
.fit_gfa <- function(V, k, burnin, iter_max, seed, verbose) {
  V_gfa <- lapply(V, t)  # samples x features
  opts <- GFA::getDefaultOpts()
  opts$iter.burnin <- burnin
  opts$iter.max    <- iter_max
  opts$verbose     <- if (verbose) 1 else 0
  norm <- GFA::normalizeData(V_gfa, type = "center")
  set.seed(seed)
  fit  <- .quietly(GFA::gfa(norm$train, K = k, opts = opts), verbose)

  Z <- if (is.list(fit$X)) fit$X[[1]] else fit$X
  Z <- as.matrix(Z)
  rownames(Z) <- rownames(V_gfa[[1]])

  view_feat <- vapply(V_gfa, ncol, integer(1))
  W <- if (is.list(fit$W)) do.call(rbind, lapply(fit$W, as.matrix)) else
    as.matrix(fit$W)
  if (nrow(W) != sum(view_feat) && ncol(W) == sum(view_feat)) W <- t(W)
  if (nrow(W) != sum(view_feat))
    stop("GFA returned loadings of unexpected shape.", call. = FALSE)
  rownames(W) <- unlist(lapply(names(V), function(v)
    paste0(v, "::", rownames(V[[v]]))), use.names = FALSE)

  list(scores = Z, loadings = W,
       meta = list(iter.max = iter_max, iter.burnin = burnin, K = k))
}


# ---- MFA ---------------------------------------------------------------------
.fit_mfa <- function(V, k, group_sizes, seed) {
  if (length(V) < 2)
    stop("MFA needs at least two views.", call. = FALSE)
  Vt <- lapply(names(V), function(v) {
    M <- t(V[[v]])
    colnames(M) <- paste0(v, "::", colnames(M))
    sds <- apply(M, 2, stats::sd)
    bad <- which(!is.finite(sds) | sds == 0)
    if (length(bad))
      M[, bad] <- M[, bad, drop = FALSE] +
        matrix(stats::rnorm(nrow(M) * length(bad), 0, 1e-8), nrow(M))
    M
  })
  X_mfa <- do.call(cbind, Vt)
  if (is.null(group_sizes)) group_sizes <- vapply(Vt, ncol, 1L)
  if (sum(group_sizes) != ncol(X_mfa))
    stop("sum(mfa_group_sizes) must equal the total number of features.",
         call. = FALSE)

  set.seed(seed)
  fit <- FactoMineR::MFA(as.data.frame(X_mfa, check.names = FALSE),
                         group = group_sizes,
                         type  = rep("c", length(group_sizes)),
                         name.group = names(V), ncp = k, graph = FALSE)
  Z <- as.matrix(fit$ind$coord)
  L <- as.matrix(fit$quanti.var$coord)
  L <- L[colnames(X_mfa), , drop = FALSE]
  eig <- fit$eig
  list(scores = Z, loadings = L,
       meta = list(eigenvalues = eig,
                   variance_percent = eig[seq_len(ncol(Z)), 2]))
}


# ---- PCA (block-weighted, concatenated) -------------------------------------
.fit_pca <- function(V, k, block_weight) {
  Vw <- if (block_weight) lapply(V, function(m) m / sqrt(nrow(m))) else V
  M  <- t(.stack_views(Vw))                 # samples x features
  k  <- min(k, nrow(M) - 1, ncol(M))
  pc <- stats::prcomp(M, center = TRUE, scale. = FALSE, rank. = k)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  list(scores = pc$x[, seq_len(k), drop = FALSE],
       loadings = pc$rotation[, seq_len(k), drop = FALSE],
       meta = list(variance_percent = 100 * ve[seq_len(k)],
                   block_weight = block_weight))
}


# ---- SGCCA (mixOmics) --------------------------------------------------------
.fit_sgcca <- function(V, k, keep) {
  if (length(V) < 2)
    stop("SGCCA needs at least two views.", call. = FALSE)
  Xs <- lapply(V, t)
  k  <- min(k, nrow(Xs[[1]]) - 1)
  keepX <- if (is.null(keep)) lapply(Xs, function(x) rep(ncol(x), k)) else
    stats::setNames(lapply(seq_along(Xs), function(i) rep(keep[i], k)),
                    names(Xs))
  fit <- mixOmics::wrapper.sgcca(X = Xs, ncomp = k, keepX = keepX,
                                 scale = FALSE)
  # one variate matrix per block; the consensus score is their average
  Z <- Reduce(`+`, lapply(fit$variates, as.matrix)) / length(fit$variates)
  rownames(Z) <- rownames(Xs[[1]])
  L <- .stack_views(lapply(fit$loadings, function(l) t(t(as.matrix(l)))))
  L <- L[unlist(lapply(names(V), function(v)
    paste0(v, "::", rownames(V[[v]]))), use.names = FALSE), , drop = FALSE]
  list(scores = Z, loadings = L,
       meta = list(keepX = keepX,
                   block_variates = fit$variates))
}


# ---- JIVE (r.jive) -----------------------------------------------------------
.fit_jive <- function(V, k, rank_individual, seed) {
  set.seed(seed)
  utils::capture.output(
    fit <- r.jive::jive(V, rankJ = k,
                        rankA = rep(rank_individual, length(V)),
                        method = "given", scale = TRUE, center = TRUE,
                        showProgress = FALSE)
  )
  J  <- do.call(rbind, fit$joint)            # features x samples
  sv <- svd(J, nu = k, nv = k)
  Z  <- sv$v[, seq_len(k), drop = FALSE] %*% diag(sv$d[seq_len(k)], k)
  rownames(Z) <- colnames(V[[1]])
  L  <- sv$u[, seq_len(k), drop = FALSE]
  rownames(L) <- rownames(.stack_views(V))
  list(scores = Z, loadings = L,
       meta = list(rank_joint = k, rank_individual = fit$rankA,
                   individual = fit$individual))
}


#' Build a loading comparison object aligned across methods
#'
#' Takes the `results` list returned by [run_factor_pipeline()] and
#' returns a named list of loading matrices with a common (union) feature set,
#' padded with zeros where a method does not have a particular feature.
#'
#' @param results Named list as returned by [run_factor_pipeline()].
#'
#' @return Named list (one entry per method) of **features x k** matrices
#'   sharing the same row-name order.
#'
#' @export
build_loading_comparison <- function(results) {
  stopifnot(length(results) >= 1)
  results  <- unclass(results)
  feat_all <- unique(unlist(lapply(results, function(r) rownames(r$loadings))))
  lapply(results, function(r) {
    L <- r$loadings
    M <- matrix(0, nrow = length(feat_all), ncol = ncol(L),
                dimnames = list(feat_all, colnames(L)))
    M[rownames(L), ] <- L
    M
  })
}
