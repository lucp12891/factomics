# ============================================================
#  factomics — factor pipeline
# ============================================================

#' Run multi-omics factor analysis pipeline
#'
#' Runs one or more latent factor methods on a (list of) numeric matrix/matrices
#' and returns scores and loadings in a unified format.
#'
#' @param X Either a single numeric matrix (**features × samples**) or a named
#'   list of such matrices.  All matrices must share the same sample names
#'   (column names).
#' @param methods Character vector of method names to run.  Any subset of
#'   `c("FABIA", "MOFA", "GFA", "MFA")` (case-insensitive).  Defaults to all
#'   four.
#' @param k Integer.  Target number of latent factors (default `10`).
#' @param scale_features Logical.  If `TRUE` (default), each feature is
#'   z-scored before fitting.
#' @param mfa_group_sizes Integer vector.  Column group sizes passed to
#'   [FactoMineR::MFA()].  If `NULL` (default), derived from view sizes.
#' @param fabia_alpha Numeric.  Sparseness parameter for FABIA (default `0.01`).
#' @param fabia_cyc Integer.  Number of FABIA cycles (default `1000`).
#' @param mofa_outfile Character.  Path for the MOFA2 HDF5 output file.
#'   Defaults to `"mofa_model.hdf5"` in the current working directory.
#' @param mofa_convergence Character. Convergence mode for MOFA2: `"fast"`,
#'   `"medium"`, or `"slow"` (default `"slow"`).
#' @param mofa_maxiter Integer.  Maximum MOFA2 iterations (default `2000`).
#' @param gfa_iter_burnin Integer.  GFA burn-in iterations (default `500`).
#' @param gfa_iter_max Integer.  GFA maximum iterations (default `1000`).
#' @param seed Integer.  Random seed (default `21`).
#' @param verbose Logical.  Print progress messages (default `TRUE`).
#'
#' @return A named list with one entry per fitted method.  Each entry contains:
#' \describe{
#'   \item{`scores`}{Matrix of sample factor scores (**samples × k**).}
#'   \item{`loadings`}{Matrix of feature loadings (**features × k**).
#'     Row names use `"view::feature"` notation for multi-view inputs.}
#'   \item{`meta`}{Method-specific metadata (number of factors, settings, etc.).}
#'   \item{`plots`}{(MOFA only) List of diagnostic ggplot objects.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' mats <- list(
#'   omic1 = matrix(rnorm(200), nrow = 20,
#'                  dimnames = list(paste0("f", 1:20), paste0("s", 1:10))),
#'   omic2 = matrix(rnorm(150), nrow = 15,
#'                  dimnames = list(paste0("g", 1:15), paste0("s", 1:10)))
#' )
#' res <- run_factor_pipeline(mats, methods = c("MFA", "FABIA"), k = 3)
#' str(res, max.level = 2)
#' }
#'
#' @export
run_factor_pipeline <- function(
    X,
    methods           = c("FABIA", "MOFA", "GFA", "MFA"),
    k                 = 10,
    scale_features    = TRUE,
    mfa_group_sizes   = NULL,
    fabia_alpha       = 0.01,
    fabia_cyc         = 1000,
    mofa_outfile      = file.path(getwd(), "mofa_model.hdf5"),
    mofa_convergence  = "slow",
    mofa_maxiter      = 2000,
    gfa_iter_burnin   = 500,
    gfa_iter_max      = 1000,
    seed              = 21,
    verbose           = TRUE
) {
  methods <- toupper(trimws(methods))

  # ---- input validation ----
  if (is.data.frame(X))
    stop("X must be a matrix or a named list of matrices, not a data.frame. ",
         "Convert with as.matrix() first.", call. = FALSE)
  if (!(is.matrix(X) || is.list(X)))
    stop("X must be a matrix or a named list of matrices.", call. = FALSE)
  V <- .as_list_views(X)
  if (is.null(names(V)) || any(!nzchar(names(V))))
    names(V) <- paste0("omic", seq_along(V))

  sample_names <- colnames(V[[1]])
  if (is.null(sample_names))
    stop("All matrices must have column (sample) names.", call. = FALSE)
  for (nm in names(V)) {
    if (is.null(colnames(V[[nm]])))
      stop("Matrix '", nm, "' has no column names.", call. = FALSE)
    if (!setequal(colnames(V[[nm]]), sample_names))
      stop("All views must share the same sample names.", call. = FALSE)
  }
  # align all views to same sample order
  V <- lapply(V, function(m) m[, sample_names, drop = FALSE])

  # ---- optional feature scaling ----
  if (scale_features) {
    V <- lapply(V, function(m) {
      m[] <- t(scale(t(m)))
      m[is.na(m)] <- 0
      m
    })
  }

  results <- list()

  # ====================================================
  # FABIA
  # ====================================================
  if ("FABIA" %in% methods) {
    .stop_pkg("fabia")
    if (verbose) message("[factomics] Running FABIA ...")
    M_fabia <- .stack_views(V)
    p_fabia  <- min(nrow(M_fabia), k)
    set.seed(seed)
    fit <- fabia::fabia(
      X     = M_fabia,
      p     = p_fabia,
      alpha = fabia_alpha,
      cyc   = fabia_cyc,
      spz   = 0.5, spl = 0.1, lap = 1.0, nL = 2
    )
    L <- tryCatch(fit@L, error = function(e) NULL)
    Z <- tryCatch(fit@Z, error = function(e) NULL)

    # ensure (features x k)
    if (!is.null(L)) {
      if (nrow(L) != nrow(M_fabia) && ncol(L) == nrow(M_fabia)) L <- t(L)
      rownames(L) <- rownames(M_fabia)
      if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
    } else {
      L <- matrix(0, nrow = nrow(M_fabia), ncol = 0,
                  dimnames = list(rownames(M_fabia), NULL))
    }
    # ensure (samples x k)
    if (!is.null(Z)) {
      if (nrow(Z) != ncol(M_fabia) && ncol(Z) == ncol(M_fabia)) Z <- t(Z)
      rownames(Z) <- colnames(M_fabia)
      if (ncol(Z) > 0) colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    } else {
      Z <- matrix(0, nrow = ncol(M_fabia), ncol = 0,
                  dimnames = list(colnames(M_fabia), NULL))
    }
    results$FABIA <- list(
      scores   = Z,
      loadings = L,
      meta     = list(n_factors = ncol(Z), alpha = fabia_alpha, cyc = fabia_cyc)
    )
  }

  # ====================================================
  # MOFA2
  # ====================================================
  if ("MOFA" %in% methods) {
    .stop_pkg("MOFA2")
    if (verbose) message("[factomics] Running MOFA2 ...")
    V_full <- lapply(V, function(m) {
      M <- as.matrix(m)
      storage.mode(M) <- "double"
      if (is.null(rownames(M))) rownames(M) <- paste0("f", seq_len(nrow(M)))
      M
    })
    mofa_obj <- MOFA2::create_mofa(V_full)

    data_opts <- MOFA2::get_default_data_options(mofa_obj)
    data_opts$scale_views  <- TRUE
    data_opts$scale_groups <- TRUE

    model_opts <- MOFA2::get_default_model_options(mofa_obj)
    model_opts$num_factors <- k

    train_opts <- MOFA2::get_default_training_options(mofa_obj)
    train_opts$seed             <- seed
    train_opts$convergence_mode <- mofa_convergence
    train_opts$maxiter          <- mofa_maxiter

    mofa_obj <- MOFA2::prepare_mofa(
      object           = mofa_obj,
      data_options     = data_opts,
      model_options    = model_opts,
      training_options = train_opts
    )

    fit_trained <- MOFA2::run_mofa(mofa_obj, mofa_outfile, use_basilisk = FALSE)

    fit <- tryCatch(
      MOFA2::load_model(mofa_outfile, remove_inactive_factors = FALSE),
      error = function(e) {
        if (methods::is(fit_trained, "MOFA")) {
          warning("MOFA2: load_model() failed; using in-memory model.")
          fit_trained
        } else {
          stop("MOFA2: load_model() failed: ", conditionMessage(e), call. = FALSE)
        }
      }
    )

    Z_list <- MOFA2::get_factors(fit, factors = "all", as.data.frame = FALSE)
    Z      <- as.matrix(Z_list[[1]])
    if (!is.null(Z) && ncol(Z) > 0) {
      colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    } else {
      Z <- matrix(0, nrow = length(sample_names), ncol = 0,
                  dimnames = list(sample_names, NULL))
    }

    W_list <- MOFA2::get_weights(fit, views = "all", factors = "all",
                                  as.data.frame = FALSE)
    L <- .stack_views(W_list)
    if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))

    p_varExp <- tryCatch(
      MOFA2::plot_variance_explained(fit, max_r2 = 15),
      error = function(e) NULL
    )
    varExp_tbl <- tryCatch(
      MOFA2::calculate_variance_explained(fit),
      error = function(e) NULL
    )

    # Factor-by-group boxplot
    p_factorByGroup <- tryCatch({
      smd <- tryCatch(MOFA2::samples_metadata(fit), error = function(e) NULL)
      if (is.null(smd) || nrow(smd) == 0) {
        smd <- data.frame(sample = rownames(Z), Group = "All",
                          row.names = rownames(Z))
      } else {
        if (is.null(rownames(smd)) && "sample" %in% names(smd))
          rownames(smd) <- smd$sample
        smd <- smd[rownames(Z), , drop = FALSE]
        if (!("Group" %in% colnames(smd))) {
          cat_cols <- names(which(vapply(smd,
                                        function(x) is.factor(x) || is.character(x),
                                        logical(1))))
          smd$Group <- if (length(cat_cols) > 0) smd[[cat_cols[1]]] else "All"
        }
      }
      df_long        <- utils::stack(as.data.frame(Z))
      df_long$sample <- rep(rownames(Z), times = ncol(Z))
      names(df_long) <- c("score", "Factor", "sample")
      df_long$Group  <- smd$Group[match(df_long$sample, rownames(smd))]
      ggplot2::ggplot(df_long,
                      ggplot2::aes(x = .data[["Factor"]],
                                   y = .data[["score"]],
                                   fill = .data[["Group"]])) +
        ggplot2::geom_boxplot(outlier.size = 0.5) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "MOFA factor scores by group",
                      x = "Factor", y = "Score", fill = "Group")
    }, error = function(e) NULL)

    results$MOFA <- list(
      scores   = Z,
      loadings = L,
      meta     = list(
        n_factors          = ncol(Z),
        trained_hdf5       = mofa_outfile,
        variance_explained = varExp_tbl
      ),
      plots = list(
        variance_by_factor = p_varExp,
        factors_by_group   = p_factorByGroup
      )
    )
  }

  # ====================================================
  # GFA
  # ====================================================
  if ("GFA" %in% methods) {
    .stop_pkg("GFA")
    if (verbose) message("[factomics] Running GFA ...")

    V_gfa <- lapply(V, t)  # samples x features

    model_option             <- GFA::getDefaultOpts()
    model_option$iter.burnin <- gfa_iter_burnin
    model_option$iter.max    <- gfa_iter_max

    norm <- GFA::normalizeData(V_gfa, type = "center")
    set.seed(seed)
    fit  <- GFA::gfa(norm$train, K = k, opts = model_option)

    # ---- scores ----
    Z <- {
      Zcand <- if (is.list(fit$X) && length(fit$X) >= 1) fit$X[[1]] else fit$X
      if (is.null(Zcand)) {
        matrix(0, nrow = nrow(V_gfa[[1]]), ncol = 0,
               dimnames = list(rownames(V_gfa[[1]]), NULL))
      } else {
        Zcand <- as.matrix(Zcand)
        if (is.null(rownames(Zcand))) rownames(Zcand) <- rownames(V_gfa[[1]])
        if (NCOL(Zcand) > 0 && is.null(colnames(Zcand)))
          colnames(Zcand) <- paste0("F", seq_len(NCOL(Zcand)))
        Zcand
      }
    }

    # ---- loadings ----
    .make_empty_gfa <- function(i) {
      nf    <- ncol(V_gfa[[i]])
      feats <- colnames(V_gfa[[i]])
      if (is.null(feats)) feats <- paste0("feature_", seq_len(nf))
      M <- matrix(0, nrow = nf, ncol = 0)
      rownames(M) <- paste0(names(V_gfa)[i], "::", feats)
      M
    }
    W_list <- list()

    if (is.list(fit$W) && length(fit$W) == length(V_gfa)) {
      W_list <- lapply(seq_along(fit$W), function(i) {
        Wi    <- as.matrix(fit$W[[i]])
        if (is.null(dim(Wi))) Wi <- matrix(Wi, nrow = 1)
        n_feat <- ncol(V_gfa[[i]])
        feats  <- colnames(V_gfa[[i]])
        if (is.null(feats)) feats <- paste0("feature_", seq_len(n_feat))
        if      (nrow(Wi) == n_feat) {}
        else if (ncol(Wi) == n_feat) { Wi <- t(Wi) }
        else if (n_feat == 1 && (nrow(Wi) == 1 || ncol(Wi) == 1)) {
          if (nrow(Wi) != 1) Wi <- t(Wi)
        } else {
          return(.make_empty_gfa(i))
        }
        rownames(Wi) <- paste0(names(V_gfa)[i], "::", feats)
        if (NCOL(Wi) > 0 && is.null(colnames(Wi)))
          colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
        Wi
      })
    } else if (is.matrix(fit$W)) {
      Wi_all     <- as.matrix(fit$W)
      view_feat  <- vapply(V_gfa, ncol, integer(1))
      total_feat <- sum(view_feat)
      if      (nrow(Wi_all) == total_feat) {}
      else if (ncol(Wi_all) == total_feat) { Wi_all <- t(Wi_all) }
      else {
        warning("GFA returned unexpected W shape; returning empty loadings.")
        W_list <- lapply(seq_along(V_gfa), .make_empty_gfa)
      }
      if (!length(W_list)) {
        cuts   <- c(0, cumsum(view_feat))
        W_list <- vector("list", length(V_gfa))
        for (i in seq_along(V_gfa)) {
          rows  <- (cuts[i] + 1):cuts[i + 1]
          Wi    <- Wi_all[rows, , drop = FALSE]
          feats <- colnames(V_gfa[[i]])
          if (is.null(feats)) feats <- paste0("feature_", seq_len(nrow(Wi)))
          rownames(Wi) <- paste0(names(V_gfa)[i], "::", feats)
          if (NCOL(Wi) > 0 && is.null(colnames(Wi)))
            colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
          W_list[[i]] <- Wi
        }
      }
    } else {
      W_list <- lapply(seq_along(V_gfa), .make_empty_gfa)
    }

    L <- if (length(W_list)) do.call(rbind, W_list) else NULL

    results$GFA <- list(
      scores   = Z,
      loadings = L,
      meta     = list(
        n_factors   = NCOL(Z),
        iter.max    = gfa_iter_max,
        iter.burnin = gfa_iter_burnin,
        K           = k
      )
    )
  }

  # ====================================================
  # MFA
  # ====================================================
  if ("MFA" %in% methods) {
    .stop_pkg("FactoMineR")
    if (verbose) message("[factomics] Running MFA ...")

    Vt <- lapply(V, t)
    if (is.null(names(Vt)) || any(!nzchar(names(Vt))))
      names(Vt) <- paste0("view", seq_along(Vt))

    Vt <- lapply(seq_along(Vt), function(i) {
      M <- as.matrix(Vt[[i]])
      if (is.null(colnames(M))) colnames(M) <- paste0("f", seq_len(ncol(M)))
      colnames(M) <- paste0(names(Vt)[i], "::", colnames(M))
      M
    })
    names(Vt) <- names(V)

    # jitter zero-variance columns
    jitter_cols <- function(M, eps = 1e-8) {
      if (ncol(M) == 0) return(M)
      sds <- apply(M, 2, stats::sd, na.rm = TRUE)
      bad <- which(!is.finite(sds) | sds == 0)
      if (length(bad))
        M[, bad] <- M[, bad, drop = FALSE] +
          matrix(stats::rnorm(nrow(M) * length(bad), 0, eps), nrow(M))
      M
    }
    Vt <- lapply(Vt, jitter_cols)

    X_mfa <- do.call(cbind, Vt)
    if (is.null(mfa_group_sizes)) mfa_group_sizes <- vapply(Vt, ncol, 1L)
    stopifnot(sum(mfa_group_sizes) == ncol(X_mfa))

    X_df <- as.data.frame(X_mfa, check.names = FALSE)
    if (!is.null(rownames(X_mfa))) rownames(X_df) <- rownames(X_mfa)

    set.seed(seed)
    fit <- FactoMineR::MFA(
      X_df,
      group = mfa_group_sizes,
      type  = rep("c", length(mfa_group_sizes)),
      ncp   = k,
      graph = FALSE
    )

    Z <- fit$ind$coord
    if (is.null(Z)) {
      Z <- matrix(0, nrow = nrow(X_df), ncol = 0,
                  dimnames = list(rownames(X_df), NULL))
    } else {
      Z <- as.matrix(Z)
      colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    }

    Lraw <- fit$var$coord
    if (is.null(Lraw)) Lraw <- fit$quanti.var$coord

    if (!is.null(Lraw)) {
      dim_cols <- grep("^Dim\\.", colnames(Lraw))
      L        <- as.matrix(Lraw[, dim_cols, drop = FALSE])
      common   <- intersect(colnames(X_mfa), rownames(L))
      L        <- L[common, , drop = FALSE]
      colnames(L) <- paste0("F", seq_len(ncol(L)))
      if (length(common) < ncol(X_mfa)) {
        miss <- setdiff(colnames(X_mfa), common)
        Lpad <- matrix(0, nrow = length(miss), ncol = ncol(L),
                       dimnames = list(miss, colnames(L)))
        L <- rbind(L, Lpad)
        L <- L[colnames(X_mfa), , drop = FALSE]
      }
    } else {
      L <- matrix(0, nrow = ncol(X_mfa), ncol = 0,
                  dimnames = list(colnames(X_mfa), NULL))
    }

    results$MFA <- list(
      scores   = Z,
      loadings = L,
      meta     = list(n_factors = ncol(Z))
    )
  }

  unknown <- setdiff(methods, c("FABIA", "MOFA", "GFA", "MFA"))
  if (length(unknown))
    warning("Unknown method(s) ignored: ", paste(unknown, collapse = ", "),
            call. = FALSE)

  if (verbose)
    message("[factomics] Done. Methods fit: ",
            paste(names(results), collapse = ", "))
  results
}


#' Build a loading comparison object aligned across methods
#'
#' Takes the `results` list returned by [run_factor_pipeline()] and
#' returns a named list of loading matrices with a common (union) feature set,
#' padded with zeros where a method does not have a particular feature.
#'
#' @param results Named list as returned by [run_factor_pipeline()].
#'
#' @return Named list (one entry per method) of **features × k** matrices
#'   sharing the same row-name order.
#'
#' @export
build_loading_comparison <- function(results) {
  stopifnot(length(results) >= 1)
  feat_all <- unique(unlist(lapply(results, function(r) rownames(r$loadings))))
  lapply(results, function(r) {
    L <- r$loadings
    M <- matrix(0, nrow = length(feat_all), ncol = ncol(L),
                dimnames = list(feat_all, colnames(L)))
    M[rownames(L), ] <- L
    M
  })
}

set.seed(1235)
benchmark_res_lena2.0 <- run_factor_pipeline(
  #list(microbiome = as.matrix(dat_16s_clean), metabolome = as.matrix(dat_met_clean),lipidome = as.matrix(dat_lip_clean)),
  list(microbiome = as.matrix(t(X_16s)), metabolome = as.matrix(t(X_met)),lipidome = as.matrix(t(X_lip))),
  methods = c("fabia","mofa","gfa","mfa"),
  k = 2,
  scale_features = FALSE
)

metadata <- meta; metadata$sample <- metadata$Mouse
mergedx_lena_res2.0 <- factomics::attach_scores_metadata(benchmark_res_lena2.0, metadata, sample_col = "sample")

p_nonex <- factomics::plot_scores_scatter_grid(
  mergedx_lena_res2.0,
  factor_idx = 1,
  group_col  = "Treatment",
  methods    = c("MOFA","FABIA","MFA","GFA"),
  panel_tags = NULL,
  tag_position = "none",
  align_to = "MOFA"
)
print(p_nonex)

p_lolli_method_mean <- factomics::plot_scores_lollipop_grid(
  mergedx_lena_res2.0,
  factor_idx = 1,
  group_col  = "Treatment",
  methods    = c("MOFA","FABIA","MFA","GFA"),
  baseline   = "method_mean",
  align_to   = "FABIA"#,
  #palette    = c("10AOM"="#00BFC4", "12AOM"="#F8766D", "NA"="#999999") # optional
)
print(p_lolli_method_mean)

p_lolli_method_mean <- factomics::plot_scores_lollipop_grid(
  mergedx_lena_res2.0,
  factor_idx = 2,
  group_col  = "Treatment",
  methods    = c("MOFA","FABIA","MFA","GFA"),
  baseline   = "method_mean",
  align_to   = "FABIA"#,
  #palette    = c("10AOM"="#00BFC4", "12AOM"="#F8766D", "NA"="#999999") # optional
)
print(p_lolli_method_mean)

benchmark_res_lena_res = mergedx_lena_res2.0
bechmarked <- factomics::benchmark_factor_methods(benchmark_res_lena_res)
