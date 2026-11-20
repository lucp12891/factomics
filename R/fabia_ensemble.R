# ================================================================
#  Ensemble FABIA: robust biclusters from many FABIA runs
#
#  Follows the ensemble procedure for biclustering described in
#  Kasim, Shkedy, Kaiser, Hochreiter and Talloen (2016), "Applied
#  Biclustering Methods for Big and High-Dimensional Data Using R",
#  Chapter 10 (Khamiakova et al.), as implemented in the superbiclust
#  and biclust packages:
#    1. run FABIA many times (different seeds, optionally on a random
#       subsample of the samples);
#    2. extract the bicluster of every factor of every run;
#    3. compute the Jaccard similarity between all biclusters
#       (row-column combinations);
#    4. complete-linkage hierarchical clustering, cut at a similarity
#       threshold (default 0.5);
#    5. score each group by its reliability = share of runs in which it
#       was found, and keep the most reliable groups;
#    6. form a robust bicluster from each group: the features and samples
#       found in at least a given share of the group's biclusters
#       (1 = core / intersection, close to 0 = super / union).
#  For factomics, each robust group also becomes a factor: its loadings
#  and scores are the sign-aligned averages over the group's members.
# ================================================================

#' Settings for ensemble FABIA
#'
#' Collects the settings of the ensemble FABIA procedure, used by
#' [fabia_ensemble()] and by `run_factor_pipeline(methods = "FABIA_ENS")`.
#'
#' @param runs Number of FABIA runs (default `30`). The book used 100.
#' @param sample_fraction Share of samples used in each run. `1` (default)
#'   reruns FABIA on all samples with a different seed; a value below 1 draws a
#'   random subsample without replacement for each run. With subsampling a
#'   bicluster varies more between runs; use a lower `similarity` (about 0.3).
#' @param similarity Jaccard similarity at which the bicluster tree is cut
#'   (default `0.5`: biclusters in the same group share at least half of their
#'   row-column combinations).
#' @param min_reliability Minimum reliability (share of runs in which a group
#'   of biclusters was found) for a group to be kept (default `0.5`).
#' @param membership Share of a group's biclusters a feature or sample must
#'   belong to, to be part of the robust bicluster. `1` gives the core
#'   (intersection) bicluster, a small value such as `0.01` the "super"
#'   (union) bicluster, and `0.2` the threshold of the biclust package.
#'   Default `0.5`.
#' @param thresZ,thresL Thresholds on factor scores and loadings used to turn
#'   each FABIA factor into a bicluster; passed to [fabia::extractBic()]
#'   (defaults `0.5` and `NULL`, i.e. FABIA's own choice).
#'   `thresZ` matters most: below 0.5 FABIA's automatic `thresL` becomes so
#'   strict that biclusters keep very few features; 0.5 to 0.8 works well.
#'
#' @return A named list.
#' @references Khamiakova T, Kaiser S, Shkedy Z (2016). Ensemble methods and
#'   robust solutions. In Kasim A, Shkedy Z, Kaiser S, Hochreiter S, Talloen W
#'   (eds), *Applied Biclustering Methods for Big and High-Dimensional Data
#'   Using R*, Chapter 10. Chapman and Hall/CRC.
#' @examples
#' fabia_ensemble_control(runs = 50, sample_fraction = 0.8)
#' @export
fabia_ensemble_control <- function(runs = 30, sample_fraction = 1,
                                   similarity = 0.5, min_reliability = 0.5,
                                   membership = 0.5, thresZ = 0.5,
                                   thresL = NULL) {
  stopifnot(runs >= 2, sample_fraction > 0, sample_fraction <= 1,
            similarity > 0, similarity <= 1, min_reliability >= 0,
            min_reliability <= 1, membership > 0, membership <= 1)
  list(runs = as.integer(runs), sample_fraction = sample_fraction,
       similarity = similarity, min_reliability = min_reliability,
       membership = membership, thresZ = thresZ, thresL = thresL)
}


#' Ensemble FABIA: robust biclusters and factors from many FABIA runs
#'
#' FABIA results depend on the random starting values: two runs on the same
#' data can return different biclusters. The ensemble approach runs FABIA many
#' times and keeps only the biclusters found again and again, which are robust
#' to the starting values (and, with `sample_fraction < 1`, to the samples
#' used). This is the procedure of Chapter 10 of Kasim et al. (2016),
#' implemented in the superbiclust and biclust packages, applied to FABIA.
#'
#' Each robust group of biclusters also gives a factor, so the result can be
#' used like any other factomics method: its loadings and scores are the
#' sign-aligned averages of the loadings and scores of the FABIA factors in
#' the group. Factors are ordered by reliability.
#'
#' This is a shortcut for `run_factor_pipeline(X, methods = "FABIA_ENS", ...)`.
#'
#' @inheritParams run_factor_pipeline
#' @param k Number of factors (biclusters) per FABIA run, and the maximum
#'   number of robust factors returned.
#' @param control Settings, see [fabia_ensemble_control()].
#' @param alpha,cyc FABIA sparseness parameter and number of cycles, as for
#'   `"FABIA"` in [run_factor_pipeline()].
#'
#' @return A list with `scores` (samples x factors), `loadings`
#'   (features x factors) and `meta`. `meta` contains:
#' \describe{
#'   \item{`reliability`}{Reliability of each returned factor (share of runs).}
#'   \item{`biclusters`}{For each factor, the robust bicluster: its
#'     `features`, `samples` and `reliability`.}
#'   \item{`groups`}{All groups of biclusters found: size, number of runs,
#'     reliability and whether they were kept.}
#'   \item{`similarity`, `tree`}{The Jaccard similarity matrix of all
#'     biclusters and the complete-linkage tree (plot it with `plot()`).}
#'   \item{`runs`, `control`}{Settings used.}
#' }
#' @references Khamiakova T, Kaiser S, Shkedy Z (2016). Ensemble methods and
#'   robust solutions. In Kasim A et al. (eds), *Applied Biclustering Methods
#'   for Big and High-Dimensional Data Using R*, Chapter 10. CRC Press.
#'
#'   Hochreiter S et al. (2010). FABIA: factor analysis for bicluster
#'   acquisition. *Bioinformatics* 26(12), 1520-1527.
#' @seealso [plot_fabia_ensemble()], [factor_stability()]
#' @examples
#' \donttest{
#' if (requireNamespace("fabia", quietly = TRUE)) {
#'   sim <- simulate_multiomics(n_samples = 40, seed = 1)
#'   fit <- fabia_ensemble(sim$data, k = 3,
#'                         control = fabia_ensemble_control(runs = 10),
#'                         cyc = 200, verbose = FALSE)
#'   fit$meta$reliability
#' }
#' }
#' @export
fabia_ensemble <- function(X, k = 5, control = fabia_ensemble_control(),
                           alpha = 0.01, cyc = 1000, scale_features = TRUE,
                           seed = 21, verbose = TRUE) {
  .stop_pkg("fabia")
  V <- .prepare_views(X, scale_features)
  fit <- .fit_fabia_ens(V, k, alpha, cyc, control, seed, verbose)
  .finalise_fit(fit, colnames(V[[1]]))
}


# ---- one FABIA run --------------------------------------------------------
.fabia_run <- function(M, k, alpha, cyc, seed, cols, thresZ, thresL) {
  set.seed(seed)
  Ms  <- M[, cols, drop = FALSE]
  fit <- .quietly(fabia::fabia(X = Ms, p = min(nrow(Ms), k), alpha = alpha,
                               cyc = cyc, spz = 0.5, spl = 0.1, lap = 1.0,
                               nL = 2), FALSE)
  L <- as.matrix(fit@L)
  Z <- as.matrix(fit@Z)
  if (nrow(L) != nrow(Ms) && ncol(L) == nrow(Ms)) L <- t(L)
  if (nrow(Z) != ncol(Ms) && ncol(Z) == ncol(Ms)) Z <- t(Z)
  Zfull <- matrix(NA_real_, ncol(M), ncol(Z))
  Zfull[cols, ] <- Z
  bic <- .quietly(fabia::extractBic(fit, thresZ = thresZ, thresL = thresL),
                  FALSE)
  numn <- bic$numn
  bics <- lapply(seq_len(ncol(L)), function(j) {
    if (j > nrow(numn)) return(list(rows = integer(0), cols = integer(0)))
    list(rows = as.integer(numn[j, "numng"][[1]]),
         cols = cols[as.integer(numn[j, "numnp"][[1]])])
  })
  list(L = L, Z = Zfull, bics = bics, cols = cols)
}


# ---- Jaccard similarity of biclusters (row-column combinations) ----------
# R: features x B, C: samples x B (0/1 memberships). P: samples x B, 1 if the
# sample was used in the run the bicluster comes from (all 1 without
# subsampling). Two biclusters are compared on the samples used by both runs,
# so that a sample left out of a run does not count as a disagreement; without
# subsampling this is the usual Jaccard index of Kasim et al. (2016, Ch. 10).
.bicluster_jaccard <- function(R, C, P = NULL) {
  if (is.null(P)) P <- matrix(1, nrow(C), ncol(C))
  inter <- crossprod(R) * crossprod(C)
  nrow_b <- colSums(R)
  CP <- crossprod(C, P)              # [i, j] = |cols_i in samples of run j|
  size_i <- nrow_b * CP              # |A restricted to run j's samples|
  size_j <- t(size_i)
  union <- size_i + size_j - inter
  J <- ifelse(union > 0, inter / union, 0)
  diag(J) <- 1
  J
}


# ---- fitter ----------------------------------------------------------------
.fit_fabia_ens <- function(V, k, alpha, cyc, control, seed, verbose) {
  control <- utils::modifyList(fabia_ensemble_control(), as.list(control))
  M <- .stack_views(V)
  n <- ncol(M); p <- nrow(M)
  n_sub <- max(3L, round(control$sample_fraction * n))

  set.seed(seed)
  run_seeds <- sample.int(.Machine$integer.max, control$runs)
  runs <- vector("list", control$runs)
  for (r in seq_len(control$runs)) {
    set.seed(run_seeds[r])
    cols <- if (n_sub >= n) seq_len(n) else sort(sample.int(n, n_sub))
    runs[[r]] <- .fabia_run(M, k, alpha, cyc, run_seeds[r], cols,
                            control$thresZ, control$thresL)
    if (verbose && (r %% 10 == 0 || r == control$runs))
      message("[factomics] FABIA_ENS: ", r, "/", control$runs, " runs done")
  }

  # collect non-empty biclusters
  info <- do.call(rbind, lapply(seq_along(runs), function(r)
    data.frame(run = r, factor = seq_along(runs[[r]]$bics),
               n_rows = vapply(runs[[r]]$bics, function(b) length(b$rows), 1L),
               n_cols = vapply(runs[[r]]$bics, function(b) length(b$cols), 1L))))
  info <- info[info$n_rows > 0 & info$n_cols > 0, , drop = FALSE]
  if (nrow(info) < 2)
    stop("FABIA found fewer than two non-empty biclusters over all runs; ",
         "try a smaller thresZ or a larger k.", call. = FALSE)
  B <- nrow(info)
  labels <- paste0("run", info$run, "_F", info$factor)
  R <- matrix(0, p, B); C <- matrix(0, n, B); P <- matrix(0, n, B)
  for (b in seq_len(B)) {
    run <- runs[[info$run[b]]]
    bc <- run$bics[[info$factor[b]]]
    R[bc$rows, b] <- 1
    C[bc$cols, b] <- 1
    P[run$cols, b] <- 1
  }
  J <- .bicluster_jaccard(R, C, P)
  dimnames(J) <- list(labels, labels)

  tree   <- stats::hclust(stats::as.dist(1 - J), method = "complete")
  groups <- stats::cutree(tree, h = 1 - control$similarity)

  gtab <- do.call(rbind, lapply(sort(unique(groups)), function(g) {
    idx <- which(groups == g)
    data.frame(group = g, size = length(idx),
               n_runs = length(unique(info$run[idx])),
               reliability = length(unique(info$run[idx])) / control$runs)
  }))
  gtab <- gtab[order(-gtab$reliability, -gtab$size), , drop = FALSE]
  keep <- gtab$group[gtab$reliability >= control$min_reliability]
  keep <- utils::head(keep, k)
  gtab$kept <- gtab$group %in% keep
  rownames(gtab) <- NULL
  if (!length(keep))
    stop("No bicluster was found in at least ",
         round(100 * control$min_reliability), "% of the runs (most reliable: ",
         round(100 * max(gtab$reliability)), "%). Lower ",
         "min_reliability or similarity in fabia_ensemble_control().",
         call. = FALSE)
  if (verbose)
    message("[factomics] FABIA_ENS: ", length(keep), " robust bicluster(s) ",
            "out of ", nrow(gtab), " groups")

  feat <- rownames(M); samp <- colnames(M)
  L_out <- matrix(0, p, length(keep), dimnames = list(feat, NULL))
  Z_out <- matrix(NA_real_, n, length(keep), dimnames = list(samp, NULL))
  bic_out <- vector("list", length(keep))
  for (i in seq_along(keep)) {
    idx <- which(groups == keep[i])
    # medoid = member most similar to the others, used as sign reference
    med <- idx[which.max(colSums(J[idx, idx, drop = FALSE]))]
    ref <- runs[[info$run[med]]]$L[, info$factor[med]]
    Ls <- Zs <- vector("list", length(idx))
    for (t in seq_along(idx)) {
      run <- runs[[info$run[idx[t]]]]
      l <- run$L[, info$factor[idx[t]]]
      s <- if (sum(l * ref) < 0) -1 else 1
      Ls[[t]] <- s * l
      Zs[[t]] <- s * run$Z[, info$factor[idx[t]]]
    }
    L_out[, i] <- rowMeans(do.call(cbind, Ls))
    Zm <- do.call(cbind, Zs)
    Z_out[, i] <- rowMeans(Zm, na.rm = TRUE)
    row_prob <- rowMeans(R[, idx, drop = FALSE])
    used <- rowSums(P[, idx, drop = FALSE])
    col_prob <- ifelse(used > 0, rowSums(C[, idx, drop = FALSE]) / used, 0)
    rel <- gtab$reliability[gtab$group == keep[i]]
    bic_out[[i]] <- list(features = feat[row_prob >= control$membership],
                         samples  = samp[col_prob >= control$membership],
                         reliability = rel, n_biclusters = length(idx))
  }
  # samples never seen in a group's runs (only with subsampling): project
  if (anyNA(Z_out)) {
    proj <- t(M) %*% L_out %*%
      solve(crossprod(L_out) + diag(1e-6, ncol(L_out)))
    Z_out[is.na(Z_out)] <- proj[is.na(Z_out)]
  }
  names(bic_out) <- paste0("F", seq_along(keep))
  reliability <- vapply(bic_out, `[[`, 1, "reliability")

  list(scores = Z_out, loadings = L_out,
       meta = list(reliability = reliability, biclusters = bic_out,
                   groups = gtab, similarity = J, tree = tree,
                   runs = control$runs, alpha = alpha, cyc = cyc,
                   control = control,
                   reference = paste("Khamiakova, Kaiser & Shkedy (2016),",
                                     "Applied Biclustering Methods, Ch. 10")))
}


#' Plot the reliability of ensemble FABIA biclusters
#'
#' Shows every group of biclusters found by ensemble FABIA, ordered by
#' reliability (the share of runs in which it was found), with the groups
#' kept as robust factors highlighted and the `min_reliability` cut-off as a
#' dashed line. Groups found in only one or a few runs are the
#' start-dependent, unstable biclusters a single FABIA run may report.
#'
#' @param x The output of [fabia_ensemble()], or a `factomics_results` object
#'   containing a `"FABIA_ENS"` fit.
#' @param top Number of groups to show (default `25`).
#' @return A ggplot object.
#' @examples
#' \donttest{
#' if (requireNamespace("fabia", quietly = TRUE)) {
#'   sim <- simulate_multiomics(n_samples = 40, seed = 1)
#'   fit <- fabia_ensemble(sim$data, k = 3, cyc = 200, verbose = FALSE,
#'                         control = fabia_ensemble_control(runs = 10))
#'   plot_fabia_ensemble(fit)
#' }
#' }
#' @export
plot_fabia_ensemble <- function(x, top = 25) {
  if (inherits(x, "factomics_results")) {
    if (!"FABIA_ENS" %in% names(x))
      stop("No 'FABIA_ENS' fit in these results.", call. = FALSE)
    x <- x[["FABIA_ENS"]]
  }
  g <- x$meta$groups
  if (is.null(g)) stop("x is not an ensemble FABIA fit.", call. = FALSE)
  g <- utils::head(g, top)
  g$label <- factor(paste0("G", g$group), levels = rev(paste0("G", g$group)))
  g$status <- ifelse(g$kept, "Robust (kept)", "Not kept")
  ggplot2::ggplot(g, ggplot2::aes(x = .data$reliability, y = .data$label,
                                  fill = .data$status)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = x$meta$control$min_reliability,
                        linetype = "dashed", colour = "grey30") +
    ggplot2::scale_fill_manual(values = c("Robust (kept)" = "#1A2A4F",
                                          "Not kept" = "#C9CED8")) +
    ggplot2::scale_x_continuous(limits = c(0, 1),
                                labels = function(v) paste0(100 * v, "%")) +
    ggplot2::labs(x = sprintf("Reliability (share of %d FABIA runs)",
                              x$meta$runs),
                  y = "Bicluster group", fill = NULL,
                  title = "Ensemble FABIA: how often each bicluster is found") +
    theme_factomics()
}


#' Stability of factor methods across random seeds
#'
#' Fits the chosen methods once per seed and measures how similar the results
#' are: for every pair of seeds the factors are matched one-to-one (Hungarian
#' algorithm on the absolute correlation of loadings) and the matched
#' correlations are averaged. A method whose result does not depend on the
#' seed has stability 1. If the true factor scores are known (for example
#' from [simulate_multiomics()]), the recovery of the true factors is reported
#' for every seed as well.
#'
#' This is how plain FABIA and ensemble FABIA (`"FABIA_ENS"`) can be compared:
#' a single FABIA run depends on its starting values, the ensemble much less.
#'
#' @inheritParams run_factor_pipeline
#' @param seeds Integer vector of seeds (at least two).
#' @param truth Optional matrix of true factor scores (samples x factors).
#' @param ... Further arguments passed to [run_factor_pipeline()].
#' @return A list of class `factomics_stability` with `pairs` (one row per
#'   method and pair of seeds), `recovery` (one row per method and seed, if
#'   `truth` is given), `summary` (per method) and `fits` (the results per
#'   seed).
#' @seealso [plot_factor_stability()]
#' @examples
#' \donttest{
#' sim <- simulate_multiomics(n_samples = 40, seed = 1)
#' st <- factor_stability(sim$data, methods = c("PCA", "MFA"), k = 3,
#'                        seeds = 1:3, truth = sim$truth$Z, verbose = FALSE)
#' st$summary
#' }
#' @export
factor_stability <- function(X, methods = c("FABIA", "FABIA_ENS"), k = 5,
                             seeds = 1:5, truth = NULL, verbose = FALSE, ...) {
  seeds <- unique(as.integer(seeds))
  if (length(seeds) < 2) stop("Give at least two seeds.", call. = FALSE)
  methods <- unique(toupper(methods))
  fits <- lapply(seeds, function(s) {
    if (verbose) message("[factomics] factor_stability: seed ", s)
    run_factor_pipeline(X, methods = methods, k = k, seed = s,
                        verbose = FALSE, ...)
  })
  names(fits) <- paste0("seed", seeds)

  match_mean <- function(A, B) {
    common <- intersect(rownames(A), rownames(B))
    C <- suppressWarnings(stats::cor(A[common, , drop = FALSE],
                                     B[common, , drop = FALSE]))
    mm <- .match_factors(C)
    if (is.null(mm)) NA_real_ else mean(abs(mm$corr), na.rm = TRUE)
  }

  pairs <- list()
  rec <- list()
  for (m in methods) {
    ok <- which(vapply(fits, function(f) m %in% names(f), logical(1)))
    if (length(ok) >= 2) {
      cmb <- utils::combn(ok, 2)
      pairs[[m]] <- data.frame(
        method = m, seed_a = seeds[cmb[1, ]], seed_b = seeds[cmb[2, ]],
        stability = apply(cmb, 2, function(ij)
          match_mean(fits[[ij[1]]][[m]]$loadings, fits[[ij[2]]][[m]]$loadings)),
        stringsAsFactors = FALSE)
    }
    if (!is.null(truth) && length(ok)) {
      rec[[m]] <- data.frame(
        method = m, seed = seeds[ok],
        n_factors = vapply(ok, function(i) ncol(fits[[i]][[m]]$scores), 1L),
        recovery = vapply(ok, function(i)
          match_mean(as.matrix(truth), fits[[i]][[m]]$scores), 1),
        stringsAsFactors = FALSE)
    }
  }
  pairs <- if (length(pairs)) do.call(rbind, pairs) else NULL
  rec   <- if (length(rec)) do.call(rbind, rec) else NULL
  if (!is.null(pairs)) rownames(pairs) <- NULL
  if (!is.null(rec)) rownames(rec) <- NULL

  summ <- if (is.null(pairs)) NULL else do.call(rbind, lapply(
    split(pairs, factor(pairs$method, levels = unique(pairs$method))),
    function(d) {
      out <- data.frame(method = d$method[1],
                        mean_stability = mean(d$stability, na.rm = TRUE),
                        min_stability  = min(d$stability, na.rm = TRUE))
      if (!is.null(rec) && d$method[1] %in% rec$method) {
        r <- rec[rec$method == d$method[1], ]
        out$mean_recovery <- mean(r$recovery, na.rm = TRUE)
        out$sd_recovery   <- stats::sd(r$recovery, na.rm = TRUE)
      }
      out
    }))
  if (!is.null(summ)) rownames(summ) <- NULL
  structure(list(pairs = pairs, recovery = rec, summary = summ, fits = fits,
                 seeds = seeds),
            class = "factomics_stability")
}

#' @export
print.factomics_stability <- function(x, ...) {
  cat("<factomics_stability> ", length(x$seeds), " seeds\n", sep = "")
  if (!is.null(x$summary)) print(x$summary, digits = 3, row.names = FALSE)
  invisible(x)
}


#' Plot the stability of factor methods across seeds
#'
#' @param x The output of [factor_stability()].
#' @param type `"stability"` (default): agreement between every pair of
#'   seeds. `"recovery"`: recovery of the true factors per seed (needs `truth`
#'   in [factor_stability()]).
#' @return A ggplot object.
#' @export
plot_factor_stability <- function(x, type = c("stability", "recovery")) {
  type <- match.arg(type)
  d <- if (type == "stability") x$pairs else x$recovery
  if (is.null(d)) stop("No ", type, " results to plot.", call. = FALSE)
  d$value <- if (type == "stability") d$stability else d$recovery
  d$method <- factor(d$method, levels = unique(d$method))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$method, y = .data$value,
                                  colour = .data$method)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5,
                         colour = "grey40") +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = factomics_palette(nlevels(d$method)),
                                 guide = "none") +
    ggplot2::coord_cartesian(ylim = c(min(0.5, min(d$value, na.rm = TRUE)), 1)) +
    ggplot2::labs(x = NULL,
                  y = if (type == "stability")
                    "Agreement between seeds (matched |r| of loadings)" else
                    "Recovery of true factors (matched |r| of scores)",
                  title = if (type == "stability")
                    "How much do the factors change with the seed?" else
                    "How well are the true factors recovered?") +
    theme_factomics()
}
