# ============================================================
#  factomics — RABFA: Robust and Adaptive-shrinkage Bayesian
#  Factor Analysis with graph information (Zhang et al., 2024)
# ============================================================
#
# Model (Gaussian data, features x samples X, p features, n samples):
#   x_ji      ~ N(m_j + w_j' z_i, 1 / rho_j)
#   z_i       ~ N(0, I_L)
#   m_j       ~ N(0, sigma_m^2)
#   rho_j     ~ Gamma(phi / 2, phi / 2)
#   w_jl      ~ Laplace(lambda_jl)          (as N(0, tau_jl), tau_jl ~ Exp(lambda_jl^2 / 2))
#   alpha_l   = log(lambda_.l) ~ N(nu1 1, nu2 Omega^{-1})
#   Omega     ~ Wishart-type prior, density proportional to
#               |Omega|^{eta (1 + eps) / 2} exp(-eta / 2 tr((11' + eps I) Omega)),
#               with omega_jk = 0 unless (j, k) is an edge of the feature graph.
#
# Sampled by Gibbs, with Metropolis-Hastings steps for alpha and a
# column-wise block Gibbs update of Omega on each connected component of
# the graph.
#
# Reference: Zhang Q, Chang C, Shen L, Long Q (2024). Incorporating graph
# information in Bayesian factor analysis with robust and adaptive shrinkage
# priors. Biometrics 80(1), ujad014. doi:10.1093/biomtc/ujad014
# This is an independent implementation written from the published model.

#' Default control settings for RABFA
#'
#' @param nu1 Prior mean of the log shrinkage parameters; `exp(nu1)` is the
#'   typical Laplace rate. Larger values shrink loadings more (default `0`).
#' @param nu2 Prior variance scale of the log shrinkage parameters. It sets how
#'   much each loading's shrinkage can adapt to the data (default `1`).
#' @param eta Concentration of the graph prior: confidence in the graph
#'   (default `10`, as in the paper).
#' @param eps Controls the strength of the prior correlation between
#'   neighbouring shrinkage parameters (default `0.2`, as in the paper).
#' @param sigma_m2 Prior variance of the feature means (default `100`).
#' @param phi Prior parameter for the noise precisions (default `1`).
#' @param mh_sd Standard deviation of the random-walk proposal for the log
#'   shrinkage parameters (default `1`).
#' @param thin Keep every `thin`-th draw after burn-in (default `1`).
#'
#' @return A named list.
#' @export
rabfa_control <- function(nu1 = 0, nu2 = 1, eta = 10, eps = 0.2,
                          sigma_m2 = 100, phi = 1, mh_sd = 1, thin = 1) {
  list(nu1 = nu1, nu2 = nu2, eta = eta, eps = eps, sigma_m2 = sigma_m2,
       phi = phi, mh_sd = mh_sd, thin = thin)
}


#' Build a feature graph for RABFA from the data
#'
#' RABFA uses a graph of features (genes, proteins, metabolites) that are
#' expected to act together, usually from a pathway database such as KEGG.
#' When no such knowledge is available, this helper builds a data-driven graph
#' by connecting features within the same view whose absolute correlation is
#' at least `cutoff`, keeping at most `max_neighbors` strongest edges per
#' feature. Prefer a knowledge-based graph when you have one: a graph learnt
#' from the same data partly re-uses the signal.
#'
#' @param X Named list of **features x samples** matrices (or one matrix).
#' @param cutoff Minimum absolute Pearson correlation for an edge
#'   (default `0.5`).
#' @param max_neighbors Maximum number of edges kept per feature
#'   (default `5`).
#'
#' @return A data frame of edges with columns `from` and `to`, using
#'   `"view::feature"` names, ready for `rabfa_graph`.
#' @examples
#' sim <- simulate_multiomics(n_samples = 30, seed = 1)
#' g <- build_feature_graph(sim$data, cutoff = 0.6)
#' head(g)
#' @export
build_feature_graph <- function(X, cutoff = 0.5, max_neighbors = 5) {
  V <- .prepare_views(X, scale_features = TRUE)
  out <- list()
  for (v in names(V)) {
    M <- V[[v]]
    if (nrow(M) < 2) next
    C <- abs(stats::cor(t(M)))
    diag(C) <- 0
    C[!is.finite(C)] <- 0
    keep <- matrix(FALSE, nrow(C), ncol(C))
    for (j in seq_len(nrow(C))) {
      cand <- which(C[j, ] >= cutoff)
      if (length(cand) > max_neighbors)
        cand <- cand[order(-C[j, cand])][seq_len(max_neighbors)]
      keep[j, cand] <- TRUE
    }
    keep <- keep | t(keep)
    idx <- which(keep & upper.tri(keep), arr.ind = TRUE)
    if (nrow(idx))
      out[[v]] <- data.frame(
        from = paste0(v, "::", rownames(M)[idx[, 1]]),
        to   = paste0(v, "::", rownames(M)[idx[, 2]]),
        stringsAsFactors = FALSE)
  }
  if (!length(out))
    return(data.frame(from = character(0), to = character(0)))
  res <- do.call(rbind, unname(out))
  rownames(res) <- NULL
  res
}


# ---- graph handling -----------------------------------------------------------

# Turn any accepted graph input into a list of neighbour index vectors over the
# stacked features `feats` ("view::feature"). Edges across views are dropped.
.rabfa_neighbours <- function(graph, feats, views) {
  p <- length(feats)
  nb <- vector("list", p)
  if (is.null(graph)) return(lapply(nb, function(x) integer(0)))

  short <- sub("^.*::", "", feats)
  view_of <- sub("::.*$", "", feats)
  resolve <- function(x, view = NULL) {
    x <- as.character(x)
    i <- match(x, feats)
    miss <- is.na(i)
    if (any(miss)) {
      cand <- if (is.null(view)) x[miss] else paste0(view, "::", x[miss])
      i[miss] <- match(cand, feats)
    }
    still <- is.na(i)
    if (any(still) && is.null(view)) {
      # plain feature names: accept only if unique across views
      dup <- short[duplicated(short)]
      ok  <- !(x[still] %in% dup)
      i[still][ok] <- match(x[still][ok], short)
    }
    i
  }

  edges <- NULL
  if (is.data.frame(graph)) {
    if (!all(c("from", "to") %in% names(graph)))
      stop("A graph data frame needs columns 'from' and 'to'.", call. = FALSE)
    edges <- cbind(resolve(graph$from), resolve(graph$to))
  } else if (is.matrix(graph) || inherits(graph, "Matrix")) {
    graph <- as.matrix(graph)
    if (is.null(rownames(graph)) || is.null(colnames(graph))) {
      if (nrow(graph) != p || ncol(graph) != p)
        stop("An unnamed graph matrix must be p x p, with p = ", p,
             " (all features of all views, stacked in order).", call. = FALSE)
      dimnames(graph) <- list(feats, feats)
    }
    idx <- which(graph != 0 & upper.tri(graph), arr.ind = TRUE)
    edges <- cbind(resolve(rownames(graph)[idx[, 1]]),
                   resolve(colnames(graph)[idx[, 2]]))
  } else if (is.list(graph)) {
    # one adjacency matrix per view
    el <- list()
    for (v in intersect(names(graph), views)) {
      G <- as.matrix(graph[[v]])
      if (is.null(rownames(G))) stop("Graph for view '", v,
                                     "' needs row and column names.",
                                     call. = FALSE)
      idx <- which(G != 0 & upper.tri(G), arr.ind = TRUE)
      el[[v]] <- cbind(resolve(rownames(G)[idx[, 1]], v),
                       resolve(colnames(G)[idx[, 2]], v))
    }
    edges <- if (length(el)) do.call(rbind, el) else matrix(0L, 0, 2)
  } else {
    stop("`rabfa_graph` must be NULL, an edge data frame (from, to), a ",
         "feature x feature adjacency matrix, or a named list of per-view ",
         "adjacency matrices.", call. = FALSE)
  }

  edges <- edges[stats::complete.cases(edges), , drop = FALSE]
  edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]
  cross <- view_of[edges[, 1]] != view_of[edges[, 2]]
  if (any(cross)) {
    message("[factomics] RABFA: ", sum(cross),
            " edge(s) between different views ignored (graphs are used ",
            "within each view).")
    edges <- edges[!cross, , drop = FALSE]
  }
  for (e in seq_len(nrow(edges))) {
    a <- edges[e, 1]; b <- edges[e, 2]
    nb[[a]] <- c(nb[[a]], b)
    nb[[b]] <- c(nb[[b]], a)
  }
  lapply(nb, function(x) sort(unique(as.integer(x))))
}

# connected components of the neighbour list
.components <- function(nb) {
  p <- length(nb)
  comp <- integer(p)
  cur <- 0L
  for (s in seq_len(p)) {
    if (comp[s]) next
    cur <- cur + 1L
    queue <- s
    comp[s] <- cur
    while (length(queue)) {
      v <- queue[1]; queue <- queue[-1]
      new <- nb[[v]][comp[nb[[v]]] == 0L]
      comp[new] <- cur
      queue <- c(queue, new)
    }
  }
  comp
}

# inverse-Gaussian sampler (Michael, Schucany and Haas, 1976), vectorised
.rinvgauss <- function(mu, shape) {
  n  <- length(mu)
  y  <- stats::rnorm(n)^2
  x  <- mu + mu^2 * y / (2 * shape) -
    mu / (2 * shape) * sqrt(4 * mu * shape * y + mu^2 * y^2)
  u  <- stats::runif(n)
  out <- ifelse(u <= mu / (mu + x), x, mu^2 / x)
  out[!is.finite(out) | out <= 0] <- .Machine$double.eps
  out
}


# ---- sampler -------------------------------------------------------------------

.rabfa_gibbs <- function(X, L, nb, iter, burnin, ctrl, seed, verbose) {
  set.seed(seed)
  p <- nrow(X); n <- ncol(X)
  nu1 <- ctrl$nu1; nu2 <- ctrl$nu2; eta <- ctrl$eta; eps <- ctrl$eps
  shape_xi <- 1 + (eta * (1 + eps) + L) / 2

  # ---- initial values from a truncated SVD ----
  m   <- rowMeans(X)
  Xc  <- X - m
  sv  <- svd(Xc, nu = L, nv = L)
  W   <- sv$u[, seq_len(L), drop = FALSE] %*%
    diag(sv$d[seq_len(L)] / sqrt(n), L)
  Z   <- t(sv$v[, seq_len(L), drop = FALSE]) * sqrt(n)        # L x n
  rho <- 1 / pmax(apply(Xc - W %*% Z, 1, stats::var), 1e-3)
  alpha <- matrix(nu1, p, L)
  tau   <- matrix(1, p, L)

  # ---- graph structure ----
  comp <- .components(nb)
  sizes <- tabulate(comp)
  iso  <- which(sizes[comp] == 1L)               # features with no edges
  big  <- which(sizes > 1L)
  comps <- lapply(big, function(cc) {
    idx <- which(comp == cc)
    loc <- match(seq_len(p), idx)
    list(idx = idx,
         nb  = lapply(idx, function(j) loc[nb[[j]]]),
         Om  = diag(length(idx)),
         Sig = diag(length(idx)))
  })
  omega_diag <- rep(1, p)          # diagonal of Omega, for every feature

  keep_it <- seq.int(burnin + 1L, iter, by = ctrl$thin)
  nk <- length(keep_it)
  W_draws <- array(0, c(p, L, nk))
  Z_draws <- array(0, c(L, n, nk))
  lam_sum <- matrix(0, p, L)
  loglik  <- numeric(iter)
  acc     <- 0
  k <- 0L
  t0 <- Sys.time()

  for (it in seq_len(iter)) {
    WZ <- W %*% Z

    # m_j
    prec_m <- n * rho + 1 / ctrl$sigma_m2
    mean_m <- rho * rowSums(X - WZ) / prec_m
    m <- stats::rnorm(p, mean_m, sqrt(1 / prec_m))

    # rho_j
    R   <- X - m - WZ
    rss <- rowSums(R^2)
    rho <- stats::rgamma(p, shape = (ctrl$phi + n) / 2,
                         rate = (ctrl$phi + rss) / 2)
    loglik[it] <- sum(0.5 * n * log(rho) - 0.5 * n * log(2 * pi) -
                        0.5 * rho * rss)

    # alpha (log shrinkage), random-walk Metropolis-Hastings
    logr_fun <- function(a_new, a_old, tau_j, lin) {
      2 * (a_new - a_old) - (exp(2 * a_new) - exp(2 * a_old)) * tau_j / 2 -
        (a_new - a_old) * lin / (2 * nu2)
    }
    if (length(iso)) {
      a_old <- alpha[iso, , drop = FALSE]
      a_new <- a_old + stats::rnorm(length(a_old), 0, ctrl$mh_sd)
      lin <- omega_diag[iso] * (2 * (a_old - nu1) + (a_new - a_old))
      ok  <- log(stats::runif(length(a_old))) <
        logr_fun(a_new, a_old, tau[iso, , drop = FALSE], lin)
      ok[is.na(ok)] <- FALSE
      a_old[ok] <- a_new[ok]
      alpha[iso, ] <- a_old
      acc <- acc + sum(ok)
    }
    for (cc in comps) {
      for (jj in seq_along(cc$idx)) {
        j  <- cc$idx[jj]
        nbj <- cc$nb[[jj]]
        a_old <- alpha[j, ]
        a_new <- a_old + stats::rnorm(L, 0, ctrl$mh_sd)
        off <- if (length(nbj))
          colSums(cc$Om[jj, nbj] * (alpha[cc$idx[nbj], , drop = FALSE] - nu1))
        else 0
        lin <- 2 * (cc$Om[jj, jj] * (a_old - nu1) + off) +
          cc$Om[jj, jj] * (a_new - a_old)
        ok <- log(stats::runif(L)) < logr_fun(a_new, a_old, tau[j, ], lin)
        ok[is.na(ok)] <- FALSE
        a_old[ok] <- a_new[ok]
        alpha[j, ] <- a_old
        acc <- acc + sum(ok)
      }
    }

    # Omega: isolated features (diagonal only)
    D <- alpha - nu1
    if (length(iso)) {
      A_jj <- eta * (1 + eps) + rowSums(D[iso, , drop = FALSE]^2) / nu2
      omega_diag[iso] <- stats::rgamma(length(iso), shape = shape_xi,
                                       rate = A_jj / 2)
    }
    # Omega: connected components, column-wise block Gibbs
    for (ci in seq_along(comps)) {
      cc <- comps[[ci]]
      Dc <- D[cc$idx, , drop = FALSE]
      A  <- eta * (matrix(1, length(cc$idx), length(cc$idx)) +
                     diag(eps, length(cc$idx))) + tcrossprod(Dc) / nu2
      Om <- cc$Om; Sig <- cc$Sig
      for (jj in seq_along(cc$idx)) {
        nbj <- cc$nb[[jj]]
        rest <- seq_along(cc$idx)[-jj]
        O11 <- Sig[rest, rest, drop = FALSE] -
          tcrossprod(Sig[rest, jj]) / Sig[jj, jj]       # Omega_11^{-1}
        xi <- stats::rgamma(1, shape = shape_xi, rate = A[jj, jj] / 2)
        u  <- numeric(length(rest))
        if (length(nbj)) {
          pos <- match(nbj, rest)
          C   <- O11[pos, pos, drop = FALSE]
          Ch  <- chol(A[jj, jj] * C)
          mu  <- -backsolve(Ch, forwardsolve(t(Ch), A[nbj, jj]))
          u[pos] <- mu + backsolve(Ch, stats::rnorm(length(nbj)))
        }
        Om[rest, jj] <- u
        Om[jj, rest] <- u
        pp <- O11 %*% u
        Om[jj, jj] <- xi + sum(u * pp)
        Sig[rest, rest] <- O11 + tcrossprod(pp) / xi
        Sig[rest, jj] <- -pp / xi
        Sig[jj, rest] <- -pp / xi
        Sig[jj, jj] <- 1 / xi
      }
      comps[[ci]]$Om <- Om
      comps[[ci]]$Sig <- Sig
      omega_diag[cc$idx] <- diag(Om)
    }

    # Z: all samples share the same posterior precision
    Xm <- X - m
    P  <- crossprod(W, W * rho) + diag(L)
    Pc <- chol(P)
    B  <- crossprod(W * rho, Xm)                         # L x n
    Zm <- backsolve(Pc, forwardsolve(t(Pc), B))
    Z  <- Zm + backsolve(Pc, matrix(stats::rnorm(L * n), L, n))

    # W and tau
    ZZ <- tcrossprod(Z)
    ZX <- Z %*% t(Xm)                                    # L x p
    for (j in seq_len(p)) {
      Pj <- rho[j] * ZZ
      diag(Pj) <- diag(Pj) + 1 / tau[j, ]
      Cj <- chol(Pj)
      W[j, ] <- backsolve(Cj, forwardsolve(t(Cj), rho[j] * ZX[, j]) +
                            stats::rnorm(L))
    }
    lam <- exp(alpha)
    tau <- 1 / .rinvgauss(lam / pmax(abs(W), 1e-10), lam^2)

    if (it %in% keep_it) {
      k <- k + 1L
      W_draws[, , k] <- W
      Z_draws[, , k] <- Z
      lam_sum <- lam_sum + lam
    }
    if (verbose && (it %% max(1, iter %/% 10) == 0))
      message(sprintf("[factomics] RABFA iteration %d/%d (%.0fs)", it, iter,
                      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  list(W = W_draws, Z = Z_draws, lambda = lam_sum / nk, loglik = loglik,
       accept = acc / (iter * p * L), n_edges = sum(lengths(nb)) / 2)
}


# Align the factors of every kept draw to a reference draw (permutation and
# sign), so that posterior summaries are not blurred by label switching.
.align_draws <- function(W, Z) {
  nk <- dim(Z)[3]
  ref <- Z[, , 1]
  L <- dim(Z)[1]
  for (s in seq_len(nk)) {
    C <- stats::cor(t(ref), t(Z[, , s]))
    C[!is.finite(C)] <- 0
    perm <- if (L > 1) as.integer(clue::solve_LSAP(abs(C), maximum = TRUE)) else 1L
    sg <- sign(C[cbind(seq_len(L), perm)])
    sg[sg == 0] <- 1
    Z[, , s] <- Z[perm, , s] * sg
    W[, , s] <- t(t(W[, perm, s]) * sg)
  }
  list(W = W, Z = Z)
}


.fit_rabfa <- function(V, k, graph, iter, burnin, control, seed, verbose) {
  ctrl <- utils::modifyList(rabfa_control(), control)
  if (burnin >= iter) stop("rabfa_burnin must be smaller than rabfa_iter.",
                           call. = FALSE)
  X <- .stack_views(V)                         # p x n, "view::feature" rows
  nb <- .rabfa_neighbours(graph, rownames(X), names(V))
  fit <- .rabfa_gibbs(X, L = k, nb = nb, iter = iter, burnin = burnin,
                      ctrl = ctrl, seed = seed, verbose = verbose)
  al <- .align_draws(fit$W, fit$Z)
  W_mean <- apply(al$W, c(1, 2), mean)
  Z_mean <- t(apply(al$Z, c(1, 2), mean))
  lo <- apply(al$W, c(1, 2), stats::quantile, 0.025)
  hi <- apply(al$W, c(1, 2), stats::quantile, 0.975)
  selected <- lo * hi > 0
  dimnames(W_mean) <- dimnames(selected) <- dimnames(fit$lambda) <-
    list(rownames(X), paste0("F", seq_len(k)))
  rownames(Z_mean) <- colnames(X)

  # BIC-type criterion from the paper (degrees of freedom = selected loadings)
  mu  <- W_mean %*% t(Z_mean)
  mu  <- mu + rowMeans(X - mu)
  rss <- rowSums((X - mu)^2)
  s2  <- pmax(rss / ncol(X), 1e-12)
  ll  <- sum(-0.5 * ncol(X) * (log(2 * pi * s2) + 1))
  bic <- -2 * ll + log(ncol(X)) * sum(selected)

  list(scores = Z_mean, loadings = W_mean,
       meta = list(
         lambda = fit$lambda, selected = selected,
         n_selected = sum(selected), bic = bic,
         n_edges = fit$n_edges, iter = iter, burnin = burnin,
         acceptance = fit$accept, loglik_trace = fit$loglik,
         control = ctrl,
         reference = paste("Zhang, Chang, Shen and Long (2024) Biometrics",
                           "80(1):ujad014")))
}


#' Robust and adaptive-shrinkage Bayesian factor analysis (RABFA)
#'
#' Fits the graph-guided Bayesian factor model of Zhang et al. (2024) to one or
#' more omics views. Each loading has its own Laplace shrinkage parameter; the
#' log shrinkage parameters of features that are connected in a feature graph
#' (for example genes in the same pathway) are linked through a sparse
#' precision matrix that is itself learnt, which makes the model robust to
#' noisy or wrong edges. Without a graph it reduces to adaptive-shrinkage
#' Bayesian factor analysis.
#'
#' This is the same model as `run_factor_pipeline(X, methods = "RABFA")`;
#' `rabfa()` is a convenience entry point. Inference is by MCMC: use more
#' iterations (the paper used 5,000-10,000) for final analyses.
#'
#' @param X Named list of **features x samples** matrices (or one matrix).
#' @param k Number of factors.
#' @param graph Optional feature graph: an edge data frame with columns `from`
#'   and `to` (see [build_feature_graph()]), a feature x feature adjacency
#'   matrix, or a named list with one adjacency matrix per view. Feature names
#'   can be plain (`"TP53"`) or view-prefixed (`"rna::TP53"`). Edges between
#'   different views are ignored.
#' @param iter Total MCMC iterations (default `2000`).
#' @param burnin Burn-in iterations discarded (default `1000`).
#' @param control List of hyperparameters; see [rabfa_control()].
#' @param scale_features Z-score each feature first (default `TRUE`).
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `factomics_results` object with one entry, `RABFA`, holding
#'   `scores`, `loadings` (posterior means) and `meta`, which includes the
#'   posterior mean shrinkage parameters (`lambda`), the loadings whose 95%
#'   credible interval excludes zero (`selected`), a BIC-type criterion
#'   (`bic`) for choosing `k` or hyperparameters, and the log-likelihood trace
#'   (`loglik_trace`) for checking convergence.
#'
#' @references Zhang Q, Chang C, Shen L, Long Q (2024). Incorporating graph
#'   information in Bayesian factor analysis with robust and adaptive
#'   shrinkage priors. *Biometrics* 80(1), ujad014.
#'   \doi{10.1093/biomtc/ujad014}
#' @examples
#' sim <- simulate_multiomics(n_samples = 30,
#'                            views = c(rna = 30, met = 15), seed = 1)
#' g   <- build_feature_graph(sim$data, cutoff = 0.6)
#' fit <- rabfa(sim$data, k = 3, graph = g, iter = 300, burnin = 150,
#'              verbose = FALSE)
#' fit$RABFA$meta$n_selected
#' @export
rabfa <- function(X, k = 5, graph = NULL, iter = 2000, burnin = 1000,
                  control = list(), scale_features = TRUE, seed = 21,
                  verbose = TRUE) {
  run_factor_pipeline(X, methods = "RABFA", k = k,
                      scale_features = scale_features, rabfa_graph = graph,
                      rabfa_iter = iter, rabfa_burnin = burnin,
                      rabfa_control = control, seed = seed, verbose = verbose,
                      stop_on_error = TRUE)
}
