# ============================================================
#  factomics — simulated multi-omics data
# ============================================================

#' Simulate a multi-omics data set with known latent factors
#'
#' Generates several omics views that share a small number of latent factors.
#' The first factors separate the sample groups, so a good factor method
#' should recover the groups, and only a sparse subset of features in each
#' view carries signal, so feature-level recovery can also be scored.
#'
#' @param n_samples Integer. Number of samples (default `60`).
#' @param views Named integer vector giving the number of features per view.
#' @param n_groups Integer. Number of sample groups (default `3`).
#' @param n_factors Integer. Number of true latent factors (default `3`).
#' @param n_group_factors Integer. How many of the factors are shifted by group
#'   (default `2`). The remaining factors are continuous nuisance variation.
#' @param effect Numeric. Distance of each group centroid from the origin on
#'   the group factors (default `2.5`).
#' @param sparsity Numeric in (0, 1]. Fraction of features in each view loaded
#'   on each factor (default `0.2`).
#' @param noise Numeric. Standard deviation of the Gaussian noise (default `1`).
#' @param seed Integer random seed, or `NULL` (default `21`).
#'
#' @return A list with
#' \describe{
#'   \item{`data`}{Named list of **features x samples** matrices.}
#'   \item{`metadata`}{Data frame with `sample`, `group`, `batch` and `age`.}
#'   \item{`truth`}{List with the true scores `Z`, loadings `W` (per view),
#'     and `feature_labels`, a named 0/1 vector flagging signal features
#'     (`"view::feature"` names).}
#' }
#' @examples
#' sim <- simulate_multiomics(n_samples = 30, seed = 1)
#' str(sim$data)
#' table(sim$metadata$group)
#' @export
simulate_multiomics <- function(
    n_samples       = 60,
    views           = c(transcriptome = 200, metabolome = 80, microbiome = 50),
    n_groups        = 3,
    n_factors       = 3,
    n_group_factors = 2,
    effect          = 2.5,
    sparsity        = 0.2,
    noise           = 1,
    seed            = 21
) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(names(views))) names(views) <- paste0("omic", seq_along(views))
  n_group_factors <- min(n_group_factors, n_factors)

  samples <- sprintf("S%03d", seq_len(n_samples))
  group   <- factor(paste0("G", rep_len(seq_len(n_groups), n_samples)))

  Z <- matrix(stats::rnorm(n_samples * n_factors), n_samples, n_factors,
              dimnames = list(samples, paste0("LF", seq_len(n_factors))))
  # group centroids on the group factors: evenly spaced on a circle of
  # radius `effect` in the first two group factors (well separated groups)
  centroids <- matrix(0, n_groups, n_group_factors)
  ang <- 2 * pi * (seq_len(n_groups) - 1) / n_groups
  centroids[, 1] <- effect * cos(ang)
  if (n_group_factors >= 2) centroids[, 2] <- effect * sin(ang)
  for (g in seq_len(n_groups)) {
    idx <- which(as.integer(group) == g)
    Z[idx, seq_len(n_group_factors)] <-
      Z[idx, seq_len(n_group_factors), drop = FALSE] +
      matrix(centroids[g, ], length(idx), n_group_factors, byrow = TRUE)
  }

  W <- list(); data <- list(); labels <- list()
  prefix <- c(transcriptome = "gene", metabolome = "met", microbiome = "otu",
              proteome = "prot", lipidome = "lip", methylome = "cpg")
  for (v in names(views)) {
    p   <- views[[v]]
    pre <- if (v %in% names(prefix)) prefix[[v]] else v
    feats <- sprintf("%s_%03d", pre, seq_len(p))
    Wv <- matrix(0, p, n_factors, dimnames = list(feats, colnames(Z)))
    for (f in seq_len(n_factors)) {
      on <- sample.int(p, max(1, round(sparsity * p)))
      Wv[on, f] <- stats::rnorm(length(on), mean = 0, sd = 1) +
        sign(stats::rnorm(length(on))) * 1
    }
    Xv <- Wv %*% t(Z) +
      matrix(stats::rnorm(p * n_samples, sd = noise), p, n_samples)
    dimnames(Xv) <- list(feats, samples)
    W[[v]]    <- Wv
    data[[v]] <- Xv
    labels[[v]] <- stats::setNames(as.integer(rowSums(Wv != 0) > 0),
                                   paste0(v, "::", feats))
  }

  metadata <- data.frame(
    sample = samples,
    group  = group,
    batch  = factor(rep_len(c("B1", "B2"), n_samples)),
    age    = round(stats::runif(n_samples, 30, 70)),
    stringsAsFactors = FALSE
  )

  list(data = data, metadata = metadata,
       truth = list(Z = Z, W = W,
                    feature_labels = unlist(unname(labels))))
}
