# Sensitivity of ensemble FABIA to its parameters, on 3 simulated data sets.
# Run from the package root: Rscript data-raw/fabia-ensemble-sensitivity.R
# Writes vignettes/articles/fabia-ensemble-sensitivity.csv
suppressMessages(library(factomics))
eval_fit <- function(fit, sim) {
  Z <- fit$scores
  C <- abs(cor(sim$truth$Z, Z))
  found <- sum(apply(C, 1, max) > 0.8)
  rec <- mean(apply(C, 1, max))
  truth <- sim$truth$feature_labels
  if (!is.null(fit$meta$biclusters)) {
    feats <- unique(unlist(lapply(fit$meta$biclusters, `[[`, "features")))
  } else {               # plain FABIA: features of its own biclusters
    feats <- fit$bic_features
  }
  sel <- names(truth) %in% feats
  c(n_factors = ncol(Z), true_found = found, recovery = rec,
    precision = if (any(sel)) mean(truth[sel] == 1) else NA,
    recall = sum(sel & truth == 1) / sum(truth == 1))
}
base <- list(runs = 30, sample_fraction = 1, similarity = 0.5,
             min_reliability = 0.5, membership = 0.5, thresZ = 0.5)
grid <- list(
  list(par = "default"),
  list(par = "runs", val = 10), list(par = "runs", val = 100),
  list(par = "sample_fraction", val = 0.8), list(par = "sample_fraction", val = 0.6),
  list(par = "sample_fraction 0.8 + similarity", val = 0.3),
  list(par = "similarity", val = 0.3), list(par = "similarity", val = 0.7),
  list(par = "min_reliability", val = 0.3), list(par = "min_reliability", val = 0.8),
  list(par = "membership", val = 0.2), list(par = "membership", val = 1),
  list(par = "thresZ", val = 0.3), list(par = "thresZ", val = 0.8),
  list(par = "alpha", val = 0.05), list(par = "alpha", val = 0.1),
  list(par = "k", val = 3), list(par = "k", val = 6))
out <- list()
for (ds in 1:3) {
  sim <- simulate_multiomics(n_samples = 60, seed = ds)
  for (g in grid) {
    ctl <- base; alpha <- 0.01; k <- 4
    if (g$par %in% names(ctl)) ctl[[g$par]] <- g$val
    if (grepl("^sample_fraction 0.8", g$par)) {
      ctl$sample_fraction <- 0.8; ctl$similarity <- g$val
    }
    if (g$par == "alpha") alpha <- g$val
    if (g$par == "k") k <- g$val
    for (m in c("FABIA", "FABIA_ENS")) {
      if (m == "FABIA" && !(g$par %in% c("default", "alpha", "k"))) next
      t0 <- Sys.time()
      res <- tryCatch(run_factor_pipeline(sim$data, methods = m, k = k,
                                          fabia_alpha = alpha,
                                          fabia_ens_control = ctl, seed = 21,
                                          verbose = FALSE, stop_on_error = TRUE),
                      error = function(e) e)
      tt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (m == "FABIA" && !inherits(res, "error")) {
        # the biclusters of this single run, extracted as in the ensemble
        M <- factomics:::.stack_views(attr(res, "views"))
        one <- factomics:::.fabia_run(M, k, alpha, 1000, 21, seq_len(ncol(M)),
                                      0.5, NULL)
        res$FABIA$bic_features <- rownames(M)[unique(unlist(
          lapply(one$bics, `[[`, "rows")))]
      }
      v <- if (inherits(res, "error")) c(n_factors = 0, true_found = 0,
        recovery = NA, precision = NA, recall = NA) else eval_fit(res[[m]], sim)
      out[[length(out) + 1]] <- data.frame(dataset = ds, method = m,
        param = g$par, value = if (is.null(g$val)) "" else as.character(g$val),
        t(v), seconds = tt)
    }
  }
  message("dataset ", ds, " done")
}
d <- do.call(rbind, out)
write.csv(d, "vignettes/articles/fabia-ensemble-sensitivity.csv",
          row.names = FALSE)
agg <- aggregate(cbind(n_factors, true_found, recovery, precision, recall, seconds) ~
                 method + param + value, d, mean, na.action = na.pass)
agg <- agg[order(agg$method, match(agg$param, sapply(grid, `[[`, "par"))), ]
print(format(agg, digits = 3), row.names = FALSE)
