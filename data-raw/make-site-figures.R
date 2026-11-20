# Generates the home-page gallery figures for the website.
# Run from the package root:  Rscript data-raw/make-site-figures.R
library(factomics)
library(ggplot2)

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
gal <- "pkgdown/assets/img/gallery"
dir.create(gal, recursive = TRUE, showWarnings = FALSE)

# The hex logo (man/figures/logo.png) is a supplied image, not generated here.

# ---- gallery ---------------------------------------------------------------
sim <- simulate_multiomics(n_samples = 60, seed = 1)
res <- run_factor_pipeline(sim$data,
                           methods = c("PCA", "MFA", "FABIA", "GFA", "SGCCA",
                                       "JIVE"),
                           k = 4, verbose = FALSE, gfa_iter_burnin = 300,
                           gfa_iter_max = 600)
save_g <- function(p, f, w = 6.5, h = 4.6)
  ggsave(file.path(gal, f), p, width = w, height = h, dpi = 130, bg = "white")

emb <- embed_results(res, type = "umap", methods = c("MFA", "GFA", "JIVE"),
                     include_data = FALSE)
save_g(plot_embedding(emb, sim$metadata, color_by = "group", ncol = 3) +
         labs(title = "UMAP of factor scores"), "embedding.png", 9, 3.6)
save_g(plot_variance_explained(res, methods = c("MFA", "FABIA", "GFA")),
       "variance.png", 8, 3.8)
cl <- cluster_results(res, k = 2:6)
truth <- setNames(sim$metadata$group, sim$metadata$sample)
save_g(plot_clusters(cl, "agreement", truth = truth), "agreement.png", 6.5, 5.5)
save_g(plot_clusters(cl), "silhouette.png")
bench <- benchmark_factor_methods(res)
save_g(bench$all_scores_heatmap + theme(axis.text = element_text(size = 6)),
       "benchmark.png", 6.5, 6.5)
save_g(plot_factor_associations(res, sim$metadata,
                                methods = c("PCA", "MFA", "GFA", "JIVE")),
       "associations.png", 9, 3.6)
# ---- ensemble FABIA ---------------------------------------------------------
res_ens <- run_factor_pipeline(sim$data, methods = "FABIA_ENS", k = 4,
                               fabia_ens_control = list(runs = 30),
                               verbose = FALSE)
save_g(plot_fabia_ensemble(res_ens) +
         labs(title = "Ensemble FABIA: bicluster reliability"),
       "ensemble.png", 6.5, 3.4)
st <- factor_stability(sim$data, methods = c("FABIA", "FABIA_ENS"), k = 4,
                       seeds = 1:5, truth = sim$truth$Z)
save_g(plot_factor_stability(st) +
         labs(title = "FABIA vs ensemble FABIA across seeds"),
       "stability.png", 6.5, 4)
message("Figures written.")
