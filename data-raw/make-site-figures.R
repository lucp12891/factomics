# Generates the hex logo and the home-page gallery figures for the website.
# Run from the package root:  Rscript data-raw/make-site-figures.R
library(factomics)
library(ggplot2)

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
gal <- "pkgdown/assets/img/gallery"
dir.create(gal, recursive = TRUE, showWarnings = FALSE)

# ---- hex logo ---------------------------------------------------------------
hex <- data.frame(x = cos(seq(0, 2 * pi, length.out = 7) + pi / 6),
                  y = sin(seq(0, 2 * pi, length.out = 7) + pi / 6))
set.seed(3)
pts <- data.frame(
  x = c(rnorm(14, -0.32, 0.09), rnorm(14, 0.30, 0.09), rnorm(14, 0, 0.08)),
  y = c(rnorm(14, -0.02, 0.08), rnorm(14, -0.02, 0.08), rnorm(14, 0.26, 0.07)),
  g = rep(c("a", "b", "c"), each = 14))
bars <- data.frame(x = seq(-0.42, 0.42, length.out = 7),
                   h = c(0.30, 0.22, 0.16, 0.11, 0.08, 0.05, 0.04))
logo <- ggplot() +
  geom_polygon(data = hex, aes(x, y), fill = "#1A2A4F", colour = "#F59F00",
               linewidth = 3) +
  geom_rect(data = bars, aes(xmin = x - 0.045, xmax = x + 0.045,
                             ymin = -0.62, ymax = -0.62 + h),
            fill = "#F4EDDE", alpha = 0.9) +
  geom_point(data = pts, aes(x, y, colour = g), size = 3.2, alpha = 0.95) +
  scale_colour_manual(values = c("#E8366F", "#12B886", "#F59F00"),
                      guide = "none") +
  annotate("text", x = 0, y = 0.60, label = "factomics", colour = "#F4EDDE",
           size = 8, fontface = "bold") +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))
ggsave("man/figures/logo.png", logo, width = 3.2, height = 3.6, dpi = 200,
       bg = "transparent")

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
message("Figures written.")
