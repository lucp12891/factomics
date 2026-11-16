#!/usr/bin/env Rscript
# ================================================================
#  factomics command-line interface
#
#  Runs the full factomics workflow on CSV files and writes scores,
#  loadings, clusters, embeddings and figures to an output folder.
#
#  Each view is a CSV with features in rows and samples in columns
#  (first column = feature ID). Metadata is a CSV with one row per sample.
#
#  Example:
#    factomics-cli --view rna=data/rna.csv --view met=data/met.csv \
#                  --metadata data/meta.csv --sample-col sample \
#                  --methods PCA,MFA,FABIA,GFA --k 5 --color-by group \
#                  --out output
# ================================================================

usage <- function() {
  cat("
Usage: factomics-cli --view NAME=FILE [--view NAME=FILE ...] [options]

Options:
  --view NAME=FILE      A view (features x samples CSV). Repeat per omic.
  --metadata FILE       Sample metadata CSV (optional).
  --sample-col NAME     Sample-ID column in metadata (default: sample).
  --methods LIST        Comma-separated methods (default: PCA,MFA,FABIA,GFA).
                        Available: PCA,MFA,FABIA,GFA,MOFA,SGCCA,JIVE
  --k N                 Number of factors (default: 5).
  --embeddings LIST     Comma-separated: umap,tsne,pca (default: umap,tsne).
  --cluster METHOD      kmeans | hclust | pam | none (default: kmeans).
  --n-clusters A:B      Candidate cluster numbers (default: 2:6).
  --color-by NAME       Metadata column used to colour plots.
  --no-scale            Do not z-score features before fitting.
  --seed N              Random seed (default: 21).
  --out DIR             Output directory (default: factomics_output).
  --demo                Ignore --view and run on simulated data.
  --help                Show this help.
")
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || "--help" %in% args) { usage(); quit(status = 0) }

opt <- list(view = character(0), metadata = NULL, sample_col = "sample",
            methods = "PCA,MFA,FABIA,GFA", k = "5", embeddings = "umap,tsne",
            cluster = "kmeans", n_clusters = "2:6", color_by = NULL,
            scale = TRUE, seed = "21", out = "factomics_output", demo = FALSE)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  val <- if (i < length(args)) args[i + 1] else NA
  switch(a,
    "--view"       = { opt$view <- c(opt$view, val); i <- i + 1 },
    "--metadata"   = { opt$metadata <- val; i <- i + 1 },
    "--sample-col" = { opt$sample_col <- val; i <- i + 1 },
    "--methods"    = { opt$methods <- val; i <- i + 1 },
    "--k"          = { opt$k <- val; i <- i + 1 },
    "--embeddings" = { opt$embeddings <- val; i <- i + 1 },
    "--cluster"    = { opt$cluster <- val; i <- i + 1 },
    "--n-clusters" = { opt$n_clusters <- val; i <- i + 1 },
    "--color-by"   = { opt$color_by <- val; i <- i + 1 },
    "--seed"       = { opt$seed <- val; i <- i + 1 },
    "--out"        = { opt$out <- val; i <- i + 1 },
    "--no-scale"   = { opt$scale <- FALSE },
    "--demo"       = { opt$demo <- TRUE },
    stop("Unknown argument: ", a, call. = FALSE)
  )
  i <- i + 1
}

suppressPackageStartupMessages(library(factomics))
split <- function(x) trimws(strsplit(x, ",")[[1]])

if (opt$demo) {
  sim <- simulate_multiomics(seed = as.integer(opt$seed))
  views <- sim$data
  metadata <- sim$metadata
  if (is.null(opt$color_by)) opt$color_by <- "group"
} else {
  if (!length(opt$view)) stop("Give at least one --view NAME=FILE (or --demo).",
                               call. = FALSE)
  views <- list()
  for (v in opt$view) {
    kv <- strsplit(v, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("--view must look like NAME=FILE", call. = FALSE)
    df <- utils::read.csv(kv[2], row.names = 1, check.names = FALSE)
    views[[kv[1]]] <- as.matrix(df)
  }
  common <- Reduce(intersect, lapply(views, colnames))
  if (length(common) < 3) stop("Views share fewer than 3 samples.",
                               call. = FALSE)
  views <- lapply(views, function(m) m[, common, drop = FALSE])
  metadata <- if (!is.null(opt$metadata))
    utils::read.csv(opt$metadata, check.names = FALSE) else NULL
}

nc <- as.integer(strsplit(opt$n_clusters, ":")[[1]])
nc <- if (length(nc) == 2) seq(nc[1], nc[2]) else nc
emb_types <- if (opt$embeddings %in% c("", "none")) NULL else split(opt$embeddings)

fx <- factomics(views, metadata = metadata, sample_col = opt$sample_col,
                methods = toupper(split(opt$methods)), k = as.integer(opt$k),
                embeddings = emb_types,
                cluster_method = if (opt$cluster == "none") NULL else opt$cluster,
                n_clusters = nc, scale_features = opt$scale,
                seed = as.integer(opt$seed))
print(fx)

files <- export_results(fx, opt$out)
fig <- file.path(opt$out, "figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, w = 10, h = 6) {
  if (is.null(p)) return(invisible())
  tryCatch(ggplot2::ggsave(file.path(fig, name), p, width = w, height = h,
                           dpi = 150),
           error = function(e) message("Could not save ", name, ": ",
                                       conditionMessage(e)))
}
cb <- if (!is.null(metadata)) opt$color_by else NULL
res <- fx$results
save_plot(plot_variance_explained(res), "variance_explained.png", 12, 5)
save_plot(plot_factor_scores(res, metadata = metadata, color_by = cb,
                             sample_col = opt$sample_col),
          "factor_scores_F1_F2.png", 12, 7)
if (!is.null(fx$embeddings))
  save_plot(plot_embedding(fx$embeddings, metadata = metadata, color_by = cb,
                           sample_col = opt$sample_col),
            "embeddings.png", 14, 7)
if (!is.null(fx$clusters)) {
  save_plot(plot_clusters(fx$clusters), "cluster_silhouette.png", 8, 5)
  save_plot(plot_clusters(fx$clusters, "agreement"), "cluster_agreement.png",
            8, 7)
}
if (!is.null(metadata))
  save_plot(plot_factor_associations(res, metadata,
                                     sample_col = opt$sample_col),
            "factor_associations.png", 12, 5)
for (m in names(res))
  save_plot(plot_top_loadings(res, m, 1), paste0("top_loadings_", m, "_F1.png"),
            7, 6)
if (!is.null(fx$benchmark))
  save_plot(fx$benchmark$all_scores_heatmap, "benchmark_heatmap.png", 9, 8)

saveRDS(fx, file.path(opt$out, "factomics_result.rds"))
cat("\nResults written to", normalizePath(opt$out), "\n")
