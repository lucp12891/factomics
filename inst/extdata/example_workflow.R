#!/usr/bin/env Rscript
# ================================================================
#  factomics — example workflow (real data)
#  Mirrors the colorectal-cancer integromics study (PTD38 samples).
#
#  Prerequisites:
#    library(factomics)
#    library(readxl)    # to read the xlsx data files
# ================================================================

library(factomics)
library(ggplot2)
library(readxl)
library(dplyr)

# ── 0. Load raw data ─────────────────────────────────────────────

read_omic <- function(path, sheet_data = 1) {
  df  <- as.data.frame(read_excel(path, sheet = sheet_data))
  # Assume first column is feature ID
  rownames(df) <- df[[1]]
  df  <- df[, -1, drop = FALSE]
  as.matrix(df)
}

dat_16s  <- read_omic("data/WP1_16S_integromics.xlsx")
dat_met  <- read_omic("data/WP1_metabolomics_integromics_raw.xlsx")
dat_lip  <- read_omic("data/WP1_lipidomics_integromics_raw.xlsx")
metadata <- as.data.frame(read_excel("data/WP1_metadata_integromics.xlsx"))
metadata$sample <- metadata$Mouse

# ── 1. Build input list (features × samples) ────────────────────

data_list <- list(
  microbiome = dat_16s,
  metabolome = dat_met,
  lipidome   = dat_lip
)

# ── 2.A  Run all four methods ────────────────────────────────────

results_all <- run_factor_pipeline(
  X              = data_list,
  methods        = c("FABIA", "MOFA", "GFA", "MFA"),
  k              = 3,
  scale_features = TRUE,
  seed           = 21,
  verbose        = TRUE
)

# ── 2.B  Run a single method ─────────────────────────────────────

results_mofa <- run_factor_pipeline(
  X       = data_list,
  methods = "MOFA",
  k       = 3
)

# ── 2.C  Run a custom subset ─────────────────────────────────────

results_sub <- run_factor_pipeline(
  X       = data_list,
  methods = c("FABIA", "MFA"),
  k       = 5
)

# ── 3. Inspect results ───────────────────────────────────────────

str(results_all, max.level = 2)

# Scores matrix for MOFA: samples × k
head(results_all$MOFA$scores)

# Loadings matrix: view::feature rows
head(results_all$FABIA$loadings)

# MOFA variance explained plot
print(results_all$MOFA$plots$variance_by_factor)

# ── 4. Attach sample metadata ────────────────────────────────────

results_with_meta <- attach_scores_metadata(
  results    = results_all,
  metadata   = metadata,
  sample_col = "sample"
)

# Now each method has a $scores_with_meta slot
head(results_with_meta$MOFA$scores_with_meta)

# ── 5. Scatter grid: F1 coloured by Treatment ────────────────────

p_scatter <- plot_scores_scatter_grid(
  bench      = results_with_meta,
  factor_idx = 1,
  group_col  = "Treatment",
  methods    = c("MOFA", "FABIA", "MFA", "GFA"),
  align_to   = "MOFA",
  title_text = "Factor 1 — Treatment groups"
)
print(p_scatter)

# ── 6. Lollipop grid ─────────────────────────────────────────────

p_lolli <- plot_scores_lollipop_grid(
  bench      = results_with_meta,
  factor_idx = 2,
  group_col  = "Treatment",
  methods    = c("MOFA", "FABIA", "MFA", "GFA"),
  baseline   = "method_mean",
  align_to   = "MOFA",
  title_text = "Factor 2 — Treatment groups (lollipop)"
)
print(p_lolli)

# ── 7. Benchmark: pairwise method comparison ─────────────────────

bench <- benchmark_factor_methods(
  results = results_all,
  methods = c("FABIA", "MOFA", "GFA", "MFA")
)

# Summary tables
print(bench$summary_table_scores)
print(bench$summary_table_loadings)

# Global cross-method heatmap
print(bench$all_scores_heatmap)

# Pairwise score heatmap: FABIA vs MOFA
print(bench$pair_heatmaps_scores[["FABIA|MOFA"]])

# Loading heatmap for microbiome: FABIA vs MOFA
print(bench$pair_heatmaps_loadings[["FABIA|MOFA|microbiome"]])

# Pairs plot (patchwork) for Factor 1
print(bench$pairs_scores_by_factor[["F1"]])

# ── 8. Split loadings by omic ─────────────────────────────────────

loading_df          <- data.frame(results_all$FABIA$loadings,
                                   check.names = FALSE)
loading_df$feature  <- rownames(loading_df)

split_loads <- split_loadings_by_omic(
  df       = loading_df,
  column   = "feature",
  patterns = c(microbiome = "^microbiome::",
               metabolome = "^metabolome::",
               lipidome   = "^lipidome::")
)
head(split_loads$microbiome)

# ── 9. Merge scores with metadata (long format) ──────────────────

long <- merge_scores_with_metadata(
  results    = results_all,
  metadata   = metadata,
  sample_col = "sample"
)
# long$by_method$MOFA  — one row per sample
# long$all_long        — all methods stacked, 'method' column added

# ── 10. Save outputs ─────────────────────────────────────────────

dir.create("output", showWarnings = FALSE)
ggsave("output/scatter_F1_treatment.pdf", p_scatter, width = 12, height = 10)
ggsave("output/lollipop_F2_treatment.pdf", p_lolli,  width = 12, height = 10)
ggsave("output/benchmark_heatmap.pdf",    bench$all_scores_heatmap,
       width = 10, height = 8)
write.csv(bench$summary_table_scores,
          "output/benchmark_summary_scores.csv", row.names = FALSE)
write.csv(bench$summary_table_loadings,
          "output/benchmark_summary_loadings.csv", row.names = FALSE)

# ── 11. Embeddings, clustering and interpretation ────────────────

emb <- embed_results(results_all, type = c("umap", "tsne"))
print(plot_embedding(emb, metadata = metadata, color_by = "Treatment"))

cl <- cluster_results(results_all, cluster_method = "kmeans", k = 2:6)
print(cl)
print(plot_clusters(cl, "agreement"))

print(plot_variance_explained(results_all))
print(plot_factor_associations(results_all, metadata))
print(plot_top_loadings(results_all, "MOFA", factor = 1, n = 20))

export_results(results_all, "output/tables")

message("Done — outputs written to ./output/")
