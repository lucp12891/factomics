# factomics <img src="man/figures/logo.png" align="right" height="139" alt="factomics logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/lucp12891/factomics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lucp12891/factomics/actions/workflows/R-CMD-check.yaml)
[![Docker image](https://github.com/lucp12891/factomics/actions/workflows/docker.yaml/badge.svg)](https://github.com/lucp12891/factomics/actions/workflows/docker.yaml)
[![Website](https://img.shields.io/badge/docs-website-1A2A4F)](https://lucp12891.github.io/factomics/)
<!-- badges: end -->

**One toolbox for multi-omics factor analysis.** Instead of learning, installing
and reformatting data for many separate packages, `factomics` runs nine latent
factor methods through a single function, returns every result in the same
format, and adds the downstream steps you usually need: embeddings (PCA, UMAP,
t-SNE), clustering, method benchmarking, interpretation and plots. A Docker
image ships everything pre-installed, including MOFA2's Python backend.

Website and tutorials: **<https://lucp12891.github.io/factomics/>**

## Methods

| Method | Package | Family |
|--------|---------|--------|
| `PCA`   | stats      | Block-weighted PCA on concatenated views |
| `MFA`   | FactoMineR | Multiple Factor Analysis |
| `FABIA` | fabia      | Factor Analysis for Bicluster Acquisition |
| `FABIA_ENS` | fabia  | Ensemble FABIA: robust biclusters kept over many runs (Kasim et al. 2016, Ch. 10) |
| `GFA`   | GFA        | Group Factor Analysis (Bayesian) |
| `MOFA`  | MOFA2      | Multi-Omics Factor Analysis v2 |
| `SGCCA` | mixOmics   | Sparse Generalised Canonical Correlation Analysis |
| `JIVE`  | r.jive     | Joint and Individual Variation Explained |
| `RABFA` | built in   | Graph-guided Bayesian factor analysis with robust, adaptive shrinkage (Zhang et al. 2024) |

Downstream: `embed_results()` (PCA / UMAP / t-SNE), `cluster_results()`
(k-means / hierarchical / PAM, k chosen by silhouette), `compare_clusterings()`
(adjusted Rand index), `benchmark_factor_methods()`, `variance_explained()`,
`factor_associations()`, `plot_top_loadings()`, `plot_factor_heatmap()` and more.

### Robust biclusters with ensemble FABIA

A single FABIA run depends on its random start. `FABIA_ENS` runs FABIA many
times (`fabia_ens_control = list(runs = 30)`), groups the biclusters by Jaccard
similarity and keeps those found in most runs, following the ensemble method of
Chapter 10 of *Applied Biclustering Methods* (Kasim et al. 2016).
`factor_stability()` compares plain and ensemble FABIA across seeds:

```r
st <- factor_stability(sim$data, methods = c("FABIA", "FABIA_ENS"), k = 4,
                       seeds = 1:5, truth = sim$truth$Z)
plot_factor_stability(st)
```

### Graph-guided factor analysis (RABFA)

`RABFA` lets prior knowledge of which features act together (for example genes
in the same KEGG pathway) guide the factor loadings, and stays robust when that
graph contains wrong edges:

```r
# edges: data frame with columns from, to (feature names, within a view)
fit <- rabfa(sim$data, k = 4, graph = edges, iter = 5000, burnin = 2500)
# or, with no pathway database, a data-driven graph:
fit <- rabfa(sim$data, k = 4, graph = build_feature_graph(sim$data))
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("lucp12891/factomics")

# Optional: the factor methods you want to use
install.packages(c("FactoMineR", "GFA", "r.jive", "uwot", "Rtsne"))
BiocManager::install(c("fabia", "MOFA2", "mixOmics"))
```

MOFA2 also needs the Python package `mofapy2` (`pip install mofapy2`); point to
that Python with `run_factor_pipeline(..., mofa_python = "/path/to/python")`.

## Quick start

```r
library(factomics)

sim <- simulate_multiomics(n_samples = 60, seed = 1)   # or your own data
# sim$data: named list of features x samples matrices

fx <- factomics(sim$data, metadata = sim$metadata,
                methods = c("PCA", "MFA", "FABIA", "GFA", "SGCCA", "JIVE"),
                k = 4, embeddings = c("umap", "tsne"))
fx

plot_embedding(fx$embeddings, sim$metadata, color_by = "group")
plot_variance_explained(fx$results)
plot_clusters(fx$clusters, "agreement")
fx$benchmark$all_scores_heatmap
```

Step by step, the same workflow is:

```r
res <- run_factor_pipeline(sim$data, methods = c("MFA", "GFA"), k = 4)
res <- attach_scores_metadata(res, sim$metadata)
emb <- embed_results(res, type = "umap")
cl  <- cluster_results(res, cluster_method = "kmeans", k = 2:6)
bm  <- benchmark_factor_methods(res)
export_results(res, "results/")
```

## Docker

```bash
docker pull ghcr.io/lucp12891/factomics:latest

# Full analysis from CSV files (features x samples, first column = feature ID)
docker run --rm -v "$PWD":/work ghcr.io/lucp12891/factomics:latest \
  factomics-cli --view rna=data/rna.csv --view met=data/met.csv \
                --metadata data/meta.csv --color-by group \
                --methods PCA,MFA,FABIA,GFA,MOFA --k 5 --out output

# Interactive R, or RStudio Server on http://localhost:8787
docker run --rm -it -v "$PWD":/work ghcr.io/lucp12891/factomics:latest R
docker run --rm -p 8787:8787 -e PASSWORD=choose-one ghcr.io/lucp12891/factomics:latest
```

See the [Docker article](https://lucp12891.github.io/factomics/articles/docker.html)
for details.

## Citation

Please also cite the method(s) you use:

- **FABIA**: Hochreiter et al. (2010) *Bioinformatics* 26:1520-1527.
- **Ensemble FABIA**: Khamiakova, Kaiser & Shkedy (2016), Ensemble methods and
  robust solutions, in Kasim et al. (eds) *Applied Biclustering Methods for Big
  and High-Dimensional Data Using R*, Ch. 10, CRC Press.
- **MOFA2**: Argelaguet et al. (2020) *Genome Biology* 21:111.
- **GFA**: Virtanen et al. (2012) *AISTATS*; Leppaaho et al. (2017) *JMLR* 18:1-5.
- **MFA**: Escofier & Pages (1994) *Computational Statistics & Data Analysis* 18:121-140.
- **SGCCA**: Tenenhaus et al. (2014) *Biostatistics* 15:569-583.
- **JIVE**: Lock et al. (2013) *Annals of Applied Statistics* 7:523-542.
- **RABFA**: Zhang, Chang, Shen & Long (2024) *Biometrics* 80(1):ujad014.

## License

MIT © Bernard Isekah Osang'ir
