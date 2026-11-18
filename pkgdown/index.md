<div class="mc-hero">
<img src="logo.png" class="mc-hero-logo" alt="factomics logo">
<div class="mc-hero-text">
<h1>factomics</h1>
<p>One toolbox for multi-omics factor analysis in R. Eight factor methods behind one function, one result format, and the whole downstream workflow: embeddings, clustering, benchmarking and interpretation. Also available as a ready-to-run Docker image.</p>
<p class="mc-cta">
<a class="mc-btn mc-btn-primary" href="articles/factomics.html">Get started</a>
<a class="mc-btn" href="articles/methods.html">Methods</a>
<a class="mc-btn" href="articles/docker.html">Docker</a>
<a class="mc-btn" href="reference/index.html">Reference</a>
<a class="mc-btn" href="https://github.com/lucp12891/factomics">GitHub</a>
</p>
</div>
</div>

## Why factomics

Multi-omics studies measure the same samples through several layers (transcriptome, metabolome, microbiome, ...). Latent factor methods summarise these layers into a few shared factors, but each method lives in its own package, with its own input format, output format and quirks. **factomics** gives them **one input** (a named list of *features x samples* matrices) and **one output** (`scores`, `loadings`, `meta` per method), so that every downstream function works with every method and methods can be compared side by side.

<div class="mc-cards">
<div class="mc-card"><h3>Matrix factorisation</h3><p>Block-weighted PCA, Multiple Factor Analysis and JIVE. <code>PCA</code>, <code>MFA</code>, <code>JIVE</code>.</p></div>
<div class="mc-card"><h3>Bayesian factor models</h3><p>Sparse, view-aware factors with automatic relevance. <code>MOFA</code>, <code>GFA</code>, <code>FABIA</code>.</p></div>
<div class="mc-card"><h3>Graph-guided</h3><p>Pathway knowledge guides the loadings, robust to wrong edges. <code>RABFA</code> (Zhang et al. 2024).</p></div>
<div class="mc-card"><h3>Correlation-based</h3><p>Components that maximise covariance between views. <code>SGCCA</code>.</p></div>
<div class="mc-card"><h3>Embeddings</h3><p>PCA, UMAP and t-SNE of each method's factor space. <code>embed_results()</code>.</p></div>
<div class="mc-card"><h3>Clustering</h3><p>k-means, hierarchical or PAM, with k chosen by silhouette. <code>cluster_results()</code>.</p></div>
<div class="mc-card"><h3>Benchmark and interpret</h3><p>Factor matching, variance explained, top features and covariate tests. <code>benchmark_factor_methods()</code>.</p></div>
</div>

## See it in action

Every figure below comes from the simulated example data, in a few lines of factomics.

<div class="mc-gallery">
<figure><img src="img/gallery/embedding.png" alt="UMAP of factor scores coloured by group"><figcaption>UMAP of each method's factors</figcaption></figure>
<figure><img src="img/gallery/variance.png" alt="Variance explained by each factor in each view"><figcaption>Variance explained per factor and view</figcaption></figure>
<figure><img src="img/gallery/agreement.png" alt="Adjusted Rand index between clusterings"><figcaption>Cluster agreement across methods</figcaption></figure>
<figure><img src="img/gallery/silhouette.png" alt="Silhouette width against number of clusters"><figcaption>Choosing the number of clusters</figcaption></figure>
<figure><img src="img/gallery/benchmark.png" alt="Correlation of factor scores across all methods"><figcaption>Benchmark: factor correlation across methods</figcaption></figure>
<figure><img src="img/gallery/associations.png" alt="Association of factors with sample covariates"><figcaption>Which factors track the covariates</figcaption></figure>
</div>

## 60-second quick start

```r
# install.packages("remotes")
remotes::install_github("lucp12891/factomics")

library(factomics)
sim <- simulate_multiomics(n_samples = 60, seed = 1)   # or a list of your own matrices

fx <- factomics(sim$data, metadata = sim$metadata,
                methods = c("PCA", "MFA", "FABIA", "GFA", "SGCCA", "JIVE"),
                k = 4, embeddings = c("umap", "tsne"))

plot_embedding(fx$embeddings, sim$metadata, color_by = "group")
plot_variance_explained(fx$results)
plot_clusters(fx$clusters, "agreement")
```

Or, with nothing installed but Docker:

```bash
docker run --rm -v "$PWD":/work ghcr.io/lucp12891/factomics:latest \
  factomics-cli --demo --out output
```

## Where to next

- [Get started](articles/factomics.html): the full workflow on simulated data.
- [Methods](articles/methods.html): what each method does, and how to run and tune it.
- [Embeddings and visualisation](articles/visualization.html): PCA, UMAP, t-SNE, heatmaps, loadings.
- [Clustering](articles/clustering.html): finding sample subgroups from the factors.
- [Benchmarking methods](articles/benchmarking.html): do the methods agree, and which recovers the truth?
- [Your own data](articles/own-data.html): reading files, preprocessing and metadata.
- [Docker and command line](articles/docker.html): the container and `factomics-cli`.
