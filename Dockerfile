# ================================================================
#  factomics Docker image
#  One container with every supported factor method (PCA, MFA, FABIA,
#  GFA, MOFA2, SGCCA, JIVE), embeddings (UMAP, t-SNE) and clustering.
#
#  Base: Bioconductor image (R + RStudio Server + Bioconductor tooling)
# ================================================================
FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL org.opencontainers.image.title="factomics" \
      org.opencontainers.image.description="One toolbox for multi-omics factor analysis" \
      org.opencontainers.image.source="https://github.com/lucp12891/factomics" \
      org.opencontainers.image.licenses="MIT" \
      maintainer="Bernard Isekah Osang'ir <bernard.osangir@sckcen.be>"

# ---- system libraries ----
RUN apt-get update && apt-get install -y --no-install-recommends \
      libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
      libfontconfig1-dev libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
      libglpk-dev python3-dev python3-venv python3-pip \
    && rm -rf /var/lib/apt/lists/*

# ---- Python environment for MOFA2 (mofapy2) ----
RUN python3 -m venv /opt/mofa \
    && /opt/mofa/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/mofa/bin/pip install --no-cache-dir mofapy2
ENV RETICULATE_PYTHON=/opt/mofa/bin/python
RUN echo "RETICULATE_PYTHON=/opt/mofa/bin/python" >> /usr/local/lib/R/etc/Renviron.site

# ---- R / Bioconductor packages ----
RUN Rscript -e "BiocManager::install(c( \
      'clue', 'cluster', 'dplyr', 'ggplot2', 'patchwork', 'reshape2', \
      'scales', 'tibble', 'tidyr', \
      'fabia', 'GFA', 'FactoMineR', 'MOFA2', 'mixOmics', 'r.jive', \
      'uwot', 'Rtsne', 'pROC', 'reticulate', \
      'knitr', 'rmarkdown', 'testthat', 'remotes' \
    ), ask = FALSE, update = FALSE, Ncpus = parallel::detectCores())" \
 && Rscript -e "pk <- c('fabia','GFA','FactoMineR','MOFA2','mixOmics','r.jive','uwot','Rtsne'); \
      miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]; \
      if (length(miss)) stop('Missing packages: ', paste(miss, collapse = ', '))"

# ---- install factomics from this repository ----
COPY . /tmp/factomics
RUN R CMD INSTALL /tmp/factomics \
 && install -m 0755 /tmp/factomics/inst/scripts/factomics-cli.R /usr/local/bin/factomics-cli \
 && rm -rf /tmp/factomics

WORKDIR /work

# Default: the base image's init, which serves RStudio Server on port 8787.
# Override with `R`, `bash` or `factomics-cli ...` for command-line use.
