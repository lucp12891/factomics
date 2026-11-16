# ================================================================
#  OmicsFactoR Docker image
#  Base: Bioconductor (includes R + many Bioc packages)
# ================================================================
FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL maintainer="your@email.com" \
      description="OmicsFactoR — multi-omics factor analysis pipeline" \
      version="0.1.0"

# ---- system libraries ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    python3-dev \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# ---- Python / mofapy2 (needed for MOFA2) ----
RUN pip3 install --no-cache-dir mofapy2

# ---- CRAN packages ----
RUN Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c( \
      'fabia', \
      'GFA', \
      'FactoMineR', \
      'clue', \
      'ggplot2', \
      'dplyr', \
      'tidyr', \
      'reshape2', \
      'patchwork', \
      'scales', \
      'tibble', \
      'pROC', \
      'reticulate', \
      'remotes', \
      'devtools', \
      'knitr', \
      'rmarkdown' \
    ), dependencies = TRUE)"

# ---- Bioconductor packages ----
RUN Rscript -e "BiocManager::install(c('MOFA2'), ask = FALSE, update = FALSE)"

# ---- Configure reticulate to use system python ----
RUN Rscript -e "reticulate::use_python(Sys.which('python3'), required = TRUE)"

# ---- Copy and install OmicsFactoR ----
COPY . /OmicsFactoR
RUN Rscript -e "devtools::install('/OmicsFactoR', dependencies = FALSE, upgrade = 'never')"

# ---- Working directory ----
WORKDIR /work

# ---- Default: launch R ----
CMD ["R"]
