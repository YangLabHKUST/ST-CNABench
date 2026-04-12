FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL org.opencontainers.image.title="ST-CNVBench-SCEVAN" \
      org.opencontainers.image.description="Docker image for SCEVAN (v1.0.3) and yaGST" \
      org.opencontainers.image.authors="HANShi <shanav@connect.ust.hk>"

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('doParallel', 'ggplot2', 'parallelDist', 'pheatmap', 'forcats', 'dplyr', 'cluster', 'Rtsne', 'ape', 'tidytree', 'ggrepel', 'BiocManager', 'R.utils'), repos='https://cran.rstudio.com')"

RUN git clone https://github.com/miccec/yaGST.git /opt/yaGST && \
    cd /opt/yaGST && \
    R CMD INSTALL .

RUN Rscript -e 'BiocManager::install(c("XVector", "SummarizedExperiment", "DelayedArray", "SingleCellExperiment", "SparseArray", "ScaledMatrix", "beachmat", "BiocSingular", "scuttle", "ggiraph", "fgsea", "scran", "ggtree"), ask=FALSE, force=TRUE)'

RUN git clone https://github.com/AntonioDeFalco/SCEVAN.git /opt/SCEVAN && \
    cd /opt/SCEVAN && \
    R CMD INSTALL .
