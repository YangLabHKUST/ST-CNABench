FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL org.opencontainers.image.title="ST-CNABench-Clonalscope" \
      org.opencontainers.image.description="Docker image for Clonalscope and its dependencies (Seurat, rtracklayer)" \
      org.opencontainers.image.authors="HAN SHI"

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    git \
    pkg-config \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

RUN Rscript -e "options(Ncpus = parallel::detectCores()); \
    install.packages(c('BiocManager', 'devtools', 'amap')); \
    BiocManager::install(c('XVector', 'Biostrings', 'Rhtslib', 'SparseArray', 'DelayedArray', \
    'SummarizedExperiment', 'GenomicAlignments', 'rtracklayer', 'AUCell', \
    'BiocFileCache'), ask=FALSE)"

RUN Rscript -e "options(Ncpus = parallel::detectCores()); \
    install.packages(c('Matrix', 'data.table', 'matrixStats', 'ggplot2', \
    'stringr', 'caTools', 'pheatmap', 'cowplot', 'HiddenMarkov', 'RColorBrewer', \
    'Seurat', 'XML', 'restfulr'))"

RUN git clone https://github.com/seasoncloud/Clonalscope.git /opt/Clonalscope && \
    cd /opt/Clonalscope && \
    R CMD INSTALL . && \
    rm -rf /opt/Clonalscope
