FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL org.opencontainers.image.title="ST-CNVBench-CopyKAT" \
      org.opencontainers.image.description="Docker image for CopyKat caller in CNV benchmarking workflow" \
      org.opencontainers.image.authors="HANShi <shanav@connect.ust.hk>"

RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('data.table', 'parallelDist', 'dlm', 'gplots', 'RColorBrewer', 'mixtools', 'cluster', 'MCMCpack', 'plotly', 'Seurat'), repos='https://cran.rstudio.com')"

RUN git clone https://github.com/navinlabcode/copykat.git && \
    cd copykat && \
    R CMD INSTALL . && \
    cd .. && rm -rf copykat

CMD ["Rscript", "-e", "library(copykat); cat('CopyKat version:', as.character(packageVersion('copykat')), '\n')"]
