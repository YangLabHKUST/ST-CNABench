FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL org.opencontainers.image.title="ST-CNVBench-Numbat" \
      org.opencontainers.image.authors="HAN SHI"

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libx11-dev libxt-dev libcairo2-dev \
    libgsl-dev zlib1g-dev \
    git pkg-config curl bzip2 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xjO bin/micromamba > /usr/bin/micromamba && \
    chmod +x /usr/bin/micromamba && \
    /usr/bin/micromamba create -y -p /opt/conda -c conda-forge && \
    /usr/bin/micromamba install -y -p /opt/conda \
        -c bioconda -c conda-forge \
        cellsnp-lite \
        samtools \
        bcftools && \
    /usr/bin/micromamba clean --all -y

ENV PATH="/opt/conda/bin:${PATH}"

RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com'), download.file.method = 'libcurl', Ncpus = parallel::detectCores())" >> /usr/local/lib/R/etc/Rprofile.site

RUN Rscript -e "install.packages(c('BiocManager', 'textshaping', 'ragg', 'Cairo', 'png', 'data.table', 'Seurat'))"

RUN Rscript -e "BiocManager::install(c('ggtree', 'GenomicRanges', 'IRanges'), ask=FALSE)"

RUN Rscript -e "install.packages('numbat', dependencies = TRUE)"

WORKDIR /data
CMD ["R"]
