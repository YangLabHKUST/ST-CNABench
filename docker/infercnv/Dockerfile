FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL org.opencontainers.image.title="ST-CNVBench-inferCNV" \
      org.opencontainers.image.description="Docker image for inferCNV caller in CNV benchmarking workflow" \
      org.opencontainers.image.authors="HANShi <shanav@connect.ust.hk>"

RUN apt-get update && apt-get -y install curl libssl-dev libcurl4-openssl-dev \
                                        libxml2-dev git python3 jags \
                                        r-cran-rjags time && \
                      apt-get clean && rm -rf /var/tmp/* \
                                          /tmp/* /var/lib/apt/lists/*

RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com'))" >.Rprofile
RUN R -e "BiocManager::install('infercnv')"
RUN R -e "install.packages(c('cluster', 'Seurat', 'optparse', 'igraph', 'reshape2'), repos = 'http://cran.us.r-project.org')"

RUN git clone https://github.com/broadinstitute/infercnv && cd infercnv && \
      R CMD INSTALL . 

ENV PATH=${PATH}:/infercnv/scripts

CMD inferCNV.R --help
