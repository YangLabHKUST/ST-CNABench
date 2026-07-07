FROM continuumio/miniconda3:latest

LABEL org.opencontainers.image.title="ST-CNABench-XClone"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libz-dev \
    git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN conda create -y -n xclone_env \
    -c conda-forge -c bioconda \
    python=3.9 bcftools cellsnp-lite && \
    conda clean -afy

RUN conda run -n xclone_env pip install --no-cache-dir \
    xclone \
    "xcltk>=0.5.2" \
    scanpy scipy numpy pandas matplotlib

ENV PATH="/opt/conda/envs/xclone_env/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=xclone_env
