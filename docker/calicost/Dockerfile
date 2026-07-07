FROM continuumio/miniconda3:latest

LABEL org.opencontainers.image.title="ST-CNABench-CalicoST" \
      org.opencontainers.image.authors="HAN SHI"

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app
RUN git clone https://github.com/raphael-group/CalicoST.git /app

RUN conda install mamba -n base -c conda-forge -y && \
    mamba env create -f environment.yml --name calicost_env && \
    conda clean -afy

RUN conda run -n calicost_env pip install .

ENV PATH="/opt/conda/envs/calicost_env/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=calicost_env
