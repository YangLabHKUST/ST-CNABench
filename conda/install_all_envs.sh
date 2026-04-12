#!/usr/bin/env bash
set -euo pipefail

if command -v mamba >/dev/null 2>&1; then
    BIN="mamba"
else
    BIN="conda"
fi

echo "Package manager: $BIN"

check_env() {
    conda env list | awk '{print $1}' | grep -qx "$1"
}

TEMP_DIR="temp_install_src"

cleanup() {
    if [ -d "$TEMP_DIR" ]; then
        echo ">>> [Cleanup] Removing temporary source directory..."
        rm -rf "$TEMP_DIR"
    fi
}

trap cleanup EXIT
mkdir -p "$TEMP_DIR"

ENV_NAME="calicost_env"
if check_env "$ENV_NAME"; then
    echo ">>> [$ENV_NAME] Environment already exists."
else
    echo ">>> [$ENV_NAME] Installing..."
    git clone https://github.com/raphael-group/CalicoST.git "$TEMP_DIR/CalicoST"
    conda env create -f "$TEMP_DIR/CalicoST/environment.yml" --name "$ENV_NAME"
    conda run -n "$ENV_NAME" pip install "$TEMP_DIR/CalicoST"
fi

ENV_NAME="xclone_env"
if check_env "$ENV_NAME"; then
    echo ">>> [$ENV_NAME] Environment already exists."
else
    echo ">>> [$ENV_NAME] Installing..."
    $BIN create -n "$ENV_NAME" python=3.9 -y
    conda run -n "$ENV_NAME" pip install xclone "xcltk>=0.5.2"
    $BIN install -n "$ENV_NAME" -c conda-forge -c bioconda bcftools cellsnp-lite -y
fi

ENV_NAME="infercnv_env"
if check_env "$ENV_NAME"; then
    echo ">>> [$ENV_NAME] Environment already exists."
else
    echo ">>> [$ENV_NAME] Installing..."
    $BIN create -n "$ENV_NAME" -c conda-forge r-base=4.4 -y
    $BIN install -n "$ENV_NAME" -c conda-forge jags zlib r-seurat r-base=4.4 -y
    conda run -n "$ENV_NAME" Rscript -e "install.packages('data.table', repos='https://cran.rstudio.com')"
    conda run -n "$ENV_NAME" Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cran.rstudio.com'); BiocManager::install('infercnv', ask=FALSE)"
fi

ENV_NAME="copykat_env"
if check_env "$ENV_NAME"; then
    echo ">>> [$ENV_NAME] Environment already exists."
else
    echo ">>> [$ENV_NAME] Installing..."
    $BIN create -n "$ENV_NAME" -c conda-forge r-base=4.4 -y
    $BIN install -n "$ENV_NAME" -c conda-forge zlib r-seurat -y
    conda run -n "$ENV_NAME" Rscript -e "install.packages(c('data.table','parallelDist','dlm','gplots','RColorBrewer','mixtools','cluster','MCMCpack'), repos='https://cran.rstudio.com')"
    git clone https://github.com/navinlabcode/copykat.git "$TEMP_DIR/copykat"
    conda run -n "$ENV_NAME" R CMD INSTALL "$TEMP_DIR/copykat"
fi

ENV_NAME="scevan_env"
if check_env "$ENV_NAME"; then
    echo ">>> [$ENV_NAME] Environment already exists."
else
    echo ">>> [$ENV_NAME] Installing..."
    $BIN create -n "$ENV_NAME" -c conda-forge r-base=4.4 -y
    $BIN install -n "$ENV_NAME" -c conda-forge -c bioconda -y \
        bioconductor-ggtree git zlib cellsnp-lite \
        r-igraph r-gdtools libxml2 glpk gmp cairo harfbuzz fribidi \
        bioconductor-scran bioconductor-bluster bioconductor-summarizedexperiment \
        bioconductor-delayedarray bioconductor-singlecellexperiment \
        bioconductor-sparsearray bioconductor-scaledmatrix \
        bioconductor-beachmat bioconductor-biocsingular bioconductor-scuttle
    conda run -n "$ENV_NAME" Rscript -e "install.packages(c('doParallel','ggplot2','parallelDist','pheatmap','forcats','dplyr','cluster','Rtsne','ape','tidytree','ggrepel','BiocManager','R.utils'), repos='https://cran.rstudio.com')"
    git clone https://github.com/miccec/yaGST.git "$TEMP_DIR/yaGST"
    conda run -n "$ENV_NAME" R CMD INSTALL "$TEMP_DIR/yaGST"
    conda run -n "$ENV_NAME" Rscript -e "BiocManager::install(c('ggiraph','fgsea'), ask=FALSE, update=FALSE)"
    git clone https://github.com/AntonioDeFalco/SCEVAN.git "$TEMP_DIR/SCEVAN"
    conda run -n "$ENV_NAME" R CMD INSTALL "$TEMP_DIR/SCEVAN"
fi

ENV_NAME="clonalscope_env"
if check_env "$ENV_NAME"; then
    echo ">>> [$ENV_NAME] Environment already exists."
else
    echo ">>> [$ENV_NAME] Installing..."
    $BIN create -n "$ENV_NAME" -c conda-forge r-base=4.4 -y
    $BIN install -n "$ENV_NAME" -c conda-forge zlib bzip2 xz libxml2 curl pkg-config -y
    $BIN install -n "$ENV_NAME" -c conda-forge -c bioconda -y \
        r-seurat bioconductor-aucell bioconductor-xvector \
        bioconductor-rtracklayer bioconductor-biocfilecache \
        r-pheatmap r-hiddenmarkov r-amap
    git clone https://github.com/seasoncloud/Clonalscope.git "$TEMP_DIR/Clonalscope"
    conda run -n "$ENV_NAME" R CMD INSTALL "$TEMP_DIR/Clonalscope"
fi

ENV_NAME="numbat_env"
if check_env "$ENV_NAME"; then
    echo ">>> [$ENV_NAME] Environment already exists."
else
    echo ">>> [$ENV_NAME] Installing..."
    $BIN create -n "$ENV_NAME" -c conda-forge r-base=4.4 -y
    $BIN install -n "$ENV_NAME" -c conda-forge -c bioconda -y \
        zlib harfbuzz fribidi freetype pkg-config xorg-xproto xorg-libx11 xorg-libxt libpng cairo \
        r-cairo r-textshaping r-ragg r-vcfr r-igraph r-dqrng r-rcpp r-rcppeigen \
        r-catools r-dendextend r-ggraph r-logger r-optparse r-pryr r-r.utils \
        bioconductor-ggtree bioconductor-genomicranges bioconductor-iranges \
        cellsnp-lite bioconductor-genomeinfodbdata
    conda run -n "$ENV_NAME" Rscript -e "if (!require('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cran.rstudio.com'); BiocManager::install(update=FALSE, ask=FALSE); install.packages(c('data.table','png','remotes','hahmmr','scistreer'), repos='https://cran.rstudio.com', update=FALSE)"
    git clone https://github.com/kharchenkolab/numbat.git "$TEMP_DIR/numbat"
    conda run -n "$ENV_NAME" R CMD INSTALL "$TEMP_DIR/numbat"
fi

ENV_NAME="starch_env"
if check_env "$ENV_NAME"; then
    echo ">>> [$ENV_NAME] Environment already exists."
else
    echo ">>> [$ENV_NAME] Installing..."
    $BIN create -n "$ENV_NAME" -c conda-forge python=3.9 -y
    $BIN install -n "$ENV_NAME" -c conda-forge -y numpy pandas scipy scikit-learn hmmlearn scanpy
fi

echo "===================================================="
echo ">>> All public conda environments checked and installed."
echo "===================================================="
