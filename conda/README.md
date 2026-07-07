# Conda or Mamba Environment Setup

This directory provides environment setup for the public methods in ST-CNABench.

If `mamba` is available, prefer it because dependency resolution is usually faster.
In the commands below, set `PKG_MGR` once to either `mamba` or `conda`.
Keep `conda run` for running commands inside the created environments.

## Batch Install

```bash
bash conda/install_all_envs.sh
```

The batch installer already prefers `mamba` automatically when it is available.

## Individual Install

Set the package manager once before running the commands below:

```bash
export PKG_MGR=mamba
# or
# export PKG_MGR=conda

mkdir -p temp_install_src
```

### InferCNV

```bash
$PKG_MGR create -n infercnv_env -c conda-forge r-base=4.4 -y
$PKG_MGR install -n infercnv_env -c conda-forge jags zlib r-seurat r-base=4.4 -y
conda run -n infercnv_env Rscript -e "install.packages('data.table', repos='https://cran.rstudio.com')"
conda run -n infercnv_env Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cran.rstudio.com'); BiocManager::install('infercnv', ask=FALSE)"
```

### CopyKAT

```bash
$PKG_MGR create -n copykat_env -c conda-forge r-base=4.4 -y
$PKG_MGR install -n copykat_env -c conda-forge zlib r-seurat -y
conda run -n copykat_env Rscript -e "install.packages(c('data.table','parallelDist','dlm','gplots','RColorBrewer','mixtools','cluster','MCMCpack'), repos='https://cran.rstudio.com')"
git clone https://github.com/navinlabcode/copykat.git temp_install_src/copykat
conda run -n copykat_env R CMD INSTALL temp_install_src/copykat
```

### SCEVAN

```bash
$PKG_MGR create -n scevan_env -c conda-forge r-base=4.4 -y
$PKG_MGR install -n scevan_env -c conda-forge -c bioconda -y \
  bioconductor-ggtree git zlib cellsnp-lite \
  r-igraph r-gdtools libxml2 glpk gmp cairo harfbuzz fribidi \
  bioconductor-scran bioconductor-bluster bioconductor-summarizedexperiment \
  bioconductor-delayedarray bioconductor-singlecellexperiment \
  bioconductor-sparsearray bioconductor-scaledmatrix \
  bioconductor-beachmat bioconductor-biocsingular bioconductor-scuttle
conda run -n scevan_env Rscript -e "install.packages(c('doParallel','ggplot2','parallelDist','pheatmap','forcats','dplyr','cluster','Rtsne','ape','tidytree','ggrepel','BiocManager','R.utils'), repos='https://cran.rstudio.com')"
git clone https://github.com/miccec/yaGST.git temp_install_src/yaGST
conda run -n scevan_env R CMD INSTALL temp_install_src/yaGST
conda run -n scevan_env Rscript -e "BiocManager::install(c('ggiraph','fgsea'), ask=FALSE, update=FALSE)"
git clone https://github.com/AntonioDeFalco/SCEVAN.git temp_install_src/SCEVAN
conda run -n scevan_env R CMD INSTALL temp_install_src/SCEVAN
```

### Clonalscope

```bash
$PKG_MGR create -n clonalscope_env -c conda-forge r-base=4.4 -y
$PKG_MGR install -n clonalscope_env -c conda-forge zlib bzip2 xz libxml2 curl pkg-config -y
$PKG_MGR install -n clonalscope_env -c conda-forge -c bioconda -y \
  r-seurat bioconductor-aucell bioconductor-xvector \
  bioconductor-rtracklayer bioconductor-biocfilecache \
  r-pheatmap r-hiddenmarkov r-amap
git clone https://github.com/seasoncloud/Clonalscope.git temp_install_src/Clonalscope
conda run -n clonalscope_env R CMD INSTALL temp_install_src/Clonalscope
```

### Numbat

```bash
$PKG_MGR create -n numbat_env -c conda-forge r-base=4.4 -y
$PKG_MGR install -n numbat_env -c conda-forge -c bioconda -y \
  zlib harfbuzz fribidi freetype pkg-config xorg-xproto xorg-libx11 xorg-libxt libpng cairo \
  r-cairo r-textshaping r-ragg r-vcfr r-igraph r-dqrng r-rcpp r-rcppeigen \
  r-catools r-dendextend r-ggraph r-logger r-optparse r-pryr r-r.utils \
  bioconductor-ggtree bioconductor-genomicranges bioconductor-iranges \
  cellsnp-lite bioconductor-genomeinfodbdata
conda run -n numbat_env Rscript -e "if (!require('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cran.rstudio.com'); BiocManager::install(update=FALSE, ask=FALSE); install.packages(c('data.table','png','remotes','hahmmr','scistreer'), repos='https://cran.rstudio.com', update=FALSE)"
git clone https://github.com/kharchenkolab/numbat.git temp_install_src/numbat
conda run -n numbat_env R CMD INSTALL temp_install_src/numbat
```

### Xclone

```bash
$PKG_MGR create -n xclone_env python=3.9 -y
conda run -n xclone_env pip install xclone "xcltk>=0.5.2"
$PKG_MGR install -n xclone_env -c conda-forge -c bioconda bcftools cellsnp-lite -y
```

### CalicoST

```bash
git clone https://github.com/raphael-group/CalicoST.git temp_install_src/CalicoST
$PKG_MGR env create -f temp_install_src/CalicoST/environment.yml --name calicost_env
conda run -n calicost_env pip install temp_install_src/CalicoST
```

### STARCH

```bash
$PKG_MGR create -n starch_env -c conda-forge python=3.9 -y
$PKG_MGR install -n starch_env -c conda-forge -y numpy pandas scipy scikit-learn hmmlearn scanpy
```
