# Installation

ST-CNVBench separates the Python controller from the per-method runtime environments.

## Choose A Runtime Mode

| Runtime     | Best for                                         |
| ----------- | ------------------------------------------------ |
| `conda`     | Local use and method debugging                   |
| `docker`    | Reproducible runs on systems with Docker support |
| `apptainer` | HPC systems without Docker daemon access         |

### Conda Or Mamba

If `mamba` is available, prefer it for faster dependency solving.

```bash
bash conda/install_all_envs.sh
```

Detailed per-method environment setup is in `conda/README.md` in the main repository.

### Docker

Pull only the images for methods you enable.

```bash
docker pull hans0410/cnv-benchmark-copykat:1.1.0
```

Image names are listed in `configs/templates/models.template.yaml`. More runtime notes are in [External Tools And Runtime Notes](external_tools.md).

### Apptainer

Create `.sif` files from the same Docker images and set `apptainer_sif` in `models.yaml`.

```bash
mkdir -p apptainer_sif
apptainer pull apptainer_sif/copykat.sif docker://hans0410/cnv-benchmark-copykat:1.1.0
```

## Install The Controller

Clone the repository, then create the controller environment from the repository root.

```bash
git clone https://github.com/YangLabHKUST/ST-CNVBench.git
cd STCNV-Bench
conda create -n benchmark_env python=3.10 -y
conda activate benchmark_env
pip install -e .
st-cnvbench --help
```

The package requires Python `>=3.10`. The initial `pip install -e .` typically takes about `10-20 min`, depending on your machine and network.

## Install External Tools

Only a subset of wrappers need external source trees:

- `CalicoST`
- `Clonalscope_NoWGS`
- `Clonalscope_WGS`
- `Numbat`
- `Xclone`

From the repository root, run:

```bash
mkdir -p external_tools
cd external_tools

wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar -xzf Eagle_v2.4.1.tar.gz
rm Eagle_v2.4.1.tar.gz

git clone https://github.com/raphael-group/CalicoST.git CalicoST
git clone https://github.com/seasoncloud/Clonalscope.git clonalscope
git clone https://github.com/kharchenkolab/numbat.git numbat
```

If you keep the default public layout under `external_tools/`, the public config templates already point to the expected locations.
Detailed path mapping and model-specific notes are in [External Tools And Runtime Notes](external_tools.md).

## Install Reference Data

Small hg38 annotation files are already bundled in git under `refs/hg38_genome_info/`.

Large population phasing references are required only for allele-aware wrappers:

- `CalicoST`
- `Numbat`
- `Xclone`

Download the bundle from:

- [Population phasing reference bundle (Google Drive)](https://drive.google.com/file/d/12-hEUoDdaXTdap-Ro4Ekx4XCSlbYGP_X/view?usp=drive_link)

After download, extract it under:

```text
refs/
└── population_phasing/
```

Detailed file layout and usage notes are in [Reference Data](reference_data.md).

## Demo Bundle

Download the public cSCC demo bundle from:

- [cSCC demo bundle (Google Drive)](https://drive.google.com/file/d/1WvsmmnYYDYwAq87ZdIw50NHOyNNwXhEp/view?usp=drive_link)

After extraction, the expected example outputs should appear under `demo_runs/cscc_demo/`.
See [Quickstart Demo And Expected Outputs](tutorials/quickstart_demo.md) for the expected layout and demo commands.

## After Installation

- To run the public demo, go to [Quickstart Demo And Expected Outputs](tutorials/quickstart_demo.md)
- To configure your own benchmark run, go to [Dataset Preparation](data_preparation.md), [Model Run](model_run.md), and [Evaluation](evaluation.md)
