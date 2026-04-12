# ST-CNVBench

ST-CNVBench is a public benchmark framework for copy number variation inference on spatial transcriptomics.
It provides one controller for dataset preparation, model execution, and evaluation through a unified `prep -> run -> eval` workflow.

The current public release includes 8 CNV inference methods.

## Included Methods

- [`CalicoST`](https://github.com/raphael-group/CalicoST)
- [`CopyKAT`](https://github.com/navinlabcode/copykat)
- [`InferCNV`](https://github.com/broadinstitute/infercnv)
- [`Clonalscope`](https://github.com/seasoncloud/Clonalscope) (`Clonalscope_NoWGS`, `Clonalscope_WGS`)
- [`Numbat`](https://github.com/kharchenkolab/numbat)
- [`Xclone`](https://github.com/single-cell-genetics/XClone)
- [`SCEVAN`](https://github.com/AntonioDeFalco/SCEVAN)
- [`STARCH`](https://github.com/raphael-group/STARCH)

## Installation

```bash
git clone https://github.com/Hans-0410/STCNV-Bench.git
cd STCNV-Bench
conda create -n benchmark_env python=3.10 -y
conda activate benchmark_env
pip install -e .
st-cnvbench --help
```

For method-specific runtime environments, external tools, and reference data, see [docs/installation.md](docs/installation.md).

## Run

Prepare data:

```bash
st-cnvbench --steps prep \
  --data-config <DATA_CONFIG> \
  --prep-ids <DATASET_ID>
```

Run methods:

```bash
st-cnvbench --steps run \
  --data-config <DATA_CONFIG> \
  --model-config <MODEL_CONFIG> \
  --prep-ids <DATASET_ID> \
  --exec-mode <conda|docker|apptainer> \
  --models <METHOD_1> <METHOD_2>
```

Evaluate results:

```bash
st-cnvbench --steps eval \
  --data-config <DATA_CONFIG> \
  --eval-config <EVAL_CONFIG> \
  --prep-ids <DATASET_ID> \
  --models <METHOD_1> <METHOD_2> \
  --eval-tasks <TASK_NAME>
```

## Tutorials

We provide demo data and tutorial workflows for:

- the packaged cSCC demo pipeline
- `cnv_profile` task example
- `tumor_normal` task example
- subclone identification task example

Documentation site:
[https://cnvdocs.readthedocs.io/en/latest/](https://cnvdocs.readthedocs.io/en/latest/)

Quickstart tutorial:
[https://cnvdocs.readthedocs.io/en/latest/tutorials/quickstart_demo/](https://cnvdocs.readthedocs.io/en/latest/tutorials/quickstart_demo/)

## Contact

HAN Shi  
shanav@connect.ust.hk
