# ST-CNABench

ST-CNABench is a public benchmark framework for copy number alteration inference on spatial transcriptomics.
It provides one controller for dataset preparation, model execution, and evaluation through a unified `prep -> run -> eval` workflow.

The current public release includes 9 CNA inference methods.

## Included Methods

- [`CalicoST`](https://github.com/raphael-group/CalicoST)
- [`CopyKAT`](https://github.com/navinlabcode/copykat)
- [`InferCNV`](https://github.com/broadinstitute/infercnv)
- [`Clonalscope`](https://github.com/seasoncloud/Clonalscope) (`Clonalscope_NoWGS`, `Clonalscope_WGS`)
- [`Numbat`](https://github.com/kharchenkolab/numbat)
- [`Xclone`](https://github.com/single-cell-genetics/XClone)
- [`SCEVAN`](https://github.com/AntonioDeFalco/SCEVAN)
- [`STARCH`](https://github.com/raphael-group/STARCH)
- [`SlideCNA`](https://github.com/dkzhang777/SlideCNA)

## Installation

```bash
git clone https://github.com/YangLabHKUST/ST-CNABench.git
cd ST-CNABench
conda create -n benchmark_env python=3.10 -y
conda activate benchmark_env
pip install -e .
st-cnabench --help
```

For method-specific runtime environments, external tools, and reference data, see the installation guide:
[https://cnadocs.readthedocs.io/en/latest/installation/](https://cnadocs.readthedocs.io/en/latest/installation/)

## Run

Prepare data:

```bash
st-cnabench --steps prep \
  --data-config <DATA_CONFIG> \
  --prep-ids <DATASET_ID>
```

Run methods:

```bash
st-cnabench --steps run \
  --data-config <DATA_CONFIG> \
  --model-config <MODEL_CONFIG> \
  --prep-ids <DATASET_ID> \
  --exec-mode <conda|docker|apptainer> \
  --models <METHOD_1> <METHOD_2>
```

Evaluate results:

```bash
st-cnabench --steps eval \
  --data-config <DATA_CONFIG> \
  --eval-config <EVAL_CONFIG> \
  --prep-ids <DATASET_ID> \
  --models <METHOD_1> <METHOD_2> \
  --eval-tasks <TASK_NAME>
```

## Tutorials

We provide demo data and tutorial workflows for:

- the packaged cSCC demo pipeline
- `cna_profile` task example
- `tumor_normal` task example
- subclone identification task example

Documentation site:
[https://cnadocs.readthedocs.io/en/latest/](https://cnadocs.readthedocs.io/en/latest/)

Quickstart tutorial:
[https://cnadocs.readthedocs.io/en/latest/tutorials/quickstart_demo/](https://cnadocs.readthedocs.io/en/latest/tutorials/quickstart_demo/)

## Contact

HAN Shi  
shanav@connect.ust.hk
