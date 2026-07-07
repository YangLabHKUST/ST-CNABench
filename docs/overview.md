# Overview

ST-CNABench provides a public benchmark controller for CNA inference on spatial transcriptomics datasets.
The workflow is organized into three main stages:

```text
prep -> run -> eval
```

## Design Goals

- standardize heterogeneous spatial transcriptomics inputs into one benchmark-ready layout
- run multiple CNA inference methods through one shared controller
- evaluate results with task-specific loaders and metrics
- keep the public surface config-driven and method-extensible

## Pipeline Stages

### Dataset Preparation

The `prep` stage standardizes raw spatial transcriptomics inputs into a common output bundle used by downstream wrappers.
Supported public input layouts are documented in [Dataset Preparation](data_preparation.md).

### Model Run

The `run` stage executes one or more CNA inference methods against prepared datasets using `conda`, `docker`, or `apptainer` runtime modes.
See [Model Run](model_run.md) for the public model config structure.

### Evaluation

The `eval` stage loads method outputs through public adapters and runs selected benchmark tasks such as `cna_profile` or `tumor_normal`.
See [Evaluation](evaluation.md) for available tasks and GT requirements.

## Public Examples

The docs include three task-focused examples:

- [CNA Profile Task Example](tutorials/cna_profile_hcc2t.md)
- [Tumor-Normal Classification Task Example](tutorials/tumor_normal_hcc2t.md)
- [Subclone Identification Task Example](tutorials/subclone_identification_slidednaseq.md)

For the packaged cSCC walkthrough, go to [Quickstart Demo And Expected Outputs](tutorials/quickstart_demo.md).
