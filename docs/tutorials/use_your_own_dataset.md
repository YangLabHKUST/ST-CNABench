# Use Your Own Dataset

This tutorial outlines the minimum path for running ST-CNVBench on your own dataset.

## 1. Decide The Raw Input Layout

The public `prep` step currently supports:

- `SpaceRanger`
- `STpipeline`

See [Dataset Preparation](../data_preparation.md) for the exact input contract.

## 2. Create `data.yaml`

Start from:

```text
configs/templates/data.template.yaml
```

Fill in:

- dataset identity fields such as `dataset_id`, `platform`, `genome`, and `species`
- `raw.*` input paths
- `output.root`
- GT paths only for the evaluation tasks you actually plan to run

## 3. Prepare Runtime Requirements

Before running models, make sure the required runtime and references are available.

- base installation: [Installation](../installation.md)
- external tools: [External Tools And Runtime Notes](../external_tools.md)
- reference data: [Reference Data](../reference_data.md)

## 4. Run Step By Step

Start with `prep` only:

```bash
st-cnvbench --steps prep --data-config path/to/data.yaml
```

Then run selected methods:

```bash
st-cnvbench --steps run \
  --data-config path/to/data.yaml \
  --model-config path/to/models.yaml \
  --exec-mode conda \
  --models CopyKAT InferCNV
```

Then run only the evaluation tasks supported by your GT files:

```bash
st-cnvbench --steps eval \
  --data-config path/to/data.yaml \
  --eval-config path/to/eval.yaml \
  --models CopyKAT InferCNV \
  --eval-tasks cnv_profile
```

## 5. Keep The First Run Small

For a first pass, it is usually better to:

- start with one dataset
- enable one or two methods
- run one evaluation task
- verify outputs before scaling up

## Try Next

- For the packaged cSCC demo, go to [Quickstart Demo And Expected Outputs](quickstart_demo.md)
- For the CNV profile task example, go to [CNV Profile Task Example](cnv_profile_hcc2t.md)
- For the tumor-normal task example, go to [Tumor-Normal Classification Task Example](tumor_normal_hcc2t.md)
- For the subclone task example, go to [Subclone Identification Task Example](subclone_identification_slidednaseq.md)
- To configure or add a method, go to [Use Your Own Model](use_your_own_model.md)
