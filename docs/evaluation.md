# Evaluation Guide

This page describes the public `eval.yaml` contract, available evaluation tasks, and the ground-truth files required by each task.

Evaluation assumes that model outputs already exist under `project_settings.results_dir`. See [model_run.md](model_run.md) for model-running configuration.

## Run Command

Run one evaluation task for selected models:

```bash
st-cnabench --steps eval \
  --data-config configs/examples/cscc_demo/data.yaml \
  --eval-config configs/examples/cscc_demo/eval.yaml \
  --models CopyKAT \
  --eval-tasks cna_profile
```

Run evaluation for selected datasets, models and tasks:

```bash
st-cnabench --steps eval \
  --data-config data.yaml \
  --eval-config eval.yaml \
  --prep-ids sample_1 sample_2 \
  --models InferCNV CopyKAT STARCH \
  --eval-tasks cna_profile tumor_normal
```

If `--eval-tasks` is omitted, the controller attempts all registered tasks. For public examples, specify only tasks whose required GT files are present.

## Global Settings

`project_settings` defines where raw model outputs are read and where evaluation outputs are written.

| Field         | Meaning                                          |
| ------------- | ------------------------------------------------ |
| `root_dir`    | Project root used to resolve relative paths.     |
| `results_dir` | Raw model result root produced by `--steps run`. |
| `eval_dir`    | Evaluation output root.                          |

`global_params` provides shared genomic settings.

| Field             | Meaning                                                             |
| ----------------- | ------------------------------------------------------------------- |
| `bin_size`        | Genomic bin size used by CNA profile and clone-level mapping tasks. |
| `genome_version`  | Reference genome label, currently used by clonal mapping utilities. |
| `gene_annot_path` | Gene annotation table used by expression-to-genomic-bin conversion. |

Example:

```yaml
project_settings:
  root_dir: "."
  results_dir: "${root_dir}/outputs/raw_results"
  eval_dir: "${root_dir}/outputs/evaluation"

global_params:
  bin_size: 100000
  genome_version: "hg38"
  gene_annot_path: "${root_dir}/refs/hg38_genome_info/hg38_genes_annot.txt"
```

## Model Loader Settings

`eval_list` maps each model result directory to one or more evaluation loaders.

```yaml
eval_list:
  InferCNV:
    eval_name: ["InferCNV_expr", "InferCNV_cna"]
  CopyKAT:
    eval_name: ["CopyKAT"]
```

The key, for example `InferCNV`, must match the model result directory under:

```text
<results_dir>/<dataset_id>/<model_name>/
```

The `eval_name` values select loader adapters for the output format produced by that model.

| Model key           | Loader names                    |
| ------------------- | ------------------------------- |
| `InferCNV`          | `InferCNV_expr`, `InferCNV_cna` |
| `CopyKAT`           | `CopyKAT`                       |
| `SCEVAN`            | `SCEVAN_expr`, `SCEVAN_cna`     |
| `Clonalscope_WGS`   | `Clonalscope_WGS`               |
| `Clonalscope_NoWGS` | `Clonalscope_NoWGS`             |
| `Numbat`            | `Numbat_expr`, `Numbat_cna`     |
| `CalicoST`          | `CalicoST`                      |
| `STARCH`            | `STARCH`                        |
| `Xclone`            | `Xclone_expr`, `Xclone_cna`     |

Use the loader that matches the model output type you want to evaluate. Expression-derived loaders usually need `gene_annot_path` to map genes to genomic bins.

## Evaluation Tasks

Available task names for `--eval-tasks`:

| Task                          | Purpose                                                                   | Required GT or inputs                                                                           |
| ----------------------------- | ------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `efficiency`                  | Runtime and memory summary.                                               | No biological GT. Requires conda-mode `.perf` files generated during model execution.           |
| `cna_resolution`              | CNA resolution comparison across model outputs.                           | No GT. Requires model CNA outputs and `global_params.gene_annot_path`.                          |
| `tumor_normal`                | Tumor/normal prediction evaluation.                                       | `raw.tumor_normal_gt`; subset mode also requires `raw.tumor_normal`.                            |
| `cna_profile`                 | CNA profile concordance against sample-level CNA GT.                      | `raw.cna_gt` as a FACETS/VCF-like segment file.                                                 |
| `subclone_detection_in_slice` | Spot-level subclone assignment and clone-profile matching within a slice. | `raw.subclone_gt`; clone-profile metrics also need clone-level CNA profiles from `raw.cna_gt`.  |
| `subclone_detection_organ`    | Organ-level subclone assignment across slices or merged samples.          | `raw.subclone_gt`.                                                                              |
| `clonal_evolution`            | Clone-level CNA tree and spatial phylogeography plots.                    | Model clone labels/profiles and spatial coordinates; no separate GT tree is currently required. |

## GT Files In `data.yaml`

Evaluation GT paths are defined in the dataset config, not in `eval.yaml`.

### `raw.tumor_normal_gt`

Used only by `tumor_normal` evaluation.

Expected format: tab-delimited spot labels with two or four columns. The first two columns are interpreted as barcode and label.

```text
Barcode  tumor_normal
AAAC...  tumor
AAAG...  normal
```

Accepted labels are `tumor` and `normal`.

### `raw.tumor_normal`

This is not GT. It is the reference-normal annotation used during model running and subset-mode evaluation.

In `tumor_normal_mode: subset`, evaluation removes the reference normal spots from the comparison set, so both `raw.tumor_normal` and `raw.tumor_normal_gt` are required.

### `raw.cna_gt` For `cna_profile`

For `cna_profile`, `raw.cna_gt` should point to a FACETS/VCF-like segment file.

Required content:

| Field    | Meaning                                                            |
| -------- | ------------------------------------------------------------------ |
| `#CHROM` | Chromosome.                                                        |
| `POS`    | Segment start.                                                     |
| `INFO`   | Must include `END`; `TCN_EM` and `LCN_EM` are used when available. |

The evaluator derives:

| Derived field | Source                                                            |
| ------------- | ----------------------------------------------------------------- |
| `GT_Score`    | `log2(TCN_EM / 2)` with a small pseudocount for zero-copy events. |
| `GT_Event`    | loss, neutral, or gain from total copy number.                    |
| `LOH_Status`  | CN-LOH from `TCN_EM == 2` and `LCN_EM == 0`.                      |

### `raw.subclone_gt`

Used by subclone detection tasks.

Expected format: tab-delimited spot-to-subclone labels with two columns.

```text
Barcode  subclone
AAAC...  clone_1
AAAG...  clone_2
```

The evaluator renames these columns internally to `Barcodes` and `Label_preds`.

### Clone-Level CNA GT For Subclone Tasks

`subclone_detection_in_slice` can also compare predicted clone CNA profiles against GT clone CNA profiles.

Current loader expectation: a directory containing files named:

```text
<clone_id>_GT_Profile.txt
```

Each profile file should contain genomic bins and a CNA score column such as `CN_Score` or `CN_Score_Continuous`.

This is separate from the FACETS/VCF-like file used by `cna_profile`. If a dataset needs both sample-level `cna_profile` GT and clone-level subclone CNA GT, keep the paths clear in the dataset config used for that run.

### `raw.beads_mapping`

Only needed for Slide-DNA-seq style subclone evaluation when model predictions are made on pseudo-barcodes and must be mapped back to original bead barcodes.

Expected columns:

```text
pseudo_barcode,original_barcode
```

## Spatial Inputs

Evaluation uses standardized spatial files from `prep`:

```text
<output.root>/spatial/tissue_positions.csv
<output.root>/spatial/scalefactors_json.json
<output.root>/spatial/tissue_hires_image.png  # optional
```

Spatial coherence mode is selected from `platform`:

| `platform`   | Spatial mode                                   |
| ------------ | ---------------------------------------------- |
| `ST`         | KNN-based spatial coherence.                   |
| other values | Visium-style distance-based spatial coherence. |

## Output Layout

Evaluation outputs are written under:

```text
<eval_dir>/<dataset_id>/
```

Common task directories include:

```text
computational_efficiency/
cna_resolution/
tumor_normal/
cna_profile/
subclone_detection/
GT/
```

Task outputs include formatted intermediate tables, metrics summaries, and plots. Missing required GT or missing model result directories are reported explicitly rather than silently substituted.
