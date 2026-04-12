# Dataset Preparation

This page describes the public `data.yaml` contract and the standardized outputs produced by `--steps prep`.

## Run Command

Prepare one or more datasets:

```bash
st-cnvbench --steps prep \
  --data-config data.yaml \
  --prep-ids sample_1 sample_2
```

For the packaged cSCC walkthrough, see [Quickstart Demo And Expected Outputs](tutorials/quickstart_demo.md).

## Dataset Entries

Each dataset entry in `data.yaml` represents one sample or sample group.

Required top-level fields:

| Field | Meaning |
| --- | --- |
| `dataset_id` | Unique dataset name used in output and result paths. |
| `platform` | Spatial platform name, for example `Visium`, `ST`, or `SlideDNAseq`. |
| `format` | Raw input layout currently supported by `prep`. |
| `genome` | Genome label used by model wrappers, for example `hg38`. |
| `species` | Species label, for example `human` or `mouse`. |
| `ref_norm` | Whether model running should use reference normal spots when supported. |
| `tumor_normal_mode` | How tumor-normal annotations are used in the run/eval flow. |
| `raw` | Raw input paths and GT paths. |
| `output.root` | Standardized output directory for this dataset. |

Allowed values:

- `format`: `SpaceRanger` or `STpipeline`
- `tumor_normal_mode`: `subset`, `de_novo`, or `off`

## `raw.*` Fields

| Field | Meaning |
| --- | --- |
| `root` | Base raw-data directory. Used for path interpolation and optional auto-discovery. |
| `counts` | Optional explicit count matrix path when the input is not discovered from `raw.root`. |
| `barcodes` | Optional explicit barcode file path. |
| `features` | Optional explicit feature file path. |
| `scalefactors` | Spatial scale factors. Required for `STpipeline`; optional for `SpaceRanger` if present under `raw.root/spatial/`. |
| `tissue_positions` | Spatial coordinate table when not discovered automatically. |
| `tissue_image` | Optional tissue image used by prep and plotting. |
| `tumor_normal` | Model-run annotation for reference-normal selection. In `de_novo` mode this can be `null`. |
| `tumor_normal_gt` | Ground-truth tumor-normal annotation used only by evaluation. |
| `subclone_gt` | Ground-truth spot-level subclone labels for subclone tasks. |
| `cnv_gt` | Ground-truth CNV profile input for `cnv_profile`, and clone-level CNV GT for subclone tasks when applicable. |
| `bam`, `bai` | Alignment files needed by allele-aware wrappers such as `Numbat` and `Xclone`. |
| `wgs_wes_tumor_bedg`, `wgs_wes_normal_bedg` | WGS/WES count tracks for `Clonalscope_WGS`. |
| `beads_mapping` | Slide-DNA-seq style pseudo-barcode to original-barcode mapping for subclone evaluation. |

## Tumor-Normal Modes

`tumor_normal_mode` controls how the pipeline treats tumor-normal annotations.

| Mode | Meaning |
| --- | --- |
| `subset` | Use `raw.tumor_normal` to provide reference-normal spots during model run, and remove those reference spots from tumor-normal evaluation. |
| `de_novo` | Run without reference-normal labels. Models that support de novo operation can still be benchmarked. |
| `off` | Disable tumor-normal evaluation for this dataset. |

Conceptually:

- `raw.tumor_normal` is a model-run input
- `raw.tumor_normal_gt` is an evaluation GT file

## Minimal Example

```yaml
datasets:
  sample_1:
    dataset_id: sample_1
    platform: Visium
    format: SpaceRanger
    genome: hg38
    species: human
    ref_norm: false
    tumor_normal_mode: de_novo
    raw:
      root: /path/to/raw/sample_1
      tissue_positions: null
      scalefactors: null
      tissue_image: null
      tumor_normal: null
      tumor_normal_gt: null
      subclone_gt: null
      cnv_gt: null
      bam: null
      bai: null
      wgs_wes_tumor_bedg: null
      wgs_wes_normal_bedg: null
      beads_mapping: null
    output:
      root: /path/to/processed/sample_1
```

## Standardized Outputs Produced By `prep`

The public `prep` step creates a benchmark-ready dataset bundle under:

```text
<output.root>/
```

Common outputs include:

```text
filtered_feature_bc_matrix/
filtered_feature_bc_matrix.h5ad
spatial/
metadata_<dataset_id>_tumor_normal.tsv   # when a model-run tumor-normal annotation is available
```

Expected `spatial/` contents:

```text
spatial/tissue_positions.csv
spatial/scalefactors_json.json
spatial/tissue_hires_image.png   # optional when an image is available
```

The `.h5ad` file is assembled during preparation and is part of the standardized output bundle used by downstream benchmarking steps.

## Notes

- Missing required inputs should fail loudly rather than being silently substituted.
- `de_novo` datasets can still run, but only with wrappers that support reference-free operation.
- Keep GT paths `null` for evaluation tasks you do not plan to run.
- Use one dataset entry per logical sample or sample group you want the benchmark controller to manage.
