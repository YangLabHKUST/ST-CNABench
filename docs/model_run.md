# Model Run Guide

This page describes the public `models.yaml` contract used by `st-cnabench --steps run`.

Model running assumes that `prep` has already produced a standardized dataset bundle. See [data_preparation.md](data_preparation.md) for the dataset input contract.

## Included Methods

The current public release includes 8 CNA inference methods:

- [`CalicoST`](https://github.com/raphael-group/CalicoST)
- [`CopyKAT`](https://github.com/navinlabcode/copykat)
- [`InferCNV`](https://github.com/broadinstitute/infercnv)
- [`Clonalscope`](https://github.com/seasoncloud/Clonalscope) (`Clonalscope_NoWGS`, `Clonalscope_WGS`)
- [`Numbat`](https://github.com/kharchenkolab/numbat)
- [`Xclone`](https://github.com/single-cell-genetics/XClone)
- [`SCEVAN`](https://github.com/AntonioDeFalco/SCEVAN)
- [`STARCH`](https://github.com/raphael-group/STARCH)

In the public config surface, `Clonalscope` is exposed as two wrappers so users can run the no-WGS and WGS-assisted modes separately.

## Run Command

Run all enabled models on all datasets:

```bash
st-cnabench --steps run \
  --data-config configs/examples/cscc_demo/data.yaml \
  --model-config configs/examples/cscc_demo/models.yaml
```

Run selected models or datasets:

```bash
st-cnabench --steps run \
  --data-config configs/examples/cscc_demo/data.yaml \
  --model-config configs/templates/models.template.yaml \
  --models CopyKAT InferCNV STARCH \
  --prep-ids P6_vis_rep1
```

Override the execution backend from the command line:

```bash
st-cnabench --steps run --exec-mode conda --data-config data.yaml --model-config models.yaml
```

## Global Settings

`project_settings` defines shared paths and the default execution backend.

| Field               | Meaning                                             |
| ------------------- | --------------------------------------------------- |
| `root_dir`          | Project root used to resolve relative paths.        |
| `results_dir`       | Output root for raw model results.                  |
| `refs`              | Reference data directory.                           |
| `exts`              | External tool source directory.                     |
| `default_exec_mode` | Default backend: `conda`, `docker`, or `apptainer`. |

Example:

```yaml
project_settings:
  root_dir: "."
  results_dir: "${root_dir}/outputs/raw_results"
  refs: "${root_dir}/refs"
  exts: "${root_dir}/external_tools"
  default_exec_mode: "conda"
```

## Shared Model Fields

Each model section is keyed by the public model name in the wrapper registry.

| Field           | Meaning                                                                        |
| --------------- | ------------------------------------------------------------------------------ |
| `enabled`       | Whether the wrapper is run by default.                                         |
| `model_name`    | Wrapper name expected by ST-CNABench. Keep this aligned with the section name. |
| `env_name`      | Conda environment used when `exec_mode` is `conda`.                            |
| `docker_image`  | Docker image used when `exec_mode` is `docker`.                                |
| `apptainer_sif` | SIF image path used when `exec_mode` is `apptainer`.                           |
| `per_dataset`   | Optional dataset-specific overrides.                                           |

Disable a model globally:

```yaml
CopyKAT:
  enabled: false
```

Override selected parameters for one dataset:

```yaml
STARCH:
  enabled: true
  n_clusters: 3
  per_dataset:
    P6_vis_rep1:
      n_clusters: 2
```

Disable one model for one dataset:

```yaml
InferCNV:
  enabled: true
  per_dataset:
    P6_vis_rep1:
      enabled: false
```

## Reference-Normal Mode

Reference-normal behavior is controlled by `ref_norm` in `data.yaml`, not by a separate model-level flag.

Use reference spots:

```yaml
ref_norm: true
tumor_normal_mode: subset
raw:
  tumor_normal: /path/to/metadata_tumor_normal.tsv
```

Run without provided reference spots:

```yaml
ref_norm: false
tumor_normal_mode: de_novo
raw:
  tumor_normal: null
```

Current wrapper behavior when `raw.tumor_normal` is not provided:

| Status        | Models                                                                                                  |
| ------------- | ------------------------------------------------------------------------------------------------------- |
| Supported     | `CalicoST`, `CopyKAT`, `InferCNV`, `SCEVAN`, `Numbat`, `STARCH`, `SlideCNA`, `Clonalscope_NoWGS`, `Clonalscope_WGS` |
| Not supported | `Xclone`                                                                                                |

## Key Parameters By Model

| Model               | Key fields                                                                                                                                |
| ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `CalicoST`          | `calicost_dir`, `eagle_dir`, `region_vcf`, `phasing_panel`, `use_tumor_purity`, `n_clones`, `n_threads`, `UMItag`, `cellTAG`              |
| `CopyKAT`           | `genome`, `win_size`, `ks_cut`, `distance`, `n_cores`, `n_clones`                                                                         |
| `InferCNV`          | `gene_order_file`, `cutoff`, `n_threads`, `k_obs_groups`                                                                                  |
| `Clonalscope_NoWGS` | `gene_coords_file`, `aux_data_dir`, `mincell`                                                                                             |
| `Clonalscope_WGS`   | `gene_coords_file`, `aux_data_dir`, `hmm_states`, `mincell`                                                                               |
| `Numbat`            | `pileup_script`, `eagle_path`, `genetic_map`, `snp_vcf_path`, `panel_dir`, `genome_version`, `n_threads`, `UMItag`, `cellTAG`, `n_clones` |
| `Xclone`            | `snp_vcf`, `gene_region`, `eagle_path`, `genetic_map`, `panel_dir`, `n_threads`, `UMItag`, `cellTAG`, `minCOUNT`, `minMAF`, `n_clusters`  |
| `SCEVAN`            | `n_threads`                                                                                                                               |
| `STARCH`            | `gene_mapping_file`, `n_clusters`, `beta_spot`, `platform`, `returnnormal`                                                                |
| `SlideCNA`          | `allowed_platforms`, `gene_annot_file`, `chrom_ord`, `spatial`, `roll_mean_window`, `avg_bead_per_bin`, `pos`, `pos_k`, `ex_k`, `max_k_silhouette`, `use_GO_terms` |

## Related Setup

- [Installation](installation.md)
- [External Tools And Runtime Notes](external_tools.md)
- [Reference Data](reference_data.md)

## Output Layout

Each enabled wrapper writes to:

```text
<results_dir>/<dataset_id>/<model_name>/
```

The wrapper also writes `<model_name>_run.log` in the same directory. If execution fails, ST-CNABench raises an error and reports the log path plus the last log lines.
