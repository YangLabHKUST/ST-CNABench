# Quickstart Demo And Expected Outputs

This page shows the full `prep -> run -> eval` command flow for the public cSCC demo and the expected output layout.
For the public release, the prepared demo bundle and example outputs should be distributed separately from git through a Google Drive download link.
After extracting that bundle, the expected layout should start under:

```text
demo_runs/cscc_demo/
```

The shipped demo uses the real `P6_vis_rep1` fixture and the example config set under `configs/examples/cscc_demo/`.
The demo dataset ID used in the commands below is `P6_vis_rep1`.

## Before You Start

1. Install the controller environment
2. Prepare method runtime environments
3. Download the external tools and reference data
4. Download the public demo bundle when that link is published
5. Use the provided demo configs

The default public demo enables `CopyKAT`, so it does not require allele-aware reference bundles or external tool source trees.

## Demo Bundle

Add your public Google Drive link for the packaged cSCC demo bundle here before release.
The unpacked bundle should provide the `demo_runs/cscc_demo/` layout referenced below.

## Included Example Outputs

The packaged example run can be used as a reference for the expected directory layout.
Key subdirectories include:

- `demo_runs/cscc_demo/processed_data/P6_vis_rep1/`
- `demo_runs/cscc_demo/raw_results/P6_vis_rep1/CopyKAT/`
- `demo_runs/cscc_demo/evaluation/P6_vis_rep1/`

## Step 1: Prepare Data

Run:

```bash
st-cnvbench --steps prep \
  --data-config configs/examples/cscc_demo/data.yaml \
  --prep-ids P6_vis_rep1
```

Check the prepared dataset under:

```text
demo_runs/cscc_demo/processed_data/P6_vis_rep1/
```

Expected outputs include:

- `filtered_feature_bc_matrix/`
- `filtered_feature_bc_matrix.h5ad`
- `spatial/tissue_positions.csv`
- `metadata_P6_vis_rep1_tumor_normal.tsv`

## Step 2: Run A Model

Example with `CopyKAT` in conda mode:

```bash
st-cnvbench --steps run \
  --data-config configs/examples/cscc_demo/data.yaml \
  --model-config configs/examples/cscc_demo/models.yaml \
  --prep-ids P6_vis_rep1 \
  --exec-mode conda \
  --models CopyKAT
```

Check raw model outputs under:

```text
demo_runs/cscc_demo/raw_results/P6_vis_rep1/CopyKAT/
```

Important files to inspect:

- `CopyKAT_run.log`
- `CopyKAT_run.perf`
- `P6_vis_rep1_copykat_prediction.txt`
- `P6_vis_rep1_copykat_CNA_results.txt`
- `copykat_subcluster_results.txt`

## Step 3: Evaluate

Run:

```bash
st-cnvbench --steps eval \
  --data-config configs/examples/cscc_demo/data.yaml \
  --eval-config configs/examples/cscc_demo/eval.yaml \
  --prep-ids P6_vis_rep1 \
  --models CopyKAT \
  --eval-tasks cnv_profile
```

Check evaluation outputs under:

```text
demo_runs/cscc_demo/evaluation/P6_vis_rep1/
```

For the default demo command above, the main task output is under:

```text
demo_runs/cscc_demo/evaluation/P6_vis_rep1/cnv_profile/
```

Typical outputs include:

- metrics summary tables
- CNV profile plots
- karyogram-level comparison plots
- per-method exported CNV profile tables

## Full Demo Command

If you want to run all three stages in one command, use:

```bash
st-cnvbench \
  --steps prep run eval \
  --data-config configs/examples/cscc_demo/data.yaml \
  --model-config configs/examples/cscc_demo/models.yaml \
  --eval-config configs/examples/cscc_demo/eval.yaml \
  --prep-ids P6_vis_rep1 \
  --exec-mode conda \
  --models CopyKAT \
  --eval-tasks cnv_profile
```

## Try Next

- For the CNV profile task example, go to [CNV Profile Task Example](cnv_profile_hcc2t.md)
- For the tumor-normal task example, go to [Tumor-Normal Classification Task Example](tumor_normal_hcc2t.md)
- For the subclone task example, go to [Subclone Identification Task Example](subclone_identification_slidednaseq.md)
- To use your own dataset, go to [Use Your Own Dataset](use_your_own_dataset.md)
- To configure or add a method, go to [Use Your Own Model](use_your_own_model.md)
