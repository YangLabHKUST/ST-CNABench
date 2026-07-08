# High-Resolution ST Exploration

This tutorial outlines a practical first-pass workflow for applying ST-CNABench to high-resolution spatial transcriptomics datasets.
High-resolution platforms can contain many small spots, beads, or cell bins, so the first decision is whether each unit has enough gene coverage for stable CNA inference.
For broader pre-run checks, review [Dataset Properties Quick Check](dataset_properties_quick_check.md) first.

## 1. Prepare The Dataset

Start by preparing a dataset entry in `data.yaml` following the public input contract in [Dataset Preparation](../data_preparation.md).

For high-resolution data, make sure each analysis unit has:

- a count matrix
- a barcode or bin identifier
- spatial coordinates
- genome and species metadata
- tumor-normal annotation if reference-normal mode is used

Run only the preparation step first:

```bash
st-cnabench --steps prep \
  --data-config data.yaml \
  --prep-ids <DATASET_ID>
```

Check that the standardized output bundle is created under:

```text
<output.root>/
```

Expected files include:

- `filtered_feature_bc_matrix/`
- `filtered_feature_bc_matrix.h5ad`
- `spatial/tissue_positions.csv`

## 2. Check Per-Unit Gene Coverage

Before running CNA inference, check whether the median detected genes per cell, spot, bead, or bin is around 1,000.
This is a practical QC checkpoint for high-resolution data because very sparse units can produce unstable CNA profiles.

Example QC check:

```python
import scanpy as sc
import numpy as np

adata = sc.read_h5ad("<output.root>/filtered_feature_bc_matrix.h5ad")
detected_genes = np.asarray((adata.X > 0).sum(axis=1)).ravel()
median_genes = float(np.median(detected_genes))
print(f"Median detected genes per unit: {median_genes:.1f}")
```

Interpretation:

- If the median is near or above ~1,000 genes per unit, proceed to model running.
- If the median is clearly below ~1,000 genes per unit, consider aggregating neighboring units before CNA inference.
- After aggregation, rerun the same QC check on the aggregated input.

Example QC figure:

- [High-resolution platform gene coverage QC violin plot](../assets/high_resolution_st_exploration/highres_platform_resolution_gene_qc_violin.pdf)

## 3. Aggregate Sparse High-Resolution Units If Needed

If the per-unit gene coverage is too sparse, aggregate nearby high-resolution units before CNA inference.
One practical strategy is spatial grid aggregation, where cell-level or bead-level units are grouped into larger grid bins.

For example, a nested grid can aggregate cell-level CosMx units into larger 16um or 32um spatial bins:

- [Nested 16um and 32um spatial grid aggregation schematic](../assets/high_resolution_st_exploration/cell_level_CosMx_portrait_nested_16um_32um_grid.pdf)

Grid aggregation can increase the median genes per spot or bin, which improves the stability of expression-based CNA inference.
The tradeoff is lower effective spatial resolution.
Record the grid size or neighborhood rule used for aggregation, and use the aggregated dataset as a separate prepared input.

Aggregation changes the effective spatial resolution and should be documented with the chosen neighborhood size or binning rule.
Do not compare aggregated and non-aggregated results as if they were generated at the same resolution.

## 4. Choose The Analysis Path

High-resolution datasets can be used in two common modes.

### CNA Inference

Use this path when the goal is to infer genome-wide CNA profiles from spatial expression.

Recommended starting methods:

- `CopyKAT`
- `SCEVAN`
- `STARCH`

These methods are good first choices for exploratory high-resolution ST analysis because they are computationally efficient, reliable in typical public runs, and practical for application-scale datasets.

Run selected methods:

```bash
st-cnabench --steps run \
  --data-config data.yaml \
  --model-config models.yaml \
  --prep-ids <DATASET_ID> \
  --exec-mode conda \
  --models CopyKAT SCEVAN STARCH
```

### Tumor Classification

Use this path when the immediate goal is to classify tumor and normal units rather than compare full CNA profiles.

The same efficient starting methods are recommended:

- `CopyKAT`
- `SCEVAN`
- `STARCH`

If reference-normal labels are available, set `ref_norm: true` and provide `raw.tumor_normal`.
If no reference-normal labels are available, use `tumor_normal_mode: de_novo` and run only methods that support reference-free operation.

## 5. Evaluate The Results

For CNA profile evaluation, run:

```bash
st-cnabench --steps eval \
  --data-config data.yaml \
  --eval-config eval.yaml \
  --prep-ids <DATASET_ID> \
  --models CopyKAT SCEVAN STARCH \
  --eval-tasks cna_profile
```

For tumor-normal classification evaluation, run:

```bash
st-cnabench --steps eval \
  --data-config data.yaml \
  --eval-config eval.yaml \
  --prep-ids <DATASET_ID> \
  --models CopyKAT SCEVAN STARCH \
  --eval-tasks tumor_normal
```

Only run evaluation tasks for which the dataset has the required GT files.
For example, `cna_profile` needs `raw.cna_gt`, while `tumor_normal` needs `raw.tumor_normal_gt`.

## 6. Recommended First-Pass Checklist

- Run `prep` and verify the standardized output bundle.
- Check median detected genes per unit.
- Aggregate upstream with a spatial grid if the median gene coverage is too low.
- Start with `CopyKAT`, `SCEVAN`, and `STARCH`.
- Run either `cna_profile` or `tumor_normal` evaluation depending on available GT.
- Inspect output logs before scaling to all datasets.

## Additional Figure Suggestions

For a complete tutorial page, useful figures include:

- a histogram of detected genes per unit before aggregation
- a histogram of detected genes per unit after aggregation
- a tumor-normal prediction map
- a CNA karyogram or CNA profile summary plot

Place figures under:

```text
docs/assets/high_resolution_st_exploration/
```

Then reference them from this page with:

```markdown
![Median genes QC](../assets/high_resolution_st_exploration/median_genes_qc.png)
```

## Try Next

- For the packaged cSCC demo, go to [Quickstart Demo And Expected Outputs](quickstart_demo.md)
- For dataset-level pre-run checks, go to [Dataset Properties Quick Check](dataset_properties_quick_check.md)
- For the CNA profile task example, go to [CNA Profile Task Example](cna_profile_hcc2t.md)
- For the tumor-normal task example, go to [Tumor-Normal Classification Task Example](tumor_normal_hcc2t.md)
- To adapt the workflow to your own data, go to [Use Your Own Dataset](use_your_own_dataset.md)
