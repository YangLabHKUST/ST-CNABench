# Dataset Properties Quick Check

This page summarizes dataset-level checks that should be reviewed before running CNA inference or tumor-normal classification.
The goal is to identify datasets that are likely to support reliable inference and to document the assumptions used for each run.

## 1. Expression Coverage

Expression-based CNA inference depends on enough detected genes per analysis unit.
For spatial transcriptomics data, first check the median detected genes per spot, bead, cell, or aggregated bin.

We recommend using approximately 1,000 detected genes per analysis unit as a practical lower target for first-pass CNA inference.
Datasets below this level can still be useful, but sparse expression coverage makes genome-wide CNA profiles less stable and should be handled carefully.

Example check after `prep`:

```python
import numpy as np
import scanpy as sc

adata = sc.read_h5ad("<output.root>/filtered_feature_bc_matrix.h5ad")
detected_genes = np.asarray((adata.X > 0).sum(axis=1)).ravel()

print(f"Median detected genes per unit: {np.median(detected_genes):.1f}")
print(f"Mean detected genes per unit: {np.mean(detected_genes):.1f}")
```

Recommended interpretation:

- Median detected genes near or above 1,000: proceed to model running.
- Median detected genes clearly below 1,000: consider upstream spatial aggregation or a tumor-normal classification task before full CNA profile evaluation.
- After aggregation, rerun the same expression coverage check on the aggregated input.

## 2. Allelic Coverage For cnLOH Analysis

If the analysis includes copy-neutral loss of heterozygosity, run `cellSNP` or `cellSNP-lite` before model inference and inspect the allelic evidence.
cnLOH inference requires enough heterozygous SNPs and enough allele-informative UMI coverage.

Before running cnLOH-aware analysis, check:

- The number of retained heterozygous SNPs is greater than 10,000.
- The allele-informative UMI coverage is sufficient for the selected analysis units.
- As a practical lower target, use more than 200 allele-informative UMI counts after SNP and barcode filtering.

Suggested outputs to inspect from `cellSNP` or `cellSNP-lite` include:

- `cellSNP.base.vcf`
- `cellSNP.samples.tsv`
- `cellSNP.tag.AD.mtx`
- `cellSNP.tag.DP.mtx`

If the heterozygous SNP count or allele-informative UMI coverage is too low, cnLOH calls should be treated as low confidence.
In that case, document the limitation and prioritize expression-based CNA tasks unless additional sequencing depth or matched genotyping information is available.

## 3. Cancer Type And Expected Genomic Instability

Before interpreting CNA results, check whether the cancer type is expected to show strong genomic instability.
Tumors with frequent arm-level or chromosome-level alterations are usually better suited for expression-based CNA benchmarking than tumors with mostly focal, rare, or low-burden alterations.

Recommended checks:

- Review whether the cancer type is known to have recurrent broad CNA events.
- Check whether public references, cohort studies, or matched genomic data support high CNA burden in this tumor type.
- If a cancer type is expected to be genomically quiet, interpret weak CNA signals cautiously.

This check does not exclude a dataset from analysis, but it should be recorded because it affects the expected signal strength and the interpretation of negative results.

## 4. Reference Set Quality

For reference-normal methods, a high-quality reference set is critical.
The reference set should be small enough to be manually reviewed, but clean enough to represent non-malignant expression without obvious tumor contamination.

We recommend starting from cell-type annotation to separate malignant and benign compartments, then double-checking the candidate reference cells or spots with additional evidence.
Use the evidence that is available for the dataset.

Useful checks include:

- Known tumor marker expression: verify that tumor components express expected tumor markers and that candidate reference cells do not show strong tumor marker signal.
- Tissue morphology: if H&E or matched histology is available, map spots or cells back to the tissue image and check whether candidate reference regions are morphologically consistent with benign tissue.
- Spatial context: inspect whether candidate reference cells or spots are spatially separated from obvious tumor regions when the tissue structure supports this.
- Metadata consistency: confirm that sample labels, tissue regions, and cell-type annotations agree with the proposed reference set.

After review, save the selected high-quality reference group as the reference input for model inference.
Avoid using a broad unreviewed normal pool if it may contain tumor-adjacent or tumor-contaminated units.

## Recommended Pre-Run Checklist

- Check median detected genes per analysis unit and record the value.
- Use approximately 1,000 detected genes per unit as a practical target for reliable first-pass CNA inference.
- If cnLOH is needed, run `cellSNP` or `cellSNP-lite` and confirm more than 10,000 heterozygous SNPs with sufficient allele-informative UMI coverage.
- Check whether the cancer type is expected to carry strong genomic instability.
- Build a high-quality reference set using cell type, tumor markers, histology, spatial context, and metadata when available.
- Document any limitation before comparing methods or interpreting negative results.

## Try Next

- For high-resolution spatial data, go to [High-Resolution ST Exploration](high_resolution_st_exploration.md)
- To adapt the workflow to your own data, go to [Use Your Own Dataset](use_your_own_dataset.md)
- For model execution details, go to [Model Run](../model_run.md)
