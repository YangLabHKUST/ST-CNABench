import logging
import pandas as pd
import numpy as np
from .tasks.subclone_detection_in_slice.extract_data import run_extract_data
from .tasks.subclone_detection_in_slice.metrics import run_metrics
from .tasks.subclone_detection_in_slice.mapping import run_clonal_mapping 

from .tasks.subclone_detection_in_slice.plot import (
    run_spatial_plot_no_img, 
    run_plot_karyotype_composition, 
    run_plot_metric_heatmap, 
    run_plot_metric_barplot,
    run_plot_spatial_coherence,
    run_plot_clone_matching_heatmap,
)


def _sanitize_subclone_name(name: str) -> str:
    """Convert GT subclone labels to filesystem-friendly keys."""
    clean = ''.join(ch for ch in str(name) if ch.isalnum())
    return clean if clean else str(name)


def _resolve_spatial_method(ds_platform: str | None) -> str:
    platform = str(ds_platform).strip().upper()
    return "knn" if platform == "ST" else "distance"


def _export_best_scores_by_gt_subclone(matrix_dict, metric_label, save_dir):
    """
    Export per-GT-subclone best matching scores from clone-matching matrices.

    For each model, each GT subclone row takes max across predicted clones (axis=1),
    matching the same extraction logic used for global max score aggregation.
    """
    if not matrix_dict:
        return

    records = []
    for model_name, score_matrix in matrix_dict.items():
        if score_matrix is None or score_matrix.empty:
            continue

        row_best = score_matrix.max(axis=1)
        for gt_subclone, best_score in row_best.items():
            records.append({
                "Model": model_name,
                "GT_Subclone": str(gt_subclone),
                "Best_Score": float(best_score)
            })

    if not records:
        return

    df_long = pd.DataFrame(records)
    df_wide = df_long.pivot(index="Model", columns="GT_Subclone", values="Best_Score").sort_index()

    metric_key = metric_label.lower().replace("-", "_")
    df_wide.to_csv(save_dir / f"best_{metric_key}_per_subclone_scores.csv")

    for subclone in df_wide.columns:
        safe_subclone = _sanitize_subclone_name(subclone).lower()
        series_name = f"{_sanitize_subclone_name(subclone)}_Best_{metric_label.replace('-', '_')}_Score"
        df_wide[subclone].to_csv(
            save_dir / f"{safe_subclone}_best_{metric_key}_scores.csv",
            header=[series_name]
        )

def run_subclone_detection_in_slice(
    eval_list,
    gt_loader,
    cnv_gt_path,
    gene_annot_path,
    ds_subcluster_annot_path=None,
    tn_annot_path_model_run=None,
    ds_spatial_coord=None,
    threshold=2.1,
    ds_HE_img=None,
    ds_scalefactors=None,
    beads_mapping_path=None,
    ds_platform=None,
):
    if not eval_list:
        logging.warning("Eval list is empty, skipping subclone detection task.")
        return

    logging.info("--- Starting Subclone Detection (In Slice) Benchmarking Pipeline ---")
    base_save_dir = eval_list[0].save_dir.parent / "subclone_detection"
    base_save_dir.mkdir(parents=True, exist_ok=True)
    spatial_method = _resolve_spatial_method(ds_platform)

    logging.info("Extracting subclone prediction data and clone CNV profiles...")
    results, clone_cnv_profiles = run_extract_data(
        eval_list, 
        gt_loader, 
        gene_annot_path=gene_annot_path, 
        ds_subcluster_annot_path=ds_subcluster_annot_path, 
        tn_annot_path_model_run=tn_annot_path_model_run, 
        cnv_gt_path=cnv_gt_path,
        beads_mapping_path=beads_mapping_path
    )

    if results is not None:
        results_csv_path = base_save_dir / "subclone_predictions.tsv"
        results.to_csv(results_csv_path, sep='\t', index=False)
        logging.info(f"Subclone prediction results saved to {results_csv_path}")
    else:
        logging.error("No subclone prediction results to save!")
        return

    logging.info("Aligning clone CNV profiles to unified 100kb genomic bins...")
    aligned_clone_profiles = run_clonal_mapping(
        clone_cnv_profiles, 
        ref_genome='hg38', 
        bin_size=100000, 
        fill_value=np.nan
    )

    logging.info("Calculating comprehensive subclone detection metrics...")
    logging.info("Using spatial coherence mode: %s", spatial_method)
    metrics_results = run_metrics(
        results,
        aligned_clone_profiles=aligned_clone_profiles,
        spatial_method=spatial_method,
        threshold=threshold,
    )
    
    if metrics_results[0] is not None:
        (ari_matrix, ari_scores, v_measure_matrix, v_measure_scores, 
         nmi_matrix, nmi_scores, spatial_coherence_scores, 
         pcc_matrices, max_pcc_scores, f1_matrices, max_f1_scores) = metrics_results

        ari_matrix.to_csv(base_save_dir / "ari_matrix.csv")
        v_measure_matrix.to_csv(base_save_dir / "v_measure_matrix.csv")
        nmi_matrix.to_csv(base_save_dir / "nmi_matrix.csv")
        
        pd.Series(ari_scores).to_csv(base_save_dir / "ari_scores.csv", header=['ARI_Score'])
        pd.Series(v_measure_scores).to_csv(base_save_dir / "v_measure_scores.csv", header=['V_Measure_Score'])
        pd.Series(nmi_scores).to_csv(base_save_dir / "nmi_scores.csv", header=['NMI_Score'])

        if max_pcc_scores:
            pd.Series(max_pcc_scores).to_csv(base_save_dir / "max_pcc_scores.csv", header=['Max_PCC_Score'])
        if max_f1_scores:
            pd.Series(max_f1_scores).to_csv(base_save_dir / "max_f1_scores.csv", header=['Max_F1_Score'])
        if pcc_matrices:
            _export_best_scores_by_gt_subclone(pcc_matrices, "PCC", base_save_dir)
        if f1_matrices:
            _export_best_scores_by_gt_subclone(f1_matrices, "F1", base_save_dir)
        if spatial_coherence_scores:
            sc_df = pd.DataFrame(spatial_coherence_scores).T
            sc_df.to_csv(base_save_dir / "spatial_coherence_scores.csv")

    spatial_plot_path = base_save_dir / "subclone_prediction_comparison.png"
    logging.info("Generating subclone prediction spatial plot...")
    run_spatial_plot_no_img(results, spatial_plot_path)

    logging.info("Generating unified clustered karyogram and aligned purity plots...")
    run_plot_karyotype_composition(results, aligned_clone_profiles, base_save_dir)
    
    logging.info("Generating quantitative metric plots...")
    if metrics_results[0] is not None:
        metrics_data = [
            (ari_matrix, ari_scores, "ARI"),
            (v_measure_matrix, v_measure_scores, "V-Measure"),
            (nmi_matrix, nmi_scores, "NMI")
        ]
        for matrix, scores, name in metrics_data:
            run_plot_metric_heatmap(matrix, name, base_save_dir)
            run_plot_metric_barplot(scores, name, base_save_dir)

        if max_pcc_scores:
            run_plot_metric_barplot(max_pcc_scores, "Max-PCC", base_save_dir)
        if max_f1_scores:
            run_plot_metric_barplot(max_f1_scores, "Max-F1", base_save_dir)

        if pcc_matrices:
            run_plot_clone_matching_heatmap(pcc_matrices, "PCC", base_save_dir)
        if f1_matrices:
            run_plot_clone_matching_heatmap(f1_matrices, "Max-F1", base_save_dir)

        if spatial_coherence_scores:
            run_plot_spatial_coherence(spatial_coherence_scores, base_save_dir)

    logging.info("Subclone detection (in slice) task completed.")
