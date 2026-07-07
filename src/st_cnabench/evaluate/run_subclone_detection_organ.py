import logging
import pandas as pd
from .tasks.subclone_detection_organ.extract_data import run_extract_data
from .tasks.subclone_detection_organ.metrics import run_metrics
from .tasks.subclone_detection_organ.mapping import run_clonal_mapping
from .tasks.subclone_detection_organ.plot import (
    run_plot_karyotype_composition,
    run_plot_metric_barplot,
    run_plot_metric_heatmap,
    run_spatial_plot_no_img,
)

def run_subclone_detection_organ(eval_list, gt_loader, gene_annot_path, ds_tn_annot_path, ds_subcluster_annot_path=None, tn_annot_path_model_run=None, ds_spatial_coord=None, threshold=2.1, ds_HE_img=None, ds_scalefactors=None, beads_mapping_path=None):
    if not eval_list:
        logging.warning("Eval list is empty, skipping subclone detection task.")
        return

    logging.info("--- Starting Subclone Detection Benchmarking Pipeline ---")
    base_save_dir = eval_list[0].save_dir.parent / "subclone_detection"
    base_save_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Extracting subclone prediction data...")
    results, clone_cna_profiles = run_extract_data(
        eval_list,
        gt_loader,
        gene_annot_path,
        ds_tn_annot_path,
        ds_subcluster_annot_path,
        tn_annot_path_model_run,
        beads_mapping_path=beads_mapping_path,
        ds_spatial_coord=ds_spatial_coord,
    )
    if results is not None:
        results_csv_path = base_save_dir / "subclone_predictions.tsv"
        results.to_csv(results_csv_path, sep="\t", index=False)
        logging.info(f"Subclone prediction results saved to {results_csv_path}")
    else:
        logging.error("No subclone prediction results to save!")
        return

    aligned_clone_profiles = run_clonal_mapping(
        clone_cna_profiles,
        ref_genome='hg38',
        bin_size=100000,
    )

    spatial_plot_path = base_save_dir / "subclone_prediction_comparison.png"
    
    logging.info("Generating subclone prediction spatial plot...")
    run_spatial_plot_no_img(results, spatial_plot_path)

    logging.info("Generating unified clustered karyogram and aligned purity plots...")
    run_plot_karyotype_composition(results, aligned_clone_profiles, base_save_dir)

    logging.info("Calculating subclone detection metrics...")
    metrics_results = run_metrics(results)
    if metrics_results[0] is None:
        logging.error("Metrics calculation failed (GT labels missing or invalid).")
        return

    ari_matrix, ari_scores, v_measure_matrix, v_measure_scores, nmi_matrix, nmi_scores = metrics_results

    ari_matrix.to_csv(base_save_dir / "ari_matrix.csv")
    v_measure_matrix.to_csv(base_save_dir / "v_measure_matrix.csv")
    nmi_matrix.to_csv(base_save_dir / "nmi_matrix.csv")
    logging.info("Saved metric matrices to CSV.")

    pd.Series(ari_scores).to_csv(base_save_dir / "ari_scores.csv", header=['ARI_Score'])
    pd.Series(v_measure_scores).to_csv(base_save_dir / "v_measure_scores.csv", header=['V_Measure_Score'])
    pd.Series(nmi_scores).to_csv(base_save_dir / "nmi_scores.csv", header=['NMI_Score'])
    logging.info("Saved metric scores against GT to CSV.")

    logging.info("Generating quantitative metric plots...")
    metrics_data = [
        (ari_matrix, ari_scores, "ARI"),
        (v_measure_matrix, v_measure_scores, "V-Measure"),
        (nmi_matrix, nmi_scores, "NMI")
    ]
    
    for matrix, scores, name in metrics_data:
        run_plot_metric_heatmap(matrix, name, base_save_dir)
        run_plot_metric_barplot(scores, name, base_save_dir)

    logging.info("Subclone detection task completed.")
