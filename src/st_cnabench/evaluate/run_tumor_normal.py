import logging
from functools import reduce

import pandas as pd

from .tasks.tumor_normal.extract_data import run_extract_data
from .tasks.tumor_normal.metrics import run_metrics
from .tasks.tumor_normal.plot import run_spatial_plot, run_bar_plot
from .tasks.tumor_normal.plot import run_confusion_map


def _resolve_spatial_method(ds_platform: str | None) -> str:
    platform = str(ds_platform).strip().upper()
    return "knn" if platform == "ST" else "distance"


def _build_aligned_spot_predictions(results, ds_spatial_coord):
    if ds_spatial_coord is None:
        raise ValueError("ds_spatial_coord is required to export aligned spot predictions.")

    coords_df = pd.read_csv(ds_spatial_coord, index_col=0)
    index_list = [coords_df.index]

    for model_name, df in results.items():
        if df is None or df.empty:
            raise ValueError(f"Results for model '{model_name}' are empty, cannot export aligned predictions.")
        index_list.append(df.set_index(df.columns[0]).index)

    common_index = reduce(lambda x, y: x.intersection(y), index_list)
    if len(common_index) == 0:
        raise ValueError("No common barcodes found across coords and all methods for aligned prediction export.")

    aligned_df = pd.DataFrame({"Barcode": common_index})
    for model_name, res_df in results.items():
        preds = (
            res_df.set_index(res_df.columns[0]).iloc[:, 0]
            .loc[common_index]
            .astype(float)
            .astype(int)
        )
        aligned_df[model_name] = preds.values

    return aligned_df



def run_tumor_normal(
    eval_list,
    gt_loader,
    ds_tn_annot_path=None,
    ds_tn_annot_path_model_run=None,
    ds_tn_mode=None,
    ds_spatial_coord=None,
    ds_scalefactors=None,
    ds_HE_img=None,
    ds_platform=None,
    threshold=2.1,
    spatial_k=11,
):
    if not eval_list:
        logging.warning("Eval list is empty, skipping efficiency task.")
        return

    logging.info("--- Starting Tumor Normal Prediction Benchmarking Pipeline ---")
    base_save_dir = eval_list[0].save_dir.parent / "tumor_normal"
    base_save_dir.mkdir(parents=True, exist_ok=True)
    spatial_method = _resolve_spatial_method(ds_platform)

    logging.info("Extracting tumor/normal prediction data...")
    results = run_extract_data(eval_list, gt_loader, ds_tn_annot_path, ds_tn_annot_path_model_run, ds_tn_mode)

    logging.info("Saving aligned spot-level predictions across all methods...")
    aligned_spot_path = base_save_dir / "aligned_spot_predictions.tsv"
    aligned_spot_df = _build_aligned_spot_predictions(results, ds_spatial_coord)
    aligned_spot_df.to_csv(aligned_spot_path, sep="\t", index=False)
    logging.info(f"Aligned spot-level predictions saved to {aligned_spot_path}")

    logging.info("Calculating tumor/normal detection metrics...")
    logging.info("Using spatial coherence mode: %s", spatial_method)
    metrics_summary = run_metrics(
        results,
        ds_spatial_coord,
        spatial_method=spatial_method,
        spatial_k=spatial_k,
        threshold=threshold,
    )

    if metrics_summary is not None:
        metrics_csv_path = base_save_dir / "detection_metrics_summary.csv"
        metrics_summary.to_csv(metrics_csv_path, index=False)
        logging.info(f"Metrics saved to {metrics_csv_path}")

    spatial_plot_path = base_save_dir / "tumor_normal_prediction_comparison.png"

    logging.info("Generating tumor/normal prediction plot...")
    run_spatial_plot(results, ds_spatial_coord, ds_HE_img, ds_scalefactors, spatial_plot_path)

    bar_plot_path = base_save_dir / "tumor_normal_metrics_summary.pdf"
    run_bar_plot(metrics_summary, bar_plot_path)

    confusion_matrix_plot_path = base_save_dir / "tumor_normal_confusion_matrices.pdf"
    logging.info("Generating confusion matrix heatmaps...")
    run_confusion_map(results, confusion_matrix_plot_path)
    logging.info("Tumor normal task completed.")
