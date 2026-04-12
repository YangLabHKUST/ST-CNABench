import logging
from .tasks.cnv_resolution.extract_data import run_extract_data
from .tasks.cnv_resolution.plot import run_plot_resolution

def run_cnv_resolution(eval_list, gene_annot_path=None, ds_tn_annot_path=None):
    if not eval_list:
        logging.warning("Eval list is empty, skipping efficiency task.")
        return

    logging.info("--- Starting CNV Resolution Benchmarking Pipeline ---")

    logging.info("Extracting CNV resolution data...")
    results = run_extract_data(eval_list, gene_annot_path, ds_tn_annot_path)

    if not results:
        logging.error("No data extracted. Plotting aborted.")
        return
    plot_dir = eval_list[0].save_dir.parent / "cnv_resolution"
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_dir / "cnv_resolution_comparison.png"

    logging.info("Generating CNV resolution comparison plot...")
    run_plot_resolution(results, plot_path)

    logging.info("CNV resolution task completed.")