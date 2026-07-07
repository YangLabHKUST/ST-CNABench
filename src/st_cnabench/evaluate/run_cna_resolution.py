import logging
from .tasks.cna_resolution.extract_data import run_extract_data
from .tasks.cna_resolution.plot import run_plot_resolution

def run_cna_resolution(eval_list, gene_annot_path=None, ds_tn_annot_path=None):
    if not eval_list:
        logging.warning("Eval list is empty, skipping efficiency task.")
        return

    logging.info("--- Starting CNA Resolution Benchmarking Pipeline ---")

    logging.info("Extracting CNA resolution data...")
    results = run_extract_data(eval_list, gene_annot_path, ds_tn_annot_path)

    if not results:
        logging.error("No data extracted. Plotting aborted.")
        return
    plot_dir = eval_list[0].save_dir.parent / "cna_resolution"
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_dir / "cna_resolution_comparison.png"

    logging.info("Generating CNA resolution comparison plot...")
    run_plot_resolution(results, plot_path)

    logging.info("CNA resolution task completed.")
