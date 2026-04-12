import logging
from .tasks.computational_efficiency.extract_data import run_extract_efficiency_data
from .tasks.computational_efficiency.plot import run_plot

def run_computational_efficiency(eval_list):
    """
    Eval task: Computational Efficiency
    level: RunTime and Memory Usage Comparison
    Notes: if use docker image to run, this task do not support.
    """
    if not eval_list:
        logging.warning("Eval list is empty, skipping efficiency task.")
        return
    logging.info("--- Starting Computational Efficiency Benchmarking Pipeline ---")
    
    # Extract
    logging.info("Extracting efficiency data...")
    df_efficiency = run_extract_efficiency_data(eval_list)

    if df_efficiency.empty:
        logging.error("No efficiency data extracted. Plotting aborted.")
        return

    base_save_dir = eval_list[0].save_dir.parent / "computational_efficiency"
    base_save_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = base_save_dir / "computational_efficiency_summary.tsv"
    plot_path = base_save_dir / "computational_efficiency_comparison.png"

    logging.info("Generating summary table...")
    df_efficiency.to_csv(tsv_path, index=False, sep='\t')

    logging.info("Generating efficiency comparison plot...")
    run_plot(df_efficiency, plot_path)
    
    logging.info("Computational efficiency task completed.")