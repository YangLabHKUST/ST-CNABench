import logging
import os
from pathlib import Path

from .tasks.cna_profile.extract_data import run_extract_data
from .tasks.cna_profile.mapping import run_mapping
from .tasks.cna_profile.metrics import run_metrics
from .tasks.cna_profile.plot import run_bar_plot, run_karyogram_level_plot

def run_cna_profile(eval_list, gt_loader, gt_path, gene_annot_path, ds_tn_annot_path, ref_genome='hg38', bin_size=100000):
    """
    Eval task: CNA Profile
    Level: Dual-Track Unified Profile (0-centered CN Score & LOH Status)
    """
    logging.info("--- Starting CNA Profile Benchmarking Pipeline ---")

    # 1. Extract data
    logging.info("Extracting data...")
    formatted_results = run_extract_data(
        eval_list, 
        gt_loader, 
        gt_path, 
        gene_annot_path, 
        ds_tn_annot_path
    )

    # 2. Mapping
    logging.info("Mapping data...")
    aligned_results = run_mapping(
        results=formatted_results, 
        ref_genome=ref_genome, 
        bin_size=bin_size
    )

    plot_dir = Path(eval_list[0].save_dir).parent / "cna_profile"
    os.makedirs(plot_dir, exist_ok=True)

    # 3. Metrics
    logging.info("Calculating metrics...")
    metrics_df = run_metrics(
        aligned_results=aligned_results, 
        save_dir=plot_dir
    )

    # 4. Plot - Bar plots
    logging.info("Generating comprehensive bar plots...")
    run_bar_plot(
        df_metrics=metrics_df, 
        output_path=os.path.join(plot_dir, "cna_metrics_summary.pdf")
    )

    # 5. Plot - Karyogram plots
    logging.info("Generating karyogram plots...")
    # Keep the loop even though only 'Unified_Profile' is currently emitted.
    for profile_name, profile_df in aligned_results.items():
        if profile_df is not None and not profile_df.empty:
            logging.info(f"Generating Karyogram for profile: {profile_name}")
            karyo_filename = f"karyogram_{profile_name}.pdf"
            run_karyogram_level_plot(
                aligned_level_df=profile_df,
                level_name=profile_name,
                output_path=os.path.join(plot_dir, karyo_filename)
            )

    logging.info("--- CNA Profile Benchmarking Pipeline Completed ---")
    
    return metrics_df
