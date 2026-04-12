import os
import logging

from .tasks.clonal_evolution.extract_data import run_extract_data
from .tasks.clonal_evolution.mapping import run_clonal_mapping
from .tasks.clonal_evolution.tree_inference import run_tree_inference
from .tasks.clonal_evolution.plot import run_plot_phylogeography


def run_clone_evolution(
    eval_list,
    gene_annot_path,
    ds_tn_annot_path,
    ds_spatial_coord=None,
    ref_genome='hg38',
    bin_size=100000,
    tree_mode='parsimony',
    subset_tolerance=0.05,
    root_clone=None,
):
    """
    Run clonal evolution inference and spatial phylogeography plots.
    """
    logging.info("=" * 60)
    logging.info("Starting clonal evolution and spatial phylogeography task")
    logging.info("=" * 60)

    task_out_dir = eval_list[0].save_dir.parent / "subclone_detection"
    os.makedirs(task_out_dir, exist_ok=True)

    logging.info("[Step 1/4] Extracting clone labels and clone-level CNV profiles...")
    df_all, clone_cnv_profiles = run_extract_data(
        eval_list=eval_list,
        gene_annot_path=gene_annot_path,
        ds_tn_annot_path=ds_tn_annot_path,
        ds_spatial_coord=ds_spatial_coord
    )

    if not clone_cnv_profiles:
        logging.error("No valid clone CNV profiles extracted. Aborting Clonal Evolution task.")
        return

    logging.info(f"[Step 2/4] Mapping clones to {bin_size}bp bins and aligning feature spaces...")
    aligned_results = run_clonal_mapping(
        clone_cnv_profiles=clone_cnv_profiles,
        ref_genome=ref_genome,
        bin_size=bin_size
    )

    if not aligned_results:
        logging.error("Mapping resulted in empty aligned profiles. Aborting.")
        return

    logging.info("[Step 3/4] Inferring Evolutionary Trees...")
    trees = run_tree_inference(
        aligned_results,
        tree_mode=tree_mode,
        subset_tolerance=subset_tolerance,
        root_clone=root_clone,
    )

    if not trees:
        logging.error("Tree inference failed for all models. Aborting.")
        return

    logging.info("[Step 4/4] Plotting Spatial Phylogeography on Tissue Coordinates...")
    run_plot_phylogeography(
        df_all=df_all,
        trees=trees,
        out_dir=task_out_dir
    )

    logging.info("Clonal evolution and phylogeography task completed successfully.")
    logging.info("=" * 60)
