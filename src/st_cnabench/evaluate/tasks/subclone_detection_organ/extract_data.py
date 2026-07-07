import pandas as pd
import os
import logging
from ...utils.io import save_tsv

def run_extract_data(eval_list, gt_loader, gene_annot_path, ds_tn_annot_path, ds_subcluster_annot_path=None, tn_annot_path_model_run=None, beads_mapping_path=None, ds_spatial_coord=None):
    
    results = {}
    clone_cna_profiles = {}

    # Load GT subcluster labels when available.
    if ds_subcluster_annot_path is not None:
        results['GT'] = gt_loader.extract_subcluster_preds(ds_subcluster_annot_path)

    # Collect predicted clone labels and clone-level CNA profiles from each model.
    for loader in eval_list:
        model_name = loader.model_name
        logging.info(f"Extracting subcluster preds and clone CNA profiles for model={model_name}")
        
        preds = loader.extract_subcluster_preds()
        if preds is not None and not preds.empty:
            results[model_name] = preds
            
            try:
                cna_profile = loader.extract_clone_cna_profile(gene_annot_path, ds_tn_annot_path)
                if cna_profile is not None and not cna_profile.empty:
                    clone_cna_profiles[f'{model_name}_Preds'] = cna_profile
            except Exception as e:
                logging.warning(f"Failed to extract clone CNA profile for {model_name}: {e}")
                
        else:
            logging.warning(f"Model {model_name} does not support subcluster extraction or failed, dropping from results.")
    
    # Merge all spot-level labels into one comparison table.
    df_all = None
    if len(results) > 0:
        for model_name, df in results.items():
            df = df.rename(columns={'Label_preds': f'{model_name}_Preds'})
            if df_all is None:
                df_all = df
            else:
                df_all = pd.merge(df_all, df, on='Barcodes', how='inner')
    
    # Attach spatial coordinates when provided.
    if ds_spatial_coord is not None and os.path.exists(ds_spatial_coord):
        df_coords = pd.read_csv(ds_spatial_coord, header=0)
        df_coords = df_coords.drop(columns=['in_tissue', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], errors='ignore')
        df_coords = df_coords.rename(columns={'barcode': 'Barcodes', 'array_row': 'y', 'array_col': 'x'})
        df_all = df_all.merge(df_coords, left_on='Barcodes', right_on='Barcodes', how='left')
        
    return df_all, clone_cna_profiles
