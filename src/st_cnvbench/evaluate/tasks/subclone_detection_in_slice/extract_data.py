import pandas as pd
import os
import logging
from ...utils.io import save_tsv

def run_extract_data(eval_list, gt_loader, gene_annot_path=None, cnv_gt_path=None, ds_subcluster_annot_path=None, tn_annot_path_model_run=None, beads_mapping_path=None):
    
    results = {}
    clone_cnv_profiles = {} 
    
    # 1. Extract GT subcluster label
    if ds_subcluster_annot_path is not None:
        results['GT'] = gt_loader.extract_subcluster_preds(ds_subcluster_annot_path)

    # Load GT clone-level CNV profiles when they are provided.
    if cnv_gt_path is not None:
        logging.info(f"Extracting clone CNV profiles for GT from {cnv_gt_path}")
        gt_cnv_profile = gt_loader.extract_clone_cnv_profile(cnv_gt_path=cnv_gt_path)
        if gt_cnv_profile is not None and not gt_cnv_profile.empty:
            clone_cnv_profiles['GT'] = gt_cnv_profile
        else:
            logging.warning("[GT] No clone CNV profile extracted.")

    # Collect model predictions and clone-level CNV profiles.
    for loader in eval_list:
        model_name = loader.model_name
        logging.info(f"Extracting subcluster preds and clone CNV profiles for model={model_name}")
        
        preds = loader.extract_subcluster_preds()
        if preds is None:
            logging.warning(f"Model {model_name} does not support subcluster extraction, dropping from results.")
            continue
        results[model_name] = preds
        
        cnv_profile = loader.extract_clone_cnv_profile(gene_annot_path=gene_annot_path, tn_annot_path=tn_annot_path_model_run)
        if cnv_profile is not None and not cnv_profile.empty:
            clone_cnv_profiles[model_name] = cnv_profile
        else:
            logging.warning(f"[{model_name}] No clone CNV profile extracted.")

    # Merge all spot-level labels into one table.
    df_all = None
    if len(results) > 0:
        for model_name, df in results.items():
            df = df.rename(columns={'Label_preds': f'{model_name}_Preds'})
            if df_all is None:
                df_all = df
            else:
                df_all = pd.merge(df_all, df, on='Barcodes', how='inner')
    
    # Remove reference normal cells when the model-run annotation is available.
    if tn_annot_path_model_run is not None and df_all is not None:
        df_model_run = pd.read_csv(tn_annot_path_model_run, sep='\t', header=0, comment='#')
        df_model_run.rename(columns={df_model_run.columns[0]: 'Barcodes', df_model_run.columns[1]: 'tumor_normal'}, inplace=True)
        normal_mask = df_model_run['tumor_normal'].astype(str).str.strip().str.lower() == 'normal'
        normal_barcode = set(df_model_run[normal_mask]['Barcodes'])
        df_all = df_all[~df_all['Barcodes'].isin(normal_barcode)]

    # Broadcast meta-bead predictions back to raw beads for Slide-DNA-seq style inputs.
    if beads_mapping_path is not None and df_all is not None:
        df_beads_mapping = pd.read_csv(beads_mapping_path, header=0)
        df_all = df_all.merge(df_beads_mapping, left_on='Barcodes', right_on='pseudo_barcode', how='left')
        df_all['Barcodes'] = df_all['original_barcode']
        
    return df_all, clone_cnv_profiles
