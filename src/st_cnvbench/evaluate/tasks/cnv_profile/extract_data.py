import logging
import pandas as pd

def run_extract_data(eval_list, gt_loader, gt_path, gene_annot_path, ds_tn_annot_path):
    # Keep CN and LOH tracks separate from the beginning.
    results = {'CN_Profile': {}, 'LOH_Profile': {}}

    logging.info(f"Extracting Ground Truth for {gt_loader.eval_name}")
    gt_df = gt_loader.extract_cnv_profile(gt_path=gt_path)
    
    gt_cols_to_keep = ['Chromosome', 'Start', 'End', 'ID', 'GT_Score', 'GT_Event']
    
    # Keep focal-event annotations when the GT exporter provides them.
    for focal_tag in ['Is_Focal_Type1', 'Is_Focal_Type2']:
        if focal_tag in gt_df.columns:
            gt_cols_to_keep.append(focal_tag)
            
    results['CN_Profile']['GT'] = gt_df[gt_cols_to_keep].copy()
    if 'LOH_Status' in gt_df.columns:
        loh_gt = gt_df[['Chromosome', 'Start', 'End', 'ID', 'LOH_Status']].rename(columns={'LOH_Status': 'GT_LOH_Status'})
        results['LOH_Profile']['GT'] = loh_gt

    for loader in eval_list:
        logging.info(f"Extracting data for {loader.eval_name}")
        out_info = loader.extract_cnv_profile(gene_annot_path, ds_tn_annot_path)
        eval_name = loader.eval_name
        
        if out_info is not None and not out_info['data'].empty:
            model_df = out_info['data']
            
            # Track 1: total copy-number score.
            if 'CN_Score' in model_df.columns:
                results['CN_Profile'][eval_name] = model_df[['Chromosome', 'Start', 'End', 'ID', 'CN_Score']]
                
            # Track 2: LOH status, only when the model emits usable values.
            if 'LOH_Status' in model_df.columns and model_df['LOH_Status'].notna().any():
                results['LOH_Profile'][eval_name] = model_df[['Chromosome', 'Start', 'End', 'ID', 'LOH_Status']]
                logging.info(f"Model {eval_name} has valid LOH signals, added to LOH_Profile.")
            else:
                logging.info(f"Model {eval_name} is pure expression-based (no LOH), skipped for LOH_Profile.")
                
        else:
            logging.warning(f"Model {eval_name} returned None/Empty data.")

    return results
