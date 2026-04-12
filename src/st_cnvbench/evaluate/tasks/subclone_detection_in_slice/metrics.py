import pandas as pd
import os
import logging
import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import adjusted_rand_score, v_measure_score, normalized_mutual_info_score

from ...utils.metrics import calculate_spatial_coherence, calc_max_f1_sweep

def run_metrics(
    results,
    aligned_clone_profiles=None,
    spatial_method="distance",
    spatial_k=11,
    threshold=2.1,
    gt_log2_threshold=0.09,
):
    """
    Evaluate subclone detection from cellular, spatial, and genomic views.

    The default GT log2 threshold of 0.09 is retained for the current
    SlideDNAseq-specific genomic evaluation workflow.
    """
    if 'GT_Preds' not in results.columns:
        logging.error("GT data missing! Metrics cannot be calculated.")
        return [None] * 11 
    
    model_names = [name for name in results.columns if name not in ['Barcodes','pseudo_barcode','original_barcode','group','x','y']]
    
    logging.info("Calculating standard clustering metrics...")
    ari_matrix = pd.DataFrame(index=model_names, columns=model_names, dtype=float)
    v_measure_matrix = pd.DataFrame(index=model_names, columns=model_names, dtype=float)
    nmi_matrix = pd.DataFrame(index=model_names, columns=model_names, dtype=float)
    
    for i in range(len(model_names)):
        for j in range(i, len(model_names)):
            mi, mj = model_names[i], model_names[j]
            labels_i = results[mi].astype(str).values
            labels_j = results[mj].astype(str).values
            
            ari_matrix.loc[mi, mj] = ari_matrix.loc[mj, mi] = adjusted_rand_score(labels_i, labels_j)
            v_measure_matrix.loc[mi, mj] = v_measure_matrix.loc[mj, mi] = v_measure_score(labels_i, labels_j)
            nmi_matrix.loc[mi, mj] = nmi_matrix.loc[mj, mi] = normalized_mutual_info_score(labels_i, labels_j)

    ari_scores, v_measure_scores, nmi_scores = {}, {}, {}
    gt_labels = results['GT_Preds'].astype(str).values
    for model in model_names:
        if model == 'GT_Preds': continue
        pred_labels = results[model].astype(str).values
        ari_scores[model] = adjusted_rand_score(gt_labels, pred_labels)
        v_measure_scores[model] = v_measure_score(gt_labels, pred_labels)
        nmi_scores[model] = normalized_mutual_info_score(gt_labels, pred_labels)

    spatial_coherence_scores = {}
    if 'x' in results.columns and 'y' in results.columns:
        logging.info(f"Calculating Spatial Coherence using {spatial_method.upper()}...")
        coords = results[['x', 'y']].values
        
        for model in model_names:
            valid_idx = results[model].notna()
            if not valid_idx.any(): continue
            
            pts = coords[valid_idx]
            labs = results.loc[valid_idx, model].values
            
            try:
                obs, z = calculate_spatial_coherence(
                    pts, labs, 
                    method=spatial_method, 
                    k=spatial_k, 
                    threshold=threshold,
                    calculate_z_score=True, 
                    n_permutations=100
                )
                spatial_coherence_scores[model] = {'observed_score': obs, 'z_score': z}
            except Exception as e:
                logging.warning(f"Spatial coherence failed for {model}: {e}")

    pcc_matrices = {}   
    max_pcc_scores = {} 
    f1_matrices = {}
    max_f1_scores = {}

    if aligned_clone_profiles is not None and 'GT' in aligned_clone_profiles:
        logging.info("Calculating Genomic Profile Max PCC and Max F1 Sweep...")
        gt_df = aligned_clone_profiles['GT']
        
        gt_pivot = gt_df.pivot(index=['Chromosome', 'Start', 'End'], columns='Clone_ID', values='CN_Score')
        gt_subclones = [c for c in gt_pivot.columns if 'normal' not in str(c).lower()]
        
        for model in model_names:
            mod_key = model.replace('_Preds', '') 
            if mod_key == 'GT' or mod_key not in aligned_clone_profiles: continue
            
            mod_df = aligned_clone_profiles[mod_key]
            if mod_df.empty: continue
            
            mod_pivot = mod_df.pivot(index=['Chromosome', 'Start', 'End'], columns='Clone_ID', values='CN_Score')
            
            common_bins = gt_pivot.index.intersection(mod_pivot.index)
            if len(common_bins) == 0: continue
            
            gt_mat = gt_pivot.loc[common_bins, gt_subclones]
            mod_mat = mod_pivot.loc[common_bins]
            
            pcc_matrix = pd.DataFrame(index=gt_subclones, columns=mod_mat.columns, dtype=float)
            f1_matrix = pd.DataFrame(index=gt_subclones, columns=mod_mat.columns, dtype=float)
            
            for gt_c in gt_subclones:
                gt_cont = gt_mat[gt_c].fillna(0).values
                gt_events = np.zeros_like(gt_cont)
                gt_events[gt_cont > gt_log2_threshold] = 1
                gt_events[gt_cont < -gt_log2_threshold] = -1
                
                for mod_c in mod_mat.columns:
                    mod_cont = mod_mat[mod_c].fillna(0).values
                    
                    r, _ = pearsonr(gt_cont, mod_cont)
                    pcc_matrix.loc[gt_c, mod_c] = r
                    
                    try:
                        f1_res = calc_max_f1_sweep(gt_events, mod_cont)
                        f1_matrix.loc[gt_c, mod_c] = f1_res['Max_Macro_F1']
                    except Exception as e:
                        logging.warning(f"F1 Sweep failed for {model} clone {mod_c}: {e}")
                        f1_matrix.loc[gt_c, mod_c] = 0.0
            
            pcc_matrices[model] = pcc_matrix
            max_pcc_scores[model] = pcc_matrix.max(axis=1).mean()
            
            f1_matrices[model] = f1_matrix
            max_f1_scores[model] = f1_matrix.max(axis=1).mean()

    return (ari_matrix, ari_scores, v_measure_matrix, v_measure_scores, 
            nmi_matrix, nmi_scores, spatial_coherence_scores, 
            pcc_matrices, max_pcc_scores, f1_matrices, max_f1_scores)
