import pandas as pd
import os
import logging
import numpy as np
from sklearn.metrics import adjusted_rand_score, v_measure_score, normalized_mutual_info_score
from ...utils.metrics import calculate_spatial_coherence

def run_metrics(results):
    """
    results: df_all (Barcodes, GT_Preds, model1_Preds, model2_Preds, ...)
    Calculate ARI, V-Measure, and NMI between each two columns.
    Returns heatmaps (matrices) for inter-model comparisons and dictionary scores against GT.
    """
    if 'GT_Preds' not in results.columns:
        logging.error("GT data missing!")
        return None, None, None, None, None, None
    
    # Intersect
    logging.info("Building common index for alignment...")
    coords_df = results[['Barcodes', "x", "y"]].copy()
    coords_df = coords_df.set_index('Barcodes')
    
    model_names = [name for name in results.columns if name not in ['Barcodes','pseudo_barcode','original_barcode','group','x','y']]
    
    # Initialize the three pairwise metric matrices.
    ari_matrix = pd.DataFrame(index=model_names, columns=model_names, dtype=float)
    v_measure_matrix = pd.DataFrame(index=model_names, columns=model_names, dtype=float)
    nmi_matrix = pd.DataFrame(index=model_names, columns=model_names, dtype=float)
    
    for i in range(len(model_names)):
        for j in range(i, len(model_names)):
            model_i = model_names[i]
            model_j = model_names[j]
            
            # Cast to string to keep metric inputs stable even when labels are mixed.
            labels_i = results[model_i].astype(str).values
            labels_j = results[model_j].astype(str).values
            
            # Compute the three clustering metrics.
            ari = adjusted_rand_score(labels_i, labels_j)
            v_meas = v_measure_score(labels_i, labels_j)
            nmi = normalized_mutual_info_score(labels_i, labels_j)
            
            # Fill the matrices symmetrically.
            ari_matrix.loc[model_i, model_j] = ari
            ari_matrix.loc[model_j, model_i] = ari
            
            v_measure_matrix.loc[model_i, model_j] = v_meas
            v_measure_matrix.loc[model_j, model_i] = v_meas
            
            nmi_matrix.loc[model_i, model_j] = nmi
            nmi_matrix.loc[model_j, model_i] = nmi

    # Compute model-vs-GT scores.
    ari_scores = {}
    v_measure_scores = {}
    nmi_scores = {}
    
    gt_labels = results['GT_Preds'].astype(str).values
    
    for model in model_names:
        if model == 'GT_Preds': continue
        pred_labels = results[model].astype(str).values
        
        ari_scores[model] = adjusted_rand_score(gt_labels, pred_labels)
        v_measure_scores[model] = v_measure_score(gt_labels, pred_labels)
        nmi_scores[model] = normalized_mutual_info_score(gt_labels, pred_labels)
        
    return ari_matrix, ari_scores, v_measure_matrix, v_measure_scores, nmi_matrix, nmi_scores
