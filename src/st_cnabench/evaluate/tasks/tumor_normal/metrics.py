import pandas as pd
import numpy as np
import logging
from sklearn.metrics import accuracy_score, matthews_corrcoef, precision_score, recall_score, f1_score, adjusted_rand_score
from esda.moran import Moran
from esda.join_counts import Join_Counts
import libpysal
from functools import reduce

def _build_spatial_weights(coords, spatial_method, spatial_k, threshold):
    if spatial_method == "knn":
        weights = libpysal.weights.KNN.from_array(coords, k=spatial_k)
    elif spatial_method == "distance":
        weights = libpysal.weights.DistanceBand.from_array(
            coords,
            threshold=threshold,
            binary=True,
            silence_warnings=True,
        )
    else:
        raise ValueError(f"Unsupported spatial_method: {spatial_method}")

    weights.transform = "B"
    return weights


def run_metrics(results, ds_spatial_coord, spatial_method="distance", spatial_k=11, threshold=2.1):
    """
    Accuracy: the precision of tumor/normal classification.
    MCC: Matthews correlation coefficient for tumor/normal classification.
    Spatial Coherence Z: Join Count statistic Z-score for tumor spots.
    Call Rate: Proportion of spots with defined tumor/normal calls.
    """
    summary_list = []
    
    if 'GT' not in results:
        logging.error("GT data missing!")
        return None
    
    # Intersect
    logging.info("Building common index for alignment...")
    coords_df = pd.read_csv(ds_spatial_coord, index_col=0)
    index_list = [coords_df.index]
    for name, df in results.items():
        idx = df.set_index(df.columns[0]).index
        index_list.append(idx)
    common_index = reduce(lambda x, y: x.intersection(y), index_list)
    
    if len(common_index) == 0:
        logging.error("No common barcodes found!")
        return None
    
    logging.info(f"Aligned dataset size: {len(common_index)} spots.")

    # Subset
    aligned_coords = coords_df.loc[common_index, ['array_row', 'array_col']].values
    gt_df = results['GT'].set_index(results['GT'].columns[0])
    gt_labels = gt_df.loc[common_index].iloc[:, 0].astype(float).astype(int).values

    for model_name, res_df in results.items():
        if model_name == 'GT': continue
        model_series = res_df.set_index(res_df.columns[0]).iloc[:, 0]
        y_pred_series = model_series.loc[common_index].fillna(-1)
        y_pred = y_pred_series.astype(float).astype(int).values

        valid_mask = (y_pred != -1)
        call_rate = np.mean(valid_mask)

        y_valid = y_pred[valid_mask]
        gt_valid = gt_labels[valid_mask]
        coords_valid = aligned_coords[valid_mask]
        if np.any(valid_mask):
            acc = accuracy_score(gt_valid, y_valid)
            mcc = matthews_corrcoef(gt_valid, y_valid)
            precision = precision_score(gt_valid, y_valid)
            recall = recall_score(gt_valid, y_valid)
            f1 = f1_score(gt_valid, y_valid)
            ari = adjusted_rand_score(gt_valid, y_valid)
        else:
            acc, mcc, precision, recall, f1, ari = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

        w_local = _build_spatial_weights(coords_valid, spatial_method, spatial_k, threshold)
        unique_classes = np.unique(y_valid)
        if len(unique_classes) < 2:
            jc_z = np.nan
        else:
            jc = Join_Counts(y_valid, w_local)
            sim_std = np.std(jc.sim_bb)
            jc_z = (jc.bb - jc.mean_bb) / sim_std if sim_std > 0 else 0.0

        summary_list.append({
            'Model': model_name,
            'Accuracy': acc,
            'MCC': mcc,
            'Precision': precision,
            'Recall': recall,
            'F1_Score': f1,
            'ARI': ari,
            'Spatial_Coherence_Z': jc_z,
            'Call_Rate': call_rate
        })

    return pd.DataFrame(summary_list)
