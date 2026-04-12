# src/utils/metrics.py
import numpy as np
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr, spearmanr
import libpysal
import itertools
from sklearn.metrics import roc_auc_score, f1_score, average_precision_score, roc_curve, auc, mean_absolute_error

def calc_auc(gt_binary, pred_scores):
    """Compute the standard ROC AUC."""
    if len(np.unique(gt_binary)) <= 1:
        return np.nan
    return roc_auc_score(gt_binary, pred_scores)

def calc_trunc_auc_by_sens(gt_binary, pred_scores, bio_threshold=0.0):
    """
    Compute a sensitivity-truncated partial AUC.
    This follows the pROC behavior:
    auc(roc, partial.auc=c(0, sens_threshold), partial.auc.focus="sens")
    """
    if len(np.unique(gt_binary)) <= 1:
        return np.nan
        
    # Compute sensitivity at the biological threshold.
    pred_binary = (pred_scores > bio_threshold).astype(int)
    true_positives = np.sum((pred_binary == 1) & (gt_binary == 1))
    actual_positives = np.sum(gt_binary == 1)
    
    sens_threshold = true_positives / actual_positives if actual_positives > 0 else 0.0
    
    if sens_threshold == 0.0:
        return 0.0 
        
    # Build the full ROC curve and truncate it.
    fpr, tpr, _ = roc_curve(gt_binary, pred_scores)
    valid_idx = np.where(tpr <= sens_threshold)[0]
    p_fpr = fpr[valid_idx].tolist()
    p_tpr = tpr[valid_idx].tolist()
    
    # Add the boundary point with stepwise interpolation when needed.
    if len(p_tpr) > 0 and p_tpr[-1] < sens_threshold:
        idx_over = np.where(tpr > sens_threshold)[0][0]
        p_fpr.append(fpr[idx_over])
        p_tpr.append(sens_threshold)
        
    # Calculate the partial AUC on the truncated curve.
    return auc(p_fpr, p_tpr)

def calc_max_f1_sweep(gt_events, pred_scores):
    """
    Perform a data-driven max-F1 sweep and report the best thresholds together
    with sensitivity and precision.
    """
    unique_scores = np.sort(np.unique(pred_scores))
    
    # Build candidate loss thresholds.
    loss_scores = unique_scores[unique_scores < 0]
    if len(loss_scores) > 200:
        idx = np.round(np.linspace(0, len(loss_scores) - 1, 200)).astype(int)
        loss_scores = loss_scores[idx]
    elif len(loss_scores) == 0:
        loss_scores = np.array([-1.0])
        
    # Build candidate gain thresholds.
    gain_scores = unique_scores[unique_scores > 0]
    if len(gain_scores) > 200:
        idx = np.round(np.linspace(0, len(gain_scores) - 1, 200)).astype(int)
        gain_scores = gain_scores[idx]
    elif len(gain_scores) == 0:
        gain_scores = np.array([1.0])

    loss_scores = np.append(loss_scores, 0.0)
    gain_scores = np.insert(gain_scores, 0, 0.0)

    best_f1 = -1
    best_t_loss = 0.0
    best_t_gain = 0.0
    
    # Grid-search the best loss/gain threshold pair.
    for t_loss, t_gain in itertools.product(loss_scores, gain_scores):
        pred_discrete = np.zeros_like(pred_scores)
        pred_discrete[pred_scores > t_gain] = 1
        pred_discrete[pred_scores < t_loss] = -1
        
        current_f1 = f1_score(gt_events, pred_discrete, average='macro')
        if current_f1 > best_f1:
            best_f1 = current_f1
            best_t_loss = t_loss
            best_t_gain = t_gain
            
    # Compute sensitivity and precision at the optimal thresholds.
    is_true_gain = (gt_events == 1)
    is_pred_gain = (pred_scores > best_t_gain)
    tp_gains = np.sum(is_true_gain & is_pred_gain)
    sens_gains = tp_gains / np.sum(is_true_gain) if np.sum(is_true_gain) > 0 else np.nan
    prec_gains = tp_gains / np.sum(is_pred_gain) if np.sum(is_pred_gain) > 0 else np.nan

    is_true_loss = (gt_events == -1)
    is_pred_loss = (pred_scores < best_t_loss)
    tp_losses = np.sum(is_true_loss & is_pred_loss)
    sens_losses = tp_losses / np.sum(is_true_loss) if np.sum(is_true_loss) > 0 else np.nan
    prec_losses = tp_losses / np.sum(is_pred_loss) if np.sum(is_pred_loss) > 0 else np.nan
    
    return {
        'Max_Macro_F1': best_f1,
        'Opt_Gain_Threshold': best_t_gain,
        'Gain_Sensitivity': sens_gains,
        'Gain_Precision': prec_gains,
        'Opt_Loss_Threshold': best_t_loss,
        'Loss_Sensitivity': sens_losses,
        'Loss_Precision': prec_losses
    }

def evaluate_loh(gt_loh, pred_loh_freq):
    """Evaluate CN-LOH predictions."""
    gt_loh = np.asarray(gt_loh, dtype=float)
    pred_loh_freq = np.asarray(pred_loh_freq, dtype=float)

    valid_mask = (~np.isnan(gt_loh)) & (~np.isnan(pred_loh_freq))

    if valid_mask.sum() == 0:
        return {
            'LOH_ROC_AUC': 0.5,
            'LOH_PR_AUC': 0.0,
            'LOH_Max_F1': 0.0,
            'LOH_PCC': np.nan,
            'LOH_Spearman_CC': np.nan
        }

    gt_valid = gt_loh[valid_mask]
    pred_valid = pred_loh_freq[valid_mask]

    roc_auc = calc_auc(gt_valid, pred_valid)
    pr_auc = average_precision_score(gt_valid, pred_valid) if len(np.unique(gt_valid)) > 1 else np.nan

    best_f1 = 0
    for threshold in np.linspace(0.05, 0.95, 20):
        pred_binary = (pred_valid >= threshold).astype(int)
        f1 = f1_score(gt_valid, pred_binary, zero_division=0)
        if f1 > best_f1:
            best_f1 = f1

    if len(gt_valid) < 2 or np.std(gt_valid) == 0 or np.std(pred_valid) == 0:
        loh_pcc = np.nan
        loh_scc = np.nan
    else:
        loh_pcc = pearsonr(gt_valid, pred_valid)[0]
        loh_scc = spearmanr(gt_valid, pred_valid)[0]

    return {
        'LOH_ROC_AUC': roc_auc,
        'LOH_PR_AUC': pr_auc,
        'LOH_Max_F1': best_f1,
        'LOH_PCC': loh_pcc,
        'LOH_Spearman_CC': loh_scc
    }


def calculate_pcc(pred, gt):
    # check std
    if np.std(pred) == 0 or np.std(gt) == 0:
        return 0.0
    return pearsonr(pred, gt)[0]

def calculate_cos_sim(pred, gt):
    # Cosine Similarity = 1 - Cosine Distance
    return 1 - cosine(pred, gt)

def calculate_mae(pred, gt):
    return mean_absolute_error(gt, pred)

def calculate_f1_macro(pred, gt):
    return f1_score(gt, pred, average='macro')


def calculate_spatial_coherence(coords, labels, method='distance', k=11, threshold=2.1, calculate_z_score=False, n_permutations=100):
    """
    Calculate spatial coherence of given labels based on spatial coordinates using Join Count statistic.
    
    Args:
        coords: (N, 2) array, spatial coordinates
        labels: (N,) array, corresponding cluster labels
        method: 'knn' (for Slide-seq/Stereo-seq) or 'distance' (for Visium)
        k: number of nearest neighbors (used if method='knn'). Default 11 for slide rna seq.
        threshold: distance threshold to define neighbors (used if method='distance')
        calculate_z_score: whether to calculate Z-score using permutation test (increases computation time)
        n_permutations: number of permutations for Z-score calculation
        
    Returns:
        score: spatial coherence score (0.0 to 1.0)
        z_score: (optional) Z-score if calculate_z_score=True, otherwise None
    """
    # 1. Construct spatial weights (Binary)
    if method == 'knn':
        # Use KNN graphs for irregular spatial densities.
        w = libpysal.weights.KNN.from_array(coords, k=k)
    elif method == 'distance':
        # Use absolute distance neighborhoods for Visium-style lattices.
        w = libpysal.weights.DistanceBand.from_array(
            coords, 
            threshold=threshold, 
            binary=True, 
            silence_warnings=True
        )
    else:
        raise ValueError("Invalid method. Choose 'knn' or 'distance'.")
    
    # Helper function: calculate Concordance once given labels and weight object
    def _get_concordance(lab, weight_obj):
        
        match_count = 0
        total_neighbors = 0
        
        for idx, neighbors in weight_obj.neighbors.items():
            if not neighbors: continue # Skip isolated points
            
            # Get the label of the current point
            curr_label = lab[idx]
            # Get the labels of all neighbors
            neighbor_labels = lab[neighbors]
            
            # Count how many neighbors have the same label as the current point
            match_count += np.sum(neighbor_labels == curr_label)
            total_neighbors += len(neighbors)
            
        if total_neighbors == 0:
            return 0.0
        
        return match_count / total_neighbors

    # 2. Calculate observed spatial coherence
    observed_score = _get_concordance(labels, w)
    
    # 3. If Z-score is not needed, return directly
    if not calculate_z_score:
        return observed_score

    # 4. If Z-score is needed, perform permutation test
    perm_scores = []
    temp_labels = labels.copy()
    
    for _ in range(n_permutations):
        np.random.shuffle(temp_labels) # In-place shuffle
        perm_scores.append(_get_concordance(temp_labels, w))
    
    perm_scores = np.array(perm_scores)
    mean_perm = np.mean(perm_scores)
    std_perm = np.std(perm_scores)
    
    # Prevent division by zero (extremely rare case)
    if std_perm == 0:
        z_score = 0.0
    else:
        z_score = (observed_score - mean_perm) / std_perm
        
    return observed_score, z_score
