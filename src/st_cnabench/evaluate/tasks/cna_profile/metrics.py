import logging
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from ...utils.io import save_tsv
from ...utils.metrics import calc_auc, calc_trunc_auc_by_sens, calc_max_f1_sweep, evaluate_loh

def evaluate_subset_cna(gt_event, gt_score, pred_score, opt_t_gain, opt_t_loss, prefix):
    res = {}
    
    # Gain evaluation
    is_true_gain = (gt_event == 1)
    if len(np.unique(is_true_gain)) > 1:
        res[f'{prefix}_Gain_AUC'] = calc_auc(is_true_gain.astype(int), pred_score)
        
        # Reuse the globally optimized threshold for subset sensitivity.
        is_pred_gain = (pred_score > opt_t_gain)
        tp_gains = np.sum(is_true_gain & is_pred_gain)
        res[f'{prefix}_Gain_Sensitivity'] = tp_gains / np.sum(is_true_gain) if np.sum(is_true_gain) > 0 else np.nan
    else:
        res[f'{prefix}_Gain_AUC'] = np.nan
        res[f'{prefix}_Gain_Sensitivity'] = np.nan

    # Loss evaluation
    is_true_loss = (gt_event == -1)
    if len(np.unique(is_true_loss)) > 1:
        res[f'{prefix}_Loss_AUC'] = calc_auc(is_true_loss.astype(int), -pred_score)
        
        is_pred_loss = (pred_score < opt_t_loss)
        tp_losses = np.sum(is_true_loss & is_pred_loss)
        res[f'{prefix}_Loss_Sensitivity'] = tp_losses / np.sum(is_true_loss) if np.sum(is_true_loss) > 0 else np.nan
    else:
        res[f'{prefix}_Loss_AUC'] = np.nan
        res[f'{prefix}_Loss_Sensitivity'] = np.nan
    return res

def run_metrics(aligned_results, save_dir=None):
    rows = []
    
    # Evaluate the TCN track for all models.
    if 'CN_Profile' in aligned_results and not aligned_results['CN_Profile'].empty:
        cn_df = aligned_results['CN_Profile']
        gt_event = cn_df['GT_Event'].values
        gt_score = cn_df['GT_Score'].values
        
        models = [col.replace('_CN_Score', '') for col in cn_df.columns if col.endswith('_CN_Score') and col != 'GT']
        
        for model_name in models:
            pred_score = cn_df[f"{model_name}_CN_Score"].values
            base_info = {'Method': model_name, 'Level': 'CN_Profile'}
            
            rows.append({**base_info, 'Metric': 'Gain_AUC', 'Value': calc_auc(gt_event == 1, pred_score)})
            rows.append({**base_info, 'Metric': 'Loss_AUC', 'Value': calc_auc(gt_event == -1, -pred_score)})
            #rows.append({**base_info, 'Metric': 'Gain_Trunc_AUC_Sens', 'Value': calc_trunc_auc_by_sens(gt_event == 1, pred_score, 0.0)})
            #rows.append({**base_info, 'Metric': 'Loss_Trunc_AUC_Sens', 'Value': calc_trunc_auc_by_sens(gt_event == -1, -pred_score, 0.0)})
            
            sweep_res = calc_max_f1_sweep(gt_event, pred_score)
            for m_name, m_val in sweep_res.items():
                rows.append({**base_info, 'Metric': m_name, 'Value': m_val})
                
            opt_t_gain = sweep_res['Opt_Gain_Threshold']
            opt_t_loss = sweep_res['Opt_Loss_Threshold']

            pcc, _ = pearsonr(pred_score, gt_score)
            spearman_cc, _ = spearmanr(pred_score, gt_score)
            rows.append({**base_info, 'Metric': 'PCC', 'Value': pcc})
            rows.append({**base_info, 'Metric': 'Spearman_CC', 'Value': spearman_cc})

            # Evaluate focal type 1 events (< 3 Mb).
            if 'Is_Focal_Type1' in cn_df.columns:
                # Use focal type 1 bins plus all neutral bins as the evaluation set.
                mask_t1 = (cn_df['Is_Focal_Type1'] == 1) | (gt_event == 0)
                res_t1 = evaluate_subset_cna(
                    gt_event[mask_t1], gt_score[mask_t1], pred_score[mask_t1],
                    opt_t_gain, opt_t_loss, prefix="Focal_T1"
                )
                for k, v in res_t1.items():
                    rows.append({**base_info, 'Metric': k, 'Value': v})

            # Evaluate focal type 2 events (< 25% of a chromosome arm).
            if 'Is_Focal_Type2' in cn_df.columns:
                mask_t2 = (cn_df['Is_Focal_Type2'] == 1) | (gt_event == 0)
                res_t2 = evaluate_subset_cna(
                    gt_event[mask_t2], gt_score[mask_t2], pred_score[mask_t2],
                    opt_t_gain, opt_t_loss, prefix="Focal_T2"
                )
                for k, v in res_t2.items():
                    rows.append({**base_info, 'Metric': k, 'Value': v})
    # Evaluate the LOH track for allele-aware models only.
    if 'LOH_Profile' in aligned_results and not aligned_results['LOH_Profile'].empty:
        loh_df = aligned_results['LOH_Profile']
        gt_loh = loh_df['GT_LOH_Status'].values
        
        # Due to the bin mapping, re-binary the GT LOH status at 0.5 threshold to calculate the AUC
        gt_loh_binary = (gt_loh >= 0.5).astype(int)
        models = [
            col.replace('_LOH_Status', '')
            for col in loh_df.columns
            if col.endswith('_LOH_Status') and col != 'GT_LOH_Status'
        ]
        
        for model_name in models:
            pred_loh = loh_df[f"{model_name}_LOH_Status"].values
            base_info = {'Method': model_name, 'Level': 'LOH_Profile'}
            
            loh_metrics = evaluate_loh(gt_loh_binary, pred_loh)
            for m_name, m_val in loh_metrics.items():
                rows.append({**base_info, 'Metric': m_name, 'Value': m_val})

    df_score = pd.DataFrame(rows)
    df_score['Value'] = pd.to_numeric(df_score['Value'], errors='coerce')

    if save_dir:
        save_path = f"{save_dir}/metrics_long_format.tsv"
        save_tsv(df_score, save_path)
        logging.info(f"Evaluation completed. Results saved to {save_path}")

    return df_score
