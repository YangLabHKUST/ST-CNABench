import logging
import pandas as pd
import numpy as np
from ...utils.binning import map_to_bin, intersect_bins
from ...utils.genomics import generate_ref_bins

def run_clonal_mapping(clone_cnv_profiles, ref_genome='hg38', bin_size=100000, fill_value=np.nan):
    """
    Map all clone profiles to shared reference bins and keep only bins present
    across every clone/model profile.
    """
    logging.info(f"--- Mapping Clones to {bin_size}bp Bins & Intersecting (Keeping only common bins) ---")
    
    ref_bins = generate_ref_bins(ref_genome=ref_genome, bin_size=bin_size)
    ref_has_chr = str(ref_bins['Chromosome'].iloc[0]).startswith('chr')
    
    all_binned_series = {}

    for model_pred_name, profile_df in clone_cnv_profiles.items():
        if profile_df is None or profile_df.empty:
            continue
            
        logging.info(f"Processing mapping for model: {model_pred_name}")
        df = profile_df.copy()
        
        # Normalize chromosome and score column names.
        chr_col = 'Chromosome' if 'Chromosome' in df.columns else 'Chr'
        score_col = 'CN_Score' if 'CN_Score' in df.columns else 'TCN_Score'
        
        df[chr_col] = df[chr_col].astype(str)
        if ref_has_chr:
            df[chr_col] = df[chr_col].apply(lambda x: x if x.lower().startswith('chr') else f'chr{x}')
        else:
            df[chr_col] = df[chr_col].str.replace('chr', '', case=False)

        clones = df['Clone_ID'].dropna().unique()

        for clone in clones:
            temp_df = df[df['Clone_ID'] == clone].copy()
            input_df = temp_df[[chr_col, 'Start', 'End', score_col]].rename(columns={chr_col: 'Chromosome', score_col: 'Value'})
            
            # Map the clone profile to the shared reference bins.
            binned_values = map_to_bin(input_df, ref_bins, method='mean', fill_value=fill_value)

            # Preserve clones with empty coverage as neutral profiles.
            if binned_values.isna().all():
                logging.warning(
                    f"Clone {clone} in {model_pred_name} has empty CNV profile; filling with 0 to keep normal clone."
                )
                binned_values = pd.Series(0.0, index=ref_bins['BinID'])
            
            # Use a stable separator so the wide table can be split back later.
            series_key = f"{model_pred_name}___{clone}"
            all_binned_series[series_key] = binned_values

    aligned_results = {}

    if len(all_binned_series) > 1:
        logging.info("Performing global intersection across all clones and models...")
        aligned_series_dict = intersect_bins(all_binned_series, fill_value=fill_value)
        
        merged_df = pd.DataFrame(aligned_series_dict)
        merged_df.index.name = 'BinID'
        merged_df = merged_df.reset_index()
        
        # Restore chromosome/start/end metadata after bin intersection.
        ref_meta = ref_bins[['Chromosome', 'Start', 'End', 'BinID']].drop_duplicates()
        merged_df = pd.merge(ref_meta, merged_df, on='BinID', how='inner')
        
        # Split the wide intersected table back into one long table per model.
        for model_pred_name in clone_cnv_profiles.keys():
            model_cols = [c for c in merged_df.columns if c.startswith(f"{model_pred_name}___")]
            if not model_cols:
                continue
                
            meta_cols = ['BinID', 'Chromosome', 'Start', 'End']
            model_df_wide = merged_df[meta_cols + model_cols]
            
            model_df_long = model_df_wide.melt(
                id_vars=meta_cols,
                value_vars=model_cols,
                var_name='Clone_ID',
                value_name='CN_Score'
            )
            
            # Remove the temporary model prefix from Clone_ID.
            model_df_long['Clone_ID'] = model_df_long['Clone_ID'].str.replace(f"{model_pred_name}___", "", regex=False)
            
            aligned_results[model_pred_name] = model_df_long

    else:
        logging.warning("Not enough profiles to perform intersection.")
        for model_pred_name in clone_cnv_profiles.keys():
             aligned_results[model_pred_name] = pd.DataFrame()

    return aligned_results
