import logging
import pandas as pd
import numpy as np
from ...utils.binning import map_to_bin, intersect_bins
from ...utils.genomics import generate_ref_bins

def run_mapping(results, ref_genome='hg38', bin_size=100000, fill_value=np.nan):
    ref_bins = generate_ref_bins(ref_genome=ref_genome, bin_size=bin_size)
    ref_has_chr = str(ref_bins['Chromosome'].iloc[0]).startswith('chr')
    aligned_results = {}

    for profile_key, models_dict in results.items():
        if not models_dict:
            continue
            
        logging.info(f"--- Mapping and Intersecting Profile: {profile_key} ---")
        all_binned_series = {}
        
        for name, profile_df in models_dict.items():
            profile_df = profile_df.copy()
            profile_df['Chromosome'] = profile_df['Chromosome'].astype(str)
            if ref_has_chr:
                profile_df['Chromosome'] = profile_df['Chromosome'].apply(lambda x: x if x.lower().startswith('chr') else f'chr{x}')
            else:
                profile_df['Chromosome'] = profile_df['Chromosome'].str.replace('chr', '', case=False)
            
            # CN track
            if profile_key == 'CN_Profile':
                score_col = 'GT_Score' if name == 'GT' else 'CN_Score'
                if score_col in profile_df.columns:
                    temp_df = profile_df[['Chromosome', 'Start', 'End', score_col]].rename(columns={score_col: 'Value'})
                    
                    # Keep GT column names clean and avoid GT_GT_* artifacts.
                    out_col = 'GT_Score' if name == 'GT' else f"{name}_CN_Score"
                    all_binned_series[out_col] = map_to_bin(temp_df, ref_bins, method='mean', fill_value=fill_value)
                
                if name == 'GT' and 'GT_Event' in profile_df.columns:
                    temp_df = profile_df[['Chromosome', 'Start', 'End', 'GT_Event']].rename(columns={'GT_Event': 'Value'})
                    all_binned_series["GT_Event"] = map_to_bin(temp_df, ref_bins, method='mode', fill_value=fill_value)
                    # Map focal-event flags to bins when they are available.
                    for tag in ['Is_Focal_Type1', 'Is_Focal_Type2']:
                        if tag in profile_df.columns:
                            temp_df_tag = profile_df[['Chromosome', 'Start', 'End', tag]].rename(columns={tag: 'Value'})
                            # Mark a bin as focal if any overlap is focal.
                            all_binned_series[tag] = map_to_bin(temp_df_tag, ref_bins, method='max', fill_value=0)
            # LOH track
            elif profile_key == 'LOH_Profile':
                target_col = 'GT_LOH_Status' if name == 'GT' else 'LOH_Status'
                if target_col in profile_df.columns:
                    temp_df = profile_df[['Chromosome', 'Start', 'End', target_col]].rename(columns={target_col: 'Value'})
                    
                    # Keep GT column names clean and avoid GT_GT_* artifacts.
                    out_col = 'GT_LOH_Status' if name == 'GT' else f"{name}_LOH_Status"
                    all_binned_series[out_col] = map_to_bin(temp_df, ref_bins, method='mean', fill_value=fill_value)

        # Intersect only the bins shared across all features in this track.
        if len(all_binned_series) > 1:
            logging.info(f"Performing intersection for features in {profile_key}")
            aligned_series_dict = intersect_bins(all_binned_series, fill_value=fill_value)
            
            merged_df = pd.DataFrame(aligned_series_dict)
            merged_df.index.name = 'BinID'
            merged_df = merged_df.reset_index()
            
            ref_meta = ref_bins[['Chromosome', 'Start', 'End', 'BinID']].drop_duplicates()
            merged_df = pd.merge(ref_meta, merged_df, on='BinID', how='inner')
            
            aligned_results[profile_key] = merged_df
        else:
            aligned_results[profile_key] = pd.DataFrame()
    
    return aligned_results
