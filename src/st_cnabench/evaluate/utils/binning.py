import pandas as pd
import pyranges as pr
import logging
import numpy as np
def map_to_bin(df_std, ref_bins, method='mean', fill_value=np.nan):
    """
    Map standard genomic data to reference bins using specified method.
    """
    # Normalize chromosome column naming before overlap mapping.
    if 'Chr' in df_std.columns:
        df_std = df_std.rename(columns={'Chr': 'Chromosome'})
    if 'Chr' in ref_bins.columns:
        ref_bins = ref_bins.rename(columns={'Chr': 'Chromosome'})

    gr_data = pr.PyRanges(df_std)
    gr_bins = pr.PyRanges(ref_bins)

    joined = gr_data.join(gr_bins, suffix="_bin").df
    
    if joined.empty:
        return pd.Series(fill_value, index=ref_bins['BinID'])
    
    joined['overlap_len'] = (
        np.minimum(joined['End'], joined['End_bin']) - 
        np.maximum(joined['Start'], joined['Start_bin'])
    )
    
    if method == 'mean':
        bin_values = joined.groupby('BinID')['Value'].mean()
    elif method == 'mode':
        idx = joined.groupby('BinID')['overlap_len'].idxmax()
        bin_values = joined.loc[idx].set_index('BinID')['Value']
    elif method == 'max':  
        bin_values = joined.groupby('BinID')['Value'].max()
    return bin_values.reindex(ref_bins['BinID'], fill_value=fill_value)

def intersect_bins(bin_profiles, fill_value=np.nan):
    """
    Intersect bins across multiple profiles, retaining only common bins.
    """
    common_idx = None
    
    for series in bin_profiles.values():
        # Handle both NaN-based and explicit fill-value masking.
        if pd.isna(fill_value):
            valid_idx = series[series.notna()].index
        else:
            valid_idx = series[series != fill_value].index
            
        if common_idx is None:
            common_idx = valid_idx
        else:
            common_idx = common_idx.intersection(valid_idx)
    
    if common_idx is None or common_idx.empty:
        raise ValueError("No common bins found across profiles.")

    logging.info(f"Bin alignment finished: {len(common_idx)} common bins retained across {len(bin_profiles)} features.")

    return {name: series.loc[common_idx] for name, series in bin_profiles.items()}
