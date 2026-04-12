import pandas as pd
from .constants import CHR_SIZES_HG38


def generate_ref_bins(ref_genome='hg38', bin_size=100000):
    '''
    Generate reference genomic bins for given genome and bin size.
    '''
    if ref_genome == 'hg38':
        chrom_sizes=CHR_SIZES_HG38
    bin_list = []
    for chrom, length in chrom_sizes.items():
        for start in range(0, length, bin_size):
            end = min(start + bin_size, length)
            bin_list.append({
                'Chromosome': chrom,
                'Start': start,
                'End': end,
                'BinID': f"chr{chrom}_{start // bin_size}"
            })
    return pd.DataFrame(bin_list)