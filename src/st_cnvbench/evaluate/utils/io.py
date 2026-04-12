import pandas as pd
from pathlib import Path
import scanpy as sc
import anndata as ad
import numpy as np
import scipy.stats
import pyreadr
import logging
import io
# General I/O utilities
def load_txt(txt_path: Path, gzip: bool = False):
    df = pd.read_csv(txt_path, sep=r'\s+', index_col=0, compression='gzip' if gzip else None)
    return df
def load_tsv(tsv_path: Path, gzip: bool = False, header = False):
    df = pd.read_csv(tsv_path, sep='\t', index_col=0, compression='gzip' if gzip else None, header=0 if header else None)
    return df
def load_csv(tsv_path: Path, gzip: bool = False, header = False):
    df = pd.read_csv(tsv_path, index_col=0, compression='gzip' if gzip else None, header=0 if header else None)
    return df
'''
def load_gatk_seg(file_path: str):
    """Load a GATK seg file and skip header lines starting with '@'."""
    # comment='@' skips GATK dictionary headers and read-group metadata.
    df = pd.read_csv(file_path, sep='\t', comment='@')
    return df
'''
def load_h5ad(h5ad_path: Path):
    adata = sc.read_h5ad(h5ad_path)
    return adata
def load_gene_annot(gene_annot_path, index_col = 'Gene_ID'):
    """
    Helper to load gene order file.
    """
    df_ann = pd.read_csv(
        gene_annot_path, 
        sep='\t', 
        header=None,
        usecols=[0, 1, 2, 3, 4],
        names=['Chromosome', 'Start', 'End','Gene_name', 'Gene_ID']
    )
    df_ann = df_ann.set_index(index_col)
    return df_ann

def load_RData(file_path):
    result = pyreadr.read_r(str(file_path))
    
    keys = list(result.keys())
    if not keys:
        return pd.DataFrame()
        
    return result[keys[0]]

def save_tsv(df, out_path):
    df.to_csv(out_path, sep='\t', index=False)
    
def save_adata(adata, out_fn):
    return (adata.write_h5ad(out_fn, compression = "gzip"))

# GT I/O
def load_facets(gt_path):
    # Load VCF-like file, skipping header lines
    with open(gt_path, 'r') as f:
        valid_lines = [line for line in f if not line.startswith('##')]
    df = pd.read_csv(io.StringIO(''.join(valid_lines)), sep='\t')
    df = df.rename(columns={'#CHROM': 'Chromosome'})
    
    # Parse INFO column
    def get_info_val(info_str):
        kv = {x.split('=')[0]: x.split('=')[1] for x in info_str.split(';') if '=' in x}
        return kv

    infos = df['INFO'].apply(get_info_val)

    chrom = df['Chromosome'].astype(str).str.replace('chr', '', case=False)
    start = df['POS'].astype(int)
    end = infos.apply(lambda x: int(x.get('END', 0)))
    
    # Extract total copy number (TCN) and minor-allele copy number (LCN).
    # Missing '.' values are treated as diploid defaults (TCN=2, LCN=1).
    tcn = infos.apply(lambda x: int(x.get('TCN_EM', 2)) if x.get('TCN_EM', '.') != '.' else 2).values
    lcn = infos.apply(lambda x: int(x.get('LCN_EM', 1)) if x.get('LCN_EM', '.') != '.' else 1).values

    seg_ids = chrom + "_" + start.astype(str) + "_" + end.astype(str)

    # 1. Continuous GT score: convert TCN to a zero-centered log2 ratio.
    # Add a 0.1 pseudocount to avoid log2(0) at homozygous deletions.
    gt_score = np.log2(np.maximum(tcn, 0.1) / 2.0)
    
    # 2. Discrete GT event labels for Macro-F1 and AUC evaluation.
    # >2 is gain (1), <2 is loss (-1), and =2 is neutral (0).
    gt_event = np.where(tcn > 2, 1, 
               np.where(tcn < 2, -1, 0))
               
    # 3. Heterozygous loss status for the CN-LOH evaluation track.
    # CN-LOH is defined as unchanged total copy number (tcn==2) with lost
    # minor allele copy number (lcn==0).
    loh_status = np.where((tcn == 2) & (lcn == 0), 1, 0)

    # Assemble the final unified segment table.
    df_unified = pd.DataFrame({
        "Chromosome": chrom,
        "Start": start,
        "End": end,
        "ID": seg_ids,
        "GT_Score": gt_score,      # For PCC-based continuous concordance.
        "GT_Event": gt_event,      # For AUC and Max-F1 classification metrics.
        "LOH_Status": loh_status   # For dedicated CN-LOH evaluation.
    })

    return df_unified
# Model raw result I/O

def load_CalicoST_intcn(cnv_file, clone_label_file, tn_annot_path):
    df_cnv = load_tsv(cnv_file, header=True).reset_index()
    df_label = load_tsv(clone_label_file, header=True).reset_index()
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    
    # Remove normal cells
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    initial_count = len(df_label)
    df_label = df_label[~df_label['BARCODES'].isin(normal_barcodes)]
    tumor_cells_count = len(df_label)
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                 f"{tumor_cells_count} tumor cells remain for calculation.")
    
    clone_ids = sorted(list(set([c.split(' ')[0] for c in df_cnv.columns if 'clone' in c])))

    # Build parallel TCN and LOH tracks.
    clone_tcn_dict = {}
    clone_loh_dict = {}
    
    for cid in clone_ids:
        a_val = df_cnv[f"{cid} A"]
        b_val = df_cnv[f"{cid} B"]
        
        # Track 1: total copy number.
        clone_tcn_dict[cid] = a_val + b_val
        
        # Track 2: LOH status, where one allele is 0 and the other is 2.
        # Convert the boolean mask to 0/1 integers.
        clone_loh_dict[cid] = ((a_val == 0) & (b_val == 2)) | ((a_val == 2) & (b_val == 0)).astype(int)

    # Convert the track dictionaries to DataFrames.
    clone_tcn_df = pd.DataFrame(clone_tcn_dict, index=df_cnv.index)
    clone_loh_df = pd.DataFrame(clone_loh_dict, index=df_cnv.index)
    
    def _format_clone_id(x):
        s = str(x)
        return s if s.startswith('clone') else f"clone{s}"
    
    df_label['clone_key'] = df_label['clone_label'].apply(_format_clone_id)
    cell_clone_map = df_label['clone_key'].values  # shape: n_cells
    
    # Build cell-by-segment matrices.
    matrix_X = clone_tcn_df.T.loc[cell_clone_map].values
    matrix_LOH = clone_loh_df.T.loc[cell_clone_map].values
    
    var = df_cnv[['CHR', 'START', 'END']].copy()
    var.columns = ['Chromosome', 'Start', 'End']
    var.index = (var['Chromosome'].astype(str) + "_" + var['Start'].astype(str) + "_" + var['End'].astype(str))
    # Store TCN in X and LOH in an additional layer.
    adata = ad.AnnData(
        X=matrix_X,  # Main matrix stores the TCN track.
        obs=df_label.set_index('BARCODES'),
        var=var
    )
    
    # LOH uses the same Cells x Segments layout as X.
    adata.layers['LOH_Status'] = matrix_LOH
    
    return adata
'''
def load_CalicoST_cnevent(cnv_file, clone_label_file, tn_annot_path):
    df_cnv = load_tsv(cnv_file, header = True).reset_index()
    df_label = load_tsv(clone_label_file, header = True).reset_index()
    df_tn_annot = load_tsv(tn_annot_path, header = True).reset_index()
    # Remove normal cells
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'] == 'normal']['Barcode'])
    initial_count = len(df_label)
    df_label = df_label[~df_label['BARCODES'].isin(normal_barcodes)]
    tumor_cells_count = len(df_label)
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                 f"{tumor_cells_count} tumor cells remain for calculation.")
    
    # Generate CNV event matrix
    clone_ids = sorted(list(set([c.split(' ')[0] for c in df_cnv.columns if 'clone' in c])))
    
    # Map to events
    clone_event_dict = {}    
    for cid in clone_ids:
        a_val = df_cnv[f"{cid} A"]
        b_val = df_cnv[f"{cid} B"]
        total_cn = a_val + b_val
        min_val = np.minimum(a_val, b_val)
        
        # Pay attention to the mapping relationship
        events = np.ones_like(total_cn) 
        events[total_cn < 2] = 0        
        events[(total_cn == 2) & (min_val == 0)] = 2 
        events[total_cn > 2] = 3        
        
        clone_event_dict[cid] = events
    clone_event_df = pd.DataFrame(clone_event_dict, index=df_cnv.index)
    
    def _format_clone_id(x):
        s = str(x)
        return s if s.startswith('clone') else f"clone{s}"

    df_label['clone_key'] = df_label['clone_label'].apply(_format_clone_id)
    cell_clone_map = df_label['clone_key'].values  # shape: n_cells

    matrix_X = clone_event_df.T.loc[cell_clone_map].values
    var_df = df_cnv[['CHR', 'START', 'END']].copy()
    var_df.columns = ['Chromosome', 'Start', 'End']

    adata = ad.AnnData(
        X=matrix_X.astype(np.int8),
        obs=df_label.set_index('BARCODES'),
        var=var_df
    )
    
    return adata
'''
def load_InferCNV_expr(matrix_file, gene_annot_path, tn_annot_path):
    df_data = load_txt(matrix_file)
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    # Remove norm cells (Do not need actually)
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    initial_count = df_data.shape[1]
    df_data = df_data.drop(columns=normal_barcodes, errors='ignore')
    tumor_cells_count = df_data.shape[1]
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                    f"{tumor_cells_count} tumor cells remain for calculation.")
    
    obs = pd.DataFrame(index=df_data.columns)
    df_ann = load_gene_annot(gene_annot_path)
    common_genes = df_data.index.intersection(df_ann.index)
    df_data = df_data.loc[common_genes]
    var = df_ann.loc[common_genes]

    # Log2 Copy Ratio
    raw_values = np.log2(df_data.values).T
    adata = ad.AnnData(
        X=raw_values,
        obs=obs,
        var=var
    )
    return adata

def load_InferCNV_cnv(matrix_file, gene_annot_path, tn_annot_path):
    df_data = load_txt(matrix_file)
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    # Remove norm cells (Do not need actually)
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    initial_count = df_data.shape[1]
    df_data = df_data.drop(columns=normal_barcodes, errors='ignore')
    tumor_cells_count = df_data.shape[1]
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                    f"{tumor_cells_count} tumor cells remain for calculation.")
    
    obs = pd.DataFrame(index=df_data.columns)
    df_ann = load_gene_annot(gene_annot_path)
    common_genes = df_data.index.intersection(df_ann.index)
    df_data = df_data.loc[common_genes]
    var = df_ann.loc[common_genes]

    # Integer_CN
    raw_values = (df_data.values *2).T
    adata = ad.AnnData(
        X=raw_values.astype(np.int8),
        obs=obs,
        var=var
    )
    return adata

'''
#load CopyKAT raw results gene by cell matrix txt
def load_CopyKAT_genebycell(matrix_file, tn_annot_path):

    df_all = load_txt(matrix_file).reset_index()
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    all_metadata_cols = ["abspos", "chromosome_name", "start_position", "end_position", "ensembl_gene_id", "hgnc_symbol", "band"]
    existing_meta = [c for c in all_metadata_cols if c in df_all.columns]
    barcode_cols = [c for c in df_all.columns if c not in existing_meta]

    var_subset = df_all[existing_meta].copy()
    var_subset = var_subset.rename(columns={
        "chromosome_name": "Chromosome",
        "start_position": "Start",
        "end_position": "End",
        "ensembl_gene_id": "Gene_ID",
        "hgnc_symbol": "Gene_name"
    })
    if "Gene_ID" in var_subset.columns:
        var_subset.index = var_subset["Gene_ID"].values
    # Remove normal cells
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    tumor_barcodes = [bc for bc in barcode_cols if bc not in normal_barcodes]
    logging.info(f"Filtered {len(barcode_cols) - len(tumor_barcodes)} normal cells. "
                    f"{len(tumor_barcodes)} tumor cells remain for calculation.")

    X_data = df_all[tumor_barcodes].values
    raw_values = X_data.T

    adata = ad.AnnData(
        X=raw_values,
        obs=pd.DataFrame(index=tumor_barcodes),
        var=var_subset
    )
    return adata
'''
def load_CopyKAT_binbycell(matrix_file, tn_annot_path):

    df_all = load_txt(matrix_file).reset_index()
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()


    all_metadata_cols = ["chrom", "chrompos", "abspos"]
    existing_meta = [c for c in all_metadata_cols if c in df_all.columns]
    barcode_cols = [c for c in df_all.columns if c not in existing_meta]
    # rename CB e.g. ATAC.X -> ATAC-X
    rename_dict = {bc: bc.replace('.', '-') for bc in barcode_cols}
    df_all = df_all.rename(columns=rename_dict)

    barcode_cols = list(rename_dict.values())
    # Remove normal cells
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    tumor_barcodes = [bc for bc in barcode_cols if bc not in normal_barcodes]
    logging.info(f"Filtered {len(barcode_cols) - len(tumor_barcodes)} normal cells. "
                    f"{len(tumor_barcodes)} tumor cells remain for calculation.")

    # Bin: Chr Start and End (chrompos -> Bin End, Bin Start = previous Bin End / 0)
    var_df = df_all[existing_meta].copy()
    var_df['End'] = var_df['chrompos']
    var_df['Start'] = var_df.groupby('chrom')['End'].shift(1, fill_value=0).astype(int)
    var_df['BinID'] = (var_df['chrom'].astype(str) + "_" + 
                       var_df['Start'].astype(str) + "_" + 
                       var_df['End'].astype(str))
    var_df = var_df.rename(columns={
        "chrom": "Chromosome",
    })
    X_data = df_all[tumor_barcodes].values
    raw_values = X_data.T

    adata = ad.AnnData(
        X=raw_values.astype(np.float32),
        obs=pd.DataFrame(index=tumor_barcodes),
        var=var_df.set_index('BinID')
    )
    return adata

def load_SCEVAN_expr(matrix_file, gene_ann_file, tn_annot_path):

    matrix = load_RData(matrix_file)
    annot = load_RData(gene_ann_file)
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    # Remove normal cells (Do not need actually)
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    initial_count = matrix.shape[1]
    matrix = matrix.drop(columns=normal_barcodes, errors='ignore')
    tumor_cells_count = matrix.shape[1]
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                    f"{tumor_cells_count} tumor cells remain for calculation.")
    
    annot = annot.rename(columns={
        "seqnames": "Chromosome",
        "start": "Start",
        "end": "End",
        "gene_id": "Gene_ID",
        "gene_name": "Gene_name"
    })
    annot['Chromosome'] = annot['Chromosome'].astype(int).astype(str)

    obs = pd.DataFrame(index=matrix.columns)
    var = annot
    var.index = (var['Chromosome'].astype(str) + "_" + var['Start'].astype(str) + "_" + var['End'].astype(str))
    raw_values = matrix.values.T
    adata = ad.AnnData(
        X=raw_values, 
        obs=obs,
        var=var
    )
    return adata


def load_SCEVAN_cnv(seg_file):
    seg = load_txt(seg_file).reset_index()
    SegID = seg['Chr'].astype(str) + "_" + seg['Pos'].astype(str) + "_" + seg['End'].astype(str)
    df_seg = pd.DataFrame({
        "Chromosome": seg['Chr'],
        "Start": seg['Pos'],
        "End": seg['End'],
        "ID": SegID,
        "CN_Score": seg["CN"],
        "LOH_Status": np.nan
    })
    return df_seg


def load_Numbat_expr(matrix_file, gene_annot_path, tn_annot_path):
    df_data = load_txt(matrix_file, gzip = True)
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    # Remove normal cells
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    # Cell by gene
    initial_count = df_data.shape[0]
    df_data = df_data.drop(index=normal_barcodes, errors='ignore')
    tumor_cells_count = df_data.shape[0]
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                    f"{tumor_cells_count} tumor cells remain for calculation.")
    
    obs = pd.DataFrame(index=df_data.index)
    df_ann = load_gene_annot(gene_annot_path, index_col='Gene_name')
    common_genes = df_data.columns.intersection(df_ann.index)
    if len(common_genes) == 0:
        logging.error(f"No common genes found between Numbat output and annotation!")
        raise ValueError("Gene alignment failed. Check Gene_ID format.")
    df_data = df_data[common_genes]
    var = df_ann.loc[common_genes]

    raw_values =  df_data.values
    adata = ad.AnnData(
        X=raw_values,
        obs=obs,
        var=var
    )
    # set index 
    adata.var_names = adata.var['Gene_ID'].astype(str).values
    return adata

'''
def load_Numbat_cnv(matrix_file, consensus_file, gene_annot_path, tn_annot_path):

    # See https://github.com/kharchenkolab/numbat/issues/100
    # Track 1: map TCN states to (-1, 0, 1).
    def _map_tcn_event(state):
        mapping = {
            'del': -1,  # Physical deletion.
            'neu': 0,   # Diploid neutral state.
            'loh': 0,   # CN-LOH keeps total copy number neutral here.
            'amp': 1,   # Amplification.
            'bamp': 1,  # Balanced amplification.
        }
        return mapping.get(state, 0)  # Default to neutral if unknown.

    # Track 2: map LOH states to (1, 0).
    def _map_loh_event(state):
        mapping = {
            'del': 0,   
            'neu': 0,   
            'loh': 1,   # Dedicated CN-LOH detection track.
            'amp': 0,   
            'bamp': 0,
        }
        return mapping.get(state, 0)

    # Load the required inputs.
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    df_consensus = load_tsv(consensus_file, header=True).reset_index()
    
    # Build a stable segment identifier.
    df_consensus['seg_id'] = (
        df_consensus['CHROM'].astype(str) + "_" + 
        df_consensus['seg_start'].astype(str) + "_" + 
        df_consensus['seg_end'].astype(str)
    )
    
    df_var = df_consensus[['seg_id', 'CHROM', 'seg_start', 'seg_end']].drop_duplicates().set_index('seg_id')
    df_var = df_var.rename(columns={"CHROM": "Chromosome", "seg_start": "Start", "seg_end": "End"})
    
    df_data = load_tsv(matrix_file, header=True).reset_index()
    df_data['ID'] = (
        df_data['CHROM'].astype(str) + "_" + 
        df_data['seg_start'].astype(str) + "_" + 
        df_data['seg_end'].astype(str)
    )

    # Populate both TCN and LOH tracks.
    df_data['tcn_event'] = df_data['cnv_state'].apply(_map_tcn_event)
    df_data['loh_event'] = df_data['cnv_state'].apply(_map_loh_event)

    # Cell by seg (Pivot)
    matrix_tcn = df_data.pivot(index='cell', columns='ID', values='tcn_event')
    matrix_loh = df_data.pivot(index='cell', columns='ID', values='loh_event')
    
    all_seg_ids = df_var.index.tolist()
    
    matrix_tcn_all = matrix_tcn.reindex(columns=all_seg_ids).fillna(0).astype(np.int8)
    matrix_loh_all = matrix_loh.reindex(columns=all_seg_ids).fillna(0).astype(np.int8)

    # Remove normal cells
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    
    matrix_tcn_all = matrix_tcn_all[~matrix_tcn_all.index.isin(normal_barcodes)]
    matrix_loh_all = matrix_loh_all[~matrix_loh_all.index.isin(normal_barcodes)]
    
    logging.info(f"Filtered {matrix_tcn.shape[0] - matrix_tcn_all.shape[0]} normal cells. "
                 f"{matrix_tcn_all.shape[0]} tumor cells remain for calculation.")

    # Store TCN in X and LOH in a separate layer.
    adata = ad.AnnData(
        X=matrix_tcn_all.values,                        # X stores the TCN track.
        obs=pd.DataFrame(index=matrix_tcn_all.index),
        var=df_var
    )
    
    adata.layers['LOH_Status'] = matrix_loh_all.values  # layers stores the LOH track.

    return adata
'''
def load_Numbat_cnv(bulk_clones_file, mode='bulk'):
    """
    Read Numbat bulk_clones_final.tsv(.gz) directly and return segment-level
    CNV profiles without rebuilding a cell-level AnnData object.

    Args:
    - bulk_clones_file: path to bulk_clones_final.tsv or .gz.
    - mode: 'clonal' for per-subclone profiles, or 'bulk' for the weighted
      global profile.
    """
    if not str(bulk_clones_file).endswith('.gz') and not str(bulk_clones_file).endswith('.tsv'):
        logging.warning(f"File {bulk_clones_file} does not look like a standard Numbat bulk tsv file.")
        
    df = pd.read_csv(bulk_clones_file, sep='\t')
    
    # Track 1: map TCN states to (-1, 0, 1).
    def _map_tcn_event(state):
        mapping = {
            'del': -1,  # Physical deletion.
            'neu': 0,   # Diploid neutral state.
            'loh': 0,   # CN-LOH keeps total copy number neutral here.
            'amp': 1,   # Amplification.
            'bamp': 1,  # Balanced amplification.
        }
        return mapping.get(state, 0)

    # Track 2: map LOH states to (1, 0).
    def _map_loh_event(state):
        mapping = {
            'del': 0,   
            'neu': 0,   
            'loh': 1,   # Dedicated CN-LOH detection track.
            'amp': 0,   
            'bamp': 0,
        }
        return mapping.get(state, 0)

    # Map raw CNV states into numeric tracks.
    df['CN_Score'] = df['cnv_state'].apply(_map_tcn_event)
    df['LOH_Status'] = df['cnv_state'].apply(_map_loh_event)
    
    # Standardize column names.
    rename_dict = {
        'CHROM': 'Chromosome',
        'seg_start': 'Start',
        'seg_end': 'End',
        'sample': 'Clone_ID'  # Numbat bulk outputs often store clone ids in 'sample'.
    }
    df = df.rename(columns=rename_dict)

    def _add_segment_id(df_seg):
        df_seg = df_seg.copy()
        df_seg['ID'] = (
            df_seg['Chromosome'].astype(str) + "_" +
            df_seg['Start'].astype(str) + "_" +
            df_seg['End'].astype(str)
        )
        return df_seg
    
    # Follow the expected aggregation logic for each output mode.
    if mode == 'clonal':
        # Average repeated records within each clone/segment pair.
        df_out = df.groupby(['Clone_ID', 'Chromosome', 'Start', 'End'])[['CN_Score', 'LOH_Status']].mean().reset_index()
        
        # Normalize clone ids so downstream code sees stable labels.
        df_out['Clone_ID'] = df_out['Clone_ID'].apply(lambda x: f"numbat_{x}" if not str(x).startswith("numbat") else x)
        df_out = _add_segment_id(df_out)
        return df_out
        
    elif mode == 'bulk':
        # Create a weighted pseudobulk over all cells using n_cells.
        if 'n_cells' not in df.columns:
            logging.warning("Column 'n_cells' not found in Numbat bulk file. Unweighted mean will be used.")
            df['n_cells'] = 1 
            
        df['weighted_CN'] = df['CN_Score'] * df['n_cells']
        df['weighted_LOH'] = df['LOH_Status'] * df['n_cells']
        
        grouped = df.groupby(['Chromosome', 'Start', 'End'])
        
        # Compute the weighted mean profile.
        df_out = grouped.apply(lambda x: pd.Series({
            'CN_Score': x['weighted_CN'].sum() / x['n_cells'].sum(),
            'LOH_Status': x['weighted_LOH'].sum() / x['n_cells'].sum()
        })).reset_index()

        df_out = _add_segment_id(df_out)
        return df_out
        
    else:
        raise ValueError("mode must be either 'clonal' or 'bulk'")
def load_Clonalscope(seg_file, df_subcluster=None):
    """
    Load a Clonalscope segment matrix.
    If df_subcluster is None, return a global pseudobulk profile.
    If df_subcluster is provided, return one pseudobulk profile per Clone_ID.
    """
    # load_txt wraps pd.read_csv(seg_file, sep='\t', index_col=0)
    seg = load_txt(seg_file) 

    # Parse segment coordinates.
    parsed_coords = [c.split('-') for c in seg.columns]
    chr_list = [p[0].replace('chr', '') for p in parsed_coords]
    start_list = [int(float(p[1])) for p in parsed_coords]
    end_list = [int(float(p[2])) for p in parsed_coords]
    Seg_ids = [f"{c}_{s}_{e}" for c, s, e in zip(chr_list, start_list, end_list)]

    # Branch 1: no labels provided, compute a global pseudobulk.
    if df_subcluster is None or df_subcluster.empty:
        mean_values = seg.mean(axis=0).values
        values = np.log2(np.maximum(mean_values, 0.01))
        
        df_seg = pd.DataFrame({
            "Chromosome": chr_list,
            "Start": start_list,
            "End": end_list,
            "ID": Seg_ids,
            "CN_Score": values,
            "LOH_Status": np.nan
        })
        return df_seg

    # Branch 2: labels provided, compute clonal pseudobulks.
    else:
        df_subcluster = df_subcluster.set_index('Barcodes')
        common_cells = seg.index.intersection(df_subcluster.index)
        
        if len(common_cells) == 0:
            raise ValueError("No overlapping barcodes between CNV matrix and subcluster predictions.")
            
        # Keep the overlapping cells and attach Clone_ID labels.
        seg_common = seg.loc[common_cells].copy()
        seg_common['Clone_ID'] = df_subcluster.loc[common_cells, 'Label_preds']
        
        df_list = []
        clones = seg_common['Clone_ID'].dropna().unique()
        
        for clone in clones:
            # Slice cells for the current clone and remove Clone_ID before averaging.
            subset = seg_common[seg_common['Clone_ID'] == clone].drop(columns=['Clone_ID'])
            
            # Compute the clone mean profile and convert to log2 space.
            mean_values = subset.mean(axis=0).values
            values = np.log2(np.maximum(mean_values, 0.01))
            
            df_clone = pd.DataFrame({
                "Clone_ID": clone,      # Clone-specific identifier.
                "Chromosome": chr_list,
                "Start": start_list,
                "End": end_list,
                "ID": Seg_ids,
                "CN_Score": values,
                "LOH_Status": np.nan
            })
            df_list.append(df_clone)
            
        df_out = pd.concat(df_list, ignore_index=True)
        return df_out


def load_STARCH(cnv_file, clone_file, gene_annot_path, tn_annot_path, name_mapping_file):
    def _map_cnv_event(values):
        # See: https://github.com/raphael-group/STARCH/issues/5
        # STARCH: 0=Loss, 1=Neutral, 2=Amp
        mapping = {0: -1, 1: 0, 2: 1}
        return values.map(mapping)
    df_cnv = load_csv(cnv_file, header = True)
    df_label = load_csv(clone_file, header = True).reset_index()
    df_tn_annot = load_tsv(tn_annot_path, header = True).reset_index()
    df_ann = load_gene_annot(gene_annot_path)
    name_mapping = load_tsv(name_mapping_file, header=True).reset_index()
    df_label = df_label.rename(columns={
        'index': 'barcodes',
        '0': 'clone_id'         
    })

    # transfer 0.0x16.0 to 0x16
    df_label['barcodes'] = df_label['barcodes'].apply(
        lambda x: 'x'.join([str(int(float(i))) for i in x.split('x')])
    )
    # mapping to cell barcodes using index 
    df_label = df_label.merge(name_mapping, left_on='barcodes', right_on='spot_name', how='left')

    # Remove normal cells
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    initial_count = len(df_label)
    df_label = df_label[~df_label['barcode'].isin(normal_barcodes)]
    tumor_cells_count = len(df_label)
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                 f"{tumor_cells_count} tumor cells remain for calculation.")
    
    # Mapping cnv event 
    clone_cols = [c for c in df_cnv.columns]
    df_cnv[clone_cols] = df_cnv[clone_cols].apply(_map_cnv_event)

    # Add gene annotation
    common_genes = df_cnv.index.intersection(df_ann.index)
    if len(common_genes) == 0:
        logging.error(f"No common genes found between STARCH output and annotation!")
        raise ValueError("Gene alignment failed. Check Gene_ID format.")
    df_cnv = df_cnv.loc[common_genes]
    var = df_ann.loc[common_genes]
    var.index = (var['Chromosome'].astype(str) + "_" + var['Start'].astype(str) + "_" + var['End'].astype(str))
    # Transfer cnv_file to cell by gene matrix
    df_cnv = df_cnv.reset_index()
    df_cnv = df_cnv.rename(columns={'index': 'Gene_ID'})

    clone_cols = [c for c in df_cnv.columns if c != 'Gene_ID']  
    matrix_cnv = df_cnv[clone_cols].T
    matrix_cnv.index = matrix_cnv.index.astype(int)
    clone_mapping = df_label.set_index('barcodes')['clone_id'].astype(int)  
    matrix_cnv = matrix_cnv.loc[clone_mapping].values 
    
    adata = ad.AnnData(
        X=matrix_cnv.astype(np.int8), 
        obs=df_label.set_index('barcode'),
        var=var
    )
    
    return adata

# WMA_smoothed_log_ratio_ab_dynamic log2 ratio scale
def load_Xclone_expr_mtx(adata_file, tn_annot_path):
    adata = load_h5ad(adata_file)
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    # Remove normal cells
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    initial_count = adata.n_obs
    tumor_cell_mask = adata.obs_names.isin(normal_barcodes)
    adata = adata[~tumor_cell_mask]
    tumor_cells_count = adata.n_obs
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                 f"{tumor_cells_count} tumor cells remain for calculation.")
    adata_out = ad.AnnData(
        X = adata.layers['WMA_smoothed_log_ratio_ab_dynamic'],
        obs = adata.obs.copy(),
        var = adata.var[['GeneName', 'GeneID', 'chr', 'start', 'stop']].copy()
    )
    # rename var col
    adata_out.var = adata_out.var.rename(columns={
        'GeneName': 'Gene_name',
        'GeneID': 'Gene_ID',
        'chr': 'Chromosome',
        'start': 'Start',
        'stop': 'End'
    })
    return adata_out
# posterior_mtx is the posterior likelihood of three events: deletion, LOH, neutral, amplification
def load_Xclone_cnv(adata_file, tn_annot_path):
    # See: https://github.com/single-cell-genetics/XClone/issues/15
    def assign_cna_states(cna_probs):
        """Assigns CNA states to cells based on the highest probability in a cell*gene*cnv state array.

        Args:
            cna_probs: A NumPy array of shape (num_cells, num_genes, num_cna_states)
                    representing the probabilities of each CNA state for each cell and gene.

        Returns:
            A NumPy array of shape (num_cells, num_genes) containing the assigned CNA states 
            for each cell and gene.
        """
        num_cells, num_genes, _ = cna_probs.shape
        cna_states = np.argmax(cna_probs, axis=2)  # Find the index of the highest probability state
        return cna_states

    adata = load_h5ad(adata_file)
    df_tn_annot = load_tsv(tn_annot_path, header=True).reset_index()
    
    # Remove normal cells 
    normal_barcodes = set(df_tn_annot[df_tn_annot['tumor_normal'].isin(['normal', 'Normal'])]['Barcode'])
    initial_count = adata.n_obs
    normal_cell_mask = adata.obs_names.isin(normal_barcodes)
    adata = adata[~normal_cell_mask]
    tumor_cells_count = adata.n_obs
    logging.info(f"Filtered {initial_count - tumor_cells_count} normal cells. "
                 f"{tumor_cells_count} tumor cells remain for calculation.")
                 
    adata_out = ad.AnnData(
        X = np.zeros(adata.X.shape, dtype=np.int8),  # Placeholder, will be replaced
        obs = adata.obs.copy(),
        var = adata.var[['GeneName', 'GeneID', 'chr', 'start', 'stop']].copy()
    )
    
    # rename var col
    adata_out.var = adata_out.var.rename(columns={
        'GeneName': 'Gene_name',
        'GeneID': 'Gene_ID',
        'chr': 'Chromosome',
        'start': 'Start',
        'stop': 'End'
    })
    
    # Extract posterior from posterior_mtx layers
    if 'prob1_merge' not in adata.layers:
        raise ValueError("prob1_merge layer not found in the adata.")
    
    cnv_probs = adata.layers['prob1_merge']
    
    # Derive a discrete state matrix with values in {0, 1, 2, 3}.
    cna_states = assign_cna_states(cnv_probs)
    
    # Track 1 in X: Del(-1), LOH(0), Neu(0), Amp(1).
    tcn_mapping = np.array([-1, 0, 0, 1], dtype=np.int8)
    adata_out.X = tcn_mapping[cna_states]
    
    # Track 2 in layers: only state 1 is true CN-LOH.
    loh_mapping = np.array([0, 1, 0, 0], dtype=np.int8)
    adata_out.layers['LOH_Status'] = loh_mapping[cna_states]

    return adata_out

def pseudo_bulk_adata(adata):
    
    cn_mean = np.mean(adata.X, axis=0) 
    if adata.layers.get('LOH_Status') is not None:
        loh_freq = np.mean(adata.layers['LOH_Status'], axis=0)
    else:
        loh_freq = np.nan
    df_out = pd.DataFrame({
        'Chromosome': adata.var['Chromosome'].values,
        'Start': adata.var['Start'].values,
        'End': adata.var['End'].values,
        'ID': adata.var_names.values,
        'CN_Score': cn_mean,
        'LOH_Status': loh_freq
    })
    
    return df_out

def pseudo_clonal_bulk(adata, clone_col='Clone_ID'):
    """
    Compute one pseudo-bulk CNV profile for each clone in clone_col.
    """
    if clone_col not in adata.obs.columns:
        raise ValueError(f"Column '{clone_col}' not found in adata.obs")

    df_list = []
    # Extract all valid clone ids.
    clones = adata.obs[clone_col].dropna().unique()
    
    for clone in clones:
        subset = adata[adata.obs[clone_col] == clone]
        
        # Compute CN means while handling sparse matrices safely.
        if isinstance(subset.X, np.ndarray):
            cn_mean = np.mean(subset.X, axis=0)
        else:
            cn_mean = np.asarray(subset.X.mean(axis=0)).flatten()
            
        if subset.layers.get('LOH_Status') is not None:
            if isinstance(subset.layers['LOH_Status'], np.ndarray):
                loh_freq = np.mean(subset.layers['LOH_Status'], axis=0)
            else:
                loh_freq = np.asarray(subset.layers['LOH_Status'].mean(axis=0)).flatten()
        else:
            loh_freq = np.nan
            
        df_clone = pd.DataFrame({
            'Clone_ID': clone,
            'Chromosome': subset.var['Chromosome'].values, # Keep the normalized chromosome naming for downstream plots.
            'Start': subset.var['Start'].values,
            'End': subset.var['End'].values,
            'ID': subset.var_names.values,
            'CN_Score': cn_mean,                   
            'LOH_Status': loh_freq
        })
        df_list.append(df_clone)
        
    if not df_list:
        return pd.DataFrame()
        
    df_out = pd.concat(df_list, ignore_index=True)
    return df_out

def _get_resolution_df(adata, cn_normal=2):
    '''
    Get CNV Event resolution DataFrame from AnnData
    '''

    if hasattr(adata.X, "toarray"):
        X = adata.X.toarray()
    else:
        X = adata.X
    
    mean_states = np.mean(X, axis=0)
    df = pd.DataFrame({
        'Chromosome': adata.var['Chromosome'].values,
        'Start': adata.var['Start'].values,
        'End': adata.var['End'].values,
        'State': mean_states
    })

    df = df.sort_values(['Chromosome', 'Start'])

    is_new_event = (
        (df['Chromosome'] != df['Chromosome'].shift()) | 
        (df['State'] != df['State'].shift())
    )
    df['event_id'] = is_new_event.cumsum()

    df_merged = df.groupby('event_id').agg({
        'Chromosome': 'first',
        'Start': 'min',
        'End': 'max',
        'State': 'first'
    }).reset_index(drop=True)

    # If cn_normal <= 0.99, use Log2 Copy Ratio thresholding
    if cn_normal <= 0.99:
        # Log2 Copy Ratio case
        df_cna = df_merged[np.abs(df_merged['State']) > cn_normal].copy()
    else:
        df_cna = df_merged[df_merged['State'] != cn_normal].copy()

    df_cna['Length'] = df_cna['End'] - df_cna['Start']
    
    df_final = df_cna[df_cna['Length'] > 0][['Chromosome', 'Start', 'End', 'Length']]
    
    return df_final
