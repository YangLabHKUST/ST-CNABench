import pandas as pd
import anndata as ad
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_tsv, load_CopyKAT_binbycell, pseudo_bulk_adata, pseudo_clonal_bulk, _get_resolution_df, save_tsv

class CopyKATLoader(BaseLoader):
    def __init__(self, task_name: str, eval_name: str, result_dir: str, save_dir: str):
        model_name = 'CopyKAT'
        super().__init__(model_name=model_name, eval_name=eval_name, task_name=task_name, result_dir=result_dir, save_dir=save_dir)
        
    def extract_cnv_profile(self, gene_annot_path=None, tn_annot_path=None):
        level = "Log2_Copy_Ratio"
        # raw results gene by cell is not the final output
        # use: CNA_results.txt
        pattern = "*CNA_results.txt"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            logging.warning(f"CopyKAT result not found in {self.result_dir}, returning empty DataFrame")
            return {
                "data": pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                "level": level
            }
        matrix_file = matches[0]
        adata = load_CopyKAT_binbycell(matrix_file, tn_annot_path)
        df_out = pseudo_bulk_adata(adata)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")

        return {
            "data": df_out, 
            "level": level
        }
    
    def extract_clone_cnv_profile(self, gene_annot_path=None, tn_annot_path=None):
        df_subcluster = self.extract_subcluster_preds()
        pattern = "*CNA_results.txt"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            logging.warning(f"CopyKAT result not found in {self.result_dir}, returning empty DataFrame")
            return pd.DataFrame()
        matrix_file = matches[0]
        adata = load_CopyKAT_binbycell(matrix_file, tn_annot_path)
        df_subcluster = df_subcluster.set_index('Barcodes')
        common_cells = adata.obs_names.intersection(df_subcluster.index)
        
        if len(common_cells) == 0:
            raise ValueError(f"[{self.model_name}] No overlapping barcodes between CNV matrix and subcluster predictions.")
        
        adata = adata[common_cells].copy()
        adata.obs['Clone_ID'] = df_subcluster.loc[common_cells, 'Label_preds']
        
        df_out = pseudo_clonal_bulk(adata, clone_col='Clone_ID')
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        return df_out
    def extract_cnv_resolution(self, gene_annot_path=None, tn_annot_path=None):
        pattern = "*CNA_results.txt"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            logging.warning(f"CopyKAT result not found in {self.result_dir}, returning empty resolution DataFrame")
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])
        matrix_file = matches[0]
        adata = load_CopyKAT_binbycell(matrix_file, tn_annot_path)
        df_resolution = _get_resolution_df(adata, cn_normal=0.01)
        save_tsv(df_resolution, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_resolution

    def extract_tumor_normal_preds(self):
        pattern = "*copykat_prediction.txt"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            raise FileNotFoundError(f"CopyKAT tumor or normal prediction file not found in {self.result_dir}")
        preds_file = matches[0]
        df = load_tsv(preds_file, header=True).reset_index()
        df.columns = ['Barcodes', 'Preds']
        # If CopyKAT results like c1:diploid:low.conf, c2:aneuploid:low.conf, we only keep diploid/aneuploid info
        if df['Preds'].str.contains(':').any():
            df['Preds'] = df['Preds'].apply(lambda x: x.split(':')[1] if isinstance(x, str) else x)
        mapping = {
            'aneuploid': 1,
            'diploid': 0,
            'not.defined': -1
        }
        df['Preds'] = df['Preds'].map(mapping)
        save_tsv(df, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df
    
    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        pattern = "*copykat_prediction.txt"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            raise FileNotFoundError(f"CopyKAT tumor or normal prediction file not found in {self.result_dir}")
        preds_file = matches[0]
        df = load_tsv(preds_file, header=True).reset_index()
        df.columns = ['Barcodes', 'Preds']
        mapping = {
            'aneuploid': 1,
            'diploid': 0,
            'not.defined': -1
        }
        df['Preds'] = df['Preds'].map(mapping)
        # Subset normal mode, remove ref normal cells
        df_model_run = pd.read_csv(tn_annot_path_model_run, sep='\t', header=0, comment='#')
        if df_model_run.shape[1] == 2:
            df_model_run.columns = ['Barcodes', 'tumor_normal']
        elif df_model_run.shape[1] == 4:
            df_model_run.columns = ['Barcodes', 'tumor_normal', 'leiden', 'celltype']
        normal_barcode = df_model_run[(df_model_run['tumor_normal'] == 'normal') | (df_model_run['tumor_normal'] == 'Normal')]['Barcodes'].tolist()
        df = df[~df['Barcodes'].isin(normal_barcode)]
        save_tsv(df, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df
    
    def extract_subcluster_preds(self):
        subcluster_file = self.result_dir / "copykat_subcluster_results.txt"
        if not subcluster_file.exists():
            raise FileNotFoundError(f"No CopyKAT subcluster file found in {self.result_dir}")
        # df_subcluster have processed in copykat run R script (header, label mapping)
        df_subcluster = pd.read_csv(subcluster_file, sep='\t', header=0)
        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster