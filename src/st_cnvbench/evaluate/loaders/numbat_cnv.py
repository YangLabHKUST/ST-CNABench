import pandas as pd
import anndata as ad
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_Numbat_cnv, load_tsv, pseudo_bulk_adata, pseudo_clonal_bulk, _get_resolution_df, save_tsv
from ..utils.constants import CHR_SIZES_HG38

class NumbatCNVLoader(BaseLoader):
    _NO_CNV_LOG_MARKER = "No CNV remains after filtering by LLR in pseudobulks."

    def __init__(self, task_name:str, eval_name:str, result_dir: str, save_dir: str):
        model_name = 'Numbat'
        super().__init__(model_name=model_name,eval_name=eval_name,task_name=task_name,result_dir=result_dir, save_dir=save_dir)

    def _has_no_cnv_log_marker(self):
        log_files = sorted(self.result_dir.glob("**/*.log"))
        if not log_files:
            return False

        for log_path in log_files:
            try:
                with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        if self._NO_CNV_LOG_MARKER in line:
                            logging.info(f"Numbat no-CNV marker found in log: {log_path}")
                            return True
            except OSError as e:
                logging.warning(f"Failed to read log file {log_path}: {e}")
        return False

    def _build_zero_cnv_profile(self):
        rows = []
        for chrom, chrom_len in CHR_SIZES_HG38.items():
            rows.append({
                "Chromosome": str(chrom),
                "Start": 0,
                "End": int(chrom_len),
                "ID": f"{chrom}_0_{int(chrom_len)}",
                "CN_Score": 0.0,
                "LOH_Status": 0.0,
            })
        return pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "ID", "CN_Score", "LOH_Status"])
        
    def extract_cnv_profile(self, gene_annot_path=None, tn_annot_path=None):
        level = 'CNV_Event'
        # Read the bulk_clones_final output directly when available.
        bulk_files = list(self.result_dir.glob("bulk_clones_final.tsv*"))
        if not bulk_files:
            if self._has_no_cnv_log_marker():
                logging.info(
                    f"No bulk file found in {self.result_dir}, but no-CNV marker detected in log. "
                    "Returning neutral (all-zero) CN/LOH profile."
                )
                df_out = self._build_zero_cnv_profile()
                save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")
                return {'data': df_out, 'level': level}
            logging.warning(f"No Numbat bulk clones file found in {self.result_dir}, returning empty DataFrame")
            empty_df = pd.DataFrame(columns=["Chromosome", "Start", "End", "ID", "CN_Score", "LOH_Status"])
            return {'data': empty_df, 'level': level}
            
        bulk_file = sorted(bulk_files)[-1]
        
        df_out = load_Numbat_cnv(bulk_file, mode='bulk')
        
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")
        return {'data': df_out, 'level': level}

    def extract_clone_cnv_profile(self, gene_annot_path=None, tn_annot_path=None):
        bulk_files = list(self.result_dir.glob("bulk_clones_final.tsv*"))
        if not bulk_files:
            logging.warning(f"No Numbat bulk clones file found in {self.result_dir}")
            return pd.DataFrame()
            
        bulk_file = sorted(bulk_files)[-1]
        
        df_out = load_Numbat_cnv(bulk_file, mode='clonal')
        
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        return df_out
    def extract_cnv_resolution(self, gene_annot_path, tn_annot_path):
        bulk_files = list(self.result_dir.glob("bulk_clones_final.tsv*"))
        if not bulk_files:
            logging.warning(f"No Numbat bulk clones file found in {self.result_dir}, returning empty resolution DataFrame")
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        bulk_file = sorted(bulk_files)[-1]
        df = load_Numbat_cnv(bulk_file, mode='bulk')
        df = df[np.abs(df['CN_Score']) != 0]
        df_resolution = pd.DataFrame({
            'Chromosome': df['Chromosome'],
            'Start': df['Start'],
            'End': df['End'],
            'Length': df['End'] - df['Start']
        })
        save_tsv(df_resolution, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_resolution

    def extract_tumor_normal_preds(self):
        pattern = "clone_post_*.tsv"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            logging.warning(f"Numbat tumor or normal prediction file not found in {self.result_dir}, returning empty DataFrame")
            return None
        preds_file = sorted(matches)[-1]
        df_clone = load_tsv(preds_file, header=True).reset_index()
        df = df_clone[['cell', 'compartment_opt']].copy()
        df.columns = ['Barcodes', 'Preds']
        mapping = {
            'tumor': 1,
            'normal': 0
        }
        df['Preds'] = df['Preds'].map(mapping)
        save_tsv(df, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df

    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        pattern = "clone_post_*.tsv"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            logging.warning(f"Numbat tumor or normal prediction file not found in {self.result_dir}, returning empty DataFrame")
            return None
        preds_file = sorted(matches)[-1]
        df_clone = load_tsv(preds_file, header=True).reset_index()
        df = df_clone[['cell', 'compartment_opt']].copy()
        df.columns = ['Barcodes', 'Preds']
        
        # load tn_annot_path_model_run to remove ref normal cells
        df_model_run = pd.read_csv(tn_annot_path_model_run, sep='\t', header=0, comment='#')
        if df_model_run.shape[1] == 2:
            df_model_run.columns = ['Barcodes', 'tumor_normal']
        elif df_model_run.shape[1] == 4:
            df_model_run.columns = ['Barcodes', 'tumor_normal', 'leiden', 'celltype']
        normal_barcode = df_model_run[(df_model_run['tumor_normal'] == 'normal') | (df_model_run['tumor_normal'] == 'Normal')]['Barcodes'].tolist()
        df_out = df[~df['Barcodes'].isin(normal_barcode)].copy()
        df_out['Preds'] = df_out['Preds'].map({'normal':0, 'tumor':1})
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_out
    
    def extract_subcluster_preds(self):
        subcluster_files = list(self.result_dir.glob("clone_post_*.tsv"))
        if not subcluster_files:
            raise FileNotFoundError(f"No Numbat subcluster prediction file found in {self.result_dir} matching clone_post_*.tsv")
        subcluster_file = sorted(subcluster_files)[-1] # Use the last matched file
        df_clone = load_tsv(subcluster_file, header=True).reset_index()
        df_subcluster = df_clone[['cell', 'clone_opt']].copy()
        df_subcluster.columns = ['Barcodes', 'Label_preds']
        unique_labels = sorted(df_subcluster['Label_preds'].unique())
        label_mapping = {label: f"numbat_{i}" for i, label in enumerate(unique_labels, start=1)}
        df_subcluster['Label_preds'] = df_subcluster['Label_preds'].map(label_mapping)
        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster
