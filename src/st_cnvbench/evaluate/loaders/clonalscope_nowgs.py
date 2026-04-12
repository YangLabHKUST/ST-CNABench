import pandas as pd
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_Clonalscope, save_tsv

class ClonalscopeNoWGSLoader(BaseLoader):
    def __init__(self, task_name: str, eval_name: str, result_dir: str, save_dir: str):
        model_name = 'Clonalscope_NoWGS'
        super().__init__(model_name=model_name, eval_name=eval_name, task_name=task_name, result_dir=result_dir, save_dir=save_dir)
        
    def extract_cnv_profile(self, gene_annot_path=None, tn_annot_path=None):
        level = 'Log2_Copy_Ratio'
        seg_file = self.result_dir / "clonalscope_cell_seg_matrix.tsv"
        if not seg_file.exists():
            logging.warning(f"Clonalscope_NoWGS seg matrix file not found in {self.result_dir}, returning empty DataFrame")
            return {
                'data': pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                'level': level,
            }

        df_out = load_Clonalscope(seg_file)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")

        return {
            'data': df_out,
            'level': level,
        }
    def extract_clone_cnv_profile(self, gene_annot_path=None, tn_annot_path=None):
        df_subcluster = self.extract_subcluster_preds()
        seg_file = self.result_dir / "clonalscope_cell_seg_matrix.tsv"
        if not seg_file.exists():
            logging.warning(f"Clonalscope_NoWGS seg matrix file not found in {self.result_dir}, returning empty DataFrame")
            return pd.DataFrame()
        df_out = load_Clonalscope(seg_file, df_subcluster=df_subcluster)
        
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        logging.info(f"[{self.model_name}] Successfully extracted clone CNV profiles.")
        return df_out

    
    def extract_cnv_resolution(self, gene_annot_path=None, tn_annot_path=None):
        seg_file = self.result_dir / "clonalscope_cell_seg_matrix.tsv"
        if not seg_file.exists():
            logging.warning(
                f"Clonalscope_NoWGS seg matrix file not found in {self.result_dir}, returning empty resolution DataFrame"
            )
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        df_seg = load_Clonalscope(seg_file)
        # filter normal segments (Problem: threshold)
        thres = 0.05
        df_seg = df_seg[np.abs(df_seg['CN_Score']) >= thres]
        df_resolution = pd.DataFrame({
            'Chromosome': df_seg['Chromosome'],
            'Start': df_seg['Start'],
            'End': df_seg['End'],
            'Length': df_seg['End'] - df_seg['Start']
        })
        save_tsv(df_resolution, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_resolution
    
    def extract_tumor_normal_preds(self):
        annot_file = self.result_dir / "clonalscope_final_tn_assignment.tsv"
        if not annot_file.exists():
            raise FileNotFoundError(f"Clonalscope_NoWGS cell annotation file not found in {self.result_dir}")
        df_annot = pd.read_csv(annot_file, sep="\t")
        df_annot.columns = ['Barcodes', 'Preds']
        df_annot['Preds'] = df_annot['Preds'].map({'Tumor': 1, 'Normal': 0, 'Unknown': -1})
        save_tsv(df_annot, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_annot
    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        annot_file = self.result_dir / "clonalscope_final_tn_assignment.tsv"
        if not annot_file.exists():
            raise FileNotFoundError(f"Clonalscope_NoWGS cell annotation file not found in {self.result_dir}")

        df_annot = pd.read_csv(annot_file, sep="\t")
        df_annot.columns = ['Barcodes', 'Preds']
        df_model_run = pd.read_csv(tn_annot_path_model_run, sep='\t', header=0, comment='#')
        if df_model_run.shape[1] == 2:
            df_model_run.columns = ['Barcodes', 'tumor_normal']
        elif df_model_run.shape[1] == 4:
            df_model_run.columns = ['Barcodes', 'tumor_normal', 'leiden', 'celltype']
        normal_barcode = df_model_run[(df_model_run['tumor_normal'] == 'normal') | (df_model_run['tumor_normal'] == 'Normal')]['Barcodes'].tolist()
        df = df_annot[~df_annot['Barcodes'].isin(normal_barcode)].copy()
        df['Preds'] = df['Preds'].map({'Tumor': 1, 'Normal': 0, 'Unknown': -1})
        save_tsv(df, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df
    
    def extract_subcluster_preds(self):
        subcluster_file = self.result_dir / "clonalscope_zest_clusters.tsv"
        if not subcluster_file.exists():
            raise FileNotFoundError(f"Clonalscope_NoWGS subcluster annotation file not found in {self.result_dir}")
        df_subcluster = pd.read_csv(subcluster_file, sep="\t", header=0)
        df_subcluster.columns = ['Barcodes', 'Label_preds']
        unique_labels = sorted(df_subcluster['Label_preds'].unique())
        label_mapping = {label: f"{self.model_name}_{i+1}" for i, label in enumerate(unique_labels)}
        df_subcluster['Label_preds'] = df_subcluster['Label_preds'].map(label_mapping)
        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster
