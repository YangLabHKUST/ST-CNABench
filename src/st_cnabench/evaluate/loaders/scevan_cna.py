import pandas as pd
import numpy as np
import logging
import re
from .base import BaseLoader
from ..utils.io import load_SCEVAN_cna, load_tsv, save_tsv

class SCEVANCNALoader(BaseLoader):
    def __init__(self, task_name: str, eval_name: str, result_dir: str, save_dir: str):
        model_name = 'SCEVAN'
        super().__init__(model_name=model_name, eval_name=eval_name, task_name=task_name, result_dir=result_dir, save_dir=save_dir)
        
    def extract_cna_profile(self, gene_annot_path=None, tn_annot_path=None):
        level = 'Integer_CN'
        seg_matches = list(self.result_dir.glob("**/output/*_Clonal_CN.seg"))
        if not seg_matches:
            logging.warning(f"SCEVAN clonal CN seg file not found in {self.result_dir}, returning empty DataFrame")
            return {
                "data": pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                "level": level
             }        
        df_out = load_SCEVAN_cna(seg_matches[0])
        df_out['CN_Score'] = np.log2(np.maximum(df_out['CN_Score'], 0.1) / 2.0)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")
        return {
            "data": df_out,
            "level": level
        }

    def extract_clone_cna_profile(self, gene_annot_path=None, tn_annot_path=None):
        clonal_seg_matches = list(self.result_dir.glob("**/output/*_subclone*_CN.seg"))
        if not clonal_seg_matches:
            logging.warning(f"[{self.model_name}] clonal CN seg file not found in {self.result_dir}, returning empty DataFrame")
            return pd.DataFrame()
        
        clonal_segs = []
        
        for match in clonal_seg_matches:
            m = re.search(r'_subclone(\d+)_CN\.seg', match.name)
            if not m:
                continue
                
            subclone_idx = m.group(1)
            clone_id = f"{self.model_name}_{subclone_idx}"
            
            df_seg = load_SCEVAN_cna(match)
            if df_seg is None or df_seg.empty:
                continue
            df_seg['Clone_ID'] = clone_id
            clonal_segs.append(df_seg)

        if not clonal_segs:
            return pd.DataFrame()

        df_out = pd.concat(clonal_segs, ignore_index=True)
        df_out['CN_Score'] = np.log2(np.maximum(df_out['CN_Score'], 0.1) / 2.0)
        
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        return df_out
    def extract_cna_resolution(self, gene_annot_path=None, tn_annot_path=None):
        seg_matches = list(self.result_dir.glob("**/output/*_Clonal_CN.seg"))
        if not seg_matches:
            logging.warning(f"SCEVAN clonal CN seg file not found in {self.result_dir}, returning empty resolution DataFrame")
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        df_seg = load_SCEVAN_cna(seg_matches[0])
        # filter normal segments
        df_seg = df_seg[df_seg['CN_Score'] != 2]
        df_resolution = pd.DataFrame({
            'Chromosome': df_seg['Chromosome'],
            'Start': df_seg['Start'],
            'End': df_seg['End'],
            'Length': df_seg['End'] - df_seg['Start']
        })
        save_tsv(df_resolution, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_resolution

    def extract_tumor_normal_preds(self):
        pattern = "*SCEVAN_res_meta.tsv"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            raise FileNotFoundError(f"SCEVAN tumor or normal result not found in {self.result_dir}")
        preds_file = matches[0]
        df = load_tsv(preds_file, header=True).reset_index()
        if df.shape[1] == 4:
            df.columns = ['Barcodes', 'Preds', 'ConfidentNormal', 'Subclone']
        # If No subclone detect
        elif df.shape[1] == 3:
            df.columns = ['Barcodes', 'Preds', 'ConfidentNormal']
        mapping = {
            'tumor': 1,
            'normal': 0,
            'filtered': -1
        }
        df['Preds'] = df['Preds'].map(mapping)
        save_tsv(df, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df
    
    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        pattern = "*SCEVAN_res_meta.tsv"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            raise FileNotFoundError(f"SCEVAN tumor or normal result not found in {self.result_dir}")
        preds_file = matches[0]
        df = load_tsv(preds_file, header=True).reset_index()
        if df.shape[1] == 4:
            df.columns = ['Barcodes', 'Preds', 'ConfidentNormal', 'Subclone']
        # If No subclone detect
        elif df.shape[1] == 3:
            df.columns = ['Barcodes', 'Preds', 'ConfidentNormal']
        mapping = {
            'tumor': 1,
            'normal': 0,
            'filtered': -1
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
        pattern = "*SCEVAN_res_meta.tsv"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            raise FileNotFoundError(f"SCEVAN tumor or normal result not found in {self.result_dir}")
        meta_file = matches[0]
        df_meta = load_tsv(meta_file, header=True).reset_index()
        # Valid df_meta have column 'subclone' for subcluster prediction, if not exist, return empty df
        if 'subclone' not in df_meta.columns:
            logging.warning(f"[{self.model_name}] detected no significant subclone, no subclone column in {meta_file}")
            df_subcluster = pd.DataFrame(columns=['Barcodes', 'Label_preds'])
        else:
            # First drop confidentNomral == yes cells, then extract subcluster preds
            df_meta = df_meta[df_meta['confidentNormal'] != 'yes']
            df_subcluster = df_meta[['Barcode', 'subclone']].copy()
            df_subcluster.columns = ['Barcodes', 'Label_preds']
            df_subcluster['Label_preds'] = (
                pd.to_numeric(df_subcluster['Label_preds'], errors='coerce')
                .fillna(0)
                .astype(int)
            )
            # Transform subcluster label to format: {MODELNAME}_{i}
            df_subcluster['Label_preds'] = df_subcluster['Label_preds'].apply(lambda x: f"{self.model_name}_{x}")
        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster
