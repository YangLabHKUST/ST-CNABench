import pandas as pd
import anndata as ad
import logging
import numpy as np
from .base import BaseLoader
from ..utils.io import load_InferCNV_cna, pseudo_bulk_adata, pseudo_clonal_bulk, _get_resolution_df, save_tsv

class InferCNVCNALoader(BaseLoader):
    def __init__(self, task_name:str, eval_name:str, result_dir: str, save_dir: str):
        model_name = 'InferCNV'
        super().__init__(model_name=model_name,eval_name=eval_name,task_name=task_name,result_dir=result_dir, save_dir=save_dir)
        
    def extract_cna_profile(self, gene_annot_path, tn_annot_path):
        level = 'Integer_CN'
        # InferCNV file: infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt (0, 0.5, 1, 1.5, 2, 2.5, 3)
        #matrix_file = self.result_dir /"infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observations.txt"
        matrix_file = self.result_dir / "infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
        if not matrix_file.exists():
            logging.warning(f"No InferCNV observation file found in {self.result_dir}, returning empty DataFrame")
            return {
                "data": pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                "level": level
            }
        adata = load_InferCNV_cna(matrix_file, gene_annot_path, tn_annot_path)

        df_out = pseudo_bulk_adata(adata)
        df_out['CN_Score'] = np.log2(np.maximum(df_out['CN_Score'], 0.1) / 2.0)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")

        return {
            "data": df_out, 
            "level": level
        }

    def extract_clone_cna_profile(self, gene_annot_path, tn_annot_path):
        df_subcluster = self.extract_subcluster_preds()
        matrix_file = self.result_dir / "infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
        
        if not matrix_file.exists():
            logging.warning(f"No InferCNV observation file found in {self.result_dir}, returning empty DataFrame")
            return pd.DataFrame() 

        adata = load_InferCNV_cna(matrix_file, gene_annot_path, tn_annot_path)
        df_subcluster = df_subcluster.set_index('Barcodes')
        common_cells = adata.obs_names.intersection(df_subcluster.index)
        
        if len(common_cells) == 0:
            raise ValueError(f"[{self.model_name}] No overlapping barcodes between CNA matrix and subcluster predictions.")
        
        adata = adata[common_cells].copy()
        adata.obs['Clone_ID'] = df_subcluster.loc[common_cells, 'Label_preds']
        
        df_out = pseudo_clonal_bulk(adata, clone_col='Clone_ID')
        
        if not df_out.empty:
            df_out['CN_Score'] = np.log2(np.maximum(df_out['CN_Score'], 0.1) / 2.0)
        
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        return df_out
    def extract_cna_resolution(self, gene_annot_path, tn_annot_path):
        matrix_file = self.result_dir / "infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
        if not matrix_file.exists():
            logging.warning(f"No InferCNV observation file found in {self.result_dir}, returning empty resolution DataFrame")
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        adata = load_InferCNV_cna(matrix_file, gene_annot_path, tn_annot_path)
        df_resolution = _get_resolution_df(adata, cn_normal=2)
        save_tsv(df_resolution, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_resolution

    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        from ..utils.find_normal import find_normal_tumor
        # extract tumor/normal preds using model run annot to remove ref normal cells
        ref_mtx_file = self.result_dir / "infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"
        if not ref_mtx_file.exists():
            raise FileNotFoundError(f"No InferCNV reference file found in {self.result_dir}")
        obs_mtx_file = self.result_dir / "infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
        if not obs_mtx_file.exists():
            raise FileNotFoundError(f"No InferCNV observation file found in {self.result_dir}")
        df_ref = pd.read_csv(ref_mtx_file, sep=r'\s+', index_col=0)
        df_obs = pd.read_csv(obs_mtx_file, sep=r'\s+', index_col=0)
        preds_series = find_normal_tumor(df_ref, df_obs, mode='ref')
        df_preds = preds_series.reset_index()
        df_preds.columns = ['Barcodes', 'Preds']
        df_preds['Preds'] = df_preds['Preds'].map({'normal':0, 'tumor':1})
        save_tsv(df_preds, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_preds
    
    def extract_subcluster_preds(self):
        subcluster_file = self.result_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"
        if not subcluster_file.exists():
            raise FileNotFoundError(f"No InferCNV subcluster file found in {self.result_dir}")
        df_subcluster = pd.read_csv(subcluster_file, sep='\t', header=0)
        df_subcluster = df_subcluster.rename(columns={'cell': 'Barcodes', 'cell_group_name': 'Label_preds'})
        df_subcluster = df_subcluster[['Barcodes', 'Label_preds']]
        # Mapping subcluster labels to infercnv_integers such as infercnv_1, infercnv_2, etc.
        # remove normal cluster first (Normal.Normal_s1)
        df_subcluster = df_subcluster[~df_subcluster['Label_preds'].str.contains('Normal.Normal')]
        unique_labels = sorted(df_subcluster['Label_preds'].unique())
        label_mapping = {label: f"infercnv_{i+1}" for i, label in enumerate(unique_labels)}
        df_subcluster['Label_preds'] = df_subcluster['Label_preds'].map(label_mapping)
        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster
        
