import pandas as pd
import anndata as ad
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_CalicoST_intcn, load_tsv, pseudo_bulk_adata, pseudo_clonal_bulk, _get_resolution_df, save_tsv

class CalicoSTCNALoader(BaseLoader):
    def __init__(self, task_name:str, eval_name:str, result_dir: str, save_dir: str):
        model_name = 'CalicoST'
        super().__init__(model_name=model_name,eval_name=eval_name,task_name=task_name,result_dir=result_dir, save_dir=save_dir)

    def extract_cna_profile(self, gene_annot_path, tn_annot_path):
        level = "Integer_CN"
        cna_file = self.result_dir / "clone3_rectangle0_w1.0/cnv_seglevel.tsv"
        clone_label_file = self.result_dir / "clone3_rectangle0_w1.0/clone_labels.tsv"
        if not cna_file.exists():
            df_out = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Bin_ID', 'CN_Score', 'LOH_Status'])
            logging.warning(f"No CalicoST CNA segment file found at {cna_file}, returning empty DataFrame")
            return {
                "data": df_out,
                "level": level
            }
        if not clone_label_file.exists():
            raise FileNotFoundError(f"No CalicoST label file found at {clone_label_file}")

        adata = load_CalicoST_intcn(cna_file, clone_label_file, tn_annot_path)
        df_out = pseudo_bulk_adata(adata)
        df_out['CN_Score'] = np.log2(np.maximum(df_out['CN_Score'], 0.1) / 2.0)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")
        return {
            "data": df_out,
            "level": level
        }

    def extract_clone_cna_profile(self, gene_annot_path, tn_annot_path):
        df_subcluster = self.extract_subcluster_preds()
        cna_file = self.result_dir / "clone3_rectangle0_w1.0/cnv_seglevel.tsv"
        clone_label_file = self.result_dir / "clone3_rectangle0_w1.0/clone_labels.tsv"
        if not cna_file.exists():
            logging.warning(f"No CalicoST CNA segment file found at {cna_file}, returning empty DataFrame")
            return pd.DataFrame()
        if not clone_label_file.exists():
            raise FileNotFoundError(f"No CalicoST label file found at {clone_label_file}")

        adata = load_CalicoST_intcn(cna_file, clone_label_file, tn_annot_path=tn_annot_path)
        df_subcluster = df_subcluster.set_index('Barcodes')
        common_cells = adata.obs_names.intersection(df_subcluster.index)
        if len(common_cells) == 0:
            raise ValueError(f"[{self.model_name}] No overlapping barcodes between CNA matrix and subcluster predictions.")
        adata = adata[common_cells].copy()
        adata.obs['Clone_ID'] = df_subcluster.loc[common_cells, 'Label_preds']
        df_out = pseudo_clonal_bulk(adata, clone_col='Clone_ID')
        df_out['CN_Score'] = np.log2(np.maximum(df_out['CN_Score'], 0.1) / 2.0)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        return df_out

    def extract_cna_resolution(self, gene_annot_path, tn_annot_path):
        cna_file = self.result_dir / "clone3_rectangle0_w1.0/cnv_seglevel.tsv"
        clone_label_file = self.result_dir / "clone3_rectangle0_w1.0/clone_labels.tsv"
        if not cna_file.exists():
            logging.warning(f"No CalicoST CNA segment file found at {cna_file}, returning empty resolution DataFrame")
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])
        if not clone_label_file.exists():
            logging.warning(f"No CalicoST label file found at {clone_label_file}, returning empty resolution DataFrame")
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])
        adata = load_CalicoST_intcn(cna_file, clone_label_file, tn_annot_path=tn_annot_path)

        df_resolution = _get_resolution_df(adata, cn_normal=2)
        save_tsv(df_resolution, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_resolution

    def extract_tumor_normal_preds(self):
        tumor_purity_file = self.result_dir / "loh_estimator_tumor_prop.tsv"
        if not tumor_purity_file.exists():
            logging.warning(f'No CalicoST tumor purity file found at {tumor_purity_file}, returning empty DataFrame')
            return None
        df_tumor_purity = load_tsv(tumor_purity_file, header=True).reset_index()
        df_tumor_purity.columns = ['Barcodes', 'Preds']
        #df_tumor_purity['Preds'] = pd.to_numeric(df_tumor_purity['Preds'])
        df_tumor_purity['Preds'] = (df_tumor_purity['Preds'] > 0.4).astype(int)
        save_tsv(df_tumor_purity, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_tumor_purity

    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        tumor_purity_file = self.result_dir / "loh_estimator_tumor_prop.tsv"
        if not tumor_purity_file.exists():
            logging.warning(f'No CalicoST tumor purity file found at {tumor_purity_file}, returning empty DataFrame')
            return None
        df_tumor_purity = load_tsv(tumor_purity_file, header=True).reset_index()
        df_tumor_purity.columns = ['Barcodes', 'Preds']
        # Remove ref normal from model run annotation
        df_model_run = pd.read_csv(tn_annot_path_model_run, sep="\t", header=0)
        if df_model_run.shape[1] == 2:
            df_model_run.columns = ['Barcodes', 'tumor_normal']
        elif df_model_run.shape[1] == 4:
            df_model_run.columns = ['Barcodes', 'tumor_normal', 'leiden', 'celltype']
        normal_barcode = df_model_run[(df_model_run['tumor_normal'] == 'normal') | (df_model_run['tumor_normal'] == 'Normal')]['Barcodes'].tolist()
        df_tumor_purity = df_tumor_purity[~df_tumor_purity['Barcodes'].isin(normal_barcode)]
        df_tumor_purity['Preds'] = (df_tumor_purity['Preds'] > 0.4).astype(int)
        save_tsv(df_tumor_purity, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_tumor_purity

    def extract_subcluster_preds(self):
        clone_label_file = self.result_dir / "clone3_rectangle0_w1.0/clone_labels.tsv"
        if not clone_label_file.exists():
            raise FileNotFoundError(f"No CalicoST label file found at {clone_label_file}")
        df_clone = load_tsv(clone_label_file, header=True).reset_index()
        df_subcluster = df_clone[['BARCODES', 'clone_label']].copy()
        df_subcluster.columns = ['Barcodes', 'Label_preds']
        # Map clone_label 0, 1, 2 to calicost_0, calicost_1, calicost_2
        df_subcluster['Label_preds'] = df_subcluster['Label_preds'].apply(lambda x: f"calicost_{x}")
        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster
