import pandas as pd
import anndata as ad
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_Xclone_cnv, pseudo_bulk_adata, pseudo_clonal_bulk, save_tsv


class XcloneCNVLoader(BaseLoader):
    def __init__(self, task_name:str, eval_name:str, result_dir: str, save_dir: str):
        model_name = 'Xclone'
        super().__init__(model_name=model_name,eval_name=eval_name,task_name=task_name,result_dir=result_dir, save_dir=save_dir)

    def extract_cnv_profile(self, gene_annot_path, tn_annot_path):
        level = 'CNV_Event'
        adata_file = self.result_dir /"xclone_output" /"data" / "combined_final.h5ad"
        if not adata_file.exists():
            logging.warning(f"No Xclone output file found in {self.result_dir}, returning empty DataFrame")
            return {
                "data": pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                "level": level
            }
        adata = load_Xclone_cnv(adata_file, tn_annot_path)
        df_out = pseudo_bulk_adata(adata)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")
        return {
            'data': df_out,
            'level': level
        }

    def extract_clone_cnv_profile(self, gene_annot_path=None, tn_annot_path=None):
        if tn_annot_path is None:
            logging.warning(f"[{self.model_name}] tn_annot_path is required for clone CNV extraction, returning empty DataFrame")
            return pd.DataFrame()

        adata_file = self.result_dir / "xclone_output" / "data" / "combined_final.h5ad"
        if not adata_file.exists():
            logging.warning(f"No Xclone output file found in {self.result_dir}, returning empty DataFrame")
            return pd.DataFrame()

        try:
            df_subcluster = self.extract_subcluster_preds()
        except FileNotFoundError as e:
            logging.warning(f"[{self.model_name}] {e}")
            return pd.DataFrame()

        if df_subcluster is None or df_subcluster.empty:
            logging.warning(f"[{self.model_name}] Empty subcluster predictions, returning empty DataFrame")
            return pd.DataFrame()

        adata = load_Xclone_cnv(adata_file, tn_annot_path)
        df_subcluster = df_subcluster.set_index('Barcodes')
        common_cells = adata.obs_names.intersection(df_subcluster.index)

        if len(common_cells) == 0:
            logging.warning(f"[{self.model_name}] No overlapping barcodes between CNV matrix and subcluster predictions.")
            return pd.DataFrame()

        adata = adata[common_cells].copy()
        adata.obs['Clone_ID'] = df_subcluster.loc[common_cells, 'Label_preds'].values
        df_out = pseudo_clonal_bulk(adata, clone_col='Clone_ID')
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        return df_out

    def extract_cnv_resolution(self, gene_annot_path, tn_annot_path):
        adata_file = self.result_dir / "xclone_output" / "data" / "combined_final.h5ad"
        if not adata_file.exists():
            logging.warning(f"No Xclone output file found in {self.result_dir}, returning empty resolution DataFrame")
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        adata = load_Xclone_cnv(adata_file, tn_annot_path)
        # for Xclone, use cn_normal=1 (neutral)
        from ..utils.io import _get_resolution_df
        df_resolution = _get_resolution_df(adata, cn_normal=0)
        save_tsv(df_resolution, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_resolution

    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        # Use 2-clone output only for tumor/normal task.
        preds_file_pattern = "xclone_output*_2clones*clone_anno.csv"
        matches = list(self.result_dir.glob(f"**/{preds_file_pattern}"))
        if not matches:
            raise FileNotFoundError(f"No Xclone 2-clone annotation file found in {self.result_dir} matching {preds_file_pattern}")
        preds_file = sorted(matches)[-1]
        df_clone = pd.read_csv(preds_file, header=0, index_col=0)
        if 'clone(2)' not in df_clone.columns:
            raise ValueError(f"Column 'clone(2)' not found in {preds_file}")

        df = df_clone[['clone(2)']].copy().reset_index()
        df.columns = ['Barcodes', 'Preds']

        # Remove ref normal from model run annotation
        df_model_run = pd.read_csv(tn_annot_path_model_run, sep="\t", header=0)
        if df_model_run.shape[1] == 2:
            df_model_run.columns = ['Barcodes', 'tumor_normal']
        elif df_model_run.shape[1] == 4:
            df_model_run.columns = ['Barcodes', 'tumor_normal', 'leiden', 'celltype']
        normal_barcode = df_model_run[(df_model_run['tumor_normal'] == 'normal') | (df_model_run['tumor_normal'] == 'Normal')]['Barcodes'].tolist()
        df = df[~df['Barcodes'].isin(normal_barcode)]
        logging.info(f">>> [Xclone] Removed {len(normal_barcode)} normal cells from tumor/normal prediction based on model run annotation.")

        # Output is Clone1 and Clone2, map to 0/1
        mapping = {
            'Clone1': 0,
            'Clone2': 1
        }
        df['Preds'] = df['Preds'].map(mapping)
        save_tsv(df, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df

    def extract_subcluster_preds(self):
        # Use subcluster output only for subclone task.
        subcluster_pattern = "xclone_output*_subcluster*clone_anno.csv"
        matches = list(self.result_dir.glob(f"**/{subcluster_pattern}"))
        if not matches:
            raise FileNotFoundError(f"No Xclone subcluster annotation file found in {self.result_dir} matching {subcluster_pattern}")

        subcluster_file = sorted(matches)[-1]
        df_clone = pd.read_csv(subcluster_file, header=0, index_col=0)
        if 'subcluster' not in df_clone.columns:
            raise ValueError(f"Column 'subcluster' not found in {subcluster_file}")

        df_subcluster = df_clone[['subcluster']].copy().reset_index()
        df_subcluster.columns = ['Barcodes', 'Label_preds']

        raw_labels = df_subcluster['Label_preds'].astype(str)

        def _sort_key(label):
            digits = ''.join(ch for ch in label if ch.isdigit())
            return (0, int(digits)) if digits else (1, label)

        unique_labels = sorted(raw_labels.unique(), key=_sort_key)
        label_mapping = {label: f"{self.model_name}_{i+1}" for i, label in enumerate(unique_labels)}
        df_subcluster['Label_preds'] = raw_labels.map(label_mapping)

        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster
