import pandas as pd
import anndata as ad
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_Xclone_expr_mtx, pseudo_bulk_adata, _get_resolution_df, save_tsv


class XcloneExprLoader(BaseLoader):
    def __init__(self, task_name: str, eval_name: str, result_dir: str, save_dir: str):
        model_name = 'Xclone'
        super().__init__(model_name=model_name, eval_name=eval_name, task_name=task_name, result_dir=result_dir, save_dir=save_dir)

    def extract_cna_profile(self, gene_annot_path, tn_annot_path):
        level = 'Log2_Copy_Ratio'
        adata_file = self.result_dir / 'xclone_output' / 'data' / 'combined_final.h5ad'
        if not adata_file.exists():
            logging.warning(f'No Xclone output file found in {self.result_dir}, returning empty DataFrame')
            return {
                'data': pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                'level': level,
            }

        adata = load_Xclone_expr_mtx(adata_file, tn_annot_path)
        df_out = pseudo_bulk_adata(adata)
        save_tsv(df_out, self.save_dir / f'{self.eval_name}_{self.task_name}_{level}.tsv')
        return {
            'data': df_out,
            'level': level,
        }

    def extract_cna_resolution(self, gene_annot_path, tn_annot_path):
        adata_file = self.result_dir / 'xclone_output' / 'data' / 'combined_final.h5ad'
        if not adata_file.exists():
            logging.warning(f'No Xclone output file found in {self.result_dir}, returning empty resolution DataFrame')
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        adata = load_Xclone_expr_mtx(adata_file, tn_annot_path)
        df_resolution = _get_resolution_df(adata, cn_normal=0.01)
        save_tsv(df_resolution, self.save_dir / f'{self.eval_name}_{self.task_name}.tsv')
        return df_resolution
