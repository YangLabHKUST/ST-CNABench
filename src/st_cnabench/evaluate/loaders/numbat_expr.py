import pandas as pd
import anndata as ad
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_Numbat_expr, pseudo_bulk_adata, _get_resolution_df, save_tsv


class NumbatExprLoader(BaseLoader):
    def __init__(self, task_name: str, eval_name: str, result_dir: str, save_dir: str):
        model_name = 'Numbat'
        super().__init__(model_name=model_name, eval_name=eval_name, task_name=task_name, result_dir=result_dir, save_dir=save_dir)

    def extract_cna_profile(self, gene_annot_path=None, tn_annot_path=None):
        level = 'Log2_Copy_Ratio'
        matrix_file = self.result_dir / 'gexp_roll_wide.tsv.gz'
        if not matrix_file.exists():
            logging.warning(f'No Numbat expression file found in {self.result_dir}')
            return {
                'data': pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                'level': level,
            }

        adata = load_Numbat_expr(matrix_file, gene_annot_path, tn_annot_path)
        df_out = pseudo_bulk_adata(adata)
        save_tsv(df_out, self.save_dir / f'{self.eval_name}_{self.task_name}_{level}.tsv')

        return {
            'data': df_out,
            'level': level,
        }

    def extract_cna_resolution(self, gene_annot_path=None, tn_annot_path=None):
        matrix_file = self.result_dir / 'gexp_roll_wide.tsv.gz'
        if not matrix_file.exists():
            logging.warning(f'No Numbat expression file found in {self.result_dir}, returning empty resolution DataFrame')
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        adata = load_Numbat_expr(matrix_file, gene_annot_path, tn_annot_path)
        df_resolution = _get_resolution_df(adata, cn_normal=0.01)
        save_tsv(df_resolution, self.save_dir / f'{self.eval_name}_{self.task_name}.tsv')
        return df_resolution
