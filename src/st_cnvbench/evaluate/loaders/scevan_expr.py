import pandas as pd
import anndata as ad
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_SCEVAN_expr, pseudo_bulk_adata, _get_resolution_df, save_tsv


class SCEVANExprLoader(BaseLoader):
    def __init__(self, task_name: str, eval_name: str, result_dir: str, save_dir: str):
        model_name = 'SCEVAN'
        super().__init__(model_name=model_name, eval_name=eval_name, task_name=task_name, result_dir=result_dir, save_dir=save_dir)

    def _resolve_input_files(self):
        mat_matches = list(self.result_dir.glob('**/output/*_CNAmtxSubclones.RData'))
        ann_matches = list(self.result_dir.glob('**/output/*_count_mtx_annot.RData'))
        matrix_file = mat_matches[0] if mat_matches else None
        gene_ann_file = ann_matches[0] if ann_matches else None
        return matrix_file, gene_ann_file

    def extract_cnv_profile(self, gene_annot_path=None, tn_annot_path=None):
        level = 'Log2_Copy_Ratio'
        matrix_file, gene_ann_file = self._resolve_input_files()

        if matrix_file is None:
            logging.warning(f'SCEVAN CNA matrix not found in {self.result_dir}, returning empty DataFrame')
            return {
                'data': pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                'level': level,
            }
        if gene_ann_file is None:
            logging.warning(f'SCEVAN annotation not found in {self.result_dir}, returning empty DataFrame')
            return {
                'data': pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                'level': level,
            }

        adata = load_SCEVAN_expr(matrix_file, gene_ann_file, tn_annot_path)
        df_out = pseudo_bulk_adata(adata)
        save_tsv(df_out, self.save_dir / f'{self.eval_name}_{self.task_name}_{level}.tsv')

        return {
            'data': df_out,
            'level': level,
        }

    def extract_cnv_resolution(self, gene_annot_path=None, tn_annot_path=None):
        matrix_file, gene_ann_file = self._resolve_input_files()

        if matrix_file is None or gene_ann_file is None:
            logging.warning(f'SCEVAN inputs not found in {self.result_dir}, returning empty resolution DataFrame')
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        adata = load_SCEVAN_expr(matrix_file, gene_ann_file, tn_annot_path)
        df_resolution = _get_resolution_df(adata, cn_normal=0.01)
        save_tsv(df_resolution, self.save_dir / f'{self.eval_name}_{self.task_name}.tsv')
        return df_resolution
