import pandas as pd
import anndata as ad
import logging
from .base import BaseLoader
from ..utils.io import load_InferCNV_expr, pseudo_bulk_adata, _get_resolution_df, save_tsv


class InferCNVExprLoader(BaseLoader):
    def __init__(self, task_name: str, eval_name: str, result_dir: str, save_dir: str):
        model_name = 'InferCNV'
        super().__init__(model_name=model_name, eval_name=eval_name, task_name=task_name, result_dir=result_dir, save_dir=save_dir)

    def extract_cna_profile(self, gene_annot_path, tn_annot_path):
        level = 'Log2_Copy_Ratio'
        matrix_file = self.result_dir / 'infercnv.observations.txt'
        if not matrix_file.exists():
            logging.warning(f'No InferCNV observation file found in {self.result_dir}, returning empty DataFrame')
            return {
                'data': pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                'level': level,
            }

        adata = load_InferCNV_expr(matrix_file, gene_annot_path, tn_annot_path)
        df_out = pseudo_bulk_adata(adata)
        save_tsv(df_out, self.save_dir / f'{self.eval_name}_{self.task_name}_{level}.tsv')

        return {
            'data': df_out,
            'level': level,
        }

    def extract_cna_resolution(self, gene_annot_path, tn_annot_path):
        matrix_file = self.result_dir / 'infercnv.observations.txt'
        if not matrix_file.exists():
            logging.warning(f'No InferCNV observation file found in {self.result_dir}, returning empty resolution DataFrame')
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])

        adata = load_InferCNV_expr(matrix_file, gene_annot_path, tn_annot_path)
        df_resolution = _get_resolution_df(adata, cn_normal=0.01)
        save_tsv(df_resolution, self.save_dir / f'{self.eval_name}_{self.task_name}.tsv')
        return df_resolution

    def extract_tumor_normal_preds(self):
        from ..utils.find_normal import find_normal_tumor

        obs_mtx_file = self.result_dir / 'infercnv.observations.txt'
        if not obs_mtx_file.exists():
            raise FileNotFoundError(f'No InferCNV observation file found in {self.result_dir}')

        df_obs = pd.read_csv(obs_mtx_file, sep=r'\s+', index_col=0)
        preds_series = find_normal_tumor(None, df_obs, mode='no_ref', baseline=1.0)
        df_preds = preds_series.reset_index()
        df_preds.columns = ['Barcodes', 'Preds']
        df_preds['Preds'] = df_preds['Preds'].map({'normal': 0, 'tumor': 1})
        save_tsv(df_preds, self.save_dir / f'{self.eval_name}_{self.task_name}.tsv')
        return df_preds

    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        from ..utils.find_normal import find_normal_tumor

        ref_mtx_file = self.result_dir / 'infercnv.references.txt'
        if not ref_mtx_file.exists():
            raise FileNotFoundError(f'No InferCNV reference file found in {self.result_dir}')

        obs_mtx_file = self.result_dir / 'infercnv.observations.txt'
        if not obs_mtx_file.exists():
            raise FileNotFoundError(f'No InferCNV observation file found in {self.result_dir}')

        df_ref = pd.read_csv(ref_mtx_file, sep=r'\s+', index_col=0)
        df_obs = pd.read_csv(obs_mtx_file, sep=r'\s+', index_col=0)
        preds_series = find_normal_tumor(df_ref, df_obs, mode='ref')
        df_preds = preds_series.reset_index()
        df_preds.columns = ['Barcodes', 'Preds']
        df_preds['Preds'] = df_preds['Preds'].map({'normal': 0, 'tumor': 1})
        save_tsv(df_preds, self.save_dir / f'{self.eval_name}_{self.task_name}.tsv')
        return df_preds
