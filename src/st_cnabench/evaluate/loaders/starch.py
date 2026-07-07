import pandas as pd
import anndata as ad
import numpy as np
import logging
from .base import BaseLoader
from ..utils.io import load_csv, load_tsv, load_STARCH, pseudo_bulk_adata,pseudo_clonal_bulk, save_tsv, _get_resolution_df, save_adata

class STARCHLoader(BaseLoader):
    def __init__(self, task_name:str, eval_name:str, result_dir: str, save_dir: str):
        model_name = 'STARCH'
        super().__init__(model_name=model_name,eval_name=eval_name,task_name=task_name,result_dir=result_dir, save_dir=save_dir)
        
    def extract_cna_profile(self, gene_annot_path, tn_annot_path):
        level = 'CNA_Event'
        cna_files = list(self.result_dir.glob("states_*.csv"))
        if not cna_files:
            logging.warning(f"No STARCH CNA state file found in {self.result_dir} matching states_*.csv, returning empty DataFrame")
            return {
                "data": pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                "level": level
             }
        cna_file = cna_files[-1]
        clone_files = list(self.result_dir.glob("labels_*.csv"))
        if not clone_files:
            logging.warning(f"No STARCH clone file found in {self.result_dir} matching labels_*.csv, returning empty DataFrame")
            return {
                "data": pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                "level": level
            }
            
        clone_file = clone_files[-1]
        processed_input_dir = self.result_dir / "processed_input"
        name_mapping_files = list(processed_input_dir.glob("*_spot_mapping.tsv"))
        if not name_mapping_files:
            logging.warning(f"No STARCH name mapping file found in {processed_input_dir} matching name_mapping_*.tsv, returning empty DataFrame")
            return {
                "data": pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'ID', 'CN_Score', 'LOH_Status']),
                "level": level
            }
        name_mapping_file = name_mapping_files[-1]
        adata = load_STARCH(cna_file, clone_file, gene_annot_path, tn_annot_path, name_mapping_file)
        df_out = pseudo_bulk_adata(adata)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_{level}.tsv")

        return {
            'data': df_out,
            'level': level
        }
    
    def extract_clone_cna_profile(self, gene_annot_path, tn_annot_path):
        df_subcluster = self.extract_subcluster_preds()
        cna_files = list(self.result_dir.glob("states_*.csv"))
        if not cna_files:
            logging.warning(f"No STARCH CNA state file found in {self.result_dir} matching states_*.csv, returning empty DataFrame")
            return pd.DataFrame()
        cna_file = cna_files[-1]
        clone_files = list(self.result_dir.glob("labels_*.csv"))
        if not clone_files:
            logging.warning(f"No STARCH clone file found in {self.result_dir} matching labels_*.csv, returning empty DataFrame")
            return pd.DataFrame()
        clone_file = clone_files[-1]
        processed_input_dir = self.result_dir / "processed_input"
        name_mapping_files = list(processed_input_dir.glob("*_spot_mapping.tsv"))
        if not name_mapping_files:
            logging.warning(f"No STARCH name mapping file found in {processed_input_dir} matching name_mapping_*.tsv, returning empty DataFrame")
            return pd.DataFrame()
        name_mapping_file = name_mapping_files[-1]
        adata = load_STARCH(cna_file, clone_file, gene_annot_path, tn_annot_path, name_mapping_file)
        df_subcluster = df_subcluster.set_index('Barcodes')
        common_cells = adata.obs_names.intersection(df_subcluster.index)
        if len(common_cells) == 0:
            logging.warning(f"[{self.model_name}] No overlapping barcodes between CNA matrix and subcluster predictions, returning empty DataFrame.")
            return pd.DataFrame()
        adata = adata[common_cells].copy()
        adata.obs['Clone_ID'] = df_subcluster.loc[common_cells, 'Label_preds']
        df_out = pseudo_clonal_bulk(adata, clone_col='Clone_ID')
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        logging.info(f"[{self.model_name}] Successfully extracted clone CNA profiles.")
        return df_out
    def extract_cna_resolution(self, gene_annot_path, tn_annot_path):
        cna_files = list(self.result_dir.glob("states_*.csv"))
        if not cna_files:
            logging.warning(
                f"No STARCH CNA state file found in {self.result_dir} matching states_*.csv, returning empty resolution DataFrame"
            )
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])
        cna_file = cna_files[-1]
        clone_files = list(self.result_dir.glob("labels_*.csv"))
        if not clone_files:
            logging.warning(
                f"No STARCH clone file found in {self.result_dir} matching labels_*.csv, returning empty resolution DataFrame"
            )
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])
        clone_file = clone_files[-1]
        processed_input_dir = self.result_dir / "processed_input"
        name_mapping_files = list(processed_input_dir.glob("*_spot_mapping.tsv"))
        if not name_mapping_files:
            logging.warning(
                f"No STARCH name mapping file found in {processed_input_dir} matching name_mapping_*.tsv, returning empty resolution DataFrame"
            )
            return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Length'])
        name_mapping_file = name_mapping_files[-1]
        adata = load_STARCH(cna_file, clone_file, gene_annot_path, tn_annot_path, name_mapping_file)
        df_resolution = _get_resolution_df(adata, cn_normal=0)
        save_tsv(df_resolution, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_resolution

    def extract_tumor_normal_preds(self):
        pattern = "labels_*.csv"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            raise FileNotFoundError(f"STARCH tumor normal prediction file not found in {self.result_dir}")
        preds_file = matches[-1]
        processed_input_dir = self.result_dir / "processed_input"
        name_mapping_files = list(processed_input_dir.glob("*_spot_mapping.tsv"))
        if not name_mapping_files:
            raise FileNotFoundError(f"No STARCH name mapping file found in {processed_input_dir} matching name_mapping_*.tsv")
        name_mapping_file = name_mapping_files[-1]
        name_mapping = load_tsv(name_mapping_file, header=True).reset_index()
        df_label = load_csv(preds_file, header = True).reset_index()
        df_label = df_label.rename(columns={
            'index': 'barcodes',
            '0': 'clone_id'         
        })
        df_label['barcodes'] = df_label['barcodes'].apply(
        lambda x: 'x'.join([str(int(float(i))) for i in x.split('x')])
        )
        
        df_label = df_label.merge(name_mapping, left_on='barcodes', right_on='spot_name', how='left')
        df = df_label[['barcode', 'clone_id']].copy()
        df.columns = ['Barcodes', 'Preds']
        # The last cluster id is normal see: https://github.com/raphael-group/STARCH/issues/8
        df['Preds'] = df['Preds'].astype(int)
        normal_id = df['Preds'].max()
        df['Preds'] = df['Preds'].apply(lambda x: 0 if x == normal_id else 1)
        save_tsv(df, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")

        return df
    
    # Copy from extract_tumor_normal_preds and remove ref normal cells
    def extract_sub_normal_tumor_normal_preds(self, tn_annot_path_model_run):
        pattern = "labels_*.csv"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            raise FileNotFoundError(f"STARCH tumor normal prediction file not found in {self.result_dir}")
        preds_file = matches[-1]
        cna_files = list(self.result_dir.glob("states_*.csv"))
        if not cna_files:
            raise FileNotFoundError(f"No STARCH CNA state file found in {self.result_dir} matching states_*.csv")
        cna_file = cna_files[-1]
        processed_input_dir = self.result_dir / "processed_input"
        name_mapping_files = list(processed_input_dir.glob("*_spot_mapping.tsv"))
        if not name_mapping_files:
            raise FileNotFoundError(f"No STARCH name mapping file found in {processed_input_dir} matching name_mapping_*.tsv")
        name_mapping_file = name_mapping_files[-1]
        name_mapping = load_tsv(name_mapping_file, header=True).reset_index()
        cna_df = load_csv(cna_file, header=True)

        df_label = load_csv(preds_file, header = True).reset_index()
        df_label = df_label.rename(columns={
            'index': 'barcodes',
            '0': 'clone_id'         
        })
        df_label['barcodes'] = df_label['barcodes'].apply(
        lambda x: 'x'.join([str(int(float(i))) for i in x.split('x')])
        )
        
        df_label = df_label.merge(name_mapping, left_on='barcodes', right_on='spot_name', how='left')
        df = df_label[['barcode', 'clone_id']].copy()
        df.columns = ['Barcodes', 'Preds']
        # The last cluster id is normal see: https://github.com/raphael-group/STARCH/issues/8
        df['Preds'] = df['Preds'].astype(int)
        # Notes: for the CNA state table, the last column is the reference normal clone.
        normal_id_ref = df['Preds'].max()
        cna_df = cna_df.drop(columns=[cna_df.columns[-1]])
        std_per_clone = cna_df.std()
        inferred_normal_id = std_per_clone.idxmin()
        logging.info(f"STARCH inferred normal clone id: {inferred_normal_id}, ref normal clone id: {normal_id_ref}")
        normal_id = set([int(normal_id_ref), int(inferred_normal_id)])
        df['Preds'] = df['Preds'].apply(lambda x: 0 if x in normal_id else 1)
        # load tn_annot_path_model_run to remove ref normal cells
        df_model_run = pd.read_csv(tn_annot_path_model_run, sep='\t', header=0, comment='#')
        if df_model_run.shape[1] == 2:
            df_model_run.columns = ['Barcodes', 'tumor_normal']
        elif df_model_run.shape[1] == 4:
            df_model_run.columns = ['Barcodes', 'tumor_normal', 'leiden', 'celltype']
        normal_barcode = df_model_run[(df_model_run['tumor_normal'] == 'normal') | (df_model_run['tumor_normal'] == 'Normal')]['Barcodes'].tolist()
        df_out = df[~df['Barcodes'].isin(normal_barcode)]

        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_out
    
    def extract_subcluster_preds(self):
        pattern = "labels_*.csv"
        matches = list(self.result_dir.glob(f"**/{pattern}"))
        if not matches:
            raise FileNotFoundError(f"STARCH tumor normal prediction file not found in {self.result_dir}")
        subcluster_file = matches[-1]
        df_label = load_csv(subcluster_file, header = True).reset_index()
        processed_input_dir = self.result_dir / "processed_input"
        name_mapping_files = list(processed_input_dir.glob("*_spot_mapping.tsv"))
        if not name_mapping_files:
            raise FileNotFoundError(f"No STARCH name mapping file found in {processed_input_dir} matching name_mapping_*.tsv")
        name_mapping_file = name_mapping_files[-1]
        name_mapping = load_tsv(name_mapping_file, header=True).reset_index()
        df_label = df_label.rename(columns={
            'index': 'barcodes',
            '0': 'clone_id'         
        })
        df_label['barcodes'] = df_label['barcodes'].apply(
        lambda x: 'x'.join([str(int(float(i))) for i in x.split('x')])
        )
        df_label = df_label.merge(name_mapping, left_on='barcodes', right_on='spot_name', how='left')
        df_subcluster = df_label[['barcode', 'clone_id']].copy()
        df_subcluster.columns = ['Barcodes', 'Label_preds']
        # Transform subcluster label to format: {MODELNAME}_{i}
        df_subcluster['Label_preds'] = df_subcluster['Label_preds'].apply(lambda x: f"{self.model_name}_{x}")
        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster
