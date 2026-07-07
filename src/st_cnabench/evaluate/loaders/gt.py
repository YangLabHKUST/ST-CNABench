import pandas as pd
import numpy as np
from .base import BaseLoader
from ..utils.io import load_facets, load_tsv, save_tsv
from ..utils.constants import HG38_INFO
from pathlib import Path
import logging
import os 
import glob
class GTLoader(BaseLoader):
    def __init__(self, task_name: str, eval_name:str, result_dir: str, save_dir: str):
        super().__init__(
            model_name='GT', 
            eval_name=eval_name, 
            task_name=task_name, 
            result_dir=result_dir, 
            save_dir=save_dir
        )
    def extract_cna_profile(self, gt_path, FOCAL_TYPE1_THRESHOLD=3e6):
        
        if not Path(gt_path).exists():
            raise FileNotFoundError(f"GT result file not found at {gt_path}")
        df_out = load_facets(gt_path)
        df_out['Length'] = df_out['End'] - df_out['Start']
        
        df_out['Is_Focal_Type1'] = df_out['Length'] < FOCAL_TYPE1_THRESHOLD
        def check_focal_type2(row):
            chrom = str(row['Chromosome']).replace('chr', '')
            if chrom not in HG38_INFO:
                return False 
            cen = HG38_INFO[chrom]['cen']
            chrom_len = HG38_INFO[chrom]['len']
            p_arm_len = cen
            q_arm_len = chrom_len - cen
        
            midpoint = (row['Start'] + row['End']) / 2
            arm_len = p_arm_len if midpoint < cen else q_arm_len
            return row['Length'] < (0.25 * arm_len)
            
        df_out['Is_Focal_Type2'] = df_out.apply(check_focal_type2, axis=1)
        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_GT.tsv")
        return df_out

    def extract_clone_cna_profile(self, cna_gt_path=None):
        search_pattern = os.path.join(cna_gt_path, "*_GT_Profile.txt")
        files = glob.glob(search_pattern)
        
        if not files:
            logging.warning(f"No *_GT_Profile.txt found in {cna_gt_path}")
            return None

        all_profiles = []
        for file_path in files:
            file_name = os.path.basename(file_path)
            clone_id = file_name.replace("_GT_Profile.txt", "")
            
            if clone_id.lower() in ["all_tumors", "normal"]:
                continue
                
            df = pd.read_csv(file_path, sep='\t')
            if 'CN_Score_Continuous' in df.columns:
                df['CN_Score'] = df['CN_Score_Continuous']
            
            df['Clone_ID'] = clone_id
            all_profiles.append(df)

        if not all_profiles:
            logging.warning("No valid subclone GT profiles found after filtering All_Tumors.")
            return None

        gt_profile = pd.concat(all_profiles, ignore_index=True)
        
        # SlideDNAseq CNA profiles use the normal label as the reference
        # baseline, so the normal-cluster profile is defined as all-zero.
        # Reuse the first clone's genomic bin layout and add that explicit
        # zero-valued "Normal" profile for downstream clone-level evaluation.
        first_clone = gt_profile['Clone_ID'].iloc[0]
        df_normal = gt_profile[gt_profile['Clone_ID'] == first_clone].copy()
        
        df_normal['Clone_ID'] = 'Normal'
        df_normal['CN_Score'] = 0.0
        
        gt_profile = pd.concat([gt_profile, df_normal], ignore_index=True)

        logging.info(f"Loaded GT CNA profiles for {len(gt_profile['Clone_ID'].unique())} clones: {gt_profile['Clone_ID'].unique().tolist()}")
        return gt_profile
    
    def extract_tumor_normal_preds(self, tn_annot_path=None, tn_annot_path_model_run=None):
        df = pd.read_csv(tn_annot_path, sep='\t', header=0, comment='#')
        
        # Label file type
        if df.shape[1] == 2:
            df.columns = ['Barcodes', 'Preds']
        elif df.shape[1] == 4:
            df.columns = ['Barcodes', 'Preds', 'leiden', 'celltype']
        mapping = {
            'tumor': 1,
            'normal': 0
        }
        df['Preds'] = df['Preds'].map(mapping)
        if tn_annot_path_model_run is not None:
            # Subset normal mode, remove ref normal cells
            df_model_run = pd.read_csv(tn_annot_path_model_run, sep='\t', header=0, comment='#')
            if df_model_run.shape[1] == 2:
                df_model_run.columns = ['Barcodes', 'Preds']
            elif df_model_run.shape[1] == 4:
                df_model_run.columns = ['Barcodes', 'Preds', 'leiden', 'celltype']
            df = df[df['Barcodes'].isin(df_model_run['Barcodes'])]
        save_tsv(df, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df
    
    def extract_subcluster_preds(self, subcluster_annot_path=None):
        df_subcluster = pd.read_csv(subcluster_annot_path, sep='\t', header=0)
        df_subcluster.columns = ['Barcodes', 'Label_preds']
        save_tsv(df_subcluster, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_subcluster
''' old version read gatk4 seg file
def extract_cna_profile(self, gt_path: str, threshold: int=20):
    self.gt_path = gt_path
    self.threshold = threshold
    df_gt = load_gatk_seg(gt_path)

    df_gt = df_gt[df_gt['NUM_POINTS_COPY_RATIO'] >= self.threshold].copy()

    values = np.power(2, df_gt['MEAN_LOG2_COPY_RATIO']) * 2
    df_out = pd.DataFrame({
        'Chr': df_gt['CONTIG'],
        'Start': df_gt['START'],
        'End': df_gt['END'],
        'Log2_Copy_Ratio': values
    })

    out_path = self.save_dir / f"{self.eval_name}_{self.task_name}GT.tsv"
    save_tsv(df_out, out_path)

    return df_out
'''
## Return a bundle dict including log2 values Integer CN and CN Event
