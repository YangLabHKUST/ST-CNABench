from abc import ABC, abstractmethod
from pathlib import Path
import logging
import anndata as ad
import re

class BaseLoader(ABC):
    def __init__(self, model_name: str, eval_name:str, task_name:str, result_dir: str, save_dir: str):
        self.model_name = model_name
        self.eval_name = eval_name
        self.task_name = task_name
        self.result_dir = Path(result_dir)
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def extract_cna_profile(self, **kwargs):
        """
        CNA Profile Extraction (Unified Zero-Centered):
        Returns:
            pd.DataFrame with columns:
            ['Chr', 'Start', 'End', 'Bin_ID', 'TCN_Score', 'LOH_Status']
        
        TCN_Score is a zero-centered continuous score (gain > 0, neutral ~ 0, loss < 0).
        LOH_Status is 1 (LOH), 0 (non-LOH), or NaN when unavailable.
        """
        pass

    def extract_clone_cna_profile(self, **kwargs):
        """
        Clone-Specific (Pseudo-bulk) CNA Profile Extraction:
        Extract the shared CNA profile for each predicted subclone.
        
        Returns:
            pd.DataFrame with columns:
            ['Clone_ID', 'Chr', 'Start', 'End', 'TCN_Score']
            Clone_ID must match the labels returned by extract_subcluster_preds
            (for example, infercnv_1).
            TCN_Score can be a zero-centered continuous score or an absolute copy number.
        """
        pass

    def extract_cna_resolution(self, **kwargs):
        """
        CNA Resolution Extraction:
        {
            "data": pd.DataFrame, # with columns: ['Chromosome', 'Start', 'End', 'Length']
        }
        """
        pass

    def extract_tumor_normal_preds(self, **kwargs):
        """
        Tumor Normal Prediction Extraction:
        {
            "data": pd.DataFrame, # with columns: ['Barcodes', 'Preds']
        }
        Values Mapping:
        mapping = {
            'tumor' , 'Tumor', 'malignant': 1
            'normal', 'Normal': 0
            'not.defined', 'Unassigned', 'unknown': -1
        }
        """
        pass

    def extract_subcluster_preds(self):
        """
        Subcluster Prediction Extraction:
        {
            "data": pd.DataFrame, # with columns: ['Barcodes', 'Label_preds']
        }
        Subcluster labels transformed to format: {MODELNAME}_{i}
        """
        pass
    def extract_computational_efficiency_data(self):
        """
        Global Computational Efficiency Extraction:
        {
            "model": str,
            "runtime_sec": float,
            "mem_gb": float,
            "cpu_percent": float,
            "exit_status": int
        }
        """
        # Load from .perf file
        perf_files = list(self.result_dir.glob("*.perf"))
        if not perf_files:
            logging.warning(f"[{self.model_name}] No .perf file found in {self.result_dir}")
            return None

        perf_path = perf_files[0]
        with open(perf_path, 'r') as f:
            content = f.read()

        patterns = {
            'wall_time_str': r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)",
            'max_rss_kb': r"Maximum resident set size \(kbytes\): (\d+)",
            'cpu_percent': r"Percent of CPU this job got: (\d+)%",
            'exit_status': r"Exit status: (\d+)"
        }

        stats = {'model': self.model_name}
        for key, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                stats[key] = match.group(1)

        if 'wall_time_str' in stats:
            parts = list(map(float, stats['wall_time_str'].split(':')))
            if len(parts) == 3: # h:mm:ss
                stats['runtime_sec'] = parts[0] * 3600 + parts[1] * 60 + parts[2]
            else: # m:ss
                stats['runtime_sec'] = parts[0] * 60 + parts[1]
        
        if 'max_rss_kb' in stats:
            stats['mem_gb'] = round(int(stats['max_rss_kb']) / (1024 * 1024), 2)

        return stats
