import logging
import pandas as pd

from .base import BaseLoader
from ..utils.io import save_tsv


class SlideCNALoader(BaseLoader):
    """Loader for SlideCNA subclone labels."""

    def __init__(self, task_name: str, eval_name: str, result_dir: str, save_dir: str):
        model_name = "SlideCNA"
        super().__init__(
            model_name=model_name,
            eval_name=eval_name,
            task_name=task_name,
            result_dir=result_dir,
            save_dir=save_dir,
        )

    def _get_exported_label_file(self):
        """Return the exported malignant barcode label file path."""
        label_file = self.result_dir / "slidecna_malig_barcode_labels.tsv.gz"
        if not label_file.exists():
            raise FileNotFoundError(f"No exported SlideCNA malignant label file found in {self.result_dir}")
        return label_file

    def _get_exported_profile_file(self):
        """Return the exported malignant clone-profile file path."""
        profile_file = self.result_dir / "slidecna_malig_clone_profiles.tsv.gz"
        if not profile_file.exists():
            raise FileNotFoundError(f"No exported SlideCNA malignant clone profile file found in {self.result_dir}")
        return profile_file

    def extract_cna_profile(self, **kwargs):
        """Expose misuse clearly because SlideCNA is only enabled for in-slice subclone labels."""
        raise NotImplementedError("SlideCNA loader is only configured for subclone_detection_in_slice.")

    def extract_clone_cna_profile(self, **kwargs):
        """Load the exported SlideCNA malignant clone CNA profiles."""
        profile_file = self._get_exported_profile_file()
        df_out = pd.read_csv(profile_file, sep="\t", compression="gzip")
        required_cols = {"Clone_ID", "Chromosome", "Start", "End", "ID", "CN_Score", "LOH_Status"}
        missing_cols = required_cols - set(df_out.columns)
        if missing_cols:
            raise ValueError(f"Exported SlideCNA clone profile file missing columns: {sorted(missing_cols)}")

        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}_clone_profiles.tsv")
        logging.info(f"[{self.model_name}] Loaded exported malignant-only clone CNA profiles.")
        return df_out

    def extract_tumor_normal_preds(self, **kwargs):
        """Expose misuse clearly because SlideCNA does not provide independent tumor-normal calls here."""
        raise NotImplementedError("SlideCNA tumor-normal evaluation is disabled for this benchmark scope.")

    def extract_subcluster_preds(self):
        """Load the exported SlideCNA malignant barcode labels."""
        label_file = self._get_exported_label_file()
        df_out = pd.read_csv(label_file, sep="\t", compression="gzip")
        required_cols = {"Barcodes", "Label_preds"}
        missing_cols = required_cols - set(df_out.columns)
        if missing_cols:
            raise ValueError(f"Exported SlideCNA label file missing columns: {sorted(missing_cols)}")
        df_out = df_out[["Barcodes", "Label_preds"]].drop_duplicates()

        save_tsv(df_out, self.save_dir / f"{self.eval_name}_{self.task_name}.tsv")
        return df_out
