# src/model/tools_scripts/infercnv/infercnv_model.py

from pathlib import Path
from typing import Dict, Any, List
import pandas as pd
import yaml
import logging

from ...base import BaseModel
from ....modules.io_utils import get_per_dataset_param
class InferCNVModel(BaseModel):
    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="InferCNV",
            exec_mode=exec_mode
        )

    def prepare_inputs(
        self,
        dataset_cfg: Dict[str, Any],
        output_dir: Path,
        overwrite: bool = False
    ) -> Dict[str, Path]:
        
        ds_id = dataset_cfg.get("dataset_id")
        dataset_root = Path(dataset_cfg["output"]["root"]).resolve()
        ref_norm = dataset_cfg.get("ref_norm")
        # Overwrite check
        final_res = output_dir / "infercnv_obj.rdata"
        if not overwrite and final_res.exists():
            logging.warning(f"[{ds_id}] InferCNV output exists. Skipping.")
            return {"_skip_execution": True}
        
        # MTX
        matrix_path = dataset_root / "filtered_feature_bc_matrix"
        if not matrix_path.exists():
            raise FileNotFoundError(f"[{ds_id}] 10x Matrix dir not found: {matrix_path}")
        
        annot_out_path = output_dir / "infercnv_annotations.txt"

        # Annot
        if ref_norm is True:
            label_file = dataset_root / f"metadata_{ds_id}_tumor_normal.tsv"

            if not label_file.exists():
                 raise FileNotFoundError(f"[{ds_id}] Metadata file not found: {label_file}")

            logging.info(f"[{ds_id}] Processing annotations from {label_file.name}...")

            df = pd.read_csv(label_file, sep='\t', comment='#')
            annot_df = df[['Barcode', 'tumor_normal']].copy()
            annot_df.columns = ['cell_id', 'group']
            annot_df['group'] = annot_df['group'].apply(
                lambda x: "Normal" if str(x).lower() == "normal" else "Tumor"
            )
        else:
            annot_out_path = output_dir / "infercnv_observation_groups.txt"
            barcodes_file = matrix_path / "barcodes.tsv.gz"
            if not barcodes_file.exists():
                raise FileNotFoundError(f"[{ds_id}] Barcodes file not found: {barcodes_file}")

            logging.info(f"[{ds_id}] Creating technical observation groups for no-reference InferCNV run.")
            barcodes = pd.read_csv(
                barcodes_file,
                sep="\t",
                header=None,
                compression="gzip",
            ).iloc[:, 0].astype(str)
            annot_df = pd.DataFrame({"cell_id": barcodes, "group": "Observation"})
        annot_df.to_csv(annot_out_path, sep='\t', index=False, header=False)
        logging.info(f"[{ds_id}] Saved {len(annot_df)} cell annotations.")

        # Ref data
        gene_order_path = self.model_cfg.get("gene_order_file")
        if not gene_order_path:
            raise ValueError("[InferCNV] 'gene_order_file' not specified in models.yaml")

        # Dataset-specific config
        k_obs_groups = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="k_obs_groups",
            default_key="k_obs_groups",
            default_value=2
        )
        return {
            "counts_path": matrix_path,
            "annotation_file": annot_out_path,
            "gene_order_file": Path(gene_order_path),
            "ref_norm": ref_norm,
            "k_obs_groups": k_obs_groups
        }

    def build_command(
        self,
        dataset_cfg: Dict[str, Any],
        input_files: Dict[str, Path],
        output_dir: Path
    ) -> List[str]:
        
        if input_files.get("_skip_execution", False):
            return []

        # src/model/tools_scripts/infercnv/run_infercnv.R
        script_path = self.get_script_path(
            script_name="run_infercnv.R", 
            subfolder="infercnv"
        )
        
        cmd = [
            "Rscript",
            str(script_path),
            str(input_files["counts_path"]),
            str(input_files["annotation_file"]),
            str(input_files["gene_order_file"]),
            str(input_files["ref_norm"]),
            str(output_dir),
            str(self.model_cfg.get("cutoff", 0.1)),
            str(self.model_cfg.get("n_threads", 4)),
            str(input_files["k_obs_groups"])
        ]
        
        return cmd
