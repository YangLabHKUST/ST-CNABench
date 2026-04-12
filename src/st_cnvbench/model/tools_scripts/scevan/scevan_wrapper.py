from pathlib import Path
from typing import Dict, Any, List
import yaml
import logging

from ...base import BaseModel

class SCEVANModel(BaseModel):
    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="SCEVAN",
            exec_mode=exec_mode
        )

    def prepare_inputs(self, dataset_cfg, output_dir, overwrite=False) -> Dict[str, Any]:
        ds_id = dataset_cfg.get("dataset_id")
        dataset_root = Path(dataset_cfg["output"]["root"]).resolve()
        ref_norm = dataset_cfg.get("ref_norm")
        # Overwrite check
        final_res = output_dir / "output"
        if not overwrite and final_res.exists():
            logging.warning(f"[{ds_id}] SCEVAN output exists. Skipping.")
            return {"_skip_execution": True}
        # MTX
        matrix_dir = dataset_root / "filtered_feature_bc_matrix"
        if not (matrix_dir / "matrix.mtx.gz").exists():
             raise FileNotFoundError(f"[{ds_id}] matrix.mtx.gz not found in {matrix_dir}")

        # Annot
        if ref_norm is False:
            # Without ref normal model
            logging.info(f"[{ds_id}] Model run in no reference normal mode.")
            return {
                "matrix_dir": matrix_dir,
                "meta_file": None
            }
        else:
            # With ref normal model
            label_file = None
            for f in dataset_root.iterdir():
                if "tumor_normal" in f.name:
                    label_file = f
                    break
            if not label_file:
                raise FileNotFoundError(f"[{ds_id}] Annotation file not found in {dataset_root}")

            return {
                "matrix_dir": matrix_dir,
                "meta_file": label_file
            }

    def build_command(self, dataset_cfg, input_files, output_dir) -> List[str]:
        
        if input_files.get("_skip_execution", False):
            return []
        
        # src/model/tools_scripts/scevan/run_scevan.R
        script_path = self.get_script_path("run_scevan.R", subfolder="scevan")
        ds_id = dataset_cfg.get("dataset_id")
        
        cmd = [
            "Rscript",
            str(script_path),
            str(input_files["matrix_dir"]),
            str(output_dir),
            str(input_files["meta_file"]),
            str(ds_id),
            str(self.model_cfg.get("n_threads", 10))
        ]
        return cmd