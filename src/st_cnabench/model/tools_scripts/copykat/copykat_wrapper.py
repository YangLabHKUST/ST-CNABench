from pathlib import Path
from typing import Dict, Any, List
import pandas as pd
import yaml
import logging

from ...base import BaseModel
from ....modules.io_utils import get_per_dataset_param

class CopyKATModel(BaseModel):
    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="CopyKAT",
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
        final_res = output_dir / "copykat_object.rds"
        ## skip if output exists
        if not overwrite and final_res.exists():
            logging.warning(f"[{ds_id}] CopyKAT output exists. Skipping.")
            return {"_skip_execution": True}
        
        # Params
        genome = self.model_cfg.get("genome", "hg20")
        win_size = self.model_cfg.get("win_size", 25)
        n_cores = self.model_cfg.get("n_cores", 20)
        n_clones = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="n_clones",
            default_key="n_clones",
            default_value=2
        )
        # MTX
        matrix_dir = dataset_root / "filtered_feature_bc_matrix"
        if not matrix_dir.exists():
            raise FileNotFoundError(f"[{ds_id}] CopyKAT input dir not found: {matrix_dir}")
        # Annot
        if ref_norm is False:
            # create empty normal cells file
            logging.info(f"[{ds_id}] Model run in no reference normal mode.")
            return {
                "matrix_dir": matrix_dir,
                "normal_file": None,
                "genome": genome,
                "win_size": win_size,
                "n_cores": n_cores,
                "n_clones": n_clones
            }
        else: 
            normal_cells_path = output_dir / "normal_cells.txt"
            label_file = None
            for f in dataset_root.iterdir():
                if f.is_file() and "tumor_normal" in f.name:
                    label_file = f
                    break
            
            # Process Annot (extract normal cells)
            if label_file and label_file.exists():
                logging.info(f"[{ds_id}] Extracting normal cells from {label_file.name}...")
                df = pd.read_csv(label_file, sep='\t', comment='#', header=0, engine='python')
                
                df.columns = [c.strip() for c in df.columns]
                
                if "Barcode" not in df.columns or "tumor_normal" not in df.columns:
                    raise ValueError(f"Required columns 'Barcode' or 'tumor_normal' not found in {label_file.name}")

                normals = df[df["tumor_normal"].astype(str).str.lower() == "normal"]["Barcode"]
                
                normals.to_csv(normal_cells_path, index=False, header=False)
                logging.info(f"[{ds_id}] Saved {len(normals)} normal references.")
                    
            else:
                raise FileNotFoundError(f"[{ds_id}] Label file with tumor/normal info not found in {dataset_root}")

            return {
                "matrix_dir": matrix_dir,
                "normal_file": normal_cells_path,
                "genome": genome,
                "win_size": win_size,
                "n_cores": n_cores,
                "n_clones": n_clones
            }

    def build_command(
        self,
        dataset_cfg: Dict[str, Any],
        input_files: Dict[str, Path],
        output_dir: Path
    ) -> List[str]:
        
        if input_files.get("_skip_execution", False):
            return []
        ds_id = dataset_cfg.get("dataset_id", "unknown_sample")

        # src/model/tools_scripts/copykat/run_copykat.R
        script_path = self.get_script_path(
            script_name="run_copykat.R", 
            subfolder="copykat"
        )
        
        cmd = [
            "Rscript",
            str(script_path),
            str(input_files["matrix_dir"]),
            str(output_dir),
            str(input_files["normal_file"]),
            str(input_files["genome"]),
            str(input_files["win_size"]),
            str(input_files["n_cores"]),
            str(input_files["n_clones"]),
            str(ds_id)
        ]
        
        return cmd
