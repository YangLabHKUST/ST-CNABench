from pathlib import Path
from typing import Dict, Any, List
import yaml
import logging
from ...base import BaseModel
from ....modules.io_utils import get_per_dataset_param
class ClonalscopeNoWGSModel(BaseModel):
    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="Clonalscope_NoWGS",
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
        final_res = output_dir / "Cov_obj_second.rds"
        if not overwrite and final_res.exists():
            logging.warning(f"[{ds_id}] Clonalscope_NoWGS output exists. Skipping.")
            return {"_skip_execution": True}
        
        # MTX
        matrix_dir = dataset_root / "filtered_feature_bc_matrix"
        if not matrix_dir.exists():
            raise FileNotFoundError(f"[{ds_id}] Matrix dir not found: {matrix_dir}")

        # Annot
        if ref_norm is True:
            label_file = None
            if dataset_root.exists():
                for f in dataset_root.iterdir():
                    if "tumor_normal" in f.name:
                        label_file = f
                        break
            
            if not label_file:
                raise FileNotFoundError(f"[{ds_id}] Tumor/Normal annotation file not found in {dataset_root}")
        else:
            logging.info(f"[{ds_id}] Clonalscope_NoWGS running in No reference normal cell mode")
            label_file = "None"
        # Ref data
        gene_coords_file = self.model_cfg.get("gene_coords_file")
        aux_data_dir = self.model_cfg.get("aux_data_dir")
        
        if not gene_coords_file or not aux_data_dir:
            raise ValueError("[Clonalscope_NoWGS] 'gene_coords_file' or 'aux_data_dir' not set in models.yaml")

        mincell = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="mincell",
            default_key="mincell",
            default_value=50
        )
        return {
            "matrix_dir": matrix_dir,
            "meta_file": label_file,
            "gene_coords_file": Path(gene_coords_file).resolve(),
            "aux_data_dir": Path(aux_data_dir).resolve(),
            "mincell": mincell
        }

    def build_command(
        self,
        dataset_cfg: Dict[str, Any],
        input_files: Dict[str, Path],
        output_dir: Path
    ) -> List[str]:
        
        if input_files.get("_skip_execution", False):
            return []
        # src/model/tools_scripts/clonalscope_nowgs/run_clonalscope_nowgs.R
        script_path = self.get_script_path(
            script_name="run_clonalscope_nowgs.R",
            subfolder="clonalscope_nowgs"
        )
        
        cmd = [
            "Rscript",
            str(script_path),
            str(input_files["matrix_dir"]),
            str(output_dir),
            str(input_files["meta_file"]),
            str(input_files["gene_coords_file"]),
            str(input_files["aux_data_dir"]),
            str(input_files["mincell"])
        ]
        
        return cmd