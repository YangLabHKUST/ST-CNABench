from pathlib import Path
from typing import Dict, Any, List
import yaml
import logging
from ...base import BaseModel
from ....modules.io_utils import get_per_dataset_param

class ClonalscopeWGSModel(BaseModel):
    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="Clonalscope_WGS",
            exec_mode=exec_mode
        )

    def prepare_inputs(self, dataset_cfg, output_dir, overwrite=False) -> Dict[str, Path]:
        ds_id = dataset_cfg.get("dataset_id")
        dataset_root = Path(dataset_cfg["output"]["root"]).resolve()
        ref_norm = dataset_cfg.get("ref_norm")
        # Overwrite check
        final_res = output_dir / "clonalscope_wgs_obj.rds"
        if not overwrite and final_res.exists():
            logging.warning(f"[{ds_id}] Clonalscope_WGS output exists. Skipping.")
            return {"_skip_execution": True}
        
        # MTX
        matrix_dir = dataset_root / "filtered_feature_bc_matrix"
        if not matrix_dir.exists():
             raise FileNotFoundError(f"[{ds_id}] Matrix dir not found: {matrix_dir}")

        # Annot
        if ref_norm is True:
            label_file = None
            for f in dataset_root.iterdir():
                if "tumor_normal" in f.name:
                    label_file = f
                    break
            if not label_file:
                raise FileNotFoundError(f"[{ds_id}] Annotation file not found in {dataset_root}")
        else:
            logging.info(f"[{ds_id}] Clonalscope_WGS running in No reference normal cell mode")
            label_file = "None"

        # WGS/WES file
        wgs_tumor = dataset_root / "wgs_wes_tumor_data.bedg" 
        wgs_normal = dataset_root / "wgs_wes_normal_data.bedg"

        if not wgs_tumor.exists() or not wgs_normal.exists():
            raise FileNotFoundError(
                f"[{ds_id}] [Clonalscope_WGS] requires WGS bedgraph files.\n"
                f"    Expected: {wgs_tumor} AND {wgs_normal}\n"
                "    Please ensure 'wgs_wes_tumor_bedg' and 'wgs_wes_normal_bedg' are defined in data.config.yaml"
            )

        # Ref data
        aux_dir = self.model_cfg.get("aux_data_dir")
        gene_coords = self.model_cfg.get("gene_coords_file")
        
        if not aux_dir or not Path(aux_dir).exists():
            raise ValueError(f"Invalid 'aux_data_dir' in config: {aux_dir}")
        if not gene_coords or not Path(gene_coords).exists():
            raise ValueError(f"Invalid 'gene_coords_file' in config: {gene_coords}")

        # Dataset-specific params
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
            "wgs_tumor": wgs_tumor,
            "wgs_normal": wgs_normal,
            "aux_dir": Path(aux_dir).resolve(),
            "gene_coords": Path(gene_coords).resolve(),
            "mincell": mincell
        }

    def build_command(
        self, 
        dataset_cfg, 
        input_files, 
        output_dir
    ) -> List[str]:
        
        if input_files.get("_skip_execution", False):
            return []
        
        script_path = self.get_script_path("run_clonalscope_wgs.R", subfolder="clonalscope_wgs")
        ds_id = dataset_cfg.get("dataset_id")
        
        # HMM States for WGS/WES seg
        hmm_states = self.model_cfg.get("hmm_states", [0.8, 1, 1.2])
        hmm_str = ",".join(map(str, hmm_states))

        cmd = [
            "Rscript",
            str(script_path),
            str(input_files["matrix_dir"]),
            str(output_dir),
            str(input_files["meta_file"]),
            str(input_files["wgs_tumor"]),
            str(input_files["wgs_normal"]),
            str(input_files["aux_dir"]),
            str(input_files["gene_coords"]),
            str(ds_id),
            str(hmm_str),
            str(input_files["mincell"])
        ]
        return cmd