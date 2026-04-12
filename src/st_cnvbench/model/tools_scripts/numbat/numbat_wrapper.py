from pathlib import Path
from typing import Dict, Any, List
import yaml
import logging
from ...base import BaseModel
from ....modules.io_utils import get_per_dataset_param

class NumbatModel(BaseModel):
    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="Numbat",
            exec_mode=exec_mode
        )

    def prepare_inputs(self, dataset_cfg, output_dir, overwrite=False) -> Dict[str, Path]:
        ds_id = dataset_cfg.get("dataset_id")
        dataset_root = Path(dataset_cfg["output"]["root"]).resolve()
        ref_norm = dataset_cfg.get("ref_norm")
        # Overwrite check
        final_res = output_dir / "bulk_clones_final.tsv.gz"
        if not overwrite and final_res.exists():
            logging.warning(f"[{ds_id}] Numbat output exists. Skipping.")
            return {"_skip_execution": True}
        # MTX
        matrix_dir = dataset_root / "filtered_feature_bc_matrix"
        if not dataset_root.exists():
             raise FileNotFoundError(f"[{ds_id}] Dataset root not found: {dataset_root}")
        barcodes_file = matrix_dir / "barcodes.tsv.gz"
        if not barcodes_file.exists():
            raise FileNotFoundError(f"[{ds_id}] Barcodes file not found in {matrix_dir}")

        # Annot
        if ref_norm is False:
            # Without ref normal model
            logging.info(f"[{ds_id}] Model run in no reference normal mode.")
        else:
            label_file = None
            for f in dataset_root.iterdir():
                if "tumor_normal" in f.name:
                    label_file = f
                    break
            if not label_file:
                raise FileNotFoundError(f"[{ds_id}] Annotation file not found in {dataset_root}")

        # BAM File
        bam_file = None
        for f in dataset_root.iterdir():
            if f.suffix == ".bam" and ("possorted" in f.name or "genome" in f.name):
                bam_file = f
                break
        if not bam_file:
             raise FileNotFoundError(f"[{ds_id}] BAM file not found in {dataset_root}")

        # Ref data and external tools
        deps = {
            "pileup_script": Path(self.model_cfg.get("pileup_script")).resolve(),
            "eagle_path": Path(self.model_cfg.get("eagle_path")).resolve(),
            "genetic_map": Path(self.model_cfg.get("genetic_map")).resolve(),
            "snp_vcf": Path(self.model_cfg.get("snp_vcf_path")).resolve(),
            "panel_dir": Path(self.model_cfg.get("panel_dir")).resolve()
        }
        for name, path_str in deps.items():
            if not path_str:
                raise ValueError(f"[Numbat] Missing required config: '{name}' in models.yaml")
            if not Path(path_str).exists():
                raise FileNotFoundError(f"[Numbat] Configured '{name}' path does not exist: {path_str}")

        # Dataset-specific parameters
        umi_tag = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="UMItag",         
            default_key="UMItag", 
            default_value="Auto" 
        )

        cell_tag = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="cellTAG",
            default_key="cellTAG",
            default_value="CB"
        )
        n_clones = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="n_clones",
            default_key="n_clones",
            default_value=3
        )

        if ref_norm is False:
            return {
                "matrix_dir": matrix_dir,
                "meta_file": None,
                "bam_file": bam_file,
                "barcodes_file": barcodes_file,
                "deps": deps,
                "umi_tag": umi_tag,
                "cell_tag": cell_tag,
                "n_clones": n_clones
            }
        else:
            return {
                "matrix_dir": matrix_dir,
                "meta_file": label_file,
                "bam_file": bam_file,
                "barcodes_file": barcodes_file,
                "deps": deps,
                "umi_tag": umi_tag,
                "cell_tag": cell_tag,
                "n_clones": n_clones
            }

    def build_command(self, dataset_cfg, input_files, output_dir) -> List[str]:

        if input_files.get("_skip_execution", False):
            return []
        
        # src/model/tools_scripts/numbat/run_numbat.sh
        script_path = self.get_script_path(
            script_name="run_numbat.sh", 
            subfolder="numbat"
        )
        ds_id = dataset_cfg.get("dataset_id")
        deps = input_files["deps"]

        cmd = [
            "bash",
            str(script_path),
            str(ds_id),                  
            str(input_files["bam_file"]),
            str(input_files["barcodes_file"]), 
            str(output_dir),            
            str(self.model_cfg.get("n_threads", 20)), 
            str(input_files["matrix_dir"]), 
            str(input_files["meta_file"]),  
            str(self.model_cfg.get("genome_version", "hg38")), 
            str(deps["pileup_script"]),  
            str(deps["eagle_path"]),     
            str(deps["genetic_map"]),    
            str(deps["snp_vcf"]),       
            str(deps["panel_dir"]),      
            str(input_files["umi_tag"]),      
            str(input_files["cell_tag"]),
            str(input_files["n_clones"])
        ]
        return cmd