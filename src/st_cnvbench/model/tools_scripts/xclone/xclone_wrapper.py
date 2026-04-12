from pathlib import Path
from typing import Dict, Any, List
import yaml
import logging

from ...base import BaseModel
from ....modules.io_utils import get_per_dataset_param

class XcloneModel(BaseModel):
    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="Xclone",
            exec_mode=exec_mode
        )

    def prepare_inputs(self, dataset_cfg, output_dir, overwrite=False) -> Dict[str, Any]:
        ds_id = dataset_cfg.get("dataset_id")
        dataset_root = Path(dataset_cfg["output"]["root"]).resolve()
        
        # Overwrite check
        final_res = output_dir / "xclone_output"
        if not overwrite and final_res.exists():
            logging.warning(f"[{ds_id}] Xclone output exists. Skipping.")
            return {"_skip_execution": True}
        # MTX
        matrix_dir = dataset_root / "filtered_feature_bc_matrix"
        barcodes_file = matrix_dir / "barcodes.tsv.gz"
        
        if not barcodes_file.exists():
            raise FileNotFoundError(f"[{ds_id}] Barcodes file not found in {matrix_dir}")

        # Annot
        label_file = None
        for f in dataset_root.iterdir():
            if "tumor_normal" in f.name:
                label_file = f
                break
        if not label_file:
            raise FileNotFoundError(f"[{ds_id}] Annotation file not found")

        # Spatial
        spatial_file = dataset_root / "spatial" / "tissue_positions.csv"
        if not spatial_file.exists():
             raise FileNotFoundError(f"[{ds_id}] Spatial positions file not found")

        # BAM
        bam_file = None
        for f in dataset_root.iterdir():
            if f.suffix == ".bam" and ("possorted" in f.name or "genome" in f.name):
                bam_file = f
                break
        if not bam_file:
             raise FileNotFoundError(f"[{ds_id}] BAM file not found")

        # External tools and Ref data
        deps = {
            "snp_vcf": Path(self.model_cfg.get("snp_vcf")).resolve(),
            "gene_region": Path(self.model_cfg.get("gene_region")).resolve(),
            "eagle_path": Path(self.model_cfg.get("eagle_path")).resolve(),
            "genetic_map": Path(self.model_cfg.get("genetic_map")).resolve(),
            "panel_dir": Path(self.model_cfg.get("panel_dir")).resolve()
        }
        for name, path_str in deps.items():
            if not path_str or not Path(path_str).exists():
                raise FileNotFoundError(f"[Xclone] Missing or Invalid config: {name} -> {path_str}")

        # Dataset-specific params
        umi_tag = get_per_dataset_param(
            model_cfg=self.model_cfg, dataset_cfg=dataset_cfg,
            key="UMItag", default_key="UMItag", default_value="Auto"
        )
        cell_tag = get_per_dataset_param(
            model_cfg=self.model_cfg, dataset_cfg=dataset_cfg,
            key="cellTAG", default_key="cellTAG", default_value="CB"
        )
        min_COUNT = get_per_dataset_param(
            model_cfg=self.model_cfg, dataset_cfg=dataset_cfg,
            key="minCOUNT", default_key="minCOUNT", default_value="11"
        )
        min_MAF = get_per_dataset_param(
            model_cfg=self.model_cfg, dataset_cfg=dataset_cfg,
            key="minMAF", default_key="minMAF", default_value="0.1"
        )
        n_clusters = get_per_dataset_param(
            model_cfg=self.model_cfg, dataset_cfg=dataset_cfg,
            key="n_clusters", default_key="n_clusters", default_value="3"
        )
        return {
            "bam_file": bam_file,
            "barcodes_file": barcodes_file,
            "matrix_dir": matrix_dir,
            "meta_file": label_file,
            "spatial_file": spatial_file,
            "deps": deps,
            "umi_tag": umi_tag,
            "cell_tag": cell_tag,
            "min_COUNT": min_COUNT,
            "min_MAF": min_MAF,
            "n_clusters": n_clusters
        }

    def build_command(self, dataset_cfg, input_files, output_dir) -> List[str]:
        if input_files.get("_skip_execution", False):
            return []
        
        # src/model/tools_scripts/xclone/run_xclone.sh
        script_path = self.get_script_path("run_xclone.sh", subfolder="xclone")
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
            str(input_files["spatial_file"]), 
            str(deps["snp_vcf"]),      
            str(deps["gene_region"]),  
            str(deps["genetic_map"]),  
            str(deps["eagle_path"]),   
            str(deps["panel_dir"]),    
            str(input_files["umi_tag"]),
            str(input_files["cell_tag"]),
            str(input_files["min_COUNT"]),
            str(input_files["min_MAF"]),
            str(input_files["n_clusters"])
        ]
        return cmd