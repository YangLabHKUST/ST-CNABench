# src/models/calicost_model.py

from pathlib import Path
from typing import Dict, Any, List, Optional
import logging
import yaml
import anndata as ad
import pandas as pd
from ...base import BaseModel
from ....modules.io_utils import render_text_template, get_per_dataset_param


class CalicoSTModel(BaseModel):
    """
    Wrapper for running CalicoST.
    """

    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="CalicoST",
            exec_mode=exec_mode
        )

    def is_enabled_for_dataset(self, dataset_cfg: Dict[str, Any]) -> bool:
        if not self.model_cfg.get("enabled", True):
            return False
        ds_id = dataset_cfg.get("dataset_id")
        per_ds = self.model_cfg.get("per_dataset", {})
        if ds_id in per_ds:
            return per_ds[ds_id].get("enabled", True)
        return True

    def prepare_inputs(
        self,
        dataset_cfg: Dict[str, Any],
        output_dir: Path,
        overwrite: bool = False,
    ) -> Dict[str, Path]:
        """
        Prepare for model input
        args:
            dataset_cfg: Dataset configuration dictionary.
            output_dir: Directory to store outputs.
            overwrite: Whether to overwrite existing prepared inputs.
        """
        ds_id = dataset_cfg.get("dataset_id", "unknown_dataset")

        # Overwrite check
        final_res = output_dir / "clone3_rectangle0_w1.0"
        if not overwrite and final_res.exists():
            logging.warning(f"[{ds_id}] CalicoST output exists. Skipping.")
            return {"_skip_execution": True}

        # File PATH
        spaceranger_like_dir = Path(dataset_cfg["output"]["root"]).resolve()

        calicost_dir = Path(self.model_cfg["calicost_dir"]).resolve()
        eagle_dir = Path(self.model_cfg["eagle_dir"]).resolve()
        region_vcf = Path(self.model_cfg["region_vcf"]).resolve()
        phasing_panel = Path(self.model_cfg["phasing_panel"]).resolve()

        template_dir_cfg = self.model_cfg.get("template_dir")
        if template_dir_cfg:
            template_dir_cfg = template_dir_cfg.replace(
                "${project_root}", str(self.project_root)
            )
            template_dir = Path(template_dir_cfg).resolve()
        else:
            template_dir = Path(__file__).resolve().parent / "config"

        config_dir = output_dir / "configs"
        preproc_dir = output_dir / "preprocess"

        config_yaml = config_dir / "config.yaml"
        purity_cfg = config_dir / "configuration_purity"
        cna_cfg = config_dir / "configuration_cna"

        config_dir.mkdir(parents=True, exist_ok=True)
        preproc_dir.mkdir(parents=True, exist_ok=True)
        '''
        # Normal idx (Do not support now)
        label_file = None
        dataset_root = Path(dataset_cfg["output"]["root"]).resolve()
        for f in dataset_root.iterdir():
            if "tumor_normal" in f.name:
                label_file = f
                break
        if not label_file:
            raise FileNotFoundError(f"[{ds_id}] Annotation file not found in {dataset_root}")
        normal_spots_file = output_dir / f"{ds_id}_normal_spots.txt"
        df_labels = pd.read_csv(label_file, sep="\t")
        normal_barcodes = df_labels.loc[df_labels['tumor_normal'] == 'normal', 'Barcode']
        normal_indices = df_labels.index[df_labels['Barcode'].isin(normal_barcodes)].tolist()
        with open(normal_spots_file, "w") as f:
            for idx in normal_indices:
                f.write(f"{idx}\n")
        '''
        # Tumor purity file generate from step 2
        raw_use_tumor_purity = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="use_tumor_purity",
            default_key="use_tumor_purity",
            default_value=True,
        )
        if isinstance(raw_use_tumor_purity, str):
            use_tumor_purity = raw_use_tumor_purity.strip().lower() in {"1", "true", "yes", "y"}
        else:
            use_tumor_purity = bool(raw_use_tumor_purity)

        if use_tumor_purity:
            tumor_purity_file = output_dir / "loh_estimator_tumor_prop.tsv"
        else:
            tumor_purity_file = None

        # Dataset-specific params
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
            default_value=2
        )
        cores = self.model_cfg.get("n_threads", 20)

        context = {
            "{{ calicost_dir }}": str(calicost_dir),
            "{{ eagle_dir }}": str(eagle_dir),
            "{{ region_vcf }}": str(region_vcf),
            "{{ phasing_panel }}": str(phasing_panel),
            "{{ spaceranger_dir }}": str(spaceranger_like_dir),
            "{{ preprocess_dir }}": str(preproc_dir),
            "{{ output_dir }}": str(output_dir),
            "{{ tumor_purity_file }}": str(tumor_purity_file) if tumor_purity_file else "None",
            "{{ n_clones }}": str(n_clones),
            "{{ UMItag }}": str(umi_tag),
            "{{ cellTAG }}": str(cell_tag),
            "{{ cores }}": str(cores)
        }

        # Config template
        config_template = template_dir / "config.yaml.template"
        render_text_template(config_template, config_yaml, context)
        purity_template = template_dir / "configuration_purity.template"
        render_text_template(purity_template, purity_cfg, context)

        cna_template = template_dir / "configuration_cna.template"
        render_text_template(cna_template, cna_cfg, context)

        runtime_flags = {
            "dataset_id": ds_id,
            "model_name": self.name,
            "use_tumor_purity": use_tumor_purity,
            "n_clones": n_clones,
            "UMItag": umi_tag,
            "cellTAG": cell_tag,
            "n_threads": cores,
        }
        runtime_flags_path = config_dir / "model_runtime_flags.yaml"
        runtime_flags_path.write_text(yaml.safe_dump(runtime_flags, sort_keys=False))

        logging.info(f"[{ds_id}][CalicoST] Configs generated in {config_dir}")

        return {
            "config_yaml": config_yaml,
            "purity_config": purity_cfg,
            "cna_config": cna_cfg,
            "calicost_dir": calicost_dir,
            "preprocess_dir": preproc_dir,
            "use_tumor_purity": use_tumor_purity,
            "runtime_flags": runtime_flags_path,
        }

    def build_command(
        self,
        dataset_cfg: Dict[str, Any],
        input_files: Dict[str, Path],
        output_dir: Path,
    ) -> List[str]:
        """
        Build the command to run CalicoST.
        """
        if input_files.get("_skip_execution", False):
            return []
        # src/models/tools_scripts/calicost/run_calicost.sh
        script_path = self.get_script_path(
            script_name="run_calicost.sh",
            subfolder="calicost",
        )
        cores = self.model_cfg.get("n_threads", 20)

        cmd: List[str] = [
            "bash",
            str(script_path),
            str(input_files["calicost_dir"]),
            str(input_files["config_yaml"]),
            str(input_files["purity_config"]),
            str(input_files["cna_config"]),
            str(input_files["preprocess_dir"]),
            str(cores),
            str(input_files.get("use_tumor_purity", True)).lower(),
        ]
        return cmd
