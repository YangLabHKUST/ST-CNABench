from pathlib import Path
from typing import Any, Dict, List
import logging

from ...base import BaseModel
from ....modules.io_utils import get_per_dataset_param


class SlideCNAModel(BaseModel):
    """Wrapper for running SlideCNA on standardized spatial transcriptomics data."""

    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "conda"):
        if exec_mode != "conda":
            logging.warning("[SlideCNA] Forcing conda execution because SlideCNA is installed in a conda env.")
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="SlideCNA",
            exec_mode="conda",
        )

    def is_enabled_for_dataset(self, dataset_cfg: Dict[str, Any]) -> bool:
        """Limit SlideCNA execution to configured spatial platforms."""
        if not super().is_enabled_for_dataset(dataset_cfg):
            return False
        allowed_platforms = self.model_cfg.get("allowed_platforms", ["SlideRNAseq_v2"])
        platform = dataset_cfg.get("platform")
        if allowed_platforms and platform not in allowed_platforms:
            logging.warning(f"[SlideCNA] Skipped: platform {platform} is not in {allowed_platforms}")
            return False
        return True

    def prepare_inputs(
        self,
        dataset_cfg: Dict[str, Any],
        output_dir: Path,
        overwrite: bool = False,
    ) -> Dict[str, Any]:
        """Validate standardized inputs and collect SlideCNA parameters."""
        ds_id = dataset_cfg.get("dataset_id", "unknown_dataset")
        dataset_root = Path(dataset_cfg["output"]["root"]).resolve()

        final_outputs = [
            output_dir / "md_bin.txt",
            output_dir / "cluster_labels_malig.txt",
            output_dir / "slidecna_malig_barcode_labels.tsv.gz",
            output_dir / "slidecna_malig_clone_profiles.tsv.gz",
            output_dir / "SlideCNA.log",
        ]
        if not overwrite and all(path.exists() for path in final_outputs):
            logging.warning(f"[{ds_id}] SlideCNA output exists. Skipping.")
            return {"_skip_execution": True}

        if dataset_cfg.get("ref_norm") is not True:
            raise ValueError(f"[{ds_id}] SlideCNA requires ref_norm=True with Normal reference spots.")

        matrix_dir = dataset_root / "filtered_feature_bc_matrix"
        positions_file = dataset_root / "spatial" / "tissue_positions.csv"
        metadata_file = dataset_root / f"metadata_{ds_id}_tumor_normal.tsv"
        if not metadata_file.exists():
            matches = sorted(dataset_root.glob("*tumor_normal*.tsv"))
            if len(matches) != 1:
                raise FileNotFoundError(f"[{ds_id}] Expected one tumor_normal metadata file in {dataset_root}")
            metadata_file = matches[0]

        required_files = [
            matrix_dir / "matrix.mtx.gz",
            matrix_dir / "features.tsv.gz",
            matrix_dir / "barcodes.tsv.gz",
            positions_file,
            metadata_file,
        ]
        for path in required_files:
            if not path.exists():
                raise FileNotFoundError(f"[{ds_id}] Required SlideCNA input not found: {path}")

        gene_annot_file = Path(
            get_per_dataset_param(
                model_cfg=self.model_cfg,
                dataset_cfg=dataset_cfg,
                key="gene_annot_file",
                default_key="gene_annot_file",
            )
        ).resolve()
        if not gene_annot_file.exists():
            raise FileNotFoundError(f"[{ds_id}] SlideCNA gene annotation not found: {gene_annot_file}")

        chrom_ord = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="chrom_ord",
            default_key="chrom_ord",
        )
        if not chrom_ord:
            chrom_ord = self._default_chrom_order(dataset_cfg)

        plot_dir = output_dir / "plots"
        plot_dir.mkdir(parents=True, exist_ok=True)

        return {
            "dataset_root": dataset_root,
            "matrix_dir": matrix_dir,
            "metadata_file": metadata_file,
            "positions_file": positions_file,
            "gene_annot_file": gene_annot_file,
            "plot_dir": plot_dir,
            "chrom_ord": self._format_chrom_order(chrom_ord),
            "roll_mean_window": get_per_dataset_param(self.model_cfg, dataset_cfg, "roll_mean_window", "roll_mean_window", 101),
            "avg_bead_per_bin": get_per_dataset_param(self.model_cfg, dataset_cfg, "avg_bead_per_bin", "avg_bead_per_bin", 12),
            "spatial": get_per_dataset_param(self.model_cfg, dataset_cfg, "spatial", "spatial", True),
            "pos": get_per_dataset_param(self.model_cfg, dataset_cfg, "pos", "pos", True),
            "pos_k": get_per_dataset_param(self.model_cfg, dataset_cfg, "pos_k", "pos_k", 55),
            "ex_k": get_per_dataset_param(self.model_cfg, dataset_cfg, "ex_k", "ex_k", 1),
            "use_GO_terms": get_per_dataset_param(self.model_cfg, dataset_cfg, "use_GO_terms", "use_GO_terms", False),
            "max_k_silhouette": get_per_dataset_param(self.model_cfg, dataset_cfg, "max_k_silhouette", "max_k_silhouette", 10),
        }

    def build_command(
        self,
        dataset_cfg: Dict[str, Any],
        input_files: Dict[str, Any],
        output_dir: Path,
    ) -> List[str]:
        """Build the Rscript command used to run SlideCNA."""
        if input_files.get("_skip_execution", False):
            return []

        script_path = self.get_script_path("run_slidecna.R", subfolder="slidecna")
        ds_id = dataset_cfg.get("dataset_id", "unknown_dataset")

        return [
            "Rscript",
            str(script_path),
            str(input_files["matrix_dir"]),
            str(input_files["metadata_file"]),
            str(input_files["positions_file"]),
            str(input_files["gene_annot_file"]),
            str(output_dir),
            str(input_files["plot_dir"]),
            str(ds_id),
            str(input_files["chrom_ord"]),
            str(input_files["roll_mean_window"]),
            str(input_files["avg_bead_per_bin"]),
            str(input_files["spatial"]),
            str(input_files["pos"]),
            str(input_files["pos_k"]),
            str(input_files["ex_k"]),
            str(input_files["use_GO_terms"]),
            str(input_files["max_k_silhouette"]),
        ]

    def _default_chrom_order(self, dataset_cfg: Dict[str, Any]) -> List[str]:
        """Return the standard chromosome order for supported genomes."""
        genome = str(dataset_cfg.get("genome", "")).lower()
        species = str(dataset_cfg.get("species", "")).lower()
        if genome in {"hg38", "hg19"} or species == "human":
            return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        if genome in {"mm10", "mm39"} or species == "mouse":
            return [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]
        raise ValueError(f"[SlideCNA] Cannot infer chrom_ord for genome={genome}, species={species}")

    def _format_chrom_order(self, chrom_ord: Any) -> str:
        """Convert YAML list or comma-separated chromosome order to a CLI string."""
        if isinstance(chrom_ord, str):
            return chrom_ord
        if isinstance(chrom_ord, list):
            return ",".join(str(chrom) for chrom in chrom_ord)
        raise ValueError(f"[SlideCNA] chrom_ord must be a list or comma-separated string, got {type(chrom_ord)}")
