from pathlib import Path
from typing import Dict, Any, List, Optional
import logging
import pandas as pd
import scanpy as sc
import os

from ...base import BaseModel
from ....modules.io_utils import get_per_dataset_param


class STARCHModel(BaseModel):
    """
    Wrapper for running STARCH on a given dataset.
    """

    def __init__(self, project_root: str, model_cfg: str, result_dir: str = None, exec_mode: str = "docker"):
        super().__init__(
            project_root=project_root,
            model_cfg=model_cfg,
            result_dir=result_dir,
            model_name="STARCH",
            exec_mode=exec_mode
        )

    def prepare_inputs(
        self,
        dataset_cfg: Dict[str, Any],
        output_dir: Path,
        overwrite: bool = False
    ) -> Dict[str, Path]:
        """
        Convert standardized dataset to STARCH-compatible format.
        
        Returns:
            Dictionary containing paths to prepared input files
        """
        ds_id = dataset_cfg.get("dataset_id", "unknown_dataset")
        dataset_root = Path(dataset_cfg["output"]["root"]).resolve()
        ref_norm = dataset_cfg.get("ref_norm")
        # MTX and Spatial
        matrix_dir = dataset_root / "filtered_feature_bc_matrix"
        coords_file = dataset_root / "spatial" / "tissue_positions.csv"
        label_file = None
        if ref_norm is True:
            for f in dataset_root.iterdir():
                if "tumor_normal" in f.name:
                    label_file = f
                    break
        if ref_norm is True and not label_file:
            raise FileNotFoundError(f"[{ds_id}] Annotation file not found in {dataset_root}")

        # Gene Annot
        gene_mapping = self.model_cfg.get("gene_mapping_file")
        if not gene_mapping:
            raise ValueError("[STARCH] Missing required config: 'gene_mapping_file' in models.yaml")
        gene_mapping_file = Path(gene_mapping).resolve()
        if not gene_mapping_file.exists():
            raise FileNotFoundError(f"[STARCH] Gene mapping file not found: {gene_mapping_file}")
        
        # Processed dir
        processed_dir = output_dir / "processed_input"
        processed_dir.mkdir(parents=True, exist_ok=True)
        processed_count_file = processed_dir / f"{ds_id}_counts_coords.tsv"
        processed_normal_file = processed_dir / f"{ds_id}_normal_spots.txt"
        processed_gene_mapping = processed_dir / "gene_mapping_processed.txt"
        
        # Overwrite check
        result_file = output_dir / f"states_{ds_id}.csv"
        if not overwrite and processed_count_file.exists() and processed_gene_mapping.exists():
            logging.info(f"[{ds_id}][STARCH] Found existing processed inputs. Skipping conversion.")
            normal_spots_file = processed_normal_file if ref_norm is True else None
            return {
                "count_file": processed_count_file,
                "normal_spots_file": normal_spots_file,
                "gene_mapping_file": processed_gene_mapping,
                "_skip_execution": result_file.exists()
            }
        
        # Process
        logging.info(f"[{ds_id}][STARCH] Converting count matrix to STARCH format...")
        self._process_STARCH_input(matrix_dir, coords_file, label_file, processed_count_file, processed_normal_file)
        
        # Process gene mapping file
        logging.info(f"[{ds_id}][STARCH] Processing gene mapping file...")
        self._process_gene_mapping(gene_mapping_file, processed_gene_mapping)
        
        # With/Without normal spots
        if ref_norm is False:
            logging.info(f"[{ds_id}] Model run in no reference normal mode.")
            return {
                "count_file": processed_count_file,
                "normal_spots_file": None,
                "gene_mapping_file": processed_gene_mapping,
                "preprocess_dir": processed_dir,
                "_skip_execution": False
            }
        else:
            return {
                "count_file": processed_count_file,
                "normal_spots_file": processed_normal_file,
                "gene_mapping_file": processed_gene_mapping,
                "preprocess_dir": processed_dir,
                "_skip_execution": False
            }

    def build_command(
        self,
        dataset_cfg: Dict[str, Any],
        input_files: Dict[str, Path],
        output_dir: Path
    ) -> List[str]:
        """
        Build the command to run STARCH.
        
        For STARCH, we call a Python script instead of a shell script.
        """
        # src/model/tools_scripts/starch/run_starch.py
        script_path = self.get_script_path(
            script_name="run_starch.py",
            subfolder="starch"
        )
        
        ds_id = dataset_cfg.get("dataset_id", "unknown_dataset")
        ref_norm = dataset_cfg.get("ref_norm")
        # Dataset-specific parameters
        n_clusters = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="n_clusters",
            default_key="n_clusters",
            default_value=2
        )
        
        beta_spot = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="beta_spot",
            default_key="beta_spot",
            default_value=2
        )
        
        platform = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="platform",
            default_key="platform",
            default_value="Visium"
        )

        returnnormal = get_per_dataset_param(
            model_cfg=self.model_cfg,
            dataset_cfg=dataset_cfg,
            key="returnnormal",
            default_key="returnnormal",
            default_value=True
        )

        # Build command
        cmd = [
            "python",
            str(script_path),
            "--input", str(input_files["count_file"]),
            "--gene_mapping_file_name", str(input_files["gene_mapping_file"]),
            "--outdir", str(output_dir),
            "--n_clusters", str(n_clusters),
            "--beta_spot", str(beta_spot),
            "--platform", str(platform),
            "--output", str(ds_id)
        ]
        
        if input_files["normal_spots_file"] is not None:
            cmd.extend(["--normal_spots", str(input_files["normal_spots_file"])])
        
        cmd.extend(["--returnnormal", str(int(bool(returnnormal)))])
        
        return cmd
    # Helper methods
    def _process_STARCH_input(self, count_file_path: Path, coord_file_path: Path, label_file: Optional[Path], output_file_path: Path, normal_spots_file: Path):
        """
        Convert spatial transcriptomics count matrix and coordinate file to STARCH-compatible TSV format.
        """
        adata = sc.read_10x_mtx(
            count_file_path,
            var_names="gene_ids",
            make_unique=True
        )
        original_barcodes = adata.obs_names.to_list()
        barcode_to_index = {bc: i for i, bc in enumerate(original_barcodes)}
        if label_file is not None:
            df_labels = pd.read_csv(label_file, sep="\t")
            normal_barcodes = df_labels.loc[
                (df_labels['tumor_normal'] == 'normal') | (df_labels['tumor_normal'] == 'Normal'),
                'Barcode'
            ]

            normal_indices = [
                barcode_to_index[bc]
                for bc in normal_barcodes
                if bc in barcode_to_index
            ]

            with open(normal_spots_file, "w") as f:
                for idx in normal_indices:
                    f.write(f"{idx}\n")

        coords = pd.read_csv(coord_file_path)
        coords = coords.set_index("barcode").loc[original_barcodes]
        adata.obsm["spatial"] = coords[["array_row", "array_col"]].values
        new_obs_names = [
            f"{int(r)}x{int(c)}"
            for r, c in coords[["array_row", "array_col"]].values
        ]
        adata.obs_names = new_obs_names
        mapping_df = pd.DataFrame({
            "barcode": original_barcodes,
            "spot_name": new_obs_names,
            "spot_index": range(len(original_barcodes))
        })

        mapping_df.to_csv(
            output_file_path.parent / f"{output_file_path.stem}_spot_mapping.tsv",
            sep="\t",
            index=False
        )
        x_data = adata.X
        from scipy.sparse import issparse
        if issparse(x_data):
            x_data = x_data.toarray()

        # Build a genes-by-spots table for the STARCH text input.
        expr_df = pd.DataFrame(
            x_data.T,
            index=adata.var_names,
            columns=adata.obs_names
        )

        # Write the matrix as a tab-delimited file.
        expr_df.to_csv(output_file_path, sep="\t")
        logging.info(f"[STARCH] Converted count matrix: {expr_df.shape[0]} genes, {expr_df.shape[1]} spots")

    def _process_gene_mapping(self, file_path: Path, output_file: Path):
        """
        Convert gene data file to DataFrame with index and column names.
        """
        df = pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1, 2, 3])
        
        df.columns = ['name2', 'chrom', 'cdsStart', 'cdsEnd']
        df.reset_index(drop=True, inplace=True)
        # add "chr"
        df['chrom'] = df['chrom'].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else str(x))
        df.to_csv(output_file, index=True, sep='\t')
        
        logging.info(f"[STARCH] Processed gene mapping: {len(df)} genes")
