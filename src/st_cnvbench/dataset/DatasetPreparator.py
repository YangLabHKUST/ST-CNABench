import os
import shutil
import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Any, List
import scanpy as sc
import time

from ..config import load_data_config
from ..modules import io_utils
from ..modules.utils import create_hard_link
class DatasetPreparator:
    """
    Dataset Preparator & Checker.
    Make various raw data formats into a standardized structure (Spaceranger output like).
    """

    def __init__(self, config_path: str):
        self.config_path = Path(config_path)
        self.full_config = load_data_config(str(self.config_path))
        settings = self.full_config.get("project_settings", {})
        self.dataset_root = Path(settings.get("dataset_root"))
        self.output_root = Path(settings.get("output_root"))
        self.dataset_cfgs = {
            k: v for k, v in self.full_config.items()
            if k != "project_settings" and isinstance(v, dict)
        }
    def list_datasets(self) -> List[str]:
        return list(self.dataset_cfgs.keys())

    def prepare_dataset(self, dataset_key: str, overwrite: bool = False):
        # load config
        cfg = self.dataset_cfgs.get(dataset_key)
        if not cfg: raise KeyError(f"Dataset '{dataset_key}' not found.")
        ds_id = cfg.get("dataset_id", dataset_key)
        out_root = Path(cfg["output"]["root"])
        if out_root.exists() and not overwrite:
            logging.info(f"[{ds_id}] Output exists, skipping.")
            return

        # Prepare dirs
        (out_root / "filtered_feature_bc_matrix").mkdir(parents=True, exist_ok=True)
        (out_root / "spatial").mkdir(parents=True, exist_ok=True)
        (out_root / "metadata").mkdir(parents=True, exist_ok=True)

        logging.info(f"[{ds_id}] Checking & Preparing...")

        # data type
        fmt = cfg.get("format", "").lower()
        if fmt == "stpipeline":
            data_components = self._load_stpipeline(cfg)
        elif fmt == "spaceranger":
            data_components = self._load_spaceranger(cfg)
        else:
            raise ValueError(f"[{ds_id}] Unsupported format: {fmt}")

        # Valid step
        validated_data = self._validate_and_align(ds_id, *data_components)

        # Save 
        self._save_standardized(cfg, out_root, *validated_data)
        
        logging.info(f"[{ds_id}] Success. Data validated and saved.")
    
    # stpileline loader
    def _load_stpipeline(self, cfg: Dict[str, Any]):
        """
        Load STpipeline raw data files from user custom path
        """
        ds_id = cfg.get("dataset_id")
        raw = cfg["raw"]
        logging.info(f"[{ds_id}] Loading raw files (STpipeline)...")
        
        mtx = io_utils.read_counts_matrix(raw["counts"])
        barcodes = io_utils.read_barcodes(raw["barcodes"])
        features = io_utils.read_features(raw["features"])
        coords = io_utils.read_coords(raw["coords"])
        labels = io_utils.read_labels(raw["tumor_normal"]) if raw.get("tumor_normal") else None
        return mtx, barcodes, features, coords, labels

    # SpaceRanger loader
    def _load_spaceranger(self, cfg: Dict[str, Any]):
        """
        Load SpaceRanger raw data files from dataset root
        """
        ds_id = cfg.get("dataset_id")
        raw_root = Path(cfg["raw"]["root"])
        logging.info(f"[{ds_id}] Loading raw files (SpaceRanger)...")

        # Matrix + Barcodes + Features
        src_mtx_dir = raw_root / "filtered_feature_bc_matrix"
        if not src_mtx_dir.exists():
             raise FileNotFoundError(
                 f"[{ds_id}] Expected filtered_feature_bc_matrix under {raw_root}"
             )

        mtx = io_utils.read_counts_matrix(src_mtx_dir / "matrix.mtx.gz")
        barcodes = pd.read_csv(src_mtx_dir / "barcodes.tsv.gz", header=None, compression='gzip', sep="\t").iloc[:, 0].astype(str).values.tolist()
        features = pd.read_csv(src_mtx_dir / "features.tsv.gz", header=None, compression='gzip', sep="\t")

        src_spatial = raw_root / "spatial"
        pos_file = src_spatial / "tissue_positions_list.csv"
        if not pos_file.exists():
            pos_file = src_spatial / "tissue_positions.csv"
        if not pos_file.exists():
             raise FileNotFoundError(f"[{ds_id}] No Spatial Coordinates found in {src_spatial}")        
        coords = io_utils.read_coords(pos_file)
        labels = io_utils.read_labels(cfg["raw"]["tumor_normal"]) if cfg["raw"].get("tumor_normal") else None
        return mtx, barcodes, features, coords, labels

    def _validate_and_align(self, ds_id, mtx, barcodes, features, coords_df, label_df):
        """
        Ensures strict alignment between Matrix Columns and Spatial Rows.
        """
        # Alignment Matrix <-> barcodes
        if mtx.shape[1] == len(barcodes):
            pass 
        elif mtx.shape[0] == len(barcodes):
             logging.info(f"[{ds_id}] Notice: Matrix seems transposed (Spots x Genes). Transposing to standard.")
             mtx = mtx.T
        else:
             raise ValueError(f"[{ds_id}] Matrix shape {mtx.shape} does not match barcode list length {len(barcodes)}!")
        # Intersection
        bc_sets = [set(barcodes), set(coords_df.index)]
        if label_df is not None:
            if "Barcode" not in label_df.columns:
                raise ValueError(f"[{ds_id}] tumor_normal file must contain a 'Barcode' column.")
            bc_sets.append(set(label_df["Barcode"]))
        common_barcodes = [bc for bc in barcodes if bc in set.intersection(*bc_sets)]
        if len(common_barcodes) == 0:
            raise ValueError(f"[{ds_id}] Zero overlapping barcodes among configured inputs!")
        # Subset
        bc_to_idx = {bc: i for i, bc in enumerate(barcodes)}
        keep_indices = [bc_to_idx[bc] for bc in common_barcodes]
        mtx_subset = mtx[:, keep_indices]
        coords_subset = coords_df.loc[common_barcodes]
        label_df_subset = (
            label_df.set_index("Barcode").loc[common_barcodes].reset_index()
            if label_df is not None
            else None
        )
        valid_barcodes = list(common_barcodes)
        n_valid = len(valid_barcodes)
        logging.info(f"[{ds_id}] --- Validation Report ---")
        logging.info(f"[{ds_id}] Alignment Valid spot      : {n_valid}")
        

        return mtx_subset, valid_barcodes, features, coords_subset, label_df_subset

    def _save_standardized(self, cfg, out_root, mtx, barcodes, features, coords, labels):
        raw = cfg["raw"]
        ds_id = cfg.get("dataset_id", "Unknown_ID")
        fmt = cfg.get("format", "").lower()

        ffm_dir = out_root / "filtered_feature_bc_matrix"
        spatial_dir = out_root / "spatial"
        meta_dir = out_root  
        raw_root = Path(raw.get("root", "")) if "root" in raw else Path(cfg["raw"]["root"])

        # Write matrix + barcodes + features
        io_utils.write_matrix_bundle(ffm_dir, mtx, barcodes, features)

        # Write label file
        if labels is not None:
            label_out_path = out_root / f"metadata_{ds_id}_tumor_normal.tsv"
            labels.to_csv(label_out_path, sep="\t", index=False)
        # Save spatial info
        sf_path = raw.get("scalefactors")

        if not sf_path and fmt == "spaceranger":
            potential_sf = raw_root / "spatial" / "scalefactors_json.json"
            if potential_sf.exists():
                sf_path = str(potential_sf)

        if not sf_path or not Path(sf_path).exists():
            raise FileNotFoundError(
                f"[{ds_id}] Missing scalefactors JSON. Set raw.scalefactors or provide "
                f"{raw_root / 'spatial' / 'scalefactors_json.json'}."
            )

        with open(sf_path, "r") as f:
            sf_data = json.load(f)
        logging.info(f"[{ds_id}] Loaded scalefactors from {sf_path}")

        img_path = Path(raw.get("tissue_image")) if raw.get("tissue_image") else None
        if not img_path and fmt == "spaceranger":
            potential_img = raw_root / "spatial" / "tissue_hires_image.png"
            if not potential_img.exists():
                potential_img = raw_root / "spatial" / "tissue_lowres_image.png"
            if potential_img.exists():
                img_path = potential_img

        io_utils.write_spatial_bundle(spatial_dir, coords, sf_data, img_path)

        meta_dir.mkdir(parents=True, exist_ok=True)

        # copy annot and gt
        for key in ["cnv_gt"]:
            path_str = raw.get(key)
            if path_str is None:
                continue

            src = Path(path_str)
            if not src.exists():
                raise FileNotFoundError(f"[{ds_id}] Configured {key} not found: {src}")

            dst = meta_dir / src.name
            shutil.copy2(src, dst)
            logging.info(f"[{ds_id}] Copied metadata: {key} -> {dst.name}")
                
        # BAM and BAI files (hard link)
        bam_path_str = raw.get("bam")
        bai_path_str = raw.get("bai")

        if not bam_path_str and fmt == "spaceranger":
            potential_bam = raw_root / "possorted_genome_bam.bam"
            if potential_bam.exists():
                bam_path_str = str(potential_bam)
                potential_bai = raw_root / "possorted_genome_bam.bam.bai"
                if potential_bai.exists():
                    bai_path_str = str(potential_bai)

        if bam_path_str:
            bam_src = Path(bam_path_str)
            if not bam_src.exists():
                raise FileNotFoundError(f"[{ds_id}] BAM file not found: {bam_src}")

            if bai_path_str:
                bai_src = Path(bai_path_str)
                if not bai_src.exists():
                    raise FileNotFoundError(
                        f"[{ds_id}] BAI file specified but not found: {bai_src}"
                    )

                # Check BAI timestamp
                if bai_src.stat().st_mtime < bam_src.stat().st_mtime:
                    raise ValueError(
                        f"[{ds_id}] CRITICAL: BAI file is older than BAM file!\n"
                        f"  BAM: {bam_src} ({time.ctime(bam_src.stat().st_mtime)})\n"
                        f"  BAI: {bai_src} ({time.ctime(bai_src.stat().st_mtime)})\n"
                        f"  Please re-index using 'samtools index'."
                    )

                target_bai_name = "possorted_genome_bam.bam.bai"
                create_hard_link(bai_src, out_root / target_bai_name)
                logging.info(f"[{ds_id}] Linked BAI: {target_bai_name}")

            target_bam_name = "possorted_genome_bam.bam"
            create_hard_link(bam_src, out_root / target_bam_name)
            logging.info(f"[{ds_id}] Linked BAM: {target_bam_name}")

        # WGS/WES Counts file (for Clonalscope_WGS)
        file_mapping = {
            "wgs_wes_tumor_bedg":  "wgs_wes_tumor_data.bedg",
            "wgs_wes_normal_bedg": "wgs_wes_normal_data.bedg"
        }
        
        for cfg_key, target_name in file_mapping.items():
            path_str = raw.get(cfg_key)
            if path_str is None:
                continue

            src = Path(path_str)
            if not src.exists():
                raise FileNotFoundError(f"[{ds_id}] Configured {cfg_key} not found: {src}")

            dst = meta_dir / target_name
            shutil.copy2(src, dst)
            logging.info(f"[{ds_id}] Copied WGS/WES data: {cfg_key} -> {dst.name}")
        if ffm_dir.exists():
            logging.info(
                f"[{ds_id}] Creating filtered_feature_bc_matrix.h5ad from {ffm_dir}"
            )
            adata = sc.read_10x_mtx(
                ffm_dir,
                var_names="gene_symbols",
                make_unique=True,
            )
            adata.layers["count"] = adata.X.copy()

            out_h5ad = out_root / "filtered_feature_bc_matrix.h5ad"
            adata.write_h5ad(out_h5ad)
            logging.info(f"[{ds_id}] Wrote {out_h5ad}")
        else:
            logging.warning(
                f"[{ds_id}] filtered_feature_bc_matrix directory not found under {out_root}; "
                "skip h5ad creation."
            )
