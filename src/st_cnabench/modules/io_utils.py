import os
import gzip
import shutil
import json
import logging
import pandas as pd
from scipy import io, sparse
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Any

# -------------------------------------------------------------------------
# Data preparation utils
# -------------------------------------------------------------------------

# Readers
def read_counts_matrix(path: Path) -> sparse.csr_matrix:
    mtx = io.mmread(str(path))
    if sparse.issparse(mtx):
        return mtx.tocsr()
    return sparse.csr_matrix(mtx)

def read_barcodes(path: Path) -> List[str]:
    return pd.read_csv(path, header=None, sep="\t").iloc[:, 0].astype(str).values.tolist()

def read_features(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, header=None, sep="\t")

def read_labels(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=0, index_col=None)
    return df
def read_coords(path: Path) -> pd.DataFrame:
    """
    Read spatial coordinates file (csv/tsv).
    Auto-detects separator and header presence.
    Standardizes column names and returns DataFrame indexed by 'barcode'.
    """
    path = Path(path)
    # define sep
    sep = "," if path.suffix.lower() == ".csv" else "\t"

    # deal with header
    first_row = pd.read_csv(path, sep=sep, nrows=1, header=None).iloc[0]
    first_row_vals = [str(x).lower().strip() for x in first_row.values]

    # barcode header
    possible_bc_headers = {"barcode", "barcodes", "spot", "spot_id", "cell", "cell_id"}
    has_header = bool(set(first_row_vals) & possible_bc_headers)

    # expected col
    expected_cols = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",  
        "pxl_col_in_fullres", 
    ]

    if has_header:
        df = pd.read_csv(path, sep=sep, header=0)
        df.columns = [c.lower().strip() for c in df.columns]

        found_bc = next((c for c in df.columns if c in possible_bc_headers), None)
        if not found_bc:
            raise ValueError(
                f"Header detected but no barcode column found "
                f"(expected one of {possible_bc_headers})."
            )

        if found_bc != "barcode":
            df.rename(columns={found_bc: "barcode"}, inplace=True)

        if df.shape[1] != 6:
            logging.warning(
                f"Coordinates file has {df.shape[1]} columns (expected 6). "
                f"Proceeding based on header names."
            )

    else:
        logging.info("No header detected in coordinates file. Assigning standard column names.")
        df = pd.read_csv(path, sep=sep, header=None)

        if df.shape[1] != 6:
            raise ValueError(
                f"Coordinates file without header has {df.shape[1]} columns "
                f"(expected 6)."
            )

        df.columns = expected_cols

    # Final cleanup
    df["barcode"] = df["barcode"].astype(str)

    if df["barcode"].duplicated().any():
        logging.warning("Duplicate barcodes found in coordinates. Keeping first occurrence.")
        df = df.drop_duplicates(subset="barcode", keep="first")

    return df.set_index("barcode")

# Writers
def write_matrix_bundle(out_dir: Path, mtx: sparse.csr_matrix, barcodes: List[str], features_df: pd.DataFrame):
    """Write MTX Barcodes and Features files"""
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # MTX
    mtx_path = out_dir / "matrix.mtx"
    io.mmwrite(str(mtx_path), mtx)
    gzip_file(mtx_path)
    
    # Barcodes
    bc_path = out_dir / "barcodes.tsv.gz"
    with gzip.open(bc_path, "wt") as f:
        f.write("\n".join(barcodes) + "\n")
        
    # Features 
    clean_feats = pd.DataFrame()
    # if only 1 column, use it for id and name
    if features_df.shape[1] == 1:
        clean_feats["id"] = features_df.iloc[:, 0]
        clean_feats["name"] = features_df.iloc[:, 0]
        clean_feats["type"] = "Gene Expression"
    elif features_df.shape[1] >= 2:
        clean_feats["id"] = features_df.iloc[:, 0]
        clean_feats["name"] = features_df.iloc[:, 1]
        clean_feats["type"] = features_df.iloc[:, 2] if features_df.shape[1] >= 3 else "Gene Expression"
            
    ft_path = out_dir / "features.tsv.gz"
    clean_feats.to_csv(ft_path, sep="\t", header=False, index=False, compression="gzip")

def write_spatial_bundle(out_dir: Path, coords_df: pd.DataFrame, scalefactors: dict, image_path: Optional[Path]):
    """Write Spatial folder content (spaceranger-like)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # tissue_positions
    out_cols = ["in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    for c in out_cols:
        if c not in coords_df.columns:
            coords_df[c] = 0 
            
    coords_df[out_cols].to_csv(out_dir / "tissue_positions.csv", index=True, index_label="barcode")
    
    # scalefactors
    with open(out_dir / "scalefactors_json.json", "w") as f:
        json.dump(scalefactors, f)
        
    # Image
    if image_path and image_path.exists():
        # hires or lowres
        dst_name = "tissue_lowres_image.png" if "lowres" in image_path.name.lower() else "tissue_hires_image.png"
        shutil.copy2(image_path, out_dir / dst_name)

# -------------------------------------------------------------------------
# Utils
# -------------------------------------------------------------------------

def gzip_file(path: Path):
    with open(path, 'rb') as f_in:
        with gzip.open(str(path) + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    path.unlink()

def smart_copy_gzip(src: Path, dst: Path):
    if src.suffix == ".gz":
        shutil.copy2(src, dst)
    else:
        with open(src, 'rb') as f_in, gzip.open(dst, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)




# -------------------------------------------------------------------------
# Model Run utils
# -------------------------------------------------------------------------
# src/modules/utils.py

def render_text_template(
    template_path: Path,
    output_path: Path,
    replacements: Dict[str, str],
) -> Path:
    """
    Replace placeholders in a text template file and write to output.
    """
    text = template_path.read_text()
    for old, new in replacements.items():
        text = text.replace(old, new)
    output_path.write_text(text)
    return output_path


def get_per_dataset_param(
    model_cfg: Dict[str, Any],
    dataset_cfg: Dict[str, Any],
    key: str,
    default_key: str | None = None,
    default_value=None,
):
    """
    Get parameter value with per-dataset override from configs.
    """
    ds_id = dataset_cfg.get("dataset_id")
    per_ds = model_cfg.get("per_dataset", {})

    if ds_id in per_ds and key in per_ds[ds_id]:
        return per_ds[ds_id][key]

    if default_key and default_key in model_cfg:
        return model_cfg[default_key]
    if key in model_cfg:
        return model_cfg[key]

    return default_value