import pandas as pd
import os
import logging


def run_extract_data(eval_list, gene_annot_path, ds_tn_annot_path, ds_spatial_coord=None):
    """
    Extract the inputs needed for clonal-evolution evaluation.

    Returns:
      1. df_all: spot-level clone labels plus spatial coordinates.
      2. clone_cnv_profiles: clone-level CNV signals for tree inference.
    """
    results = {}
    clone_cnv_profiles = {}

    for loader in eval_list:
        model_name = loader.model_name
        logging.info(f"Extracting subcluster preds and clone CNV profiles for model={model_name}")

        preds = loader.extract_subcluster_preds()
        if preds is not None and not preds.empty:
            results[model_name] = preds

            try:
                cnv_profile = loader.extract_clone_cnv_profile(gene_annot_path, ds_tn_annot_path)
                if cnv_profile is not None and not cnv_profile.empty:
                    clone_cnv_profiles[f'{model_name}_Preds'] = cnv_profile
            except Exception as e:
                logging.warning(f"Failed to extract clone CNV profile for {model_name}: {e}")

        else:
            logging.warning(f"Model {model_name} does not support subcluster extraction or failed, dropping from results.")

    # Merge all spot-level prediction tables.
    df_all = None
    if len(results) > 0:
        for model_name, df in results.items():
            df = df.rename(columns={'Label_preds': f'{model_name}_Preds'})
            if df_all is None:
                df_all = df
            else:
                df_all = pd.merge(df_all, df, on='Barcodes', how='inner')

    # Attach physical array coordinates when available.
    if ds_spatial_coord is not None and os.path.exists(ds_spatial_coord):
        df_coords = pd.read_csv(ds_spatial_coord, header=0)
        df_coords = df_coords.rename(columns={
            'barcode': 'Barcodes',
            'array_row': 'y',
            'array_col': 'x'
        })

        df_coords = df_coords.drop(columns=['in_tissue', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], errors='ignore')

        if 'x' in df_coords.columns:
            df_coords['x'] = df_coords['x'].astype(float)
        if 'y' in df_coords.columns:
            df_coords['y'] = df_coords['y'].astype(float)

        if df_all is not None:
            df_all = df_all.merge(df_coords, on='Barcodes', how='left')
        else:
            df_all = df_coords

    return df_all, clone_cnv_profiles
