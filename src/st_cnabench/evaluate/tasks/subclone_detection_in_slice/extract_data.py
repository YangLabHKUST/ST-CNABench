import pandas as pd
import os
import logging

def run_extract_data(
    eval_list,
    gt_loader,
    gene_annot_path=None,
    cna_gt_path=None,
    ds_subcluster_annot_path=None,
    tn_annot_path_model_run=None,
    beads_mapping_path=None,
    ds_spatial_coord=None,
):
    """Extract in-slice subclone labels, clone CNA profiles, and optional coordinates."""

    results = {}
    clone_cna_profiles = {}

    # Load optional spot-level GT subclone labels.
    if ds_subcluster_annot_path is not None:
        results['GT'] = gt_loader.extract_subcluster_preds(ds_subcluster_annot_path)

    # Load optional GT clone-level CNA profiles.
    if cna_gt_path is not None:
        logging.info(f"Extracting clone CNA profiles for GT from {cna_gt_path}")
        gt_cna_profile = gt_loader.extract_clone_cna_profile(cna_gt_path=cna_gt_path)
        if gt_cna_profile is not None and not gt_cna_profile.empty:
            clone_cna_profiles['GT'] = gt_cna_profile
        else:
            logging.warning("[GT] No clone CNA profile extracted.")

    # Collect model predictions and clone-level CNA profiles.
    for loader in eval_list:
        model_name = loader.model_name
        logging.info(f"Extracting subcluster preds and clone CNA profiles for model={model_name}")

        preds = loader.extract_subcluster_preds()
        if preds is None:
            logging.warning(f"Model {model_name} does not support subcluster extraction, dropping from results.")
            continue
        results[model_name] = preds

        cna_profile = loader.extract_clone_cna_profile(gene_annot_path=gene_annot_path, tn_annot_path=tn_annot_path_model_run)
        if cna_profile is not None and not cna_profile.empty:
            clone_cna_profiles[model_name] = cna_profile
        else:
            logging.warning(f"[{model_name}] No clone CNA profile extracted.")

    # Merge labels onto a barcode union so uncalled spots remain visible.
    df_all = None
    if len(results) > 0:
        for model_name, df in results.items():
            required_cols = {'Barcodes', 'Label_preds'}
            missing_cols = required_cols - set(df.columns)
            if missing_cols:
                raise ValueError(f"[{model_name}] subclone prediction table missing columns: {sorted(missing_cols)}")
            df = df.rename(columns={'Label_preds': f'{model_name}_Preds'})
            if df_all is None:
                df_all = df
            else:
                df_all = pd.merge(df_all, df, on='Barcodes', how='outer')

    # Remove reference normal cells when the model-run annotation is available.
    if tn_annot_path_model_run is not None and df_all is not None:
        df_model_run = pd.read_csv(tn_annot_path_model_run, sep='\t', header=0, comment='#')
        df_model_run.rename(columns={df_model_run.columns[0]: 'Barcodes', df_model_run.columns[1]: 'tumor_normal'}, inplace=True)
        normal_mask = df_model_run['tumor_normal'].astype(str).str.strip().str.lower() == 'normal'
        normal_barcode = set(df_model_run[normal_mask]['Barcodes'])
        df_all = df_all[~df_all['Barcodes'].isin(normal_barcode)]

    # Broadcast meta-bead predictions back to raw beads for Slide-DNA-seq style inputs.
    if beads_mapping_path is not None and df_all is not None:
        df_beads_mapping = pd.read_csv(beads_mapping_path, header=0)
        required_mapping_cols = {'pseudo_barcode', 'original_barcode'}
        missing_mapping_cols = required_mapping_cols - set(df_beads_mapping.columns)
        if missing_mapping_cols:
            raise ValueError(f"Beads mapping file missing columns: {sorted(missing_mapping_cols)}")
        df_all = df_all.merge(df_beads_mapping, left_on='Barcodes', right_on='pseudo_barcode', how='left')
        mapped_mask = df_all['original_barcode'].notna()
        df_all.loc[mapped_mask, 'Barcodes'] = df_all.loc[mapped_mask, 'original_barcode']

    # Attach spatial coordinates when a standardized coordinate file is available.
    if ds_spatial_coord is not None and df_all is not None and os.path.exists(ds_spatial_coord):
        df_coords = pd.read_csv(ds_spatial_coord, header=0)
        rename_map = {}
        if 'barcode' in df_coords.columns:
            rename_map['barcode'] = 'Barcodes'
        if 'pxl_col_in_fullres' in df_coords.columns:
            rename_map['pxl_col_in_fullres'] = 'x'
        if 'pxl_row_in_fullres' in df_coords.columns:
            rename_map['pxl_row_in_fullres'] = 'y'
        df_coords = df_coords.rename(columns=rename_map)
        required_coord_cols = {'Barcodes', 'x', 'y'}
        missing_coord_cols = required_coord_cols - set(df_coords.columns)
        if missing_coord_cols:
            raise ValueError(f"Spatial coordinate file missing columns: {sorted(missing_coord_cols)}")
        df_all = df_all.merge(df_coords[['Barcodes', 'x', 'y']], on='Barcodes', how='left')

    # Normalize duplicate coordinate columns after successive merges.
    if df_all is not None:
        x_candidates = [c for c in ['x', 'x_x', 'x_y'] if c in df_all.columns]
        y_candidates = [c for c in ['y', 'y_x', 'y_y'] if c in df_all.columns]

        if x_candidates:
            x_series = None
            for col in x_candidates:
                x_series = df_all[col] if x_series is None else x_series.combine_first(df_all[col])
            df_all['x'] = x_series

        if y_candidates:
            y_series = None
            for col in y_candidates:
                y_series = df_all[col] if y_series is None else y_series.combine_first(df_all[col])
            df_all['y'] = y_series

        redundant_coord_cols = [
            col for col in ['x_x', 'x_y', 'y_x', 'y_y']
            if col in df_all.columns and col not in {'x', 'y'}
        ]
        if redundant_coord_cols:
            df_all = df_all.drop(columns=redundant_coord_cols)

    return df_all, clone_cna_profiles
