import logging

def run_extract_data(eval_list, gene_annot_path, ds_tn_annot_path):
    '''
    Extract CNV resolution data from each loader in eval_list.
    Structure:
            results = {
                'ModelA': {
                    'LoaderX': DataFrame_for_ModelA_LoaderX,
                    'LoaderY': DataFrame_for_ModelA_LoaderY,
                    ...
                },
                'ModelB': {
                    'LoaderX': DataFrame_for_ModelB_LoaderX,
                    ...
                },
                ...
            }
    '''
    results = {}
    for loader in eval_list:
        model_name = loader.model_name
        eval_name = loader.eval_name
        logging.info(f"Extracting resolution for model={model_name}, loader={eval_name}")

        try:
            df_resolution = loader.extract_cnv_resolution(gene_annot_path, ds_tn_annot_path)
        except Exception as e:
            logging.warning(
                f"Skipping resolution for model={model_name}, loader={eval_name} due to error: {e}"
            )
            continue

        if df_resolution is None:
            logging.warning(f"Skipping resolution for model={model_name}, loader={eval_name}: returned None")
            continue

        if hasattr(df_resolution, 'empty') and df_resolution.empty:
            logging.warning(f"Skipping resolution for model={model_name}, loader={eval_name}: empty DataFrame")
            continue

        results.setdefault(model_name, {})[eval_name] = df_resolution

    return results
