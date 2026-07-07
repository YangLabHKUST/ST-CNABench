import logging

def run_extract_data(eval_list, gt_loader, ds_tn_annot_path=None, ds_tn_annot_path_model_run=None, ds_tn_mode=None):
    '''
    Extract tumor/normal prediction data from each loader in eval_list.
    '''
    if ds_tn_mode not in {"subset", "de_novo"}:
        raise ValueError(f"tumor_normal_mode must be 'subset' or 'de_novo', got: {ds_tn_mode}")
    if not ds_tn_annot_path:
        raise ValueError("tumor_normal_gt path is required for tumor_normal evaluation.")
    if ds_tn_mode == "subset" and not ds_tn_annot_path_model_run:
        raise ValueError("tumor_normal path is required when tumor_normal_mode='subset'.")

    results = {}
    # Extract GT
    if ds_tn_mode == "subset":
        logging.info("Extracting GT tumor normal preds in subset mode")
        results['GT'] = gt_loader.extract_tumor_normal_preds(ds_tn_annot_path, ds_tn_annot_path_model_run)
    else:
        logging.info("Extracting GT tumor normal preds in de novo mode")
        results['GT'] = gt_loader.extract_tumor_normal_preds(ds_tn_annot_path)

    for loader in eval_list:
        model_name = loader.model_name
        if ds_tn_mode == "subset":
            logging.info(f"Extracting sub normal tumor normal preds in subset mode for model={model_name}")
            results[model_name] = loader.extract_sub_normal_tumor_normal_preds(ds_tn_annot_path_model_run)
            if results[model_name] is None:
                logging.warning(f"Model {model_name} does not support sub normal tumor/normal extraction, dropping from results.")
                results.pop(model_name)
        else:
            logging.info(f"Extracting tumor normal preds in de novo mode for model={model_name}")
            results[model_name] = loader.extract_tumor_normal_preds()
            if results[model_name] is None:
                logging.warning(f"Model {model_name} does not support tumor/normal extraction, dropping from results.")
                results.pop(model_name)

    return results
