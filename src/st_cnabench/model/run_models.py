# src/run_models.py

import logging
from pathlib import Path
from typing import List, Optional
from ..config import load_model_config
from ..dataset.DatasetPreparator import DatasetPreparator

# Model registry
from .tools_scripts import MODEL_REGISTRY

def run_all_models(
    dataset_cfg_path: str,
    model_cfg_path: str,              
    prep_ids: Optional[List[str]] = None,
    model_names: Optional[List[str]] = None,
    overwrite: bool = False,
    exec_mode: Optional[str] = None,
    verbose: bool = False,
):
    """
    Run Step Core API
    Run specified models on specified datasets.
    Args:
        dataset_cfg_path (str): Path to dataset configuration YAML file.
        model_cfg_path (str): Path to model configuration YAML file.
        prep_ids (Optional[List[str]]): List of dataset IDs to process. If None, process all.
        model_names (Optional[List[str]]): List of model names to run. If None, run all.
        overwrite (bool): Whether to overwrite existing results.
    """

    project_root = Path(".").resolve()
    logging.info(f"Loading global model config: {model_cfg_path}")
    resolved_model_cfg = load_model_config(model_cfg_path)
    output_dir = Path(resolved_model_cfg["project_settings"].get("results_dir")).resolve()
    preparator = DatasetPreparator(dataset_cfg_path)
    all_datasets = preparator.list_datasets()
    # Run mode (docker/conda)
    if not exec_mode:
        exec_mode = resolved_model_cfg["project_settings"].get("default_exec_mode", "conda")
    logging.info(f"Execution Mode: {exec_mode}")

    if prep_ids is None:
        target_datasets = all_datasets
    else:
        target_datasets = []
        for dataset_id in prep_ids:
            if dataset_id not in all_datasets:
                logging.warning(f"Dataset ID [{dataset_id}] not found in config. Skipping.")
                continue
            target_datasets.append(dataset_id)
    if not target_datasets:
        logging.error("No valid datasets selected to run.")
        return
    logging.info(f"Datasets to process: {target_datasets}")

    # target models
    available_models = [name for name in MODEL_REGISTRY if name in resolved_model_cfg]
    
    if not model_names:
        target_models = available_models
    else:
        target_models = [m for m in model_names if m in MODEL_REGISTRY and m in resolved_model_cfg]
        invalid_models = set(model_names) - set(available_models)
        if invalid_models:
            logging.warning(f"Skipping unknown or unconfigured models: {invalid_models}")
    if not target_models:
        logging.error("No valid models selected to run.")
        return
    logging.info(f"Models to run: {target_models}")

   # Initialize models
    initialized_models = []
    for m_name in target_models:
        model_cls = MODEL_REGISTRY[m_name]
        model_instance = model_cls(
            project_root=str(project_root),
            model_cfg=resolved_model_cfg,
            result_dir=output_dir,
            exec_mode=exec_mode
        )
        initialized_models.append(model_instance)
    
    # Main execution loop
    for ds_id in target_datasets:
        ds_cfg = preparator.dataset_cfgs.get(ds_id)
        
        logging.info(f"--- Processing Dataset: {ds_id} ---")
        for model in initialized_models:
            model.run(ds_cfg, overwrite=overwrite,verbose=verbose)

    logging.info("All model executions finished.")
