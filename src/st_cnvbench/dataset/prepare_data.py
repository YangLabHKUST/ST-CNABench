import sys
import time
import logging
from pathlib import Path
from .DatasetPreparator import DatasetPreparator

def run_data_preparation(config_path, dataset_ids=None, overwrite=False):
    """
    Core API function for data preparation.
    
    Args:
        config_path (str/Path): Path to the DATASET configuration file.
        dataset_ids (list, optional): List of dataset IDs to process. None processes all.
        overwrite (bool): Whether to overwrite existing output directories.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        logging.error(f"Config file not found at {config_path}")
        return False
    # init
    preparator = DatasetPreparator(str(config_path))

    # target dataset
    available_datasets = preparator.list_datasets()
    if dataset_ids is None:
        target_datasets = available_datasets
    else:
        target_datasets = []
        for dataset_id in dataset_ids:
            if dataset_id not in available_datasets:
                logging.warning(f"Dataset ID [{dataset_id}] not found in config. Skipping.")
                continue
            target_datasets.append(dataset_id)

    logging.info("-" * 60)
    logging.info("Starting Data Preparation Pipeline")
    logging.info(f"Config: {config_path.name}")
    logging.info(f"Targets: {len(target_datasets)} dataset(s)")
    logging.info(f"Overwrite: {overwrite}")
    logging.info("-" * 60)

    success_count = 0
    fail_count = 0
    start_global = time.time()

    # Execute loop
    for ds_id in target_datasets:
        logging.info(f"Processing [{ds_id}]...")
        start_time = time.time()
        preparator.prepare_dataset(ds_id, overwrite=overwrite)
        elapsed = time.time() - start_time
        logging.info(f"  > Done in {elapsed:.2f}s")
        success_count += 1

    # Summary
    total_time = time.time() - start_global
    logging.info("-" * 60)
    logging.info(f"Completed in {total_time:.2f}s")
    logging.info(f"Success: {success_count} | Failed: {fail_count}")
    logging.info("-" * 60)
    
    return True
