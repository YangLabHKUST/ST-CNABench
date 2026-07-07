import pandas as pd
import logging

def run_extract_efficiency_data(eval_list):
    """
    Extract computational efficiency data from each loader in eval_list.
    """
    logging.info("Starting extraction of computational efficiency data...")
    
    all_efficiency_stats = []

    for loader in eval_list:
        stats = loader.extract_computational_efficiency_data()
            
        if stats is not None:
            all_efficiency_stats.append(stats)
            logging.info(f"Successfully extracted performance stats for {loader.model_name}")
        else:
            logging.warning(f"No performance data found for {loader.model_name}")

    df_efficiency = pd.DataFrame(all_efficiency_stats)
    
    return df_efficiency