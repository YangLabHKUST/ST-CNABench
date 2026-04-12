import logging
from pathlib import Path
from ..config import load_eval_config
from ..dataset.DatasetPreparator import DatasetPreparator
from .loaders import LOADER_REGISTRY
from .loaders.gt import GTLoader
from .loaders import EFFICIENCY_LOADER_LIST
from .loaders import RESOLUTION_LOADER_LIST
from .loaders import TN_PRED_LOADER_LIST
from .loaders import SUBCLONE_PRED_LOADER_LIST
from .run_cnv_profile import run_cnv_profile
from .run_computational_efficiency import run_computational_efficiency 
from .run_cnv_resolution import run_cnv_resolution
from .run_tumor_normal import run_tumor_normal
from .run_subclone_detection_in_slice import run_subclone_detection_in_slice
from .run_subclone_detection_organ import run_subclone_detection_organ
from .run_clone_evolution import run_clone_evolution
def run_evaluation(dataset_cfg_path, eval_cfg_path, target_ds_ids=None, model_names=None, tasks=None):
    """
    Evaluation Wrapper Function
    """
    # Config
    eval_cfg = load_eval_config(eval_cfg_path)
    preparator = DatasetPreparator(dataset_cfg_path)
    
    results_root = Path(eval_cfg['project_settings']['results_dir'])
    eval_root = Path(eval_cfg['project_settings']['eval_dir'])
    
    # Global params
    global_params = eval_cfg.get('global_params', {})
    ref_genome = global_params.get('genome_version', 'hg38')
    bin_size = global_params.get('bin_size', 100000)
    gene_annot_path = global_params.get('gene_annot_path')

    if tasks is None:
        tasks = ['efficiency','resolution','tumor_normal','subclone_detection_in_slice', 'subclone_detection_organ','clonal_evolution', 'cnv_profile']
    else:
        tasks = [t for t in tasks if t in ['efficiency','resolution','tumor_normal','subclone_detection_in_slice','subclone_detection_organ','clonal_evolution', 'cnv_profile']]
    logging.info(f"Tasks to run: {tasks}")

    # Target datasets and models
    all_datasets = preparator.list_datasets()
    target_datasets = []
    # target datasets
    if target_ds_ids:
        for dataset_id in target_ds_ids:
            if dataset_id not in all_datasets:
                logging.warning(f"Dataset ID [{dataset_id}] not found in config. Skipping.")
                continue
            target_datasets.append(dataset_id)
    else:
        target_datasets = all_datasets
    logging.info(f"Datasets to process: {target_datasets}")
    eval_list_cfg = eval_cfg.get('eval_list', {})
    target_models = [m for m in model_names if m in eval_list_cfg] if model_names else list(eval_list_cfg.keys())
    logging.info(f"Models to evaluate: {target_models}")
    for ds_id in target_datasets:
        ds_cfg = preparator.dataset_cfgs.get(ds_id)
        ds_platform = ds_cfg.get('platform')
        # Explicit tumor/normal annotations:
        # - tumor_normal_gt: real GT labels for tumor/normal task
        # - tumor_normal: model-run annotation (subset normal reference)
        ds_tn_annot_path = ds_cfg['raw'].get('tumor_normal_gt')
        ds_tn_annot_path_model_run = ds_cfg['raw'].get('tumor_normal')
        ds_tn_mode = ds_cfg.get('tumor_normal_mode')
        ds_ref_norm = ds_cfg.get('ref_norm')
        ds_spatial_coord = Path(ds_cfg['output']['root']) / 'spatial' / 'tissue_positions.csv'
        ds_scalefactors = Path(ds_cfg['output']['root']) / 'spatial' / 'scalefactors_json.json'
        ds_HE_img = Path(ds_cfg['output']['root']) / 'spatial' / 'tissue_hires_image.png'
        logging.info(f"=== Process Dataset: {ds_id} ===")
        
        dataset_eval_loaders = []
        
        # Init Model loaders
        for model_name in target_models:
            methods = eval_list_cfg[model_name].get('eval_name', [])
            model_res_dir = results_root / ds_id / model_name
            
            if not model_res_dir.exists():
                logging.warning(f"  [Skip] Result directory not found: {model_res_dir}")
                continue
            
            for eval_name in methods:
                loader_cls = LOADER_REGISTRY.get(eval_name)
                if not loader_cls:
                    logging.error(f"  ! Loader '{eval_name}' not found")
                    continue
                
                save_dir = eval_root / ds_id / eval_name
                save_dir.mkdir(parents=True, exist_ok=True)
                
                loader_instance = loader_cls(
                    eval_name=eval_name,
                    task_name="evaluation", 
                    result_dir=str(model_res_dir),
                    save_dir=str(save_dir)
                )
                dataset_eval_loaders.append(loader_instance)

        if not dataset_eval_loaders:
            logging.warning(f"No valid loaders for dataset {ds_id}, skipping.")
            continue
        # Eval Tasks
        if 'efficiency' in tasks:
            logging.info(f"  [Run] Computational Efficiency Task for {ds_id}...")
            # Use one loader for each model
            efficiency_loaders = [loader for loader in dataset_eval_loaders if loader.eval_name in EFFICIENCY_LOADER_LIST]
            for loader in efficiency_loaders:
                loader.task_name = "computational_efficiency"
            run_computational_efficiency(efficiency_loaders)
        if 'resolution' in tasks:
            logging.info(f"  [Run] CNV Resolution Task for {ds_id}...")
            # Use one loader for each model
            resolution_loaders = [loader for loader in dataset_eval_loaders if loader.eval_name in RESOLUTION_LOADER_LIST]
            for loader in resolution_loaders:
                loader.task_name = "cnv_resolution"
            logging.info(f"  [Run] CNV Resolution Task for {ds_id}...")
            run_cnv_resolution(resolution_loaders, gene_annot_path, ds_tn_annot_path_model_run)
        if 'tumor_normal' in tasks:
            logging.info(f"  [Run] Tumor Normal Prediction Task for {ds_id}...")          
            if ds_tn_mode not in {"subset", "de_novo", "off"}:
                raise ValueError(
                    f"[{ds_id}] tumor_normal_mode must be one of 'subset'/'de_novo'/'off', "
                    f"got: {ds_tn_mode}"
                )
            if ds_tn_mode == "off":
                logging.info(f"  [Skip] Tumor Normal Prediction Task disabled by tumor_normal_mode='off' for {ds_id}.")
            else:
                if not ds_tn_annot_path:
                    raise ValueError(f"[{ds_id}] raw.tumor_normal_gt is required for tumor_normal evaluation.")
                if ds_tn_mode == "subset" and not ds_tn_annot_path_model_run:
                    raise ValueError(f"[{ds_id}] raw.tumor_normal is required when tumor_normal_mode='subset'.")
                if ds_tn_mode == "subset" and ds_ref_norm is not True:
                    raise ValueError(f"[{ds_id}] tumor_normal_mode='subset' requires ref_norm=True.")

                gt_save_dir = eval_root / ds_id / "GT"
                gt_save_dir.mkdir(parents=True, exist_ok=True)
                gt_loader = GTLoader(
                    task_name="tumor_normal_prediction",
                    eval_name='GT',
                    result_dir=str(Path(ds_tn_annot_path_model_run).parent),
                    save_dir=str(gt_save_dir)
                )

                tn_preds_loaders = [loader for loader in dataset_eval_loaders if loader.eval_name in TN_PRED_LOADER_LIST]

                for loader in tn_preds_loaders:
                    loader.task_name = "tumor_normal_prediction"
                run_tumor_normal(
                    eval_list=tn_preds_loaders,
                    gt_loader=gt_loader,
                    ds_tn_annot_path=ds_tn_annot_path,
                    ds_tn_annot_path_model_run=ds_tn_annot_path_model_run,
                    ds_tn_mode=ds_tn_mode,
                    ds_spatial_coord=ds_spatial_coord,
                    ds_scalefactors=ds_scalefactors,
                    ds_HE_img=ds_HE_img,
                    ds_platform=ds_platform,
                    threshold=2.1,
                )
        if 'subclone_detection_in_slice' in tasks:
            logging.info(f"  [Run] Subclone Detection Task for {ds_id}...")
            gt_path = ds_cfg['raw'].get('subclone_gt')
            cnv_gt_path = ds_cfg['raw'].get('cnv_gt')                
            gt_save_dir = eval_root / ds_id / "GT"
            gt_save_dir.mkdir(parents=True, exist_ok=True)
            gt_loader = GTLoader(
                task_name="subclone_detection_in_slice",
                eval_name='GT',
                result_dir=str(Path(gt_path).parent),
                save_dir=str(gt_save_dir)
             )
            for loader in dataset_eval_loaders:
                loader.task_name = "subclone_detection_in_slice"
            subclone_preds_loaders = [loader for loader in dataset_eval_loaders if loader.eval_name in SUBCLONE_PRED_LOADER_LIST]
            run_subclone_detection_in_slice(
                eval_list=subclone_preds_loaders,
                gt_loader=gt_loader,
                cnv_gt_path=cnv_gt_path,
                gene_annot_path=gene_annot_path,
                ds_subcluster_annot_path=gt_path,
                tn_annot_path_model_run=ds_tn_annot_path_model_run,
                ds_spatial_coord=ds_spatial_coord,
                ds_scalefactors=ds_scalefactors,
                ds_HE_img=ds_HE_img,
                threshold=2.1,
                beads_mapping_path=ds_cfg['raw'].get('beads_mapping'),
                ds_platform=ds_platform,
            )
        
        if 'subclone_detection_organ' in tasks:
            logging.info(f"  [Run] Subclone Detection Task for {ds_id}...")
            gt_path = ds_cfg['raw'].get('subclone_gt')                
            gt_save_dir = eval_root / ds_id / "GT"
            gt_save_dir.mkdir(parents=True, exist_ok=True)
            gt_loader = GTLoader(
                task_name="subclone_detection_organ",
                eval_name='GT',
                result_dir=str(Path(gt_path).parent),
                save_dir=str(gt_save_dir)
             )
            for loader in dataset_eval_loaders:
                loader.task_name = "subclone_detection_organ"
            subclone_preds_loaders = [loader for loader in dataset_eval_loaders if loader.eval_name in SUBCLONE_PRED_LOADER_LIST]
            run_subclone_detection_organ(
                eval_list=subclone_preds_loaders,
                gt_loader=gt_loader,
                gene_annot_path=gene_annot_path,
                ds_tn_annot_path=ds_tn_annot_path_model_run,
                ds_subcluster_annot_path=gt_path,
                tn_annot_path_model_run=ds_tn_annot_path_model_run,
                ds_spatial_coord=ds_spatial_coord,
                ds_scalefactors=ds_scalefactors,
                ds_HE_img=ds_HE_img,
                threshold=2.1,
                beads_mapping_path=ds_cfg['raw'].get('beads_mapping')
            )
        if 'clonal_evolution' in tasks:
            logging.info(f"  [Run] Clonal Evolution Task for {ds_id}...")
            subclone_preds_loaders = [loader for loader in dataset_eval_loaders if loader.eval_name in SUBCLONE_PRED_LOADER_LIST]
            run_clone_evolution(
                eval_list=subclone_preds_loaders,
                gene_annot_path=gene_annot_path,
                ds_tn_annot_path=ds_tn_annot_path_model_run,
                ds_spatial_coord=ds_spatial_coord,
                ref_genome=ref_genome,
                bin_size=bin_size
            )
        if 'cnv_profile' in tasks:
            gt_path = ds_cfg['raw'].get('cnv_gt')                
            gt_save_dir = eval_root / ds_id / "GT"
            gt_save_dir.mkdir(parents=True, exist_ok=True)
            gt_loader = GTLoader(
                task_name="cnv_profile",
                eval_name='GT',
                result_dir=str(Path(gt_path).parent),
                save_dir=str(gt_save_dir)
            )

            for loader in dataset_eval_loaders:
                loader.task_name = "cnv_profile"
            
            # Use all loaders
            run_cnv_profile(
                eval_list=dataset_eval_loaders,
                gt_loader=gt_loader,
                gt_path=gt_path,
                gene_annot_path=gene_annot_path,
                ds_tn_annot_path=ds_tn_annot_path_model_run,
                ref_genome=ref_genome,
                bin_size=bin_size
            )

    logging.info("Evaluation Pipeline Finished.")
