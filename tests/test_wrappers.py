from pathlib import Path

from st_cnabench.config import load_model_config
from st_cnabench.dataset import DatasetPreparator, run_data_preparation
from st_cnabench.model.tools_scripts import MODEL_REGISTRY

from tests.helpers import create_support_tree, write_demo_data_config, write_full_model_config


def test_all_wrappers_prepare_inputs_and_build_commands(tmp_path: Path) -> None:
    data_cfg = write_demo_data_config(tmp_path)
    run_data_preparation(data_cfg, dataset_ids=["P6_vis_rep1"], overwrite=True)

    results_dir = tmp_path / "raw_results"
    support = create_support_tree(tmp_path)
    model_cfg_path = write_full_model_config(tmp_path, results_dir, support)
    model_cfg = load_model_config(model_cfg_path)

    preparator = DatasetPreparator(data_cfg)
    dataset_cfg = preparator.dataset_cfgs["P6_vis_rep1"]

    for model_name, model_cls in MODEL_REGISTRY.items():
        model_output = results_dir / dataset_cfg["dataset_id"] / model_name
        model_output.mkdir(parents=True, exist_ok=True)

        model = model_cls(
            project_root=str(tmp_path),
            model_cfg=model_cfg,
            result_dir=str(results_dir),
            exec_mode="docker",
        )
        input_files = model.prepare_inputs(dataset_cfg, model_output, overwrite=True)
        command = model.build_command(dataset_cfg, input_files, model_output)

        assert command, f"{model_name} returned an empty command"
        assert all(str(part) for part in command), f"{model_name} returned an invalid command"
