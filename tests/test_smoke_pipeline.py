from pathlib import Path

from st_cnvbench.dataset import run_data_preparation
from st_cnvbench.evaluate import run_evaluation
from st_cnvbench.model import run_all_models
from st_cnvbench.model.base import BaseModel

from tests.helpers import (
    write_demo_data_config,
    write_smoke_eval_config,
    write_smoke_model_config,
)


def _fake_execute_in_env(self, cmd, cwd, log_file, verbose=False):  # noqa: ANN001
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("fake wrapper execution\n")

    perf_file = log_file.with_suffix(".perf")
    perf_file.write_text(
        "\n".join(
            [
                "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.00",
                "Maximum resident set size (kbytes): 1024",
                "Percent of CPU this job got: 100%",
                "Exit status: 0",
            ]
        )
        + "\n"
    )


def test_smoke_prep_run_eval_pipeline(tmp_path: Path, monkeypatch) -> None:
    data_cfg = write_demo_data_config(tmp_path)
    run_data_preparation(data_cfg, dataset_ids=["P6_vis_rep1"], overwrite=True)

    results_dir = tmp_path / "raw_results"
    eval_dir = tmp_path / "evaluation"
    model_cfg = write_smoke_model_config(tmp_path, results_dir)
    eval_cfg = write_smoke_eval_config(tmp_path, results_dir, eval_dir)

    monkeypatch.setattr(BaseModel, "_execute_in_env", _fake_execute_in_env)

    run_all_models(
        dataset_cfg_path=str(data_cfg),
        model_cfg_path=str(model_cfg),
        prep_ids=["P6_vis_rep1"],
        model_names=["CopyKAT"],
        overwrite=True,
        exec_mode="conda",
        verbose=False,
    )

    run_evaluation(
        dataset_cfg_path=str(data_cfg),
        eval_cfg_path=str(eval_cfg),
        target_ds_ids=["P6_vis_rep1"],
        model_names=["CopyKAT"],
        tasks=["efficiency"],
    )

    summary_tsv = (
        eval_dir
        / "P6_vis_rep1"
        / "computational_efficiency"
        / "computational_efficiency_summary.tsv"
    )
    assert summary_tsv.exists()
    assert "CopyKAT" in summary_tsv.read_text()
