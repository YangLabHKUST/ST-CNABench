from pathlib import Path

import pytest

from st_cnabench.cli import build_parser
from st_cnabench.config import (
    ConfigValidationError,
    load_data_config,
    load_eval_config,
    load_model_config,
)

from tests.helpers import get_repo_root


def test_template_and_example_configs_load() -> None:
    repo_root = get_repo_root()

    load_data_config(repo_root / "configs" / "templates" / "data.template.yaml")
    load_model_config(repo_root / "configs" / "templates" / "models.template.yaml")
    load_eval_config(repo_root / "configs" / "templates" / "eval.template.yaml")

    load_data_config(repo_root / "configs" / "examples" / "cscc_demo" / "data.yaml")
    load_model_config(repo_root / "configs" / "examples" / "cscc_demo" / "models.yaml")
    load_eval_config(repo_root / "configs" / "examples" / "cscc_demo" / "eval.yaml")


def test_invalid_data_config_fails_loudly(tmp_path: Path) -> None:
    bad_cfg = tmp_path / "bad_data.yaml"
    bad_cfg.write_text(
        "\n".join(
            [
                "project_settings:",
                "  dataset_root: ./data",
                "  output_root: ./out",
                "broken_dataset:",
                "  dataset_id: broken",
                "  platform: Visium",
                "  format: SpaceRanger",
                "  genome: hg38",
                "  species: human",
                "  ref_norm: true",
                "  tumor_normal_mode: subset",
                "  raw:",
                "    root: ./data/raw",
                "    tumor_normal: ./data/meta.tsv",
                "  output: {}",
            ]
        )
        + "\n"
    )

    with pytest.raises(ConfigValidationError, match="output is missing required key"):
        load_data_config(bad_cfg)


def test_de_novo_data_config_allows_missing_tumor_normal(tmp_path: Path) -> None:
    cfg = tmp_path / "de_novo_data.yaml"
    cfg.write_text(
        "\n".join(
            [
                "project_settings:",
                "  dataset_root: ./data",
                "  output_root: ./out",
                "de_novo_dataset:",
                "  dataset_id: de_novo_dataset",
                "  platform: Visium",
                "  format: SpaceRanger",
                "  genome: hg38",
                "  species: human",
                "  ref_norm: false",
                "  tumor_normal_mode: de_novo",
                "  raw:",
                "    root: ./data/raw",
                "  output:",
                "    root: ./out/de_novo_dataset",
            ]
        )
        + "\n"
    )

    loaded = load_data_config(cfg)
    assert loaded["de_novo_dataset"]["raw"]["root"] == "./data/raw"


def test_reference_mode_requires_tumor_normal(tmp_path: Path) -> None:
    cfg = tmp_path / "missing_tn_data.yaml"
    cfg.write_text(
        "\n".join(
            [
                "project_settings:",
                "  dataset_root: ./data",
                "  output_root: ./out",
                "subset_dataset:",
                "  dataset_id: subset_dataset",
                "  platform: Visium",
                "  format: SpaceRanger",
                "  genome: hg38",
                "  species: human",
                "  ref_norm: true",
                "  tumor_normal_mode: subset",
                "  raw:",
                "    root: ./data/raw",
                "  output:",
                "    root: ./out/subset_dataset",
            ]
        )
        + "\n"
    )

    with pytest.raises(ConfigValidationError, match="raw.tumor_normal is required"):
        load_data_config(cfg)


def test_cli_exposes_public_entrypoint_shape() -> None:
    parser = build_parser()
    args = parser.parse_args(
        [
            "--steps",
            "prep",
            "run",
            "eval",
            "--data-config",
            "configs/examples/cscc_demo/data.yaml",
            "--model-config",
            "configs/examples/cscc_demo/models.yaml",
            "--eval-config",
            "configs/examples/cscc_demo/eval.yaml",
        ]
    )
    assert args.steps == ["prep", "run", "eval"]
    assert args.data_config.endswith("data.yaml")
