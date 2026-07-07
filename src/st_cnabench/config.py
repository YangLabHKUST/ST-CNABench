from __future__ import annotations

from pathlib import Path
from typing import Any

from .modules.utils import load_config


class ConfigValidationError(ValueError):
    """Raised when a public config file is incomplete or invalid."""


def _require_mapping(value: Any, path: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigValidationError(f"{path} must be a mapping.")
    return value


def _require_keys(mapping: dict[str, Any], path: str, keys: list[str]) -> None:
    missing = [key for key in keys if key not in mapping]
    if missing:
        joined = ", ".join(missing)
        raise ConfigValidationError(f"{path} is missing required key(s): {joined}.")


def _require_choice(path: str, value: Any, allowed: set[str]) -> None:
    if value not in allowed:
        choices = ", ".join(sorted(allowed))
        raise ConfigValidationError(f"{path} must be one of: {choices}. Got: {value!r}.")


def load_data_config(config_path: str | Path) -> dict[str, Any]:
    config = load_config(str(config_path))
    validate_data_config(config)
    return config


def load_model_config(config_path: str | Path) -> dict[str, Any]:
    config = load_config(str(config_path))
    validate_model_config(config)
    return config


def load_eval_config(config_path: str | Path) -> dict[str, Any]:
    config = load_config(str(config_path))
    validate_eval_config(config)
    return config


def validate_data_config(config: dict[str, Any]) -> None:
    settings = _require_mapping(config.get("project_settings"), "project_settings")
    _require_keys(settings, "project_settings", ["dataset_root", "output_root"])

    dataset_sections = {
        key: value for key, value in config.items() if key != "project_settings"
    }
    if not dataset_sections:
        raise ConfigValidationError("Data config must define at least one dataset entry.")

    for name, dataset_cfg in dataset_sections.items():
        ds_path = f"dataset[{name}]"
        ds = _require_mapping(dataset_cfg, ds_path)
        _require_keys(
            ds,
            ds_path,
            [
                "dataset_id",
                "platform",
                "format",
                "genome",
                "species",
                "ref_norm",
                "tumor_normal_mode",
                "raw",
                "output",
            ],
        )
        _require_choice(
            f"{ds_path}.format",
            str(ds["format"]).lower(),
            {"spaceranger", "stpipeline"},
        )
        _require_choice(
            f"{ds_path}.tumor_normal_mode",
            ds["tumor_normal_mode"],
            {"subset", "de_novo", "off"},
        )

        raw = _require_mapping(ds["raw"], f"{ds_path}.raw")
        output = _require_mapping(ds["output"], f"{ds_path}.output")
        _require_keys(output, f"{ds_path}.output", ["root"])

        needs_tumor_normal = (
            ds["tumor_normal_mode"] == "subset" or ds.get("ref_norm") is True
        )

        if str(ds["format"]).lower() == "spaceranger":
            _require_keys(raw, f"{ds_path}.raw", ["root"])
        else:
            _require_keys(
                raw,
                f"{ds_path}.raw",
                ["counts", "barcodes", "features", "coords"],
            )
        if needs_tumor_normal and not raw.get("tumor_normal"):
            raise ConfigValidationError(
                f"{ds_path}.raw.tumor_normal is required when "
                "tumor_normal_mode='subset' or ref_norm=True."
            )


def validate_model_config(config: dict[str, Any]) -> None:
    settings = _require_mapping(config.get("project_settings"), "project_settings")
    _require_keys(
        settings,
        "project_settings",
        ["root_dir", "results_dir", "refs", "exts", "default_exec_mode"],
    )
    _require_choice(
        "project_settings.default_exec_mode",
        settings["default_exec_mode"],
        {"conda", "docker", "apptainer"},
    )

    model_sections = {
        key: value for key, value in config.items() if key != "project_settings"
    }
    if not model_sections:
        raise ConfigValidationError("Model config must define at least one model section.")

    for name, model_cfg in model_sections.items():
        model = _require_mapping(model_cfg, f"model[{name}]")
        if model.get("enabled", True):
            _require_keys(model, f"model[{name}]", ["model_name"])


def validate_eval_config(config: dict[str, Any]) -> None:
    settings = _require_mapping(config.get("project_settings"), "project_settings")
    _require_keys(settings, "project_settings", ["root_dir", "results_dir", "eval_dir"])

    global_params = _require_mapping(config.get("global_params"), "global_params")
    _require_keys(
        global_params,
        "global_params",
        ["bin_size", "genome_version", "gene_annot_path"],
    )

    eval_list = _require_mapping(config.get("eval_list"), "eval_list")
    if not eval_list:
        raise ConfigValidationError("eval_list must define at least one model mapping.")

    for name, section in eval_list.items():
        eval_section = _require_mapping(section, f"eval_list[{name}]")
        _require_keys(eval_section, f"eval_list[{name}]", ["eval_name"])
        eval_names = eval_section["eval_name"]
        if not isinstance(eval_names, list) or not eval_names:
            raise ConfigValidationError(
                f"eval_list[{name}].eval_name must be a non-empty list."
            )
