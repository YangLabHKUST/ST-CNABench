"""ST-CNVBench public package."""

from importlib import import_module
from typing import Any

__version__ = "0.1.0"

_LAZY_ATTRS = {
    "DatasetPreparator": ("st_cnvbench.dataset", "DatasetPreparator"),
    "run_all_models": ("st_cnvbench.model", "run_all_models"),
    "run_data_preparation": ("st_cnvbench.dataset", "run_data_preparation"),
    "run_evaluation": ("st_cnvbench.evaluate", "run_evaluation"),
}

__all__ = [*_LAZY_ATTRS, "__version__"]


def __getattr__(name: str) -> Any:
    if name not in _LAZY_ATTRS:
        raise AttributeError(f"module 'st_cnvbench' has no attribute {name!r}")

    module_name, attr_name = _LAZY_ATTRS[name]
    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value
    return value
