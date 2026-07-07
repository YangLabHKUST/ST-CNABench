"""Evaluation entrypoints."""


def run_evaluation(*args, **kwargs):
    from .run_eval import run_evaluation as _run_evaluation

    return _run_evaluation(*args, **kwargs)


__all__ = ["run_evaluation"]
