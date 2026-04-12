from __future__ import annotations

import argparse
import logging
import sys


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "ST-CNVBench: benchmark CNV inference methods on spatial transcriptomics."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        choices=["prep", "run", "eval"],
        default=["prep", "run", "eval"],
        help="Pipeline steps to execute.",
    )
    parser.add_argument(
        "--exec-mode",
        choices=["conda", "docker", "apptainer"],
        default=None,
        help="Execution mode for model wrappers. Overrides config when provided.",
    )
    parser.add_argument(
        "--data-config",
        default="configs/templates/data.template.yaml",
        help="Path to the dataset config YAML.",
    )
    parser.add_argument(
        "--model-config",
        default="configs/templates/models.template.yaml",
        help="Path to the model config YAML.",
    )
    parser.add_argument(
        "--eval-config",
        default="configs/templates/eval.template.yaml",
        help="Path to the evaluation config YAML.",
    )
    parser.add_argument(
        "--prep-ids",
        nargs="+",
        default=None,
        help="Specific dataset IDs to process.",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help="Specific model names to run or evaluate.",
    )
    parser.add_argument(
        "--eval-tasks",
        nargs="+",
        default=None,
        help="Specific evaluation tasks to run.",
    )
    parser.add_argument(
        "--overwrite-prep",
        action="store_true",
        help="Overwrite prepared datasets when they already exist.",
    )
    parser.add_argument(
        "--overwrite-results",
        action="store_true",
        help="Overwrite model outputs when they already exist.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Stream wrapper logs to stdout during model execution.",
    )
    return parser


def setup_logging() -> None:
    logger = logging.getLogger()
    if logger.handlers:
        return

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    setup_logging()
    logging.info("Starting ST-CNVBench")
    logging.info("Steps to execute: %s", args.steps)

    if "prep" in args.steps:
        from .dataset import run_data_preparation

        success = run_data_preparation(
            config_path=args.data_config,
            dataset_ids=args.prep_ids,
            overwrite=args.overwrite_prep,
        )
        if not success:
            raise RuntimeError("Data preparation failed.")

    if "run" in args.steps:
        from .model import run_all_models

        run_all_models(
            dataset_cfg_path=args.data_config,
            model_cfg_path=args.model_config,
            prep_ids=args.prep_ids,
            model_names=args.models,
            overwrite=args.overwrite_results,
            exec_mode=args.exec_mode,
            verbose=args.verbose,
        )

    if "eval" in args.steps:
        from .evaluate import run_evaluation

        run_evaluation(
            dataset_cfg_path=args.data_config,
            eval_cfg_path=args.eval_config,
            target_ds_ids=args.prep_ids,
            model_names=args.models,
            tasks=args.eval_tasks,
        )

    logging.info("ST-CNVBench completed successfully.")
    return 0
