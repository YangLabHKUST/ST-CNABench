import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ...utils.constants import (
    MODEL_ORDER,
)
from ...utils.plot_style import apply_plot_style, resolve_method_color


def run_plot(df_efficiency, output_path):
    logging.info("Generating computational efficiency plots...")

    if df_efficiency is None or df_efficiency.empty:
        logging.warning("Efficiency DataFrame is empty, skipping plot.")
        return

    df_plot = df_efficiency.copy()
    df_plot["runtime_sec"] = pd.to_numeric(df_plot["runtime_sec"], errors="coerce")
    df_plot["mem_gb"] = pd.to_numeric(df_plot["mem_gb"], errors="coerce")
    df_plot = df_plot.dropna(subset=["model", "runtime_sec", "mem_gb"])

    if df_plot.empty:
        logging.warning("No valid efficiency values after numeric coercion. Skipping plot.")
        return

    existing_models = df_plot["model"].unique().tolist()
    plot_order = [m for m in MODEL_ORDER if m in existing_models]
    plot_order += [m for m in existing_models if m not in plot_order]

    df_plot["model"] = pd.Categorical(df_plot["model"], categories=plot_order, ordered=True)
    df_plot = df_plot.sort_values("model")

    apply_plot_style(dpi=300)
    sns.set_theme(style="ticks")

    n_models = max(1, len(plot_order))
    fig_w = min(20.0, max(12.0, 1.2 * n_models + 4.0))
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, 5.6), dpi=300)

    palette_map = {m: resolve_method_color(m) for m in plot_order}

    metrics = [
        ("runtime_sec", "Runtime (s)", "Total Runtime"),
        ("mem_gb", "Memory (GB)", "Peak Memory"),
    ]

    for idx, (col, ylabel, title) in enumerate(metrics):
        ax = axes[idx]
        sns.barplot(
            data=df_plot,
            x="model",
            y=col,
            order=plot_order,
            ax=ax,
            palette=palette_map,
            hue="model",
            legend=False,
            edgecolor="black",
            linewidth=1.0,
            errorbar=None,
        )

        max_val = float(df_plot[col].max()) if not df_plot[col].isna().all() else 0.0
        upper = max_val * 1.15 if max_val > 0 else 1.0
        ax.set_ylim(0, upper)
        ax.margins(y=0)

        ax.set_title(title, fontsize=14, fontweight="bold", pad=10)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlabel("")

        rotation = 45 if n_models > 8 else 30
        for label in ax.get_xticklabels():
            label.set_rotation(rotation)
            label.set_ha("right")
            label.set_fontsize(9)
        ax.tick_params(axis="y", labelsize=10)

        for p in ax.patches:
            h = p.get_height()
            if pd.isna(h) or np.isclose(h, 0.0):
                continue
            label = f"{h:.1f}" if h < 100 else f"{h:.0f}"
            ax.annotate(
                label,
                (p.get_x() + p.get_width() / 2.0, h),
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
                xytext=(0, 3),
                textcoords="offset points",
            )

        sns.despine(ax=ax)

    fig.tight_layout()

    save_path = Path(output_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
