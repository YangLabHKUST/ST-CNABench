import logging
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ...utils.constants import (
    MASTER_MODEL_ORDER,
)
from ...utils.io import save_tsv
from ...utils.plot_style import apply_plot_style, resolve_method_color

import warnings

warnings.filterwarnings("ignore")


def _to_pdf_path(path_like):
    return Path(path_like).with_suffix(".pdf")


def _infer_metric_ylim(metric, values):
    metric_lower = str(metric).lower()
    finite_vals = np.asarray(values, dtype=float)
    finite_vals = finite_vals[np.isfinite(finite_vals)]

    if metric in {"PCC", "Spearman_CC"} or "pcc" in metric_lower or "spearman" in metric_lower:
        if finite_vals.size == 0:
            return (-1.05, 1.05)
        vmin = float(np.min(finite_vals))
        vmax = float(np.max(finite_vals))
        if vmin >= 0:
            return (0.0, 1.05)
        if vmax <= 0:
            return (-1.05, 0.05)
        return (-1.05, 1.05)

    if any(
        key in metric_lower
        for key in ["auc", "f1", "precision", "sensitivity", "recall", "accuracy", "specificity"]
    ):
        return (0.0, 1.05)

    if finite_vals.size == 0:
        return None

    vmin = float(np.min(finite_vals))
    vmax = float(np.max(finite_vals))
    if np.isclose(vmin, vmax):
        pad = 0.1 * abs(vmax) if not np.isclose(vmax, 0.0) else 0.1
    else:
        pad = 0.1 * (vmax - vmin)

    lower = vmin - pad
    upper = vmax + pad
    if np.isclose(lower, upper):
        upper = lower + 1.0
    return (lower, upper)


def run_bar_plot(df_metrics, output_path):
    """
    Render one bar plot per metric under the output directory.
    """
    if df_metrics is None or df_metrics.empty:
        logging.warning("Metrics DataFrame is empty, skipping cnv_profile bar plots.")
        return

    out_path = Path(output_path)
    apply_plot_style(dpi=300)
    sns.set_style("ticks")

    metrics = [m for m in df_metrics["Metric"].dropna().unique()]

    for metric in metrics:
        sub_df = df_metrics[df_metrics["Metric"] == metric].copy()
        sub_df = sub_df.dropna(subset=["Value"])

        if sub_df.empty:
            logging.warning(f"No valid data for metric: {metric}, skipping plot.")
            continue

        local_order = [m for m in MASTER_MODEL_ORDER if m in sub_df["Method"].unique()]
        local_order += [m for m in sub_df["Method"].unique() if m not in local_order]

        current_palette = [resolve_method_color(m) for m in local_order]

        n_models = max(1, len(local_order))
        fig_w = min(12.5, max(5.2, 0.36 * n_models + 2.8))
        fig, ax = plt.subplots(figsize=(fig_w, 5.6), dpi=300)

        ordered_df = (
            sub_df[["Method", "Value"]]
            .drop_duplicates(subset=["Method"], keep="last")
            .set_index("Method")
            .reindex(local_order)
        )
        bar_values = ordered_df["Value"].to_numpy(dtype=float)

        # Compact layout: tighter method spacing + slightly narrower bars.
        center_step = 0.44
        bar_width = 0.34
        x_pos = np.arange(n_models, dtype=float) * center_step

        ax.bar(
            x_pos,
            bar_values,
            width=bar_width,
            color=current_palette,
            edgecolor="black",
            linewidth=1.0,
        )

        ax.set_xticks(x_pos)
        ax.set_xticklabels(local_order)
        if n_models == 1:
            ax.set_xlim(-0.35, 0.35)
        else:
            margin = center_step * 0.55
            ax.set_xlim(x_pos[0] - margin, x_pos[-1] + margin)
        ax.margins(x=0.01)

        y_limits = _infer_metric_ylim(metric, sub_df["Value"].to_numpy())
        if y_limits is not None:
            ax.set_ylim(*y_limits)
            if np.isclose(y_limits[0], 0.0):
                ax.margins(y=0)
            if y_limits[0] < 0 < y_limits[1]:
                ax.axhline(0, color="#555555", linewidth=1.0, linestyle="--", alpha=0.7)

        for p in ax.patches:
            h = p.get_height()
            if np.isnan(h) or np.isclose(h, 0.0):
                continue
            offset = 4 if h >= 0 else -4
            va = "bottom" if h >= 0 else "top"
            ax.annotate(
                f"{h:.3f}",
                (p.get_x() + p.get_width() / 2.0, h),
                ha="center",
                va=va,
                fontsize=9,
                fontweight="bold",
                xytext=(0, offset),
                textcoords="offset points",
            )

        ax.set_title(f"Performance: {metric}", fontsize=15, fontweight="bold", pad=12)
        ax.set_xlabel("")
        ax.set_ylabel(metric.replace("_", " "), fontsize=12, fontweight="bold")

        rotation = 45 if n_models > 8 else 30
        for label in ax.get_xticklabels():
            label.set_rotation(rotation)
            label.set_ha("right")
            label.set_fontsize(10)

        ax.tick_params(axis="y", labelsize=10)
        sns.despine(ax=ax)
        fig.tight_layout()

        metric_safe = str(metric).replace("/", "_")
        metric_out_path = out_path.parent / f"{out_path.stem}_{metric_safe}.pdf"
        fig.savefig(metric_out_path, bbox_inches="tight", format="pdf")
        plt.close(fig)

    logging.info(f"Generated {len(metrics)} separate bar plots in {out_path.parent}")


def run_karyogram_level_plot(aligned_level_df, level_name, output_path):
    """
    Karyogram plot for a specific CNV profile level.
    Adapted for Dual-Track (CN_Score + LOH_Status) with rigorous Z-score standardization.
    """
    if aligned_level_df is None or aligned_level_df.empty:
        logging.warning(f"No valid aligned data for level: {level_name}. Skipping karyogram plot.")
        return

    apply_plot_style(dpi=300)
    sns.set_style("white")

    df = aligned_level_df.copy()

    if "BinID" in df.columns:
        df = df.set_index("BinID")
    elif "ID" in df.columns:
        df = df.set_index("ID")

    meta_cols = ["Chromosome", "Start", "End"]
    meta_cols_present = [c for c in meta_cols if c in df.columns]
    plot_df = df.drop(columns=[c for c in meta_cols if c in df.columns], errors="ignore")

    def _get_sort_key(bin_id):
        parts = str(bin_id).split("_")
        ch = parts[0].replace("chr", "").replace("Chr", "")
        ch_map = {"X": 23, "Y": 24, "x": 23, "y": 24}
        c_val = ch_map.get(ch, int(ch) if ch.isdigit() else 99)
        pos_val = int(parts[1]) if len(parts) > 1 and str(parts[1]).isdigit() else 0
        return (c_val, pos_val)

    plot_df["sort_key"] = plot_df.index.map(_get_sort_key)
    plot_df = plot_df.sort_values("sort_key").drop(columns="sort_key")
    meta_df = df.loc[plot_df.index, meta_cols_present].copy() if meta_cols_present else pd.DataFrame(index=plot_df.index)

    cn_cols = [c for c in plot_df.columns if c.endswith("_CN_Score") or c == "GT_Score"]
    loh_cols = [c for c in plot_df.columns if c.endswith("_LOH_Status") or c == "GT_LOH_Status"]

    if not cn_cols and not loh_cols:
        cn_cols = list(plot_df.columns)

    out_path = Path(output_path)
    loh_cmap = mcolors.LinearSegmentedColormap.from_list(
        "loh_status",
        ["#ffffff", "#fcae91", "#fb6a4a", "#cb181d"],
    )

    def _format_model_label(model_name):
        if model_name == "GT":
            return "GT"

        lower = str(model_name).lower()
        if lower.endswith("_expr"):
            return f"{model_name[:-5]}(expr)"
        if lower.endswith("_cnv"):
            return f"{model_name[:-4]}(cnv)"
        if lower.endswith("_nowgs"):
            return f"{model_name[:-7]}(NoWGS)"
        if lower.endswith("_wgs"):
            return f"{model_name[:-4]}(WGS)"
        return str(model_name)

    def _draw_heatmap(
        sub_df,
        suffix,
        cmap,
        vmin,
        vmax,
        center,
        cbar_label,
        cbar_ticks,
        cbar_tick_labels,
        export_sub_df=None,
    ):
        if sub_df.empty or sub_df.columns.empty:
            return

        if export_sub_df is None:
            export_sub_df = sub_df.copy()

        rename_map = {}
        for c in sub_df.columns:
            clean_name = c.replace("_CN_Score", "").replace("_LOH_Status", "")
            clean_name = clean_name.replace("GT_Score", "GT").replace("GT_LOH_Status", "GT")
            rename_map[c] = clean_name
        sub_df = sub_df.rename(columns=rename_map)

        plot_data = sub_df.T
        present_models = [m for m in MASTER_MODEL_ORDER if m in plot_data.index]
        final_order = (["GT"] if "GT" in plot_data.index else []) + [m for m in present_models if m != "GT"]
        final_order += [m for m in plot_data.index if m not in final_order]
        plot_data = plot_data.reindex(final_order)

        n_models = max(1, plot_data.shape[0])
        fig_w = 10.0
        fig_h = min(6.0, max(2.4, 0.25 * n_models + 1.4))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)

        hm = sns.heatmap(
            plot_data,
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            xticklabels=False,
            yticklabels=True,
            ax=ax,
            cbar_kws={
                "label": cbar_label,
                "orientation": "horizontal",
                "pad": 0.24,
                "shrink": 0.55,
                "ticks": cbar_ticks,
            },
        )

        if cbar_tick_labels is not None:
            cbar = hm.collections[0].colorbar
            cbar.set_ticklabels(cbar_tick_labels)
            cbar.ax.tick_params(labelsize=8)
            cbar.set_label(cbar_label, fontsize=9)

        chroms = [str(bid).split("_")[0] for bid in sub_df.index]
        boundaries = [0]
        labels = []
        curr_chr = None
        for i, c in enumerate(chroms):
            if c != curr_chr:
                if curr_chr is not None:
                    boundaries.append(i)
                labels.append(c.lower().replace("chr", ""))
                curr_chr = c
        boundaries.append(len(sub_df))

        for b in boundaries:
            ax.axvline(x=b, color="#4a4a4a", linewidth=0.70, alpha=0.40)

        # Subtle outer frame improves panel definition in multi-panel figures.
        ax.add_patch(
            plt.Rectangle(
                (0, 0),
                plot_data.shape[1],
                plot_data.shape[0],
                fill=False,
                edgecolor="#5a5a5a",
                linewidth=0.55,
                alpha=0.45,
                clip_on=False,
            )
        )

        label_pos = [(boundaries[i] + boundaries[i + 1]) / 2 for i in range(len(boundaries) - 1)]
        ax.set_xticks(label_pos)
        ax.set_xticklabels(labels, fontsize=8.5, fontweight="normal")
        ax.tick_params(axis="x", which="both", length=0)

        display_model_labels = [_format_model_label(m) for m in plot_data.index]
        ax.set_yticklabels(display_model_labels, rotation=0, fontweight="normal", fontsize=9)
        ax.tick_params(axis="y", length=0)

        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_label_position("bottom")
        ax.yaxis.set_ticks_position("left")

        ax.set_ylabel("Models", fontsize=11)
        ax.set_xlabel("Chromosome", fontsize=11)

        fig.tight_layout()

        if suffix:
            save_path = out_path.parent / f"{out_path.stem}_{suffix}.pdf"
        else:
            save_path = _to_pdf_path(out_path)

        fig.savefig(save_path, bbox_inches="tight", format="pdf")
        plt.close(fig)
        logging.info(f"Karyogram plot saved to {save_path}")

        # Export the exact aligned profile used in this heatmap for downstream inspection.
        export_data = export_sub_df.rename(columns=rename_map).reindex(columns=final_order)
        export_df = pd.concat([meta_df.reindex(export_data.index), export_data], axis=1)
        export_df.insert(0, "BinID", export_df.index)
        export_df = export_df.reset_index(drop=True)

        if suffix:
            export_path = out_path.parent / f"{out_path.stem}_{suffix}_aligned_profile.tsv"
        else:
            export_path = out_path.parent / f"{out_path.stem}_aligned_profile.tsv"
        save_tsv(export_df, export_path)
        logging.info(f"Karyogram aligned profile saved to {export_path}")

    if cn_cols:
        raw_cn_df = plot_df[cn_cols].copy()
        df_cn = raw_cn_df.copy()
        for col in df_cn.columns:
            std_val = df_cn[col].std()
            mean_val = df_cn[col].mean()
            if std_val == 0 or pd.isna(std_val):
                std_val = 1
            df_cn[col] = (df_cn[col] - mean_val) / std_val

        _draw_heatmap(
            sub_df=df_cn,
            suffix="CN_Score",
            cmap="RdBu_r",
            vmin=-4,
            vmax=4,
            center=0,
            cbar_label="Z-scored Relative CN Signal (Loss/Neutral/Gain)",
            cbar_ticks=[-4, 0, 4],
            cbar_tick_labels=["Loss", "Neutral", "Gain"],
            export_sub_df=raw_cn_df,
        )

    if loh_cols:
        df_loh = plot_df[loh_cols].copy()

        _draw_heatmap(
            sub_df=df_loh,
            suffix="LOH_Status",
            cmap=loh_cmap,
            vmin=0,
            vmax=1,
            center=None,
            cbar_label="LOH Population Frequency",
            cbar_ticks=[0, 0.5, 1.0],
            cbar_tick_labels=["0 (Hetero)", "0.5", "1 (LOH)"],
        )
