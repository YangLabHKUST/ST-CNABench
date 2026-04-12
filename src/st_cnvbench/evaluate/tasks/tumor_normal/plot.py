import json
import logging
from functools import reduce
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import squidpy as sq
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from PIL import Image, ImageEnhance

from ...utils.plot_style import apply_plot_style, resolve_method_color


def _to_pdf_path(path_like):
    return Path(path_like).with_suffix(".pdf")


def _auto_ylim(values):
    finite_vals = np.asarray(values, dtype=float)
    finite_vals = finite_vals[np.isfinite(finite_vals)]
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


def run_bar_plot(metrics_df, plot_path):
    """
    Bar plot for Accuracy, MCC, Spatial Coherence Z with Call Rate in legend.
    """
    if metrics_df is None or metrics_df.empty:
        logging.warning("Metrics summary is empty, skipping bar plot.")
        return

    apply_plot_style(dpi=300)
    sns.set_theme(style="ticks")

    metrics_to_plot = [
        "Accuracy",
        "Precision",
        "Recall",
        "F1_Score",
        "MCC",
        "ARI",
        "Spatial_Coherence_Z",
    ]
    metrics_to_plot = [m for m in metrics_to_plot if m in metrics_df.columns]
    if not metrics_to_plot:
        logging.warning("No expected metrics found for tumor_normal bar plot.")
        return

    models = metrics_df["Model"].unique()
    n_models = max(1, len(models))
    fig_w = min(28.0, max(16.0, 1.6 * n_models + 8.0))

    fig, axes = plt.subplots(2, 4, figsize=(fig_w, 10), dpi=300)
    fig.suptitle("Tumor Detection Quantitative Evaluation", fontsize=18, fontweight="bold", y=1.02)

    axes_flat = axes.flatten()
    model_colors = {m: resolve_method_color(m) for m in models}

    legend_labels = {
        row["Model"]: f"{row['Model']} (Call: {row['Call_Rate']:.1%})"
        for _, row in metrics_df.iterrows()
    }

    unit_metrics = {"Accuracy", "Precision", "Recall", "F1_Score"}
    signed_unit_metrics = {"MCC", "ARI"}

    for i, metric in enumerate(metrics_to_plot):
        ax = axes_flat[i]
        sns.barplot(
            data=metrics_df,
            x="Model",
            y=metric,
            ax=ax,
            palette=model_colors,
            hue="Model",
            legend=False,
            edgecolor="black",
            linewidth=0.8,
            errorbar=None,
        )

        ax.set_title(metric.replace("_", " "), fontsize=13, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(metric.replace("_", " "), fontsize=11)

        if metric in unit_metrics:
            ax.set_ylim(0, 1.05)
        elif metric in signed_unit_metrics:
            ax.set_ylim(-1.05, 1.05)
            ax.axhline(0, color="#555555", linewidth=1.0, linestyle="--", alpha=0.7)
        else:
            limits = _auto_ylim(metrics_df[metric].to_numpy())
            if limits is not None:
                ax.set_ylim(*limits)
                if limits[0] < 0 < limits[1]:
                    ax.axhline(0, color="#555555", linewidth=1.0, linestyle="--", alpha=0.7)

        for p in ax.patches:
            height = p.get_height()
            if np.isnan(height) or np.isclose(height, 0.0):
                continue
            offset = 3 if height >= 0 else -3
            va = "bottom" if height >= 0 else "top"
            ax.annotate(
                f"{height:.2f}",
                (p.get_x() + p.get_width() / 2.0, height),
                ha="center",
                va=va,
                fontsize=8,
                fontweight="bold",
                xytext=(0, offset),
                textcoords="offset points",
            )

        rotation = 45 if n_models > 8 else 30
        for label in ax.get_xticklabels():
            label.set_rotation(rotation)
            label.set_ha("right")
            label.set_fontsize(9)
        ax.tick_params(axis="y", labelsize=9)

    for j in range(len(metrics_to_plot), len(axes_flat)):
        axes_flat[j].axis("off")

    legend_elements = [Patch(facecolor=model_colors[m], label=legend_labels[m]) for m in models]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.07),
        ncol=min(len(models), 4),
        frameon=False,
        title="Models & Data Call Rate",
        fontsize=9,
        title_fontsize=10,
    )

    fig.tight_layout(rect=[0, 0.09, 1, 0.96])
    save_path = Path(plot_path)
    if not save_path.suffix:
        save_path = save_path.with_suffix(".png")
    fig.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close(fig)


def run_confusion_map(results, plot_path):
    apply_plot_style(dpi=300)
    sns.set_theme(style="ticks")

    if "GT" not in results:
        logging.error("GT data missing in results! Cannot generate confusion map.")
        return

    gt_df = results["GT"].copy()

    models_to_plot = [m for m in results.keys() if m != "GT"]
    n_models = len(models_to_plot)

    if n_models == 0:
        logging.warning("No models found in results for comparison.")
        return

    n_cols = min(3, n_models)
    n_rows = (n_models + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows), dpi=300)

    if n_models == 1:
        axes_flat = [axes]
    else:
        axes_flat = np.array(axes).flatten()

    for i, model_name in enumerate(models_to_plot):
        ax = axes_flat[i]
        model_df = results[model_name]

        blues_custom = ListedColormap(
            [
                "#f7fbff",
                "#deebf7",
                "#c6dbef",
                "#9ecae1",
                "#6baed6",
                "#4292c6",
                "#2171b5",
                "#08519c",
                "#08306b",
            ]
        )

        model_df = model_df[model_df.iloc[:, 1] != -1]
        gt_df_filtered = gt_df[gt_df["Barcodes"].isin(model_df["Barcodes"])]
        merged = pd.merge(
            gt_df_filtered[["Barcodes", "Preds"]],
            model_df.iloc[:, [0, 1]],
            on="Barcodes",
            how="inner",
            suffixes=("_true", "_pred"),
        )

        y_true = merged.iloc[:, 1]
        y_pred = merged.iloc[:, 2]
        cm = pd.crosstab(
            y_true,
            y_pred,
            rownames=["Actual"],
            colnames=["Predicted"],
            dropna=False,
        )

        cm = cm.reindex(index=[0, 1], columns=[0, 1], fill_value=0)

        cm.index = ["Normal", "Tumor"]
        cm.columns = ["Normal", "Tumor"]

        sns.heatmap(
            cm,
            annot=True,
            fmt="d",
            cmap=blues_custom,
            cbar=False,
            ax=ax,
            annot_kws={"size": 14, "fontweight": "bold"},
        )

        ax.set_title(f"Model: {model_name}", fontsize=14, fontweight="bold", pad=8)
        ax.set_xlabel("Predicted Label", fontsize=11)
        ax.set_ylabel("Actual Label", fontsize=11)

    for j in range(i + 1, n_rows * n_cols):
        axes_flat[j].axis("off")

    fig.suptitle("Confusion Matrix Analysis", fontsize=18, fontweight="bold", y=1.02)
    fig.tight_layout()

    save_path = Path(plot_path)
    if not save_path.suffix:
        save_path = save_path.with_suffix(".png")
    fig.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    logging.info(f"Confusion matrix heatmap saved to: {save_path}")


def run_spatial_plot(results, ds_spatial_coord, ds_HE_img, ds_scalefactors, plot_path):
    """
    Spatial Visualization of tumor/normal predictions from multiple models.
    Plot with HE image background.
    """
    apply_plot_style(dpi=300)

    coords_df = pd.read_csv(ds_spatial_coord, index_col=0)
    with open(ds_scalefactors, "r") as f:
        scalefactors = json.load(f)

    img_pil = Image.open(ds_HE_img)
    img_pil = ImageEnhance.Brightness(img_pil).enhance(1.5)
    image = np.array(img_pil)

    logging.info("Aligning spatial plot indices across all sources...")
    index_list = [coords_df.index]
    for _, df in results.items():
        index_list.append(df.set_index(df.columns[0]).index)
    common_index = reduce(lambda x, y: x.intersection(y), index_list)
    logging.info(f"Common index size for plotting: {len(common_index)}")

    palette_dict = {
        "1": "#FA9D75",
        "0": "#82B0D2",
        "-1": "#BDBDBD",
    }

    adata = ad.AnnData(obs=pd.DataFrame(index=common_index))
    for model_name, res_df in results.items():
        res_series = res_df.set_index(res_df.columns[0]).iloc[:, 0]

        adata.obs[model_name] = (
            res_series.loc[common_index]
            .fillna(-1)
            .astype(float)
            .astype(int)
            .astype(str)
            .astype("category")
        )

        adata.obs[model_name] = adata.obs[model_name].cat.set_categories(["-1", "0", "1"])
        adata.uns[f"{model_name}_colors"] = [palette_dict[c] for c in adata.obs[model_name].cat.categories]

    adata.obsm["spatial"] = coords_df.loc[common_index, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].values

    library_id = "sample_id"
    adata.uns["spatial"] = {
        library_id: {
            "images": {"hires": image.astype(float) / 255.0 if image.max() > 1 else image},
            "scalefactors": scalefactors,
            "metadata": {
                "coords_columns": list(coords_df.columns),
                "image_path": ds_HE_img,
                "image_shape": image.shape[:2],
            },
        }
    }

    model_list = ["GT"] + [m for m in results.keys() if m != "GT"]
    n_plots = len(model_list)
    n_cols = min(4, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    panel_w = 5.2
    panel_h = 5.2
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(panel_w * n_cols, panel_h * n_rows),
        dpi=300,
    )
    if n_plots == 1:
        axes_flat = [axes]
    else:
        axes_flat = np.array(axes).flatten()

    for i, model in enumerate(model_list):
        sq.pl.spatial_scatter(
            adata,
            color=model,
            library_id=library_id,
            img_res_key="hires",
            ax=axes_flat[i],
            title=f"{model} Detection",
            size=1.5,
            alpha=1,
        )
        leg = axes_flat[i].get_legend()
        if leg is not None:
            leg.remove()

    for j in range(i + 1, n_rows * n_cols):
        axes_flat[j].axis("off")

    legend_elements = [
        Patch(facecolor="#FA9D75", edgecolor="none", label="Tumor"),
        Patch(facecolor="#82B0D2", edgecolor="none", label="Normal"),
        Patch(facecolor="#BDBDBD", edgecolor="gray", label="Filtered"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        title="Detection Results",
        title_fontsize=12,
        fontsize=10,
        frameon=False,
    )

    fig.tight_layout(rect=[0, 0, 0.9, 1])
    fig.subplots_adjust(hspace=0.34, wspace=0.10)
    save_path = Path(plot_path)
    if not save_path.suffix:
        save_path = save_path.with_suffix(".png")
    fig.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close(fig)

    per_model_dir = save_path.parent / "spatial_plots_by_model"
    per_model_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Saving per-model tumor/normal spatial plots to: {per_model_dir}")

    used_filenames = set()
    for model in model_list:
        fig_single, ax_single = plt.subplots(figsize=(panel_w, panel_h), dpi=300)
        sq.pl.spatial_scatter(
            adata,
            color=model,
            library_id=library_id,
            img_res_key="hires",
            ax=ax_single,
            title="",
            size=1.5,
            alpha=1,
        )

        single_leg = ax_single.get_legend()
        if single_leg is not None:
            single_leg.remove()
        for fig_leg in list(fig_single.legends):
            fig_leg.remove()

        ax_single.set_title("")
        ax_single.set_xlabel("")
        ax_single.set_ylabel("")
        ax_single.set_axis_off()

        safe_model_name = "".join(c if (c.isalnum() or c in "._-") else "_" for c in str(model)).strip("_")
        if not safe_model_name:
            raise ValueError(f"Invalid model name for output filename: {model}")
        if safe_model_name in used_filenames:
            raise ValueError(f"Duplicated output filename after sanitization: {safe_model_name}")
        used_filenames.add(safe_model_name)

        single_plot_path = per_model_dir / f"{safe_model_name}.png"
        fig_single.savefig(single_plot_path, bbox_inches="tight", dpi=300, pad_inches=0)
        plt.close(fig_single)
        logging.info(f"Saved per-model tumor/normal spatial plot: {single_plot_path}")
