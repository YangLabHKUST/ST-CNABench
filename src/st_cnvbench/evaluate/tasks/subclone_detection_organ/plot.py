import logging
import os
import re
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage

from ...utils.plot_style import apply_plot_style, resolve_method_color

_CATEGORICAL_PALETTE = "tab20"

def run_plot_metric_heatmap(matrix, metric_name, save_dir):
    """Render a generic metric heatmap."""
    if matrix is None or matrix.empty:
        return

    apply_plot_style(dpi=300)
    sns.set_style("ticks")

    fig, ax = plt.subplots(figsize=(8.2, 6.0), dpi=300)
    sns.heatmap(
        matrix.astype(float),
        annot=True,
        fmt=".2f",
        cmap="RdYlBu_r",
        vmin=0,
        vmax=1,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"shrink": 0.85, "label": metric_name},
        ax=ax,
    )
    ax.set_title(f"{metric_name} Correlation Heatmap", fontsize=15, fontweight="bold", pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=10, rotation=30)
    ax.tick_params(axis="y", labelsize=10, rotation=0)
    fig.tight_layout()

    file_name = f"{metric_name.lower().replace('-', '_')}_heatmap.png"
    save_path = save_dir / file_name
    fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    logging.info(f"Saved {metric_name} heatmap to {save_path}")
    plt.close(fig)


def run_plot_metric_barplot(scores_dict, metric_name, save_dir):
    """Render a generic score barplot against ground truth."""
    if not scores_dict:
        return

    apply_plot_style(dpi=300)
    sns.set_style("ticks")

    scores_series = pd.Series(scores_dict).sort_values(ascending=False)
    methods = scores_series.index.tolist()
    colors = [resolve_method_color(m) for m in methods]

    fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=300)
    bars = ax.bar(
        methods,
        scores_series.values,
        color=colors,
        edgecolor="black",
        linewidth=0.9,
        width=0.66,
    )

    for bar in bars:
        height = bar.get_height()
        if pd.notna(height):
            ax.annotate(
                f"{height:.2f}",
                (bar.get_x() + bar.get_width() / 2.0, height),
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
                xytext=(0, 3),
                textcoords="offset points",
            )

    ax.set_ylabel(f"{metric_name} Score", fontsize=12, fontweight="bold")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Model Performance vs Ground Truth ({metric_name})", fontsize=15, fontweight="bold", pad=12)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelsize=10, rotation=35)
    ax.tick_params(axis="y", labelsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()

    file_name = f"{metric_name.lower().replace('-', '_')}_barplot.png"
    save_path = save_dir / file_name
    fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    logging.info(f"Saved {metric_name} barplot to {save_path}")
    plt.close(fig)


def run_spatial_plot_no_img(results, save_path):
    """Render organ-wide spatial scatter plots from the raw label table."""
    save_path = Path(save_path)
    panel_dir = save_path.parent / f"{save_path.stem}_panels"
    panel_dir.mkdir(parents=True, exist_ok=True)

    if "GT_Preds" not in results.columns:
        logging.error("GT_Preds not found in results, skip organ spatial plotting.")
        return

    exclude_cols = ["Barcodes", "pseudo_barcode", "original_barcode", "group", "x", "y"]
    model_list = ["GT_Preds"] + [c for c in results.columns if c not in exclude_cols and c != "GT_Preds"]
    if not model_list:
        logging.warning("No model columns found for spatial plotting.")
        return

    n_plots = len(model_list)
    n_cols = min(2, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    apply_plot_style(dpi=300)
    sns.set_style("white")

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11.5 * n_cols, 4.8 * n_rows), dpi=300)

    if n_plots == 1:
        axes_flat = [axes]
    else:
        axes_flat = axes.flatten()

    for i, model in enumerate(model_list):
        ax = axes_flat[i]
        plot_data, class_order, color_map = _prepare_spatial_plot_data(results, model)
        if plot_data is None or plot_data.empty:
            ax.axis("off")
            continue

        _render_spatial_subplot(ax, plot_data, model, class_order, color_map)

        # Fixed two-column layout keeps the spatial panel size identical
        # across methods, independent of legend text length.
        fig_single = plt.figure(figsize=(7.2, 5.0), dpi=300)
        gs_single = fig_single.add_gridspec(1, 2, width_ratios=[1.0, 0.50], wspace=0.02)
        ax_single = fig_single.add_subplot(gs_single[0, 0])
        ax_single_legend = fig_single.add_subplot(gs_single[0, 1])
        _render_spatial_subplot(ax_single, plot_data, model, class_order, color_map, legend_ax=ax_single_legend)
        fig_single.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.07)
        safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", _strip_preds_suffix(model))
        panel_path = panel_dir / f"{i + 1:02d}_{safe_name}.png"
        fig_single.savefig(panel_path, dpi=300, facecolor="white")
        plt.close(fig_single)

    for j in range(n_plots, len(axes_flat)):
        axes_flat[j].axis("off")

    fig.tight_layout()
    fig.savefig(save_path, bbox_inches="tight", dpi=300, facecolor="white")
    logging.info(f"Saved organ-wide spatial scatter plot to {save_path}")
    logging.info(f"Saved per-panel spatial plots to {panel_dir}")
    plt.close(fig)


def _run_plot_karyotype_composition_single(df_all, clone_cnv_profiles, save_dir, output_filename):
    apply_plot_style(dpi=300)

    if "GT_Preds" not in df_all.columns:
        logging.error("GT_Preds not found in dataframe, skip unified karyogram plot.")
        return
    if not clone_cnv_profiles:
        logging.warning(f"No CNV profiles found, skipping {output_filename}.")
        return

    os.makedirs(save_dir, exist_ok=True)

    all_gt_classes = [str(x) for x in df_all["GT_Preds"].dropna().unique()]
    if not all_gt_classes:
        logging.error("GT_Preds has no valid labels, skip unified karyogram plot.")
        return

    all_cnv_list = []
    for mod, df in clone_cnv_profiles.items():
        if df is None or df.empty:
            continue
        temp = df.copy()
        if "Chr" in temp.columns and "Chromosome" not in temp.columns:
            temp = temp.rename(columns={"Chr": "Chromosome"})
        score_col = "CN_Score" if "CN_Score" in temp.columns else "TCN_Score"
        if score_col in temp.columns:
            temp["Unified_Score"] = temp[score_col]
            temp["Clone_ID"] = temp["Clone_ID"].map(lambda x: _format_clone_name(mod, x))
            temp = temp.dropna(subset=["Clone_ID"]).copy()
            req_cols = ["Clone_ID", "Chromosome", "Start", "End", "Unified_Score"]
            if all(col in temp.columns for col in req_cols):
                all_cnv_list.append(temp[req_cols])
    if not all_cnv_list:
        logging.warning(f"No valid CNV profiles found for {output_filename}.")
        return

    unified_cnv_df = pd.concat(all_cnv_list, ignore_index=True)
    unified_cnv_df["Bin_ID"] = (
        unified_cnv_df["Chromosome"].astype(str)
        + "_"
        + unified_cnv_df["Start"].astype(str)
        + "_"
        + unified_cnv_df["End"].astype(str)
    )

    def _get_sort_key(chrom):
        chrom = str(chrom).replace("chr", "").replace("Chr", "")
        chrom_map = {"X": 23, "Y": 24, "x": 23, "y": 24}
        return chrom_map.get(chrom, int(chrom) if chrom.isdigit() else 99)

    unified_cnv_df["chr_sort"] = unified_cnv_df["Chromosome"].map(_get_sort_key)
    unified_cnv_df["Start"] = unified_cnv_df["Start"].astype(float)
    unified_cnv_df = unified_cnv_df.sort_values(["chr_sort", "Start"])

    cnv_matrix = unified_cnv_df.pivot(index="Clone_ID", columns="Bin_ID", values="Unified_Score")
    ordered_bins = unified_cnv_df["Bin_ID"].unique()
    cnv_matrix = cnv_matrix[ordered_bins]

    # Export pre-Zscore clone-by-bin matrix for downstream analyses.
    output_stem = Path(output_filename).stem
    pre_zscore_path = Path(save_dir) / f"{output_stem}_CNV_Matrix_PreZScore.tsv"
    cnv_matrix.to_csv(pre_zscore_path, sep="\t")
    logging.info(f"Saved pre-Zscore CNV matrix to {pre_zscore_path}")

    for clone in cnv_matrix.index:
        std_val = cnv_matrix.loc[clone].std()
        mean_val = cnv_matrix.loc[clone].mean()
        if std_val == 0 or pd.isna(std_val):
            std_val = 1
        cnv_matrix.loc[clone] = (cnv_matrix.loc[clone] - mean_val) / std_val
    cnv_matrix = cnv_matrix.fillna(0)

    dendro_linkage = None
    if len(cnv_matrix) > 1:
        dendro_linkage = linkage(cnv_matrix.values, method="ward")
        idx = leaves_list(dendro_linkage)
        ordered_clones = cnv_matrix.index[idx].tolist()
    else:
        ordered_clones = cnv_matrix.index.tolist()
    cnv_matrix = cnv_matrix.loc[ordered_clones]

    all_cross_tabs = []
    all_gt_labels = [str(x) for x in df_all["GT_Preds"].dropna().unique()]
    for mod in clone_cnv_profiles.keys():
        if mod == "GT":
            gt_clones = clone_cnv_profiles["GT"]["Clone_ID"].unique()
            cross_data = {}
            for clone in gt_clones:
                row_name = _format_clone_name(mod, clone)
                if row_name is None:
                    continue
                cross_data[row_name] = _build_gt_identity_row(clone, all_gt_labels)
            if cross_data:
                cross = pd.DataFrame.from_dict(cross_data, orient="index")
                all_cross_tabs.append(cross)
            continue

        mod_col = mod if mod in df_all.columns else f"{mod}_Preds"
        if mod_col not in df_all.columns:
            logging.warning(f"Model column {mod_col} not found in dataframe, skipping purity for {mod}.")
            continue
        valid_df = df_all.dropna(subset=["GT_Preds", mod_col]).copy()
        if valid_df.empty:
            continue
        valid_df[mod_col] = valid_df[mod_col].astype(str)
        cross = pd.crosstab(valid_df[mod_col], valid_df["GT_Preds"])
        cross.index = [_format_clone_name(mod, clone) for clone in cross.index]
        cross = cross[cross.index.notna()]
        all_cross_tabs.append(cross)

    if not all_cross_tabs:
        logging.error(f"Failed to compute GT composition table for {output_filename}.")
        return

    master_cross_tab = pd.concat(all_cross_tabs).fillna(0)
    master_cross_tab = master_cross_tab.groupby(master_cross_tab.index).sum()
    cross_tab_norm = master_cross_tab.div(master_cross_tab.sum(axis=1), axis=0)

    cols_order = sorted(cross_tab_norm.columns.tolist(), key=_gt_class_sort_key)
    gt_color_map = _build_categorical_color_map(cols_order)
    gt_colors = [gt_color_map[col] for col in cols_order]

    comp_data = cross_tab_norm.reindex(ordered_clones).fillna(0)
    comp_data = comp_data[cols_order]

    row_height = 0.22
    base_height = 3.5
    fig_height = max(5.0, base_height + len(ordered_clones) * row_height)
    fig = plt.figure(figsize=(16, fig_height), dpi=300)
    gs = gridspec.GridSpec(1, 5, width_ratios=[1.0, 1.4, 7.3, 0.1, 1.5], wspace=0.0)

    ax_dendro = fig.add_subplot(gs[0])
    if dendro_linkage is not None:
        dendrogram(dendro_linkage, orientation="left", ax=ax_dendro, color_threshold=0, above_threshold_color="black")
        ax_dendro.set_ylim(len(ordered_clones) * 10, 0)
    ax_dendro.axis("off")

    ax_cnv = fig.add_subplot(gs[2])
    bottom_margin = max(0.12, 1.8 / fig_height)
    cbar_ax = fig.add_axes([0.35, bottom_margin / 3, 0.3, 0.015])

    sns.heatmap(
        cnv_matrix,
        cmap="RdBu_r",
        center=0,
        vmin=-5,
        vmax=5,
        ax=ax_cnv,
        cbar_ax=cbar_ax,
        xticklabels=False,
        yticklabels=True,
        cbar_kws={"label": "Standardized Log2 Ratio", "orientation": "horizontal"},
    )
    cbar_ax.set_xticks([-5, 0, 5])
    cbar_ax.set_xticklabels(["Loss", "Neutral", "Gain"])

    chroms = [str(bin_id).split("_")[0] for bin_id in cnv_matrix.columns]
    boundaries = [0]
    labels = []
    curr_chr = None
    for i, chrom in enumerate(chroms):
        if chrom != curr_chr:
            if curr_chr is not None:
                boundaries.append(i)
            labels.append(chrom.lower().replace("chr", ""))
            curr_chr = chrom
    boundaries.append(len(cnv_matrix.columns))

    for boundary in boundaries:
        ax_cnv.axvline(x=boundary, color="#4a4a4a", linewidth=0.70, alpha=0.40)

    ax_cnv.add_patch(
        plt.Rectangle(
            (0, 0),
            cnv_matrix.shape[1],
            cnv_matrix.shape[0],
            fill=False,
            edgecolor="#5a5a5a",
            linewidth=0.55,
            alpha=0.45,
            clip_on=False,
        )
    )

    label_pos = [(boundaries[i] + boundaries[i + 1]) / 2 for i in range(len(boundaries) - 1)]
    ax_cnv.set_xticks(label_pos)
    ax_cnv.set_xticklabels(labels, fontsize=10, fontweight="normal", rotation=0)
    ax_cnv.tick_params(axis="x", which="both", length=0)

    y_fontsize = 12 if len(ordered_clones) <= 40 else 8
    ax_cnv.set_yticklabels(ax_cnv.get_yticklabels(), fontweight="normal", rotation=0, fontsize=y_fontsize)
    ax_cnv.tick_params(axis="y", length=0)
    ax_cnv.xaxis.set_ticks_position("bottom")
    ax_cnv.xaxis.set_label_position("bottom")
    ax_cnv.yaxis.set_ticks_position("left")
    ax_cnv.set_title("Unified Clustered Karyogram Profiles", fontsize=16, fontweight="bold", pad=15)
    ax_cnv.set_ylabel("")
    ax_cnv.set_xlabel("Chromosome", fontsize=11)

    ax_bar = fig.add_subplot(gs[4], sharey=ax_cnv)
    y_positions = np.arange(len(ordered_clones)) + 0.5
    left_bases = np.zeros(len(ordered_clones))
    for i, col in enumerate(cols_order):
        widths = comp_data[col].values
        ax_bar.barh(
            y=y_positions,
            width=widths,
            left=left_bases,
            height=1.0,
            color=gt_colors[i],
            edgecolor="white",
            linewidth=0.5,
            label=col,
        )
        left_bases += widths

    ax_bar.set_xlim(0, 1.0)
    ax_bar.set_title("Tumor Composition", fontsize=14, fontweight="bold", pad=15)
    ax_bar.set_facecolor("white")
    ax_bar.grid(False)
    for spine_loc in ["top", "right", "left"]:
        ax_bar.spines[spine_loc].set_visible(False)
    ax_bar.spines["bottom"].set_color("black")
    ax_bar.tick_params(axis="y", left=False, labelleft=False)

    handles, labels = ax_bar.get_legend_handles_labels()
    legend_cols = 2 if len(cols_order) > 3 else len(cols_order)
    fig.legend(
        handles,
        labels,
        title="Tumor Class",
        bbox_to_anchor=(0.88, bottom_margin / 5),
        loc="lower center",
        ncol=legend_cols,
        frameon=False,
        fontsize=10,
    )
    plt.subplots_adjust(bottom=bottom_margin, top=0.92, left=0.02, right=0.96)

    save_file = os.path.join(save_dir, output_filename)
    plt.savefig(save_file, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close()
    logging.info(f"Saved unified karyogram to {save_file}")


def run_plot_karyotype_composition(df_all, clone_cnv_profiles, save_dir):
    """Save the generic unified karyogram for all available methods."""
    _run_plot_karyotype_composition_single(
        df_all=df_all,
        clone_cnv_profiles=clone_cnv_profiles,
        save_dir=save_dir,
        output_filename="Unified_All_Methods_Clustered_Karyogram.png",
    )


def _strip_preds_suffix(model_name):
    name = str(model_name)
    return name[:-6] if name.endswith("_Preds") else name


def _parse_integer_label(label_value):
    text = str(label_value).strip()
    if not text or text.lower() == "nan":
        return None
    try:
        numeric = float(text)
    except ValueError:
        return None
    if not np.isfinite(numeric):
        return None
    integer = int(numeric)
    return integer if np.isclose(numeric, integer) else None


def _subclone_label_sort_key(label_text):
    text = str(label_text).strip()
    parsed = _parse_integer_label(text)
    if parsed is not None:
        return (0, parsed, "")
    suffix = _extract_suffix_integer(text)
    if suffix is not None:
        return (1, suffix, text.lower())
    return (2, 0, text.lower())


def _normalize_gt_class_name(class_name):
    return re.sub(r"[\s_\-]+", "", str(class_name).strip().lower())


def _extract_suffix_integer(text):
    match = re.search(r"(\d+)$", str(text).strip())
    return int(match.group(1)) if match else None


def _gt_class_sort_key(label_text):
    return _subclone_label_sort_key(label_text)


def _clean_label(label_value):
    text = str(label_value).strip()
    if not text or text.lower() == "nan":
        return None
    parsed = _parse_integer_label(text)
    if parsed is not None:
        return str(parsed)
    return text


def _simplify_model_label(model_name, label_value):
    text = _clean_label(label_value)
    if text is None:
        return None

    prefixes = [str(model_name).strip(), _strip_preds_suffix(str(model_name)).strip()]
    changed = True
    while changed:
        changed = False
        for prefix in prefixes:
            if not prefix:
                continue
            token = f"{prefix}_"
            if text.lower().startswith(token.lower()):
                text = text[len(token):].strip()
                changed = True

    parsed = _parse_integer_label(text)
    if parsed is not None:
        return str(parsed)
    return text


def _build_categorical_color_map(class_order):
    if not class_order:
        return {}
    palette = sns.color_palette(_CATEGORICAL_PALETTE, n_colors=len(class_order)).as_hex()
    return {label: palette[idx] for idx, label in enumerate(class_order)}


def _prepare_spatial_plot_data(results, model):
    if model not in results.columns:
        return None, [], {}

    plot_data = results.dropna(subset=[model, "x", "y"]).copy()
    if plot_data.empty:
        return plot_data, [], {}

    plot_data["_plot_label"] = plot_data[model].map(lambda x: _simplify_model_label(model, x))
    plot_data = plot_data.dropna(subset=["_plot_label"]).copy()
    if plot_data.empty:
        return plot_data, [], {}
    class_order = sorted(plot_data["_plot_label"].dropna().unique().tolist(), key=_subclone_label_sort_key)
    color_map = _build_categorical_color_map(class_order)
    return plot_data, class_order, color_map


def _format_clone_name(model_name, clone_label):
    clone_text = _simplify_model_label(model_name, clone_label)
    if clone_text is None:
        return None
    model_label = _strip_preds_suffix(str(model_name))
    return f"{model_label}_{clone_text}"


def _build_gt_identity_row(clone_label, gt_labels):
    clone_text = _clean_label(clone_label)
    if clone_text is None:
        return {label: 0.0 for label in gt_labels}

    clone_suffix = _extract_suffix_integer(clone_text)
    row_dict = {label: 0.0 for label in gt_labels}

    for label in gt_labels:
        label_text = _clean_label(label)
        if label_text == clone_text:
            row_dict[label] = 100.0
            return row_dict

        label_suffix = _extract_suffix_integer(label_text)
        if clone_suffix is not None and label_suffix is not None and clone_suffix == label_suffix:
            row_dict[label] = 100.0
            return row_dict

    row_dict[clone_text] = 100.0
    return row_dict


def _render_spatial_subplot(ax, plot_data, model, class_order, color_map, legend_ax=None):
    sns.scatterplot(
        data=plot_data,
        x="x",
        y="y",
        hue="_plot_label",
        hue_order=class_order,
        palette=color_map,
        s=8,
        edgecolor="none",
        alpha=0.95,
        legend=False,
        ax=ax,
    )

    model_title = _strip_preds_suffix(model)
    ax.set_title(f"{model_title} Detection", fontsize=14, fontweight="bold", pad=8)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(axis="both", length=0)
    for side in ["top", "right", "left", "bottom"]:
        ax.spines[side].set_visible(False)
    ax.axis("equal")


    if class_order:
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markerfacecolor=color_map.get(cls, "#999999"),
                markeredgecolor="none",
                markersize=6.5,
                label=str(cls),
            )
            for cls in class_order
        ]

        if legend_ax is not None:
            legend_ax.axis("off")
            legend_ax.legend(
                handles=handles,
                labels=[str(cls) for cls in class_order],
                title="Subclone Label",
                loc="upper left",
                bbox_to_anchor=(0.0, 1.0),
                frameon=False,
                markerscale=1.2,
                fontsize=8.0,
                title_fontsize=9.0,
            )
        else:
            ax.legend(
                handles=handles,
                labels=[str(cls) for cls in class_order],
                title="Subclone Label",
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                frameon=False,
                markerscale=1.5,
                fontsize=8.5,
                title_fontsize=9.5,
            )
