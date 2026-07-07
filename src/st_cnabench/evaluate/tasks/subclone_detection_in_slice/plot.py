import pandas as pd
import os
import logging
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import re
from pathlib import Path
from PIL import Image, ImageEnhance
import squidpy as sq
import scanpy as sc
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram

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
        cmap='RdYlBu_r',
        vmin=0,
        vmax=1,
        linewidths=0.4,
        linecolor='white',
        cbar_kws={"shrink": 0.85, "label": metric_name},
        ax=ax
    )
    ax.set_title(f"{metric_name} Correlation Heatmap", fontsize=15, fontweight='bold', pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis='x', labelsize=10, rotation=30)
    ax.tick_params(axis='y', labelsize=10, rotation=0)
    fig.tight_layout()

    file_name = f"{metric_name.lower().replace('-', '_')}_heatmap.png"
    save_path = save_dir / file_name
    fig.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
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
        edgecolor='black',
        linewidth=0.9,
        width=0.66,
    )

    for b in bars:
        h = b.get_height()
        if pd.notna(h):
            ax.annotate(f"{h:.2f}", (b.get_x() + b.get_width() / 2.0, h),
                        ha='center', va='bottom', fontsize=9, fontweight='bold',
                        xytext=(0, 3), textcoords='offset points')

    ax.set_ylabel(f'{metric_name} Score', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Model Performance vs Ground Truth ({metric_name})", fontsize=15, fontweight='bold', pad=12)
    ax.set_xlabel("")
    ax.tick_params(axis='x', labelsize=10, rotation=35)
    ax.tick_params(axis='y', labelsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()

    file_name = f"{metric_name.lower().replace('-', '_')}_barplot.png"
    save_path = save_dir / file_name
    fig.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    logging.info(f"Saved {metric_name} barplot to {save_path}")
    plt.close(fig)

def run_plot_spatial_coherence(spatial_coherence, save_dir):
    """
    spatial_coherence: dict of model_name -> {'observed_score': observed_score, 'z_score': z_score}
    Plot barplots of spatial coherence z-score and observed score.
    """
    model_names = list(spatial_coherence.keys())
    if not model_names:
        return

    apply_plot_style(dpi=300)
    sns.set_style("ticks")
    colors = [resolve_method_color(m) for m in model_names]

    z_scores = [spatial_coherence[model]['z_score'] for model in model_names]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=300)
    ax.bar(model_names, z_scores, color=colors, edgecolor='black', linewidth=0.9)
    ax.axhline(0, color='#666666', linestyle='--', linewidth=0.9, alpha=0.8)
    ax.set_ylabel('Spatial Coherence Z-score', fontsize=11, fontweight='bold')
    ax.set_title('Spatial Coherence Comparison', fontsize=14, fontweight='bold', pad=10)
    ax.tick_params(axis='x', labelsize=10, rotation=35)
    ax.tick_params(axis='y', labelsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_path = save_dir / "spatial_coherence_barplot.png"
    fig.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    logging.info(f"Saved spatial coherence barplot to {save_path}")
    plt.close(fig)

    observed_scores = [spatial_coherence[model]['observed_score'] for model in model_names]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=300)
    ax.bar(model_names, observed_scores, color=colors, edgecolor='black', linewidth=0.9)
    ax.set_ylabel('Spatial Coherence Observed Score', fontsize=11, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.set_title('Spatial Coherence Observed Score Comparison', fontsize=14, fontweight='bold', pad=10)
    ax.tick_params(axis='x', labelsize=10, rotation=35)
    ax.tick_params(axis='y', labelsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_path = save_dir / "spatial_coherence_observed_score_barplot.png"
    fig.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    logging.info(f"Saved spatial coherence observed score barplot to {save_path}")
    plt.close(fig)

def run_spatial_plot_no_img(results, save_path):
    """Render slice-wide spatial scatter plots from the raw label table."""
    save_path = Path(save_path)
    panel_dir = save_path.parent / f"{save_path.stem}_panels"
    panel_dir.mkdir(parents=True, exist_ok=True)

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
            ax.axis('off')
            continue

        _render_spatial_subplot(ax, plot_data, model, class_order, color_map)

        fig_single, ax_single = plt.subplots(figsize=(6.0, 5.0), dpi=300)
        _render_spatial_subplot(ax_single, plot_data, model, class_order, color_map)
        fig_single.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.06)
        safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", _strip_preds_suffix(model))
        panel_path = panel_dir / f"{i + 1:02d}_{safe_name}.png"
        fig_single.savefig(panel_path, dpi=300, facecolor="white")
        plt.close(fig_single)

    for j in range(n_plots, len(axes_flat)):
        axes_flat[j].axis("off")

    fig.tight_layout()
    fig.savefig(save_path, dpi=300, facecolor="white")
    logging.info(f"Saved slice-wide spatial scatter plot to {save_path}")
    logging.info(f"Saved per-panel spatial plots to {panel_dir}")
    plt.close(fig)

def run_plot_clone_matching_heatmap(matrix_dict, metric_name, save_dir):
    """Render GT-vs-model clone matching heatmaps for genomic-view metrics."""
    if not matrix_dict:
        return

    apply_plot_style(dpi=300)
    out_dir = save_dir / f"clone_matching_{metric_name.lower()}"
    os.makedirs(out_dir, exist_ok=True)

    for model, matrix in matrix_dict.items():
        if matrix is None or matrix.empty:
            continue

        fig, ax = plt.subplots(
            figsize=(max(6.2, len(matrix.columns) * 1.15), max(4.2, len(matrix.index) * 0.8)),
            dpi=300,
        )

        if 'PCC' in metric_name:
            cmap = 'RdBu_r'
            center = 0
            vmin, vmax = -1, 1
        else:
            cmap = 'YlGnBu'
            center = None
            vmin, vmax = 0, 1

        sns.heatmap(
            matrix.astype(float),
            annot=True,
            fmt=".2f",
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            linewidths=0.4,
            linecolor='white',
            cbar_kws={"shrink": 0.82},
            ax=ax,
        )
        ax.set_title(f"{model} vs GT ({metric_name})", fontsize=14, fontweight='bold', pad=10)
        ax.set_ylabel("Ground Truth Subclones", fontsize=11, fontweight='bold')
        ax.set_xlabel(f"{model} Predicted Clones", fontsize=11, fontweight='bold')
        ax.tick_params(axis='x', labelsize=9, rotation=30)
        ax.tick_params(axis='y', labelsize=9, rotation=0)
        fig.tight_layout()

        save_path = out_dir / f"{model}_{metric_name}_heatmap.png"
        fig.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
    logging.info(f"Saved {metric_name} clone matching heatmaps to {out_dir}")

def run_plot_karyotype_composition(df_all, clone_cna_profiles, save_dir):
    """Render the unified clustered karyogram with GT composition bars."""
    apply_plot_style(dpi=300)

    if 'GT_Preds' not in df_all.columns or not clone_cna_profiles:
        logging.warning("Missing GT_Preds or CNA profiles, skipping unified karyotype composition plot.")
        return

    os.makedirs(save_dir, exist_ok=True)
    all_gt_classes = [str(x) for x in df_all["GT_Preds"].dropna().unique()]
    if not all_gt_classes:
        logging.warning("GT_Preds has no valid labels, skipping unified karyotype composition plot.")
        return

    all_cna_list = []
    for mod, df in clone_cna_profiles.items():
        if df.empty:
            continue
        df = df.copy()

        if 'Chr' in df.columns and 'Chromosome' not in df.columns:
            df = df.rename(columns={'Chr': 'Chromosome'})

        score_col = 'CN_Score' if 'CN_Score' in df.columns else 'TCN_Score'
        if score_col in df.columns:
            df['Unified_Score'] = df[score_col]
            df['Clone_ID'] = df['Clone_ID'].map(lambda x: _format_clone_name(mod, x))
            df = df.dropna(subset=['Clone_ID']).copy()

            req_cols = ['Clone_ID', 'Chromosome', 'Start', 'End', 'Unified_Score']
            if all(c in df.columns for c in req_cols):
                all_cna_list.append(df[req_cols])

    if not all_cna_list:
        logging.warning("No valid CNA profiles found for plotting.")
        return

    unified_cna_df = pd.concat(all_cna_list, ignore_index=True)

    unified_cna_df['Bin_ID'] = unified_cna_df['Chromosome'].astype(str) + '_' + unified_cna_df['Start'].astype(str) + '_' + unified_cna_df['End'].astype(str)

    def _get_sort_key(ch):
        ch = str(ch).replace('chr', '').replace('Chr', '')
        ch_map = {'X': 23, 'Y': 24, 'x': 23, 'y': 24}
        return ch_map.get(ch, int(ch) if ch.isdigit() else 99)

    unified_cna_df['chr_sort'] = unified_cna_df['Chromosome'].map(_get_sort_key)
    unified_cna_df['Start'] = unified_cna_df['Start'].astype(float)
    unified_cna_df = unified_cna_df.sort_values(['chr_sort', 'Start'])

    cna_matrix = unified_cna_df.pivot(index='Clone_ID', columns='Bin_ID', values='Unified_Score')
    ordered_bins = unified_cna_df['Bin_ID'].unique()
    cna_matrix = cna_matrix[ordered_bins]

    for clone in cna_matrix.index:
        std_val = cna_matrix.loc[clone].std()
        mean_val = cna_matrix.loc[clone].mean()
        if std_val == 0 or pd.isna(std_val):
            std_val = 1
        cna_matrix.loc[clone] = (cna_matrix.loc[clone] - mean_val) / std_val

    cna_matrix = cna_matrix.fillna(0)

    Z = None
    if len(cna_matrix) > 1:
        Z = linkage(cna_matrix.values, method='ward')
        idx = leaves_list(Z)
        ordered_clones = cna_matrix.index[idx].tolist()
    else:
        ordered_clones = cna_matrix.index.tolist()

    cna_matrix = cna_matrix.loc[ordered_clones]

    all_cross_tabs = []
    all_gt_labels = [str(x) for x in df_all['GT_Preds'].dropna().unique()]

    for mod in clone_cna_profiles.keys():
        if mod == 'GT':
            gt_clones = clone_cna_profiles['GT']['Clone_ID'].unique()
            cross_data = {}
            for clone in gt_clones:
                row_name = _format_clone_name(mod, clone)
                if row_name is None:
                    continue
                cross_data[row_name] = _build_gt_identity_row(clone, all_gt_labels)
            if cross_data:
                cross = pd.DataFrame.from_dict(cross_data, orient='index')
                all_cross_tabs.append(cross)
            continue

        mod_col = f"{mod}_Preds" if f"{mod}_Preds" in df_all.columns else mod

        if mod_col not in df_all.columns:
            logging.warning(f"Model column {mod_col} not found in dataframe, skipping purity calculation for {mod}.")
            continue

        valid_df = df_all.dropna(subset=['GT_Preds', mod_col]).copy()
        if valid_df.empty:
            continue

        valid_df[mod_col] = valid_df[mod_col].astype(str)
        cross = pd.crosstab(valid_df[mod_col], valid_df['GT_Preds'])
        cross.index = [_format_clone_name(mod, c) for c in cross.index]
        cross = cross[cross.index.notna()]
        all_cross_tabs.append(cross)

    if not all_cross_tabs:
        logging.error("Failed to generate cross-tabs for GT purity (check column names). Unified karyogram plot skipped.")
        return

    master_cross_tab = pd.concat(all_cross_tabs).fillna(0)
    master_cross_tab = master_cross_tab.groupby(master_cross_tab.index).sum()

    cross_tab_norm = master_cross_tab.div(master_cross_tab.sum(axis=1), axis=0)
    cols_order = sorted(cross_tab_norm.columns.tolist(), key=_gt_class_sort_key)
    gt_color_map = _build_categorical_color_map(cols_order)
    gt_colors = [gt_color_map[c] for c in cols_order]

    comp_data = cross_tab_norm.reindex(ordered_clones).fillna(0)
    comp_data = comp_data[cols_order]

    row_height = 0.22
    base_height = 3.5
    fig_height = max(5.0, base_height + len(ordered_clones) * row_height)

    fig = plt.figure(figsize=(16, fig_height), dpi=300)

    # Narrow the spacer so dendrogram branches sit closer to method labels.
    gs = gridspec.GridSpec(1, 5, width_ratios=[1.0, 1.4, 7.3, 0.1, 1.5], wspace=0.0)

    ax_dendro = fig.add_subplot(gs[0])
    if Z is not None:
        dendrogram(Z, orientation='left', ax=ax_dendro, color_threshold=0, above_threshold_color='black')
        ax_dendro.set_ylim(len(ordered_clones) * 10, 0)
    ax_dendro.axis('off')

    ax_cna = fig.add_subplot(gs[2])

    bottom_margin = max(0.12, 1.8 / fig_height)
    cbar_ax = fig.add_axes([0.35, bottom_margin / 3, 0.3, 0.015])

    sns.heatmap(
        cna_matrix,
        cmap='RdBu_r',
        center=0, vmin=-5, vmax=5,
        ax=ax_cna,
        cbar_ax=cbar_ax,
        xticklabels=False,
        yticklabels=True,
        cbar_kws={'label': 'Standardized Log2 Ratio', 'orientation': 'horizontal'}
    )
    cbar_ax.set_xticks([-5, 0, 5])
    cbar_ax.set_xticklabels(['Loss', 'Neutral', 'Gain'])

    chroms = [str(bid).split('_')[0] for bid in cna_matrix.columns]
    boundaries = [0]
    labels = []
    curr_chr = None
    for i, c in enumerate(chroms):
        if c != curr_chr:
            if curr_chr is not None: boundaries.append(i)
            labels.append(c.lower().replace('chr', ''))
            curr_chr = c
    boundaries.append(len(cna_matrix.columns))

    for b in boundaries:
        ax_cna.axvline(x=b, color="#4a4a4a", linewidth=0.70, alpha=0.40)

    ax_cna.add_patch(
        plt.Rectangle(
            (0, 0),
            cna_matrix.shape[1],
            cna_matrix.shape[0],
            fill=False,
            edgecolor="#5a5a5a",
            linewidth=0.55,
            alpha=0.45,
            clip_on=False,
        )
    )

    label_pos = [(boundaries[i] + boundaries[i+1]) / 2 for i in range(len(boundaries)-1)]
    ax_cna.set_xticks(label_pos)

    chr_fontsize = 10
    ax_cna.set_xticklabels(labels, fontsize=chr_fontsize, fontweight="normal", rotation=0)
    ax_cna.tick_params(axis='x', which='both', length=0)

    y_fontsize = 12 if len(ordered_clones) <= 40 else 8
    ax_cna.set_yticklabels(ax_cna.get_yticklabels(), fontweight="normal", rotation=0, fontsize=y_fontsize)
    ax_cna.tick_params(axis="y", length=0)
    ax_cna.xaxis.set_ticks_position("bottom")
    ax_cna.xaxis.set_label_position("bottom")
    ax_cna.yaxis.set_ticks_position("left")
    ax_cna.set_title("Unified Clustered Karyogram Profiles", fontsize=16, fontweight='bold', pad=15)
    ax_cna.set_ylabel("")
    ax_cna.set_xlabel("Chromosome", fontsize=11)

    ax_bar = fig.add_subplot(gs[4], sharey=ax_cna)

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
            edgecolor='white',
            linewidth=0.5,
            label=col
        )
        left_bases += widths

    ax_bar.set_xlim(0, 1.0)
    ax_bar.set_title("Tumor Composition", fontsize=14, fontweight='bold', pad=15)

    ax_bar.set_facecolor('white')
    ax_bar.grid(False)
    for spine_loc in ['top', 'right', 'left']:
        ax_bar.spines[spine_loc].set_visible(False)
    ax_bar.spines['bottom'].set_color('black')

    ax_bar.tick_params(axis='y', left=False, labelleft=False)

    handles, labels = ax_bar.get_legend_handles_labels()
    legend_cols = 2 if len(cols_order) > 3 else len(cols_order)

    fig.legend(
        handles, labels,
        title="Tumor Class",
        bbox_to_anchor=(0.88, bottom_margin / 5),
        loc='lower center',
        ncol=legend_cols,
        frameon=False,
        fontsize=10
    )

    plt.subplots_adjust(bottom=bottom_margin, top=0.92, left=0.02, right=0.96)

    save_file = os.path.join(save_dir, "Unified_All_Methods_Clustered_Karyogram.png")
    plt.savefig(save_file, bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()

    logging.info(f"Saved strictly aligned unified karyogram with absolute legends to {save_file}")

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


def _prepare_spatial_plot_data(results, model):
    if model not in results.columns:
        return None, [], {}

    plot_data = results.dropna(subset=[model, "x", "y"]).copy()
    if plot_data.empty:
        return plot_data, [], {}

    plot_data["_plot_label"] = plot_data[model].map(lambda x: _simplify_model_label(model, x))
    plot_data = plot_data.dropna(subset=["_plot_label"]).copy()
    class_order = sorted(plot_data["_plot_label"].dropna().unique().tolist(), key=_subclone_label_sort_key)
    color_map = _build_categorical_color_map(class_order)
    return plot_data, class_order, color_map


def _render_spatial_subplot(ax, plot_data, model, class_order, color_map):
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
