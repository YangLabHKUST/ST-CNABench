import logging
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from itertools import cycle


def run_plot_resolution(results, plot_path):
    plot_path = Path(plot_path)
    plot_path.parent.mkdir(parents=True, exist_ok=True)

    # Pattle choose!!
    professional_colors = [
        "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
        "#F39B7FFF", "#8491B3FF", "#91D1C2FF", "#DC0000FF",
        "#7E6148FF", "#B09C85FF", "#BC3C29FF", "#0072B5FF"
    ]
    colors = cycle(professional_colors)

    bins = np.logspace(4, 9, num=25)
    xticks = [1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
    xtick_labels = ["10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb", "1 Gb"]

    num_models = len(results)
    if num_models == 0:
        logging.warning("No resolution results found. Skip plotting.")
        return

    # two columns layout
    ncols = 2
    nrows = (num_models + 1) // 2

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(14, 3.2 * nrows),
        sharex=True,
        sharey=True
    )

    # Always normalize to a flat Axes list for both 1-model and multi-model cases.
    axes_flat = np.atleast_1d(axes).flatten()

    for (model_name, loaders_data), ax in zip(results.items(), axes_flat):
        color = next(colors)

        all_lengths = []
        for df in loaders_data.values():
            if df is not None and not df.empty:
                all_lengths.append(df["Length"].values)

        if not all_lengths:
            continue

        lengths = np.concatenate(all_lengths)
        lengths = lengths[lengths > 0]
        if lengths.size == 0:
            continue

        median_len = np.median(lengths)

        hist_counts, bin_edges = np.histogram(lengths, bins=bins)
        hist_sum = hist_counts.sum()
        if hist_sum == 0:
            continue

        hist_percent = hist_counts / hist_sum * 100
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_widths = bin_edges[1:] - bin_edges[:-1]

        ax.bar(
            bin_centers, hist_percent, width=bin_widths, align="center",
            color=color, alpha=0.3, edgecolor="black", linewidth=0.5, zorder=1
        )
        ax.plot(
            bin_centers, hist_percent, color=color, lw=2,
            marker="o", markersize=3, zorder=3
        )
        # Median line
        ax.axvline(
            median_len, linestyle="--", color=color, alpha=0.8,
            linewidth=2, zorder=2
        )
        info_title = f"{model_name} (Median: {median_len/1e6:.2f} Mb)"

        ax.set_title(
            info_title,
            #color=color,
            fontsize=12,
            fontweight="bold",
            pad=10
        )

        ax.set_xscale("log")
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels)
        ax.tick_params(labelbottom=True, which="both", direction="in", length=0)

        ax.set_ylabel("Percentage (%)", fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for ax in axes_flat[num_models:]:
        ax.set_visible(False)

    for ax in axes_flat[:num_models]:
        ax.set_xlabel("CNA Segment Length", fontsize=12)

    plt.suptitle("CNA Segment Length Distribution (Resolution Comparison)", fontsize=14, y=0.99)

    plt.tight_layout(rect=[0, 0, 1, 0.98])

    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    logging.info(f"Final resolution plot with unified labels saved at {plot_path}")
