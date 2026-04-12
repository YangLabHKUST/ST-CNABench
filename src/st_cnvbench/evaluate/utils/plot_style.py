import logging

import matplotlib.colors as mcolors
import numpy as np
import matplotlib.pyplot as plt

from .constants import METHOD_FAMILY_BASE_COLORS, METHOD_VARIANT_SHADE


def apply_plot_style(dpi=300):
    """Apply a consistent publication-style matplotlib setup."""
    # Silence verbose font subsetting logs when exporting many PDF figures.
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    logging.getLogger("fontTools.subset").setLevel(logging.WARNING)

    plt.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "Arial",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _adjust_color_lightness(color_hex, factor):
    rgb = np.array(mcolors.to_rgb(color_hex))
    if factor >= 1.0:
        adjusted = rgb + (1.0 - rgb) * (factor - 1.0)
    else:
        adjusted = rgb * factor
    adjusted = np.clip(adjusted, 0.0, 1.0)
    return mcolors.to_hex(adjusted)


def resolve_method_color(model_name):
    model_lower = str(model_name).lower()

    family = None
    for candidate in METHOD_FAMILY_BASE_COLORS:
        if candidate in model_lower:
            family = candidate
            break

    if family is None:
        raise ValueError(f"Unsupported model name for color resolution: {model_name}")

    base_color = METHOD_FAMILY_BASE_COLORS[family]
    shade_factor = 1.0

    if "_expr" in model_lower or model_lower.endswith("expr"):
        shade_factor = METHOD_VARIANT_SHADE["expr"]
    elif "_cnv" in model_lower or model_lower.endswith("cnv"):
        shade_factor = METHOD_VARIANT_SHADE["cnv"]
    elif "nowgs" in model_lower or "no_wgs" in model_lower:
        shade_factor = METHOD_VARIANT_SHADE["nowgs"]
    elif "wgs" in model_lower:
        shade_factor = METHOD_VARIANT_SHADE["wgs"]

    return _adjust_color_lightness(base_color, shade_factor)
