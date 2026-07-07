import logging
import hashlib

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
    """Adjust a hex color by multiplying or lightening its RGB channels."""
    rgb = np.array(mcolors.to_rgb(color_hex))
    if factor >= 1.0:
        adjusted = rgb + (1.0 - rgb) * (factor - 1.0)
    else:
        adjusted = rgb * factor
    adjusted = np.clip(adjusted, 0.0, 1.0)
    return mcolors.to_hex(adjusted)


def _fallback_method_color(model_name):
    """Generate a stable fallback color for methods outside the known palette."""
    digest = hashlib.sha1(str(model_name).encode("utf-8")).hexdigest()
    hue = int(digest[:8], 16) / 0xFFFFFFFF
    rgb = mcolors.hsv_to_rgb((hue, 0.62, 0.72))
    return mcolors.to_hex(rgb)


def resolve_method_color(model_name):
    """Resolve a stable plot color for a known or custom method name."""
    model_lower = str(model_name).lower()

    family = None
    for candidate in METHOD_FAMILY_BASE_COLORS:
        if candidate in model_lower:
            family = candidate
            break

    if family is None:
        logging.warning("Using fallback color for unsupported model name: %s", model_name)
        return _fallback_method_color(model_name)

    base_color = METHOD_FAMILY_BASE_COLORS[family]
    shade_factor = 1.0

    if "_expr" in model_lower or model_lower.endswith("expr"):
        shade_factor = METHOD_VARIANT_SHADE["expr"]
    elif "_cna" in model_lower or model_lower.endswith("cna"):
        shade_factor = METHOD_VARIANT_SHADE["cna"]
    elif "nowgs" in model_lower or "no_wgs" in model_lower:
        shade_factor = METHOD_VARIANT_SHADE["nowgs"]
    elif "wgs" in model_lower:
        shade_factor = METHOD_VARIANT_SHADE["wgs"]

    return _adjust_color_lightness(base_color, shade_factor)
