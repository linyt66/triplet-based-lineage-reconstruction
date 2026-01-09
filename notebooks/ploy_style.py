"""
Matplotlib global style configuration for publication-quality figures.

This module defines a unified plotting style suitable for top-tier journals
(e.g., Systematic Biology, Nature, Science), with vector-friendly output,
compact typography, and consistent axis formatting.

Usage:
    import plot_style  # noqa: F401
    # All matplotlib figures created afterward will follow this style.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# -----------------------------------------------------------------------------
# Global matplotlib configuration (journal-ready, vector-first)
# -----------------------------------------------------------------------------

plt.rcParams.update({
    # --- Output format and font embedding (vector-friendly) ---
    'savefig.format': 'pdf',
    'figure.dpi': 300,
    'savefig.dpi': 600,
    'pdf.fonttype': 42,          # TrueType fonts for Illustrator compatibility
    'ps.fonttype': 42,
    'svg.fonttype': 'none',

    # --- Figure layout and size ---
    # Single-column: (3, 2); double-column: (5.5, 4); adjust as needed
    'figure.figsize': (3, 2),
    'figure.facecolor': 'white',
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'savefig.transparent': False,
    'figure.constrained_layout.use': False,

    # --- Font settings (sans-serif preferred by most journals) ---
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.unicode_minus': True,
    'text.usetex': False,

    # --- Math text and scientific notation ---
    'mathtext.fontset': 'dejavusans',
    'axes.formatter.use_mathtext': True,
    'axes.formatter.limits': (-3, 3),

    # --- Font sizes (compact, print-friendly) ---
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5,
    'legend.title_fontsize': 7,

    # --- Axes and ticks (boxed axes, inward ticks) ---
    'axes.spines.top': True,
    'axes.spines.right': True,
    'axes.linewidth': 0.6,
    'axes.labelpad': 2.0,

    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,

    'xtick.minor.visible': True,
    'ytick.minor.visible': True,

    'xtick.major.size': 4,
    'xtick.major.width': 0.6,
    'xtick.minor.size': 2,
    'xtick.minor.width': 0.4,
    'xtick.major.pad': 2,

    'ytick.major.size': 4,
    'ytick.major.width': 0.6,
    'ytick.minor.size': 2,
    'ytick.minor.width': 0.4,
    'ytick.major.pad': 2,

    # --- Lines, markers, and legends ---
    'lines.linewidth': 0.8,
    'lines.markersize': 3.0,
    'errorbar.capsize': 2.0,

    'legend.frameon': False,
    'legend.handlelength': 2.0,
    'legend.handleheight': 0.7,
    'legend.borderaxespad': 1,
    'legend.columnspacing': 2,
})


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def log10_formatter(x, pos):
    """
    Format logarithmic tick labels as base-10 exponents.

    Example:
        1e-2 -> -2
        1e0  -> 0
        1e2  -> 2
    """
    if x > 0:
        return f"{int(np.log10(x))}"
    return ""


# -----------------------------------------------------------------------------
# Standard color palette (colorblind-friendly, publication-safe)
# -----------------------------------------------------------------------------

COLORS = {
    'blue': '#1f77b4',
    'orange': '#ff7f0e',
    'green': '#2ca02c',
    'red': '#d62728',
    'purple': '#9467bd',
    'dark_gray': '#595959',
}
