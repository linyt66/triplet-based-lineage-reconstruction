"""Publication plotting defaults used by the manuscript figures."""

import math


COLORS = {
    "blue": "#1f77b4",
    "orange": "#ff7f0e",
    "green": "#2ca02c",
    "red": "#d62728",
    "purple": "#9467bd",
    "dark_gray": "#595959",
}


def configure_journal_style() -> None:
    """Apply compact, vector-friendly matplotlib defaults."""
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "savefig.format": "pdf",
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.figsize": (3, 2),
        "figure.facecolor": "white",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "savefig.transparent": False,
        "figure.constrained_layout.use": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.unicode_minus": True,
        "text.usetex": False,
        "mathtext.fontset": "dejavusans",
        "axes.formatter.use_mathtext": True,
        "axes.formatter.limits": (-3, 3),
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "legend.title_fontsize": 7,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.linewidth": 0.6,
        "axes.labelpad": 2.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": 4,
        "xtick.major.width": 0.6,
        "xtick.minor.size": 2,
        "xtick.minor.width": 0.4,
        "xtick.major.pad": 2,
        "ytick.major.size": 4,
        "ytick.major.width": 0.6,
        "ytick.minor.size": 2,
        "ytick.minor.width": 0.4,
        "ytick.major.pad": 2,
        "lines.linewidth": 0.8,
        "lines.markersize": 3.0,
        "errorbar.capsize": 2.0,
        "legend.frameon": False,
        "legend.handlelength": 2.0,
        "legend.handleheight": 0.7,
        "legend.borderaxespad": 1,
        "legend.columnspacing": 2,
    })


def log10_formatter(x, pos=None) -> str:
    """Format positive log-scale ticks as base-10 exponents."""
    if x > 0:
        return f"{int(math.log10(x))}"
    return ""
