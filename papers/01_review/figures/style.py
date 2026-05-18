"""Paper figure style — shared across all hydroflux 01_review figures.

NHESS / Copernicus targeting:
- Single-column width 88 mm (3.46 in), two-column 170 mm (6.7 in).
- Sans-serif typography (Copernicus uses Latin Modern Sans at typeset);
  figures shipped with sans-serif fallback for consistency.
- Wong palette (Nature Methods 2011) for categorical colour,
  colourblind-safe.
- Subtle grid (alpha 0.3, lw 0.4) only where it aids reading.
- Top and right spines off; ticks pointing inward.

Per-entity colours are defined ONCE here and consumed by every gen_*.py
so that a solver or basin appearing in multiple figures keeps the same
colour identity across the paper.

Import at the top of each generator:

    from style import (
        setup, add_panel_label,
        SOLVER_COLORS, BASIN_COLORS, COLOR_WONG,
        FIG_W_SC, FIG_W_DC,
    )
    setup()
"""

from __future__ import annotations

import matplotlib as mpl

# ── Wong palette (Nature Methods 2011), colourblind-safe categorical ──
COLOR_WONG = {
    "orange":    "#E69F00",
    "sky":       "#56B4E9",
    "green":     "#009E73",
    "yellow":    "#F0E442",
    "blue":      "#0072B2",
    "vermilion": "#D55E00",
    "purple":    "#CC79A7",
    "black":     "#000000",
}

# ── Per-entity semantic colours used paper-wide for cross-fig consistency ──

# Solvers (Figure 1 intersection diagram).
SOLVER_COLORS = {
    "HEC-RAS":     COLOR_WONG["vermilion"],  # regulatory anchor
    "LISFLOOD-FP": COLOR_WONG["sky"],         # open-source CUDA
    "BASEMENT":    COLOR_WONG["orange"],      # closed academic
    "TELEMAC":     COLOR_WONG["yellow"],      # legacy LGPL
    "TUFLOW HPC":  COLOR_WONG["green"],       # commercial GPU
    "hydroflux":   COLOR_WONG["blue"],        # the target
}

# Chilean pilot basins (Figure 4 flagship + future per-basin plots).
BASIN_COLORS = {
    "Maule":  COLOR_WONG["blue"],     # Mediterranean-temperate, navy
    "Huasco": COLOR_WONG["orange"],   # semiarid Andean, warm
}

# ── Numerical / analytical solution colours (Figures 2, 3) ──
COLOR_NUMERICAL  = COLOR_WONG["blue"]        # solver output
COLOR_ANALYTICAL = COLOR_WONG["vermilion"]   # closed-form reference
COLOR_BED        = "#5c4630"                  # earth-tone for bathymetry
COLOR_REFERENCE  = "#8a8a8a"                  # grey for order-1 / critical lines

# ── Journal widths (NHESS / Copernicus, mm → inches) ──
FIG_W_SC = 88.0  / 25.4   # single-column ≈ 3.46 in
FIG_W_DC = 170.0 / 25.4   # two-column   ≈ 6.69 in


def setup(serif: bool = False) -> None:
    """Apply global matplotlib rcParams. Call once at the top of every
    generator script, before any ``plt`` calls."""
    family = "serif" if serif else "sans-serif"
    mpl.rcParams.update({
        # output
        "figure.dpi":      130,
        "savefig.dpi":     300,
        "savefig.bbox":    "tight",
        "savefig.transparent": False,

        # typography
        "font.family":       family,
        "font.size":         9,
        "axes.titlesize":    10,
        "axes.labelsize":    9,
        "xtick.labelsize":   8,
        "ytick.labelsize":   8,
        "legend.fontsize":   8,

        # spines
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.linewidth":    0.6,

        # ticks
        "xtick.direction":   "in",
        "ytick.direction":   "in",
        "xtick.major.size":  3,
        "ytick.major.size":  3,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,

        # grid — disabled by default; enable per-fig with ax.grid(True, ...)
        "axes.grid":        False,
        "grid.alpha":       0.3,
        "grid.linewidth":   0.4,
        "grid.color":       "#cccccc",

        # legend
        "legend.frameon":     False,
        "legend.handlelength": 1.6,

        # lines / markers
        "lines.linewidth": 1.4,
        "lines.markersize": 4,
    })


def add_panel_label(
    ax,
    label: str,
    x: float = -0.12,
    y: float = 1.05,
    fontsize: float = 10,
) -> None:
    """Add a (a)/(b)/(c) panel label top-left of an axes, in bold."""
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=fontsize, fontweight="bold",
        va="top", ha="left",
    )
