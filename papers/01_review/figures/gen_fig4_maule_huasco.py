"""Figure 4 — Flagship par insignia: Río Maule vs Río Huasco.

Paper-quality version of `examples/composite_figure.py`. Loads the
solver outputs from each demo's `output/` directory, renders the same
4-panel layout, applies the paper-wide style and adds an editorial
callout on the counter-intuitive Froude finding.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio

from style import (
    BASIN_COLORS,
    COLOR_BED,
    FIG_W_DC,
    add_panel_label,
    setup,
)

setup()

G = 9.81
REPO_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).parent


@dataclass
class Reach:
    name: str
    label: str
    demo_dir: str
    color: str
    subtitle: str


def load_reach(reach: Reach):
    out = REPO_ROOT / "examples" / reach.demo_dir / "output"
    with rasterio.open(out / "bed.tif") as src:
        bed = src.read(1).astype(np.float64).ravel()
        dx = float(src.transform.a)
    with rasterio.open(out / "depth.tif") as src:
        h = src.read(1).astype(np.float64).ravel()
    with rasterio.open(out / "discharge.tif") as src:
        hu = src.read(1).astype(np.float64).ravel()
    n = len(bed)
    x_km = (np.arange(n) + 0.5) * dx / 1_000.0
    u = np.where(h > 1e-6, hu / h, 0.0)
    fr = np.where(h > 1e-6, np.abs(u) / np.sqrt(G * h), 0.0)
    return x_km, bed, h, bed + h, u, fr


REACHES = [
    Reach(
        name="Maule",
        label="Río Maule (BNA #11)",
        demo_dir="maule_reach_demo",
        color=BASIN_COLORS["Maule"],
        subtitle=r"Mediterranean-temperate · slope $\approx 1\%$ · $q = 3$ m$^2$/s · $n = 0.04$",
    ),
    Reach(
        name="Huasco",
        label="Río Huasco (BNA #06)",
        demo_dir="huasco_reach_demo",
        color=BASIN_COLORS["Huasco"],
        subtitle=r"Semiarid Andean · slope $\approx 3.5\%$ · $q = 1$ m$^2$/s · $n = 0.06$",
    ),
]


def main() -> None:
    fig, axes = plt.subplots(
        2, 2,
        figsize=(FIG_W_DC, FIG_W_DC * 0.55),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 2.0],
                     "hspace": 0.16, "wspace": 0.18},
    )

    panel_labels = [("(a)", "(b)"), ("(c)", "(d)")]
    fr_max_global = 0.0

    for col, reach in enumerate(REACHES):
        x_km, bed, h, water, u, fr = load_reach(reach)
        ax_top = axes[0, col]
        ax_bot = axes[1, col]
        fr_max_global = max(fr_max_global, float(fr.max()))

        # Longitudinal profile.
        ax_top.fill_between(x_km, bed, water,
                            color=reach.color, alpha=0.22,
                            label="Water depth")
        ax_top.plot(x_km, water, color=reach.color, lw=1.3,
                    label="Water surface")
        ax_top.plot(x_km, bed, color=COLOR_BED, lw=1.3, label="Bed")
        ax_top.set_ylabel("Elevation (m)" if col == 0 else "", labelpad=2)
        # Solver setup as a low-stakes inline title (the panel label
        # carries the figure-paper identity; this just situates the run).
        ax_top.text(
            0.50, 0.96, f"{reach.label}\n{reach.subtitle}",
            transform=ax_top.transAxes,
            fontsize=8, ha="center", va="top",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=2),
        )
        if col == 0:
            ax_top.legend(loc="lower left", fontsize=7.5,
                          ncol=3, columnspacing=0.8, handlelength=1.2)
        add_panel_label(ax_top, panel_labels[0][col],
                        x=-0.13 if col == 0 else -0.08)

        # Froude.
        ax_bot.axhline(1.0, color="#a04040", lw=0.6, ls="--",
                       label="Critical $Fr = 1$")
        ax_bot.plot(x_km, fr, color=reach.color, lw=1.1)
        ax_bot.set_xlabel("Distance along channel (km)")
        ax_bot.set_ylabel("Froude number" if col == 0 else "", labelpad=2)
        ax_bot.set_ylim(0.0, 1.10)
        if col == 0:
            ax_bot.legend(loc="upper right", fontsize=7.5)
        add_panel_label(ax_bot, panel_labels[1][col],
                        x=-0.13 if col == 0 else -0.08)

        for ax in (ax_top, ax_bot):
            ax.grid(True, alpha=0.25, lw=0.4, color="#cccccc")

    # Editorial callout on the counter-intuitive Fr finding. Place it
    # in the Huasco bottom panel (col=1), pointing at the Fr plateau
    # which is visibly lower than the Maule oscillations on the left.
    ax_huasco_fr = axes[1, 1]
    ax_huasco_fr.annotate(
        "Despite a 3.5× steeper slope,\n"
        "Huasco has *lower* Froude than\n"
        r"Maule — boulder roughness ($n=0.06$)" "\n"
        "absorbs the extra gradient",
        xy=(5.0, 0.48), xytext=(2.5, 0.92),
        fontsize=7.5, color="#444444", style="italic",
        ha="left", va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2),
        arrowprops=dict(arrowstyle="-|>", color="#888888",
                        lw=0.5, mutation_scale=7,
                        connectionstyle="arc3,rad=-0.2"),
    )

    fig.savefig(OUT_DIR / "fig4_maule_huasco.png", dpi=300,
                bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig4_maule_huasco.pdf", bbox_inches="tight")
    print(f"Wrote {OUT_DIR / 'fig4_maule_huasco.png'}")
    print(f"Global max Froude observed: {fr_max_global:.3f}")


if __name__ == "__main__":
    main()
