"""Composite figure: hydroflux-solver-1d on Maule vs Huasco.

Side-by-side longitudinal profiles + Froude traces for the two pilot
basins of the postdoctorate. Intended as the "flagship figure" of the
Q4 2026 review paper: same solver, two contrasting climatic regimes,
spanning the Mediterranean-temperate / semiarid divide of the Chilean
piloto-cuencas.

Loads `maule_reach_demo/output/*.tif` and `huasco_reach_demo/output/*.tif`,
emits `figures/maule_vs_huasco.{png,pdf}` at the repository root level.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio

G = 9.81
ROOT = Path(__file__).parent
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


@dataclass
class Reach:
    name: str
    path_id: str
    color_water: str
    color_line: str
    title: str

    def load(self):
        out = ROOT / f"{self.path_id}_reach_demo" / "output"
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
        return x_km, bed, h, bed + h, u, fr, dx, n


MAULE = Reach(
    name="Río Maule",
    path_id="maule",
    color_water="#9bbcdc",
    color_line="#1a4a72",
    title="Río Maule — Mediterranean temperate (BNA #11)\n"
          "Andean foothills, slope 1.0 %, q = 3 m²/s, n = 0.04",
)
HUASCO = Reach(
    name="Río Huasco",
    path_id="huasco",
    color_water="#dcc69b",
    color_line="#8b6a1a",
    title="Río Huasco — semiarid Andean (BNA #06)\n"
          "Boulder bed, slope 3.5 %, q = 1 m²/s, n = 0.06",
)


def main() -> None:
    fig, axes = plt.subplots(
        2, 2, figsize=(12.5, 6.0), sharex=True,
        gridspec_kw={"height_ratios": [3, 2], "hspace": 0.18, "wspace": 0.18},
    )

    for col, reach in enumerate([MAULE, HUASCO]):
        x_km, bed, h, water, u, fr, dx, n = reach.load()
        ax_top = axes[0, col]
        ax_bot = axes[1, col]

        # ---- Longitudinal profile.
        ax_top.fill_between(
            x_km, bed, water, color=reach.color_water, alpha=0.85,
            label="Water depth",
        )
        ax_top.plot(x_km, water, color=reach.color_line, lw=1.4,
                    label="Water surface")
        ax_top.plot(x_km, bed, color="#5c4630", lw=1.4, label="Bed")
        ax_top.set_ylabel("Elevation (m)" if col == 0 else "")
        ax_top.set_title(reach.title, fontsize=10, loc="left")
        ax_top.legend(loc="upper right", frameon=False, fontsize=8)

        # ---- Froude.
        ax_bot.axhline(1.0, color="#a04040", lw=0.8, ls="--",
                       label="Critical (Fr = 1)")
        ax_bot.plot(x_km, fr, color=reach.color_line, lw=1.2)
        ax_bot.set_xlabel("Distance along channel (km)")
        ax_bot.set_ylabel("Froude number" if col == 0 else "")
        ax_bot.set_ylim(0.0, 1.1)
        ax_bot.legend(loc="upper right", frameon=False, fontsize=8)

        for ax in (ax_top, ax_bot):
            ax.grid(True, alpha=0.25, lw=0.5)
            ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "hydroflux-solver-1d on Chilean pilot basins — climatic contrast",
        fontsize=11, y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    png_path = OUT_DIR / "maule_vs_huasco.png"
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(OUT_DIR / "maule_vs_huasco.pdf", bbox_inches="tight")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
