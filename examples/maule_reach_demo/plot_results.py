"""Plot the hydroflux-solver-1d demo on the Río Maule reach.

Loads bed, depth, discharge GeoTIFFs and produces a two-panel longitudinal
profile figure for the Q4 review paper. Panel A: bed elevation + water
surface elevation along the channel. Panel B: Froude number along the
channel, with the critical line.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio

G = 9.81
OUT_DIR = Path(__file__).parent / "output"


def _load_row(path: Path) -> tuple[np.ndarray, float]:
    with rasterio.open(path) as src:
        data = src.read(1).astype(np.float64).ravel()
        dx = float(src.transform.a)
    return data, dx


def main() -> None:
    bed, dx = _load_row(OUT_DIR / "bed.tif")
    h, _ = _load_row(OUT_DIR / "depth.tif")
    hu, _ = _load_row(OUT_DIR / "discharge.tif")

    n = len(bed)
    x_km = (np.arange(n) + 0.5) * dx / 1_000.0
    water_surface = bed + h
    u = np.where(h > 1e-6, hu / h, 0.0)
    froude = np.where(h > 1e-6, np.abs(u) / np.sqrt(G * h), 0.0)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.5, 5.5), sharex=True,
        gridspec_kw={"height_ratios": [3, 2], "hspace": 0.12},
    )

    # ---- Panel A: longitudinal profile.
    ax_top.fill_between(x_km, bed, water_surface, color="#9bbcdc", alpha=0.85,
                        label="Water depth")
    ax_top.plot(x_km, water_surface, color="#1a4a72", lw=1.4,
                label="Water surface")
    ax_top.plot(x_km, bed, color="#5c4630", lw=1.4, label="Bed")
    ax_top.set_ylabel("Elevation (m)")
    ax_top.legend(loc="upper right", frameon=False, fontsize=9)
    ax_top.set_title(
        "Río Maule reach — hydroflux-solver-1d steady state "
        f"(n = {n}, dx = {dx:.1f} m, q = 3 m²/s, n = 0.04)",
        fontsize=10, loc="left",
    )

    # ---- Panel B: Froude.
    ax_bot.axhline(1.0, color="#a04040", lw=0.8, ls="--",
                   label="Critical (Fr = 1)")
    ax_bot.plot(x_km, froude, color="#1a4a72", lw=1.2)
    ax_bot.set_xlabel("Distance along channel (km)")
    ax_bot.set_ylabel("Froude number")
    ax_bot.set_ylim(0.0, max(1.05, froude.max() * 1.1))
    ax_bot.legend(loc="upper right", frameon=False, fontsize=9)

    # ---- Styling.
    for ax in (ax_top, ax_bot):
        ax.grid(True, alpha=0.25, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)

    fig_path = OUT_DIR / "figure.png"
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figure.pdf", bbox_inches="tight")
    print(f"Wrote {fig_path}")

    # ---- Compact numerical summary for the README / paper caption.
    print()
    print(f"Reach length:           {x_km[-1]:.2f} km")
    print(f"Bed drop:               {bed[0] - bed[-1]:.1f} m "
          f"(mean slope {(bed[0] - bed[-1]) / (x_km[-1] * 1000):.4f})")
    print(f"Depth range:            {h.min():.3f} – {h.max():.3f} m")
    print(f"Velocity range:         {u.min():.3f} – {u.max():.3f} m/s")
    print(f"Froude range:           {froude.min():.3f} – {froude.max():.3f}")
    print(f"Discharge (mean ± std): {hu.mean():.3f} ± {hu.std():.4f} m²/s")


if __name__ == "__main__":
    main()
