"""Figure 3 — MacDonald variable-depth steady state."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from style import (
    COLOR_ANALYTICAL,
    COLOR_BED,
    COLOR_NUMERICAL,
    COLOR_REFERENCE,
    FIG_W_DC,
    add_panel_label,
    setup,
)

setup()
OUT_DIR = Path(__file__).parent

G = 9.81
Q = 1.0
N_MANNING = 0.03
L = 50.0
H_BASE = 1.0
H_AMP = 0.2

L1_TABLE = [
    # (n, dx, L1(h), L1(hu))
    (50, 1.000, 0.39196, 3.21456),
    (100, 0.500, 0.18464, 1.63426),
    (200, 0.250, 0.08913, 0.82122),
    (400, 0.125, 0.04382, 0.41134),
    (800, 0.063, 0.02174, 0.20584),
]


def h_analytical(x):
    return H_BASE + H_AMP * np.sin(2.0 * np.pi * x / L)


def dh_dx_analytical(x):
    return H_AMP * (2.0 * np.pi / L) * np.cos(2.0 * np.pi * x / L)


def friction_slope(h):
    return N_MANNING ** 2 * Q * Q / h ** (10.0 / 3.0)


def dz_dx_analytical(x):
    h = h_analytical(x)
    fr_sq = Q * Q / (G * h ** 3)
    return -(1.0 - fr_sq) * dh_dx_analytical(x) - friction_slope(h)


def bed_profile(n_pts):
    x = np.linspace(0.0, L, n_pts + 1)
    dz = 0.5 * (dz_dx_analytical(x[:-1]) + dz_dx_analytical(x[1:])) * (L / n_pts)
    z = np.concatenate([[0.0], np.cumsum(dz)])
    return x, z


def main() -> None:
    fig, (ax_l, ax_r) = plt.subplots(
        1, 2, figsize=(FIG_W_DC, FIG_W_DC * 0.42),
    )

    # ── Panel (a): inverse design.
    x_fine, z = bed_profile(1000)
    h = h_analytical(x_fine)
    water_surface = z + h

    # Water column fill — light tint of the numerical-output blue.
    ax_l.fill_between(x_fine, z, water_surface,
                      color=COLOR_NUMERICAL, alpha=0.18,
                      label="Water depth $h(x)$")
    ax_l.plot(x_fine, water_surface,
              color=COLOR_NUMERICAL, lw=1.4,
              label=r"Water surface $\eta = z + h$")
    ax_l.plot(x_fine, z,
              color=COLOR_BED, lw=1.4,
              label="Derived bed $z(x)$")
    ax_l.set_xlabel(r"$x$ (m)")
    ax_l.set_ylabel("Elevation (m)")
    ax_l.set_xlim(0, L)
    ax_l.legend(loc="lower left", fontsize=8)
    add_panel_label(ax_l, "(a)")

    # Editorial inline note explaining the inverse-design construction
    # (placed in upper-right where there is whitespace).
    ax_l.text(
        0.97, 0.95,
        r"$h(x) = 1 + 0.2\sin(2\pi x/L)$ prescribed;""\n"
        r"$z(x)$ derived from $dz/dx = -(1{-}Fr^2)\,dh/dx - S_f$",
        transform=ax_l.transAxes,
        fontsize=7.5, color="#444444", style="italic",
        ha="right", va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2),
    )

    # ── Panel (b): convergence.
    dxs = np.array([row[1] for row in L1_TABLE])
    l1_h = np.array([row[2] for row in L1_TABLE])
    l1_hu = np.array([row[3] for row in L1_TABLE])

    slope_h, _ = np.polyfit(np.log(dxs), np.log(l1_h), 1)
    slope_hu, _ = np.polyfit(np.log(dxs), np.log(l1_hu), 1)

    ax_r.loglog(dxs, l1_h, marker="o", color=COLOR_NUMERICAL, lw=1.3,
                label=fr"$L^1(h)$, order = {slope_h:.2f}")
    ax_r.loglog(dxs, l1_hu, marker="s", color=COLOR_ANALYTICAL, lw=1.3,
                label=fr"$L^1(hu)$, order = {slope_hu:.2f}")
    ref_x = np.array([dxs.min() * 0.7, dxs.max() * 1.3])
    ref_y = ref_x * (l1_h[0] / dxs[0])
    ax_r.loglog(ref_x, ref_y, color=COLOR_REFERENCE, ls="--", lw=0.6,
                label="Order-1 reference")
    ax_r.set_xlabel(r"Cell size $\Delta x$ (m)")
    ax_r.set_ylabel(r"$L^1$ error")
    ax_r.legend(loc="lower right", fontsize=8)
    ax_r.grid(True, which="both", color="#dddddd", lw=0.3, alpha=0.6)
    add_panel_label(ax_r, "(b)", x=-0.14)

    # Editorial inline note explaining the contrast with Stoker.
    ax_r.text(
        0.04, 0.96,
        "Clean first-order:\nno shock to smear\n(cf. Fig. 2)",
        transform=ax_r.transAxes,
        fontsize=7.5, color="#444444", style="italic",
        ha="left", va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2),
    )

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig3_macdonald.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig3_macdonald.pdf", bbox_inches="tight")
    print(f"Wrote {OUT_DIR / 'fig3_macdonald.png'}")
    print(f"L1(h) order = {slope_h:.3f}, L1(hu) order = {slope_hu:.3f}")


if __name__ == "__main__":
    main()
