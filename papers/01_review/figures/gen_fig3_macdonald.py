"""Figure 3 — MacDonald variable-depth steady state.

Two-panel figure:
- Left: prescribed depth profile h(x) = 1 + 0.2 sin(2π x/L) overlaid on
  the derived bed z(x) = ∫ dz/dx dx', with the water surface η = z + h
  shown for context. Demonstrates the inverse-design construction.
- Right: log-log L1(h) convergence across n ∈ {50, 100, 200, 400, 800}.
  Reference order-1 slope; empirical slope 1.03 annotated.

Data: benchmarks/macdonald-variable-results.md table.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

G = 9.81
Q = 1.0
N_MANNING = 0.03
L = 50.0
H_BASE = 1.0
H_AMP = 0.2

# From benchmarks/macdonald-variable-results.md
L1_TABLE = [
    # (n, dx, L1(h), L1(hu), rel L1, ratio)
    (50, 1.000, 0.39196, 3.21456, 0.00784, None),
    (100, 0.500, 0.18464, 1.63426, 0.00369, 2.12),
    (200, 0.250, 0.08913, 0.82122, 0.00178, 2.07),
    (400, 0.125, 0.04382, 0.41134, 0.00088, 2.03),
    (800, 0.063, 0.02174, 0.20584, 0.00043, 2.02),
]

OUT_DIR = Path(__file__).parent


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


def bed_profile(n_pts: int) -> tuple[np.ndarray, np.ndarray]:
    """Trapezoid-integrate z(x) starting from z(0) = 0.
    Returns (x_nodes, z_values)."""
    x = np.linspace(0.0, L, n_pts + 1)
    dz = 0.5 * (dz_dx_analytical(x[:-1]) + dz_dx_analytical(x[1:])) * (L / n_pts)
    z = np.concatenate([[0.0], np.cumsum(dz)])
    return x, z


def main() -> None:
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(10.0, 4.0))

    # ---- Left panel: profile.
    x_fine, z = bed_profile(1000)
    h = h_analytical(x_fine)
    water_surface = z + h

    ax_l.fill_between(x_fine, z, water_surface, color="#9bbcdc", alpha=0.75,
                      label="Water depth $h(x)$")
    ax_l.plot(x_fine, water_surface, color="#1a4a72", lw=1.4,
              label=r"Water surface $\eta = z + h$")
    ax_l.plot(x_fine, z, color="#5c4630", lw=1.4,
              label="Derived bed $z(x)$")
    ax_l.set_xlabel(r"$x$ (m)")
    ax_l.set_ylabel("Elevation (m)")
    ax_l.set_xlim(0, L)
    ax_l.legend(loc="lower left", frameon=False, fontsize=9)
    ax_l.set_title(
        r"(a) Inverse design: prescribe $h(x)$, derive $z(x)$",
        loc="left", fontsize=10,
    )

    # Annotate with the analytical formula.
    ax_l.text(0.5, 0.95,
              r"$h(x) = 1.0 + 0.2 \sin(2\pi x / L)$,  $q = 1$ m$^2$/s,  $n = 0.03$",
              transform=ax_l.transAxes,
              ha="center", va="top", fontsize=9,
              bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=3))

    # ---- Right panel: convergence.
    dxs = np.array([row[1] for row in L1_TABLE])
    l1_h = np.array([row[2] for row in L1_TABLE])
    l1_hu = np.array([row[3] for row in L1_TABLE])

    slope_h, intercept_h = np.polyfit(np.log(dxs), np.log(l1_h), 1)
    slope_hu, intercept_hu = np.polyfit(np.log(dxs), np.log(l1_hu), 1)

    ax_r.loglog(dxs, l1_h, marker="o", color="#1a4a72", lw=1.4,
                label=fr"$L^1(h)$, order = {slope_h:.2f}")
    ax_r.loglog(dxs, l1_hu, marker="s", color="#5c4630", lw=1.4,
                label=fr"$L^1(hu)$, order = {slope_hu:.2f}")
    ref_x = np.array([dxs.min() * 0.7, dxs.max() * 1.3])
    ref_y = ref_x * (l1_h[0] / dxs[0])
    ax_r.loglog(ref_x, ref_y, color="#888888", ls="--", lw=0.8,
                label="Order-1 reference")
    ax_r.set_xlabel(r"Cell size $\Delta x$ (m)")
    ax_r.set_ylabel(r"$L^1$ error")
    ax_r.legend(loc="lower right", frameon=False, fontsize=9)
    ax_r.set_title(
        r"(b) Convergence on $n \in \{50, 100, 200, 400, 800\}$",
        loc="left", fontsize=10,
    )

    for ax in (ax_l, ax_r):
        ax.grid(True, alpha=0.25, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig3_macdonald.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig3_macdonald.pdf", bbox_inches="tight")
    print(f"Wrote {OUT_DIR / 'fig3_macdonald.png'}")
    print(f"L1(h) empirical order = {slope_h:.3f}")
    print(f"L1(hu) empirical order = {slope_hu:.3f}")


if __name__ == "__main__":
    main()
