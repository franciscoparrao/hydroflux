"""Figure 2 — Stoker wet-wet dam break: profile + convergence."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from style import (
    COLOR_ANALYTICAL,
    COLOR_NUMERICAL,
    COLOR_REFERENCE,
    FIG_W_DC,
    add_panel_label,
    setup,
)

setup()
OUT_DIR = Path(__file__).parent

G = 9.81
H_L = 1.0
H_R = 0.1
X_DAM = 0.5
L_DOMAIN = 1.0
T_END = 0.075

# Convergence table from benchmarks/dam-break-results.md
L1_TABLE = [
    # (n, dx, L1(h), L1(hu))
    (50, 0.02000, 0.022062, 0.050882),
    (100, 0.01000, 0.012883, 0.028212),
    (200, 0.00500, 0.007525, 0.017454),
    (400, 0.00250, 0.004217, 0.009419),
    (800, 0.00125, 0.002409, 0.005334),
]


def stoker_star(h_l, h_r):
    c_l = np.sqrt(G * h_l)
    c_r = np.sqrt(G * h_r)

    def f(h):
        return 2.0 * (np.sqrt(G * h) - c_l) + (h - h_r) * np.sqrt(
            G * (h + h_r) / (2.0 * h * h_r)
        )

    lo, hi = h_r, h_l
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    h_star = 0.5 * (lo + hi)
    c_star = np.sqrt(G * h_star)
    f_r = (h_star - h_r) * np.sqrt(G * (h_star + h_r) / (2.0 * h_star * h_r))
    u_star = f_r
    shock = c_r * np.sqrt(h_star * (h_star + h_r) / (2.0 * h_r * h_r))
    return h_star, u_star, shock, -c_l, u_star - c_star


def analytical_h(x, t, h_l, h_r):
    h_star, u_star, s_r, head_rare, tail_rare = stoker_star(h_l, h_r)
    c_l = np.sqrt(G * h_l)
    xi = (x - X_DAM) / t
    h = np.empty_like(x)
    for i, xi_i in enumerate(xi):
        if xi_i < head_rare:
            h[i] = h_l
        elif xi_i < tail_rare:
            c = (2.0 * c_l - xi_i) / 3.0
            h[i] = c * c / G
        elif xi_i < s_r:
            h[i] = h_star
        else:
            h[i] = h_r
    return h, h_star, s_r


def stoker_solver_run(n):
    """Faithful Python HLL+forward Euler reproduction of solver-1d for
    the Stoker setup. Inline so the figure is regenerable without cargo."""
    dx = L_DOMAIN / n
    cfl = 0.4
    h = np.where(np.arange(n) * dx + 0.5 * dx < X_DAM, H_L, H_R)
    hu = np.zeros(n)
    t = 0.0
    while t < T_END:
        c = np.sqrt(G * np.maximum(h, 0))
        u = np.where(h > 0, hu / h, 0)
        smax = float(np.max(np.abs(u) + c))
        dt = min(cfl * dx / smax, T_END - t)
        h_pad = np.concatenate([[h[0]], h, [h[-1]]])
        hu_pad = np.concatenate([[hu[0]], hu, [hu[-1]]])
        c_pad = np.sqrt(G * np.maximum(h_pad, 0))
        u_pad = np.where(h_pad > 0, hu_pad / h_pad, 0)
        sl = np.minimum(u_pad[:-1] - c_pad[:-1], u_pad[1:] - c_pad[1:])
        sr = np.maximum(u_pad[:-1] + c_pad[:-1], u_pad[1:] + c_pad[1:])
        f_mass_l = hu_pad[:-1]
        f_mom_l = np.where(
            h_pad[:-1] > 0,
            hu_pad[:-1] ** 2 / np.maximum(h_pad[:-1], 1e-12) + 0.5 * G * h_pad[:-1] ** 2,
            0.0,
        )
        f_mass_r = hu_pad[1:]
        f_mom_r = np.where(
            h_pad[1:] > 0,
            hu_pad[1:] ** 2 / np.maximum(h_pad[1:], 1e-12) + 0.5 * G * h_pad[1:] ** 2,
            0.0,
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            denom = sr - sl
            f_mass = np.where(
                sl >= 0, f_mass_l,
                np.where(sr <= 0, f_mass_r,
                         (sr * f_mass_l - sl * f_mass_r
                          + sl * sr * (h_pad[1:] - h_pad[:-1])) / denom)
            )
            f_mom = np.where(
                sl >= 0, f_mom_l,
                np.where(sr <= 0, f_mom_r,
                         (sr * f_mom_l - sl * f_mom_r
                          + sl * sr * (hu_pad[1:] - hu_pad[:-1])) / denom)
            )
        h = h - (dt / dx) * (f_mass[1:] - f_mass[:-1])
        hu = hu - (dt / dx) * (f_mom[1:] - f_mom[:-1])
        t += dt
    x_centres = (np.arange(n) + 0.5) * dx
    return x_centres, h


def main() -> None:
    fig, (ax_l, ax_r) = plt.subplots(
        1, 2, figsize=(FIG_W_DC, FIG_W_DC * 0.42),
    )

    # ── Panel (a): depth profile at n = 400.
    n_show = 400
    x_num, h_num = stoker_solver_run(n_show)
    x_fine = np.linspace(0.0, L_DOMAIN, 2000)
    h_exact, h_star, s_r = analytical_h(x_fine, T_END, H_L, H_R)

    ax_l.plot(x_fine, h_exact,
              color=COLOR_ANALYTICAL, lw=1.6,
              label="Analytical (Stoker)")
    ax_l.plot(x_num, h_num,
              color=COLOR_NUMERICAL, lw=0,
              marker="o", markersize=2.4, mew=0, alpha=0.75,
              label=f"Numerical, $n = {n_show}$")

    # Editorial callout: arrow + caption highlighting the shock smear.
    shock_x = X_DAM + s_r * T_END
    ax_l.annotate(
        "Shock smeared\nover 3–5 cells\n(HLL signature)",
        xy=(shock_x, h_star * 0.6),
        xytext=(shock_x + 0.18, h_star + 0.15),
        fontsize=7.5, color="#444444", style="italic", ha="left",
        arrowprops=dict(arrowstyle="-|>", color="#888888",
                        lw=0.6, mutation_scale=8,
                        connectionstyle="arc3,rad=-0.15"),
    )

    ax_l.set_xlabel("$x$ (m)")
    ax_l.set_ylabel("Depth $h$ (m)")
    ax_l.set_xlim(0, 1)
    ax_l.set_ylim(0, 1.18)
    ax_l.legend(loc="upper right", fontsize=8)
    add_panel_label(ax_l, "(a)")

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

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig2_stoker.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig2_stoker.pdf", bbox_inches="tight")
    print(f"Wrote {OUT_DIR / 'fig2_stoker.png'}")
    print(f"L1(h) order = {slope_h:.3f}, L1(hu) order = {slope_hu:.3f}")


if __name__ == "__main__":
    main()
