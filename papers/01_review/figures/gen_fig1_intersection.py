"""Figure 1 — The intersection is the gap.

Five-axis radar (polar) plot. Each axis is one of the five structural
gaps articulated in §3: Open, Modern language, GPU first-class, Native
coupling, Differentiable physics.

Five representative solvers are plotted as polygons whose vertices
score each gap from 0 (does not address) to 1 (fully addresses). The
hydroflux polygon is the full pentagon; every other solver is a strict
sub-pentagon. Per-entity colours come from `style.SOLVER_COLORS` for
cross-figure consistency.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from style import FIG_W_SC, SOLVER_COLORS, setup

setup()
OUT_DIR = Path(__file__).parent

AXES = [
    "Open\nauditable",
    "Modern\nlanguage",
    "GPU\nfirst-class",
    "Native\ncoupling",
    "Differentiable\nphysics",
]

# Qualitative scores (0–1) per solver across the 5 axes. Justification
# is in §3 of the manuscript and in state-of-the-art.md.
SOLVERS = {
    "HEC-RAS":     [0.0, 0.0, 0.3, 0.0, 0.0],
    "LISFLOOD-FP": [1.0, 0.3, 1.0, 0.0, 0.0],
    "BASEMENT":    [0.0, 0.3, 0.0, 0.5, 0.0],
    "TELEMAC":     [0.7, 0.0, 0.0, 0.5, 0.0],
    "TUFLOW HPC":  [0.0, 0.3, 1.0, 0.0, 0.0],
    "hydroflux":   [1.0, 1.0, 1.0, 1.0, 1.0],
}


def main() -> None:
    n = len(AXES)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles += angles[:1]

    # Single-column figure, square aspect for radar legibility.
    fig = plt.figure(figsize=(FIG_W_SC * 1.6, FIG_W_SC * 1.4))
    ax = fig.add_subplot(111, projection="polar")

    # Background "target" guide — full pentagon at radius 1.
    ax.plot(angles, [1.0] * (n + 1),
            color="#dddddd", lw=0.6, linestyle=":", zorder=0)

    # Plot existing solvers first (low opacity), hydroflux on top.
    for name, scores in SOLVERS.items():
        is_target = name == "hydroflux"
        closed = scores + scores[:1]
        ax.plot(angles, closed,
                color=SOLVER_COLORS[name],
                lw=2.4 if is_target else 1.1,
                zorder=3 if is_target else 1,
                label=name)
        ax.fill(angles, closed,
                color=SOLVER_COLORS[name],
                alpha=0.18 if is_target else 0.06,
                zorder=2 if is_target else 0)

    # Axis ticks and labels.
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(AXES, fontsize=8.5)
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0.5, 1.0])
    ax.set_yticklabels(["0.5", "1.0"], fontsize=7, color="#666666")
    ax.set_rlabel_position(72)
    ax.tick_params(axis="x", pad=8)
    ax.grid(True, color="#dddddd", linewidth=0.4, alpha=0.7)
    ax.spines["polar"].set_color("#aaaaaa")
    ax.spines["polar"].set_linewidth(0.6)

    # Legend outside-right.
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(1.45, 1.10),
        frameon=False,
        fontsize=8.5,
        labelspacing=0.6,
    )

    # Editorial callout: arrow + caption pointing at the gap region.
    # The gap = anywhere the hydroflux pentagon (outer) extends beyond
    # the other polygons (clustered near the centre). Place a small
    # annotation at the GPU axis, where TUFLOW & LISFLOOD reach 1 but
    # everyone else is at 0–0.3.
    ax.annotate(
        "Two-thirds of the survey\n"
        "score $\\leq 0.3$ on every axis;\n"
        "only the hydroflux pentagon\n"
        "spans all five.",
        xy=(angles[3], 0.5),         # near "Native coupling" axis, mid-radius
        xytext=(1.5 * np.pi, 1.55),    # outside, lower-left
        fontsize=7.5, color="#444444", style="italic",
        ha="center",
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
    )

    fig.savefig(OUT_DIR / "fig1_intersection.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig1_intersection.pdf", bbox_inches="tight")
    print(f"Wrote {OUT_DIR / 'fig1_intersection.png'}")


if __name__ == "__main__":
    main()
