"""Figure 1 — The intersection is the gap.

Five-axis radar (a.k.a. polar) plot. Each axis is one of the five gaps
articulated in §3 of the paper:

  1. Open auditable codebase
  2. Modern host language with ergonomic autograd
  3. GPU as first-class execution target
  4. Native coupling to non-hydraulic hazards
  5. Differentiable physics across the engine

Five representative solvers are plotted as polygons whose vertices
score each gap from 0 (does not address) to 1 (fully addresses). The
hydroflux polygon is the full pentagon; every other solver is a strict
sub-pentagon. The visual story: the intersection is exactly the
hydroflux outline, and no existing solver reaches it.

Scoring rationale documented inline. The scores are deliberately
qualitative — the goal is to render the *shape* of the gap, not to
rank solvers numerically.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path(__file__).parent

# Axes in order: Open, Modern lang, GPU, Coupled, Differentiable.
AXES = [
    "Open\nauditable",
    "Modern\nlanguage",
    "GPU\nfirst-class",
    "Native\ncoupling",
    "Differentiable\nphysics",
]

# Scoring rationale per solver across the 5 axes.
# 1.0 = fully addresses; 0.5 = partial; 0.0 = absent.
# See REVIEW_CHECKLIST.md and state-of-the-art.md for justification.
SOLVERS = {
    "HEC-RAS": {
        "scores": [0.0, 0.0, 0.3, 0.0, 0.0],
        "color": "#a04040",  # regulatory red
    },
    "LISFLOOD-FP": {
        "scores": [1.0, 0.3, 1.0, 0.0, 0.0],
        "color": "#5c8aab",  # open-source blue
    },
    "BASEMENT": {
        "scores": [0.0, 0.3, 0.0, 0.5, 0.0],
        "color": "#8b6a1a",  # closed-academic brown
    },
    "TELEMAC": {
        "scores": [0.7, 0.0, 0.0, 0.5, 0.0],
        "color": "#5c4630",  # legacy brown
    },
    "TUFLOW HPC": {
        "scores": [0.0, 0.3, 1.0, 0.0, 0.0],
        "color": "#5c5c5c",  # commercial grey
    },
    "hydroflux": {
        "scores": [1.0, 1.0, 1.0, 1.0, 1.0],
        "color": "#1a4a72",  # navy = the target
    },
}


def main() -> None:
    n = len(AXES)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles += angles[:1]  # close the polygon

    fig, ax = plt.subplots(
        figsize=(7.5, 7.5),
        subplot_kw={"projection": "polar"},
    )

    # Plot order: existing solvers first (lower alpha), hydroflux last
    # (full opacity) so it sits on top.
    for name, props in SOLVERS.items():
        scores = props["scores"] + props["scores"][:1]
        is_target = name == "hydroflux"
        ax.plot(angles, scores,
                color=props["color"],
                lw=2.6 if is_target else 1.4,
                linestyle="-" if is_target else "-",
                label=name)
        ax.fill(angles, scores,
                color=props["color"],
                alpha=0.20 if is_target else 0.08)

    # Axis ticks and labels.
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(AXES, fontsize=10)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["", "0.5", "", "1.0"], fontsize=8)
    ax.set_rlabel_position(72)

    # Legend outside the plot.
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(1.30, 1.05),
        frameon=False,
        fontsize=10,
    )

    # Subtle annotation explaining the visual claim.
    fig.text(
        0.50, 0.02,
        "Each axis: 0 = does not address, 1 = fully addresses. "
        "The hydroflux pentagon is the intersection of all five gaps.",
        ha="center", fontsize=9, style="italic", color="#444444",
    )

    fig.savefig(OUT_DIR / "fig1_intersection.png", dpi=220,
                bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig1_intersection.pdf", bbox_inches="tight")
    print(f"Wrote {OUT_DIR / 'fig1_intersection.png'}")


if __name__ == "__main__":
    main()
