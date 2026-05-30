"""ANUGA reference run of the Stoker/Ritter dam-break for head-to-head
comparison against hydroflux (§3.8 of the methods paper).

Setup matches `solver-2d/examples/gen_verification_data.rs::gen_stoker`:
- 100 m × 5 m flat-bed channel (1D-equivalent: tangential walls).
- Initial dam at x = 50 m: h_L = 1 m, downstream dry.
- Walls on north/south; transmissive on east/west.
- t_end = 4 s.

ANUGA uses a triangular unstructured mesh — we choose
`maximum_triangle_area` so the centroid spacing in x is close to 1 m
(matched to a 100-cell hydroflux re-run for a fair head-to-head).

Output:
  papers/01_review/figures/data/anuga_stoker.csv
    columns: x, h_sim
  (centreline-equivalent: depth at the cell whose centroid has the
  smallest |y − 2.5| at each x).

Run (inside the anuga venv):
  /tmp/anuga_venv/bin/python solver-2d/examples/anuga_stoker_compare.py
"""

import csv
import math
import os

import numpy as np
import anuga

LENGTH_X = 100.0
LENGTH_Y = 5.0
H_L = 1.0
X_DAM = 50.0
T_END = 4.0
MAX_TRI_AREA = 0.25  # ~1 m effective dx (two triangles per square)

OUT_DIR = "papers/01_review/figures/data"
OUT = os.path.join(OUT_DIR, "anuga_stoker.csv")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Rectangular domain with structured-like triangulation.
    points, vertices, boundary = anuga.rectangular_cross(
        int(LENGTH_X / 1.0), int(LENGTH_Y / 1.0),
        LENGTH_X, LENGTH_Y
    )
    domain = anuga.Domain(points, vertices, boundary)
    domain.set_name("stoker_anuga")
    domain.set_quantity("elevation", 0.0)
    domain.set_quantity("friction", 0.0)

    # Initial condition: stage = h_L for x < X_DAM, 0 elsewhere.
    def stage_init(x, y):
        return np.where(x < X_DAM, H_L, 0.0)
    domain.set_quantity("stage", stage_init)

    # Boundaries: tangential walls on long sides, transmissive on the
    # ends (free outflow into the dry region).
    bcs = {
        "left":   anuga.Transmissive_boundary(domain),
        "right":  anuga.Transmissive_boundary(domain),
        "top":    anuga.Reflective_boundary(domain),
        "bottom": anuga.Reflective_boundary(domain),
    }
    domain.set_boundary(bcs)

    # March to t_end.
    for t in domain.evolve(yieldstep=T_END, finaltime=T_END):
        pass

    # Extract centreline (y ≈ LENGTH_Y / 2) depth profile.
    centroids = domain.get_centroid_coordinates(absolute=True)
    stage = domain.quantities["stage"].centroid_values
    elev = domain.quantities["elevation"].centroid_values
    h = stage - elev
    h = np.maximum(h, 0.0)

    # Bin centroids into 1 m x-bins, take the one closest to y=2.5
    # per bin (a centreline-equivalent sampling).
    bins = np.arange(0.0, LENGTH_X + 1e-9, 1.0)
    centres = 0.5 * (bins[:-1] + bins[1:])
    centreline = []
    for i in range(len(centres)):
        x_lo, x_hi = bins[i], bins[i + 1]
        mask = (centroids[:, 0] >= x_lo) & (centroids[:, 0] < x_hi)
        if not mask.any():
            continue
        dy = np.abs(centroids[mask, 1] - 0.5 * LENGTH_Y)
        j = np.argmin(dy)
        idxs = np.where(mask)[0]
        k = idxs[j]
        centreline.append((centres[i], h[k]))

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x", "h_sim"])
        w.writerows(centreline)
    print(f"wrote {OUT} ({len(centreline)} rows)")


if __name__ == "__main__":
    main()
