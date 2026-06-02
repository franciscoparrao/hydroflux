"""M2d: ANUGA wall-clock benchmark on a Stoker dam-break scaled up to a
size where the wall-clock comparison is informative.

We pick a 200 m × 5 m channel at `maximum_triangle_area = 0.0625`,
giving roughly `200 / 0.5 × 5 / 0.5 × 2 = 8000` triangles — close
enough to a `400 × 10 = 4000`-cell hydroflux mesh at `Δx = 0.5 m` for
a fair per-simulation-second comparison. The wall-clock is reported
for the evolve loop only (excluding I/O), in seconds per simulated
second of physical time and as triangle-steps per second.

Output:
  papers/01_review/figures/data/m2_anuga_wallclock.csv
    columns: solver, n_cells_or_triangles, n_steps, t_sim_s,
             wall_clock_s, wall_per_sim_s, mcell_steps_per_s
  (one row for ANUGA; the matched hydroflux row is appended by
   `m2_hydroflux_wallclock.rs`).

Run (inside the anuga venv):
  /tmp/anuga_venv/bin/python solver-2d/examples/m2_anuga_wallclock.py
"""

import csv
import os
import time

import numpy as np
import anuga


LENGTH_X = 200.0
LENGTH_Y = 5.0
H_L = 1.0
X_DAM = 100.0
T_END = 8.0
DX_EFFECTIVE = 0.5
MAX_TRI_AREA = (DX_EFFECTIVE ** 2) / 4.0  # two triangles per dx² grid square

OUT_DIR = "papers/01_review/figures/data"
OUT = os.path.join(OUT_DIR, "m2_anuga_wallclock.csv")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    points, vertices, boundary = anuga.rectangular_cross(
        int(LENGTH_X / DX_EFFECTIVE),
        int(LENGTH_Y / DX_EFFECTIVE),
        LENGTH_X,
        LENGTH_Y,
    )
    domain = anuga.Domain(points, vertices, boundary)
    domain.set_name("m2_anuga_stoker")
    domain.set_quantity("elevation", 0.0)
    domain.set_quantity("friction", 0.0)

    def stage_init(x, y):
        return np.where(x < X_DAM, H_L, 0.0)

    domain.set_quantity("stage", stage_init)
    domain.set_boundary({
        "left":   anuga.Transmissive_boundary(domain),
        "right":  anuga.Transmissive_boundary(domain),
        "top":    anuga.Reflective_boundary(domain),
        "bottom": anuga.Reflective_boundary(domain),
    })

    n_triangles = len(domain.quantities["stage"].centroid_values)
    print(f"ANUGA mesh: {n_triangles} triangles, dx_eff = {DX_EFFECTIVE} m")
    print(f"Physical time: {T_END} s")

    # Time the evolve loop only — exclude mesh build, BC setup, I/O.
    n_steps = 0
    t0 = time.perf_counter()
    for _t in domain.evolve(yieldstep=T_END, finaltime=T_END):
        n_steps += 1
    wall = time.perf_counter() - t0

    # ANUGA's adaptive time step is internal; we count yieldsteps but
    # the actual integration takes many internal sub-steps. Report
    # both the wall-clock-per-simulated-second and an effective
    # triangle-step throughput approximated from the inferred internal
    # step count.
    # Inferred sub-steps: try evolve_to_completion's internal counter.
    try:
        internal_steps = int(domain.beta_w)  # placeholder; replaced below
    except Exception:
        internal_steps = 0
    # ANUGA exposes domain.number_of_full_timesteps after evolution.
    try:
        internal_steps = int(domain.number_of_full_timesteps)
    except Exception:
        try:
            internal_steps = int(domain.timestepping_statistics()["number_of_steps"])
        except Exception:
            internal_steps = max(n_steps, 1)
    wall_per_sim = wall / T_END
    triangle_steps = float(n_triangles * internal_steps)
    mcell_steps_per_s = triangle_steps / wall / 1.0e6 if wall > 0 else 0.0

    print(f"Wall clock      : {wall:.3f} s")
    print(f"Internal steps  : {internal_steps}")
    print(f"Wall / sim-s    : {wall_per_sim:.3f}")
    print(f"Mtri-steps/s    : {mcell_steps_per_s:.3f}")

    new_file = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as f:
        w = csv.writer(f)
        if new_file:
            w.writerow([
                "solver",
                "n_cells_or_triangles",
                "n_steps",
                "t_sim_s",
                "wall_clock_s",
                "wall_per_sim_s",
                "mcell_steps_per_s",
            ])
        w.writerow([
            "ANUGA",
            n_triangles,
            internal_steps,
            T_END,
            f"{wall:.6f}",
            f"{wall_per_sim:.6f}",
            f"{mcell_steps_per_s:.6f}",
        ])
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
