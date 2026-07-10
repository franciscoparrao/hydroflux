//! UK Environment Agency 2D benchmark suite — **Test 8A** (Néelz &
//! Pender 2013, EA report SC120002 §4.9: "Rainfall and point source
//! surface flow in urban areas" — Cockenzie Street, Glasgow), run
//! against the *official* input geometry (Sharifian et al. 2023,
//! LISFLOOD-FP 8.1 reproducibility package, Zenodo 10.5281/
//! zenodo.6907286, `4-Glasgow.zip/Setup/`).
//!
//! Closes review WP3 (`papers/01_review/ROADMAP_REVISION_EMS.md`) for
//! the third of three UK EA tests reproduced on official geometry.
//!
//! # Reference data caveat
//!
//! Unlike Test 4, no LISFLOOD-FP numeric reference *time series* is
//! redistributed for Test 8A (the Zenodo record ships only the model
//! **inputs**, not per-run outputs). The comparison below is
//! therefore against the **qualitative bounds reported in the text**
//! of SC120002 §4.9.3 (agreement ranges observed across the ~15
//! industry packages that ran this test), not a point-by-point RMSE
//! as in Test 4. See `benchmarks/data/uk_ea/test8a_glasgow/README.md`
//! (if present) for provenance detail.
//!
//! # Setup (from `benchmarks/data/uk_ea/test8a_glasgow/`)
//!
//! - Domain: **962 m × 398 m** (481×199 cells @ 2 m), real DEM
//!   (`ea8-2m.dem.gz`, elevation 21.1-37.4 m, no NODATA — buildings
//!   and kerbs are encoded as raised bed elevation, not masked cells).
//! - Friction: spatially varying Manning `n` (`ea8-2m.n.gz`, two
//!   values in the source raster: 0.02 roads, 0.05 elsewhere).
//! - Forcing (two independent sources, `ea8-2m.rain` / `ea8-2m.bdy`):
//!   - Uniform rainfall, **400 mm/h for a 3-minute pulse** (`t` in
//!     `[1, 4]` min) — report §4.9.1.
//!   - A point inflow at `(264894, 664750)`, trapezoidal, **peaking
//!     at 5 m³/s at t ≈ 37-39 min** — report §4.9.1 ("~35 min after
//!     the rainfall event"). LISFLOOD's point-BC `QVAR` convention is
//!     per-unit-width like the line BCs (confirmed: raw `.bdy` peak
//!     value 2.5 × cellsize 2 m = 5 m³/s, matching the report).
//! - All four domain boundaries closed (`Boundaries2D::WALLS`).
//! - Nine control points from `ea8-2m.stage`, order = report's point
//!   numbering 1-9.
//! - `sim_time = 18000` s (5 h) — `ea8-2m.par`.
//!
//! ```text
//! cargo run --release -p hydroflux-solver-2d --example uk_ea_test8a_official
//! ```

#[path = "uk_ea_common/mod.rs"]
mod uk_ea_common;

use std::time::Instant;

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2D, Mesh2D, PointSource, Simulation, SimulationConfig,
    apply_point_sources, apply_rain, read_ascii_grid,
};
use ndarray::Array2;

const DEM_PATH: &str = "benchmarks/data/uk_ea/test8a_glasgow/ea8-2m.dem.gz";
const MANNING_PATH: &str = "benchmarks/data/uk_ea/test8a_glasgow/ea8-2m.n.gz";

const CFL: f64 = 0.4;
const SIM_TIME: f64 = 18_000.0; // ea8-2m.par: sim_time
const POINT_BC_X: f64 = 264_894.0;
const POINT_BC_Y: f64 = 664_750.0;

/// Rainfall hyetograph, `(time [min], intensity [mm/h])`, from
/// `ea8-2m.rain` (raw file is `value time` per line; `.rain`'s two
/// near-duplicate breakpoints at `t = 0.9999/4.0001` approximate a
/// clean rectangular pulse without a literal vertical jump).
const RAIN_MIN: [(f64, f64); 6] = [
    (0.0, 0.0),
    (0.9999, 0.0),
    (1.0, 400.0),
    (4.0, 400.0),
    (4.0001, 0.0),
    (300.0, 0.0),
];

/// Point-source hydrograph, `(time [min], QVAR [m²/s])`, from
/// `ea8-2m.bdy`. Multiply by `cellsize` for `Q [m³/s]` (point-BC
/// `QVAR` uses the same per-unit-width convention as a line BC of
/// length one cell — see module docs).
const POINT_MIN: [(f64, f64); 23] = [
    (0.0, 0.0),
    (20.0, 0.0),
    (23.0, 0.03125),
    (25.0, 0.125),
    (27.0, 0.3125),
    (29.0, 0.65625),
    (31.0, 1.21875),
    (33.0, 1.78125),
    (35.0, 2.3125),
    (37.0, 2.5),
    (39.0, 2.5),
    (41.0, 2.3125),
    (43.0, 1.78125),
    (45.0, 1.21875),
    (47.0, 0.65625),
    (49.0, 0.3125),
    (51.0, 0.125),
    (53.0, 0.03125),
    (55.0, 0.0),
    (57.0, 0.0),
    (75.0, 0.0),
    (120.0, 0.0),
    (300.0, 0.0),
];

/// Official control points (`ea8-2m.stage`), report numbering 1-9.
const CONTROL_POINTS: [(f64, f64); 9] = [
    (264_680.0, 664_582.0),
    (264_536.0, 664_668.0),
    (264_354.0, 664_490.0),
    (264_200.0, 664_556.0),
    (264_332.0, 664_564.0),
    (264_572.0, 664_556.0),
    (264_708.0, 664_702.0),
    (264_306.0, 664_650.0),
    (264_220.0, 664_614.0),
];

fn interp_min(table: &[(f64, f64)], t_s: f64) -> f64 {
    let t_min = t_s / 60.0;
    let times: Vec<f64> = table.iter().map(|(t, _)| *t).collect();
    let vals: Vec<f64> = table.iter().map(|(_, v)| *v).collect();
    uk_ea_common::interp(&times, &vals, t_min)
}

fn rain_rate_m_s(t_s: f64) -> f64 {
    let mm_per_hour = interp_min(&RAIN_MIN, t_s);
    mm_per_hour * 1e-3 / 3600.0
}

fn point_discharge_m3_s(t_s: f64, cellsize: f64) -> f64 {
    interp_min(&POINT_MIN, t_s) * cellsize
}

fn main() {
    let t_start = Instant::now();

    let (bed, header) =
        read_ascii_grid(DEM_PATH).expect("failed to read Test 8A DEM — check the path");
    let (manning, manning_header) =
        read_ascii_grid(MANNING_PATH).expect("failed to read Test 8A Manning field");
    assert_eq!(
        header.ncols, manning_header.ncols,
        "DEM and Manning grids must share geometry (ncols)"
    );
    assert_eq!(
        header.nrows, manning_header.nrows,
        "DEM and Manning grids must share geometry (nrows)"
    );
    assert!(
        bed.iter().all(|&z| z.is_finite()),
        "Test 8A DEM has no NODATA cells (verified at acquisition) — found a non-finite elevation"
    );

    let mesh = Mesh2D::with_manning_field(bed, header.cellsize, header.cellsize, manning);
    println!(
        "Domain: {:.0} m x {:.0} m ({} x {} cells at {} m)",
        header.ncols as f64 * header.cellsize,
        header.nrows as f64 * header.cellsize,
        mesh.n_rows(),
        mesh.n_cols(),
        header.cellsize,
    );

    let point_cell = header
        .cell_at(POINT_BC_X, POINT_BC_Y)
        .expect("point BC falls outside the grid");
    println!("Point BC at ({POINT_BC_X}, {POINT_BC_Y}) -> cell {point_cell:?}");

    let control_cells: Vec<(usize, usize)> = CONTROL_POINTS
        .iter()
        .map(|&(x, y)| {
            header
                .cell_at(x, y)
                .unwrap_or_else(|| panic!("control point ({x},{y}) falls outside the grid"))
        })
        .collect();

    let states = Array2::from_elem((mesh.n_rows(), mesh.n_cols()), Conserved2D::DRY);

    // Fully dry domain, closed boundaries, rainfall does not start
    // until t = 60 s (1 min): the CFL bound sees no signal at all
    // before that and would return dt = INFINITY. Same mechanism as
    // Test 4 (`max_dt`); 10 s comfortably resolves the ~1 min rain
    // onset and the point source's 2-min ramp segments.
    let config = SimulationConfig {
        cfl: CFL,
        boundaries: Boundaries2D::WALLS,
        max_dt: 10.0,
        ..Default::default()
    };
    let mut sim = Simulation::new(mesh.clone(), states, config).expect("valid setup");

    let mut sim_times = vec![vec![0.0_f64]; CONTROL_POINTS.len()];
    let mut sim_depths: Vec<Vec<f64>> = control_cells
        .iter()
        .map(|&(i, j)| vec![sim.states()[(i, j)].h])
        .collect();

    let mut steps = 0usize;
    while sim.time() < SIM_TIME {
        let dt = sim
            .step(SIM_TIME - sim.time())
            .expect("simulation step failed");

        let t_new = sim.time();
        apply_rain(sim.states_mut(), rain_rate_m_s(t_new), dt);

        let q = point_discharge_m3_s(t_new, mesh.dx);
        let sources = [PointSource {
            row: point_cell.0,
            col: point_cell.1,
            q_mass: q,
        }];
        apply_point_sources(sim.states_mut(), &sources, dt, mesh.dx, mesh.dy);
        steps += 1;

        for (p, &(i, j)) in control_cells.iter().enumerate() {
            sim_times[p].push(t_new);
            sim_depths[p].push(sim.states()[(i, j)].h);
        }

        if steps % 5000 == 0 {
            let h_max = sim.states().iter().map(|s| s.h).fold(0.0_f64, f64::max);
            println!(
                "  t = {:>8.3} s ({:>5.1} %), dt = {dt:.5} s, h_max = {h_max:.4} m, steps = {steps}, wall = {:.1} s",
                sim.time(),
                100.0 * sim.time() / SIM_TIME,
                t_start.elapsed().as_secs_f64(),
            );
        }
    }
    println!(
        "\nDone: {steps} steps, t = {:.1} s, wall time {:.1} s",
        sim.time(),
        t_start.elapsed().as_secs_f64()
    );

    // --- Comparison against SC120002 §4.9.3's qualitative bounds ---
    println!("\n=== Peak / final depths at the 9 control points ===");
    println!("  {:>3}  {:>10}  {:>10}  {:>10}", "pt", "x", "y", "peak_h");
    let mut peaks = vec![0.0_f64; CONTROL_POINTS.len()];
    for (p, series) in sim_depths.iter().enumerate() {
        peaks[p] = series.iter().cloned().fold(0.0_f64, f64::max);
        let (x, y) = CONTROL_POINTS[p];
        println!("  {:>3}  {:>10.1}  {:>10.1}  {:>10.4}", p + 1, x, y, peaks[p]);
    }
    let finals: Vec<f64> = sim_depths.iter().map(|s| *s.last().unwrap()).collect();

    println!("\n=== Qualitative checks (report §4.9.3, no numeric reference series available) ===");
    println!(
        "  Point 1 peak depth > 0.5 m (report: models agree within ~5%): {:.4} m -> {}",
        peaks[0],
        if peaks[0] > 0.5 { "PASS" } else { "FAIL" }
    );
    for p in [1usize, 3, 6] {
        // points 2, 4, 7 (0-indexed 1, 3, 6)
        println!(
            "  Point {} peak depth <= ~0.35 m (report: models agree within ~0.04 m): {:.4} m -> {}",
            p + 1,
            peaks[p],
            if peaks[p] <= 0.40 { "PASS (within margin)" } else { "FAIL" }
        );
    }
    println!(
        "  Point 3 (downstream pond) final depth ~0.8 m (report: models agree within ~0.07 m): {:.4} m",
        finals[2]
    );
}
