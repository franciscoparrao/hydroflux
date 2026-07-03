//! UK Environment Agency 2D benchmark suite — **Test 4** (Néelz &
//! Pender 2013, EA report SC120002 §4.5: "Speed of flood propagation
//! over an extended floodplain"), run against the *official* input
//! geometry and compared quantitatively against LISFLOOD-FP reference
//! results at the six official control points.
//!
//! Closes review WP3 (`papers/01_review/ROADMAP_REVISION_EMS.md`):
//! the manuscript's §3.6 currently reports only a qualitative "passes
//! all six" on synthetic stand-ins (`solver-2d/tests/uk_ea_test*.rs`).
//! This is the first test run on the actual EA/LISFLOOD-FP geometry
//! with a numeric, citable reference.
//!
//! # Setup (from `benchmarks/data/uk_ea/test4/`, provenance in its
//! `README.md`)
//!
//! - Domain: **1000 m × 2000 m** flat floodplain (`z ≡ 0`; verified
//!   from the shipped DEM), Manning `n = 0.05` uniform (`ea4.par`
//!   `fpfric`).
//! - Inflow: a single opening **20 m wide, centred on the west edge**
//!   (`y ∈ [990, 1010]`, from `ea4.bci`'s `W 990.0 1010.0 QVAR test4`),
//!   trapezoidal hydrograph **peaking at 20 m³/s** (report §4.5.1),
//!   ramping 0→peak over `t ∈ [300, 3600]` s, holding to `t = 14400` s,
//!   then 0 by `t = 18000` s (5 h total — breakpoints from `ea4.bdy`,
//!   scaled from LISFLOOD's unit-width `QVAR` convention:
//!   `bdy_value [m²/s] × 20 m segment = Q [m³/s]`, confirmed against
//!   the report's stated 20 m³/s peak). East/South/North: `FREE`
//!   (`Boundary::Transmissive`).
//! - Six control points from `ea4.stage`: (50,1000), (100,1000),
//!   (200,1000), (300,1000), (400,1000), (300,1300).
//!
//! Run at 5 m resolution (`ea4-5m.dem.gz`) to match the resolution of
//! the DG2 reference (LISFLOOD-FP's high-order full-SWE scheme — the
//! closest conceptual match to hydroflux's HLLC+MUSCL+SSP-RK2); the
//! 1 m/ACC reference is also compared for context, at a different
//! resolution AND a different (inertial-storage-cell) scheme, so its
//! numbers are not a same-footing comparison.
//!
//! ```text
//! cargo run --release -p hydroflux-solver-2d --example uk_ea_test4_official
//! ```

#[path = "uk_ea_common/mod.rs"]
mod uk_ea_common;

use std::path::Path;
use std::time::Instant;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, PointSource, Simulation, SimulationConfig,
    apply_point_sources, read_ascii_grid,
};
use ndarray::Array2;

const DEM_PATH: &str = "benchmarks/data/uk_ea/test4/ea4-5m.dem.gz";
const REFERENCE_DG2_5M: &str = "benchmarks/data/uk_ea/test4/reference/ea4-5m-dg2.stage";
const REFERENCE_ACC_1M: &str = "benchmarks/data/uk_ea/test4/reference/ea4-1m-acc.stage";

const MANNING_N: f64 = 0.05; // ea4.par: fpfric
const CFL: f64 = 0.4;
const SIM_TIME: f64 = 18_000.0; // ea4.par: sim_time
const INLET_Y_MIN: f64 = 990.0;
const INLET_Y_MAX: f64 = 1010.0;
const PEAK_Q_M3S: f64 = 20.0; // report §4.5.1

/// Trapezoidal hydrograph breakpoints, `(time [s], fraction of peak)`,
/// from `ea4.bdy` (`0 0 / 0 300 / 1 3600 / 1 14400 / 0 18000`).
const HYDROGRAPH: [(f64, f64); 5] = [
    (0.0, 0.0),
    (300.0, 0.0),
    (3600.0, 1.0),
    (14400.0, 1.0),
    (18000.0, 0.0),
];

/// Official control points (`ea4.stage`): `(x, y)` in metres.
const CONTROL_POINTS: [(f64, f64); 6] = [
    (50.0, 1000.0),
    (100.0, 1000.0),
    (200.0, 1000.0),
    (300.0, 1000.0),
    (400.0, 1000.0),
    (300.0, 1300.0),
];

fn discharge_at(t: f64) -> f64 {
    let frac = uk_ea_common::interp(
        &HYDROGRAPH.iter().map(|(t, _)| *t).collect::<Vec<_>>(),
        &HYDROGRAPH.iter().map(|(_, q)| *q).collect::<Vec<_>>(),
        t,
    );
    frac * PEAK_Q_M3S
}

fn main() {
    let t_start = Instant::now();

    let (bed, header) =
        read_ascii_grid(DEM_PATH).expect("failed to read Test 4 DEM (5 m) — check the path");
    assert!(
        bed.iter().all(|&z| z == 0.0),
        "Test 4 bed must be perfectly flat (report §4.5.1); found non-zero elevation"
    );
    let mesh = Mesh2D::new(bed, header.cellsize, header.cellsize, MANNING_N);
    println!(
        "Domain: {:.0} m x {:.0} m ({} x {} cells at {} m), Manning n = {}",
        header.ncols as f64 * header.cellsize,
        header.nrows as f64 * header.cellsize,
        mesh.n_rows(),
        mesh.n_cols(),
        header.cellsize,
        MANNING_N
    );

    // Inlet: rows overlapping the 20 m segment, at the west edge (col 0).
    let inlet_rows = header.rows_overlapping_y_range(INLET_Y_MIN, INLET_Y_MAX);
    assert!(!inlet_rows.is_empty(), "inlet y-range fell outside the grid");
    println!(
        "Inlet: y in [{INLET_Y_MIN}, {INLET_Y_MAX}] m -> {} cell(s) at col 0, rows {:?}",
        inlet_rows.len(),
        inlet_rows
    );

    let bcs = Boundaries2D {
        north: Boundary::Transmissive,
        south: Boundary::Transmissive,
        east: Boundary::Transmissive,
        west: Boundary::Wall, // background west edge is closed; the inlet is a point-source override
    };

    let control_cells: Vec<(usize, usize)> = CONTROL_POINTS
        .iter()
        .map(|&(x, y)| {
            header
                .cell_at(x, y)
                .unwrap_or_else(|| panic!("control point ({x},{y}) falls outside the grid"))
        })
        .collect();

    let states = Array2::from_elem((mesh.n_rows(), mesh.n_cols()), Conserved2D::DRY);

    // The domain starts fully dry with a closed background boundary
    // (Wall) and no wet ghost anywhere — the CFL bound sees no signal
    // and would return dt = INFINITY, swallowing the whole 5 h event
    // (and the entire hydrograph) in a single step, since inflow only
    // enters via the point-source override below, which the CFL
    // calculation cannot see. `max_dt` caps this: once water enters,
    // the physical CFL bound governs as usual (max_dt only bites
    // during the dry/near-dry cold start). 15 s comfortably resolves
    // the hydrograph's 300 s ramp-up.
    let config = SimulationConfig {
        cfl: CFL,
        boundaries: bcs,
        max_dt: 15.0,
        ..Default::default()
    };
    let mut sim = Simulation::new(mesh.clone(), states, config).expect("valid setup");

    // Per-point (time, depth) series, sampled every step.
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

        let q_total = discharge_at(sim.time());
        let q_per_cell = q_total / inlet_rows.len() as f64;
        let sources: Vec<PointSource> = inlet_rows
            .iter()
            .map(|&row| PointSource {
                row,
                col: 0,
                q_mass: q_per_cell,
            })
            .collect();
        apply_point_sources(sim.states_mut(), &sources, dt, mesh.dx, mesh.dy);
        steps += 1;

        for (p, &(i, j)) in control_cells.iter().enumerate() {
            sim_times[p].push(sim.time());
            sim_depths[p].push(sim.states()[(i, j)].h);
        }

        if steps % 2000 == 0 {
            let h_max = sim.states().iter().map(|s| s.h).fold(0.0_f64, f64::max);
            println!(
                "  t = {:>8.1} s ({:>5.1} %), h_max = {h_max:.4} m, steps = {steps}",
                sim.time(),
                100.0 * sim.time() / SIM_TIME
            );
        }
    }
    println!(
        "\nDone: {steps} steps, t = {:.1} s, wall time {:.1} s",
        sim.time(),
        t_start.elapsed().as_secs_f64()
    );

    // --- Comparison against LISFLOOD-FP references -----------------
    for (label, path, note) in [
        (
            "DG2 @ 5 m (resolution-matched, full 2nd-order SWE)",
            REFERENCE_DG2_5M,
            "",
        ),
        (
            "ACC @ 1 m (different resolution AND scheme — context only)",
            REFERENCE_ACC_1M,
            " [not a same-footing comparison]",
        ),
    ] {
        println!("\n=== vs LISFLOOD-FP {label}{note} ===");
        let reference = uk_ea_common::parse_reference_stage(Path::new(path));
        println!(
            "  {:>3}  {:>8}  {:>8}  {:>10}  {:>10}  {:>10}  {:>10}",
            "pt", "x", "y", "ref_peak", "rmse", "peak_bias", "arrival_dt"
        );
        for (p, &(x, y)) in CONTROL_POINTS.iter().enumerate() {
            let cmp = uk_ea_common::compare_point(
                &reference,
                p,
                &sim_times[p],
                &sim_depths[p],
                0.02, // 2 cm arrival threshold
            );
            let arrival_dt = match (cmp.sim_arrival, cmp.reference_arrival) {
                (Some(s), Some(r)) => format!("{:+.1} s", s - r),
                (None, Some(_)) => "never (sim)".to_string(),
                (Some(_), None) => "never (ref)".to_string(),
                (None, None) => "never (both)".to_string(),
            };
            println!(
                "  {:>3}  {:>8.1}  {:>8.1}  {:>10.4}  {:>10.4}  {:>+10.4}  {:>10}",
                p + 1,
                x,
                y,
                cmp.reference_peak,
                cmp.rmse,
                cmp.peak_bias,
                arrival_dt
            );
        }
    }
}
