//! UK Environment Agency 2D benchmark suite — **Test 5** (Néelz &
//! Pender 2013, EA report SC120002 §4.6 / Appendix A.5: "Valley
//! flooding", dam-break wave propagation down a river valley), run
//! against a **synthetic reconstruction**, not the official geometry.
//!
//! # Why synthetic, and what that means for the result
//!
//! Unlike Test 4 and Test 8A, no official Test 5 input package
//! (`Test5DEM.asc`, `Test5BC.csv`, `Test5Output.csv`) is redistributed
//! anywhere public — it is proprietary EA data, requested by email.
//! An attempt to recover the true valley *footprint* from the
//! reference-output rasters' NODATA pattern (Zenodo 10.5281/
//! zenodo.4066824, `ea5.zip`) failed: those rasters are a plain
//! rectangle with zero NODATA cells — the valley shape lives only in
//! the bed elevation, which is not in the package. See
//! `benchmarks/data/uk_ea/README.md` for the full investigation.
//!
//! This example instead builds an **idealised straight valley** from
//! the report's *text* description (§4.6.1, Appendix A.5) — every
//! number below that is not a direct quote of the report is an
//! explicit, documented assumption:
//!
//! | Quantity | Value | Source |
//! |---|---|---|
//! | Valley length | 17,000 m | report: "~17 km" |
//! | Valley width | 800 m | report: "~0.8 km" |
//! | Manning n | 0.04 uniform | report §A5.3 |
//! | Grid resolution | 50 m | report §A5.3 (the "expected" resolution) |
//! | Upper-reach slope | 0.01 | report: "~0.01 in its upper region" |
//! | Lower-reach slope | 0.001 | report: "easing to ~0.001" |
//! | **Upper/lower transition** | **s = 5,000 m** | **ASSUMED — not in the report** |
//! | Inflow line length | 260 m | report: "~260m long line" |
//! | **Inflow hydrograph shape** | **see below** | **ASSUMED — `Test5BC.csv` not available; only "skewed trapezoidal, short early peak at 3000 m³/s" is stated (§4.6.1), no breakpoint table** |
//! | **Cross-section shape** | **flat 200 m thalweg + linear banks to +8 m** | **ASSUMED — the real valley is not straight/symmetric** |
//! | Control-point along-valley distance | 6 of 7 from report text (§4.6.4/Figs. 4.21-4.23) | point 5 (downstream pond) has no stated distance — **ASSUMED** near the domain's closed downstream end |
//!
//! **Consequence**: this is not a replica of the official test in the
//! sense Test 4/8A are. Treat the comparison against the redistributed
//! LISFLOOD-FP reference series (`benchmarks/data/uk_ea/test5/
//! reference-ea5-50m-{acc,dg2}.stage`) as an **order-of-magnitude /
//! qualitative sanity check** (does hydroflux produce peak levels and
//! arrival times in the right ballpark for a valley of this scale and
//! forcing), not a quantitative validation. A real quantitative Test 5
//! result requires the official `Test5DEM.asc` (request from
//! `fcerm.evidence@environment-agency.gov.uk`).
//!
//! ```text
//! cargo run --release -p hydroflux-solver-2d --example uk_ea_test5_synthetic_valley
//! ```

use std::path::Path;
use std::time::Instant;

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2D, Mesh2D, PointSource, Simulation, SimulationConfig,
    apply_point_sources,
};
use ndarray::Array2;

#[path = "uk_ea_common/mod.rs"]
mod uk_ea_common;

const CELLSIZE: f64 = 50.0;
const VALLEY_LENGTH_M: f64 = 17_000.0;
const VALLEY_WIDTH_M: f64 = 800.0;
const N_ROWS: usize = (VALLEY_LENGTH_M / CELLSIZE) as usize; // 340, row 0 = upstream/inflow
const N_COLS: usize = (VALLEY_WIDTH_M / CELLSIZE) as usize; // 16

const MANNING_N: f64 = 0.04;
const CFL: f64 = 0.4;
const SIM_TIME: f64 = 30.0 * 3600.0; // report: run to t = 30 h

const SLOPE_UPPER: f64 = 0.01;
const SLOPE_LOWER: f64 = 0.001;
const TRANSITION_S_M: f64 = 5_000.0; // ASSUMED — not stated in the report

const THALWEG_LO: usize = 6; // cols 6-9 (4 * 50 m = 200 m) are the flat channel bottom
const THALWEG_HI: usize = 9;
const BANK_HEIGHT_M: f64 = 8.0; // ASSUMED valley-wall relief

const INLET_LINE_M: f64 = 260.0; // report: inflow line length
const INLET_COL_LO: usize = 5; // ceil(260/50) = 6 cols, straddling the thalweg
const INLET_COL_HI: usize = 10;

const POINT_CENTER_COL: usize = 7; // all points placed on the thalweg centerline (no cross-valley position is known)

/// Along-valley distance from the inflow [m] for each control point,
/// report numbering 1-7. Points 1, 2, 3, 4, 6, 7 are stated in the
/// report text (§4.6.4 notes to Figs. 4.21-4.23); point 5 (the
/// downstream pond) has no stated distance and is placed, as an
/// ASSUMPTION, near the domain's closed downstream end.
const POINT_DISTANCE_M: [f64; 7] = [
    3_240.0,  // point 1
    5_290.0,  // point 2
    7_080.0,  // point 3
    10_460.0, // point 4
    16_000.0, // point 5 — ASSUMED (downstream pond, no report value)
    3_670.0,  // point 6
    7_330.0,  // point 7
];

/// Inflow hydrograph, `(time [s], Q [m³/s])` — ASSUMED shape: fast
/// rise to the report's stated peak (3000 m³/s, "short early peak"),
/// long recession tail over the 30 h window ("skewed trapezoidal").
/// `Test5BC.csv` (the official breakpoint table) is not available;
/// no numeric hydrograph table appears in the SC120002 report text.
const HYDROGRAPH: [(f64, f64); 9] = [
    (0.0, 0.0),
    (900.0, 0.0),
    (1_800.0, 3_000.0), // peak — the one number the report actually states
    (3_600.0, 2_400.0),
    (7_200.0, 1_200.0),
    (14_400.0, 500.0),
    (28_800.0, 150.0),
    (54_000.0, 50.0),
    (108_000.0, 0.0),
];

fn discharge_at(t: f64) -> f64 {
    let times: Vec<f64> = HYDROGRAPH.iter().map(|(t, _)| *t).collect();
    let vals: Vec<f64> = HYDROGRAPH.iter().map(|(_, q)| *q).collect();
    uk_ea_common::interp(&times, &vals, t)
}

fn z_cross(col: usize) -> f64 {
    if (THALWEG_LO..=THALWEG_HI).contains(&col) {
        0.0
    } else if col < THALWEG_LO {
        let d = (THALWEG_LO - col) as f64;
        (d / THALWEG_LO as f64) * BANK_HEIGHT_M
    } else {
        let d = (col - THALWEG_HI) as f64;
        let bank_cells = (N_COLS - 1 - THALWEG_HI) as f64;
        (d / bank_cells) * BANK_HEIGHT_M
    }
}

fn z_longitudinal(row: usize) -> f64 {
    let s = row as f64 * CELLSIZE;
    let z_transition = SLOPE_LOWER * (VALLEY_LENGTH_M - TRANSITION_S_M);
    if s <= TRANSITION_S_M {
        z_transition + SLOPE_UPPER * (TRANSITION_S_M - s)
    } else {
        SLOPE_LOWER * (VALLEY_LENGTH_M - s)
    }
}

/// Parse the redistributed LISFLOOD-FP Test 5 reference `.stage`
/// files: unlike Test 4's, these have **no header at all** (no
/// coordinate block) — just `time depth_1 .. depth_7` per line.
fn parse_headerless_series(path: &Path, n_points: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {path:?}: {e}"));
    let mut times = Vec::new();
    let mut depths = Vec::new();
    for line in text.lines() {
        let vals: Result<Vec<f64>, _> = line.split_whitespace().map(str::parse::<f64>).collect();
        let Ok(vals) = vals else { continue };
        if vals.len() != n_points + 1 {
            continue;
        }
        times.push(vals[0]);
        depths.push(vals[1..].to_vec());
    }
    assert!(!times.is_empty(), "no data rows parsed from {path:?}");
    (times, depths)
}

fn main() {
    let t_start = Instant::now();

    let inlet_cols = INLET_COL_HI - INLET_COL_LO + 1;
    assert_eq!(
        inlet_cols,
        (INLET_LINE_M / CELLSIZE).ceil() as usize,
        "inlet column span must match the report's ~260 m inflow line"
    );

    println!(
        "Synthetic Test 5 valley: {N_ROWS} x {N_COLS} cells at {CELLSIZE} m ({:.1} km x {:.1} km)",
        N_ROWS as f64 * CELLSIZE / 1000.0,
        N_COLS as f64 * CELLSIZE / 1000.0,
    );

    let mut bed = Array2::<f64>::zeros((N_ROWS, N_COLS));
    for i in 0..N_ROWS {
        for j in 0..N_COLS {
            bed[(i, j)] = z_longitudinal(i) + z_cross(j);
        }
    }
    println!(
        "Thalweg elevation: {:.1} m (inflow) -> {:.1} m (downstream end)",
        bed[(0, POINT_CENTER_COL)],
        bed[(N_ROWS - 1, POINT_CENTER_COL)]
    );

    let mesh = Mesh2D::new(bed, CELLSIZE, CELLSIZE, MANNING_N);
    let states = Array2::from_elem((mesh.n_rows(), mesh.n_cols()), Conserved2D::DRY);

    let config = SimulationConfig {
        cfl: CFL,
        boundaries: Boundaries2D::WALLS,
        max_dt: 30.0,
        ..Default::default()
    };
    let mut sim = Simulation::new(mesh.clone(), states, config).expect("valid setup");

    let control_rows: Vec<usize> = POINT_DISTANCE_M
        .iter()
        .map(|&d| ((d / CELLSIZE).round() as usize).min(N_ROWS - 1))
        .collect();

    let mut sim_times = vec![vec![0.0_f64]; POINT_DISTANCE_M.len()];
    let mut sim_depths: Vec<Vec<f64>> = control_rows
        .iter()
        .map(|&row| vec![sim.states()[(row, POINT_CENTER_COL)].h])
        .collect();

    let mut steps = 0usize;
    while sim.time() < SIM_TIME {
        let dt = sim
            .step(SIM_TIME - sim.time())
            .expect("simulation step failed");

        let q_total = discharge_at(sim.time());
        let q_per_cell = q_total / inlet_cols as f64;
        let sources: Vec<PointSource> = (INLET_COL_LO..=INLET_COL_HI)
            .map(|col| PointSource {
                row: 0,
                col,
                q_mass: q_per_cell,
            })
            .collect();
        apply_point_sources(sim.states_mut(), &sources, dt, mesh.dx, mesh.dy);
        steps += 1;

        for (p, &row) in control_rows.iter().enumerate() {
            sim_times[p].push(sim.time());
            sim_depths[p].push(sim.states()[(row, POINT_CENTER_COL)].h);
        }

        if steps % 2000 == 0 {
            let h_max = sim.states().iter().map(|s| s.h).fold(0.0_f64, f64::max);
            println!(
                "  t = {:>8.1} s ({:>5.1} %), dt = {dt:.4} s, h_max = {h_max:.4} m, steps = {steps}, wall = {:.1} s",
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

    println!("\n=== Peak depths at the 7 control points (thalweg centerline) ===");
    let mut peaks = vec![0.0_f64; POINT_DISTANCE_M.len()];
    for (p, series) in sim_depths.iter().enumerate() {
        peaks[p] = series.iter().cloned().fold(0.0_f64, f64::max);
        println!(
            "  pt {}  s={:>7.0} m  peak_h = {:.4} m  (final = {:.4} m)",
            p + 1,
            POINT_DISTANCE_M[p],
            peaks[p],
            series.last().unwrap()
        );
    }

    println!("\n=== Order-of-magnitude comparison vs LISFLOOD-FP references (50 m, official resolution) — see module docs: this is a SYNTHETIC geometry, not the official one ===");
    for (label, path) in [
        (
            "DG2 @ 50 m",
            "benchmarks/data/uk_ea/test5/reference-ea5-50m-dg2.stage",
        ),
        (
            "ACC @ 50 m",
            "benchmarks/data/uk_ea/test5/reference-ea5-50m-acc.stage",
        ),
    ] {
        println!("\n--- vs {label} ---");
        let (ref_times, ref_depths) = parse_headerless_series(Path::new(path), 7);
        println!(
            "  {:>3}  {:>10}  {:>10}  {:>10}",
            "pt", "ref_peak", "sim_peak", "peak_diff"
        );
        for p in 0..7 {
            let ref_series: Vec<f64> = ref_depths.iter().map(|row| row[p]).collect();
            let ref_peak = ref_series.iter().cloned().fold(0.0_f64, f64::max);
            let _ = &ref_times;
            println!(
                "  {:>3}  {:>10.4}  {:>10.4}  {:>+10.4}",
                p + 1,
                ref_peak,
                peaks[p],
                peaks[p] - ref_peak
            );
        }
    }
}
