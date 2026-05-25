//! Phase 2 iter 2: full Atacama 2017 event over 21 days at Santa
//! Juana, with time-varying daily Q from the DGA + Manning-normal-depth
//! warm-start on channel cells (reduces initial transient versus the
//! thin-film start of iter 1).
//!
//! # Setup vs iter 1
//!
//! Same geometry / BCs / Manning as `huasco_2d_steady`:
//! - Subset 200×67 portrait DEM (6 km N-S × 2 km E-W, gauge-centred).
//! - W = Transmissive (outflow), N/S/E = Wall.
//! - Single `PointSource` at the E-edge channel cell (135, 66).
//! - Manning `n = 0.04` (gravel-bed Andean, Chow 1959).
//!
//! Differences:
//! - **Time-varying Q**: 21 daily values from DGA station 03820003,
//!   window 2017-02-20 → 2017-03-12 (Aluvión Atacama 2017). PointSource
//!   `q_mass` is updated once per day; intra-day Q is held constant.
//! - **Manning warm-start**: channel cells (acc > 1e6) start with
//!   `h_initial = h_n(Q_day1, slope_mean) = ((n · Q_day1 / W_eff) / √S₀)^(3/5)`
//!   instead of thin film `h = 1 cm`. The mean reach slope `S₀ = 0.0074`
//!   comes from the 1D longitudinal-profile extraction.
//! - **Snapshot writer**: one depth GeoTIFF per simulated day, named
//!   `huasco_2d_depth_day_NN.tif`, plus a mass-balance log row.
//!
//! Two run modes:
//! - `--days N` (default 21): how many days of the event to simulate.
//!   Use 5 for a quick rising-limb-to-peak smoke test (~1 h wall time).
//! - Default `--days 21` runs the full event (~4-5 h wall time, plan
//!   for overnight or background).
//!
//! Reproducir (full event):
//! ```text
//! cargo run --release -p hydroflux-solver-2d --example huasco_2d_event
//! ```
//! Quick test (5 days):
//! ```text
//! cargo run --release -p hydroflux-solver-2d --example huasco_2d_event -- --days 5
//! ```

use std::env;
use std::path::PathBuf;
use std::time::Instant;

use ndarray::Array2;
use surtgis_core::io::read_geotiff;
use surtgis_core::raster::Raster;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, PointSource, apply_point_sources, cfl_time_step_with_bcs,
    manning_friction_step, mesh_from_geotiff, ssprk2_step, write_depth_geotiff,
};

const SUBSET_DEM: &str = "examples/huasco_2d_phase2/output/huasco_subset_dem.tif";
const SUBSET_ACC: &str = "examples/huasco_2d_phase2/output/huasco_subset_acc.tif";
const OUTPUT_DIR: &str = "examples/huasco_2d_phase2/output";

const MANNING_N: f64 = 0.04;
const ACC_THRESHOLD: f64 = 1_000_000.0;
const SLOPE_MEAN: f64 = 0.0074; // from 1D longitudinal-profile extraction
const CFL: f64 = 0.4;
const SECONDS_PER_DAY: f64 = 86_400.0;
const INFLOW_ROW: usize = 135;
const INFLOW_COL: usize = 66;

/// Atacama 2017 daily Q at Santa Juana (DGA station 03820003),
/// window 2017-02-20 → 2017-03-12. Same series as the 1D paper
/// (`autograd/examples/calibrate_manning_huasco_2017.rs`).
const Q_DAILY_M3S: [f64; 21] = [
    17.5, 18.7, 18.4, 18.5, 20.5, 31.9, 34.8, 35.5, 37.8, 38.8, 38.9, 38.1, 37.5, 37.5, 36.0, 36.0,
    35.2, 34.8, 34.9, 33.9, 33.6,
];

fn manning_normal_depth(q_m3s: f64, n: f64, slope: f64, cell_width_m: f64) -> f64 {
    let q_per_w = q_m3s / cell_width_m;
    (n * q_per_w / slope.sqrt()).powf(3.0 / 5.0)
}

fn parse_days_arg() -> usize {
    let mut iter = env::args().skip(1);
    while let Some(arg) = iter.next() {
        if arg == "--days" {
            if let Some(v) = iter.next() {
                return v.parse().unwrap_or(Q_DAILY_M3S.len());
            }
        }
    }
    Q_DAILY_M3S.len()
}

fn main() {
    let n_days = parse_days_arg().min(Q_DAILY_M3S.len());
    let t_start = Instant::now();

    let (mesh, transform) = mesh_from_geotiff(SUBSET_DEM, MANNING_N)
        .expect("failed to load DEM subset; run extract_subset.py first");
    println!(
        "DEM: {}×{} cells, dx={} m, dy={} m, manning={}",
        mesh.n_rows(),
        mesh.n_cols(),
        mesh.dx,
        mesh.dy,
        mesh.manning
    );

    let acc: Raster<f64> =
        read_geotiff(SUBSET_ACC, None).expect("failed to load flow_accumulation subset");
    let acc_data = acc.data();
    assert_eq!(acc_data.dim(), mesh.bed.dim());
    let n_channel: usize = acc_data.iter().filter(|&&v| v > ACC_THRESHOLD).count();
    println!(
        "Channel cells (acc > {:.0e}): {} ({:.1}% of domain)",
        ACC_THRESHOLD,
        n_channel,
        100.0 * n_channel as f64 / mesh.n_cells() as f64
    );

    // Manning warm-start: pre-fill channel cells with h_n at day-1 Q.
    let h_warm = manning_normal_depth(Q_DAILY_M3S[0], MANNING_N, SLOPE_MEAN, mesh.dx);
    println!(
        "\nManning warm-start: Q_day1={} m³/s, S₀={}, dx={} m → h_n={:.3} m",
        Q_DAILY_M3S[0], SLOPE_MEAN, mesh.dx, h_warm
    );
    let mut states = Array2::from_shape_fn((mesh.n_rows(), mesh.n_cols()), |(i, j)| {
        if acc_data[(i, j)] > ACC_THRESHOLD {
            Conserved2D::new(h_warm, 0.0, 0.0)
        } else {
            Conserved2D::DRY
        }
    });

    let bcs = Boundaries2D {
        north: Boundary::Wall,
        south: Boundary::Wall,
        east: Boundary::Wall,
        west: Boundary::Transmissive,
    };

    let initial_mass: f64 = states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    println!("Initial wet volume (warm start): {:.2e} m³", initial_mass);

    // --- Main event loop ----------------------------------------
    println!(
        "\nSimulating {} days of Atacama 2017 (peak day = 11, Q = {} m³/s)",
        n_days, Q_DAILY_M3S[10]
    );
    println!(
        "{:>3}  {:>5}  {:>7}  {:>11}  {:>9}  {:>6}  {:>8}  {:>8}",
        "day", "Q", "h_max", "mass [m³]", "Δm/Δt", "n_wet", "t_day", "t_total"
    );
    println!(
        "{:>3}  {:>5}  {:>7}  {:>11}  {:>9}  {:>6}  {:>8}  {:>8}",
        "", "m³/s", "[m]", "", "[m³/s]", "", "[s]", "[s]"
    );

    let mut t_total_sim = 0.0;
    let mut cumulative_inflow = 0.0;
    let mut prev_mass = initial_mass;
    let t_loop_start = Instant::now();

    for day in 0..n_days {
        let q_day = Q_DAILY_M3S[day];
        let sources = vec![PointSource {
            row: INFLOW_ROW,
            col: INFLOW_COL,
            q_mass: q_day,
        }];

        let t_day_start = Instant::now();
        let mut t = 0.0;
        let mut steps_day = 0usize;
        while t < SECONDS_PER_DAY {
            let dt = cfl_time_step_with_bcs(&states, &mesh, bcs, CFL).min(SECONDS_PER_DAY - t);
            ssprk2_step(&mut states, &mesh, bcs, dt);
            manning_friction_step(&mut states, MANNING_N, dt, 1.0e-9);
            apply_point_sources(&mut states, &sources, dt, mesh.dx, mesh.dy);
            t += dt;
            steps_day += 1;
        }
        t_total_sim += SECONDS_PER_DAY;
        cumulative_inflow += q_day * SECONDS_PER_DAY;

        let m: f64 = states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
        let max_h = states.iter().map(|s| s.h).fold(0.0_f64, f64::max);
        let n_wet = states.iter().filter(|s| s.h > 0.01).count();
        let dm = m - prev_mass;
        prev_mass = m;

        let t_day = t_day_start.elapsed().as_secs_f64();
        let t_total = t_loop_start.elapsed().as_secs_f64();
        println!(
            "{:>3}  {:>5.1}  {:>7.3}  {:>11.3e}  {:>+9.2}  {:>6}  {:>8.1}  {:>8.1}",
            day + 1,
            q_day,
            max_h,
            m,
            dm / SECONDS_PER_DAY,
            n_wet,
            t_day,
            t_total
        );

        // Write per-day depth snapshot.
        let out_path = PathBuf::from(OUTPUT_DIR).join(format!(
            "huasco_2d_depth_day_{:02}.tif",
            day + 1
        ));
        write_depth_geotiff(&out_path, &states, transform, Some(-9999.0))
            .expect("failed to write per-day depth GeoTIFF");

        let _ = steps_day;
        if !max_h.is_finite() {
            eprintln!("NaN/Inf at day {} — aborting", day + 1);
            break;
        }
    }

    // --- Final mass balance ------------------------------------
    let final_mass: f64 = states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    let net = final_mass - initial_mass;
    let outflow = cumulative_inflow - net;
    println!(
        "\nMass balance over {} days ({:.1} h simulated):",
        n_days, t_total_sim / 3600.0
    );
    println!("  initial wet volume     : {:>12.3e} m³", initial_mass);
    println!("  cumulative inflow      : {:>12.3e} m³", cumulative_inflow);
    println!("  final wet volume       : {:>12.3e} m³", final_mass);
    println!(
        "  net storage change     : {:>12.3e} m³ ({:+.1}% of cumulative)",
        net,
        100.0 * net / cumulative_inflow
    );
    println!(
        "  implied outflow (W)    : {:>12.3e} m³  (= {:.2} m³/s mean)",
        outflow,
        outflow / t_total_sim
    );

    println!("\nTotal wall time: {:.1} min", t_start.elapsed().as_secs_f64() / 60.0);
    println!(
        "Per-day snapshots written to {}/huasco_2d_depth_day_*.tif",
        OUTPUT_DIR
    );
}
