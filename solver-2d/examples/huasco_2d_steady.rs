//! Phase 2 first real-data 2D run: Atacama 2017 peak Q = 38.9 m³/s
//! sustained for 1 day on a 6 km × 2 km subset of the Huasco DEM
//! around Santa Juana.
//!
//! # Setup (see `examples/huasco_2d_phase2/extract_subset.py`)
//!
//! - Domain: 200 rows × 67 cols at 30 m → 6 km (N-S) × 2 km (E-W).
//! - DEM: SRTM 30 m pit-filled.
//! - River geometry within window: enters at E edge (row 135, col 66),
//!   exits at W edge (row 25, col 0). River runs roughly diagonally
//!   SE → NW inside the subset.
//! - BCs: West = Transmissive (outflow), all other edges = Wall.
//! - Inflow: single `PointSource` at (135, 66) with mass rate
//!   `Q = 38.9 m³/s` (peak Atacama 2017 at Santa Juana, DGA station
//!   03820003).
//! - Warm start: thin film `h = 1 cm` on channel cells (identified
//!   from the co-extracted flow-accumulation raster with threshold
//!   1×10⁶ cells ≈ 900 km² catchment); dry on banks (`h = 0`).
//! - Manning: uniform `n = 0.04` (gravel-bed Andean, Chow 1959).
//! - Simulation time: 86 400 s (1 day → equilibrium; channel residence
//!   time at u ≈ 2 m/s is ~50 min, so 1 day is ~30× equilibration).
//! - Output: depth GeoTIFF + mass-balance + timing report.
//!
//! Run from the repo root:
//! ```text
//! cargo run --release -p hydroflux-solver-2d --example huasco_2d_steady
//! ```

use std::path::Path;
use std::time::Instant;

use ndarray::Array2;
use surtgis_core::io::read_geotiff;
use surtgis_core::raster::Raster;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, PointSource, apply_point_sources, cfl_time_step_with_bcs,
    manning_friction_step, mesh_from_geotiff, ssprk2_step, write_depth_geotiff,
};

const SUBSET_DEM: &str = "examples/huasco_2d_phase2/data/huasco_subset_dem.tif";
const SUBSET_ACC: &str = "examples/huasco_2d_phase2/data/huasco_subset_acc.tif";
const OUTPUT_DEPTH: &str = "examples/huasco_2d_phase2/output/huasco_2d_depth.tif";

const MANNING_N: f64 = 0.04;
const Q_PEAK: f64 = 38.9; // m³/s — Atacama 2017 peak at Santa Juana
const ACC_THRESHOLD: f64 = 1_000_000.0; // ≈ 900 km² catchment
const WARM_START_H: f64 = 0.01; // 1 cm thin film on channel
const T_END: f64 = 86_400.0; // 1 day = 86400 s
const CFL: f64 = 0.4;

// Inflow cell (single PointSource at the channel cell on E edge).
// Identified by extract_subset.py: (row 135, col 66).
const INFLOW_ROW: usize = 135;
const INFLOW_COL: usize = 66;

fn main() {
    let t_start = Instant::now();

    // --- Load DEM and build mesh ----------------------------------
    let (mesh, transform) = mesh_from_geotiff(SUBSET_DEM, MANNING_N)
        .expect("failed to load DEM subset; run extract_subset.py first");
    println!(
        "Loaded DEM: {} rows × {} cols, dx={} m, dy={} m, manning={}",
        mesh.n_rows(),
        mesh.n_cols(),
        mesh.dx,
        mesh.dy,
        MANNING_N
    );
    println!(
        "  bed elev range: [{:.2}, {:.2}] m",
        mesh.bed.iter().copied().fold(f64::INFINITY, f64::min),
        mesh.bed.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    );

    // --- Load flow_accumulation to identify channel cells ---------
    let acc: Raster<f64> =
        read_geotiff(SUBSET_ACC, None).expect("failed to load flow_accumulation subset");
    let acc_data = acc.data();
    assert_eq!(acc_data.dim(), mesh.bed.dim(), "acc and DEM shapes differ");
    let mut n_channel = 0usize;
    for &v in acc_data {
        if v > ACC_THRESHOLD {
            n_channel += 1;
        }
    }
    println!(
        "Channel cells (acc > {:.0e}): {} / {} ({:.1}%)",
        ACC_THRESHOLD,
        n_channel,
        mesh.n_cells(),
        100.0 * n_channel as f64 / mesh.n_cells() as f64
    );

    // --- Warm-start: thin film on channel cells, dry elsewhere ----
    let mut states = Array2::from_shape_fn((mesh.n_rows(), mesh.n_cols()), |(i, j)| {
        if acc_data[(i, j)] > ACC_THRESHOLD {
            Conserved2D::new(WARM_START_H, 0.0, 0.0)
        } else {
            Conserved2D::DRY
        }
    });

    // --- Boundaries: W = Transmissive (outflow), rest = Wall ------
    let bcs = Boundaries2D {
        north: Boundary::Wall,
        south: Boundary::Wall,
        east: Boundary::Wall,
        west: Boundary::Transmissive,
    };

    // --- Single PointSource at the E-edge channel cell ------------
    let sources = vec![PointSource {
        row: INFLOW_ROW,
        col: INFLOW_COL,
        q_mass: Q_PEAK,
    }];
    println!(
        "Inflow PointSource at ({}, {}): Q = {} m³/s",
        INFLOW_ROW, INFLOW_COL, Q_PEAK
    );

    let initial_mass: f64 = states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    println!(
        "Initial wet volume: {:.1} m³ ({} channel cells × {:.2} m × {} m²)",
        initial_mass,
        n_channel,
        WARM_START_H,
        mesh.dx * mesh.dy
    );

    // --- Time integration -----------------------------------------
    println!(
        "\nSimulating {} s ({:.1} h) with CFL = {}, ssprk2 + Manning + PointSource",
        T_END, T_END / 3600.0, CFL
    );
    let mut t = 0.0;
    let mut steps = 0usize;
    let mut next_report = 600.0; // report every 10 sim-minutes
    let t_sim_start = Instant::now();
    while t < T_END {
        let dt = cfl_time_step_with_bcs(&states, &mesh, bcs, CFL).min(T_END - t);
        ssprk2_step(&mut states, &mesh, bcs, dt);
        manning_friction_step(&mut states, &mesh, dt, 1.0e-9);
        apply_point_sources(&mut states, &sources, dt, mesh.dx, mesh.dy);
        t += dt;
        steps += 1;

        if t >= next_report || t >= T_END {
            let m: f64 = states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
            let max_h = states
                .iter()
                .map(|s| s.h)
                .fold(0.0_f64, f64::max);
            let n_wet = states.iter().filter(|s| s.h > 0.01).count();
            println!(
                "t = {:>6.0} s ({:>5.2} h)  steps = {:>6}  mass = {:>10.2e} m³  max h = {:>5.2} m  wet cells = {}",
                t,
                t / 3600.0,
                steps,
                m,
                max_h,
                n_wet
            );
            next_report = (t / 600.0).floor() * 600.0 + 600.0;
        }

        if !states[(INFLOW_ROW, INFLOW_COL)].h.is_finite() {
            eprintln!(
                "NaN/Inf at inflow cell after {} steps, t = {} — aborting",
                steps, t
            );
            break;
        }
    }
    let sim_elapsed = t_sim_start.elapsed();
    println!(
        "\nSimulation complete: {} steps in {:.1} s wall time ({:.0} steps/s)",
        steps,
        sim_elapsed.as_secs_f64(),
        steps as f64 / sim_elapsed.as_secs_f64()
    );

    // --- Mass balance ---------------------------------------------
    let final_mass: f64 = states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    let cumulative_in = Q_PEAK * T_END;
    let net = final_mass - initial_mass;
    let outflow = cumulative_in - net;
    println!("\nMass balance:");
    println!("  initial wet volume     : {:>10.2e} m³", initial_mass);
    println!("  cumulative inflow      : {:>10.2e} m³  (= {} m³/s × {} s)", cumulative_in, Q_PEAK, T_END);
    println!("  final wet volume       : {:>10.2e} m³", final_mass);
    println!(
        "  net storage change     : {:>10.2e} m³  ({:+.1}% of cumulative in)",
        net,
        100.0 * net / cumulative_in
    );
    println!(
        "  implied outflow (W)    : {:>10.2e} m³  (≡ {:.2} m³/s mean)",
        outflow,
        outflow / T_END
    );

    // --- Write depth raster ---------------------------------------
    write_depth_geotiff(Path::new(OUTPUT_DEPTH), &states, transform, Some(-9999.0))
        .expect("failed to write depth GeoTIFF");
    println!("\nDepth raster written: {}", OUTPUT_DEPTH);
    println!("Total wall time: {:.1} s", t_start.elapsed().as_secs_f64());
}
