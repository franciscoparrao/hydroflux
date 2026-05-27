//! Atacama 2017 event with **spatially varying** Manning roughness
//! derived from an ESA WorldCover landcover raster.
//!
//! Mirrors [`huasco_2d_event`] except for the Manning field:
//!
//! - `huasco_2d_event.rs`        — uniform `n = 0.04`, single
//!   calibrated value from the 1D autograd run.
//! - `huasco_2d_event_landcover.rs` (this example) — `n(x, y)` from
//!   `huasco_subset_landcover.tif` (200×67 ESA WorldCover codes
//!   resampled to the DEM grid) mapped through
//!   [`esa_worldcover_to_manning`]. Channel cells (bare ground,
//!   code 60) get `n = 0.025`, while riparian tree cover (code 10)
//!   gets `n = 0.10` — a 4× contrast that the uniform run misses.
//!
//! All other physics is identical: same DEM, same BCs (W = Transmissive,
//! N/S/E = Wall), same time-varying daily Q from DGA 03820003, same
//! Manning-normal-depth warm-start (using a mean reference `n_ref`).
//! The PointSource is at the same channel cell at the E edge.
//!
//! # How to reproduce the landcover raster
//!
//! ```text
//! gdalwarp \
//!   -t_srs EPSG:32719 \
//!   -te 333620.42 6826925.80 335630.42 6832925.80 \
//!   -tr 30 30 -r mode -ot Byte -overwrite \
//!   /vsicurl/https://esa-worldcover.s3.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_S30W072_Map.tif \
//!   examples/huasco_2d_phase2/output/huasco_subset_landcover.tif
//! ```
//!
//! # Reproducir
//!
//! ```text
//! cargo run --release -p hydroflux-solver-2d --example huasco_2d_event_landcover -- --days 1
//! cargo run --release -p hydroflux-solver-2d --example huasco_2d_event_landcover            # full 21 days
//! ```

use std::env;
use std::path::PathBuf;
use std::time::Instant;

use ndarray::Array2;
use surtgis_core::io::read_geotiff;
use surtgis_core::raster::Raster;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, PointSource, apply_point_sources, cfl_time_step_with_bcs,
    esa_worldcover_to_manning, manning_friction_step, mesh_from_geotiff_with_landcover,
    ssprk2_step, write_depth_geotiff,
};

const SUBSET_DEM: &str = "examples/huasco_2d_phase2/output/huasco_subset_dem.tif";
const SUBSET_ACC: &str = "examples/huasco_2d_phase2/output/huasco_subset_acc.tif";
const SUBSET_LC: &str = "examples/huasco_2d_phase2/output/huasco_subset_landcover.tif";
const OUTPUT_DIR: &str = "examples/huasco_2d_phase2/output";

const N_REF: f64 = 0.04; // reference Manning for warm-start depth only
const ACC_THRESHOLD: f64 = 1_000_000.0;
const SLOPE_MEAN: f64 = 0.0074;
const CFL: f64 = 0.4;
const SECONDS_PER_DAY: f64 = 86_400.0;
const INFLOW_ROW: usize = 135;
const INFLOW_COL: usize = 66;

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

    let (mesh, transform) =
        mesh_from_geotiff_with_landcover(SUBSET_DEM, SUBSET_LC, esa_worldcover_to_manning)
            .expect("failed to load DEM + landcover; check huasco_subset_landcover.tif exists");

    let n_min = mesh
        .manning
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let n_max = mesh
        .manning
        .iter()
        .copied()
        .fold(0.0_f64, f64::max);
    let n_mean = mesh.manning.iter().sum::<f64>() / mesh.manning.len() as f64;
    println!(
        "DEM: {}×{} cells, dx={} m, dy={} m",
        mesh.n_rows(),
        mesh.n_cols(),
        mesh.dx,
        mesh.dy
    );
    println!(
        "Manning field: n_min={:.4}, n_mean={:.4}, n_max={:.4} (variable from landcover)",
        n_min, n_mean, n_max
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

    // Warm-start uses N_REF (uniform mean) for the depth estimate;
    // the per-cell Manning field still governs the per-step friction.
    let h_warm = manning_normal_depth(Q_DAILY_M3S[0], N_REF, SLOPE_MEAN, mesh.dx);
    println!(
        "\nManning warm-start: Q_day1={} m³/s, n_ref={}, S₀={}, dx={} m → h_n={:.3} m",
        Q_DAILY_M3S[0], N_REF, SLOPE_MEAN, mesh.dx, h_warm
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
            manning_friction_step(&mut states, &mesh, dt, 1.0e-9);
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

        let out_path = PathBuf::from(OUTPUT_DIR).join(format!(
            "huasco_2d_depth_day_{:02}_landcover.tif",
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

    let final_mass: f64 = states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    let net = final_mass - initial_mass;
    let outflow = cumulative_inflow - net;
    println!(
        "\nMass balance over {} days ({:.1} h simulated):",
        n_days,
        t_total_sim / 3600.0
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
        "Per-day snapshots written to {}/huasco_2d_depth_day_*_landcover.tif",
        OUTPUT_DIR
    );
}
