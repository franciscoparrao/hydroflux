//! Closed-domain variant of the Huasco day-1 run, for the SynxFlow
//! cross-validation diagnostic (`docs/xval-synxflow-huasco-results.md`).
//!
//! The open-domain comparison showed excellent agreement on inundation
//! extent (CSI 0.935, zero false negatives) but a systematic depth
//! bias: SynxFlow stored ~31 % more water. A mass-balance inference
//! pointed at the outlet — SynxFlow's `open` boundary appearing to
//! evacuate less than hydroflux's `Transmissive` — rather than at the
//! interior scheme. That is a hypothesis, and this example is how we
//! test it.
//!
//! Sealing every boundary removes the outlet as a degree of freedom:
//! total stored mass is then pinned by the inflow alone and must agree
//! between the two solvers by conservation. Whatever depth difference
//! survives is purely how each scheme *distributes* that water. If the
//! bias collapses here, the open-domain discrepancy was the boundary
//! condition and the interior comparison is clean.
//!
//! Identical to `huasco_2d_event_landcover --days 1` except that the
//! western edge is Wall rather than Transmissive.
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example huasco_closed_domain

use std::path::PathBuf;

use ndarray::Array2;
use surtgis_core::io::read_geotiff;
use surtgis_core::raster::Raster;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, PointSource, apply_point_sources,
    cfl_time_step_with_bcs, esa_worldcover_to_manning, manning_friction_step,
    mesh_from_geotiff_with_landcover, ssprk2_step, write_depth_geotiff,
};

const SUBSET_DEM: &str = "examples/huasco_2d_phase2/data/huasco_subset_dem.tif";
const SUBSET_ACC: &str = "examples/huasco_2d_phase2/data/huasco_subset_acc.tif";
const SUBSET_LC: &str = "examples/huasco_2d_phase2/data/huasco_subset_landcover.tif";
const OUTPUT_DIR: &str = "examples/huasco_2d_phase2/output";

const N_REF: f64 = 0.04;
const ACC_THRESHOLD: f64 = 1_000_000.0;
const SLOPE_MEAN: f64 = 0.0074;
const CFL: f64 = 0.4;
const SECONDS_PER_DAY: f64 = 86_400.0;
const INFLOW_ROW: usize = 135;
const INFLOW_COL: usize = 66;
const Q_DAY1: f64 = 17.5;

fn manning_normal_depth(q_m3s: f64, n: f64, slope: f64, cell_width_m: f64) -> f64 {
    let q_per_w = q_m3s / cell_width_m;
    (n * q_per_w / slope.sqrt()).powf(3.0 / 5.0)
}

fn main() {
    // `--cfl X` exists to separate two things the cross-validation of §4.5
    // would otherwise conflate: how much of the disagreement with an
    // independent solver is the spatial scheme, and how much is our choice
    // of time step. Both codes converge as the CFL number falls, so if the
    // residual shrinks with CFL it is temporal discretisation; if it
    // plateaus, it is the scheme.
    let cfl = {
        let a: Vec<String> = std::env::args().collect();
        a.iter().position(|x| x == "--cfl")
            .and_then(|i| a.get(i + 1))
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(CFL)
    };

    let (mesh, transform) =
        mesh_from_geotiff_with_landcover(SUBSET_DEM, SUBSET_LC, esa_worldcover_to_manning)
            .expect("failed to load DEM + landcover");

    let acc: Raster<f64> =
        read_geotiff(SUBSET_ACC, None).expect("failed to load flow accumulation");
    let acc_data = acc.data();

    let h_warm = manning_normal_depth(Q_DAY1, N_REF, SLOPE_MEAN, mesh.dx);
    let mut states = Array2::from_shape_fn((mesh.n_rows(), mesh.n_cols()), |(i, j)| {
        if acc_data[(i, j)] > ACC_THRESHOLD {
            Conserved2D::new(h_warm, 0.0, 0.0)
        } else {
            Conserved2D::DRY
        }
    });

    // The one difference from the open run: every edge is a wall.
    let bcs = Boundaries2D::WALLS;

    let cell_area = mesh.dx * mesh.dy;
    let initial_mass: f64 = states.iter().map(|s| s.h * cell_area).sum();
    println!("Closed-domain Huasco, day 1 (all boundaries Wall), CFL = {cfl}");
    println!("  warm start h_n = {h_warm:.3} m, initial volume {initial_mass:.4e} m³");

    // `--no-inflow` gives the cleanest possible scheme-to-scheme
    // comparison: a sealed domain, an identical initial water body, and
    // no source at all, so a full day of integration is pure
    // redistribution. It exists because the two codes realise an
    // upstream discharge differently — hydroflux as an exact volumetric
    // point source, SynxFlow by converting the discharge series to
    // boundary velocities — and that difference alone moves the
    // delivered volume by ~8 %, which would otherwise contaminate the
    // comparison.
    let no_inflow = std::env::args().any(|a| a == "--no-inflow");

    let q = if no_inflow { 0.0 } else { Q_DAY1 };
    let sources = vec![PointSource {
        row: INFLOW_ROW,
        col: INFLOW_COL,
        q_mass: q,
    }];

    let mut t = 0.0;
    let mut steps = 0usize;
    while t < SECONDS_PER_DAY {
        let dt = cfl_time_step_with_bcs(&states, &mesh, bcs, cfl).min(SECONDS_PER_DAY - t);
        ssprk2_step(&mut states, &mesh, bcs, dt);
        manning_friction_step(&mut states, &mesh, dt, 1.0e-9);
        apply_point_sources(&mut states, &sources, dt, mesh.dx, mesh.dy);
        t += dt;
        steps += 1;
    }

    let final_mass: f64 = states.iter().map(|s| s.h * cell_area).sum();
    let cumulative_inflow = q * SECONDS_PER_DAY;
    let expected = initial_mass + cumulative_inflow;
    let h_max = states.iter().map(|s| s.h).fold(0.0_f64, f64::max);
    let n_wet = states.iter().filter(|s| s.h > 0.01).count();

    println!("  steps {steps}, t = {t:.1} s");
    println!("  final volume    {final_mass:.6e} m³");
    println!("  expected        {expected:.6e} m³  (initial + inflow)");
    println!(
        "  closure error   {:+.3e} m³ ({:+.2e} relative)",
        final_mass - expected,
        (final_mass - expected) / expected
    );
    println!("  h_max {h_max:.4} m, wet cells (h > 0.01) {n_wet}");

    let tag = format!("{:.0}", cfl * 100.0);
    let name = if no_inflow {
        format!("huasco_closed_noflow_cfl{tag}.tif")
    } else {
        format!("huasco_closed_day_01_cfl{tag}.tif")
    };
    let out = PathBuf::from(OUTPUT_DIR).join(name);
    write_depth_geotiff(&out, &states, transform, Some(-9999.0)).expect("write depth");
    println!("  wrote {}", out.display());
}
