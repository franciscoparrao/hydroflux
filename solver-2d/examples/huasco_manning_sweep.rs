//! WP7: one-at-a-time sensitivity of the §4.3 headline result to the
//! land-cover → Manning lookup.
//!
//! §4.2 of the manuscript asserts that "the qualitative direction
//! (riparian vegetation retains water) is robust to plausible lookup
//! choices, the quantitative magnitude is not". This example measures
//! that claim instead of asserting it.
//!
//! The three land-cover classes that actually occupy the channel and
//! its banks are swept one at a time over the plausible range each
//! carries in the standard compilations, with the other two held at
//! the baseline used in §4.2:
//!
//! | class | ESA code | baseline | swept over      |
//! |-------|----------|----------|-----------------|
//! | tree  | 10       | 0.100    | 0.06, 0.10, 0.15 |
//! | shrub | 20       | 0.060    | 0.04, 0.06, 0.08 |
//! | bare  | 60       | 0.025    | 0.02, 0.025, 0.03 |
//!
//! Every configuration is compared against the SAME uniform-`n = 0.04`
//! reference, re-run here so the comparison is internally consistent
//! rather than read off a previous session's log.
//!
//! Reported per configuration: retained wet volume (and its change vs
//! the uniform reference — the "+22 %" headline), mean outflow through
//! the western boundary, mean depth change over channel cells
//! (`acc > 1e6`), wetted-cell count and peak depth.
//!
//! Output:
//!   papers/01_review/figures/data/m_manning_sweep.csv
//!     columns: class, n_value, retained_m3, retained_vs_uniform_pct,
//!              outflow_m3s, outflow_vs_uniform_pct, dh_mean_channel_m,
//!              n_wet, h_max_m
//!
//! Run (~8 one-day simulations, a few minutes each):
//!   cargo run --release -p hydroflux-solver-2d --example huasco_manning_sweep

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use ndarray::Array2;
use surtgis_core::io::read_geotiff;
use surtgis_core::raster::Raster;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, PointSource, apply_point_sources,
    cfl_time_step_with_bcs, manning_friction_step, mesh_from_geotiff,
    mesh_from_geotiff_with_landcover, ssprk2_step,
};

const SUBSET_DEM: &str = "examples/huasco_2d_phase2/data/huasco_subset_dem.tif";
const SUBSET_ACC: &str = "examples/huasco_2d_phase2/data/huasco_subset_acc.tif";
const SUBSET_LC: &str = "examples/huasco_2d_phase2/data/huasco_subset_landcover.tif";
const OUT_CSV: &str = "papers/01_review/figures/data/m_manning_sweep.csv";

const N_UNIFORM: f64 = 0.04;
const N_REF_WARMSTART: f64 = 0.04;
const ACC_THRESHOLD: f64 = 1_000_000.0;
const SLOPE_MEAN: f64 = 0.0074;
const CFL: f64 = 0.4;
const SECONDS_PER_DAY: f64 = 86_400.0;
const INFLOW_ROW: usize = 135;
const INFLOW_COL: usize = 66;
const Q_DAY1: f64 = 17.5;

/// Baseline lookup values for the three channel-adjacent classes.
const TREE_BASE: f64 = 0.100;
const SHRUB_BASE: f64 = 0.060;
const BARE_BASE: f64 = 0.025;

/// Outcome of one one-day simulation.
struct RunResult {
    retained_m3: f64,
    outflow_m3s: f64,
    n_wet: usize,
    h_max: f64,
    depth: Array2<f64>,
}

fn manning_normal_depth(q_m3s: f64, n: f64, slope: f64, cell_width_m: f64) -> f64 {
    let q_per_w = q_m3s / cell_width_m;
    (n * q_per_w / slope.sqrt()).powf(3.0 / 5.0)
}

/// Integrate one day with the given mesh, returning the summary plus
/// the final depth field (needed for the channel-mean depth change).
fn run_one_day(mesh: &Mesh2D, acc: &Array2<f64>) -> RunResult {
    let h_warm = manning_normal_depth(Q_DAY1, N_REF_WARMSTART, SLOPE_MEAN, mesh.dx);
    let mut states = Array2::from_shape_fn((mesh.n_rows(), mesh.n_cols()), |(i, j)| {
        if acc[(i, j)] > ACC_THRESHOLD {
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

    let cell_area = mesh.dx * mesh.dy;
    let initial_mass: f64 = states.iter().map(|s| s.h * cell_area).sum();

    let sources = vec![PointSource {
        row: INFLOW_ROW,
        col: INFLOW_COL,
        q_mass: Q_DAY1,
    }];

    let mut t = 0.0;
    while t < SECONDS_PER_DAY {
        let dt = cfl_time_step_with_bcs(&states, mesh, bcs, CFL).min(SECONDS_PER_DAY - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        manning_friction_step(&mut states, mesh, dt, 1.0e-9);
        apply_point_sources(&mut states, &sources, dt, mesh.dx, mesh.dy);
        t += dt;
    }

    let final_mass: f64 = states.iter().map(|s| s.h * cell_area).sum();
    let cumulative_inflow = Q_DAY1 * SECONDS_PER_DAY;
    let outflow = cumulative_inflow - (final_mass - initial_mass);

    RunResult {
        retained_m3: final_mass,
        outflow_m3s: outflow / SECONDS_PER_DAY,
        n_wet: states.iter().filter(|s| s.h > 0.01).count(),
        h_max: states.iter().map(|s| s.h).fold(0.0_f64, f64::max),
        depth: Array2::from_shape_fn(states.dim(), |(i, j)| states[(i, j)].h),
    }
}

/// Mean depth difference over channel cells only.
fn channel_mean_delta(a: &Array2<f64>, b: &Array2<f64>, acc: &Array2<f64>) -> f64 {
    let mut sum = 0.0;
    let mut n = 0usize;
    for ((i, j), &acc_v) in acc.indexed_iter() {
        if acc_v > ACC_THRESHOLD {
            sum += a[(i, j)] - b[(i, j)];
            n += 1;
        }
    }
    if n == 0 { 0.0 } else { sum / n as f64 }
}

fn main() {
    fs::create_dir_all(Path::new(OUT_CSV).parent().unwrap()).expect("create data dir");

    let acc_raster: Raster<f64> =
        read_geotiff(SUBSET_ACC, None).expect("failed to load flow-accumulation subset");
    let acc = acc_raster.data().to_owned();

    // Uniform reference, re-run here so every percentage below is
    // relative to a number produced in this same session.
    let (mesh_u, _) =
        mesh_from_geotiff(SUBSET_DEM, N_UNIFORM).expect("failed to load DEM subset");
    println!("Uniform reference (n = {N_UNIFORM}) …");
    let uniform = run_one_day(&mesh_u, &acc);
    println!(
        "  retained {:.4e} m³, outflow {:.2} m³/s, n_wet {}, h_max {:.3} m\n",
        uniform.retained_m3, uniform.outflow_m3s, uniform.n_wet, uniform.h_max
    );

    // One-at-a-time configurations. The baseline appears once, not
    // three times, so the CSV has no duplicated rows.
    let mut configs: Vec<(&str, f64, f64, f64, f64)> = vec![
        ("baseline", f64::NAN, TREE_BASE, SHRUB_BASE, BARE_BASE),
    ];
    for v in [0.06, 0.15] {
        configs.push(("tree", v, v, SHRUB_BASE, BARE_BASE));
    }
    for v in [0.04, 0.08] {
        configs.push(("shrub", v, TREE_BASE, v, BARE_BASE));
    }
    for v in [0.02, 0.03] {
        configs.push(("bare", v, TREE_BASE, SHRUB_BASE, v));
    }
    // One-at-a-time never varies two classes together, so it cannot by
    // itself support a claim about the whole range. The two corners do
    // the adversarial work: "all minimum" is the configuration most
    // likely to flip the sign, since it drives the channel classes
    // toward the uniform n = 0.04 the comparison is against.
    configs.push(("corner-min", f64::NAN, 0.06, 0.04, 0.02));
    configs.push(("corner-max", f64::NAN, 0.15, 0.08, 0.03));

    println!(
        "{:>9} {:>7} {:>12} {:>10} {:>10} {:>9} {:>7} {:>8}",
        "class", "n", "retained", "vs unif", "outflow", "vs unif", "n_wet", "dh_chan"
    );

    let mut rows: Vec<String> = Vec::new();
    for (class, n_swept, tree, shrub, bare) in configs {
        let lookup = move |code: u8| -> f64 {
            match code {
                10 => tree,
                20 => shrub,
                60 => bare,
                30 => 0.040,
                40 => 0.035,
                50 => 0.015,
                70 | 80 => 0.030,
                90 => 0.050,
                95 => 0.100,
                100 => 0.045,
                _ => 0.040,
            }
        };
        let (mesh, _) = mesh_from_geotiff_with_landcover(SUBSET_DEM, SUBSET_LC, lookup)
            .expect("failed to load DEM + landcover");

        let r = run_one_day(&mesh, &acc);
        let d_ret = 100.0 * (r.retained_m3 / uniform.retained_m3 - 1.0);
        let d_out = 100.0 * (r.outflow_m3s / uniform.outflow_m3s - 1.0);
        let dh = channel_mean_delta(&r.depth, &uniform.depth, &acc);

        let n_label = if n_swept.is_nan() {
            "—".to_string()
        } else {
            format!("{n_swept:.3}")
        };
        println!(
            "{class:>9} {n_label:>7} {:>12.4e} {d_ret:>+9.1}% {:>10.2} {d_out:>+8.1}% {:>7} {dh:>+8.3}",
            r.retained_m3, r.outflow_m3s, r.n_wet
        );
        rows.push(format!(
            "{class},{n_label},{:.6e},{d_ret:.3},{:.4},{d_out:.3},{dh:.5},{},{:.4}",
            r.retained_m3, r.outflow_m3s, r.n_wet, r.h_max
        ));
    }

    let mut f = File::create(OUT_CSV).expect("create CSV");
    writeln!(
        f,
        "class,n_value,retained_m3,retained_vs_uniform_pct,outflow_m3s,\
         outflow_vs_uniform_pct,dh_mean_channel_m,n_wet,h_max_m"
    )
    .unwrap();
    writeln!(
        f,
        "uniform,{N_UNIFORM},{:.6e},0.000,{:.4},0.000,0.00000,{},{:.4}",
        uniform.retained_m3, uniform.outflow_m3s, uniform.n_wet, uniform.h_max
    )
    .unwrap();
    for r in &rows {
        writeln!(f, "{r}").unwrap();
    }
    println!("\n  wrote {OUT_CSV}");
}
