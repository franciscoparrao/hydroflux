//! Diagnostic: does the forward-mode tangent survive a long integration
//! with active wet/dry fronts on real terrain?
//!
//! Motivation. The gradient-based calibration of the tree-class Manning
//! (`huasco_calibrate_tree`) returned dJ/dn ≈ -4e133 and a Gauss-Newton
//! curvature ≈ 9e270 after a sealed 1-day solve — not a gradient, an
//! overflow. The AD-vs-FD locking suite of §2.5 agrees to better than
//! 1e-6, but those tests integrate O(100) steps on smooth synthetic
//! problems. This run takes ~78,000 steps over a 30 m DEM with a
//! moving shoreline.
//!
//! This example measures where the two regimes part company: it
//! integrates the same sealed Huasco configuration in `Dual`, reporting
//! the maximum |dval| across the field at intervals, so exponential
//! growth (if that is what happens) is visible as a straight line on a
//! log axis against step count.
//!
//! Output:
//!   papers/01_review/figures/data/m_dual_growth.csv
//!     columns: step, t_seconds, max_abs_dval, max_h, n_wet
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example diag_dual_growth

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use hydroflux_autograd::{Dual, Real};
use ndarray::Array2;
use surtgis_core::io::read_geotiff;
use surtgis_core::raster::Raster;

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2DG, Mesh2DG, cfl_time_step_with_bcs, manning_friction_step,
    ssprk2_step,
};

const SUBSET_DEM: &str = "examples/huasco_2d_phase2/data/huasco_subset_dem.tif";
const SUBSET_ACC: &str = "examples/huasco_2d_phase2/data/huasco_subset_acc.tif";
const SUBSET_LC: &str = "examples/huasco_2d_phase2/data/huasco_subset_landcover.tif";
const OUT_CSV: &str = "papers/01_review/figures/data/m_dual_growth.csv";

const ACC_THRESHOLD: f64 = 1_000_000.0;
const H_WARM: f64 = 0.457;
const CFL: f64 = 0.4;
const T_END: f64 = 86_400.0;
const N_TREE: f64 = 0.100;
const REPORT_EVERY: usize = 250;

fn manning_for(code: u8, n_tree: Dual) -> Dual {
    match code {
        10 => n_tree,
        20 => Dual::constant(0.060),
        30 => Dual::constant(0.040),
        40 => Dual::constant(0.035),
        50 => Dual::constant(0.015),
        60 => Dual::constant(0.025),
        70 | 80 => Dual::constant(0.030),
        90 => Dual::constant(0.050),
        95 => Dual::constant(0.100),
        100 => Dual::constant(0.045),
        _ => Dual::constant(0.040),
    }
}

fn main() {
    fs::create_dir_all(Path::new(OUT_CSV).parent().unwrap()).expect("data dir");

    let dem_r: Raster<f64> = read_geotiff(SUBSET_DEM, None).expect("DEM");
    let lc_r: Raster<f64> = read_geotiff(SUBSET_LC, None).expect("landcover");
    let acc_r: Raster<f64> = read_geotiff(SUBSET_ACC, None).expect("acc");
    let t = dem_r.transform();
    let (dx, dy) = (t.pixel_width.abs(), t.pixel_height.abs());
    let dem = dem_r.data();
    let codes = lc_r.data().mapv(|v| v as u8);
    let acc = acc_r.data();
    let (nr, nc) = dem.dim();

    let n_tree = Dual::variable(N_TREE);
    let bed = Array2::from_shape_fn((nr, nc), |(i, j)| Dual::constant(dem[(i, j)]));
    let manning = Array2::from_shape_fn((nr, nc), |(i, j)| manning_for(codes[(i, j)], n_tree));
    let mesh = Mesh2DG::<Dual>::with_manning_field(bed, dx, dy, manning);

    let mut states = Array2::from_shape_fn((nr, nc), |(i, j)| {
        if acc[(i, j)] > ACC_THRESHOLD {
            Conserved2DG::<Dual>::new_generic(Dual::constant(H_WARM), Dual::zero(), Dual::zero())
        } else {
            Conserved2DG::<Dual>::new_generic(Dual::zero(), Dual::zero(), Dual::zero())
        }
    });

    let bcs = Boundaries2D::WALLS;
    let mut f = File::create(OUT_CSV).expect("csv");
    writeln!(f, "step,t_seconds,max_abs_dval,max_h,n_wet").unwrap();

    println!("{:>8}  {:>10}  {:>14}  {:>8}  {:>6}", "step", "t [s]", "max |dh/dn|", "max h", "n_wet");

    let mut t_sim = 0.0;
    let mut step = 0usize;
    while t_sim < T_END {
        let dt = cfl_time_step_with_bcs(&states, &mesh, bcs, CFL).min(T_END - t_sim);
        ssprk2_step(&mut states, &mesh, bcs, dt);
        manning_friction_step(&mut states, &mesh, dt, 1.0e-9);
        t_sim += dt;
        step += 1;

        if step % REPORT_EVERY == 0 || step == 1 {
            let mut max_d = 0.0_f64;
            let mut max_h = 0.0_f64;
            let mut n_wet = 0usize;
            for s in states.iter() {
                if s.h.dval.abs() > max_d {
                    max_d = s.h.dval.abs();
                }
                if s.h.val > max_h {
                    max_h = s.h.val;
                }
                if s.h.val > 0.01 {
                    n_wet += 1;
                }
            }
            println!("{step:>8}  {t_sim:>10.1}  {max_d:>14.4e}  {max_h:>8.4}  {n_wet:>6}");
            writeln!(f, "{step},{t_sim:.3},{max_d:.6e},{max_h:.6},{n_wet}").unwrap();
            if !max_d.is_finite() || max_d > 1.0e100 {
                println!("\nAborting: |dh/dn| exceeded 1e100 — the tangent has diverged.");
                break;
            }
        }
    }
    println!("\nwrote {OUT_CSV}");
}
