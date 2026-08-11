//! Gradient-based calibration of the riparian (tree-class) Manning
//! coefficient against an independently computed depth field.
//!
//! §4.4 shows that the roughness problem on this reach is
//! low-dimensional: of the three land-cover classes occupying the
//! channel and its banks, one carries essentially all of the
//! sensitivity (elasticity 0.261, against 0.022 and 0.001). That is
//! exactly the regime in which the forward-mode gradient of §2.5 is the
//! right tool — the break-even against reverse mode sits at P ≈ 2, and
//! here P = 1.
//!
//! This example closes that loop. The target is the SynxFlow depth
//! field from the sealed-domain cross-validation of §4.5, shipped with
//! the repository. The objective is the sum of squared depth
//! differences over cells wet in either field:
//!
//!   J(n) = Σ_i ( h_i(n) − h_i^target )²
//!
//! A single forward pass with `T = Dual` yields, for every cell,
//! both `h_i` and `∂h_i/∂n`. That gives the gradient and the
//! Gauss-Newton curvature from the *same* evaluation:
//!
//!   dJ/dn  = 2 Σ_i r_i (∂h_i/∂n)
//!   d²J/dn² ≈ 2 Σ_i (∂h_i/∂n)²        (Gauss-Newton)
//!
//! so each iteration costs one solve, not two. This is the property
//! that makes forward-mode worthwhile at low parameter count: the
//! per-cell sensitivities come out of the primal pass for free, and a
//! spatially distributed misfit turns them into a curvature estimate
//! without a second-order model.
//!
//! Honest scope: the target is another *model*, not an observation, so
//! what this demonstrates is the calibration machinery operating on a
//! spatially distributed target — not physical identification of a
//! roughness value. Recovering the value that makes two independent
//! solvers agree is a statement about the codes, not about the river.
//!
//! Output:
//!   papers/01_review/figures/data/m_calibrate_tree.csv
//!     columns: iter, n_tree, loss, dJ_dn, gn_curvature, step
//!
//! Run (release; each iteration is one sealed 1-day solve in Dual):
//!   cargo run --release -p hydroflux-solver-2d --example huasco_calibrate_tree

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use hydroflux_autograd::{Dual, Real};
use ndarray::Array2;
use surtgis_core::io::read_geotiff;
use surtgis_core::raster::Raster;

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2DG, MaybeSendSync, Mesh2DG, cfl_time_step_with_bcs,
    manning_friction_step, read_ascii_grid, ssprk2_step,
};

const SUBSET_DEM: &str = "examples/huasco_2d_phase2/data/huasco_subset_dem.tif";
const SUBSET_ACC: &str = "examples/huasco_2d_phase2/data/huasco_subset_acc.tif";
const SUBSET_LC: &str = "examples/huasco_2d_phase2/data/huasco_subset_landcover.tif";
const TARGET: &str = "examples/huasco_2d_phase2/data/synxflow_noflow_h_86400.asc.gz";
const OUT_CSV: &str = "papers/01_review/figures/data/m_calibrate_tree.csv";

const ACC_THRESHOLD: f64 = 1_000_000.0;
const H_WARM: f64 = 0.457;
const CFL: f64 = 0.4;
const T_END: f64 = 86_400.0;
const H_THR: f64 = 0.01;

/// Deliberately wrong starting point: the low end of the plausible tree
/// range from §4.4, well away from the 0.100 baseline.
const N_INIT: f64 = 0.060;
const MAX_ITER: usize = 6;
const TOL: f64 = 1.0e-5;

/// Land-cover → Manning, with the tree class (ESA code 10) left free.
/// Every other class is held at the §4.2 lookup as a constant, so the
/// derivative that propagates is exactly ∂/∂n_tree.
fn manning_for<T: Real>(code: u8, n_tree: T) -> T {
    match code {
        10 => n_tree,
        20 => T::from_f64(0.060),
        30 => T::from_f64(0.040),
        40 => T::from_f64(0.035),
        50 => T::from_f64(0.015),
        60 => T::from_f64(0.025),
        70 | 80 => T::from_f64(0.030),
        90 => T::from_f64(0.050),
        95 => T::from_f64(0.100),
        100 => T::from_f64(0.045),
        _ => T::from_f64(0.040),
    }
}

struct Inputs {
    dem: Array2<f64>,
    codes: Array2<u8>,
    acc: Array2<f64>,
    dx: f64,
    dy: f64,
}

fn load() -> Inputs {
    let dem_r: Raster<f64> = read_geotiff(SUBSET_DEM, None).expect("DEM");
    let lc_r: Raster<f64> = read_geotiff(SUBSET_LC, None).expect("landcover");
    let acc_r: Raster<f64> = read_geotiff(SUBSET_ACC, None).expect("acc");
    let t = dem_r.transform();
    let codes = lc_r.data().mapv(|v| v as u8);
    Inputs {
        dem: dem_r.data().to_owned(),
        codes,
        acc: acc_r.data().to_owned(),
        dx: t.pixel_width.abs(),
        dy: t.pixel_height.abs(),
    }
}

/// One sealed, source-free day. Returns the final depth field in `T`,
/// so a Dual run carries ∂h_i/∂n_tree in every cell.
fn forward<T: Real + MaybeSendSync>(inp: &Inputs, n_tree: T) -> Array2<T> {
    let (nr, nc) = inp.dem.dim();
    let bed = Array2::from_shape_fn((nr, nc), |(i, j)| T::from_f64(inp.dem[(i, j)]));
    let manning =
        Array2::from_shape_fn((nr, nc), |(i, j)| manning_for::<T>(inp.codes[(i, j)], n_tree));
    let mesh = Mesh2DG::<T>::with_manning_field(bed, inp.dx, inp.dy, manning);

    let mut states = Array2::from_shape_fn((nr, nc), |(i, j)| {
        if inp.acc[(i, j)] > ACC_THRESHOLD {
            Conserved2DG::<T>::new_generic(T::from_f64(H_WARM), T::zero(), T::zero())
        } else {
            Conserved2DG::<T>::new_generic(T::zero(), T::zero(), T::zero())
        }
    });

    let bcs = Boundaries2D::WALLS;
    let mut t = 0.0;
    while t < T_END {
        let dt = cfl_time_step_with_bcs(&states, &mesh, bcs, CFL).min(T_END - t);
        ssprk2_step(&mut states, &mesh, bcs, dt);
        manning_friction_step(&mut states, &mesh, dt, 1.0e-9);
        t += dt;
    }
    Array2::from_shape_fn((nr, nc), |(i, j)| states[(i, j)].h)
}

fn main() {
    fs::create_dir_all(Path::new(OUT_CSV).parent().unwrap()).expect("data dir");
    let inp = load();
    let (target, _hdr) = read_ascii_grid(TARGET).expect("SynxFlow target field");
    assert_eq!(target.dim(), inp.dem.dim(), "target grid mismatch");

    println!(
        "Calibrating the tree-class Manning against the SynxFlow field\n\
         Domain {}×{} at {} m, sealed, no source, t_end = {} s\n\
         Start n_tree = {N_INIT} (deliberately wrong; §4.2 baseline is 0.100)\n",
        inp.dem.dim().0,
        inp.dem.dim().1,
        inp.dx,
        T_END
    );
    println!(
        "{:>4}  {:>9}  {:>12}  {:>13}  {:>13}  {:>10}",
        "iter", "n_tree", "loss", "dJ/dn", "GN curvature", "step"
    );

    let mut n = N_INIT;
    let mut rows: Vec<String> = Vec::new();
    let t0 = Instant::now();

    for it in 0..MAX_ITER {
        // One Dual pass: value and per-cell sensitivity together.
        let h = forward::<Dual>(&inp, Dual::variable(n));

        let mut loss = 0.0;
        let mut grad = 0.0;
        let mut curv = 0.0;
        for ((i, j), hd) in h.indexed_iter() {
            let tgt = target[(i, j)];
            // Cells dry in both fields carry no information and would
            // otherwise dilute the misfit with a large block of zeros.
            if hd.val <= H_THR && tgt <= H_THR {
                continue;
            }
            let r = hd.val - tgt;
            loss += r * r;
            grad += 2.0 * r * hd.dval;
            curv += 2.0 * hd.dval * hd.dval;
        }

        // Gauss-Newton step. The curvature is positive semi-definite by
        // construction, so this cannot step uphill; guard the
        // degenerate case where the parameter has no local influence.
        let step = if curv > 1.0e-12 { -grad / curv } else { 0.0 };
        println!(
            "{it:>4}  {n:>9.6}  {loss:>12.6e}  {grad:>+13.5e}  {curv:>13.5e}  {step:>+10.6}"
        );
        rows.push(format!("{it},{n:.8},{loss:.8e},{grad:.8e},{curv:.8e},{step:.8e}"));

        if step.abs() < TOL {
            println!("\nconverged: |step| = {:.2e} < {TOL:.0e}", step.abs());
            n += step;
            break;
        }
        n += step;
    }

    println!("\n  final n_tree = {n:.6}");
    println!("  §4.2 lookup value = 0.100");
    println!("  wall time = {:.1} min", t0.elapsed().as_secs_f64() / 60.0);

    let mut f = File::create(OUT_CSV).expect("create CSV");
    writeln!(f, "iter,n_tree,loss,dJ_dn,gn_curvature,step").unwrap();
    for r in &rows {
        writeln!(f, "{r}").unwrap();
    }
    println!("  wrote {OUT_CSV}");
}
