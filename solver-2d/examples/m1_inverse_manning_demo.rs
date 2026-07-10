//! M1b: single-variable inverse demo for Manning n on the full 2D
//! solver via forward-mode AD.
//!
//! Setup: a `(3 × 60)`-cell flat-bed channel with a depth jump at
//! the centre (`h_L = 1 m`, `h_R = 0.1 m`), walls on the long sides,
//! transmissive ends. The dam-break wave propagates outward; Manning
//! friction damps it. We pick `n_true = 0.040`, run the forward
//! solver for `t = 20 s`, measure the depth-weighted centroid of the
//! wet region as the "observation", then start from `n_init = 0.080`
//! and recover `n_true` by Newton iteration on the squared error. Each
//! iteration evaluates a single forward pass with `T = Dual<f64>`;
//! the derivative of the loss w.r.t. `n` is read directly from `.dval`.
//!
//! This is the direct evidence the EMS reviewer M1 asked for: the
//! differentiability wedge is delivered end-to-end (forward solver +
//! gradient) on the 2D code, not just at the per-kernel level. The
//! observable (depth-weighted centroid) is sensitive to friction
//! through wave damping; a smaller n means a less-damped wave with a
//! larger spread, so the centroid sits closer to the leading edge.
//!
//! Output:
//!   papers/01_review/figures/data/m1_inverse_manning.csv
//!     columns: iter, n_estimate, centroid, loss, dcentroid_dn, step
//!
//! Run (release for speed; Dual carries ~2× cost over f64):
//!   cargo run --release -p hydroflux-solver-2d --example m1_inverse_manning_demo

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use hydroflux_autograd::{Dual, Real};
use ndarray::Array2;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2DG, MaybeSendSync, Mesh2DG, cfl_time_step,
    manning_friction_step, ssprk2_step,
};

const OUT_CSV: &str = "papers/01_review/figures/data/m1_inverse_manning.csv";

const N_ROWS: usize = 3;
const N_COLS: usize = 60;
const DX: f64 = 1.0; // 60 m channel
const DY: f64 = 1.0;
const H_LEFT: f64 = 1.0;
const H_RIGHT: f64 = 0.1;
const T_END: f64 = 20.0;
const N_TRUE: f64 = 0.040;
const N_INIT: f64 = 0.080;
const ITERATIONS: usize = 12;

/// Build a flat-bed mesh in the chosen scalar type.
fn build_mesh<T: Real>(n: T) -> Mesh2DG<T> {
    let bed = Array2::<T>::from_elem((N_ROWS, N_COLS), T::zero());
    Mesh2DG::<T>::new(bed, DX, DY, n)
}

/// Initial dam-break condition: water of depth `h_L` for `j < N_COLS/2`,
/// `h_R` elsewhere; everything at rest.
fn initial_states<T: Real>() -> Array2<Conserved2DG<T>> {
    let h_l = T::from_f64(H_LEFT);
    let h_r = T::from_f64(H_RIGHT);
    Array2::from_shape_fn((N_ROWS, N_COLS), |(_i, j)| {
        let h = if j < N_COLS / 2 { h_l } else { h_r };
        Conserved2DG::<T>::new_generic(h, T::zero(), T::zero())
    })
}

/// Run the dam-break for `T_END` seconds and return the depth-weighted
/// centroid of the wet region — `Σ x · h / Σ h` — along the long axis.
/// Generic over T so the same code evaluates with f64 (for the
/// synthetic observation) and with Dual (for AD-based inversion).
fn depth_centroid<T: Real + MaybeSendSync>(n: T) -> T {
    let mesh = build_mesh(n);
    let mut states = initial_states::<T>();
    let bcs = Boundaries2D {
        north: Boundary::Wall,
        south: Boundary::Wall,
        west: Boundary::Transmissive,
        east: Boundary::Transmissive,
    };
    let mut t = 0.0_f64;
    let mut n_steps = 0usize;
    while t < T_END && n_steps < 5000 {
        let dt = cfl_time_step(&states, &mesh, 0.4).min(T_END - t);
        if !dt.is_finite() || dt <= 0.0 {
            break;
        }
        ssprk2_step(&mut states, &mesh, bcs, dt);
        manning_friction_step(&mut states, &mesh, dt, 1.0e-6);
        t += dt;
        n_steps += 1;
    }
    // Depth-weighted centroid along x. Only the middle row contributes
    // — that pins the test to the channel centreline and reduces the
    // contribution of wet/dry transitions near the walls.
    let i_mid = N_ROWS / 2;
    let mut num = T::zero();
    let mut den = T::zero();
    for j in 0..N_COLS {
        let x = T::from_f64((j as f64 + 0.5) * DX);
        let h = states[(i_mid, j)].h;
        num = num + h * x;
        den = den + h;
    }
    num / den
}

fn main() {
    fs::create_dir_all(Path::new(OUT_CSV).parent().unwrap()).expect("create data dir");

    let c_target = depth_centroid::<f64>(N_TRUE);
    println!("Synthetic observation");
    println!("  n_true               = {N_TRUE:.6}");
    println!("  depth centroid x [m] = {c_target:.6}");

    let mut n = N_INIT;
    let mut log = Vec::<(usize, f64, f64, f64, f64, f64)>::new();
    println!("\nInversion (Newton on depth centroid)");
    println!(
        "  iter  n_estimate     centroid       loss            dc/dn        step"
    );
    for iter in 0..ITERATIONS {
        let n_dual = Dual::variable(n);
        let c_pred = depth_centroid::<Dual>(n_dual);
        let err = c_pred.val - c_target;
        let loss = err * err;
        let dc_dn = c_pred.dval;
        let step = if dc_dn.abs() > 1e-12 {
            -err / dc_dn
        } else {
            0.0
        };
        println!(
            "  {:>4}  {:>10.6}  {:>10.6}  {:>14.3e}  {:>10.3}  {:>10.6}",
            iter, n, c_pred.val, loss, dc_dn, step
        );
        log.push((iter, n, c_pred.val, loss, dc_dn, step));
        n += step;
        n = n.clamp(1e-4, 1.0);
    }

    let mut f = File::create(OUT_CSV).expect("create CSV");
    writeln!(f, "iter,n_estimate,centroid,loss,dc_dn,step").unwrap();
    for (it, n_e, c, l, d, s) in &log {
        writeln!(f, "{it},{n_e:.10},{c:.10},{l:.10e},{d:.10e},{s:.10e}").unwrap();
    }

    let (_, n_final, _, _, _, _) = *log.last().unwrap();
    let rel_err = (n_final - N_TRUE).abs() / N_TRUE;
    println!("\nResult");
    println!("  n_final  = {n_final:.6}");
    println!("  n_true   = {N_TRUE:.6}");
    println!("  rel err  = {:.3e}", rel_err);
    println!("  wrote {OUT_CSV}");
}
