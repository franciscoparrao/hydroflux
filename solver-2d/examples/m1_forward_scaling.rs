//! WP6: how the cost of a full gradient scales with the number of
//! parameters under forward-mode AD, and where reverse-mode overtakes it.
//!
//! `Dual {val, dval}` carries a single derivative direction, so the
//! gradient of an objective with respect to `P` parameters costs `P`
//! independent forward passes — one per seeded parameter. Reverse-mode,
//! by contrast, recovers the whole gradient in a single reverse sweep
//! whose cost is a small constant multiple of the primal, independent
//! of `P` [Griewank & Walther 2008].
//!
//! This example measures the forward-mode side of that trade-off
//! directly: the objective is the total water volume after a fixed
//! number of SSP-RK2 + Manning steps on a dam-break, and the parameters
//! are `P` zonal Manning coefficients tiling the domain in vertical
//! stripes. For each `P` we time the full gradient (all `P` passes) and
//! normalise by the `f64` primal cost, which is the quantity that can
//! be compared against a reverse-mode implementation's own multiple.
//!
//! The measured slope pins the constant in `cost(P) / cost_primal ≈ r·P`;
//! the break-even against reverse-mode is then `P* = k / r`, where `k`
//! is reverse-mode's cost multiple. We report `P*` across the standard
//! `k ∈ [3, 5]` band rather than asserting a single figure.
//!
//! Output:
//!   papers/01_review/figures/data/m1_forward_scaling.csv
//!     columns: n_params, wall_clock_s, cost_vs_primal, per_param_cost
//!
//! Run (quiet machine — this is a timing measurement):
//!   cargo run --release -p hydroflux-solver-2d --example m1_forward_scaling

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use hydroflux_autograd::{Dual, Real};
use ndarray::Array2;

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2DG, MaybeSendSync, Mesh2DG, manning_friction_step, ssprk2_step,
};

const OUT_CSV: &str = "papers/01_review/figures/data/m1_forward_scaling.csv";

const N_ROWS: usize = 64;
const N_COLS: usize = 64;
const DX: f64 = 1.0;
const DY: f64 = 1.0;
const N_STEPS: usize = 200;
const DT_FIXED: f64 = 0.01;

const N_BASE: f64 = 0.03;

/// Parameter counts to sweep. Powers of two so the linear trend is
/// visible over a decade without spending time on redundant points.
const PARAM_COUNTS: [usize; 5] = [1, 2, 4, 8, 16];

/// Reverse-mode cost multiples (× primal) bracketing the range usually
/// quoted for a well-implemented adjoint of an explicit stencil code.
const REVERSE_BAND: [f64; 2] = [3.0, 5.0];

/// Build the initial dam-break state. Identical geometry for every
/// backend and every `P`, so the only thing that varies is the
/// arithmetic type and the number of passes.
fn build_states<T: Real>() -> Array2<Conserved2DG<T>> {
    Array2::from_shape_fn((N_ROWS, N_COLS), |(_, j)| {
        let h = if j < N_COLS / 2 { 1.0 } else { 0.1 };
        Conserved2DG::<T>::new_generic(T::from_f64(h), T::zero(), T::zero())
    })
}

/// Which vertical stripe a column belongs to, under a tiling of the
/// domain into `n_params` zones.
fn zone_of(col: usize, n_params: usize) -> usize {
    (col * n_params) / N_COLS
}

/// Uniform Manning field for the `f64` primal run.
fn zonal_manning_f64() -> Array2<f64> {
    Array2::from_elem((N_ROWS, N_COLS), N_BASE)
}

/// Zonal Manning field for one gradient pass: the domain is tiled into
/// `n_params` vertical stripes, and `seed` selects the single stripe
/// carrying `dval = 1`. Every other stripe holds the same value as a
/// constant, so consecutive passes differ only in which partial
/// derivative they propagate — the values, and therefore the primal
/// trajectory, are identical across all `P` passes.
fn zonal_manning_dual(n_params: usize, seed: usize) -> Array2<Dual> {
    Array2::from_shape_fn((N_ROWS, N_COLS), |(_, j)| {
        if zone_of(j, n_params) == seed {
            Dual::variable(N_BASE)
        } else {
            Dual::constant(N_BASE)
        }
    })
}

/// One forward pass: integrate and return the objective (total volume).
fn forward_pass<T: Real + MaybeSendSync>(manning: Array2<T>) -> T {
    let mut states = build_states::<T>();
    let bed = Array2::<T>::from_elem((N_ROWS, N_COLS), T::zero());
    let mesh = Mesh2DG::<T>::with_manning_field(bed, DX, DY, manning);
    let bcs = Boundaries2D::TRANSMISSIVE;

    for _ in 0..N_STEPS {
        ssprk2_step(&mut states, &mesh, bcs, DT_FIXED);
        manning_friction_step(&mut states, &mesh, DT_FIXED, 1.0e-6);
    }

    let cell_area = T::from_f64(DX * DY);
    states
        .iter()
        .fold(T::zero(), |acc, s| acc + s.h * cell_area)
}

fn main() {
    fs::create_dir_all(Path::new(OUT_CSV).parent().unwrap()).expect("create data dir");

    println!(
        "Forward-mode gradient scaling: {N_ROWS}×{N_COLS} cells, \
         {N_STEPS} SSP-RK2 + Manning steps\n"
    );

    // Warm-up both backends so the allocator and page-fault paths are
    // hot before anything is recorded.
    let _ = forward_pass::<f64>(zonal_manning_f64());
    let _ = forward_pass::<Dual>(zonal_manning_dual(1, 0));

    // Primal reference: the f64 run, no derivatives at all.
    let t0 = Instant::now();
    let volume = forward_pass::<f64>(zonal_manning_f64());
    let t_primal = t0.elapsed().as_secs_f64();
    println!("  primal (f64)      {t_primal:>7.3} s   objective = {volume:.6} m³\n");

    println!(
        "  {:>8}  {:>10}  {:>14}  {:>16}",
        "P", "wall [s]", "cost / primal", "per-param cost"
    );

    let mut rows: Vec<(usize, f64, f64, f64)> = Vec::new();
    for &p in PARAM_COUNTS.iter() {
        let t0 = Instant::now();
        let mut grad = Vec::with_capacity(p);
        for seed in 0..p {
            let out = forward_pass::<Dual>(zonal_manning_dual(p, seed));
            grad.push(out.dval);
        }
        let wall = t0.elapsed().as_secs_f64();
        let ratio = wall / t_primal;
        let per_param = ratio / p as f64;
        println!(
            "  {p:>8}  {wall:>10.3}  {ratio:>14.2}×  {per_param:>15.2}×"
        );
        rows.push((p, wall, ratio, per_param));

        assert!(
            grad.iter().all(|g| g.is_finite()),
            "non-finite gradient component at P = {p}"
        );
    }

    // Fit the per-parameter constant `r` as the mean of the measured
    // per-parameter costs — with a scalar Dual the relation is linear
    // by construction, so the mean is the estimator, not a regression.
    let r = rows.iter().map(|r| r.3).sum::<f64>() / rows.len() as f64;
    println!("\n  per-parameter cost r = {r:.2}× primal (mean over P)");
    println!(
        "  break-even vs reverse-mode: P* = {:.1} (k = {:.0}×)  …  {:.1} (k = {:.0}×)",
        REVERSE_BAND[0] / r,
        REVERSE_BAND[0],
        REVERSE_BAND[1] / r,
        REVERSE_BAND[1]
    );

    let mut f = File::create(OUT_CSV).expect("create CSV");
    writeln!(f, "n_params,wall_clock_s,cost_vs_primal,per_param_cost").unwrap();
    for (p, wall, ratio, per_param) in &rows {
        writeln!(f, "{p},{wall:.6},{ratio:.4},{per_param:.4}").unwrap();
    }
    println!("\n  wrote {OUT_CSV}");
}
