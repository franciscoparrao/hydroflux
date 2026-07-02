//! Measured convergence order of the 2D MUSCL + SSP-RK2 scheme.
//!
//! The 1D suite verifies convergence rates (Stoker, MacDonald); in 2D
//! the audit found only fixed L1 bounds — nothing asserted that the
//! second-order machinery actually delivers super-linear convergence.
//! This test measures self-convergence on a smooth problem: a gentle
//! free-surface Gaussian bump at rest on a flat bed, evolved for a
//! short time (well before any shock forms), on 32², 64² and 128²
//! grids. Successive solutions are compared after block-average
//! restriction (the L1 difference between grid N and grid 2N is the
//! standard Richardson error surrogate), and the observed order
//!
//!   p = log2( e_32 / e_64 )
//!
//! must clear 1.3: comfortably above first order, while allowing for
//! the minmod limiter clipping at the bump extremum (which is what
//! keeps a TVD scheme below the clean p = 2 of unlimited
//! reconstruction).

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2D, Mesh2D, StepWorkspace2D, cfl_time_step, ssprk2_step_with,
};
use ndarray::Array2;

const L_DOMAIN: f64 = 10.0;
const T_END: f64 = 0.25;
const H_BASE: f64 = 1.0;
const H_AMP: f64 = 0.1;

/// Run the smooth-bump problem on an `n × n` grid to `T_END`.
fn run_bump(n: usize) -> Array2<Conserved2D> {
    let d = L_DOMAIN / n as f64;
    let mesh = Mesh2D::new(Array2::<f64>::zeros((n, n)), d, d, 0.0);
    let mut states = Array2::from_shape_fn((n, n), |(i, j)| {
        let x = (j as f64 + 0.5) * d - 0.5 * L_DOMAIN;
        let y = (i as f64 + 0.5) * d - 0.5 * L_DOMAIN;
        let h = H_BASE + H_AMP * (-(x * x + y * y) / 2.0).exp();
        Conserved2D::new(h, 0.0, 0.0)
    });
    let mut ws = StepWorkspace2D::for_mesh(&mesh);
    let mut t = 0.0;
    while t < T_END {
        let dt = cfl_time_step(&states, &mesh, 0.4).min(T_END - t);
        ssprk2_step_with(&mut states, &mesh, Boundaries2D::WALLS, dt, &mut ws);
        t += dt;
    }
    states
}

/// Restrict a `2n × 2n` solution to `n × n` by block averaging (the
/// conservative restriction: each coarse cell is the mean of its four
/// fine children).
fn restrict(fine: &Array2<Conserved2D>) -> Array2<f64> {
    let (rows, cols) = fine.dim();
    Array2::from_shape_fn((rows / 2, cols / 2), |(i, j)| {
        0.25 * (fine[(2 * i, 2 * j)].h
            + fine[(2 * i, 2 * j + 1)].h
            + fine[(2 * i + 1, 2 * j)].h
            + fine[(2 * i + 1, 2 * j + 1)].h)
    })
}

/// Mean absolute difference between a coarse solution and the
/// restriction of the next-finer one.
fn l1_vs_restricted(coarse: &Array2<Conserved2D>, fine: &Array2<Conserved2D>) -> f64 {
    let restricted = restrict(fine);
    let n = coarse.len() as f64;
    coarse
        .indexed_iter()
        .map(|((i, j), s)| (s.h - restricted[(i, j)]).abs())
        .sum::<f64>()
        / n
}

#[test]
fn smooth_bump_self_convergence_is_superlinear() {
    let u32 = run_bump(32);
    let u64 = run_bump(64);
    let u128 = run_bump(128);

    let e_coarse = l1_vs_restricted(&u32, &u64);
    let e_fine = l1_vs_restricted(&u64, &u128);
    assert!(
        e_coarse > 0.0 && e_fine > 0.0,
        "degenerate errors: {e_coarse:.3e}, {e_fine:.3e}"
    );

    let order = (e_coarse / e_fine).log2();
    assert!(
        order > 1.3,
        "MUSCL + SSP-RK2 self-convergence order {order:.3} \
         (e_32/64 = {e_coarse:.4e}, e_64/128 = {e_fine:.4e}) — expected > 1.3; \
         a drop toward 1.0 means the second-order machinery is not engaging"
    );
    // Guard against silently comparing garbage: the absolute error
    // scale must be small compared to the bump amplitude.
    assert!(
        e_fine < 0.05 * H_AMP,
        "fine-grid Richardson error {e_fine:.3e} is out of scale"
    );
}
