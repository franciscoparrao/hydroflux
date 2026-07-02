//! End-to-end forward-mode AD through the production 1D solver.
//!
//! This is the point of making `hydroflux-solver-1d` generic over
//! `Real`: the calibration workflow no longer needs the first-order
//! Lax-Friedrichs stand-ins in `hydroflux-autograd` to obtain
//! gradients — the HLL + Audusse well-balanced production scheme is
//! differentiable as-is. The test drives a steady Manning flow to
//! equilibrium and compares `∂h/∂n` from a seeded `Dual` run against
//! a central finite difference of two `f64` runs.

use hydroflux_autograd::{Dual, Real, physics::manning_normal_depth};
use hydroflux_solver_1d::{
    Boundaries, Boundary, Channel1D, Conserved, cfl_time_step, forward_euler_step,
    manning_friction_step,
};
use ndarray::Array1;

const CFL: f64 = 0.4;

/// Run the production solver (FV step + operator-split friction) to
/// `t_end` and return the depth at `probe`.
fn steady_depth_at<T: Real>(
    n_cells: usize,
    dx: f64,
    slope: f64,
    manning: T,
    q_in: f64,
    h_init: f64,
    t_end: f64,
    probe: usize,
) -> T {
    let bed = Array1::from_iter((0..n_cells).map(|i| -slope * (i as f64 + 0.5) * dx));
    let channel = Channel1D::new(bed, dx, manning);
    let mut states: Vec<Conserved<T>> = (0..n_cells)
        .map(|_| Conserved::new(T::from_f64(h_init), T::from_f64(q_in)))
        .collect();

    // Discharge upstream, transmissive downstream: the interior relaxes
    // to the Manning normal depth of the *run's own* roughness rather
    // than being anchored by a fixed downstream stage (which would damp
    // ∂h/∂n via the backwater profile).
    let bcs = Boundaries {
        left: Boundary::Discharge { q: q_in },
        right: Boundary::Transmissive,
    };

    let mut t = 0.0;
    let mut steps = 0usize;
    while t < t_end {
        let dt = cfl_time_step(&states, dx, CFL).min(t_end - t);
        assert!(dt.is_finite() && dt > 0.0, "degenerate dt = {dt}");
        forward_euler_step(&mut states, &channel, bcs, dt);
        manning_friction_step(&mut states, channel.manning, dt, 1e-9);
        t += dt;
        steps += 1;
        assert!(steps < 500_000, "run did not finish in 500k steps");
    }
    states[probe].h
}

#[test]
fn ad_gradient_of_steady_depth_wrt_manning_matches_finite_difference() {
    let n_cells = 80;
    let dx = 2.0;
    let slope = 1.0e-3;
    let n0 = 0.04_f64;
    let q_in = 1.5_f64;
    let t_end = 2000.0;
    let probe = n_cells / 2;

    // Initial (and downstream-BC) depth: analytical Manning normal
    // depth at the BASELINE n, held fixed across perturbations so FD
    // and AD measure the same quantity (BC and IC enter the Dual run
    // as constants).
    let h_n = manning_normal_depth(q_in, n0, slope);

    let eps = 1.0e-5;
    let h_plus = steady_depth_at::<f64>(n_cells, dx, slope, n0 + eps, q_in, h_n, t_end, probe);
    let h_minus = steady_depth_at::<f64>(n_cells, dx, slope, n0 - eps, q_in, h_n, t_end, probe);
    let grad_fd = (h_plus - h_minus) / (2.0 * eps);

    let h_dual = steady_depth_at::<Dual>(
        n_cells,
        dx,
        slope,
        Dual::variable(n0),
        q_in,
        h_n,
        t_end,
        probe,
    );
    let grad_ad = h_dual.dval;

    // Sanity: the steady interior depth tracks the analytical normal
    // depth, so the gradient must be positive and of the order of
    // dh_n/dn = (3/5)·h_n/n.
    let grad_scale = 0.6 * h_n / n0;
    assert!(
        grad_fd > 0.2 * grad_scale && grad_fd < 5.0 * grad_scale,
        "FD gradient implausible: {grad_fd:.4e} vs analytical scale {grad_scale:.4e}"
    );

    // Tolerance 5e-3, not machine-level: the CFL dt is state-dependent
    // and deliberately not differentiated, and the semi-implicit
    // friction factor depends on dt through |q*| — the same O(dt)
    // frozen-dt residual quantified in
    // `hydroflux-autograd::power_law_swe1d` (test
    // `ad_gradient_matches_central_finite_difference_for_n_c_p`).
    // Genuine Dual-rule or threading bugs show up at O(1).
    let rel_err = (grad_ad - grad_fd).abs() / grad_fd.abs();
    assert!(
        rel_err < 5.0e-3,
        "AD vs FD mismatch through the production solver: \
         ad = {grad_ad:.6e}, fd = {grad_fd:.6e}, rel_err = {rel_err:.3e}"
    );
}

#[test]
fn dual_constant_run_matches_f64_bitwise() {
    // The T = Dual code path with constant seeds must reproduce the
    // f64 primal exactly over a full production-solver run — same
    // invariant the autograd LF steppers assert via proptest, now on
    // the HLL well-balanced scheme.
    let n_cells = 40;
    let dx = 2.0;
    let slope = 1.0e-3;
    let n0 = 0.04;
    let q_in = 1.5;
    let h_n = manning_normal_depth(q_in, n0, slope);
    let t_end = 100.0;

    let h_f = steady_depth_at::<f64>(n_cells, dx, slope, n0, q_in, h_n, t_end, n_cells / 2);
    let h_d = steady_depth_at::<Dual>(
        n_cells,
        dx,
        slope,
        Dual::constant(n0),
        q_in,
        h_n,
        t_end,
        n_cells / 2,
    );
    assert_eq!(h_f, h_d.val, "Dual::constant diverged from f64 primal");
    assert_eq!(h_d.dval, 0.0, "constant seed produced non-zero derivative");
}
