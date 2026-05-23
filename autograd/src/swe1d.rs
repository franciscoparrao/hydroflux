//! Minimal generic 1D shallow-water step.
//!
//! Explicit Lax-Friedrichs scheme on a uniform cell-centred mesh,
//! written over [`Real`] so it evaluates end-to-end with both `f64`
//! and [`Dual`](crate::Dual). The point is to exercise full
//! Saint-Venant arithmetic (flux + bed-slope + Manning friction) in
//! a setting small enough to differentiate in seconds.
//!
//! This is NOT a substitute for `hydroflux-solver-1d`: that crate
//! has the HLL Riemann solver, MUSCL reconstruction, well-balanced
//! Audusse source, and is tuned for production-grade benchmarks.
//! The module here exists to demonstrate gradient-based calibration
//! and to feed the 2028 Q1 paper application case; if the demo
//! reveals a need for higher-fidelity gradients we migrate (or copy)
//! the production solver to be generic.
//!
//! # Discretisation
//!
//! State `U_i = (h_i, q_i)` at cell centres. Flux at face `i+½`:
//!
//! ```text
//! F*_{i+½} = ½·(F(U_i) + F(U_{i+1})) − ½·α·(U_{i+1} − U_i)
//! ```
//!
//! with `α = max_i(|q_i/h_i| + √(g·h_i))` (global Lax-Friedrichs).
//! The cell update is
//!
//! ```text
//! U_i^{n+1} = U_i^n − (dt/dx)·(F*_{i+½} − F*_{i−½}) + dt·S_i
//! ```
//!
//! with source `S_i = (0, −g·h·dz/dx − g·n²·q·|q|/h^{7/3})`. Bed
//! slope is centred-difference on the interior, one-sided at the
//! boundaries. Friction is applied semi-implicitly to avoid stiffness
//! at small `h` (point-implicit, single Newton step is exact for the
//! linearised Manning term).
//!
//! Boundary conditions for the demo are minimal: ghost cells mirror
//! the boundary state with optional Dirichlet override on the
//! upstream side (`h_left = h_inflow`, `q_left = q_inflow`).

use crate::Real;

/// Upstream boundary specification.
#[derive(Debug, Clone, Copy)]
pub enum LeftBc<T: Real> {
    /// Dirichlet on both `h` and `q`. Use this to drive steady inflow.
    Dirichlet { h: T, q: T },
    /// Transmissive (zero-gradient).
    Transmissive,
}

/// Downstream boundary specification. Currently transmissive only;
/// the demo runs long-enough channels that the steady solution is
/// established before the wave reaches the outlet.
#[derive(Debug, Clone, Copy)]
pub enum RightBc {
    /// Transmissive (zero-gradient).
    Transmissive,
}

/// One explicit step of length `dt` on `(h, q)`.
///
/// The buffers `h_next` and `q_next` are written in place; `bed`
/// stores cell-centre bed elevations.
pub fn lax_friedrichs_step<T: Real>(
    h: &[T],
    q: &[T],
    bed: &[f64],
    dx: f64,
    dt: f64,
    manning_n: T,
    gravity: f64,
    left_bc: LeftBc<T>,
    _right_bc: RightBc,
    h_next: &mut [T],
    q_next: &mut [T],
) {
    let n = h.len();
    assert_eq!(q.len(), n);
    assert_eq!(bed.len(), n);
    assert_eq!(h_next.len(), n);
    assert_eq!(q_next.len(), n);
    assert!(n >= 2);

    // Largest physical wave speed (global LF dissipation).
    let mut alpha = T::zero();
    for i in 0..n {
        let c = (h[i].max(T::zero()) * gravity).sqrt();
        let u = q[i] / h[i].max(T::from_f64(1.0e-9));
        let a = u.abs() + c;
        if a.value() > alpha.value() {
            alpha = a;
        }
    }

    // Ghost states: dirichlet (or transmissive) on the left,
    // transmissive on the right.
    let (h_l_ghost, q_l_ghost) = match left_bc {
        LeftBc::Dirichlet { h: hg, q: qg } => (hg, qg),
        LeftBc::Transmissive => (h[0], q[0]),
    };
    let (h_r_ghost, q_r_ghost) = (h[n - 1], q[n - 1]);

    // Helper: physical flux F(U) = (q, q²/h + ½·g·h²).
    let flux = |hi: T, qi: T| -> (T, T) {
        let h_safe = hi.max(T::from_f64(1.0e-12));
        let mass = qi;
        let mom = (qi * qi) / h_safe + (hi * hi) * (0.5 * gravity);
        (mass, mom)
    };

    // LF flux at the face between (hL, qL) and (hR, qR).
    let lf_face = |hl: T, ql: T, hr: T, qr: T| -> (T, T) {
        let (fl_h, fl_q) = flux(hl, ql);
        let (fr_h, fr_q) = flux(hr, qr);
        let half_alpha = alpha * 0.5;
        let f_h = (fl_h + fr_h) * 0.5 - half_alpha * (hr - hl);
        let f_q = (fl_q + fr_q) * 0.5 - half_alpha * (qr - ql);
        (f_h, f_q)
    };

    // Faces 0 .. n (n+1 of them); face k separates cell k-1 and cell k,
    // where cell -1 and cell n are the ghosts.
    let mut f_h_left = lf_face(h_l_ghost, q_l_ghost, h[0], q[0]).0;
    let mut f_q_left = lf_face(h_l_ghost, q_l_ghost, h[0], q[0]).1;

    let dt_over_dx = dt / dx;
    let g = gravity;

    for i in 0..n {
        // Right face of cell i: between cell i and cell i+1 (or ghost).
        let (h_r, q_r) = if i + 1 < n {
            (h[i + 1], q[i + 1])
        } else {
            (h_r_ghost, q_r_ghost)
        };
        let (f_h_right, f_q_right) = lf_face(h[i], q[i], h_r, q_r);

        // Bed slope (central difference, one-sided at boundaries).
        let dz_dx = if i == 0 {
            (bed[1] - bed[0]) / dx
        } else if i == n - 1 {
            (bed[n - 1] - bed[n - 2]) / dx
        } else {
            (bed[i + 1] - bed[i - 1]) / (2.0 * dx)
        };

        // Conservative update (mass + momentum) without friction yet.
        let h_new = h[i] - (f_h_right - f_h_left) * dt_over_dx;
        let q_star = q[i] - (f_q_right - f_q_left) * dt_over_dx + h[i] * (-g * dz_dx) * dt;

        // Point-implicit Manning friction:
        //   q^{n+1} = q*  − dt · g · n² · q^{n+1} · |q*| / h^{7/3}
        // ⇒ q^{n+1} · (1 + dt·g·n²·|q*|/h^{7/3}) = q*.
        // h^{7/3} is taken from the post-mass-update h_new (clamped).
        let h_clamp = h_new.max(T::from_f64(1.0e-9));
        let coeff = manning_n * manning_n * (dt * g) * q_star.abs() / h_clamp.powf(7.0 / 3.0);
        let q_new = q_star / (coeff + 1.0);

        h_next[i] = h_new.max(T::zero());
        q_next[i] = if h_next[i].value() > 1.0e-9 {
            q_new
        } else {
            T::zero()
        };

        // Slide faces.
        f_h_left = f_h_right;
        f_q_left = f_q_right;
    }
}

/// CFL-bounded time step for a state.
pub fn cfl_dt<T: Real>(h: &[T], q: &[T], dx: f64, gravity: f64, cfl: f64) -> f64 {
    let mut max_lambda = 0.0_f64;
    for i in 0..h.len() {
        let h_v = h[i].value().max(0.0);
        if h_v < 1.0e-9 {
            continue;
        }
        let c = (gravity * h_v).sqrt();
        let u = q[i].value() / h_v;
        let a = u.abs() + c;
        if a > max_lambda {
            max_lambda = a;
        }
    }
    if max_lambda < 1.0e-12 {
        // No physical waves; advance by a default that does not blow up
        // any reasonable initial perturbation.
        return cfl * dx / 1.0;
    }
    cfl * dx / max_lambda
}

/// Run the explicit scheme until `t_end`. Returns final `(h, q)` and
/// number of steps taken.
#[allow(clippy::too_many_arguments)]
pub fn run<T: Real>(
    h0: Vec<T>,
    q0: Vec<T>,
    bed: &[f64],
    dx: f64,
    t_end: f64,
    manning_n: T,
    gravity: f64,
    cfl: f64,
    left_bc: LeftBc<T>,
    right_bc: RightBc,
) -> (Vec<T>, Vec<T>, usize) {
    let mut h = h0;
    let mut q = q0;
    let mut h_next = vec![T::zero(); h.len()];
    let mut q_next = vec![T::zero(); q.len()];
    let mut t = 0.0;
    let mut steps = 0;
    while t < t_end {
        let dt = cfl_dt(&h, &q, dx, gravity, cfl).min(t_end - t);
        lax_friedrichs_step(
            &h,
            &q,
            bed,
            dx,
            dt,
            manning_n,
            gravity,
            left_bc,
            right_bc,
            &mut h_next,
            &mut q_next,
        );
        std::mem::swap(&mut h, &mut h_next);
        std::mem::swap(&mut q, &mut q_next);
        t += dt;
        steps += 1;
        if steps > 500_000 {
            panic!("swe1d::run did not finish in 500k steps");
        }
    }
    (h, q, steps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Dual, physics::manning_normal_depth};

    const G: f64 = 9.81;

    fn linspace_bed(n: usize, slope: f64, dx: f64) -> Vec<f64> {
        // Bed descending in +x at constant slope. Cell-centre x at
        // (i + 0.5)·dx so the leftmost cell is highest.
        (0..n).map(|i| -slope * (i as f64 + 0.5) * dx).collect()
    }

    #[test]
    fn lake_at_rest_stays_at_rest_with_f64() {
        // Flat bed, zero discharge, walls (transmissive ghosts mirror
        // the boundary). The scheme must preserve the still state to
        // round-off.
        let n_cells = 50;
        let bed = vec![0.0; n_cells];
        let h0: Vec<f64> = vec![1.0; n_cells];
        let q0: Vec<f64> = vec![0.0; n_cells];
        let (h, q, _) = run(
            h0,
            q0,
            &bed,
            1.0,
            10.0,
            0.03,
            G,
            0.4,
            LeftBc::Transmissive,
            RightBc::Transmissive,
        );
        for i in 0..n_cells {
            assert!((h[i] - 1.0).abs() < 1.0e-9, "h[{i}] = {}", h[i]);
            assert!(q[i].abs() < 1.0e-9, "q[{i}] = {}", q[i]);
        }
    }

    #[test]
    fn steady_inflow_relaxes_toward_manning_normal_depth_f64() {
        // Sloping channel + constant inflow → after enough time, the
        // depth in the interior should be close to Manning normal
        // depth. Tolerance is loose: LF is dissipative and the demo
        // mesh is coarse, but the value should be in the right ball.
        let n_cells = 80;
        let dx = 2.0;
        let slope = 0.001;
        let bed = linspace_bed(n_cells, slope, dx);
        let n = 0.04;
        let q_in = 1.5_f64;
        let h_n = manning_normal_depth(q_in, n, slope);
        let h0: Vec<f64> = vec![h_n; n_cells];
        let q0: Vec<f64> = vec![q_in; n_cells];
        let (h, q, _) = run(
            h0,
            q0,
            &bed,
            dx,
            300.0,
            n,
            G,
            0.4,
            LeftBc::Dirichlet { h: h_n, q: q_in },
            RightBc::Transmissive,
        );
        // Interior (skip ghost-influenced ends) should be near h_n.
        let mid = n_cells / 2;
        let h_mid = h[mid];
        let q_mid = q[mid];
        assert!(
            (h_mid / h_n - 1.0).abs() < 0.10,
            "h_mid = {h_mid:.4}, h_n = {h_n:.4}, ratio = {:.3}",
            h_mid / h_n
        );
        assert!(
            (q_mid / q_in - 1.0).abs() < 0.10,
            "q_mid = {q_mid:.4}, q_in = {q_in:.4}"
        );
    }

    #[test]
    fn gradient_d_h_steady_d_manning_n_matches_normal_depth_derivative() {
        // h_n(n) = (n·q/√S₀)^(3/5). Analytical: dh_n/dn = (3/5)·(q/√S₀)^(3/5)·n^(-2/5).
        // Run the solver with `n` as Dual::variable; the interior
        // depth gradient w.r.t. n should match analytical to within
        // the LF dissipation tolerance.
        let n_cells = 80;
        let dx = 2.0;
        let slope = 0.001;
        let bed = linspace_bed(n_cells, slope, dx);
        let n_val = 0.04_f64;
        let q_in_val = 1.5_f64;
        let h_n_val = manning_normal_depth(q_in_val, n_val, slope);

        let n_dual = Dual::variable(n_val);
        let q_in_dual = Dual::constant(q_in_val);
        // Initial condition: at the normal depth (so the system stays
        // near steady; we mostly check that the scheme keeps it there
        // and that the gradient threads through).
        let h0: Vec<Dual> = vec![Dual::constant(h_n_val); n_cells];
        let q0: Vec<Dual> = vec![Dual::constant(q_in_val); n_cells];
        // Boundary depth depends on n (Manning normal depth as a Dual).
        let h_bc = (n_dual * Dual::constant(q_in_val) / Dual::constant(slope.sqrt()))
            .powf(3.0 / 5.0);
        let (h, _q, _) = run(
            h0,
            q0,
            &bed,
            dx,
            500.0,
            n_dual,
            G,
            0.4,
            LeftBc::Dirichlet {
                h: h_bc,
                q: q_in_dual,
            },
            RightBc::Transmissive,
        );

        // Analytical derivative of h_n w.r.t. n.
        let analytic_grad =
            (3.0 / 5.0) * (q_in_val / slope.sqrt()).powf(3.0 / 5.0) * n_val.powf(-2.0 / 5.0);

        // Interior cell gradient.
        let mid = n_cells / 2;
        let grad_sim = h[mid].dval;

        // Tolerance: LF dissipation + finite time → loose 25% match
        // is plenty to demonstrate that AD threads through the
        // time-stepping. Production solver should hit <1%.
        let rel_err = (grad_sim / analytic_grad - 1.0).abs();
        assert!(
            rel_err < 0.25,
            "dh/dn sim = {grad_sim:.4}, analytic = {analytic_grad:.4}, rel_err = {rel_err:.3}"
        );
    }
}
