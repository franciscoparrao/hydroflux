//! MacDonald inverse-design benchmark with a variable depth profile.
//!
//! Prescribe a smooth depth profile `h(x) = h_base + amp · sin(2π x / L)`
//! and a constant unit discharge `q`. The steady 1D SWE with Manning
//! friction yields the bed slope that makes the prescribed `h(x)` a
//! steady state:
//!
//! ```text
//!   dz/dx = −(1 − Fr²) · dh/dx − Sf
//!   with  Fr² = q²/(g h³),  Sf = n² q² / h^(10/3)
//! ```
//!
//! Integrating `dz/dx` over the domain gives `z(x)`. The solver, run with
//! this bed, Manning `n`, `Discharge(q)` upstream and `Depth(h(L))`
//! downstream, and initial state equal to the analytical profile, should
//! preserve the profile up to first-order discretisation error.
//!
//! Validates the full Q3 pipeline together: Audusse well-balanced source,
//! Manning friction, physical inflow/outflow BCs.

use hydroflux_solver_1d::{
    Boundaries, Boundary, Channel1D, Conserved, cfl_time_step, forward_euler_step,
    manning_friction_step,
};
use ndarray::Array1;
use std::f64::consts::PI;

const G: f64 = 9.81;
const Q: f64 = 1.0; // unit discharge [m²/s]
const N_MANNING: f64 = 0.03;
const L_DOMAIN: f64 = 50.0; // [m]
const H_BASE: f64 = 1.0; // [m]
const H_AMP: f64 = 0.2; // [m], so h ∈ [0.8, 1.2]

fn h_analytical(x: f64) -> f64 {
    H_BASE + H_AMP * (2.0 * PI * x / L_DOMAIN).sin()
}

fn dh_dx_analytical(x: f64) -> f64 {
    H_AMP * (2.0 * PI / L_DOMAIN) * (2.0 * PI * x / L_DOMAIN).cos()
}

/// Manning friction slope at depth `h` for the prescribed `q` and `n`.
fn friction_slope(h: f64) -> f64 {
    N_MANNING * N_MANNING * Q * Q / h.powf(10.0 / 3.0)
}

/// Bed slope `dz/dx` required for `h(x)` to be a steady state of the SWE.
fn dz_dx_analytical(x: f64) -> f64 {
    let h = h_analytical(x);
    let fr_sq = Q * Q / (G * h * h * h);
    -(1.0 - fr_sq) * dh_dx_analytical(x) - friction_slope(h)
}

/// Integrate `z(x)` from `z(0) = 0` via the trapezoidal rule over
/// `n_pts` sub-intervals. Returns `z` at the `n_pts + 1` nodes
/// `x_k = k · (L / n_pts)`.
fn analytical_bed_profile(n_pts: usize) -> Vec<f64> {
    let dx = L_DOMAIN / n_pts as f64;
    let mut z = vec![0.0; n_pts + 1];
    for k in 0..n_pts {
        let x_k = k as f64 * dx;
        let x_k1 = (k + 1) as f64 * dx;
        z[k + 1] = z[k] + 0.5 * (dz_dx_analytical(x_k) + dz_dx_analytical(x_k1)) * dx;
    }
    z
}

/// Build a Channel1D for the solver grid by sampling the analytical bed
/// at cell centres `x_i = (i + 0.5)·dx`. The bed is integrated on a 10×
/// finer mesh and linearly interpolated to the cell centres, so the
/// integration error is well below the solver's discretisation error.
fn build_channel(n_cells: usize) -> Channel1D {
    let dx = L_DOMAIN / n_cells as f64;
    let n_fine = n_cells * 10;
    let z_fine = analytical_bed_profile(n_fine);
    let dx_fine = L_DOMAIN / n_fine as f64;
    let bed = Array1::from_iter((0..n_cells).map(|i| {
        let x_c = (i as f64 + 0.5) * dx;
        let k_f = (x_c / dx_fine).floor() as usize;
        let k = k_f.min(n_fine - 1);
        let t = (x_c - k as f64 * dx_fine) / dx_fine;
        z_fine[k] * (1.0 - t) + z_fine[k + 1] * t
    }));
    Channel1D::new(bed, dx, N_MANNING)
}

fn analytical_state_at_cells(n_cells: usize) -> Vec<Conserved> {
    let dx = L_DOMAIN / n_cells as f64;
    (0..n_cells)
        .map(|i| {
            let x_c = (i as f64 + 0.5) * dx;
            Conserved::new(h_analytical(x_c), Q)
        })
        .collect()
}

/// Run the solver from the analytical initial condition for `t_end` and
/// return `(L1(h), L1(hu))` errors against the analytical profile.
fn run_and_measure(n_cells: usize, t_end: f64) -> (f64, f64) {
    let dx = L_DOMAIN / n_cells as f64;
    let cfl = 0.4;
    let channel = build_channel(n_cells);
    let mut states = analytical_state_at_cells(n_cells);

    let bcs = Boundaries {
        left: Boundary::Discharge { q: Q },
        right: Boundary::Depth {
            h: h_analytical(L_DOMAIN),
        },
    };

    let mut t = 0.0;
    while t < t_end {
        let dt = cfl_time_step(&states, dx, cfl).min(t_end - t);
        forward_euler_step(&mut states, &channel, bcs, dt);
        manning_friction_step(&mut states, N_MANNING, dt, 1e-9);
        t += dt;
    }

    let mut l1_h = 0.0;
    let mut l1_hu = 0.0;
    for (i, s) in states.iter().enumerate() {
        let x_c = (i as f64 + 0.5) * dx;
        let h_exact = h_analytical(x_c);
        l1_h += (s.h - h_exact).abs() * dx;
        l1_hu += (s.hu - Q).abs() * dx;
    }
    (l1_h, l1_hu)
}

#[test]
fn analytical_profile_is_subcritical_everywhere() {
    // For the prescribed (Q, h_base, h_amp), the Froude number must stay
    // below 1 across the whole domain. If it goes critical the BC choices
    // (sub-critical Discharge / Depth) become inconsistent.
    let mut fr_max = 0.0_f64;
    for i in 0..1000 {
        let x = i as f64 * L_DOMAIN / 1000.0;
        let h = h_analytical(x);
        let u = Q / h;
        let fr = u / (G * h).sqrt();
        fr_max = fr_max.max(fr);
    }
    assert!(
        fr_max < 0.9,
        "Fr_max = {fr_max} too close to or above 1 — pick a deeper h_base"
    );
}

#[test]
fn analytical_bed_integration_converges() {
    // Trapezoid integration of dz/dx should agree between a coarse and a
    // fine mesh at common nodes — bounds the integration error so it does
    // not contaminate the solver-level assertions. Trapezoid is O(dx²),
    // so 10× refinement gives ~100× smaller error; the bound below is
    // 1000× the typical z magnitude (≈ 0.1), well above the floor and
    // well below the solver's discretisation error (≈ 0.01).
    let z_coarse = analytical_bed_profile(200);
    let z_fine = analytical_bed_profile(2000);
    for k in 0..=200 {
        let diff = (z_coarse[k] - z_fine[k * 10]).abs();
        assert!(
            diff < 1e-4,
            "bed integration not converged at node {k}: {diff:.2e}"
        );
    }
}

#[test]
fn variable_macdonald_l1_error_under_bound() {
    // n = 200 cells over L = 50 m → dx = 0.25 m. After a few wave
    // transits (t_end = 30 s, transit ≈ 16 s) the system has settled
    // around the discrete steady state, which differs from the analytical
    // by O(dx) (HLL + 1st-order time). The L1 errors below are normalised
    // by the integral of the reference variable across the domain
    // (h_base · L for h, Q · L for hu) so the bounds are dimensionless.
    let (l1_h, l1_hu) = run_and_measure(200, 30.0);
    let rel_l1_h = l1_h / (H_BASE * L_DOMAIN);
    let rel_l1_hu = l1_hu / (Q * L_DOMAIN);
    assert!(
        rel_l1_h < 0.05,
        "relative L1(h) = {rel_l1_h:.5} above 5 % bound"
    );
    assert!(
        rel_l1_hu < 0.05,
        "relative L1(hu) = {rel_l1_hu:.5} above 5 % bound"
    );
}

#[test]
fn variable_macdonald_converges_at_first_order() {
    // 4× refinement should drop L1(h) by a factor of ~4 for a 1st-order
    // FV scheme on a smooth steady state. Loose [2, 6] bounds catch
    // gross order regressions without flaking on the discrete error
    // landscape.
    let (l1_h_coarse, _) = run_and_measure(100, 30.0);
    let (l1_h_fine, _) = run_and_measure(400, 30.0);
    let ratio = l1_h_coarse / l1_h_fine;
    assert!(
        (2.0..=6.0).contains(&ratio),
        "convergence ratio {ratio} not in [2, 6]; coarse L1={l1_h_coarse:.5}, fine L1={l1_h_fine:.5}"
    );
}
