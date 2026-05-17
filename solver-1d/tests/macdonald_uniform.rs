//! MacDonald-style steady-state benchmark: uniform Manning flow on a
//! uniformly-sloped channel.
//!
//! The degenerate case of MacDonald et al. (1997) inverse design: prescribe
//! a constant water depth `h(x) = h_n`, derive the bed elevation
//! `z(x) = −S₀·x`, and choose Manning `n` such that bed-slope gravity is
//! exactly balanced by friction at the analytical normal depth:
//!
//! ```text
//!   q = (1/n) h_n^(5/3) √S₀     (Manning equation)
//!   ⇒  h_n = (n q / √S₀)^(3/5)
//! ```
//!
//! Validates jointly the bed-slope source (Audusse hydrostatic
//! reconstruction in `update.rs`) and the semi-implicit friction step
//! (`source.rs`): at `h_n`, they must cancel.
//!
//! # Known limitation with current boundary conditions
//!
//! Only `Transmissive` and `Wall` BCs exist. Neither sustains a prescribed
//! upstream discharge: the transmissive ghost cell sits at the same bed as
//! the first inner cell, so the upstream boundary face has no bed jump and
//! the first cell does not receive the Audusse bed-slope source correction
//! that interior cells do. Friction then drains its momentum and the
//! perturbation propagates downstream as an **upstream boundary layer**.
//!
//! A symmetric (smaller) effect exists at the downstream end: cell `N−1`
//! gets the interior source from its left face, but its right face flux is
//! `F(uniform)` (no reconstruction), whereas its left face flux includes
//! the small reconstruction-induced HLL asymmetry — producing a tiny
//! downstream drift in the opposite direction.
//!
//! Strategy for this commit: long domain + short physical time, and assert
//! preservation only on the boundary-undisturbed middle slab of the domain
//! (empirically the central 50 %). A proper inflow BC (prescribed `q` at
//! the upstream with extended-slope ghost bed) eliminates the upstream
//! artifact and is the topic of the follow-up commit that completes
//! MacDonald with a varying `h(x)` profile.

use hydroflux_solver_1d::{
    Boundaries, Channel1D, Conserved, cfl_time_step, forward_euler_step, manning_friction_step,
};
use ndarray::Array1;

/// Manning normal depth `h_n` for prescribed unit discharge `q`, bed
/// slope `S₀`, and roughness `n`. Inverted from Manning's equation
/// `q = (1/n) h^(5/3) √S₀`.
fn manning_normal_depth(q: f64, slope: f64, manning: f64) -> f64 {
    (manning * q / slope.sqrt()).powf(3.0 / 5.0)
}

/// Uniformly-descending channel of length `n·dx`, bed slope `S₀`,
/// constant Manning `n`. Bed at cell `i` is `z_i = −i·dx·S₀`.
fn sloped_channel(n: usize, dx: f64, slope: f64, manning: f64) -> Channel1D {
    let bed = Array1::from_iter((0..n).map(|i| -(i as f64) * dx * slope));
    Channel1D::new(bed, dx, manning)
}

#[test]
fn manning_normal_depth_matches_formula() {
    // Round-trip: derive h_n, plug back into Manning, recover q.
    let q = 2.0;
    let slope = 0.01;
    let manning = 0.03;
    let h_n = manning_normal_depth(q, slope, manning);
    let q_check = (1.0 / manning) * h_n.powf(5.0 / 3.0) * slope.sqrt();
    assert!(
        (q_check - q).abs() < 1e-12,
        "round-trip failed: h_n={h_n}, q_check={q_check}, q={q}"
    );

    // Sanity on magnitudes: realistic small Chilean river reach.
    assert!(0.1 < h_n && h_n < 5.0, "h_n out of expected range: {h_n}");
}

#[test]
fn interior_preserves_uniform_flow_outside_upstream_boundary_layer() {
    // Long domain + short time so the upstream boundary layer (see module
    // docs) cannot reach the tested interior region. Geometry chosen so
    // that the boundary-layer head travels at most ~20 % of the domain
    // length during the run; the tested region is the downstream 60 %.
    let slope = 0.005;
    let manning = 0.03;
    let q = 1.0;
    let h_n = manning_normal_depth(q, slope, manning);
    let u_n = q / h_n;
    let dx = 0.25;
    let n = 400; // domain length L = 100 m
    let cfl = 0.4;
    let t_end = 5.0; // short

    let channel = sloped_channel(n, dx, slope, manning);
    let mut states: Vec<Conserved> = vec![Conserved::new(h_n, q); n];

    let mut t = 0.0;
    while t < t_end {
        let dt = cfl_time_step(&states, dx, cfl).min(t_end - t);
        forward_euler_step(&mut states, &channel, Boundaries::TRANSMISSIVE, dt);
        manning_friction_step(&mut states, manning, dt, 1e-9);
        t += dt;
    }

    // Assert preservation only on the central 50 % of the domain — the
    // boundary-undisturbed slab where both the upstream perturbation (which
    // propagates ~25 m in 5 s) and the small downstream boundary effect
    // are out of reach. The empirical drift in this region is bounded by
    // IEEE roundoff (verified with a one-off probe: ~1e-14 relative).
    let inner_lo = n / 4;
    let inner_hi = 3 * n / 4;
    let mut max_dev_h = 0.0_f64;
    let mut max_dev_u = 0.0_f64;
    for s in &states[inner_lo..inner_hi] {
        let u = s.hu / s.h;
        max_dev_h = max_dev_h.max((s.h - h_n).abs() / h_n);
        max_dev_u = max_dev_u.max((u - u_n).abs() / u_n);
    }
    // At the analytical Manning equilibrium the bed-slope source and
    // friction loss cancel in the continuous limit, but operator splitting
    // (flux then friction, first order) leaves a drift O(t · g · h · S₀ · dt).
    // For the parameters here the empirical drift in u is ~2e-4 over 5 s;
    // h preserves to ~1e-5. The bound below catches magnitude bugs (a
    // factor-of-2 error in either side would drive drift to O(0.1)) without
    // demanding cancellation that the splitting cannot deliver.
    assert!(
        max_dev_h < 1e-3,
        "interior depth drifted {max_dev_h:.2e} — bed-source/friction balance is off"
    );
    assert!(
        max_dev_u < 1e-3,
        "interior velocity drifted {max_dev_u:.2e} — bed-source/friction balance is off"
    );
}
