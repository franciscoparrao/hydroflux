//! Completion of the canonical Toro (2001/2009) SWE Riemann-problem
//! suite for the 1D solver. `tests/dam_break.rs` covers the wet-wet
//! Stoker case; this file adds the three configurations the audit
//! flagged as missing, each of which stresses a different branch of
//! the wave-speed logic:
//!
//! 1. **Two-rarefaction with near-dry star region** (diverging flow):
//!    robustness of the Davis estimate when the depression between
//!    the rarefactions almost dries out.
//! 2. **Transcritical (sonic) rarefaction**: the left fan straddles
//!    the sonic point `u − c = 0`; entropy-violating "expansion
//!    shocks" show up exactly here in under-dissipative solvers.
//! 3. **Dry bed on the LEFT** (mirrored Ritter): exercises the
//!    `(dry, wet)` wave-speed branch, the mirror image of the
//!    dry-right case covered in the unit tests.
//!
//! Exact solutions: closed forms for the two-rarefaction and Ritter
//! cases; bisection on the depth function (Toro §5.3 adapted to SWE)
//! for the transcritical dam break.

use hydroflux_solver_1d::{
    Boundaries, Channel1D, Conserved, GRAVITY, cfl_time_step, forward_euler_step,
};
use ndarray::Array1;

const CFL: f64 = 0.4;

fn flat_channel(n: usize, dx: f64) -> Channel1D {
    Channel1D::new(Array1::zeros(n), dx, 0.0)
}

/// March the frictionless solver to `t_end` on a flat bed.
fn run_to(states: &mut [Conserved], channel: &Channel1D, t_end: f64) {
    let mut t = 0.0;
    while t < t_end {
        let dt = cfl_time_step(states, channel.dx, CFL).min(t_end - t);
        assert!(dt.is_finite() && dt > 0.0, "degenerate dt = {dt}");
        forward_euler_step(states, channel, Boundaries::TRANSMISSIVE, dt);
        t += dt;
    }
}

#[test]
fn two_rarefaction_near_dry_star_region() {
    // Symmetric diverging flow: h = 1 on both sides, u_L = −v, u_R = +v.
    // Both waves are rarefactions; the star state has the closed form
    //   c* = (c_L + c_R)/2 − (u_R − u_L)/4,   h* = c*²/g,   u* = 0.
    // v = 5 m/s gives c* ≈ 0.63 m/s → h* ≈ 0.041 m: a 96 % depression,
    // the near-dry configuration that breaks naive wave-speed choices.
    //
    // Resolution note: the star plateau is only ~0.6 m wide at t_end,
    // and first-order HLL dissipation over-deepens it badly on coarse
    // grids (h_mid = 0.012 at dx = 0.05 vs 0.035 at dx = 0.00625,
    // converging monotonically toward h* = 0.041 — checked 2026-07-02).
    // The fine grid is what makes the 30 % band meaningful; the
    // unconditional assertions are positivity and finiteness.
    let n = 3200;
    let dx = 0.00625;
    let l_domain = n as f64 * dx; // 20 m
    let h0 = 1.0;
    let v = 5.0;
    let t_end = 0.5;

    let c0 = (GRAVITY * h0).sqrt();
    let c_star = c0 - v / 2.0;
    assert!(c_star > 0.0, "configuration must not fully dry out");
    let h_star = c_star * c_star / GRAVITY;

    let channel = flat_channel(n, dx);
    let mut states: Vec<Conserved> = (0..n)
        .map(|i| {
            let x = -0.5 * l_domain + (i as f64 + 0.5) * dx;
            if x < 0.0 {
                Conserved::new(h0, -h0 * v)
            } else {
                Conserved::new(h0, h0 * v)
            }
        })
        .collect();
    run_to(&mut states, &channel, t_end);

    // Positivity and finiteness everywhere — the actual robustness claim.
    for (i, s) in states.iter().enumerate() {
        assert!(
            s.h.is_finite() && s.h >= 0.0 && s.hu.is_finite(),
            "cell {i}: h = {}, hu = {}",
            s.h,
            s.hu
        );
    }

    // The star plateau sits between the two fans: sample the centre.
    let mid = n / 2;
    let h_mid = states[mid].h;
    let u_mid = states[mid].hu / states[mid].h.max(1e-12);
    // First-order + near-dry: generous 30 % band on the tiny h*, and
    // |u*| = 0 by symmetry (checked loosely against the star scales).
    assert!(
        (h_mid - h_star).abs() < 0.3 * h_star + 5.0e-3,
        "centre depth {h_mid:.4} vs analytical h* = {h_star:.4}"
    );
    assert!(
        u_mid.abs() < 0.15 * v,
        "centre velocity {u_mid:.3} should vanish by symmetry"
    );
}

/// Exact star depth for a wet-wet dam break at rest via bisection on
/// the SWE depth function (left rarefaction + right shock assumed and
/// verified by the caller's configuration h_L > h_R).
fn stoker_star(h_l: f64, h_r: f64) -> (f64, f64) {
    let c_l = (GRAVITY * h_l).sqrt();
    let f = |h: f64| -> f64 {
        // u* from the left rarefaction invariant:
        let u_from_left = 2.0 * (c_l - (GRAVITY * h).sqrt());
        // u* from the right shock jump:
        let u_from_right = (h - h_r) * (0.5 * GRAVITY * (h + h_r) / (h * h_r)).sqrt();
        u_from_left - u_from_right
    };
    let (mut lo, mut hi) = (h_r * 1.0000001, h_l);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let h_star = 0.5 * (lo + hi);
    let u_star = 2.0 * (c_l - (GRAVITY * h_star).sqrt());
    (h_star, u_star)
}

#[test]
fn transcritical_rarefaction_has_no_expansion_shock() {
    // Dam break h_L = 1, h_R = 0.05 at rest: the left rarefaction
    // spans u − c from −c_L < 0 to u* − c* > 0 — it contains the sonic
    // point, which sits exactly at x = 0 for all t. An
    // entropy-violating solver parks a stationary "expansion shock"
    // there; HLL must resolve a smooth fan through it.
    let n = 800;
    let dx = 0.025;
    let l_domain = n as f64 * dx; // 20 m
    let h_l = 1.0;
    let h_r = 0.05;
    let t_end = 0.8;

    let (h_star, u_star) = stoker_star(h_l, h_r);
    let c_star = (GRAVITY * h_star).sqrt();
    assert!(
        u_star - c_star > 0.0,
        "configuration must be transcritical: u* − c* = {:.4}",
        u_star - c_star
    );

    let channel = flat_channel(n, dx);
    let mut states: Vec<Conserved> = (0..n)
        .map(|i| {
            let x = -0.5 * l_domain + (i as f64 + 0.5) * dx;
            if x < 0.0 {
                Conserved::new(h_l, 0.0)
            } else {
                Conserved::new(h_r, 0.0)
            }
        })
        .collect();
    run_to(&mut states, &channel, t_end);

    // Exact solution inside the fan: h(x,t) = (2c_L − x/t)²/9g between
    // the head (−c_L t) and the tail ((u* − c*) t); star plateau to
    // the shock; undisturbed beyond.
    let c_l = (GRAVITY * h_l).sqrt();
    let shock_speed = u_star * h_star / (h_star - h_r);
    let exact_h = |x: f64| -> f64 {
        if x <= -c_l * t_end {
            h_l
        } else if x <= (u_star - c_star) * t_end {
            let val = (2.0 * c_l - x / t_end) / 3.0;
            val * val / GRAVITY
        } else if x <= shock_speed * t_end {
            h_star
        } else {
            h_r
        }
    };

    // (a) The sonic point x = 0 lies inside the fan: the numerical
    // solution must be locally smooth there — the jump between the
    // two cells straddling x = 0 must be of O(dx·|dh/dx|), not a
    // discontinuity. |dh/dx| at the sonic point is 4c_L/(9g t)·3 —
    // bound the jump by a few cell-slopes.
    let mid = n / 2;
    let jump = (states[mid].h - states[mid - 1].h).abs();
    let fan_slope = 2.0 * (2.0 * c_l) / (9.0 * GRAVITY * t_end) * 3.0; // |dh/dx| scale in the fan
    assert!(
        jump < 5.0 * fan_slope * dx,
        "expansion shock at the sonic point: |Δh| = {jump:.4e} vs fan-slope scale {:.4e}",
        fan_slope * dx
    );

    // (b) L1 error against the exact solution, away from the shock
    // (first-order smearing there is expected).
    let mut l1 = 0.0;
    let mut count = 0usize;
    for (i, s) in states.iter().enumerate() {
        let x = -0.5 * l_domain + (i as f64 + 0.5) * dx;
        if x.abs() < 0.45 * l_domain && (x - shock_speed * t_end).abs() > 8.0 * dx {
            l1 += (s.h - exact_h(x)).abs();
            count += 1;
        }
    }
    l1 /= count as f64;
    assert!(l1 < 8.0e-3, "L1 error vs exact transcritical solution: {l1:.4e}");
}

#[test]
fn ritter_dry_bed_on_the_left() {
    // Mirror of the classical dry-bed dam break: water at rest on the
    // RIGHT of the dam, dry bed on the LEFT. Exact (mirrored Ritter):
    //   h(x,t) = (2c_R + x/t)²/9g   for −2c_R t < x < c_R t,
    // undisturbed h_R beyond, dry before the front. Exercises the
    // (dry, wet) wave-speed branch end-to-end.
    let n = 800;
    let dx = 0.025;
    let l_domain = n as f64 * dx;
    let h_r = 1.0;
    let t_end = 0.6;
    let c_r = (GRAVITY * h_r).sqrt();

    let channel = flat_channel(n, dx);
    let mut states: Vec<Conserved> = (0..n)
        .map(|i| {
            let x = -0.5 * l_domain + (i as f64 + 0.5) * dx;
            if x < 0.0 {
                Conserved::DRY
            } else {
                Conserved::new(h_r, 0.0)
            }
        })
        .collect();
    let m0: f64 = states.iter().map(|s| s.h * dx).sum();
    run_to(&mut states, &channel, t_end);

    let exact_h = |x: f64| -> f64 {
        if x <= -2.0 * c_r * t_end {
            0.0
        } else if x <= c_r * t_end {
            let val = (2.0 * c_r + x / t_end) / 3.0;
            val * val / GRAVITY
        } else {
            h_r
        }
    };

    let mut l1 = 0.0;
    let mut count = 0usize;
    for (i, s) in states.iter().enumerate() {
        assert!(s.h >= 0.0 && s.h.is_finite(), "cell {i}: h = {}", s.h);
        let x = -0.5 * l_domain + (i as f64 + 0.5) * dx;
        if x.abs() < 0.45 * l_domain {
            l1 += (s.h - exact_h(x)).abs();
            count += 1;
        }
    }
    l1 /= count as f64;
    assert!(l1 < 1.0e-2, "L1 error vs mirrored Ritter: {l1:.4e}");

    // Mass conservation on the (transmissive but untouched) domain:
    // nothing has reached the boundaries at t_end, so mass is exact.
    let m1: f64 = states.iter().map(|s| s.h * dx).sum();
    assert!(
        ((m1 - m0) / m0).abs() < 1.0e-10,
        "mass drifted: {m0} → {m1}"
    );

    // Front position: the wet front must track −2c_R t within the
    // first-order smearing (a lagging front is the classical symptom
    // of Davis speeds on dry beds — the bug fixed in this iteration).
    let front_exact = -2.0 * c_r * t_end;
    let front_num = states
        .iter()
        .enumerate()
        .find(|(_, s)| s.h > 1.0e-4)
        .map(|(i, _)| -0.5 * l_domain + (i as f64 + 0.5) * dx)
        .expect("domain cannot be fully dry");
    assert!(
        (front_num - front_exact).abs() < 0.15 * (2.0 * c_r * t_end),
        "front at {front_num:.3} m vs exact {front_exact:.3} m"
    );
}
