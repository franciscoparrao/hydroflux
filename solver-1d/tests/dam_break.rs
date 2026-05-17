//! Wet-wet dam break benchmark against the Stoker (1957) analytical
//! solution. First benchmark with a closed-form reference, intended as a
//! regression guard for the well-balanced FV update + HLL Riemann solver.
//!
//! Problem: flat frictionless channel [0, L], initial dam at x_dam with
//! still water on both sides, h_L > h_R > 0. Released at t = 0. The exact
//! solution is a left rarefaction + constant star region + right shock
//! (Toro 2009, §6.2).
//!
//! Numerical setup: cells with centres x_i = (i + 0.5) · dx, transmissive
//! boundary conditions, CFL-limited forward Euler.

use approx::assert_relative_eq;
use hydroflux_solver_1d::{Boundaries, Channel1D, Conserved, cfl_time_step, forward_euler_step};
use ndarray::Array1;

const G: f64 = 9.81;

/// Star-region depth, velocity, right-shock speed, and the head/tail
/// similarity coordinates of the left rarefaction, for the wet-wet
/// dam break with u_L = u_R = 0 and h_L > h_R > 0.
///
/// Solves f(h*) = f_L(h*, h_L) + f_R(h*, h_R) = 0 by bisection over
/// (h_R, h_L), where the left wave is a rarefaction (h* < h_L) and the
/// right wave is a shock (h* > h_R). See Toro (2009) §6.2 eqs 6.5–6.12.
fn stoker_star_state(h_l: f64, h_r: f64) -> StokerSolution {
    assert!(
        h_l > h_r && h_r > 0.0,
        "wet-wet dam break needs h_L > h_R > 0, got h_L={h_l}, h_R={h_r}"
    );
    let c_l = (G * h_l).sqrt();
    let c_r = (G * h_r).sqrt();

    // f_L for the rarefaction branch (h ≤ h_L).
    let f_l = |h: f64| 2.0 * ((G * h).sqrt() - c_l);
    // f_R for the shock branch (h ≥ h_R).
    let f_r = |h: f64| (h - h_r) * (G * (h + h_r) / (2.0 * h * h_r)).sqrt();
    let f = |h: f64| f_l(h) + f_r(h);

    // Bisection: f is strictly increasing on (h_R, h_L); f(h_R) < 0 < f(h_L).
    let (mut lo, mut hi) = (h_r, h_l);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
        if (hi - lo) < 1e-14 * h_l {
            break;
        }
    }
    let h_star = 0.5 * (lo + hi);
    let c_star = (G * h_star).sqrt();

    // u* = ½(u_L + u_R) + ½(f_R(h*) − f_L(h*))  (Toro eq 6.4); u_L=u_R=0.
    let u_star = 0.5 * (f_r(h_star) - f_l(h_star));

    // Right shock speed (Toro eq 6.12 with u_R = 0).
    let shock_speed = c_r * (h_star * (h_star + h_r) / (2.0 * h_r * h_r)).sqrt();

    // Rarefaction fan extends from head (xi = u_L − c_L = −c_L) to tail
    // (xi = u* − c*).
    let head_rare = -c_l;
    let tail_rare = u_star - c_star;

    StokerSolution {
        h_l,
        h_r,
        c_l,
        h_star,
        u_star,
        shock_speed,
        head_rare,
        tail_rare,
    }
}

#[derive(Debug, Clone, Copy)]
struct StokerSolution {
    h_l: f64,
    h_r: f64,
    c_l: f64,
    h_star: f64,
    u_star: f64,
    shock_speed: f64,
    head_rare: f64,
    tail_rare: f64,
}

impl StokerSolution {
    /// Exact (h, u) at position `x_rel` (relative to the dam) and time `t`.
    fn at(&self, x_rel: f64, t: f64) -> (f64, f64) {
        if t == 0.0 {
            return if x_rel < 0.0 {
                (self.h_l, 0.0)
            } else {
                (self.h_r, 0.0)
            };
        }
        let xi = x_rel / t;
        if xi < self.head_rare {
            (self.h_l, 0.0)
        } else if xi < self.tail_rare {
            // Inside the left rarefaction fan: Riemann invariants
            // u + 2c = u_L + 2 c_L (constant)  and  xi = u − c
            //   ⇒ c = (u_L + 2 c_L − xi) / 3 = (2 c_L − xi) / 3
            //   ⇒ u = (2/3)(c_L + xi)
            let c = (2.0 * self.c_l - xi) / 3.0;
            let u = (2.0 / 3.0) * (self.c_l + xi);
            (c * c / G, u)
        } else if xi < self.shock_speed {
            (self.h_star, self.u_star)
        } else {
            (self.h_r, 0.0)
        }
    }
}

/// Run the solver on a uniform grid and return (L1 error of h, L1 error of hu).
fn run_and_measure(n: usize, h_l: f64, h_r: f64, x_dam: f64, l: f64, t_end: f64) -> (f64, f64) {
    let dx = l / n as f64;
    let cfl = 0.4;
    let channel = Channel1D::new(Array1::zeros(n), dx, 0.0);

    let mut states: Vec<Conserved> = (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) * dx;
            let h = if x < x_dam { h_l } else { h_r };
            Conserved::new(h, 0.0)
        })
        .collect();

    let mut t = 0.0;
    while t < t_end {
        let dt = cfl_time_step(&states, dx, cfl).min(t_end - t);
        forward_euler_step(&mut states, &channel, Boundaries::TRANSMISSIVE, dt);
        t += dt;
    }

    let exact = stoker_star_state(h_l, h_r);
    let mut l1_h = 0.0;
    let mut l1_hu = 0.0;
    for (i, s) in states.iter().enumerate() {
        let x = (i as f64 + 0.5) * dx;
        let (h_exact, u_exact) = exact.at(x - x_dam, t_end);
        l1_h += (s.h - h_exact).abs() * dx;
        l1_hu += (s.hu - h_exact * u_exact).abs() * dx;
    }
    (l1_h, l1_hu)
}

#[test]
fn analytical_solution_sanity() {
    // h_L = 1, h_R = 0.1: classic wet-wet dam break.
    // Sanity-check the star state lies in the right range and matches the
    // implicit Stoker equation to roundoff.
    let s = stoker_star_state(1.0, 0.1);
    assert!(s.h_r < s.h_star && s.h_star < s.h_l);
    assert!(s.u_star > 0.0);
    assert!(s.shock_speed > 0.0);
    assert!(s.head_rare < 0.0);

    // Residual of the Stoker equation at the computed h_star.
    let c_l = (G * s.h_l).sqrt();
    let c_star = (G * s.h_star).sqrt();
    let f_l = 2.0 * (c_star - c_l);
    let f_r = (s.h_star - s.h_r) * (G * (s.h_star + s.h_r) / (2.0 * s.h_star * s.h_r)).sqrt();
    assert_relative_eq!(f_l + f_r, 0.0, epsilon = 1e-12);
    // Star-region velocity matches both wave-jump formulas.
    assert_relative_eq!(s.u_star, -f_l, epsilon = 1e-12);
    assert_relative_eq!(s.u_star, f_r, epsilon = 1e-12);
}

#[test]
fn stoker_wet_wet_l1_error_under_bound() {
    // Classic wet-wet dam break: h_L=1, h_R=0.1, dam at x=0.5, domain [0,1],
    // t_end = 0.075 s. At t_end the shock has travelled ~0.23 m so it sits
    // well inside the domain and transmissive BC is harmless.
    //
    // With HLL + forward Euler at n=400, expected L1 errors are dominated
    // by shock smearing over 3-5 cells. The bound below leaves ~2× margin
    // for safe regression detection.
    let (l1_h, l1_hu) = run_and_measure(400, 1.0, 0.1, 0.5, 1.0, 0.075);
    assert!(l1_h < 0.010, "L1 error in h: {l1_h} (bound 0.010)");
    assert!(l1_hu < 0.020, "L1 error in hu: {l1_hu} (bound 0.020)");
}

#[test]
fn stoker_l1_error_converges_at_first_order() {
    // 4x refinement should reduce L1 error by ~4 for a 1st-order FV scheme
    // on a problem with a shock. We bound the ratio loosely to [2, 6] —
    // tight enough to catch order regressions (e.g. accidentally 0th order),
    // loose enough to avoid flakiness from the discrete error landscape.
    let (l1_h_coarse, _) = run_and_measure(100, 1.0, 0.1, 0.5, 1.0, 0.075);
    let (l1_h_fine, _) = run_and_measure(400, 1.0, 0.1, 0.5, 1.0, 0.075);
    let ratio = l1_h_coarse / l1_h_fine;
    assert!(
        (2.0..=6.0).contains(&ratio),
        "convergence ratio {ratio} not in [2, 6]; coarse L1={l1_h_coarse}, fine L1={l1_h_fine}"
    );
}
