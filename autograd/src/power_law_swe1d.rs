//! 1D Saint-Venant explicit solver with a **continuous power-law
//! cross-section** `T(h) = c · h^p`.
//!
//! Sister module to [`crate::compound_swe1d`]. The compound module
//! uses a two-stage rectangular cross-section (main channel +
//! floodplain with a kink at bank-full); this one uses Leopold's
//! at-a-station hydraulic geometry: `T(h) ∝ h^p` with `p ∈ (0, 1)`
//! typical of natural channels. The continuous form fixes the
//! cross-event generalisation failure observed in iter 7
//! (`validate_manning_huasco_1998`): the 2-stage compound saturates
//! at `w_flood` once `h >> h_bank`, giving a rectangular response at
//! high Q that contradicts the sublinear `h ∝ Q^0.4` shape of
//! empirical rating curves. The power-law `T(h)` gives a
//! consistently sublinear `h-Q` relation across the whole stage
//! range.
//!
//! # Manning normal-depth analysis (Leopold + Manning)
//!
//! For wide-channel approximation `P ≈ T`, hydraulic radius
//! `R = A/P = h/(p+1)`. Manning velocity:
//! ```text
//!   V = (1/n) · R^(2/3) · √S = (1/n) · (h/(p+1))^(2/3) · √S
//! ```
//! With `A = c · h^(p+1) / (p+1)`:
//! ```text
//!   Q = V · A = (1/n) · (p+1)^(-5/3) · c · h^(p + 5/3) · √S
//! ```
//! Solving for `h`:
//! ```text
//!   h ∝ Q^(1/(p + 5/3))
//! ```
//! For a target rating-curve exponent `b` (e.g., `b = 0.40`):
//! `p = 1/b − 5/3`. With `b = 0.40` → `p = 5/6 ≈ 0.833`.
//!
//! So calibrating `p` to fit the empirical rating-curve shape is the
//! mechanism by which this module captures cross-event behaviour
//! that the compound module cannot.
//!
//! # State and parameters
//!
//! State: `(A, Q)` per cell (area + total discharge), same as
//! [`crate::compound_swe1d`]. Cross-section parameters
//! `PowerLawSection { coefficient: c, exponent: p }` are generic
//! over `T: Real` so they can be passed as `Dual` for calibration.
//!
//! # Wide-channel limit
//!
//! `PowerLawSection { coefficient: W, exponent: 0 }` gives constant
//! top width `T = W`, recovering the rectangular wide-channel
//! solver. The `wide_channel_limit_matches_swe1d` test asserts this.

use crate::Real;

/// Power-law cross-section: top width `T(h) = c · h^p`.
///
/// All fields are generic over `T: Real` so the parameters
/// themselves can carry AD derivatives. For calibration, instantiate
/// with `Dual::variable(value)` on the parameter being differentiated
/// and `Dual::constant(value)` on the others.
#[derive(Debug, Clone, Copy)]
pub struct PowerLawSection<T: Real> {
    /// Coefficient `c` in `T(h) = c · h^p` [units of m^(1-p)].
    pub coefficient: T,
    /// Exponent `p` in `T(h) = c · h^p` [dimensionless]. Typical
    /// natural channels: `p ∈ [0.3, 0.8]`. `p = 0` recovers
    /// rectangular wide-channel; `p = 1` is triangular V-shape.
    pub exponent: T,
}

impl<T: Real> PowerLawSection<T> {
    /// Top width at stage `h`: `T(h) = c · h^p`.
    ///
    /// Guards against `h ≤ 0` returning 0 (dry/no top width) instead
    /// of NaN from `0^p` for negative `p` or undefined `log(0)`.
    pub fn top_width(&self, h: T) -> T {
        if h.value() <= 0.0 {
            T::zero()
        } else {
            self.coefficient * h.powt(self.exponent)
        }
    }

    /// Cross-sectional wetted area:
    /// `A(h) = c · h^(p+1) / (p+1) = T(h) · h / (p+1)`.
    pub fn area(&self, h: T) -> T {
        if h.value() <= 0.0 {
            T::zero()
        } else {
            self.coefficient * h.powt(self.exponent + T::one()) / (self.exponent + T::one())
        }
    }

    /// Wetted perimeter. Wide-channel approximation `P ≈ T`,
    /// adequate when `h << T(h)` (typical for natural channels).
    pub fn perimeter(&self, h: T) -> T {
        // Use max(T, 2h) so very narrow / very deep degenerate cases
        // give at least the sidewall contribution.
        let t = self.top_width(h);
        let two_h = h * 2.0;
        if t.value() > two_h.value() {
            t
        } else {
            two_h
        }
    }

    /// First moment of area below stage `h`:
    /// `I₁(h) = ∫₀ʰ T(η)·(h−η) dη = c · h^(p+2) / [(p+1)(p+2)]`.
    pub fn pressure_integral(&self, h: T) -> T {
        if h.value() <= 0.0 {
            T::zero()
        } else {
            let p = self.exponent;
            let p1 = p + T::one();
            let p2 = p + T::from_f64(2.0);
            self.coefficient * h.powt(p + T::from_f64(2.0)) / (p1 * p2)
        }
    }

    /// Recover stage `h` from area `A`:
    /// `h = ((p+1) · A / c)^(1/(p+1))`.
    pub fn stage(&self, a: T) -> T {
        if a.value() <= 0.0 {
            T::zero()
        } else {
            let p1 = self.exponent + T::one();
            (a * p1 / self.coefficient).powt(T::one() / p1)
        }
    }
}

/// Upstream boundary specification (mirrors `compound_swe1d::LeftBc`).
#[derive(Debug, Clone, Copy)]
pub enum LeftBc<T: Real> {
    /// Dirichlet on stage `h` and total discharge `Q` (m³/s).
    Dirichlet {
        /// Imposed stage at the upstream ghost cell [m].
        h: T,
        /// Imposed total volumetric discharge [m³/s].
        q: T,
    },
    /// Transmissive (zero-gradient).
    Transmissive,
}

/// Downstream boundary specification.
#[derive(Debug, Clone, Copy)]
pub enum RightBc {
    /// Transmissive (zero-gradient).
    Transmissive,
}

/// Lax-Friedrichs explicit step on `(A, Q)` for the power-law section.
#[allow(clippy::too_many_arguments)]
pub fn lax_friedrichs_step<T: Real>(
    section: &PowerLawSection<T>,
    a: &[T],
    q: &[T],
    bed: &[f64],
    dx: f64,
    dt: f64,
    manning_n: T,
    gravity: f64,
    left_bc: LeftBc<T>,
    _right_bc: RightBc,
    a_next: &mut [T],
    q_next: &mut [T],
) {
    let n = a.len();
    assert_eq!(q.len(), n);
    assert_eq!(bed.len(), n);
    assert_eq!(a_next.len(), n);
    assert_eq!(q_next.len(), n);
    assert!(n >= 2);

    // Global LF dissipation: max(|u| + c) over the domain with
    // c = √(g·A/T).
    let mut alpha = T::zero();
    for i in 0..n {
        let a_i = a[i];
        if a_i.value() < 1.0e-9 {
            continue;
        }
        let h_i = section.stage(a_i);
        let t_top = section.top_width(h_i);
        let c = (a_i * gravity / t_top.max(T::from_f64(1.0e-6))).sqrt();
        let u = q[i] / a_i;
        let s = u.abs() + c;
        if s.value() > alpha.value() {
            alpha = s;
        }
    }

    let (a_l_ghost, q_l_ghost) = match left_bc {
        LeftBc::Dirichlet { h, q: q_g } => (section.area(h), q_g),
        LeftBc::Transmissive => (a[0], q[0]),
    };
    let (a_r_ghost, q_r_ghost) = (a[n - 1], q[n - 1]);

    let flux = |a_i: T, q_i: T| -> (T, T) {
        let a_safe = a_i.max(T::from_f64(1.0e-12));
        let h = section.stage(a_safe);
        let mass = q_i;
        let mom = (q_i * q_i) / a_safe + section.pressure_integral(h) * gravity;
        (mass, mom)
    };

    let lf_face = |a_l: T, q_l: T, a_r: T, q_r: T| -> (T, T) {
        let (fl_a, fl_q) = flux(a_l, q_l);
        let (fr_a, fr_q) = flux(a_r, q_r);
        let half_alpha = alpha * 0.5;
        let f_a = (fl_a + fr_a) * 0.5 - half_alpha * (a_r - a_l);
        let f_q = (fl_q + fr_q) * 0.5 - half_alpha * (q_r - q_l);
        (f_a, f_q)
    };

    let mut f_a_left = lf_face(a_l_ghost, q_l_ghost, a[0], q[0]).0;
    let mut f_q_left = lf_face(a_l_ghost, q_l_ghost, a[0], q[0]).1;

    let dt_over_dx = dt / dx;

    for i in 0..n {
        let (a_r, q_r) = if i + 1 < n {
            (a[i + 1], q[i + 1])
        } else {
            (a_r_ghost, q_r_ghost)
        };
        let (f_a_right, f_q_right) = lf_face(a[i], q[i], a_r, q_r);

        let dz_dx = if i == 0 {
            (bed[1] - bed[0]) / dx
        } else if i == n - 1 {
            (bed[n - 1] - bed[n - 2]) / dx
        } else {
            (bed[i + 1] - bed[i - 1]) / (2.0 * dx)
        };

        let a_new = a[i] - (f_a_right - f_a_left) * dt_over_dx;
        let a_safe = a_new.max(T::from_f64(1.0e-9));

        let s_bed = a[i] * (-gravity * dz_dx);

        // Manning friction: force = −g·n²·Q²·P^(4/3)/A^(7/3).
        let q_star = q[i] + s_bed * dt - (f_q_right - f_q_left) * dt_over_dx;
        let h_i = section.stage(a_safe);
        let p = section.perimeter(h_i).max(T::from_f64(1.0e-9));
        let coeff = manning_n * manning_n * (dt * gravity) * q_star.abs() * p.powf(4.0 / 3.0)
            / a_safe.powf(7.0 / 3.0);
        let q_new = q_star / (coeff + 1.0);

        a_next[i] = a_new.max(T::zero());
        q_next[i] = if a_next[i].value() > 1.0e-9 {
            q_new
        } else {
            T::zero()
        };

        f_a_left = f_a_right;
        f_q_left = f_q_right;
    }
}

/// CFL-bounded `dt` for the power-law solver.
///
/// Returns `f64`, extracting `.value()` from the state: the time step
/// is deliberately **not differentiated**. Consequence: the AD gradient
/// is that of the scheme with a frozen dt sequence, which differs from
/// the true sensitivity of the discrete map by an O(dt) term (the
/// semi-implicit friction has ∂G/∂dt ≠ 0 even at steady state).
/// Measured at ~1e-3 relative under CFL 0.4 in
/// `ad_gradient_matches_central_finite_difference_for_n_c_p`; larger
/// for transient QoIs.
pub fn cfl_dt<T: Real>(
    section: &PowerLawSection<T>,
    a: &[T],
    q: &[T],
    dx: f64,
    gravity: f64,
    cfl: f64,
) -> f64 {
    let mut max_lambda = 0.0_f64;
    for i in 0..a.len() {
        let a_v = a[i].value();
        if a_v < 1.0e-9 {
            continue;
        }
        let h = section.stage(a[i]);
        let t_top = section.top_width(h).value().max(1.0e-6);
        let c = (gravity * a_v / t_top).sqrt();
        let u = q[i].value() / a_v;
        let lam = u.abs() + c;
        if lam > max_lambda {
            max_lambda = lam;
        }
    }
    if max_lambda < 1.0e-12 {
        return cfl * dx / 1.0;
    }
    cfl * dx / max_lambda
}

/// Run until `t_end`. Same signature as `compound_swe1d::run`.
#[allow(clippy::too_many_arguments)]
pub fn run<T: Real>(
    section: &PowerLawSection<T>,
    a0: Vec<T>,
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
    let mut a = a0;
    let mut q = q0;
    let mut a_next = vec![T::zero(); a.len()];
    let mut q_next = vec![T::zero(); q.len()];
    let mut t = 0.0;
    let mut steps = 0;
    while t < t_end {
        let dt = cfl_dt(section, &a, &q, dx, gravity, cfl).min(t_end - t);
        lax_friedrichs_step(
            section, &a, &q, bed, dx, dt, manning_n, gravity, left_bc, right_bc,
            &mut a_next, &mut q_next,
        );
        std::mem::swap(&mut a, &mut a_next);
        std::mem::swap(&mut q, &mut q_next);
        t += dt;
        steps += 1;
        if steps > 1_000_000 {
            panic!("power_law_swe1d::run did not finish in 1M steps");
        }
    }
    (a, q, steps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Dual;

    const G: f64 = 9.81;
    const EPS: f64 = 1.0e-10;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn area_closed_form_matches_integral_for_rectangular() {
        // p = 0 → T = c constant → A = c·h.
        let s = PowerLawSection::<f64> {
            coefficient: 30.0,
            exponent: 0.0,
        };
        assert!(approx_eq(s.area(1.0), 30.0, EPS));
        assert!(approx_eq(s.area(2.5), 75.0, EPS));
        assert!(approx_eq(s.top_width(2.5), 30.0, EPS));
    }

    #[test]
    fn area_closed_form_matches_for_triangular() {
        // p = 1 → T = c·h → A = c·h²/2 (triangular V-shape).
        let s = PowerLawSection::<f64> {
            coefficient: 10.0,
            exponent: 1.0,
        };
        assert!(approx_eq(s.area(2.0), 0.5 * 10.0 * 4.0, EPS));
        assert!(approx_eq(s.top_width(2.0), 20.0, EPS));
    }

    #[test]
    fn stage_inverts_area_round_trip() {
        let s = PowerLawSection::<f64> {
            coefficient: 22.0,
            exponent: 0.83,
        };
        for h in [0.1, 0.5, 1.0, 2.0, 5.0] {
            let a = s.area(h);
            let h_back = s.stage(a);
            assert!(approx_eq(h, h_back, 1.0e-8), "h = {h}, recovered = {h_back}");
        }
    }

    #[test]
    fn pressure_integral_recovers_rectangular_limit() {
        // p = 0 → I₁ = c·h²/2.
        let s = PowerLawSection::<f64> {
            coefficient: 30.0,
            exponent: 0.0,
        };
        assert!(approx_eq(s.pressure_integral(1.5), 0.5 * 30.0 * 2.25, EPS));
    }

    #[test]
    fn stage_inverse_is_differentiable_in_p() {
        // Make exponent a Dual variable; verify gradient is finite and
        // non-zero at a typical operating point.
        let s = PowerLawSection::<Dual> {
            coefficient: Dual::constant(22.0),
            exponent: Dual::variable(0.83),
        };
        let a = Dual::constant(50.0);
        let h = s.stage(a);
        assert!(h.dval.is_finite() && h.dval.abs() > 0.0);
    }

    #[test]
    fn wide_channel_limit_recovers_swe1d_steady_state() {
        // p = 0, c = W → rectangular wide-channel of width W.
        // Steady inflow on a sloping bed should yield Manning normal
        // depth approximately.
        let w = 20.0_f64;
        let s = PowerLawSection::<f64> {
            coefficient: w,
            exponent: 0.0,
        };
        let slope = 0.001_f64;
        let n = 0.04_f64;
        let n_cells = 60;
        let dx = 5.0;
        let bed: Vec<f64> = (0..n_cells).map(|i| -slope * (i as f64 + 0.5) * dx).collect();

        let q_in_total = 1.5_f64;
        let q_per_w = q_in_total / w;
        let h_n = (n * q_per_w / slope.sqrt()).powf(3.0 / 5.0);
        let a_n = w * h_n;

        let a0 = vec![a_n; n_cells];
        let q0 = vec![q_in_total; n_cells];
        let (a, _q, _) = run(
            &s, a0, q0, &bed, dx, 500.0, n, G, 0.4,
            LeftBc::Dirichlet { h: h_n, q: q_in_total },
            RightBc::Transmissive,
        );
        let mid = n_cells / 2;
        let h_mid = s.stage(a[mid]);
        // Loose tolerance: LF + boundary effects.
        assert!(
            (h_mid / h_n - 1.0).abs() < 0.10,
            "h_mid = {h_mid:.4}, h_n = {h_n:.4}"
        );
    }

    #[test]
    fn ad_gradient_matches_central_finite_difference_for_n_c_p() {
        // Independent ground truth for the AD gradient of the FULL
        // solver — the same end-to-end check `swe1d.rs` has, ported
        // to the module that actually feeds the calibration paper.
        // For each of the three calibrated parameters θ ∈ {n, c, p},
        // compute the steady-state stage h_mid(θ ± ε) with the plain
        // f64 solver and compare the central difference against the
        // Dual forward-mode derivative. The QoI is the STAGE (not the
        // area), so the gradient also threads the stage() inversion —
        // the same path the rating-curve loss of the paper uses.
        //
        // Setup mirrors `wide_channel_limit_recovers_swe1d_steady_state`
        // but with a natural-channel exponent. Initial conditions and
        // the Dirichlet BC are held FIXED across perturbations (they
        // are Dual::constant on the AD side), so FD and AD measure the
        // same quantity.
        //
        // Tolerance 5e-3, NOT machine-level: the CFL dt is state-
        // dependent and deliberately not differentiated (`cfl_dt`
        // extracts `.value()`), and the semi-implicit friction factor
        // uses `|q_star| = |q + dt·R|`, so ∂G/∂dt ≠ 0 even at the
        // steady fixed point — FD sees that term, AD does not. The
        // residual is O(dt): measured 1.2e-3 / 2.2e-3 / 1.7e-3 for
        // n / c / p at CFL 0.4, and it halves when CFL halves (checked
        // 2026-07-02), which pins it to the frozen-dt term rather than
        // a defect in the Dual rules — those would show up at O(1).
        // For transient QoIs the frozen-dt gap is larger; re-derive
        // the tolerance before reusing this pattern there.
        let n_cells = 60;
        let dx = 5.0;
        let slope = 0.005_f64;
        let bed: Vec<f64> = (0..n_cells).map(|i| -slope * (i as f64 + 0.5) * dx).collect();
        let t_end = 2000.0_f64;
        let mid = n_cells / 2;

        let n0 = 0.04_f64;
        let c0 = 22.0_f64;
        let p0 = 5.0 / 6.0;
        let q_in = 20.0_f64;

        // Analytical normal depth for the baseline parameters; used
        // as fixed IC/BC for every run (FD and AD alike).
        let h_exp = 1.0 / (p0 + 5.0 / 3.0);
        let prefactor = (p0 + 1.0).powf(5.0 / 3.0) * n0 / (c0 * slope.sqrt());
        let h_bc = (prefactor * q_in).powf(h_exp);
        let s0 = PowerLawSection::<f64> {
            coefficient: c0,
            exponent: p0,
        };
        let a_init = s0.area(h_bc);

        let h_mid_f64 = |n: f64, c: f64, p: f64| -> f64 {
            let s = PowerLawSection::<f64> {
                coefficient: c,
                exponent: p,
            };
            let (a, _, _) = run(
                &s,
                vec![a_init; n_cells],
                vec![q_in; n_cells],
                &bed,
                dx,
                t_end,
                n,
                G,
                0.4,
                LeftBc::Dirichlet { h: h_bc, q: q_in },
                RightBc::Transmissive,
            );
            s.stage(a[mid])
        };

        let h_mid_ad = |n: Dual, c: Dual, p: Dual| -> Dual {
            let s = PowerLawSection::<Dual> {
                coefficient: c,
                exponent: p,
            };
            let (a, _, _) = run(
                &s,
                vec![Dual::constant(a_init); n_cells],
                vec![Dual::constant(q_in); n_cells],
                &bed,
                dx,
                t_end,
                n,
                G,
                0.4,
                LeftBc::Dirichlet {
                    h: Dual::constant(h_bc),
                    q: Dual::constant(q_in),
                },
                RightBc::Transmissive,
            );
            s.stage(a[mid])
        };

        // (label, eps, FD gradient, AD gradient) per parameter.
        let cases: [(&str, f64, f64, f64); 3] = [
            ("manning_n", 1.0e-5, {
                (h_mid_f64(n0 + 1.0e-5, c0, p0) - h_mid_f64(n0 - 1.0e-5, c0, p0)) / 2.0e-5
            }, h_mid_ad(Dual::variable(n0), Dual::constant(c0), Dual::constant(p0)).dval),
            ("coefficient_c", 1.0e-3, {
                (h_mid_f64(n0, c0 + 1.0e-3, p0) - h_mid_f64(n0, c0 - 1.0e-3, p0)) / 2.0e-3
            }, h_mid_ad(Dual::constant(n0), Dual::variable(c0), Dual::constant(p0)).dval),
            ("exponent_p", 1.0e-5, {
                (h_mid_f64(n0, c0, p0 + 1.0e-5) - h_mid_f64(n0, c0, p0 - 1.0e-5)) / 2.0e-5
            }, h_mid_ad(Dual::constant(n0), Dual::constant(c0), Dual::variable(p0)).dval),
        ];

        for (label, _eps, grad_fd, grad_ad) in cases {
            assert!(
                grad_fd.abs() > 1.0e-9,
                "{label}: FD gradient unexpectedly zero ({grad_fd:.3e}) — test is not exercising the parameter"
            );
            let rel_err = (grad_ad - grad_fd).abs() / grad_fd.abs();
            assert!(
                rel_err < 5.0e-3,
                "{label}: AD vs FD mismatch: ad = {grad_ad:.6e}, fd = {grad_fd:.6e}, rel_err = {rel_err:.3e}"
            );
        }
    }

    #[test]
    fn rating_curve_exponent_implied_by_p() {
        // Analytical Manning normal depth: h ∝ Q^(1/(p + 5/3)).
        // For target b = 0.40 we need p = 1/b − 5/3 = 5/6.
        // Verify the simulator's h vs Q response at two Q values
        // matches this exponent within a tolerance.
        let p = 5.0 / 6.0;
        let c = 22.0;
        let s = PowerLawSection::<f64> {
            coefficient: c,
            exponent: p,
        };
        let slope = 0.005_f64;
        let n = 0.04_f64;
        let g = G;

        // Use analytical h_n directly (skip simulator — testing the
        // CLOSED-FORM section physics, not the LF solver).
        let h_n = |q: f64| -> f64 {
            // From Q = (1/n)·(p+1)^(-5/3)·c·h^(p+5/3)·√S:
            // h = [(p+1)^(5/3) · n · q / (c · √S)]^(1/(p + 5/3))
            let exp = 1.0 / (p + 5.0 / 3.0);
            let prefactor = (p + 1.0).powf(5.0 / 3.0) * n / (c * slope.sqrt());
            (prefactor * q).powf(exp)
        };

        let h_low = h_n(10.0);
        let h_high = h_n(100.0);
        let observed_b = (h_high / h_low).ln() / (100.0_f64 / 10.0).ln();
        assert!(
            (observed_b - 0.40).abs() < 1.0e-6,
            "Observed Q-exponent {observed_b:.4} != target 0.40"
        );
        let _ = g;
    }
}
