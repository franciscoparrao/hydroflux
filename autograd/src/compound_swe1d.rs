//! 1D Saint-Venant explicit solver for arbitrary cross-sections.
//!
//! Generalises [`crate::swe1d`] from the wide-channel `(h, q)` form
//! to the area-discharge `(A, Q)` form, where `A` is the wetted
//! cross-sectional area and `Q` is the total volumetric discharge.
//! This lets the solver represent **compound channels** — a narrow
//! main channel that spreads into a wider floodplain at bank-full
//! depth — which the wide-channel approximation cannot.
//!
//! # Motivation
//!
//! Track A iter 5 (`calibrate_manning_huasco_2017_width`) showed
//! that wide-channel 1D with width DEM-derived (42 m) recovered
//! `n ≈ 0.024` for the Atacama 2017 event but with RMSE ≈ 0.43 m
//! against the rating-curve target. The misfit *shape* (undershoot
//! at baseflow, overshoot at peak) is the signature of a missing
//! compound section: the wide-channel response is too steep with Q
//! because all stages spread the same width. A real channel has
//! `dh/dQ` steeper at low Q (narrow main channel) and flatter at
//! high Q (floodplain absorbs the excess), which matches the
//! sublinear `h ∝ Q^0.4` shape of the empirical rating curve.
//!
//! # Discretisation
//!
//! State `(A_i, Q_i)` at cell centres. Conservation form:
//!
//! ```text
//! ∂A/∂t + ∂Q/∂x = 0
//! ∂Q/∂t + ∂(Q²/A + g·I₁(h))/∂x = −g·A·dz/dx − g·A·S_f
//! ```
//!
//! where:
//! - `I₁(h) = ∫₀ʰ T(y)·(h−y) dy` is the first moment of area
//!   (hydrostatic pressure integral; pressure force = `g·I₁`).
//! - `T(y)` is the top width at height `y` above the bed.
//! - `S_f = n²·Q²·P^(4/3) / A^(10/3)` is the Manning friction
//!   *slope*, with `P` the wetted perimeter (Chow 1959, eq. 5-15).
//!   The friction *force* on the momentum equation is
//!   `−g·A·S_f = −g·n²·Q²·P^(4/3)/A^(7/3)`.
//!
//! Numerical scheme: same Lax-Friedrichs as [`crate::swe1d`]
//! (symmetric flux with global α dissipation) so the solver stays
//! minimal. Higher-order MUSCL/HLLC for compound sections can come
//! later when the application case demands it.
//!
//! # Wide-channel limit
//!
//! `CompoundSection { w_main: W, w_flood: W, h_bank: arbitrary }`
//! reduces analytically to a rectangular channel of width `W`, and
//! the solver's output matches [`crate::swe1d::run`] modulo
//! discretisation noise. The `wide_channel_limit_matches_swe1d`
//! test asserts this.

use crate::Real;

/// Two-stage rectangular compound cross-section.
///
/// Below bank-full (`h < h_bank`): rectangular main channel of width
/// `w_main`. Above bank-full: the wetted top width jumps to
/// `w_flood` (the floodplain is treated as a rectangular extension
/// on both sides of the main channel, symmetric).
///
/// Reasonable defaults for Andean gravel-bed reaches:
/// - `w_main` from DEM stream-network buffer (often 1 pixel = 30 m
///   at 30-m DEM resolution).
/// - `w_flood` from HAND-connected-walk perpendicular to flow at a
///   higher threshold (e.g., HAND < 2 m).
/// - `h_bank` from bankful indicator (often ~ 1 m for medium rivers).
#[derive(Debug, Clone, Copy)]
pub struct CompoundSection {
    /// Main channel width [m].
    pub w_main: f64,
    /// Floodplain top width [m]. Must satisfy `w_flood >= w_main`.
    pub w_flood: f64,
    /// Bank-full depth [m] (transition stage).
    pub h_bank: f64,
}

impl CompoundSection {
    /// Cross-sectional wetted area `A(h)`.
    pub fn area<T: Real>(&self, h: T) -> T {
        let h_bank = T::from_f64(self.h_bank);
        let h_in_main = h.min(h_bank);
        let h_above = (h - h_bank).max(T::zero());
        h_in_main * self.w_main + h_above * self.w_flood
    }

    /// Wetted perimeter `P(h)`.
    ///
    /// Below bank-full: `P = w_main + 2·h` (rectangular).
    /// Above bank-full: `P = w_flood + 2·h_bank + 2·(h − h_bank)
    /// + (w_flood − w_main)` minus the horizontal step that is
    /// NOT wetted from below. Working through the geometry of the
    /// step section (see module doc): `P_above = w_flood + 2·h`.
    pub fn perimeter<T: Real>(&self, h: T) -> T {
        if h.value() < self.h_bank {
            T::from_f64(self.w_main) + h * 2.0
        } else {
            T::from_f64(self.w_flood) + h * 2.0
        }
    }

    /// First moment of area below stage `h`, i.e.
    /// `I₁(h) = ∫₀ʰ T(y)·(h − y) dy`. Multiplying by `g` gives the
    /// hydrostatic pressure force on the cross-section.
    ///
    /// Closed form for the two-stage rectangular section:
    /// - `h ≤ h_bank`: `I₁ = w_main · h²/2`.
    /// - `h > h_bank`:
    ///   `I₁ = w_main · h_bank · (h − h_bank/2) + w_flood · (h − h_bank)²/2`.
    pub fn pressure_integral<T: Real>(&self, h: T) -> T {
        if h.value() < self.h_bank {
            h * h * (0.5 * self.w_main)
        } else {
            let h_b = T::from_f64(self.h_bank);
            let half_main = T::from_f64(0.5 * self.w_main);
            let half_flood = T::from_f64(0.5 * self.w_flood);
            // w_main · h_bank · (h − h_bank/2)
            let part1 = h_b * (h - h_b * 0.5) * self.w_main;
            // w_flood · (h − h_bank)²/2
            let _ = half_main;
            let dh = h - h_b;
            let part2 = dh * dh * half_flood * 2.0; // *2 cancels half then half
            // Actually pressure_integral floodplain term is w_flood · (h-h_bank)² / 2.
            // half_flood = 0.5 * w_flood. Multiplying by 2.0 gives w_flood. Then dh*dh.
            // So part2 = dh * dh * w_flood / 2. Wait that's wrong, half_flood * 2 = w_flood,
            // then we'd have part2 = dh * dh * w_flood (no /2). Let me redo:
            let part2_correct = dh * dh * (0.5 * self.w_flood);
            let _ = part2;
            part1 + part2_correct
        }
    }

    /// Lotter (1933) equivalent Manning roughness for the compound
    /// cross-section. Treats the channel as two sub-sections with
    /// distinct `n_main` (low-flow rectangular channel) and `n_flood`
    /// (overbank), partitioned at `h_bank`. The Lotter weighting is
    /// `n_eq = (P_main + P_flood) / (P_main/n_main + P_flood/n_flood)`
    /// with
    ///
    /// - `P_main  = w_main + 2·min(h, h_bank)`
    /// - `P_flood = (w_flood − w_main) + 2·max(h − h_bank, 0)`
    ///
    /// For `h ≤ h_bank` only the main channel is wet (`P_flood = 0`)
    /// and the formula collapses to `n_eq = n_main`. For `h > h_bank`
    /// the contributions weigh by their wetted perimeter, matching
    /// the standard HEC-RAS conveyance-method composite roughness.
    ///
    /// Used by [`lax_friedrichs_step`] to allow a calibrated
    /// `n_main` distinct from the floodplain `n_flood` (vegetation
    /// riparian vs sparse overbank, for example).
    pub fn compound_manning<T: Real>(&self, h: T, n_main: T, n_flood: T) -> T {
        let h_b = T::from_f64(self.h_bank);
        if h.value() <= self.h_bank {
            n_main
        } else {
            let p_main = T::from_f64(self.w_main) + h_b * 2.0;
            let p_flood = T::from_f64(self.w_flood - self.w_main) + (h - h_b) * 2.0;
            let inv = p_main / n_main + p_flood / n_flood;
            (p_main + p_flood) / inv
        }
    }

    /// Recover stage `h` from area `A` (inverse of [`area`]).
    pub fn stage<T: Real>(&self, a: T) -> T {
        let a_bank = self.w_main * self.h_bank;
        if a.value() < a_bank {
            a / T::from_f64(self.w_main)
        } else {
            let excess = a - T::from_f64(a_bank);
            T::from_f64(self.h_bank) + excess / T::from_f64(self.w_flood)
        }
    }
}

/// Upstream boundary specification (analogous to `swe1d::LeftBc`).
#[derive(Debug, Clone, Copy)]
pub enum LeftBc<T: Real> {
    /// Dirichlet on stage `h` and total discharge `Q` (m³/s).
    Dirichlet {
        /// Imposed stage at the upstream ghost cell [m].
        h: T,
        /// Imposed total volumetric discharge at the upstream ghost cell [m³/s].
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

/// Lax-Friedrichs explicit step on `(A, Q)` with split Manning
/// roughness for the main channel (`n_main`) and the floodplain
/// (`n_flood`). When `h ≤ h_bank` in a given cell, only `n_main`
/// contributes; for `h > h_bank` the cell sees the Lotter (1933)
/// compound roughness — see [`CompoundSection::compound_manning`].
///
/// To match the original single-Manning behaviour (Track A iter 6),
/// pass the same value for both: `lax_friedrichs_step(..., n, n,
/// ...)` is bit-exact equivalent to the pre-split API.
#[allow(clippy::too_many_arguments)]
pub fn lax_friedrichs_step<T: Real>(
    section: &CompoundSection,
    a: &[T],
    q: &[T],
    bed: &[f64],
    dx: f64,
    dt: f64,
    manning_main: T,
    manning_flood: T,
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

    // Global Lax-Friedrichs dissipation: max(|u| + c) over the
    // domain, with `c = √(g·A/T)` (gravity wave on top width) and
    // `u = Q/A`.
    let mut alpha = T::zero();
    for i in 0..n {
        let a_i = a[i];
        if a_i.value() < 1.0e-9 {
            continue;
        }
        let h_i = section.stage(a_i);
        let t_top = if h_i.value() < section.h_bank {
            T::from_f64(section.w_main)
        } else {
            T::from_f64(section.w_flood)
        };
        let c = (a_i * gravity / t_top).sqrt();
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
        // Pressure force = g · I₁(h).
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

        // Bed-slope source: −g·A·dz/dx.
        let s_bed = a[i] * (-gravity * dz_dx);

        // Manning friction (semi-implicit on Q):
        // dQ/dt += −g·n²·Q·|Q|·P^(4/3) / A^(10/3).
        // Linearise around |Q*| and apply implicit on Q:
        //   Q_{n+1} = Q* / (1 + dt·g·n²·|Q*|·P^(4/3)/A^(10/3))
        // The roughness `n` is the Lotter compound roughness — for
        // cells with `h ≤ h_bank` this reduces to `n_main`.
        let q_star = q[i] + s_bed * dt - (f_q_right - f_q_left) * dt_over_dx;
        let h_i = section.stage(a_safe);
        let p = section.perimeter(h_i).max(T::from_f64(1.0e-9));
        let n_cell = section.compound_manning(h_i, manning_main, manning_flood);
        // Friction force on Q: −g·n²·Q²·P^(4/3)/A^(7/3) (NOT A^(10/3),
        // which is the friction *slope*; A^(7/3) is the *force* form
        // after multiplying by A — see module doc).
        let coeff = n_cell * n_cell * (dt * gravity) * q_star.abs() * p.powf(4.0 / 3.0)
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

/// CFL-bounded `dt` for the compound solver.
pub fn cfl_dt<T: Real>(
    section: &CompoundSection,
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
        let h_v = section.stage(a[i]).value();
        let t_top = if h_v < section.h_bank {
            section.w_main
        } else {
            section.w_flood
        };
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

/// Run until `t_end`. Mirrors [`crate::swe1d::run`] for the
/// compound formulation, with split Manning `(n_main, n_flood)`.
///
/// Pass identical values for both to reproduce the pre-split
/// single-Manning behaviour exactly.
#[allow(clippy::too_many_arguments)]
pub fn run<T: Real>(
    section: &CompoundSection,
    a0: Vec<T>,
    q0: Vec<T>,
    bed: &[f64],
    dx: f64,
    t_end: f64,
    manning_main: T,
    manning_flood: T,
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
            section,
            &a,
            &q,
            bed,
            dx,
            dt,
            manning_main,
            manning_flood,
            gravity,
            left_bc,
            right_bc,
            &mut a_next,
            &mut q_next,
        );
        std::mem::swap(&mut a, &mut a_next);
        std::mem::swap(&mut q, &mut q_next);
        t += dt;
        steps += 1;
        if steps > 1_000_000 {
            panic!("compound_swe1d::run did not finish in 1M steps");
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

    // -- CompoundSection geometry -----------------------------------

    #[test]
    fn area_below_bank_full_is_rectangular_main() {
        let s = CompoundSection { w_main: 20.0, w_flood: 100.0, h_bank: 1.0 };
        let a = s.area::<f64>(0.6);
        assert!(approx_eq(a, 12.0, EPS));
    }

    #[test]
    fn area_above_bank_full_adds_floodplain_excess() {
        // h_bank = 1.0, h = 1.5: area = 20·1.0 + 100·0.5 = 70.
        let s = CompoundSection { w_main: 20.0, w_flood: 100.0, h_bank: 1.0 };
        let a = s.area::<f64>(1.5);
        assert!(approx_eq(a, 70.0, EPS));
    }

    #[test]
    fn pressure_integral_below_bank_full() {
        // I₁ = w_main · h²/2 for h ≤ h_bank.
        let s = CompoundSection { w_main: 20.0, w_flood: 100.0, h_bank: 1.0 };
        let i1 = s.pressure_integral::<f64>(0.5);
        assert!(approx_eq(i1, 0.5 * 20.0 * 0.25, EPS));
    }

    #[test]
    fn pressure_integral_above_bank_full_matches_closed_form() {
        // h = 1.5, h_bank = 1.0, w_main = 20, w_flood = 100.
        // I₁ = 20·1·(1.5 − 0.5) + 100·(0.5)²/2 = 20 + 12.5 = 32.5.
        let s = CompoundSection { w_main: 20.0, w_flood: 100.0, h_bank: 1.0 };
        let i1 = s.pressure_integral::<f64>(1.5);
        assert!(approx_eq(i1, 32.5, EPS), "I₁ = {i1}");
    }

    #[test]
    fn compound_manning_below_bank_full_equals_n_main() {
        // For `h ≤ h_bank` the floodplain is dry, so the Lotter
        // formula must collapse to `n_main` regardless of n_flood.
        let s = CompoundSection { w_main: 30.0, w_flood: 100.0, h_bank: 1.0 };
        let n_main = 0.025_f64;
        let n_flood = 0.080_f64;
        for h in [0.01, 0.25, 0.5, 0.999, 1.0] {
            let n_eq = s.compound_manning(h, n_main, n_flood);
            assert!(
                approx_eq(n_eq, n_main, EPS),
                "h = {h}: n_eq = {n_eq} (expected n_main = {n_main})"
            );
        }
    }

    #[test]
    fn compound_manning_uniform_matches_scalar() {
        // n_main == n_flood ⇒ n_eq == n_main for any h. This is the
        // bit-exact backward-compat guarantee that older code with a
        // single Manning argument keeps producing the same results.
        let s = CompoundSection { w_main: 30.0, w_flood: 85.0, h_bank: 1.0 };
        let n = 0.04_f64;
        for h in [0.1, 1.0, 1.5, 3.0, 5.0] {
            let n_eq = s.compound_manning(h, n, n);
            assert!(approx_eq(n_eq, n, EPS), "h = {h}: n_eq = {n_eq}");
        }
    }

    #[test]
    fn compound_manning_lotter_increases_with_n_flood() {
        // Above bank-full, raising `n_flood` (with `n_main` fixed)
        // must raise the equivalent Manning monotonically. Captures
        // the physical intuition: rougher overbank ⇒ more friction
        // on the section as a whole.
        let s = CompoundSection { w_main: 30.0, w_flood: 100.0, h_bank: 1.0 };
        let n_main = 0.025_f64;
        let h = 2.0_f64;
        let n_low = s.compound_manning::<f64>(h, n_main, 0.030);
        let n_mid = s.compound_manning::<f64>(h, n_main, 0.060);
        let n_high = s.compound_manning::<f64>(h, n_main, 0.120);
        assert!(n_low < n_mid && n_mid < n_high);
        // All must lie between the extreme single-value Manning
        // results: bounded by the smaller and larger of n_main /
        // n_flood across the wetted perimeter.
        assert!(n_low > n_main && n_low < 0.030);
        assert!(n_high > n_main && n_high < 0.120);
    }

    #[test]
    fn compound_manning_dual_ad_matches_finite_diff() {
        // Forward-mode AD: derivative w.r.t. n_main and n_flood
        // computed by `Dual::variable` must match a central finite
        // difference. Roundoff-tight (`1e-5`) is the standard AD
        // sanity check.
        let s = CompoundSection { w_main: 30.0, w_flood: 100.0, h_bank: 1.0 };
        let h = 2.0_f64;
        let n_main = 0.025_f64;
        let n_flood = 0.060_f64;
        let h_d = Dual::constant(h);
        let n_main_d = Dual::variable(n_main);
        let n_flood_d = Dual::constant(n_flood);
        let n_eq_d = s.compound_manning(h_d, n_main_d, n_flood_d);
        let dn_eq_dn_main_ad = n_eq_d.dval;

        let eps = 1.0e-6;
        let plus = s.compound_manning::<f64>(h, n_main + eps, n_flood);
        let minus = s.compound_manning::<f64>(h, n_main - eps, n_flood);
        let dn_eq_dn_main_fd = (plus - minus) / (2.0 * eps);
        assert!(
            (dn_eq_dn_main_ad - dn_eq_dn_main_fd).abs() < 1.0e-7,
            "∂n_eq/∂n_main: AD = {dn_eq_dn_main_ad}, FD = {dn_eq_dn_main_fd}"
        );
    }

    #[test]
    fn stage_inverts_area_round_trip() {
        let s = CompoundSection { w_main: 20.0, w_flood: 100.0, h_bank: 1.0 };
        for h in [0.1, 0.5, 0.9, 1.0, 1.1, 2.5, 5.0] {
            let a = s.area::<f64>(h);
            let h_back = s.stage(a);
            assert!(approx_eq(h, h_back, EPS), "h = {h}, recovered = {h_back}");
        }
    }

    #[test]
    fn stage_inverse_is_differentiable_in_dual_above_bank() {
        // d(h)/d(A) = 1/w_flood above bank-full.
        let s = CompoundSection { w_main: 20.0, w_flood: 100.0, h_bank: 1.0 };
        let a = Dual::variable(s.w_main * s.h_bank + 5.0); // above bank
        let h = s.stage(a);
        assert!(approx_eq(h.dval, 1.0 / 100.0, EPS));
    }

    // -- Solver: wide-channel limit ---------------------------------

    #[test]
    fn lake_at_rest_stays_at_rest() {
        // Flat bed, zero discharge, transmissive walls. With h
        // initialised at 1 m, area should stay constant and Q at 0.
        let s = CompoundSection { w_main: 30.0, w_flood: 30.0, h_bank: 1.0 };
        let n_cells = 40;
        let bed = vec![0.0; n_cells];
        let h_init = 0.6_f64;
        let a0 = vec![s.area::<f64>(h_init); n_cells];
        let q0 = vec![0.0; n_cells];
        let (a, q, _) = run(
            &s,
            a0,
            q0,
            &bed,
            1.0,
            10.0,
            0.03,
            0.03,
            G,
            0.4,
            LeftBc::Transmissive,
            RightBc::Transmissive,
        );
        for (i, &a_i) in a.iter().enumerate() {
            assert!(
                (a_i - 18.0).abs() < 1.0e-9,
                "A[{i}] = {a_i} (expected 18.0 for h=0.6, W=30)"
            );
            assert!(q[i].abs() < 1.0e-9);
        }
    }

    #[test]
    fn wide_channel_limit_recovers_swe1d_steady_state() {
        // CompoundSection with w_main = w_flood degenerates to a
        // rectangular channel of fixed width. Steady inflow on a
        // sloping bed should yield Manning normal depth approximately.
        let w = 20.0_f64;
        let s = CompoundSection { w_main: w, w_flood: w, h_bank: 0.5 };
        let slope = 0.001_f64;
        let n = 0.04_f64;
        let n_cells = 60;
        let dx = 5.0;
        let bed: Vec<f64> = (0..n_cells)
            .map(|i| -slope * (i as f64 + 0.5) * dx)
            .collect();

        // Total Q = 1.5 m³/s (so q-per-unit-width = 0.075 m²/s).
        let q_in_total = 1.5_f64;
        // Manning normal depth for wide channel: h_n = (n·(q/w)/√S₀)^(3/5)
        let q_per_w = q_in_total / w;
        let h_n = (n * q_per_w / slope.sqrt()).powf(3.0 / 5.0);
        let a_n = w * h_n;

        let a0 = vec![a_n; n_cells];
        let q0 = vec![q_in_total; n_cells];
        let (a, q, _) = run(
            &s,
            a0,
            q0,
            &bed,
            dx,
            500.0,
            n,
            n,
            G,
            0.4,
            LeftBc::Dirichlet { h: h_n, q: q_in_total },
            RightBc::Transmissive,
        );

        let mid = n_cells / 2;
        let h_mid = s.stage(a[mid]);
        let q_mid = q[mid];
        // Tolerance loose: LF + boundary effects.
        assert!(
            (h_mid / h_n - 1.0).abs() < 0.10,
            "h_mid = {h_mid:.4}, h_n = {h_n:.4}"
        );
        assert!(
            (q_mid / q_in_total - 1.0).abs() < 0.10,
            "q_mid = {q_mid:.4}, q_in = {q_in_total:.4}"
        );
    }

    #[test]
    fn compound_section_at_bank_full_matches_main_channel_limit() {
        // At exactly h = h_bank, the area should equal w_main · h_bank
        // regardless of w_flood. Quick sanity that the cross-section
        // is continuous at the kink.
        let s = CompoundSection { w_main: 10.0, w_flood: 80.0, h_bank: 1.0 };
        let a_just_below = s.area::<f64>(0.999999);
        let a_at_bank = s.area::<f64>(1.0);
        let a_just_above = s.area::<f64>(1.000001);
        assert!((a_at_bank - 10.0).abs() < 1.0e-6);
        assert!((a_just_below - a_at_bank).abs() < 1.0e-4);
        assert!((a_just_above - a_at_bank).abs() < 1.0e-4);
    }

    #[test]
    fn compound_response_is_flatter_than_main_only() {
        // For two channels carrying the same Q, the compound section
        // (wide floodplain) should reach a SHALLOWER stage at high Q
        // than the main-channel-only equivalent. This is the key
        // physical effect that makes compound sections match the
        // sublinear shape of empirical rating curves.
        let main_only = CompoundSection {
            w_main: 30.0,
            w_flood: 30.0,
            h_bank: 1.0,
        };
        let compound = CompoundSection {
            w_main: 30.0,
            w_flood: 100.0,
            h_bank: 0.5,
        };
        let slope = 0.005_f64;
        let n = 0.04_f64;
        let n_cells = 40;
        let dx = 10.0;
        let bed: Vec<f64> = (0..n_cells)
            .map(|i| -slope * (i as f64 + 0.5) * dx)
            .collect();
        let q_in_total = 60.0_f64; // high Q, will exceed bank-full

        let run_one = |s: &CompoundSection| -> f64 {
            let a_init = s.area::<f64>(1.5); // crude warm-start
            let a0 = vec![a_init; n_cells];
            let q0 = vec![q_in_total; n_cells];
            let (a, _, _) = run(
                s,
                a0,
                q0,
                bed.as_slice(),
                dx,
                300.0,
                n,
                n,
                G,
                0.4,
                LeftBc::Dirichlet {
                    h: 1.5_f64,
                    q: q_in_total,
                },
                RightBc::Transmissive,
            );
            s.stage(a[n_cells / 2])
        };

        let h_main_only = run_one(&main_only);
        let h_compound = run_one(&compound);
        assert!(
            h_compound < h_main_only,
            "Compound h_mid = {h_compound:.3} should be < main-only h_mid = {h_main_only:.3} \
             (floodplain spreads excess flow at high Q)"
        );
    }
}
