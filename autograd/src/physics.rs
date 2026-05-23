//! Generic shallow-water primitives written over [`Real`].
//!
//! These are the SWE building blocks (celerity, friction slope, flux
//! function) re-expressed without committing to `f64`. They evaluate
//! identically to their concrete `f64` counterparts in
//! `hydroflux-solver-2d`, but accept any `T: Real` — including
//! [`crate::Dual`] — so the gradient of any quantity derived from them
//! is available by simply seeding the input.
//!
//! Physical constants (`gravity`, `manning_n`) are passed as `f64`
//! parameters rather than baked in. That keeps the call site explicit
//! about which inputs are being differentiated: only the `T`-typed
//! arguments carry derivatives.

use crate::Real;

/// Gravity-wave celerity `c = √(g · max(h, 0))`.
///
/// The `max` with zero handles cells that round below zero in a
/// finite-precision step; for `Dual` inputs the derivative at the
/// kink is zero (see [`crate::Dual::max`]).
pub fn celerity<T: Real>(h: T, gravity: f64) -> T {
    (h.max(T::zero()) * gravity).sqrt()
}

/// Manning friction slope in 1D:
/// `S_f = n² · q · |q| / h^(10/3)`.
///
/// The `q · |q|` form preserves sign so the slope opposes flow
/// direction when `q < 0`. Caller must guarantee `h > 0`; the
/// function does not guard the division.
pub fn manning_friction_slope_1d<T: Real>(h: T, q: T, manning_n: f64) -> T {
    let abs_q = q.abs();
    (q * abs_q) * (manning_n * manning_n) / h.powf(10.0 / 3.0)
}

/// Manning friction slope magnitude in 2D:
/// `S_f = n² · |q| · ‖q‖ / h^(10/3)` where `q = (q_x, q_y)`.
///
/// Returns the scalar slope; direction is recovered by the caller as
/// `-S_f · q / ‖q‖` (i.e., friction opposes velocity). Caller must
/// guarantee `h > 0`.
pub fn manning_friction_slope_2d<T: Real>(h: T, qx: T, qy: T, manning_n: f64) -> T {
    let speed_q = (qx * qx + qy * qy).sqrt();
    speed_q * speed_q * (manning_n * manning_n) / h.powf(10.0 / 3.0)
}

/// 1D Saint-Venant flux `F(U) = (q, q²/h + ½ g h²)`.
///
/// Caller must guarantee `h > 0`.
pub fn flux_swe_1d<T: Real>(h: T, q: T, gravity: f64) -> (T, T) {
    let mass = q;
    let momentum = (q * q) / h + (h * h) * (0.5 * gravity);
    (mass, momentum)
}

/// 2D Saint-Venant `x`-direction flux:
/// `F_x = (h u, h u² + ½ g h², h u v)`.
///
/// Inputs are conserved `(h, hu, hv)`. Caller must guarantee `h > 0`.
pub fn flux_swe_2d_x<T: Real>(h: T, hu: T, hv: T, gravity: f64) -> (T, T, T) {
    let u = hu / h;
    let mass = hu;
    let x_momentum = hu * u + (h * h) * (0.5 * gravity);
    let y_momentum = hu * (hv / h);
    (mass, x_momentum, y_momentum)
}

/// 2D Saint-Venant `y`-direction flux:
/// `G_y = (h v, h u v, h v² + ½ g h²)`.
pub fn flux_swe_2d_y<T: Real>(h: T, hu: T, hv: T, gravity: f64) -> (T, T, T) {
    let v = hv / h;
    let mass = hv;
    let x_momentum = hv * (hu / h);
    let y_momentum = hv * v + (h * h) * (0.5 * gravity);
    (mass, x_momentum, y_momentum)
}

/// Manning normal depth `h_n = (n · |q| / √S₀)^(3/5)` for a wide
/// rectangular channel under steady uniform flow.
///
/// Differentiable with respect to `h`-typed quantities only via the
/// caller (everything inside is f64-parameterised). To differentiate
/// w.r.t. `q`, pass `q: T` and `slope: f64` and use
/// [`manning_normal_depth_t`].
pub fn manning_normal_depth(q: f64, manning_n: f64, slope: f64) -> f64 {
    (manning_n * q.abs() / slope.sqrt()).powf(3.0 / 5.0)
}

/// Manning normal depth, differentiable in `q`.
///
/// `h_n(q) = (n · |q| / √S₀)^(3/5)`. Useful as a ghost-state value
/// for Discharge boundary conditions over dry beds, where we want
/// the gradient of the inflow profile w.r.t. the imposed discharge.
pub fn manning_normal_depth_t<T: Real>(q: T, manning_n: f64, slope: f64) -> T {
    (q.abs() * (manning_n / slope.sqrt())).powf(3.0 / 5.0)
}

/// Critical depth `h_c = (q² / g)^(1/3)`.
pub fn critical_depth<T: Real>(q: T, gravity: f64) -> T {
    (q * q / gravity).powf(1.0 / 3.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Dual;
    use approx::assert_relative_eq;

    const G: f64 = 9.81;
    const EPS: f64 = 1.0e-10;

    // -- celerity ----------------------------------------------------

    #[test]
    fn celerity_matches_concrete_solver() {
        // h = 4 m → c = √(9.81 · 4) = √39.24.
        let c_f = celerity::<f64>(4.0, G);
        assert_relative_eq!(c_f, (G * 4.0).sqrt(), epsilon = EPS);
    }

    #[test]
    fn celerity_dual_gradient_d_c_d_h_at_unit_depth() {
        // d/dh √(g h) = √g / (2 √h). At h = 1: √g / 2.
        let h = Dual::variable(1.0);
        let c = celerity(h, G);
        assert_relative_eq!(c.val, G.sqrt(), epsilon = EPS);
        assert_relative_eq!(c.dval, G.sqrt() / 2.0, epsilon = EPS);
    }

    #[test]
    fn celerity_clamps_dry_to_zero() {
        // Negative h must clamp before sqrt.
        let c = celerity::<f64>(-0.1, G);
        assert_eq!(c, 0.0);
        // Dual gradient at the kink is taken as 0.
        let cd = celerity(Dual::variable(-0.1), G);
        assert_eq!(cd.val, 0.0);
        assert_eq!(cd.dval, 0.0);
    }

    // -- Manning friction --------------------------------------------

    #[test]
    fn manning_1d_value_and_gradient_w_r_t_h() {
        // S_f = n²·q·|q| / h^(10/3); d/dh S_f = -(10/3) · n²·q·|q| / h^(13/3).
        let q = 1.5_f64;
        let n = 0.04_f64;
        let h_val = 0.8_f64;

        let s_concrete = manning_friction_slope_1d::<f64>(h_val, q, n);
        let analytic_val = n * n * q * q.abs() / h_val.powf(10.0 / 3.0);
        assert_relative_eq!(s_concrete, analytic_val, epsilon = EPS);

        let h = Dual::variable(h_val);
        let s = manning_friction_slope_1d(h, Dual::constant(q), n);
        let analytic_grad = -(10.0 / 3.0) * n * n * q * q.abs() / h_val.powf(13.0 / 3.0);
        assert_relative_eq!(s.val, analytic_val, epsilon = EPS);
        assert_relative_eq!(s.dval, analytic_grad, epsilon = EPS);
    }

    #[test]
    fn manning_2d_reduces_to_1d_when_qy_is_zero() {
        // With q_y = 0, the 2D formula must give the same magnitude as
        // the 1D formula with q = q_x.
        let h = 0.5_f64;
        let qx = 1.2_f64;
        let n = 0.035_f64;
        let s_1d = manning_friction_slope_1d::<f64>(h, qx, n);
        let s_2d = manning_friction_slope_2d::<f64>(h, qx, 0.0, n);
        assert_relative_eq!(s_1d, s_2d, epsilon = EPS);
    }

    // -- Flux --------------------------------------------------------

    #[test]
    fn flux_1d_matches_concrete_value_with_static_water() {
        // q = 0 → F = (0, ½·g·h²).
        let (mass, momentum) = flux_swe_1d::<f64>(2.0, 0.0, G);
        assert_eq!(mass, 0.0);
        assert_relative_eq!(momentum, 0.5 * G * 4.0, epsilon = EPS);
    }

    #[test]
    fn flux_1d_gradient_d_momentum_d_h_under_static_water() {
        // d/dh (q²/h + ½·g·h²) at q = 0, h = h₀ → g·h₀.
        let h = Dual::variable(2.0);
        let q = Dual::constant(0.0);
        let (_, momentum) = flux_swe_1d(h, q, G);
        assert_relative_eq!(momentum.val, 0.5 * G * 4.0, epsilon = EPS);
        assert_relative_eq!(momentum.dval, G * 2.0, epsilon = EPS);
    }

    #[test]
    fn flux_2d_x_consistency_with_static_water() {
        // hu = hv = 0; F_x = (0, ½·g·h², 0).
        let h = 1.5;
        let (m, x, y) = flux_swe_2d_x::<f64>(h, 0.0, 0.0, G);
        assert_eq!(m, 0.0);
        assert_relative_eq!(x, 0.5 * G * h * h, epsilon = EPS);
        assert_eq!(y, 0.0);
    }

    #[test]
    fn flux_2d_y_consistency_with_static_water() {
        let h = 1.5;
        let (m, x, y) = flux_swe_2d_y::<f64>(h, 0.0, 0.0, G);
        assert_eq!(m, 0.0);
        assert_eq!(x, 0.0);
        assert_relative_eq!(y, 0.5 * G * h * h, epsilon = EPS);
    }

    // -- Manning normal depth ----------------------------------------

    #[test]
    fn manning_normal_depth_matches_closed_form() {
        // q = 2 m²/s, n = 0.04, S₀ = 0.001 → h_n = (0.04·2/√0.001)^(3/5).
        let h_n = manning_normal_depth(2.0, 0.04, 0.001);
        let analytic = (0.04 * 2.0 / 0.001_f64.sqrt()).powf(3.0 / 5.0);
        assert_relative_eq!(h_n, analytic, epsilon = EPS);
    }

    #[test]
    fn manning_normal_depth_gradient_in_q_matches_analytic() {
        // h_n(q) = (n · q / √S₀)^(3/5) for q > 0.
        // dh_n/dq = (3/5) · (n/√S₀)^(3/5) · q^(-2/5).
        let n = 0.04_f64;
        let s = 0.001_f64;
        let q_val = 2.0_f64;
        let q = Dual::variable(q_val);
        let h_n = manning_normal_depth_t(q, n, s);
        let analytic_grad = (3.0 / 5.0) * (n / s.sqrt()).powf(3.0 / 5.0) * q_val.powf(-2.0 / 5.0);
        assert_relative_eq!(h_n.dval, analytic_grad, epsilon = EPS);
    }

    // -- Critical depth ----------------------------------------------

    #[test]
    fn critical_depth_value_and_gradient() {
        // h_c(q) = (q²/g)^(1/3). At q = 2: (4/9.81)^(1/3) ≈ 0.7416.
        // dh_c/dq = (2/3) · q^(-1/3) · g^(-1/3).
        let q = Dual::variable(2.0);
        let h_c = critical_depth(q, G);
        let val_analytic = (4.0_f64 / G).powf(1.0 / 3.0);
        let grad_analytic = (2.0 / 3.0) * 2.0_f64.powf(-1.0 / 3.0) * G.powf(-1.0 / 3.0);
        assert_relative_eq!(h_c.val, val_analytic, epsilon = EPS);
        assert_relative_eq!(h_c.dval, grad_analytic, epsilon = EPS);
    }

    // -- Cross-backend round-trip ------------------------------------

    #[test]
    fn concrete_and_dual_value_agree_for_full_2d_x_flux() {
        // Same physical input through f64 and Dual must give the same value.
        let (h_v, hu_v, hv_v) = (1.2, 2.5, -0.6);
        let (mf, xf, yf) = flux_swe_2d_x::<f64>(h_v, hu_v, hv_v, G);
        let (md, xd, yd) = flux_swe_2d_x::<Dual>(
            Dual::variable(h_v),
            Dual::constant(hu_v),
            Dual::constant(hv_v),
            G,
        );
        assert_relative_eq!(mf, md.val, epsilon = EPS);
        assert_relative_eq!(xf, xd.val, epsilon = EPS);
        assert_relative_eq!(yf, yd.val, epsilon = EPS);
    }
}
