//! HLLC Riemann solver for the 2D Saint-Venant system.
//!
//! Wave-speed estimate follows Toro (2009, §10.5.1) using the
//! Davis (1988) bound on the outer waves, with the HLLC contact-wave
//! velocity computed from the algebraic identity in Toro §10.3.
//!
//! # Rotational invariance
//!
//! On an axis-aligned face with outward normal `n`, the Riemann
//! problem reduces to a 1D problem in the face-normal direction
//! (Toro §16.4). We decompose each cell state into:
//!
//! - `q_n_K = h_K · u_n_K`: momentum per unit area in the **normal**
//!   direction;
//! - `q_t_K = h_K · u_t_K`: momentum per unit area in the **tangential**
//!   direction.
//!
//! The implementation lives in a single private function
//! [`hllc_normal_flux`]; the public [`hllc_flux_x`] and [`hllc_flux_y`]
//! wrappers simply pass the right components.
//!
//! # Tangential momentum
//!
//! In the HLLC star state for SWE the *tangential* velocity is
//! constant across the contact wave (it advects with the contact
//! rather than jumping). This is the key property that distinguishes
//! HLLC from HLL in 2D: HLL would smear an artificial jump in the
//! tangential momentum.
//!
//! # Dry-bed handling
//!
//! The current implementation uses the same approximate dry treatment
//! as the 1D solver (cells with `h ≤ 0` carry `c = 0`). For
//! dam-break-on-dry benchmarks (Toro §10.5.4) a two-rarefaction
//! wave-speed estimate is needed; deferred to a follow-up iteration.

use crate::GRAVITY;
use crate::flux::{FluxX, FluxY};
use crate::state::Conserved2D;

/// HLLC Riemann flux through an `x`-face between left state `ul`
/// and right state `ur`. Mass and `x`-momentum (normal) flow are
/// computed by HLLC; the `y`-momentum (tangential) flux advects with
/// the contact wave.
pub fn hllc_flux_x(ul: Conserved2D, ur: Conserved2D) -> FluxX {
    // x-face: normal direction is x, so normal momentum = hu and
    // tangential momentum = hv.
    let (mass, normal_mom, tangential_mom) =
        hllc_normal_flux(ul.h, ul.hu, ul.hv, ur.h, ur.hu, ur.hv);
    FluxX {
        mass,
        x_momentum: normal_mom,
        y_momentum: tangential_mom,
    }
}

/// HLLC Riemann flux through a `y`-face between left state `ul`
/// (lower row index, `+y` side after orientation) and right state `ur`
/// (higher row index). Mass and `y`-momentum (normal) flow are
/// computed by HLLC; the `x`-momentum (tangential) flux advects with
/// the contact wave.
pub fn hllc_flux_y(ul: Conserved2D, ur: Conserved2D) -> FluxY {
    // y-face: normal direction is y, so normal momentum = hv and
    // tangential momentum = hu.
    let (mass, normal_mom, tangential_mom) =
        hllc_normal_flux(ul.h, ul.hv, ul.hu, ur.h, ur.hv, ur.hu);
    FluxY {
        mass,
        x_momentum: tangential_mom,
        y_momentum: normal_mom,
    }
}

/// Core HLLC solver in the face-normal direction.
///
/// Inputs are decomposed into normal / tangential momentum components.
/// Returns `(mass_flux, normal_momentum_flux, tangential_momentum_flux)`
/// at the face.
fn hllc_normal_flux(
    h_l: f64,
    qn_l: f64,
    qt_l: f64,
    h_r: f64,
    qn_r: f64,
    qt_r: f64,
) -> (f64, f64, f64) {
    // Wave speeds (Davis 1988 bound).
    let cl = (GRAVITY * h_l.max(0.0)).sqrt();
    let cr = (GRAVITY * h_r.max(0.0)).sqrt();
    let un_l = if h_l > 0.0 { qn_l / h_l } else { 0.0 };
    let un_r = if h_r > 0.0 { qn_r / h_r } else { 0.0 };
    let ut_l = if h_l > 0.0 { qt_l / h_l } else { 0.0 };
    let ut_r = if h_r > 0.0 { qt_r / h_r } else { 0.0 };

    let sl = (un_l - cl).min(un_r - cr);
    let sr = (un_l + cl).max(un_r + cr);

    // Physical fluxes in the normal direction. Tangential flux is
    // qn · u_t (advected by the normal velocity).
    let (fm_l, fn_l, ft_l) = if h_l > 0.0 {
        (qn_l, qn_l * un_l + 0.5 * GRAVITY * h_l * h_l, qn_l * ut_l)
    } else {
        (0.0, 0.0, 0.0)
    };
    let (fm_r, fn_r, ft_r) = if h_r > 0.0 {
        (qn_r, qn_r * un_r + 0.5 * GRAVITY * h_r * h_r, qn_r * ut_r)
    } else {
        (0.0, 0.0, 0.0)
    };

    if sl >= 0.0 {
        return (fm_l, fn_l, ft_l);
    }
    if sr <= 0.0 {
        return (fm_r, fn_r, ft_r);
    }

    // Contact wave speed s* (Toro §10.3). For wet-wet states with
    // sl < 0 < sr the denominator is bounded away from zero.
    let denom = h_r * (un_r - sr) - h_l * (un_l - sl);
    let s_star = if denom.abs() > 1.0e-12 {
        (sl * h_r * (un_r - sr) - sr * h_l * (un_l - sl)) / denom
    } else {
        0.5 * (sl + sr)
    };

    // Star-state depths.
    let h_star_l = h_l * (sl - un_l) / (sl - s_star);
    let h_star_r = h_r * (sr - un_r) / (sr - s_star);

    // HLLC sampling at the face (xi = 0 in the similarity variable).
    if s_star >= 0.0 {
        // In the F_L* region: F* = F_L + sL · (U_L* - U_L).
        let mass = fm_l + sl * (h_star_l - h_l);
        let normal = fn_l + sl * (h_star_l * s_star - qn_l);
        // Tangential velocity advects through the contact: u_t* = u_t_L.
        let tangential = ft_l + sl * (h_star_l * ut_l - qt_l);
        (mass, normal, tangential)
    } else {
        let mass = fm_r + sr * (h_star_r - h_r);
        let normal = fn_r + sr * (h_star_r * s_star - qn_r);
        let tangential = ft_r + sr * (h_star_r * ut_r - qt_r);
        (mass, normal, tangential)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn x_consistency_identical_states() {
        // F*(U, U) must equal F_x(U).
        let u = Conserved2D::new(1.5, 0.6, -0.3);
        let star = hllc_flux_x(u, u);
        let physical = FluxX::from_state(u);
        assert_relative_eq!(star.mass, physical.mass, epsilon = 1e-12);
        assert_relative_eq!(star.x_momentum, physical.x_momentum, epsilon = 1e-12);
        assert_relative_eq!(star.y_momentum, physical.y_momentum, epsilon = 1e-12);
    }

    #[test]
    fn y_consistency_identical_states() {
        let u = Conserved2D::new(1.5, 0.6, -0.3);
        let star = hllc_flux_y(u, u);
        let physical = FluxY::from_state(u);
        assert_relative_eq!(star.mass, physical.mass, epsilon = 1e-12);
        assert_relative_eq!(star.x_momentum, physical.x_momentum, epsilon = 1e-12);
        assert_relative_eq!(star.y_momentum, physical.y_momentum, epsilon = 1e-12);
    }

    #[test]
    fn dry_dry_interface_is_zero_flux() {
        let dry = Conserved2D::DRY;
        assert_eq!(hllc_flux_x(dry, dry), FluxX::ZERO);
        assert_eq!(hllc_flux_y(dry, dry), FluxY::ZERO);
    }

    #[test]
    fn x_dam_break_wet_bed_positive_mass_flux() {
        // Symmetric water depth jump, both at rest. Water moves from
        // deep (left) to shallow (right) — mass flux must be positive.
        let ul = Conserved2D::new(2.0, 0.0, 0.0);
        let ur = Conserved2D::new(0.5, 0.0, 0.0);
        let f = hllc_flux_x(ul, ur);
        assert!(
            f.mass > 0.0,
            "expected positive mass flux from deep to shallow, got {}",
            f.mass
        );
    }

    #[test]
    fn y_dam_break_wet_bed_positive_mass_flux() {
        // Analogous in y direction (jump along y, both at rest).
        let ul = Conserved2D::new(2.0, 0.0, 0.0);
        let ur = Conserved2D::new(0.5, 0.0, 0.0);
        let f = hllc_flux_y(ul, ur);
        assert!(
            f.mass > 0.0,
            "expected positive y mass flux from deep to shallow, got {}",
            f.mass
        );
    }

    #[test]
    fn x_supercritical_flow_returns_upwind_flux() {
        // |u| > c on both sides: both wave speeds positive, HLLC
        // collapses to upwind = F_L.
        let h = 1.0;
        let u = 10.0;
        let v = 1.0;
        let state = Conserved2D::new(h, h * u, h * v);
        let f = hllc_flux_x(state, state);
        let physical = FluxX::from_state(state);
        assert_relative_eq!(f.mass, physical.mass, epsilon = 1e-12);
        assert_relative_eq!(f.x_momentum, physical.x_momentum, epsilon = 1e-12);
        assert_relative_eq!(f.y_momentum, physical.y_momentum, epsilon = 1e-12);
    }

    #[test]
    fn x_face_preserves_tangential_velocity_through_contact() {
        // Asymmetric dam break in the x direction with a *constant*
        // tangential velocity v on both sides. HLLC should propagate
        // the same v through the contact and the y-momentum flux must
        // equal (mass flux) · v.
        let v = 2.5;
        let ul = Conserved2D::new(2.0, 0.0, 2.0 * v);
        let ur = Conserved2D::new(0.5, 0.0, 0.5 * v);
        let f = hllc_flux_x(ul, ur);
        assert!(f.mass > 0.0);
        assert_relative_eq!(f.y_momentum, f.mass * v, epsilon = 1e-10);
    }

    #[test]
    fn y_face_preserves_tangential_velocity_through_contact() {
        // Same as above but tangential is u along y-face.
        let u = -1.3;
        let ul = Conserved2D::new(2.0, 2.0 * u, 0.0);
        let ur = Conserved2D::new(0.5, 0.5 * u, 0.0);
        let f = hllc_flux_y(ul, ur);
        assert!(f.mass > 0.0);
        assert_relative_eq!(f.x_momentum, f.mass * u, epsilon = 1e-10);
    }

    #[test]
    fn lake_at_rest_x_face_no_mass_flux() {
        // Two cells at the same depth, no velocity, on an x-face.
        // Mass flux must be exactly zero (the Riemann problem is
        // trivial and HLLC reduces to the consistent flux).
        let h = 2.0;
        let state = Conserved2D::new(h, 0.0, 0.0);
        let f = hllc_flux_x(state, state);
        assert_relative_eq!(f.mass, 0.0, epsilon = 1e-12);
        assert_relative_eq!(f.x_momentum, 0.5 * GRAVITY * h * h, epsilon = 1e-12);
        assert_relative_eq!(f.y_momentum, 0.0, epsilon = 1e-12);
    }
}
