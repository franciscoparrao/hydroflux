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
//! Wet–wet interfaces use the Davis (1988) bound
//! `s_L = min(u_n_L − c_L, u_n_R − c_R)`,
//! `s_R = max(u_n_L + c_L, u_n_R + c_R)`. This is well known to
//! under-estimate the dry-bed wave speed by up to a factor of two: at
//! a wet/dry front the leading edge of the rarefaction into the dry
//! region propagates at `u_n_W ± 2·c_W` (Toro 2009 §10.5.4,
//! "Two-Rarefaction Riemann Solver"), not at `u_n_W ± c_W`.
//!
//! We therefore branch on the wet/dry pattern at the face:
//!
//! - Both sides dry → zero flux (trivial).
//! - Left dry, right wet → `s_L = u_n_R − 2·c_R`, `s_R = u_n_R + c_R`.
//!   The rarefaction propagates leftward into the dry region.
//! - Left wet, right dry → `s_L = u_n_L − c_L`, `s_R = u_n_L + 2·c_L`.
//!   Symmetric: rarefaction propagates rightward into the dry region.
//! - Both wet → Davis bound (the original wet–wet branch).
//!
//! The wet/dry threshold is fixed to [`DRY_TOL`] inside this module
//! (independent of the higher-level `H_DRY` constant used by the FV
//! update), because the wave-speed branch only needs to detect cells
//! with depth at or below the numerical noise floor — at this
//! threshold the celerity `√(g·h)` falls below 10⁻⁴ m/s and the
//! distinction between "dry" and "wet" wave speeds is meaningful.

use crate::GRAVITY;
use crate::flux::{FluxX, FluxY};
use crate::state::Conserved2D;

/// Wet/dry threshold used internally by the Riemann solver to detect
/// dry cells when picking wave-speed estimates. Tighter than the
/// user-visible `H_DRY` constant: this is the numerical-noise floor
/// for the wave-speed branch, not a physical wet/dry definition.
const DRY_TOL: f64 = 1.0e-12;

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
    let h_l_wet = h_l > DRY_TOL;
    let h_r_wet = h_r > DRY_TOL;

    // Both sides dry: no flux, no Riemann problem to solve.
    if !h_l_wet && !h_r_wet {
        return (0.0, 0.0, 0.0);
    }

    let cl = if h_l_wet { (GRAVITY * h_l).sqrt() } else { 0.0 };
    let cr = if h_r_wet { (GRAVITY * h_r).sqrt() } else { 0.0 };
    let un_l = if h_l_wet { qn_l / h_l } else { 0.0 };
    let un_r = if h_r_wet { qn_r / h_r } else { 0.0 };
    let ut_l = if h_l_wet { qt_l / h_l } else { 0.0 };
    let ut_r = if h_r_wet { qt_r / h_r } else { 0.0 };

    // Wave-speed estimate. The Davis (1988) bound is correct for two
    // wet states; at a wet/dry interface the leading edge of the
    // rarefaction into the dry region propagates at `u_n_W ± 2·c_W`
    // (Toro §10.5.4), not at `u_n_W ± c_W`. We branch on the wet/dry
    // pattern to use the right estimate.
    let (sl, sr) = match (h_l_wet, h_r_wet) {
        (true, true) => ((un_l - cl).min(un_r - cr), (un_l + cl).max(un_r + cr)),
        (false, true) => (un_r - 2.0 * cr, un_r + cr),
        (true, false) => (un_l - cl, un_l + 2.0 * cl),
        (false, false) => unreachable!("handled by the both-dry early return above"),
    };

    // Physical fluxes in the normal direction. Tangential flux is
    // qn · u_t (advected by the normal velocity).
    let (fm_l, fn_l, ft_l) = if h_l_wet {
        (qn_l, qn_l * un_l + 0.5 * GRAVITY * h_l * h_l, qn_l * ut_l)
    } else {
        (0.0, 0.0, 0.0)
    };
    let (fm_r, fn_r, ft_r) = if h_r_wet {
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
    fn right_dry_face_produces_outflow_into_dry_region() {
        // Water at rest on the left, fully dry on the right. The
        // analytical solution is a left-going rarefaction whose
        // leading edge (the dry front) propagates at +2·c_L. HLLC with
        // the correct wave speed must produce strictly positive mass
        // flux through the face — water moves from L into the dry R.
        let h_l = 1.0;
        let ul = Conserved2D::new(h_l, 0.0, 0.0);
        let ur = Conserved2D::DRY;
        let f = hllc_flux_x(ul, ur);
        let c_l = (GRAVITY * h_l).sqrt();
        // Mass flux on the dry-bed Stoker problem must be positive and
        // bounded above by h_L · 2·c_L (the dry-front speed times the
        // depth on the wet side, which is the upper bound from
        // conservation through the leading edge).
        assert!(
            f.mass > 0.0,
            "mass flux into dry region must be positive, got {}",
            f.mass
        );
        assert!(
            f.mass < h_l * 2.0 * c_l,
            "mass flux must be bounded by h_L · 2·c_L = {}, got {}",
            h_l * 2.0 * c_l,
            f.mass
        );
    }

    #[test]
    fn left_dry_face_produces_inflow_from_right_region() {
        // Symmetric mirror image: dry on the left, water at rest on the
        // right. The dry front propagates leftward; mass flux at the
        // face must be strictly negative (water moves from R into the
        // dry L).
        let h_r = 1.0;
        let ul = Conserved2D::DRY;
        let ur = Conserved2D::new(h_r, 0.0, 0.0);
        let f = hllc_flux_x(ul, ur);
        let c_r = (GRAVITY * h_r).sqrt();
        assert!(
            f.mass < 0.0,
            "mass flux from wet R into dry L must be negative, got {}",
            f.mass
        );
        assert!(
            f.mass.abs() < h_r * 2.0 * c_r,
            "|mass flux| bounded by h_R · 2·c_R = {}, got {}",
            h_r * 2.0 * c_r,
            f.mass.abs()
        );
    }

    #[test]
    fn dry_bed_mass_flux_matches_two_rarefaction_closed_form() {
        // Closed-form check. For Stoker dry-bed (u_L = 0, dry on the
        // right, wave speeds s_L = −c_L, s_R = +2·c_L), the HLLC
        // sampling of F_L* at ξ = 0 gives:
        //   h*_L  = h_L · (s_L − u_L) / (s_L − s*) = h_L / 3
        //   F.mass = s_L · (h*_L − h_L) = (−c_L) · (−2·h_L/3) = (2/3)·c_L·h_L
        // This must hold to machine precision because everything in
        // the path is exact algebra once the wave speeds are picked.
        let h_l = 1.0;
        let ul = Conserved2D::new(h_l, 0.0, 0.0);
        let f = hllc_flux_x(ul, Conserved2D::DRY);
        let c_l = (GRAVITY * h_l).sqrt();
        let expected = (2.0 / 3.0) * c_l * h_l;
        assert_relative_eq!(f.mass, expected, epsilon = 1e-12);

        // Cross-check: this exceeds what the Davis bound would yield
        // for the same problem (Davis: s_R = +c_L, giving F.mass =
        // c_L · h_L / 2). The ratio (2/3) / (1/2) = 4/3 quantifies the
        // wave-speed correction at a wet/dry interface.
        let davis_flux = 0.5 * c_l * h_l;
        assert!(
            f.mass > davis_flux,
            "two-rarefaction flux {} must exceed Davis-bound flux {}",
            f.mass,
            davis_flux
        );
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
