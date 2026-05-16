//! HLL Riemann solver for the 1D Saint-Venant system.
//!
//! Wave-speed estimate follows Toro (2009, §10.5.1), Davis (1988) bound:
//!
//! ```text
//!   S_L = min(u_L - c_L, u_R - c_R)
//!   S_R = max(u_L + c_L, u_R + c_R)
//! ```
//!
//! Dry-bed handling: cells with `h = 0` carry `c = 0` here, giving an
//! approximate but consistent flux. The exact two-rarefaction estimate
//! (Toro 2009 §10.5.4) is deferred to the iteration that introduces a
//! dam-break-on-dry-bed benchmark.

use crate::GRAVITY;
use crate::flux::Flux;
use crate::state::Conserved;

/// HLL Riemann flux at the interface between left state `ul` and right
/// state `ur`. Returns the numerical flux `F*` used by the FV update.
pub fn hll_flux(ul: Conserved, ur: Conserved) -> Flux {
    let cl = (GRAVITY * ul.h.max(0.0)).sqrt();
    let cr = (GRAVITY * ur.h.max(0.0)).sqrt();
    let u_left = if ul.h > 0.0 { ul.hu / ul.h } else { 0.0 };
    let u_right = if ur.h > 0.0 { ur.hu / ur.h } else { 0.0 };

    let sl = (u_left - cl).min(u_right - cr);
    let sr = (u_left + cl).max(u_right + cr);

    let fl = Flux::from_state(ul);
    let fr = Flux::from_state(ur);

    if sl >= 0.0 {
        fl
    } else if sr <= 0.0 {
        fr
    } else {
        let denom = sr - sl;
        Flux {
            mass: (sr * fl.mass - sl * fr.mass + sl * sr * (ur.h - ul.h)) / denom,
            momentum: (sr * fl.momentum - sl * fr.momentum + sl * sr * (ur.hu - ul.hu)) / denom,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn consistency_identical_states() {
        // F*(U, U) must equal F(U).
        let u = Conserved::new(1.5, 0.6);
        let star = hll_flux(u, u);
        let physical = Flux::from_state(u);
        assert_relative_eq!(star.mass, physical.mass, epsilon = 1e-12);
        assert_relative_eq!(star.momentum, physical.momentum, epsilon = 1e-12);
    }

    #[test]
    fn dry_dry_interface_is_zero_flux() {
        let f = hll_flux(Conserved::DRY, Conserved::DRY);
        // Both states dry → c_L = c_R = 0, u_L = u_R = 0, S_L = S_R = 0,
        // takes the `sl >= 0` branch and returns F(U_L) = ZERO.
        assert_eq!(f, Flux::ZERO);
    }

    #[test]
    fn symmetric_dam_break_has_zero_mass_flux_at_center() {
        // Symmetric Riemann problem about u = 0: left and right depths
        // mirrored, both at rest. The HLL star state must have zero mass
        // flux at the interface (lake-at-rest symmetry around the diaphragm
        // is broken only by the depth jump, but mass flux must integrate to
        // zero by symmetry of the wave fan).
        //
        // For an *asymmetric* dam break (h_L > h_R, both at rest), HLL gives
        // a non-zero star flux. Here both sides equal → trivial case.
        let h = 1.0;
        let state = Conserved::new(h, 0.0);
        let f = hll_flux(state, state);
        assert_relative_eq!(f.mass, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn dam_break_wet_bed_has_positive_mass_flux() {
        // Classic dry-vs-wet at rest: left deeper, right shallower, both
        // initially at rest. HLL must produce a positive mass flux from
        // left to right (water moves downhill).
        let ul = Conserved::new(2.0, 0.0);
        let ur = Conserved::new(0.5, 0.0);
        let f = hll_flux(ul, ur);
        assert!(
            f.mass > 0.0,
            "expected positive mass flux from deep to shallow, got {}",
            f.mass
        );
    }

    #[test]
    fn supercritical_flow_returns_upwind_flux() {
        // Supercritical: |u| > c. Both wave speeds have the same sign,
        // HLL collapses to upwind. Here u_L = u_R = 10 m/s, c ≈ 3.13 m/s.
        let h = 1.0;
        let u = 10.0;
        let ul = Conserved::new(h, h * u);
        let ur = Conserved::new(h, h * u);
        let f = hll_flux(ul, ur);
        let physical = Flux::from_state(ul);
        assert_relative_eq!(f.mass, physical.mass, epsilon = 1e-12);
        assert_relative_eq!(f.momentum, physical.momentum, epsilon = 1e-12);
    }
}
