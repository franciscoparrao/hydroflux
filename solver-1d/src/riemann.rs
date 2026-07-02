//! HLL Riemann solver for the 1D Saint-Venant system.
//!
//! Wave-speed estimate follows Toro (2009, §10.5.1), Davis (1988) bound
//! for two wet states:
//!
//! ```text
//!   S_L = min(u_L - c_L, u_R - c_R)
//!   S_R = max(u_L + c_L, u_R + c_R)
//! ```
//!
//! At a wet/dry interface the Davis bound underestimates the front: the
//! leading edge of the rarefaction into the dry region propagates at
//! `u_W ± 2·c_W` (two-rarefaction estimate, Toro 2009 §10.5.4). The
//! solver branches on the wet/dry pattern to use the right estimate,
//! mirroring `hydroflux-solver-2d`.

use crate::GRAVITY;
use crate::flux::Flux;
use crate::state::Conserved;

/// Wet/dry threshold used internally by the Riemann solver to detect
/// dry cells when picking wave-speed estimates. Tighter than the
/// user-visible `H_DRY` constant: this is the numerical-noise floor
/// for the wave-speed branch, not a physical wet/dry definition.
const DRY_TOL: f64 = 1.0e-12;

/// HLL Riemann flux at the interface between left state `ul` and right
/// state `ur`. Returns the numerical flux `F*` used by the FV update.
pub fn hll_flux(ul: Conserved, ur: Conserved) -> Flux {
    let l_wet = ul.h > DRY_TOL;
    let r_wet = ur.h > DRY_TOL;

    // Both sides dry: no flux, no Riemann problem to solve.
    if !l_wet && !r_wet {
        return Flux::ZERO;
    }

    let cl = if l_wet { (GRAVITY * ul.h).sqrt() } else { 0.0 };
    let cr = if r_wet { (GRAVITY * ur.h).sqrt() } else { 0.0 };
    let u_left = if l_wet { ul.hu / ul.h } else { 0.0 };
    let u_right = if r_wet { ur.hu / ur.h } else { 0.0 };

    let (sl, sr) = match (l_wet, r_wet) {
        (true, true) => (
            (u_left - cl).min(u_right - cr),
            (u_left + cl).max(u_right + cr),
        ),
        (false, true) => (u_right - 2.0 * cr, u_right + cr),
        (true, false) => (u_left - cl, u_left + 2.0 * cl),
        (false, false) => unreachable!("handled by the both-dry early return above"),
    };

    // Physical fluxes, written from the reconstructed velocities so a
    // near-dry state (h ≤ DRY_TOL with residual hu) cannot blow up in
    // `hu²/h`; the dry side carries exactly zero flux.
    let fl = if l_wet {
        Flux {
            mass: ul.hu,
            momentum: ul.hu * u_left + 0.5 * GRAVITY * ul.h * ul.h,
        }
    } else {
        Flux::ZERO
    };
    let fr = if r_wet {
        Flux {
            mass: ur.hu,
            momentum: ur.hu * u_right + 0.5 * GRAVITY * ur.h * ur.h,
        }
    } else {
        Flux::ZERO
    };

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
    fn wet_dry_interface_pushes_mass_into_dry_side() {
        // Dam break onto dry bed at the initial instant: left wet at
        // rest, right dry. With the two-rarefaction estimate S_L = -c,
        // S_R = +2c the HLL star flux is F_mass = 2·c·h/3 > 0 (mass
        // moves into the dry region) and must be finite.
        let h = 1.0;
        let ul = Conserved::new(h, 0.0);
        let f = hll_flux(ul, Conserved::DRY);
        let c = (GRAVITY * h).sqrt();
        assert_relative_eq!(f.mass, 2.0 * c * h / 3.0, epsilon = 1e-12);
        assert!(f.momentum.is_finite());
    }

    #[test]
    fn near_dry_state_with_residual_momentum_does_not_blow_up() {
        // h below DRY_TOL with leftover hu is a numerical pathology the
        // update clamp normally removes; the Riemann solver must not
        // amplify it through hu²/h if it ever sees one.
        let pathological = Conserved::new(1.0e-13, 1.0e-13);
        let wet = Conserved::new(1.0, 0.0);
        for f in [hll_flux(wet, pathological), hll_flux(pathological, wet)] {
            assert!(f.mass.is_finite());
            assert!(f.momentum.is_finite());
        }
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
