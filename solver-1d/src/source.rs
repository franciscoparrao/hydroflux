//! Source-term steps applied as operator-split fractional updates.
//!
//! The bed-slope source is folded into the well-balanced flux inside
//! [`crate::update`] (Audusse hydrostatic reconstruction). This module
//! covers the remaining cell-local sources — for now, just Manning friction.

use crate::GRAVITY;
use crate::state::Conserved;

/// Semi-implicit Manning friction step. In-place point-implicit update:
///
/// ```text
///   hu^{n+1} = hu^n / (1 + dt · g n² |u^n| / h^{4/3})
/// ```
///
/// Properties:
/// - Unconditionally stable (no extra CFL beyond the FV step).
/// - Strictly decreases `|hu|`; preserves `hu = 0` exactly.
/// - Cells with `h ≤ dry_tol` are skipped (friction is undefined as `h → 0`).
/// - `manning == 0` short-circuits to a no-op.
///
/// `manning` is the Manning roughness coefficient `n` (units s/m^(1/3));
/// `dry_tol` is the wet/dry threshold in metres.
pub fn manning_friction_step(states: &mut [Conserved], manning: f64, dt: f64, dry_tol: f64) {
    if manning == 0.0 {
        return;
    }
    let n_sq = manning * manning;
    for s in states.iter_mut() {
        if s.h <= dry_tol {
            continue;
        }
        let u = s.hu / s.h;
        let factor = 1.0 + dt * GRAVITY * n_sq * u.abs() / s.h.powf(4.0 / 3.0);
        s.hu /= factor;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn zero_manning_is_identity() {
        let mut states = vec![Conserved::new(1.0, 2.0), Conserved::new(2.0, -1.5)];
        let before = states.clone();
        manning_friction_step(&mut states, 0.0, 0.1, 1e-9);
        assert_eq!(states, before);
    }

    #[test]
    fn rest_state_is_preserved() {
        let mut states = vec![Conserved::new(2.0, 0.0); 5];
        manning_friction_step(&mut states, 0.03, 0.1, 1e-9);
        for s in &states {
            assert_eq!(s.hu, 0.0);
            assert_eq!(s.h, 2.0);
        }
    }

    #[test]
    fn dry_cells_are_untouched() {
        let mut states = vec![Conserved::DRY, Conserved::new(0.0, 1e-30)];
        let before = states.clone();
        manning_friction_step(&mut states, 0.03, 0.1, 1e-9);
        assert_eq!(states, before);
    }

    #[test]
    fn moving_water_decelerates_monotonically() {
        // No driver, only friction. |hu| must decrease strictly each step
        // and never overshoot through zero. The semi-implicit factor is
        // unconditionally > 1 for |u| > 0, so this is exact algebra, not
        // a numerical coincidence — checking it once is enough to detect
        // sign-flip or formula bugs.
        let mut states = vec![Conserved::new(1.0, 2.0)];
        let initial = states[0].hu;
        let mut previous = initial.abs();
        for _ in 0..50 {
            manning_friction_step(&mut states, 0.05, 0.1, 1e-9);
            let current = states[0].hu.abs();
            assert!(
                current < previous,
                "|hu| did not decrease: {previous} → {current}"
            );
            assert!(
                states[0].hu * initial >= 0.0,
                "hu overshot through zero: {}",
                states[0].hu
            );
            previous = current;
        }
    }

    #[test]
    fn one_step_matches_analytic_for_known_state() {
        // Verify the closed-form: hu^{n+1} = hu / (1 + dt g n² |u| / h^{4/3}).
        let mut s = vec![Conserved::new(1.0, 2.0)];
        let manning = 0.04_f64;
        let dt = 0.05_f64;
        let h = 1.0_f64;
        let hu = 2.0_f64;
        let u = hu / h;
        let expected = hu / (1.0 + dt * GRAVITY * manning * manning * u.abs() / h.powf(4.0 / 3.0));
        manning_friction_step(&mut s, manning, dt, 1e-9);
        assert_relative_eq!(s[0].hu, expected, epsilon = 1e-12);
        assert_relative_eq!(s[0].h, h, epsilon = 1e-12);
    }
}
