//! Source-term steps applied as operator-split fractional updates.
//!
//! The bed-slope source is folded into the well-balanced flux inside
//! [`crate::update`] (Audusse hydrostatic reconstruction in 2D). This
//! module covers the remaining cell-local sources — for now, just
//! Manning friction.

use crate::GRAVITY;
use crate::state::Conserved2D;
use ndarray::Array2;

/// Semi-implicit Manning friction step in 2D. In-place point-implicit
/// update applied independently in each cell:
///
/// ```text
///   |U|^n      = √((hu^n)² + (hv^n)²) / h^n
///   α          = dt · g · n² · |U|^n / h^{4/3}
///   hu^{n+1}   = hu^n / (1 + α)
///   hv^{n+1}   = hv^n / (1 + α)
/// ```
///
/// The factor `1 + α` is shared by both momentum components because
/// the friction vector is `g n² |U| U / h^{4/3}` — parallel to the
/// velocity. As a consequence the **flow direction is preserved
/// exactly** by this step; only the magnitude decays.
///
/// Properties:
/// - Unconditionally stable: `1 + α ≥ 1` for any `dt > 0`, so
///   `|hu|`, `|hv|` and `|U|` are non-increasing. No extra CFL beyond
///   the FV step.
/// - Preserves rest (`hu = hv = 0`) and dry cells exactly.
/// - `manning == 0` short-circuits to a no-op.
///
/// `manning` is the Manning roughness coefficient `n` (units
/// s/m^(1/3)); `dry_tol` is the wet/dry threshold in metres.
pub fn manning_friction_step(
    states: &mut Array2<Conserved2D>,
    manning: f64,
    dt: f64,
    dry_tol: f64,
) {
    if manning == 0.0 {
        return;
    }
    let n_sq = manning * manning;
    for s in states.iter_mut() {
        if s.h <= dry_tol {
            continue;
        }
        // |U| = √(u² + v²) = √((hu)² + (hv)²) / h. Compute directly to
        // avoid two divisions in the hot path.
        let speed = (s.hu * s.hu + s.hv * s.hv).sqrt() / s.h;
        let factor = 1.0 + dt * GRAVITY * n_sq * speed / s.h.powf(4.0 / 3.0);
        s.hu /= factor;
        s.hv /= factor;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn zero_manning_is_identity() {
        let mut states = array![[
            Conserved2D::new(1.0, 2.0, -0.5),
            Conserved2D::new(2.0, -1.5, 0.7)
        ]];
        let before = states.clone();
        manning_friction_step(&mut states, 0.0, 0.1, 1e-9);
        assert_eq!(states, before);
    }

    #[test]
    fn rest_state_is_preserved_exactly() {
        let mut states = Array2::from_elem((3, 4), Conserved2D::new(2.0, 0.0, 0.0));
        manning_friction_step(&mut states, 0.03, 0.1, 1e-9);
        for s in &states {
            assert_eq!(s.hu, 0.0);
            assert_eq!(s.hv, 0.0);
            assert_eq!(s.h, 2.0);
        }
    }

    #[test]
    fn dry_cells_are_untouched() {
        let mut states = array![[Conserved2D::DRY, Conserved2D::new(0.0, 1e-30, -1e-30)]];
        let before = states.clone();
        manning_friction_step(&mut states, 0.03, 0.1, 1e-9);
        assert_eq!(states, before);
    }

    #[test]
    fn moving_water_decelerates_monotonically_in_speed() {
        // No driver, only friction. The Euclidean speed |U| must
        // decrease strictly each step and never overshoot through zero.
        // Algebraic property of the semi-implicit factor (1 + α) ≥ 1;
        // catches sign-flip and formula bugs.
        let mut states = array![[Conserved2D::new(1.0, 2.0, -1.5)]];
        let initial_speed = {
            let s = states[(0, 0)];
            (s.hu * s.hu + s.hv * s.hv).sqrt() / s.h
        };
        let mut previous = initial_speed;
        for _ in 0..50 {
            manning_friction_step(&mut states, 0.05, 0.1, 1e-9);
            let s = states[(0, 0)];
            let current = (s.hu * s.hu + s.hv * s.hv).sqrt() / s.h;
            assert!(
                current < previous,
                "|U| did not decrease: {previous} → {current}"
            );
            assert!(current >= 0.0, "speed went negative: {current}");
            previous = current;
        }
    }

    #[test]
    fn flow_direction_is_preserved_exactly() {
        // Both momentum components are divided by the SAME factor, so
        // the angle atan2(hv, hu) must be invariant under the friction
        // step. This is the key property that distinguishes a proper
        // 2D Manning step (vector friction) from component-wise scalar
        // friction (would skew the direction toward the smaller axis).
        let mut states = array![[Conserved2D::new(1.0, 2.0, 1.0)]];
        let initial = states[(0, 0)];
        let initial_angle = initial.hv.atan2(initial.hu);
        for _ in 0..50 {
            manning_friction_step(&mut states, 0.05, 0.1, 1e-9);
            let s = states[(0, 0)];
            let angle = s.hv.atan2(s.hu);
            assert_relative_eq!(angle, initial_angle, epsilon = 1e-12);
        }
    }

    #[test]
    fn one_step_matches_analytic_for_known_state() {
        // Verify the closed-form against hand-computed values for a
        // single cell with non-trivial (u, v).
        let mut s = array![[Conserved2D::new(1.0, 2.0, 1.0)]];
        let manning = 0.04_f64;
        let dt = 0.05_f64;
        let h = 1.0_f64;
        let hu = 2.0_f64;
        let hv = 1.0_f64;
        let speed = (hu * hu + hv * hv).sqrt() / h;
        let factor = 1.0 + dt * GRAVITY * manning * manning * speed / h.powf(4.0 / 3.0);
        let expected_hu = hu / factor;
        let expected_hv = hv / factor;
        manning_friction_step(&mut s, manning, dt, 1e-9);
        assert_relative_eq!(s[(0, 0)].hu, expected_hu, epsilon = 1e-12);
        assert_relative_eq!(s[(0, 0)].hv, expected_hv, epsilon = 1e-12);
        assert_relative_eq!(s[(0, 0)].h, h, epsilon = 1e-12);
    }

    #[test]
    fn pure_x_flow_keeps_hv_at_zero() {
        // If hv starts at 0, friction must keep it at 0 exactly (no
        // spurious cross-coupling). The shared factor ensures
        // hv^{n+1} = 0 / (1 + α) = 0.
        let mut states = array![[Conserved2D::new(1.5, 3.0, 0.0)]];
        for _ in 0..20 {
            manning_friction_step(&mut states, 0.04, 0.05, 1e-9);
            assert_eq!(states[(0, 0)].hv, 0.0);
            assert!(states[(0, 0)].hu > 0.0);
        }
    }

    #[test]
    fn pure_y_flow_keeps_hu_at_zero() {
        let mut states = array![[Conserved2D::new(1.5, 0.0, -3.0)]];
        for _ in 0..20 {
            manning_friction_step(&mut states, 0.04, 0.05, 1e-9);
            assert_eq!(states[(0, 0)].hu, 0.0);
            assert!(states[(0, 0)].hv < 0.0);
        }
    }

    #[test]
    fn depth_is_unchanged_by_friction() {
        // Friction is a momentum-only source — mass must be invariant
        // cell-by-cell.
        let mut states = array![
            [
                Conserved2D::new(1.0, 2.0, 0.5),
                Conserved2D::new(2.0, -1.0, 1.5)
            ],
            [
                Conserved2D::new(0.5, 0.3, -0.3),
                Conserved2D::new(3.0, 0.0, 0.0)
            ],
        ];
        let depths_before: Vec<f64> = states.iter().map(|s| s.h).collect();
        manning_friction_step(&mut states, 0.035, 0.05, 1e-9);
        let depths_after: Vec<f64> = states.iter().map(|s| s.h).collect();
        for (b, a) in depths_before.iter().zip(depths_after.iter()) {
            assert_eq!(b, a);
        }
    }
}
