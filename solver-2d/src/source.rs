//! Source-term steps applied as operator-split fractional updates.
//!
//! The bed-slope source is folded into the well-balanced flux inside
//! [`crate::update`] (Audusse hydrostatic reconstruction in 2D). This
//! module covers the remaining cell-local sources:
//!
//! - [`manning_friction_step`]: semi-implicit Manning friction.
//! - [`apply_point_sources`]: point inflows (e.g. UK EA Test 1).

use crate::GRAVITY;
use crate::geometry::Mesh2D;
use crate::state::Conserved2D;
use ndarray::Array2;

/// A point inflow source: adds water mass to a single cell at a
/// prescribed volumetric rate. Used for UK EA-style tests with
/// localised inflow (e.g. "flooding a disconnected water body").
///
/// The source does **not** add momentum: the injected water enters
/// at rest and is accelerated downhill by the bed-slope source +
/// pressure gradient. This matches the physical picture of a
/// vertical drop (pipe outlet, rainfall focused at a point) and is
/// what the UK EA tests assume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointSource {
    /// Row index `i` (along `y`) of the cell receiving the inflow.
    pub row: usize,
    /// Column index `j` (along `x`).
    pub col: usize,
    /// Mass inflow rate `Q` [m³/s]. Positive injects mass; negative
    /// withdraws it (useful for sinks like infiltration patches).
    pub q_mass: f64,
}

/// Apply all point sources to the state for a single timestep.
///
/// For each source, `h_cell += Q · dt / (Δx · Δy)`. Momentum is
/// unchanged. Cells that go below zero from a negative source
/// (sink) are clamped to zero — sinks cannot draw more water than
/// the cell holds.
pub fn apply_point_sources(
    states: &mut Array2<Conserved2D>,
    sources: &[PointSource],
    dt: f64,
    dx: f64,
    dy: f64,
) {
    let cell_area = dx * dy;
    for src in sources {
        let dh = src.q_mass * dt / cell_area;
        let cell = &mut states[(src.row, src.col)];
        cell.h = (cell.h + dh).max(0.0);
        if cell.h == 0.0 {
            cell.hu = 0.0;
            cell.hv = 0.0;
        }
    }
}

/// Apply uniform rainfall to every cell for one timestep.
///
/// `rate` is the precipitation intensity in **metres of water per
/// second** (i.e. depth per unit time, NOT per unit area). To
/// convert from common units:
/// - mm/hour:    `rate = mm_per_hour · 1e-3 / 3600`
/// - mm/minute:  `rate = mm_per_minute · 1e-3 / 60`
/// - m³/(s·m²):  same value (depth per area-time is the same number)
///
/// Negative `rate` represents evaporation / infiltration: a depth
/// per unit time being removed. Cells whose depth would go below
/// zero are clamped to dry (depth + both momentum components set
/// to zero).
///
/// The rainfall does NOT add momentum (the convention from UK EA
/// Test 2 and most flood solvers): rain falls vertically, hits
/// stationary, and is then accelerated by the bed-slope source.
pub fn apply_rain(states: &mut Array2<Conserved2D>, rate: f64, dt: f64) {
    if rate == 0.0 {
        return;
    }
    let dh = rate * dt;
    for cell in states.iter_mut() {
        let new_h = cell.h + dh;
        if new_h <= 0.0 {
            *cell = Conserved2D::DRY;
        } else {
            cell.h = new_h;
        }
    }
}

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
/// - Cells with `n = 0` (Manning field zero locally) are a no-op.
///
/// The Manning roughness `n` is read per cell from `mesh.manning`;
/// `dry_tol` is the wet/dry threshold in metres.
pub fn manning_friction_step(
    states: &mut Array2<Conserved2D>,
    mesh: &Mesh2D,
    dt: f64,
    dry_tol: f64,
) {
    assert_eq!(
        states.dim(),
        mesh.manning.dim(),
        "states shape {:?} must match mesh.manning shape {:?}",
        states.dim(),
        mesh.manning.dim(),
    );
    for ((i, j), s) in states.indexed_iter_mut() {
        if s.h <= dry_tol {
            continue;
        }
        let n = mesh.manning[(i, j)];
        if n == 0.0 {
            continue;
        }
        // |U| = √(u² + v²) = √((hu)² + (hv)²) / h. Compute directly to
        // avoid two divisions in the hot path.
        let speed = (s.hu * s.hu + s.hv * s.hv).sqrt() / s.h;
        let factor = 1.0 + dt * GRAVITY * n * n * speed / s.h.powf(4.0 / 3.0);
        s.hu /= factor;
        s.hv /= factor;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    /// Flat-bed mesh with uniform Manning matching `states.dim()` —
    /// the simplest fixture for friction-step tests.
    fn mesh_for(states: &Array2<Conserved2D>, manning: f64) -> Mesh2D {
        let bed = Array2::<f64>::zeros(states.dim());
        Mesh2D::new(bed, 1.0, 1.0, manning)
    }

    #[test]
    fn zero_manning_is_identity() {
        let mut states = array![[
            Conserved2D::new(1.0, 2.0, -0.5),
            Conserved2D::new(2.0, -1.5, 0.7)
        ]];
        let before = states.clone();
        let mesh = mesh_for(&states, 0.0);
        manning_friction_step(&mut states, &mesh, 0.1, 1e-9);
        assert_eq!(states, before);
    }

    #[test]
    fn rest_state_is_preserved_exactly() {
        let mut states = Array2::from_elem((3, 4), Conserved2D::new(2.0, 0.0, 0.0));
        let mesh = mesh_for(&states, 0.03);
        manning_friction_step(&mut states, &mesh, 0.1, 1e-9);
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
        let mesh = mesh_for(&states, 0.03);
        manning_friction_step(&mut states, &mesh, 0.1, 1e-9);
        assert_eq!(states, before);
    }

    #[test]
    fn moving_water_decelerates_monotonically_in_speed() {
        // No driver, only friction. The Euclidean speed |U| must
        // decrease strictly each step and never overshoot through zero.
        // Algebraic property of the semi-implicit factor (1 + α) ≥ 1;
        // catches sign-flip and formula bugs.
        let mut states = array![[Conserved2D::new(1.0, 2.0, -1.5)]];
        let mesh = mesh_for(&states, 0.05);
        let initial_speed = {
            let s = states[(0, 0)];
            (s.hu * s.hu + s.hv * s.hv).sqrt() / s.h
        };
        let mut previous = initial_speed;
        for _ in 0..50 {
            manning_friction_step(&mut states, &mesh, 0.1, 1e-9);
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
        let mesh = mesh_for(&states, 0.05);
        let initial = states[(0, 0)];
        let initial_angle = initial.hv.atan2(initial.hu);
        for _ in 0..50 {
            manning_friction_step(&mut states, &mesh, 0.1, 1e-9);
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
        let mesh = mesh_for(&s, manning);
        let dt = 0.05_f64;
        let h = 1.0_f64;
        let hu = 2.0_f64;
        let hv = 1.0_f64;
        let speed = (hu * hu + hv * hv).sqrt() / h;
        let factor = 1.0 + dt * GRAVITY * manning * manning * speed / h.powf(4.0 / 3.0);
        let expected_hu = hu / factor;
        let expected_hv = hv / factor;
        manning_friction_step(&mut s, &mesh, dt, 1e-9);
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
        let mesh = mesh_for(&states, 0.04);
        for _ in 0..20 {
            manning_friction_step(&mut states, &mesh, 0.05, 1e-9);
            assert_eq!(states[(0, 0)].hv, 0.0);
            assert!(states[(0, 0)].hu > 0.0);
        }
    }

    #[test]
    fn pure_y_flow_keeps_hu_at_zero() {
        let mut states = array![[Conserved2D::new(1.5, 0.0, -3.0)]];
        let mesh = mesh_for(&states, 0.04);
        for _ in 0..20 {
            manning_friction_step(&mut states, &mesh, 0.05, 1e-9);
            assert_eq!(states[(0, 0)].hu, 0.0);
            assert!(states[(0, 0)].hv < 0.0);
        }
    }

    #[test]
    fn uniform_field_matches_scalar_bit_exact() {
        // `Mesh2D::with_manning_field(bed, dx, dy, uniform_array)` must
        // produce the SAME numerical result as `Mesh2D::new(bed, dx, dy,
        // scalar)` when the field is filled with the scalar. Roundoff-
        // exact equality after many steps catches any branching or
        // ordering divergence between the two construction paths.
        let n = 0.045_f64;
        let mut states_scalar = array![[
            Conserved2D::new(1.0, 2.0, -1.0),
            Conserved2D::new(2.0, 1.5, 0.3),
            Conserved2D::new(0.5, -0.7, 1.2),
        ]];
        let mut states_field = states_scalar.clone();
        let mesh_scalar = mesh_for(&states_scalar, n);
        let n_field = Array2::from_elem(states_field.dim(), n);
        let mesh_field = Mesh2D::with_manning_field(
            Array2::<f64>::zeros(states_field.dim()),
            1.0,
            1.0,
            n_field,
        );
        for _ in 0..50 {
            manning_friction_step(&mut states_scalar, &mesh_scalar, 0.05, 1e-9);
            manning_friction_step(&mut states_field, &mesh_field, 0.05, 1e-9);
        }
        for (a, b) in states_scalar.iter().zip(states_field.iter()) {
            assert_eq!(a.h, b.h);
            assert_eq!(a.hu, b.hu);
            assert_eq!(a.hv, b.hv);
        }
    }

    #[test]
    fn variable_field_decelerates_high_n_cell_faster() {
        // Two adjacent cells with identical initial state. Cell A has
        // Manning n_A = 0.10 (high roughness, like dense forest), cell
        // B has n_B = 0.02 (smooth concrete channel). After identical
        // dt the cell with higher n must have lost more momentum.
        // The relation is monotonic in n²·|U|, so a 5× difference in
        // n gives a 25× difference in deceleration ratio — easy to
        // detect even on a single step.
        let mut states = array![[
            Conserved2D::new(1.0, 3.0, 0.0),
            Conserved2D::new(1.0, 3.0, 0.0),
        ]];
        let n_field = ndarray::arr2(&[[0.10_f64, 0.02_f64]]);
        let mesh = Mesh2D::with_manning_field(
            Array2::<f64>::zeros((1, 2)),
            1.0,
            1.0,
            n_field,
        );
        manning_friction_step(&mut states, &mesh, 0.5, 1e-9);
        let hu_high_n = states[(0, 0)].hu;
        let hu_low_n = states[(0, 1)].hu;
        assert!(
            hu_high_n < hu_low_n,
            "high-n cell ({}) did not decelerate more than low-n cell ({})",
            hu_high_n,
            hu_low_n
        );
        // Sanity: both still positive (semi-implicit cannot overshoot
        // zero), and the low-n cell has barely changed.
        assert!(hu_high_n > 0.0);
        assert!(hu_low_n > 2.9);
    }

    #[test]
    #[should_panic(expected = "manning field shape")]
    fn with_manning_field_panics_on_shape_mismatch() {
        // `Mesh2D::with_manning_field` requires bed and manning to
        // have the same shape — programming error to mismatch.
        let _ = Mesh2D::with_manning_field(
            Array2::<f64>::zeros((3, 4)),
            1.0,
            1.0,
            Array2::<f64>::from_elem((3, 5), 0.03),
        );
    }

    #[test]
    #[should_panic(expected = "Manning n must be non-negative")]
    fn with_manning_field_panics_on_negative_value() {
        let mut n = Array2::<f64>::from_elem((2, 2), 0.03);
        n[(0, 0)] = -0.01;
        let _ = Mesh2D::with_manning_field(Array2::<f64>::zeros((2, 2)), 1.0, 1.0, n);
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
        let mesh = mesh_for(&states, 0.035);
        let depths_before: Vec<f64> = states.iter().map(|s| s.h).collect();
        manning_friction_step(&mut states, &mesh, 0.05, 1e-9);
        let depths_after: Vec<f64> = states.iter().map(|s| s.h).collect();
        for (b, a) in depths_before.iter().zip(depths_after.iter()) {
            assert_eq!(b, a);
        }
    }

    // -----------------------------------------------------------------
    // Point source tests
    // -----------------------------------------------------------------

    #[test]
    fn point_source_adds_mass_at_expected_rate() {
        // Q · dt / area is the depth increment per step.
        let mut states = Array2::from_elem((3, 3), Conserved2D::DRY);
        let sources = [PointSource {
            row: 1,
            col: 1,
            q_mass: 2.0,
        }];
        let dt = 0.5_f64;
        let dx = 1.0_f64;
        let dy = 1.0_f64;
        apply_point_sources(&mut states, &sources, dt, dx, dy);
        assert_relative_eq!(states[(1, 1)].h, 1.0, epsilon = 1e-12); // 2·0.5/1 = 1
        assert_eq!(states[(1, 1)].hu, 0.0);
        assert_eq!(states[(1, 1)].hv, 0.0);
        // Other cells untouched.
        for ((i, j), s) in states.indexed_iter() {
            if (i, j) != (1, 1) {
                assert_eq!(s.h, 0.0);
            }
        }
    }

    #[test]
    fn point_source_is_cumulative_across_steps() {
        let mut states = Array2::from_elem((1, 1), Conserved2D::DRY);
        let sources = [PointSource {
            row: 0,
            col: 0,
            q_mass: 1.0,
        }];
        for _ in 0..10 {
            apply_point_sources(&mut states, &sources, 0.1, 1.0, 1.0);
        }
        // 10 × (1 · 0.1 / 1) = 1.0 m
        assert_relative_eq!(states[(0, 0)].h, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn point_source_scales_inversely_with_cell_area() {
        let mut a = Array2::from_elem((1, 1), Conserved2D::DRY);
        let mut b = Array2::from_elem((1, 1), Conserved2D::DRY);
        let sources = [PointSource {
            row: 0,
            col: 0,
            q_mass: 1.0,
        }];
        apply_point_sources(&mut a, &sources, 1.0, 1.0, 1.0);
        apply_point_sources(&mut b, &sources, 1.0, 2.0, 2.0);
        // 4× the area → 1/4 the depth.
        assert_relative_eq!(a[(0, 0)].h, 4.0 * b[(0, 0)].h, epsilon = 1e-12);
    }

    #[test]
    fn negative_q_drains_to_zero_clamping_momentum() {
        // A sink (Q < 0) draws water out. If the sink takes more
        // than the cell holds, depth clamps to zero and momentum is
        // zeroed (cannot have moving water without water).
        let mut states = Array2::from_elem((1, 1), Conserved2D::new(1.0, 5.0, 0.0));
        let sources = [PointSource {
            row: 0,
            col: 0,
            q_mass: -10.0, // wants to remove 10 m³/s
        }];
        // dt = 1, area = 1 → would drain 10 m, but cell only has 1.
        apply_point_sources(&mut states, &sources, 1.0, 1.0, 1.0);
        assert_eq!(states[(0, 0)].h, 0.0);
        assert_eq!(states[(0, 0)].hu, 0.0);
        assert_eq!(states[(0, 0)].hv, 0.0);
    }

    #[test]
    fn multiple_sources_apply_independently() {
        let mut states = Array2::from_elem((2, 2), Conserved2D::DRY);
        let sources = [
            PointSource {
                row: 0,
                col: 0,
                q_mass: 1.0,
            },
            PointSource {
                row: 1,
                col: 1,
                q_mass: 2.0,
            },
        ];
        apply_point_sources(&mut states, &sources, 1.0, 1.0, 1.0);
        assert_relative_eq!(states[(0, 0)].h, 1.0, epsilon = 1e-12);
        assert_relative_eq!(states[(1, 1)].h, 2.0, epsilon = 1e-12);
        assert_eq!(states[(0, 1)].h, 0.0);
        assert_eq!(states[(1, 0)].h, 0.0);
    }

    // -----------------------------------------------------------------
    // Rain-on-grid tests
    // -----------------------------------------------------------------

    #[test]
    fn rain_adds_uniform_depth() {
        // rate · dt added to every cell. For dry cells, this is the
        // exact new depth.
        let mut states = Array2::from_elem((3, 4), Conserved2D::DRY);
        apply_rain(&mut states, 0.001, 60.0); // 1 mm/s · 60 s = 60 mm
        for s in &states {
            assert_relative_eq!(s.h, 0.06, epsilon = 1e-12);
            assert_eq!(s.hu, 0.0);
            assert_eq!(s.hv, 0.0);
        }
    }

    #[test]
    fn rain_preserves_existing_momentum() {
        // Rain on a wet cell adds depth but does not change `hu, hv`.
        // The new water enters at rest, so the cell's mean velocity
        // (`hu / h`) drops slightly, but the conserved momentum stays.
        let mut states = Array2::from_elem((1, 1), Conserved2D::new(1.0, 3.0, -1.5));
        apply_rain(&mut states, 0.005, 4.0); // adds 0.02 m
        assert_relative_eq!(states[(0, 0)].h, 1.02, epsilon = 1e-12);
        assert_eq!(states[(0, 0)].hu, 3.0);
        assert_eq!(states[(0, 0)].hv, -1.5);
    }

    #[test]
    fn negative_rain_evaporates_to_dry() {
        // Negative rate (evaporation) removes depth. A cell that drains
        // past zero is clamped DRY (depth + momentum both zero).
        let mut states = array![[Conserved2D::new(0.01, 0.5, 0.0)]];
        apply_rain(&mut states, -1.0, 1.0); // wants to remove 1 m
        assert_eq!(states[(0, 0)].h, 0.0);
        assert_eq!(states[(0, 0)].hu, 0.0);
        assert_eq!(states[(0, 0)].hv, 0.0);
    }

    #[test]
    fn zero_rain_is_identity() {
        let mut states = array![[Conserved2D::new(1.0, 0.5, -0.3)]];
        let before = states.clone();
        apply_rain(&mut states, 0.0, 100.0);
        assert_eq!(states, before);
    }

    #[test]
    fn rain_is_cumulative() {
        let mut states = Array2::from_elem((1, 1), Conserved2D::DRY);
        // 10 mm/h falling for 6 minutes = 1 mm total.
        let rate = 0.01_f64 / 3600.0;
        for _ in 0..360 {
            // 360 × 1 s = 6 min
            apply_rain(&mut states, rate, 1.0);
        }
        assert_relative_eq!(states[(0, 0)].h, 0.001, epsilon = 1e-12);
    }

    #[test]
    fn empty_sources_list_is_no_op() {
        let mut states = Array2::from_elem((2, 2), Conserved2D::new(1.0, 0.5, -0.3));
        let before = states.clone();
        apply_point_sources(&mut states, &[], 0.1, 1.0, 1.0);
        assert_eq!(states, before);
    }
}
