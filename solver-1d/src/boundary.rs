//! Boundary condition types and ghost-cell construction.
//!
//! Boundary conditions are imposed via a single ghost cell at each end of
//! the domain. The numerical flux at the boundary face is then the standard
//! Audusse well-balanced HLL flux between the ghost cell and the first
//! inner cell.
//!
//! # Bed elevation at the ghost
//!
//! "Computational" BCs (`Transmissive`, `Wall`) extend the bed elevation
//! as zero-gradient (`z_ghost = z_inner`): no bed jump at the boundary
//! face, so no source correction is generated there. This is benign for
//! lake-at-rest tests but does **not** preserve uniform flow on a slope —
//! cell 0 misses the Audusse bed-slope source it would receive from an
//! interior left neighbour, and a boundary-layer perturbation propagates
//! into the domain (see `benchmarks/macdonald-uniform-results.md`).
//!
//! "Physical" BCs (`Discharge`, `Depth`) extend the bed linearly across
//! the boundary, `z_ghost = 2·z_inner − z_next_inner`, so the boundary
//! face carries the same bed jump as an interior face. Combined with the
//! prescribed flow variable, this preserves analytical steady states
//! across the **entire** domain rather than only the interior slab.

use crate::geometry::Channel1D;
use crate::state::Conserved;
use hydroflux_autograd::Real;

/// Boundary kind at one end of the 1D domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Boundary {
    /// Zero-gradient outflow: ghost equals the adjacent inner cell.
    /// Waves leave the domain (approximately) without reflecting.
    Transmissive,
    /// Reflective wall: ghost mirrors the inner cell with reversed
    /// velocity. Mass flux through the boundary is exactly zero by
    /// symmetry of the HLL flux.
    Wall,
    /// Prescribed unit discharge `hu = q` at the boundary. Depth is
    /// extrapolated as zero-gradient (`h_ghost = h_inner`). Intended for
    /// sub-critical inflow at the upstream end; the sign convention is
    /// `q > 0` for flow in the `+x` direction regardless of side. Bed is
    /// extended linearly so the boundary face carries the interior bed
    /// jump.
    Discharge {
        /// Prescribed unit discharge `hu` [m²/s] imposed in the ghost cell.
        q: f64,
    },
    /// Prescribed depth `h` at the boundary. Discharge is extrapolated
    /// as zero-gradient (`hu_ghost = hu_inner`). Intended for sub-critical
    /// outflow at the downstream end. Bed is extended linearly.
    Depth {
        /// Prescribed water depth [m] imposed in the ghost cell.
        h: f64,
    },
}

/// Which end of the domain a boundary lives on. Used by [`ghost_cell`] to
/// pick the correct inner / next-inner cells for state and bed
/// extrapolation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    /// Upstream side, face 0.
    Left,
    /// Downstream side, face `n`.
    Right,
}

/// Pair of boundary conditions for the 1D domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boundaries {
    /// Left (upstream-side) boundary.
    pub left: Boundary,
    /// Right (downstream-side) boundary.
    pub right: Boundary,
}

impl Boundaries {
    /// Both ends transmissive (outflow). Default for "open" prototype runs.
    pub const TRANSMISSIVE: Self = Self {
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
    };

    /// Both ends walls (closed box). Useful for mass-conservation tests.
    pub const WALLS: Self = Self {
        left: Boundary::Wall,
        right: Boundary::Wall,
    };
}

/// Build the ghost cell at one side of the domain: returns the conservative
/// state `(h, hu)` and the bed elevation `z` of the ghost. Both are needed
/// by the Audusse hydrostatic reconstruction at the boundary face.
///
/// The `inner` parameter is the adjacent inner cell state (`states[0]` for
/// `Side::Left`, `states[n-1]` for `Side::Right`).
///
/// The prescribed BC values (`q`, `h`) are `f64` and enter the generic
/// state as constants (`T::from_f64`): under AD they carry zero
/// derivative, which is the "gradient with the BC held fixed" semantics
/// the calibration workflow (and its FD verification) relies on.
pub fn ghost_cell<T: Real>(
    channel: &Channel1D<T>,
    inner: Conserved<T>,
    kind: Boundary,
    side: Side,
) -> (Conserved<T>, f64) {
    let state = match kind {
        Boundary::Transmissive => inner,
        Boundary::Wall => Conserved {
            h: inner.h,
            hu: -inner.hu,
        },
        Boundary::Discharge { q } => Conserved {
            h: inner.h,
            hu: T::from_f64(q),
        },
        Boundary::Depth { h } => Conserved {
            h: T::from_f64(h),
            hu: inner.hu,
        },
    };
    let bed = ghost_bed(channel, kind, side);
    (state, bed)
}

/// Bed elevation of the ghost cell. Computational BCs use zero-gradient
/// extrapolation; physical BCs extend the bed linearly so the boundary
/// face carries the same bed jump as interior faces.
fn ghost_bed<T: Real>(channel: &Channel1D<T>, kind: Boundary, side: Side) -> f64 {
    let n = channel.n_cells();
    match kind {
        Boundary::Transmissive | Boundary::Wall => match side {
            Side::Left => channel.bed[0],
            Side::Right => channel.bed[n - 1],
        },
        Boundary::Discharge { .. } | Boundary::Depth { .. } => match side {
            // Linear extrapolation: z_ghost = z_inner + (z_inner − z_next).
            Side::Left => 2.0 * channel.bed[0] - channel.bed[1],
            Side::Right => 2.0 * channel.bed[n - 1] - channel.bed[n - 2],
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    fn flat_channel() -> Channel1D {
        Channel1D::new(array![10.0, 10.0, 10.0, 10.0], 1.0, 0.03)
    }

    fn sloped_channel() -> Channel1D {
        // Descending bed: z = -i·dx·S₀ with dx=1, S₀=0.01.
        Channel1D::new(array![0.0, -0.01, -0.02, -0.03], 1.0, 0.03)
    }

    #[test]
    fn transmissive_ghost_equals_inner() {
        let ch = flat_channel();
        let u = Conserved::new(1.5, 0.6);
        let (g, _) = ghost_cell(&ch, u, Boundary::Transmissive, Side::Left);
        assert_eq!(g, u);
    }

    #[test]
    fn wall_ghost_mirrors_velocity_and_keeps_depth() {
        let ch = flat_channel();
        let u = Conserved::new(1.5, 0.6);
        let (g, _) = ghost_cell(&ch, u, Boundary::Wall, Side::Left);
        assert_eq!(g.h, u.h);
        assert_eq!(g.hu, -u.hu);
    }

    #[test]
    fn wall_at_rest_coincides_with_transmissive() {
        let ch = flat_channel();
        let u = Conserved::new(2.0, 0.0);
        let (g_wall, _) = ghost_cell(&ch, u, Boundary::Wall, Side::Left);
        let (g_trans, _) = ghost_cell(&ch, u, Boundary::Transmissive, Side::Left);
        assert_eq!(g_wall, g_trans);
    }

    #[test]
    fn discharge_prescribes_hu_keeps_h_zero_gradient() {
        let ch = sloped_channel();
        let inner = Conserved::new(0.6, 1.2); // arbitrary inner state
        let q_in = 0.8;
        let (g, _) = ghost_cell(&ch, inner, Boundary::Discharge { q: q_in }, Side::Left);
        assert_relative_eq!(g.h, inner.h, epsilon = 1e-12);
        assert_relative_eq!(g.hu, q_in, epsilon = 1e-12);
    }

    #[test]
    fn depth_prescribes_h_keeps_hu_zero_gradient() {
        let ch = sloped_channel();
        let inner = Conserved::new(0.6, 1.2);
        let h_out = 0.45;
        let (g, _) = ghost_cell(&ch, inner, Boundary::Depth { h: h_out }, Side::Right);
        assert_relative_eq!(g.h, h_out, epsilon = 1e-12);
        assert_relative_eq!(g.hu, inner.hu, epsilon = 1e-12);
    }

    #[test]
    fn computational_bcs_use_zero_gradient_bed() {
        let ch = sloped_channel();
        let inner = Conserved::new(0.5, 0.0);
        let (_, z_left) = ghost_cell(&ch, inner, Boundary::Transmissive, Side::Left);
        let (_, z_right) = ghost_cell(&ch, inner, Boundary::Wall, Side::Right);
        assert_eq!(z_left, ch.bed[0]);
        assert_eq!(z_right, ch.bed[ch.n_cells() - 1]);
    }

    #[test]
    fn physical_bcs_extend_bed_linearly() {
        let ch = sloped_channel(); // bed = [0, -0.01, -0.02, -0.03], slope 0.01
        let inner_left = Conserved::new(0.5, 0.0);
        let inner_right = Conserved::new(0.5, 0.0);

        let (_, z_left_disch) =
            ghost_cell(&ch, inner_left, Boundary::Discharge { q: 1.0 }, Side::Left);
        // Linear extrapolation upstream: 2·bed[0] − bed[1] = 2·0 − (−0.01) = +0.01.
        assert_relative_eq!(z_left_disch, 0.01, epsilon = 1e-12);

        let (_, z_right_depth) =
            ghost_cell(&ch, inner_right, Boundary::Depth { h: 0.5 }, Side::Right);
        // Linear extrapolation downstream: 2·bed[3] − bed[2] = 2·(−0.03) − (−0.02) = −0.04.
        assert_relative_eq!(z_right_depth, -0.04, epsilon = 1e-12);
    }
}
