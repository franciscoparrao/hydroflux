//! Boundary condition types and ghost-cell construction.
//!
//! Boundary conditions are imposed via a single ghost cell at each end of
//! the domain. The numerical flux at the boundary face is then the standard
//! HLL flux between the ghost cell and the first inner cell.

use crate::state::Conserved;

/// Boundary kind at one end of the 1D domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Zero-gradient outflow: ghost cell equals the adjacent inner cell.
    /// Waves leave the domain (approximately) without reflecting.
    Transmissive,
    /// Reflective wall: ghost mirrors the inner cell with reversed velocity.
    /// Mass flux through the boundary is exactly zero by symmetry of the
    /// HLL solver.
    Wall,
}

/// Pair of boundary conditions for the 1D domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

/// Construct the ghost cell adjacent to an inner cell according to a
/// boundary kind. The ghost is placed *outside* the domain; flux is computed
/// between the ghost and the inner cell in the usual Riemann fashion.
pub fn ghost_state(inner: Conserved, kind: Boundary) -> Conserved {
    match kind {
        Boundary::Transmissive => inner,
        Boundary::Wall => Conserved {
            h: inner.h,
            hu: -inner.hu,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transmissive_ghost_equals_inner() {
        let u = Conserved::new(1.5, 0.6);
        assert_eq!(ghost_state(u, Boundary::Transmissive), u);
    }

    #[test]
    fn wall_ghost_mirrors_velocity_and_keeps_depth() {
        let u = Conserved::new(1.5, 0.6);
        let g = ghost_state(u, Boundary::Wall);
        assert_eq!(g.h, u.h);
        assert_eq!(g.hu, -u.hu);
    }

    #[test]
    fn wall_at_rest_is_self_ghost() {
        // At rest (hu = 0) the wall and transmissive ghosts coincide.
        let u = Conserved::new(2.0, 0.0);
        assert_eq!(
            ghost_state(u, Boundary::Wall),
            ghost_state(u, Boundary::Transmissive)
        );
    }
}
