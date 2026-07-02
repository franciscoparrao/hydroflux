//! Physical flux for the 1D Saint-Venant system in conservation form.
//!
//! `F(U) = (hu, hu² + g h²/2)`.

use crate::GRAVITY;
use crate::state::Conserved;
use hydroflux_autograd::Real;

/// 1D Saint-Venant flux vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Flux<T = f64> {
    /// Mass flux component `hu`.
    pub mass: T,
    /// Momentum flux component `hu² + g h²/2`.
    pub momentum: T,
}

impl<T: Real> Flux<T> {
    /// Compute `F(U)` for a conservative state. Returns the zero flux when
    /// the cell is dry; the momentum term is otherwise written to avoid
    /// `0/0` if the caller passes `h = 0` with non-zero `hu` (a numerical
    /// pathology, not a physical state).
    pub fn from_state(u: Conserved<T>) -> Self {
        let h = u.h.max(T::zero());
        if h.value() == 0.0 {
            return Self::zero();
        }
        let hu = u.hu;
        Self {
            mass: hu,
            momentum: hu * hu / h + h * h * (0.5 * GRAVITY),
        }
    }

    /// Zero flux. Used as the trivial value for dry-dry interfaces.
    pub fn zero() -> Self {
        Self {
            mass: T::zero(),
            momentum: T::zero(),
        }
    }
}

impl Flux {
    /// Zero flux constant for the `f64` configuration (generic code uses
    /// [`Flux::zero`]).
    pub const ZERO: Self = Self {
        mass: 0.0,
        momentum: 0.0,
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn zero_velocity_gives_hydrostatic_momentum() {
        let h = 2.0;
        let f = Flux::from_state(Conserved::new(h, 0.0));
        assert_relative_eq!(f.mass, 0.0, epsilon = 1e-12);
        assert_relative_eq!(f.momentum, 0.5 * GRAVITY * h * h, epsilon = 1e-12);
    }

    #[test]
    fn dry_cell_yields_zero_flux() {
        let f = Flux::from_state(Conserved::DRY);
        assert_eq!(f, Flux::ZERO);
    }

    #[test]
    fn flux_matches_textbook_formula() {
        // U = (h, hu) with u = 1.5, h = 2.0
        let h = 2.0;
        let u = 1.5;
        let hu = h * u;
        let f = Flux::from_state(Conserved::new(h, hu));
        assert_relative_eq!(f.mass, hu, epsilon = 1e-12);
        assert_relative_eq!(f.momentum, hu * u + 0.5 * GRAVITY * h * h, epsilon = 1e-12);
    }
}
