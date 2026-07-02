//! State variables for the 1D Saint-Venant system.
//!
//! The conserved variables are `U = (h, hu)`, with `h` the water depth and
//! `hu` the discharge per unit width. Primitives `(h, u)` are kept for the
//! API boundary where users reason in physical terms.
//!
//! All state types are generic over `T: Real` (defaulting to `f64`), so
//! the same solver code path runs the `f64` production configuration and
//! the `Dual` forward-mode AD configuration — the pattern proven in
//! `hydroflux-solver-2d`.

use crate::GRAVITY;
use hydroflux_autograd::Real;

/// Conservative state of a single finite-volume cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conserved<T = f64> {
    /// Water depth [m]. Non-negative; the solver enforces the wet/dry
    /// threshold at the FV update level, not at the type level.
    pub h: T,
    /// Discharge per unit width [m²/s], equal to `h * u`.
    pub hu: T,
}

impl<T: Real> Conserved<T> {
    /// New conserved state. No validation: invariants are the solver's job.
    pub fn new(h: T, hu: T) -> Self {
        Self { h, hu }
    }

    /// Dry cell: zero depth, zero discharge.
    pub fn dry() -> Self {
        Self {
            h: T::zero(),
            hu: T::zero(),
        }
    }

    /// Convert to primitives. Returns zero velocity for dry cells (depth
    /// below `dry_tol`) to avoid `0/0`.
    pub fn to_primitive(self, dry_tol: f64) -> Primitive<T> {
        let u = if self.h.value() > dry_tol {
            self.hu / self.h
        } else {
            T::zero()
        };
        Primitive { h: self.h, u }
    }

    /// Gravity wave celerity `sqrt(g h)` for this cell.
    pub fn celerity(self) -> T {
        (self.h.max(T::zero()) * GRAVITY).sqrt()
    }
}

impl Conserved {
    /// Dry cell constant for the `f64` configuration (generic code uses
    /// [`Conserved::dry`]).
    pub const DRY: Self = Self { h: 0.0, hu: 0.0 };
}

/// Primitive state: depth and velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Primitive<T = f64> {
    /// Water depth [m].
    pub h: T,
    /// Depth-averaged velocity [m/s].
    pub u: T,
}

impl<T: Real> Primitive<T> {
    /// New primitive state.
    pub fn new(h: T, u: T) -> Self {
        Self { h, u }
    }

    /// Convert to conserved variables.
    pub fn to_conserved(self) -> Conserved<T> {
        Conserved {
            h: self.h,
            hu: self.h * self.u,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn primitive_conserved_roundtrip() {
        let p = Primitive::new(1.5, 2.3);
        let c = p.to_conserved();
        let p2 = c.to_primitive(1e-9);
        assert_relative_eq!(p.h, p2.h, epsilon = 1e-12);
        assert_relative_eq!(p.u, p2.u, epsilon = 1e-12);
    }

    #[test]
    fn dry_cell_has_zero_velocity_under_tol() {
        let c = Conserved::new(1e-10, 1e-20);
        let p = c.to_primitive(1e-9);
        assert_eq!(p.u, 0.0);
    }

    #[test]
    fn celerity_matches_sqrt_g_h() {
        let c = Conserved::new(4.0, 0.0);
        assert_relative_eq!(c.celerity(), (GRAVITY * 4.0).sqrt(), epsilon = 1e-12);
    }

    #[test]
    fn dual_constant_state_matches_f64() {
        use hydroflux_autograd::Dual;
        let cf = Conserved::new(2.0, 3.0);
        let cd = Conserved::new(Dual::constant(2.0), Dual::constant(3.0));
        assert_eq!(cf.celerity(), cd.celerity().val);
        assert_eq!(cd.celerity().dval, 0.0);
    }
}
