//! State variables for the 1D Saint-Venant system.
//!
//! The conserved variables are `U = (h, hu)`, with `h` the water depth and
//! `hu` the discharge per unit width. Primitives `(h, u)` are kept for the
//! API boundary where users reason in physical terms.

use crate::GRAVITY;

/// Conservative state of a single finite-volume cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conserved {
    /// Water depth [m]. Non-negative; the solver enforces the wet/dry
    /// threshold at the FV update level, not at the type level.
    pub h: f64,
    /// Discharge per unit width [m²/s], equal to `h * u`.
    pub hu: f64,
}

impl Conserved {
    /// New conserved state. No validation: invariants are the solver's job.
    pub const fn new(h: f64, hu: f64) -> Self {
        Self { h, hu }
    }

    /// Dry cell: zero depth, zero discharge.
    pub const DRY: Self = Self { h: 0.0, hu: 0.0 };

    /// Convert to primitives. Returns zero velocity for dry cells (depth
    /// below `dry_tol`) to avoid `0/0`.
    pub fn to_primitive(self, dry_tol: f64) -> Primitive {
        let u = if self.h > dry_tol {
            self.hu / self.h
        } else {
            0.0
        };
        Primitive { h: self.h, u }
    }

    /// Gravity wave celerity `sqrt(g h)` for this cell.
    pub fn celerity(self) -> f64 {
        (GRAVITY * self.h.max(0.0)).sqrt()
    }
}

/// Primitive state: depth and velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Primitive {
    /// Water depth [m].
    pub h: f64,
    /// Depth-averaged velocity [m/s].
    pub u: f64,
}

impl Primitive {
    /// New primitive state.
    pub const fn new(h: f64, u: f64) -> Self {
        Self { h, u }
    }

    /// Convert to conserved variables.
    pub fn to_conserved(self) -> Conserved {
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
}
