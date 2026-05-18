//! State variables for the 2D Saint-Venant system.
//!
//! Conserved variables are `U = (h, hu, hv)` with `h` the water depth,
//! `hu` the discharge per unit area in the `x` direction, and `hv` the
//! discharge per unit area in the `y` direction. Primitives `(h, u, v)`
//! are kept for the API boundary where users reason in physical terms.

use crate::GRAVITY;

/// Conservative state of a single finite-volume cell in 2D.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conserved2D {
    /// Water depth [m]. Non-negative; the solver enforces the wet/dry
    /// threshold at the FV update level, not at the type level.
    pub h: f64,
    /// Discharge per unit area in the `x` direction [m²/s], `h · u`.
    pub hu: f64,
    /// Discharge per unit area in the `y` direction [m²/s], `h · v`.
    pub hv: f64,
}

impl Conserved2D {
    /// New conserved state. No validation: invariants are the solver's job.
    pub const fn new(h: f64, hu: f64, hv: f64) -> Self {
        Self { h, hu, hv }
    }

    /// Dry cell: zero depth, zero discharge in both directions.
    pub const DRY: Self = Self {
        h: 0.0,
        hu: 0.0,
        hv: 0.0,
    };

    /// Convert to primitives. Returns zero velocities for cells with
    /// depth below `dry_tol` to avoid `0/0`.
    pub fn to_primitive(self, dry_tol: f64) -> Primitive2D {
        if self.h > dry_tol {
            Primitive2D {
                h: self.h,
                u: self.hu / self.h,
                v: self.hv / self.h,
            }
        } else {
            Primitive2D {
                h: self.h,
                u: 0.0,
                v: 0.0,
            }
        }
    }

    /// Gravity wave celerity `√(g h)` for this cell.
    pub fn celerity(self) -> f64 {
        (GRAVITY * self.h.max(0.0)).sqrt()
    }

    /// Velocity magnitude `√(u² + v²)` for cells with non-zero depth.
    pub fn speed(self, dry_tol: f64) -> f64 {
        if self.h > dry_tol {
            (self.hu * self.hu + self.hv * self.hv).sqrt() / self.h
        } else {
            0.0
        }
    }
}

/// Primitive state: depth and (u, v) velocity components.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Primitive2D {
    /// Water depth [m].
    pub h: f64,
    /// Depth-averaged velocity in the `x` direction [m/s].
    pub u: f64,
    /// Depth-averaged velocity in the `y` direction [m/s].
    pub v: f64,
}

impl Primitive2D {
    /// New primitive state.
    pub const fn new(h: f64, u: f64, v: f64) -> Self {
        Self { h, u, v }
    }

    /// Convert to conserved variables.
    pub fn to_conserved(self) -> Conserved2D {
        Conserved2D {
            h: self.h,
            hu: self.h * self.u,
            hv: self.h * self.v,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn primitive_conserved_roundtrip() {
        let p = Primitive2D::new(1.5, 2.3, -0.7);
        let c = p.to_conserved();
        let p2 = c.to_primitive(1e-9);
        assert_relative_eq!(p.h, p2.h, epsilon = 1e-12);
        assert_relative_eq!(p.u, p2.u, epsilon = 1e-12);
        assert_relative_eq!(p.v, p2.v, epsilon = 1e-12);
    }

    #[test]
    fn dry_cell_has_zero_velocity_components() {
        let c = Conserved2D::new(1e-10, 1.0, -1.0);
        let p = c.to_primitive(1e-9);
        assert_eq!(p.u, 0.0);
        assert_eq!(p.v, 0.0);
    }

    #[test]
    fn celerity_matches_sqrt_g_h() {
        let c = Conserved2D::new(4.0, 0.0, 0.0);
        assert_relative_eq!(c.celerity(), (GRAVITY * 4.0).sqrt(), epsilon = 1e-12);
    }

    #[test]
    fn speed_is_velocity_magnitude() {
        // u = 3, v = 4 → speed = 5.
        let c = Conserved2D::new(1.0, 3.0, 4.0);
        assert_relative_eq!(c.speed(1e-9), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn dry_cell_has_zero_speed() {
        let c = Conserved2D::new(1e-10, 1.0, 1.0);
        assert_eq!(c.speed(1e-9), 0.0);
    }
}
