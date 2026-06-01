//! State variables for the 2D Saint-Venant system.
//!
//! Conserved variables are `U = (h, hu, hv)` with `h` the water depth,
//! `hu` the discharge per unit area in the `x` direction, and `hv` the
//! discharge per unit area in the `y` direction. Primitives `(h, u, v)`
//! are kept for the API boundary where users reason in physical terms.
//!
//! The types are generic over a [`Real`] scalar so the identical code
//! evaluates with `T = f64` for production and with `T = Dual` for
//! forward-mode automatic differentiation. The legacy `f64`-only
//! aliases `Conserved2D` and `Primitive2D` are re-exported so existing
//! call sites compile unchanged.

use hydroflux_autograd::Real;

use crate::GRAVITY;

/// Conservative state of a single finite-volume cell in 2D, generic
/// over the scalar type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conserved2DG<T> {
    /// Water depth [m]. Non-negative; the solver enforces the wet/dry
    /// threshold at the FV update level, not at the type level.
    pub h: T,
    /// Discharge per unit area in the `x` direction [m²/s], `h · u`.
    pub hu: T,
    /// Discharge per unit area in the `y` direction [m²/s], `h · v`.
    pub hv: T,
}

/// Conservative state with `f64` storage — the default for production
/// runs. Existing call sites use this name unchanged.
pub type Conserved2D = Conserved2DG<f64>;

impl<T: Real> Conserved2DG<T> {
    /// Generic constructor. No validation: invariants are the solver's
    /// job.
    pub fn new_generic(h: T, hu: T, hv: T) -> Self {
        Self { h, hu, hv }
    }

    /// Dry cell: zero depth, zero discharge in both directions.
    pub fn dry() -> Self {
        Self {
            h: T::zero(),
            hu: T::zero(),
            hv: T::zero(),
        }
    }

    /// Convert to primitives. Returns zero velocities for cells with
    /// depth below `dry_tol` to avoid `0/0`. The branch decision uses
    /// the scalar value of `h` so derivative carry-over is well-defined
    /// on each side of the threshold.
    pub fn to_primitive(self, dry_tol: f64) -> Primitive2DG<T> {
        if self.h.value() > dry_tol {
            Primitive2DG {
                h: self.h,
                u: self.hu / self.h,
                v: self.hv / self.h,
            }
        } else {
            Primitive2DG {
                h: self.h,
                u: T::zero(),
                v: T::zero(),
            }
        }
    }

    /// Gravity wave celerity `√(g h)` for this cell. Clamped to
    /// non-negative depth so a transiently negative draft from the
    /// flux divergence does not produce a `NaN` here.
    pub fn celerity(self) -> T {
        (self.h.max(T::zero()) * GRAVITY).sqrt()
    }

    /// Velocity magnitude `√(u² + v²)` for cells with non-zero depth.
    pub fn speed(self, dry_tol: f64) -> T {
        if self.h.value() > dry_tol {
            (self.hu * self.hu + self.hv * self.hv).sqrt() / self.h
        } else {
            T::zero()
        }
    }
}

// f64-specific niceties preserved for back-compat.
impl Conserved2D {
    /// New conserved state with `f64` storage. Convenience constructor
    /// preserved for back-compat (callers that don't need generics).
    pub const fn new(h: f64, hu: f64, hv: f64) -> Self {
        Self { h, hu, hv }
    }

    /// Dry cell as a compile-time constant — only available for the
    /// `f64` alias because `const fn` cannot call `T::zero()`.
    pub const DRY: Self = Self {
        h: 0.0,
        hu: 0.0,
        hv: 0.0,
    };
}

/// Primitive state: depth and `(u, v)` velocity components, generic
/// over the scalar type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Primitive2DG<T> {
    /// Water depth [m].
    pub h: T,
    /// Depth-averaged velocity in the `x` direction [m/s].
    pub u: T,
    /// Depth-averaged velocity in the `y` direction [m/s].
    pub v: T,
}

/// Primitive state with `f64` storage.
pub type Primitive2D = Primitive2DG<f64>;

impl<T: Real> Primitive2DG<T> {
    /// Generic constructor.
    pub fn new_generic(h: T, u: T, v: T) -> Self {
        Self { h, u, v }
    }

    /// Convert to conserved variables.
    pub fn to_conserved(self) -> Conserved2DG<T> {
        Conserved2DG {
            h: self.h,
            hu: self.h * self.u,
            hv: self.h * self.v,
        }
    }
}

impl Primitive2D {
    /// `f64`-only convenience constructor preserved for back-compat.
    pub const fn new(h: f64, u: f64, v: f64) -> Self {
        Self { h, u, v }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use hydroflux_autograd::Dual;

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

    // ----- Generic-over-Real instantiations: the wedge of the paper. -----

    #[test]
    fn celerity_with_dual_carries_derivative() {
        // d/dh √(g h) = (1/2) · √(g / h).
        let h = Dual::variable(4.0);
        let c = Conserved2DG::<Dual>::new_generic(h, Dual::constant(0.0), Dual::constant(0.0));
        let cel = c.celerity();
        assert_relative_eq!(cel.val, (GRAVITY * 4.0).sqrt(), epsilon = 1e-12);
        // d/dh √(g h) = √g / (2 √h) = √(g / (4 h)).
        let expected = 0.5 * (GRAVITY / 4.0).sqrt();
        assert_relative_eq!(cel.dval, expected, epsilon = 1e-12);
    }

    #[test]
    fn speed_with_dual_carries_derivative_through_quotient() {
        // For h held variable, hu, hv constant: speed = √(hu² + hv²) / h,
        // so d(speed)/dh = − √(hu² + hv²) / h².
        let h = Dual::variable(2.0);
        let hu = Dual::constant(6.0);
        let hv = Dual::constant(8.0); // mom mag = 10
        let c = Conserved2DG::<Dual>::new_generic(h, hu, hv);
        let s = c.speed(1e-9);
        assert_relative_eq!(s.val, 10.0 / 2.0, epsilon = 1e-12);
        assert_relative_eq!(s.dval, -10.0 / (2.0 * 2.0), epsilon = 1e-12);
    }

    #[test]
    fn f64_and_dual_value_match_on_the_same_inputs() {
        // The wedge claim: identical code path, f64 and Dual produce
        // bit-identical .val for any operation that doesn't differentiate
        // (constants, no seeds).
        let f = Conserved2D::new(1.5, 0.6, -0.3);
        let d = Conserved2DG::<Dual>::new_generic(
            Dual::constant(1.5),
            Dual::constant(0.6),
            Dual::constant(-0.3),
        );
        assert_eq!(f.celerity(), d.celerity().val);
        assert_eq!(f.speed(1e-9), d.speed(1e-9).val);
    }

    #[test]
    fn dry_below_tol_branch_uses_value_not_dual() {
        // h.val = 1e-10 < dry_tol = 1e-9 — branch must trigger
        // regardless of dval. Result is T::zero() with .dval = 0.
        let h = Dual { val: 1e-10, dval: 1.0 };
        let c = Conserved2DG::<Dual>::new_generic(h, Dual::constant(1.0), Dual::constant(-1.0));
        let p = c.to_primitive(1e-9);
        assert_eq!(p.u.val, 0.0);
        assert_eq!(p.u.dval, 0.0);
        assert_eq!(p.v.val, 0.0);
        assert_eq!(p.v.dval, 0.0);
    }
}
