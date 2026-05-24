//! Abstract real-valued scalar.
//!
//! The [`Real`] trait captures the surface area of `f64` that the
//! hydroflux solvers actually use: arithmetic, scaling by `f64`,
//! square root, absolute value, branch-by-branch selection (max/min),
//! integer and real powers. Both [`f64`] and [`crate::Dual`] implement
//! it, so functions written generically over `T: Real` evaluate as
//! normal arithmetic when `T = f64` and propagate forward-mode
//! derivatives when `T = Dual`.
//!
//! Comparisons go through [`Real::value`], not `PartialOrd`. The point
//! is to make explicit that branching decisions are taken on the
//! scalar value alone; the derivative carried by `Dual` does not — and
//! must not — influence control flow.

use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::Dual;

/// Scalar type supporting the arithmetic the SWE primitives need.
pub trait Real:
    Copy
    + Default
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + Add<f64, Output = Self>
    + Sub<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    /// Additive identity.
    fn zero() -> Self;

    /// Multiplicative identity.
    fn one() -> Self;

    /// Lift a `f64` constant. For `Dual` this becomes [`Dual::constant`].
    fn from_f64(x: f64) -> Self;

    /// The scalar value. For `f64` returns `self`; for `Dual` returns
    /// `self.val`. Use this for branch decisions and for emitting
    /// results back to non-generic code.
    fn value(self) -> f64;

    /// Square root.
    fn sqrt(self) -> Self;

    /// Absolute value.
    fn abs(self) -> Self;

    /// Maximum of two values.
    fn max(self, other: Self) -> Self;

    /// Minimum of two values.
    fn min(self, other: Self) -> Self;

    /// Integer power.
    fn powi(self, n: i32) -> Self;

    /// Real power with a non-differentiable exponent.
    fn powf(self, n: f64) -> Self;

    /// Real power with a *differentiable* exponent (`Self^Self`).
    /// Enables calibrating both base and exponent of a power-law
    /// relation (e.g., cross-section top width `T(h) = c · h^p`
    /// where both `c` and `p` are inferred parameters).
    fn powt(self, exponent: Self) -> Self;
}

impl Real for f64 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn from_f64(x: f64) -> Self {
        x
    }
    fn value(self) -> f64 {
        self
    }
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    fn abs(self) -> Self {
        f64::abs(self)
    }
    fn max(self, other: Self) -> Self {
        f64::max(self, other)
    }
    fn min(self, other: Self) -> Self {
        f64::min(self, other)
    }
    fn powi(self, n: i32) -> Self {
        f64::powi(self, n)
    }
    fn powf(self, n: f64) -> Self {
        f64::powf(self, n)
    }
    fn powt(self, exponent: f64) -> f64 {
        f64::powf(self, exponent)
    }
}

impl Default for Dual {
    fn default() -> Self {
        Dual::constant(0.0)
    }
}

impl Real for Dual {
    fn zero() -> Self {
        Dual::constant(0.0)
    }
    fn one() -> Self {
        Dual::constant(1.0)
    }
    fn from_f64(x: f64) -> Self {
        Dual::constant(x)
    }
    fn value(self) -> f64 {
        self.val
    }
    fn sqrt(self) -> Self {
        Dual::sqrt(self)
    }
    fn abs(self) -> Self {
        Dual::abs(self)
    }
    fn max(self, other: Self) -> Self {
        Dual::max(self, other)
    }
    fn min(self, other: Self) -> Self {
        Dual::min(self, other)
    }
    fn powi(self, n: i32) -> Self {
        Dual::powi(self, n)
    }
    fn powf(self, n: f64) -> Self {
        Dual::powf(self, n)
    }
    fn powt(self, exponent: Dual) -> Dual {
        Dual::powd(self, exponent)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const EPS: f64 = 1.0e-12;

    /// Same generic function consumed by both backends.
    fn quadratic<T: Real>(x: T) -> T {
        x * x + x * 3.0 + T::from_f64(2.0)
    }

    #[test]
    fn quadratic_runs_with_f64_and_dual_with_matching_value() {
        let y_f = quadratic::<f64>(4.0);
        let y_d = quadratic::<Dual>(Dual::variable(4.0));
        assert_relative_eq!(y_f, y_d.val, epsilon = EPS);
        // f'(x) = 2x + 3, f'(4) = 11.
        assert_relative_eq!(y_d.dval, 11.0, epsilon = EPS);
    }

    /// Sqrt over Real picks the right branch for both backends.
    fn norm2<T: Real>(x: T, y: T) -> T {
        (x * x + y * y).sqrt()
    }

    #[test]
    fn norm2_is_consistent_across_backends() {
        let r_f = norm2::<f64>(3.0, 4.0);
        let r_d = norm2::<Dual>(Dual::variable(3.0), Dual::constant(4.0));
        assert_relative_eq!(r_f, 5.0, epsilon = EPS);
        assert_relative_eq!(r_d.val, 5.0, epsilon = EPS);
        // d/dx √(x² + y²) at x=3, y=4 → x / √(x²+y²) = 3/5.
        assert_relative_eq!(r_d.dval, 0.6, epsilon = EPS);
    }

    #[test]
    fn value_extracts_scalar_for_branch_decisions() {
        let a = 1.5_f64;
        let b = Dual::variable(1.5);
        assert_eq!(a.value(), 1.5);
        assert_eq!(b.value(), 1.5);
        // Branch on .value(), not on the dual itself.
        assert!(a.value() > 1.0);
        assert!(b.value() > 1.0);
    }

    #[test]
    fn powf_with_fractional_exponent_carries_chain_rule() {
        // f(x) = x^(5/3). Manning normal-depth shape (after inversion).
        // f'(2) = (5/3) · 2^(2/3).
        let x = Dual::variable(2.0);
        let f = x.powf(5.0 / 3.0);
        assert_relative_eq!(f.val, 2.0_f64.powf(5.0 / 3.0), epsilon = EPS);
        assert_relative_eq!(f.dval, (5.0 / 3.0) * 2.0_f64.powf(2.0 / 3.0), epsilon = EPS);
    }
}
