//! Vector-mode dual number: one value, `N` derivative components.
//!
//! The scalar [`Dual`](crate::Dual) carries a single derivative, so a
//! gradient over `P` parameters costs `P` independent forward passes.
//! Each of those passes recomputes the same values — the same square
//! roots, the same powers, the same mesh traversal — and differs only in
//! which partial derivative rides along. `DualN<N>` carries `N`
//! derivatives through one pass, so that shared work is paid once.
//!
//! The saving is large because of where the cost sits. Propagating a
//! tangent is multiply-add; computing the value it rides on involves
//! `sqrt`, `powf` and division. Measured on a kernel with the solver's
//! arithmetic mix (`examples/vector_mode_spike.rs`), each additional
//! component costs about 0.06x a scalar pass rather than the 1.04x a
//! separate pass costs, which moves the break-even against reverse-mode
//! AD from roughly two parameters to a few dozen.
//!
//! The non-smooth points follow the same conventions as `Dual` and for
//! the same reasons: `sqrt` at zero returns the subdifferential element
//! 0 so a clamp-to-dry composition yields a finite gradient at the
//! shoreline rather than NaN, `abs` at zero returns 0, and `max`/`min`
//! at a tie average the incoming derivatives so an optimiser sitting on
//! the kink sees a bounded, symmetric value.

use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A dual number `val + Σ dval[i]·εᵢ` with `εᵢεⱼ = 0`.
///
/// Component `i` is the partial derivative with respect to the `i`-th
/// seed variable. Seeding is by [`DualN::seeded`]; everything else is a
/// constant and carries a zero row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DualN<const N: usize> {
    /// Function value at the current point.
    pub val: f64,
    /// Partial derivatives with respect to each of the `N` seeds.
    pub dval: [f64; N],
}

impl<const N: usize> DualN<N> {
    /// Explicit value and derivative row.
    pub const fn new(val: f64, dval: [f64; N]) -> Self {
        Self { val, dval }
    }

    /// A constant: all derivative components zero.
    pub const fn constant(val: f64) -> Self {
        Self { val, dval: [0.0; N] }
    }

    /// The `k`-th seed variable: `dval[k] = 1`, the rest zero.
    ///
    /// Panics if `k >= N`, which is a programming error rather than a
    /// runtime condition — the seed index is fixed when the gradient is
    /// set up.
    pub fn seeded(val: f64, k: usize) -> Self {
        assert!(k < N, "seed index {k} out of range for DualN<{N}>");
        let mut dval = [0.0; N];
        dval[k] = 1.0;
        Self { val, dval }
    }

    /// Apply the chain rule with a scalar local derivative: given
    /// `f(x)` and `f'(x)`, produce the dual for `f` at this point.
    #[inline]
    fn chain(self, val: f64, dfdx: f64) -> Self {
        let mut dval = [0.0; N];
        for i in 0..N {
            dval[i] = self.dval[i] * dfdx;
        }
        Self { val, dval }
    }

    /// Combine two duals componentwise under a bilinear rule.
    #[inline]
    fn combine(self, other: Self, val: f64, a: f64, b: f64) -> Self {
        let mut dval = [0.0; N];
        for i in 0..N {
            dval[i] = a * self.dval[i] + b * other.dval[i];
        }
        Self { val, dval }
    }

    /// Square root. At `val = 0` the analytical derivative is
    /// unbounded; the subdifferential element 0 is returned instead so
    /// that `h.max(0).sqrt()` at a dry cell gives a finite, zero
    /// gradient rather than NaN.
    pub fn sqrt(self) -> Self {
        let s = self.val.sqrt();
        self.chain(s, if s == 0.0 { 0.0 } else { 0.5 / s })
    }

    /// Natural exponential.
    pub fn exp(self) -> Self {
        let e = self.val.exp();
        self.chain(e, e)
    }

    /// Natural logarithm.
    pub fn ln(self) -> Self {
        self.chain(self.val.ln(), 1.0 / self.val)
    }

    /// Sine.
    pub fn sin(self) -> Self {
        self.chain(self.val.sin(), self.val.cos())
    }

    /// Cosine.
    pub fn cos(self) -> Self {
        self.chain(self.val.cos(), -self.val.sin())
    }

    /// Absolute value. At the origin the derivative is taken as 0, the
    /// symmetric choice, so an optimiser on the kink sees no gradient
    /// rather than an arbitrary sign.
    pub fn abs(self) -> Self {
        let sign = if self.val > 0.0 {
            1.0
        } else if self.val < 0.0 {
            -1.0
        } else {
            0.0
        };
        self.chain(self.val.abs(), sign)
    }

    /// Integer power. `n` is a non-differentiable constant.
    pub fn powi(self, n: i32) -> Self {
        self.chain(self.val.powi(n), (n as f64) * self.val.powi(n - 1))
    }

    /// Real power with a constant exponent.
    pub fn powf(self, n: f64) -> Self {
        self.chain(self.val.powf(n), n * self.val.powf(n - 1.0))
    }

    /// Real power where both base and exponent carry derivatives:
    /// `d/dx [a(x)^b(x)] = a^b · (b'·ln a + b·a'/a)`.
    pub fn powd(self, other: Self) -> Self {
        let val = self.val.powf(other.val);
        let ln_a = self.val.ln();
        let inv_a = 1.0 / self.val;
        let mut dval = [0.0; N];
        for i in 0..N {
            dval[i] = val * (other.dval[i] * ln_a + other.val * self.dval[i] * inv_a);
        }
        Self { val, dval }
    }

    /// Maximum. On a tie the incoming derivatives are averaged, which
    /// keeps the result bounded and symmetric on the kink.
    pub fn max(self, other: Self) -> Self {
        if self.val > other.val {
            self
        } else if other.val > self.val {
            other
        } else {
            self.combine(other, self.val, 0.5, 0.5)
        }
    }

    /// Minimum. Symmetric counterpart to [`DualN::max`].
    pub fn min(self, other: Self) -> Self {
        if self.val < other.val {
            self
        } else if other.val < self.val {
            other
        } else {
            self.combine(other, self.val, 0.5, 0.5)
        }
    }
}

impl<const N: usize> fmt::Display for DualN<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} + {:?}ε", self.val, self.dval)
    }
}

impl<const N: usize> Add for DualN<N> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        self.combine(o, self.val + o.val, 1.0, 1.0)
    }
}

impl<const N: usize> Sub for DualN<N> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self.combine(o, self.val - o.val, 1.0, -1.0)
    }
}

impl<const N: usize> Mul for DualN<N> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        self.combine(o, self.val * o.val, o.val, self.val)
    }
}

impl<const N: usize> Div for DualN<N> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let inv = 1.0 / o.val;
        self.combine(o, self.val * inv, inv, -self.val * inv * inv)
    }
}

impl<const N: usize> Neg for DualN<N> {
    type Output = Self;
    fn neg(self) -> Self {
        self.chain(-self.val, -1.0)
    }
}

macro_rules! scalar_ops {
    ($($tr:ident, $m:ident, $op:tt;)*) => {$(
        impl<const N: usize> $tr<f64> for DualN<N> {
            type Output = Self;
            fn $m(self, rhs: f64) -> Self {
                self $op DualN::<N>::constant(rhs)
            }
        }
        impl<const N: usize> $tr<DualN<N>> for f64 {
            type Output = DualN<N>;
            fn $m(self, rhs: DualN<N>) -> DualN<N> {
                DualN::<N>::constant(self) $op rhs
            }
        }
    )*};
}
scalar_ops! { Add, add, +; Sub, sub, -; Mul, mul, *; Div, div, /; }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Dual;

    /// The whole point of the type is that it agrees with N separate
    /// scalar passes. This exercises the full arithmetic surface the
    /// solver uses and checks every component against the scalar Dual
    /// seeded on the same variable.
    #[test]
    fn components_match_independent_scalar_passes() {
        // f(a, b, c) mixes the operations the solver's inner loop uses.
        fn f_scalar(a: Dual, b: Dual, c: Dual) -> Dual {
            let wave = (a * 9.81).sqrt();
            let fric = (b * b) / a.powf(4.0 / 3.0);
            (wave * fric + c.abs()).max(Dual::constant(1e-9)) - c / (a + 1.0)
        }
        fn f_vec(a: DualN<3>, b: DualN<3>, c: DualN<3>) -> DualN<3> {
            let wave = (a * 9.81).sqrt();
            let fric = (b * b) / a.powf(4.0 / 3.0);
            (wave * fric + c.abs()).max(DualN::<3>::constant(1e-9)) - c / (a + 1.0)
        }

        let (av, bv, cv) = (1.7_f64, 0.035_f64, -0.42_f64);
        let vec = f_vec(
            DualN::<3>::seeded(av, 0),
            DualN::<3>::seeded(bv, 1),
            DualN::<3>::seeded(cv, 2),
        );

        for k in 0..3 {
            let mk = |i: usize, v: f64| {
                if i == k { Dual::variable(v) } else { Dual::constant(v) }
            };
            let sc = f_scalar(mk(0, av), mk(1, bv), mk(2, cv));
            assert!(
                (vec.val - sc.val).abs() < 1e-14,
                "value mismatch: {} vs {}",
                vec.val,
                sc.val
            );
            assert!(
                (vec.dval[k] - sc.dval).abs() < 1e-12,
                "component {k}: vector {} vs scalar {}",
                vec.dval[k],
                sc.dval
            );
        }
    }

    #[test]
    fn constants_carry_no_derivative() {
        let x = DualN::<4>::constant(3.0);
        assert_eq!(x.dval, [0.0; 4]);
    }

    #[test]
    fn seeding_sets_exactly_one_component() {
        let x = DualN::<4>::seeded(2.0, 2);
        assert_eq!(x.dval, [0.0, 0.0, 1.0, 0.0]);
    }

    /// The wet/dry conventions are the reason the solver can take a
    /// gradient through a moving shoreline at all, so they are pinned
    /// here as well as in the scalar type.
    #[test]
    fn sqrt_at_zero_is_finite_and_zero() {
        let r = DualN::<2>::seeded(0.0, 0).sqrt();
        assert_eq!(r.val, 0.0);
        assert!(r.dval.iter().all(|d| d.is_finite() && *d == 0.0));
    }

    #[test]
    fn max_on_a_tie_averages_components() {
        let a = DualN::<2>::seeded(1.0, 0);
        let b = DualN::<2>::seeded(1.0, 1);
        assert_eq!(a.max(b).dval, [0.5, 0.5]);
    }

    #[test]
    fn division_matches_quotient_rule_per_component() {
        let a = DualN::<2>::seeded(3.0, 0);
        let b = DualN::<2>::seeded(0.5, 1);
        let q = a / b;
        assert!((q.val - 6.0).abs() < 1e-14);
        assert!((q.dval[0] - 1.0 / 0.5).abs() < 1e-14);
        assert!((q.dval[1] - (-3.0 / 0.25)).abs() < 1e-12);
    }
}
