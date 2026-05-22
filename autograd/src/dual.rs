//! Dual number for forward-mode automatic differentiation.

use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A dual number `val + dval·ε` with `ε² = 0`.
///
/// `val` is the function value, `dval` is the partial derivative
/// with respect to the seed variable. Arithmetic on `Dual` propagates
/// derivatives via the standard rules (sum, product, quotient, chain).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual {
    /// Function value at the current point.
    pub val: f64,
    /// Derivative of the value with respect to the seed variable.
    pub dval: f64,
}

impl Dual {
    /// A dual number with explicit value and derivative parts.
    pub const fn new(val: f64, dval: f64) -> Self {
        Self { val, dval }
    }

    /// The seed variable: `val = x`, `dval = 1`.
    ///
    /// Use this for the input being differentiated.
    pub const fn variable(x: f64) -> Self {
        Self { val: x, dval: 1.0 }
    }

    /// A constant: `val = x`, `dval = 0`.
    ///
    /// Use this for parameters that are NOT being differentiated.
    pub const fn constant(x: f64) -> Self {
        Self { val: x, dval: 0.0 }
    }

    /// Square root. Undefined for negative `val`; returns NaN.
    pub fn sqrt(self) -> Self {
        let s = self.val.sqrt();
        Self {
            val: s,
            dval: self.dval / (2.0 * s),
        }
    }

    /// Natural exponential.
    pub fn exp(self) -> Self {
        let e = self.val.exp();
        Self {
            val: e,
            dval: self.dval * e,
        }
    }

    /// Natural logarithm. Undefined for non-positive `val`.
    pub fn ln(self) -> Self {
        Self {
            val: self.val.ln(),
            dval: self.dval / self.val,
        }
    }

    /// Sine.
    pub fn sin(self) -> Self {
        Self {
            val: self.val.sin(),
            dval: self.dval * self.val.cos(),
        }
    }

    /// Cosine.
    pub fn cos(self) -> Self {
        Self {
            val: self.val.cos(),
            dval: -self.dval * self.val.sin(),
        }
    }

    /// Absolute value. Derivative at `val = 0` is taken as 0 (the
    /// subdifferential 0 lies in [-1, 1]; picking 0 avoids spurious
    /// gradient signals when the optimiser sits on the kink).
    pub fn abs(self) -> Self {
        let sign = if self.val > 0.0 {
            1.0
        } else if self.val < 0.0 {
            -1.0
        } else {
            0.0
        };
        Self {
            val: self.val.abs(),
            dval: self.dval * sign,
        }
    }

    /// Integer power. `n` is treated as a non-differentiable constant.
    pub fn powi(self, n: i32) -> Self {
        Self {
            val: self.val.powi(n),
            dval: (n as f64) * self.val.powi(n - 1) * self.dval,
        }
    }

    /// Real power with a constant exponent. The exponent is NOT
    /// differentiated; use [`Dual::powd`] for the case where it is.
    pub fn powf(self, n: f64) -> Self {
        Self {
            val: self.val.powf(n),
            dval: n * self.val.powf(n - 1.0) * self.dval,
        }
    }

    /// Real power where both base and exponent carry derivatives.
    /// `d/dx [a(x)^b(x)] = a^b · (b'·ln a + b·a'/a)`.
    pub fn powd(self, other: Dual) -> Self {
        let val = self.val.powf(other.val);
        let dval = val * (other.dval * self.val.ln() + other.val * self.dval / self.val);
        Self { val, dval }
    }

    /// Maximum of two duals. The derivative jumps at the crossing
    /// point; on the kink we return the average to keep the gradient
    /// bounded and symmetric.
    pub fn max(self, other: Dual) -> Self {
        if self.val > other.val {
            self
        } else if other.val > self.val {
            other
        } else {
            Self {
                val: self.val,
                dval: 0.5 * (self.dval + other.dval),
            }
        }
    }

    /// Minimum of two duals. Symmetric counterpart to [`Dual::max`].
    pub fn min(self, other: Dual) -> Self {
        if self.val < other.val {
            self
        } else if other.val < self.val {
            other
        } else {
            Self {
                val: self.val,
                dval: 0.5 * (self.dval + other.dval),
            }
        }
    }
}

impl fmt::Display for Dual {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} + {}ε", self.val, self.dval)
    }
}

// -- Dual op Dual -----------------------------------------------------

impl Add for Dual {
    type Output = Dual;
    fn add(self, rhs: Dual) -> Dual {
        Dual {
            val: self.val + rhs.val,
            dval: self.dval + rhs.dval,
        }
    }
}

impl Sub for Dual {
    type Output = Dual;
    fn sub(self, rhs: Dual) -> Dual {
        Dual {
            val: self.val - rhs.val,
            dval: self.dval - rhs.dval,
        }
    }
}

impl Mul for Dual {
    type Output = Dual;
    fn mul(self, rhs: Dual) -> Dual {
        Dual {
            val: self.val * rhs.val,
            dval: self.dval * rhs.val + self.val * rhs.dval,
        }
    }
}

impl Div for Dual {
    type Output = Dual;
    fn div(self, rhs: Dual) -> Dual {
        Dual {
            val: self.val / rhs.val,
            dval: (self.dval * rhs.val - self.val * rhs.dval) / (rhs.val * rhs.val),
        }
    }
}

impl Neg for Dual {
    type Output = Dual;
    fn neg(self) -> Dual {
        Dual {
            val: -self.val,
            dval: -self.dval,
        }
    }
}

// -- Dual op f64 / f64 op Dual ---------------------------------------
//
// `f64` is treated as a constant (`dval = 0`). Providing both sides
// makes call sites read naturally (`x * 2.0` and `2.0 * x` both work).

impl Add<f64> for Dual {
    type Output = Dual;
    fn add(self, rhs: f64) -> Dual {
        Dual {
            val: self.val + rhs,
            dval: self.dval,
        }
    }
}

impl Add<Dual> for f64 {
    type Output = Dual;
    fn add(self, rhs: Dual) -> Dual {
        Dual {
            val: self + rhs.val,
            dval: rhs.dval,
        }
    }
}

impl Sub<f64> for Dual {
    type Output = Dual;
    fn sub(self, rhs: f64) -> Dual {
        Dual {
            val: self.val - rhs,
            dval: self.dval,
        }
    }
}

impl Sub<Dual> for f64 {
    type Output = Dual;
    fn sub(self, rhs: Dual) -> Dual {
        Dual {
            val: self - rhs.val,
            dval: -rhs.dval,
        }
    }
}

impl Mul<f64> for Dual {
    type Output = Dual;
    fn mul(self, rhs: f64) -> Dual {
        Dual {
            val: self.val * rhs,
            dval: self.dval * rhs,
        }
    }
}

impl Mul<Dual> for f64 {
    type Output = Dual;
    fn mul(self, rhs: Dual) -> Dual {
        Dual {
            val: self * rhs.val,
            dval: self * rhs.dval,
        }
    }
}

impl Div<f64> for Dual {
    type Output = Dual;
    fn div(self, rhs: f64) -> Dual {
        Dual {
            val: self.val / rhs,
            dval: self.dval / rhs,
        }
    }
}

impl Div<Dual> for f64 {
    type Output = Dual;
    fn div(self, rhs: Dual) -> Dual {
        Dual {
            val: self / rhs.val,
            dval: -self * rhs.dval / (rhs.val * rhs.val),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    const EPS: f64 = 1.0e-12;

    #[test]
    fn variable_seed_carries_unit_derivative() {
        let x = Dual::variable(3.0);
        assert_eq!(x.val, 3.0);
        assert_eq!(x.dval, 1.0);
    }

    #[test]
    fn constant_has_zero_derivative() {
        let c = Dual::constant(7.5);
        assert_eq!(c.val, 7.5);
        assert_eq!(c.dval, 0.0);
    }

    #[test]
    fn add_sub_propagate_derivatives() {
        let x = Dual::variable(2.0);
        let f = x + x - Dual::constant(1.0);
        // f(x) = 2x - 1; f'(x) = 2.
        assert_eq!(f.val, 3.0);
        assert_eq!(f.dval, 2.0);
    }

    #[test]
    fn power_rule_via_mul() {
        // f(x) = x² ⇒ f'(2) = 4.
        let x = Dual::variable(2.0);
        let f = x * x;
        assert_eq!(f.val, 4.0);
        assert_eq!(f.dval, 4.0);
    }

    #[test]
    fn quotient_rule() {
        // f(x) = x / (x + 1); f'(x) = 1/(x+1)².
        let x = Dual::variable(3.0);
        let f = x / (x + 1.0);
        assert_relative_eq!(f.val, 0.75, epsilon = EPS);
        assert_relative_eq!(f.dval, 1.0 / 16.0, epsilon = EPS);
    }

    #[test]
    fn neg_flips_value_and_derivative() {
        let x = Dual::variable(5.0);
        let f = -x;
        assert_eq!(f.val, -5.0);
        assert_eq!(f.dval, -1.0);
    }

    #[test]
    fn scalar_left_and_right_multiplication_are_consistent() {
        let x = Dual::variable(4.0);
        let a = 2.0 * x;
        let b = x * 2.0;
        assert_eq!(a, b);
        assert_eq!(a.dval, 2.0);
    }

    #[test]
    fn scalar_division_constant_over_variable() {
        // f(x) = 1 / x ⇒ f'(2) = -1/4.
        let x = Dual::variable(2.0);
        let f = 1.0 / x;
        assert_relative_eq!(f.val, 0.5, epsilon = EPS);
        assert_relative_eq!(f.dval, -0.25, epsilon = EPS);
    }

    #[test]
    fn sqrt_derivative() {
        // f(x) = sqrt(x) ⇒ f'(4) = 1/4.
        let x = Dual::variable(4.0);
        let f = x.sqrt();
        assert_relative_eq!(f.val, 2.0, epsilon = EPS);
        assert_relative_eq!(f.dval, 0.25, epsilon = EPS);
    }

    #[test]
    fn exp_is_its_own_derivative() {
        let x = Dual::variable(1.5);
        let f = x.exp();
        assert_relative_eq!(f.val, 1.5_f64.exp(), epsilon = EPS);
        assert_relative_eq!(f.dval, 1.5_f64.exp(), epsilon = EPS);
    }

    #[test]
    fn ln_derivative() {
        // f(x) = ln(x) ⇒ f'(e) = 1/e.
        let x = Dual::variable(std::f64::consts::E);
        let f = x.ln();
        assert_relative_eq!(f.val, 1.0, epsilon = EPS);
        assert_relative_eq!(f.dval, 1.0 / std::f64::consts::E, epsilon = EPS);
    }

    #[test]
    fn sin_cos_derivatives_at_pi_over_4() {
        let x = Dual::variable(PI / 4.0);
        let s = x.sin();
        let c = x.cos();
        let sqrt2_2 = 2.0_f64.sqrt() / 2.0;
        // sin(π/4) = √2/2; d/dx sin = cos ⇒ cos(π/4) = √2/2.
        assert_relative_eq!(s.val, sqrt2_2, epsilon = EPS);
        assert_relative_eq!(s.dval, sqrt2_2, epsilon = EPS);
        // cos(π/4) = √2/2; d/dx cos = -sin ⇒ -√2/2.
        assert_relative_eq!(c.val, sqrt2_2, epsilon = EPS);
        assert_relative_eq!(c.dval, -sqrt2_2, epsilon = EPS);
    }

    #[test]
    fn chain_rule_sin_of_x_squared() {
        // f(x) = sin(x²) ⇒ f'(x) = 2x · cos(x²); at x = 1: 2·cos(1).
        let x = Dual::variable(1.0);
        let f = (x * x).sin();
        assert_relative_eq!(f.val, 1.0_f64.sin(), epsilon = EPS);
        assert_relative_eq!(f.dval, 2.0 * 1.0_f64.cos(), epsilon = EPS);
    }

    #[test]
    fn powi_matches_repeated_multiplication() {
        let x = Dual::variable(2.0);
        let by_powi = x.powi(4);
        let by_mul = x * x * x * x;
        assert_relative_eq!(by_powi.val, by_mul.val, epsilon = EPS);
        assert_relative_eq!(by_powi.dval, by_mul.dval, epsilon = EPS);
    }

    #[test]
    fn powf_derivative_with_constant_exponent() {
        // f(x) = x^(3/5); f'(x) = (3/5) x^(-2/5). Manning normal-depth
        // shape (h_n = (n·q/√S₀)^(3/5)): exponent 3/5 over the bracket.
        let x = Dual::variable(2.0);
        let f = x.powf(0.6);
        assert_relative_eq!(f.val, 2.0_f64.powf(0.6), epsilon = EPS);
        assert_relative_eq!(f.dval, 0.6 * 2.0_f64.powf(-0.4), epsilon = EPS);
    }

    #[test]
    fn powd_with_both_sides_dual() {
        // f(x) = x^x at x = 2: val = 4; f'(x) = x^x · (ln x + 1)
        // ⇒ f'(2) = 4·(ln 2 + 1).
        let x = Dual::variable(2.0);
        let f = x.powd(x);
        assert_relative_eq!(f.val, 4.0, epsilon = EPS);
        assert_relative_eq!(f.dval, 4.0 * (2.0_f64.ln() + 1.0), epsilon = EPS);
    }

    #[test]
    fn abs_pos_neg_and_zero() {
        let pos = Dual::variable(2.0).abs();
        assert_eq!(pos.val, 2.0);
        assert_eq!(pos.dval, 1.0);
        let neg = Dual::variable(-3.0).abs();
        assert_eq!(neg.val, 3.0);
        assert_eq!(neg.dval, -1.0);
        let zero = Dual::variable(0.0).abs();
        assert_eq!(zero.val, 0.0);
        assert_eq!(zero.dval, 0.0);
    }

    #[test]
    fn max_min_pick_branch_and_carry_its_derivative() {
        let x = Dual::variable(2.0);
        let y = Dual::constant(5.0);
        let mx = x.max(y);
        assert_eq!(mx.val, 5.0);
        assert_eq!(mx.dval, 0.0); // y was the constant branch
        let mn = x.min(y);
        assert_eq!(mn.val, 2.0);
        assert_eq!(mn.dval, 1.0); // x was the variable branch
    }

    #[test]
    fn manning_friction_term_derivative_matches_finite_diff() {
        // S_f(h, q) = n² q |q| / h^(10/3). Differentiate w.r.t. h
        // at h = 0.5, q = 1.2, n = 0.035.
        // Analytical d/dh S_f = -(10/3) n² q |q| / h^(13/3).
        let q = 1.2_f64;
        let n = 0.035_f64;
        let h = Dual::variable(0.5);
        let sf = Dual::constant(n * n * q * q.abs()) / h.powf(10.0 / 3.0);
        let analytic = -(10.0 / 3.0) * n * n * q * q.abs() / 0.5_f64.powf(13.0 / 3.0);
        assert_relative_eq!(sf.dval, analytic, epsilon = 1.0e-10);
    }

    #[test]
    fn critical_depth_derivative_w_r_t_discharge() {
        // h_c(q) = (q² / g)^(1/3); dh_c/dq = (2/3)·(q/g)^(2/3)/q^(1/3)
        // simpler: dh_c/dq = (2/3)·q^(-1/3)·g^(-1/3) = (2/3)·(g·q)^(-1/3)·(q^0)
        // Direct: d/dq [(q²/g)^(1/3)] = (2/3) q^(-1/3) g^(-1/3).
        let g = 9.81;
        let q = Dual::variable(2.0);
        let hc = (q * q / g).powf(1.0 / 3.0);
        let analytic = (2.0 / 3.0) * 2.0_f64.powf(-1.0 / 3.0) * g.powf(-1.0 / 3.0);
        assert_relative_eq!(hc.dval, analytic, epsilon = 1.0e-12);
    }
}
