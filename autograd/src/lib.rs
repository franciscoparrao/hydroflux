//! Forward-mode automatic differentiation primitives.
//!
//! The crate ships a single number type, [`Dual`], that propagates a
//! scalar value together with its derivative with respect to one
//! seed variable. Any function written generically over operations
//! supported by [`Dual`] computes its derivative automatically when
//! evaluated on a [`Dual`] input.
//!
//! Usage pattern: pick the variable being differentiated, build a
//! [`Dual`] with [`Dual::variable`] for that input and
//! [`Dual::constant`] for all others; evaluate the function; read
//! `.dval` from the result.
//!
//! ```
//! use hydroflux_autograd::Dual;
//!
//! // f(x) = x² + 3x; expected f'(2) = 2·2 + 3 = 7.
//! let x = Dual::variable(2.0);
//! let f = x * x + Dual::constant(3.0) * x;
//! assert_eq!(f.val, 10.0);
//! assert_eq!(f.dval, 7.0);
//! ```
//!
//! Forward-mode is the right fit for the immediate hydroflux use
//! case (low-dimensional parameter sweeps over Manning n, friction,
//! a small handful of BC coefficients). Reverse-mode (tape-based)
//! is a separate, later crate.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod dual;
pub mod physics;
mod real;

pub use dual::Dual;
pub use real::Real;
