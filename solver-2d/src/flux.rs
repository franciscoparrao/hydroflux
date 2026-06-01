//! 2D Saint-Venant flux vectors.
//!
//! The shallow-water system in conservation form has two flux
//! components, one per spatial direction:
//!
//! ```text
//!   F_x(U) = (hu, hu² + g h²/2, huv)
//!   F_y(U) = (hv, huv, hv² + g h²/2)
//! ```
//!
//! The first component is the mass flux, the second is the
//! `x`-momentum flux, the third is the `y`-momentum flux. The
//! hydrostatic-pressure term `g h²/2` appears only in the diagonal
//! components — the `xx` for `F_x` and the `yy` for `F_y`.
//!
//! Types are generic over a [`Real`] scalar; the `f64`-only aliases
//! `FluxX` and `FluxY` are re-exported so existing call sites compile
//! unchanged.

use hydroflux_autograd::Real;

use crate::GRAVITY;
use crate::state::Conserved2DG;

/// Flux through a face whose outward normal is in the `+x` direction.
/// Generic over the scalar type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluxXG<T> {
    /// Mass flux `h · u`.
    pub mass: T,
    /// `x`-momentum flux `h u² + g h²/2`.
    pub x_momentum: T,
    /// `y`-momentum flux `h u v`.
    pub y_momentum: T,
}

/// Flux through a face whose outward normal is in the `+y` direction.
/// Generic over the scalar type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluxYG<T> {
    /// Mass flux `h · v`.
    pub mass: T,
    /// `x`-momentum flux `h u v`.
    pub x_momentum: T,
    /// `y`-momentum flux `h v² + g h²/2`.
    pub y_momentum: T,
}

/// `f64`-storage face flux in the `+x` direction.
pub type FluxX = FluxXG<f64>;
/// `f64`-storage face flux in the `+y` direction.
pub type FluxY = FluxYG<f64>;

impl<T: Real> FluxXG<T> {
    /// Zero flux. Used as the trivial value for dry interfaces.
    pub fn zero() -> Self {
        Self {
            mass: T::zero(),
            x_momentum: T::zero(),
            y_momentum: T::zero(),
        }
    }

    /// Compute `F_x(U)` for a conservative state. Returns zero flux
    /// when the cell is dry. The branch decision uses the scalar value
    /// of `h`; if `h.value() == 0`, the entire flux is identically zero
    /// and no derivative carry-over is meaningful at that operating point.
    pub fn from_state(u: Conserved2DG<T>) -> Self {
        let h = u.h.max(T::zero());
        if h.value() == 0.0 {
            return Self::zero();
        }
        Self {
            mass: u.hu,
            x_momentum: u.hu * u.hu / h + h * h * (0.5 * GRAVITY),
            y_momentum: u.hu * u.hv / h,
        }
    }
}

impl<T: Real> FluxYG<T> {
    /// Zero flux.
    pub fn zero() -> Self {
        Self {
            mass: T::zero(),
            x_momentum: T::zero(),
            y_momentum: T::zero(),
        }
    }

    /// Compute `F_y(U)` for a conservative state. Returns zero flux
    /// when the cell is dry.
    pub fn from_state(u: Conserved2DG<T>) -> Self {
        let h = u.h.max(T::zero());
        if h.value() == 0.0 {
            return Self::zero();
        }
        Self {
            mass: u.hv,
            x_momentum: u.hv * u.hu / h,
            y_momentum: u.hv * u.hv / h + h * h * (0.5 * GRAVITY),
        }
    }
}

// Back-compat constants — only available on the f64 aliases because
// `const fn` cannot call `T::zero()`.
impl FluxX {
    /// Zero flux at compile time. Preserved for existing call sites.
    pub const ZERO: Self = Self {
        mass: 0.0,
        x_momentum: 0.0,
        y_momentum: 0.0,
    };
}

impl FluxY {
    /// Zero flux at compile time. Preserved for existing call sites.
    pub const ZERO: Self = Self {
        mass: 0.0,
        x_momentum: 0.0,
        y_momentum: 0.0,
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::Conserved2D;
    use approx::assert_relative_eq;
    use hydroflux_autograd::Dual;

    #[test]
    fn fx_zero_velocity_gives_hydrostatic_x_momentum() {
        let h = 2.0;
        let f = FluxX::from_state(Conserved2D::new(h, 0.0, 0.0));
        assert_relative_eq!(f.mass, 0.0, epsilon = 1e-12);
        assert_relative_eq!(f.x_momentum, 0.5 * GRAVITY * h * h, epsilon = 1e-12);
        assert_relative_eq!(f.y_momentum, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn fy_zero_velocity_gives_hydrostatic_y_momentum() {
        let h = 2.0;
        let f = FluxY::from_state(Conserved2D::new(h, 0.0, 0.0));
        assert_relative_eq!(f.mass, 0.0, epsilon = 1e-12);
        assert_relative_eq!(f.x_momentum, 0.0, epsilon = 1e-12);
        assert_relative_eq!(f.y_momentum, 0.5 * GRAVITY * h * h, epsilon = 1e-12);
    }

    #[test]
    fn cross_momentum_appears_only_when_both_velocities_nonzero() {
        // u = 1.5, v = 2.0, h = 1.0 → cross term huv = 1.0 * 1.5 * 2.0 = 3.0
        let h = 1.0;
        let u = 1.5;
        let v = 2.0;
        let state = Conserved2D::new(h, h * u, h * v);
        let fx = FluxX::from_state(state);
        let fy = FluxY::from_state(state);
        assert_relative_eq!(fx.y_momentum, h * u * v, epsilon = 1e-12);
        assert_relative_eq!(fy.x_momentum, h * u * v, epsilon = 1e-12);
        // Cross-component is symmetric: F_x.y_mom == F_y.x_mom == huv.
        assert_relative_eq!(fx.y_momentum, fy.x_momentum, epsilon = 1e-12);
    }

    #[test]
    fn fx_matches_textbook_formula() {
        let h = 2.0;
        let u = 1.5;
        let v = -0.7;
        let f = FluxX::from_state(Conserved2D::new(h, h * u, h * v));
        assert_relative_eq!(f.mass, h * u, epsilon = 1e-12);
        assert_relative_eq!(
            f.x_momentum,
            h * u * u + 0.5 * GRAVITY * h * h,
            epsilon = 1e-12
        );
        assert_relative_eq!(f.y_momentum, h * u * v, epsilon = 1e-12);
    }

    #[test]
    fn dry_cell_yields_zero_flux_in_both_directions() {
        assert_eq!(FluxX::from_state(Conserved2D::DRY), FluxX::ZERO);
        assert_eq!(FluxY::from_state(Conserved2D::DRY), FluxY::ZERO);
    }

    // ----- Generic-over-Real instantiations: AD-ready flux. -----

    #[test]
    fn fx_with_dual_seed_on_h_gives_correct_derivative() {
        // F_x.x_momentum(h, hu, hv) = hu²/h + g h²/2 with hu, hv fixed.
        // ∂/∂h = − hu²/h² + g h. At h=2, hu=3, hv=0:
        //   value     = 9/2 + g · 2 = 4.5 + 2g
        //   derivative = −9/4 + 2g
        let h = Dual::variable(2.0);
        let hu = Dual::constant(3.0);
        let hv = Dual::constant(0.0);
        let f = FluxXG::<Dual>::from_state(Conserved2DG::new_generic(h, hu, hv));
        assert_relative_eq!(f.x_momentum.val, 4.5 + 2.0 * GRAVITY, epsilon = 1e-12);
        assert_relative_eq!(f.x_momentum.dval, -9.0 / 4.0 + 2.0 * GRAVITY, epsilon = 1e-12);
        // Mass is hu (held constant in this test), so dval = 0.
        assert_eq!(f.mass.val, 3.0);
        assert_eq!(f.mass.dval, 0.0);
    }

    #[test]
    fn fy_with_dual_seed_on_v_gives_correct_y_momentum_derivative() {
        // F_y.y_momentum = hv²/h + g h²/2. Held h, hu constant; seed hv.
        // ∂/∂hv at fixed h = 2·hv/h. At h=2, hv=4: deriv = 4.
        let h = Dual::constant(2.0);
        let hu = Dual::constant(0.0);
        let hv = Dual::variable(4.0);
        let f = FluxYG::<Dual>::from_state(Conserved2DG::new_generic(h, hu, hv));
        // value: 16/2 + g · 2 = 8 + 2g
        assert_relative_eq!(f.y_momentum.val, 8.0 + 2.0 * GRAVITY, epsilon = 1e-12);
        assert_relative_eq!(f.y_momentum.dval, 4.0, epsilon = 1e-12);
    }

    #[test]
    fn fx_dual_value_matches_f64_value_on_same_inputs() {
        // Bit-identical .val between f64 and Dual when no seed is active.
        let state_f = Conserved2D::new(1.5, 0.6, -0.3);
        let state_d = Conserved2DG::<Dual>::new_generic(
            Dual::constant(1.5),
            Dual::constant(0.6),
            Dual::constant(-0.3),
        );
        let f = FluxX::from_state(state_f);
        let d = FluxXG::<Dual>::from_state(state_d);
        assert_eq!(f.mass, d.mass.val);
        assert_eq!(f.x_momentum, d.x_momentum.val);
        assert_eq!(f.y_momentum, d.y_momentum.val);
    }

    #[test]
    fn dry_branch_under_dual_returns_zero_flux_with_zero_derivative() {
        // h = 0 exactly — the dry branch must fire and produce
        // FluxXG::zero() with dval = 0 on every component.
        let h = Dual { val: 0.0, dval: 1.0 };
        let f = FluxXG::<Dual>::from_state(Conserved2DG::new_generic(
            h,
            Dual::constant(0.5),
            Dual::constant(0.5),
        ));
        assert_eq!(f.mass.val, 0.0);
        assert_eq!(f.mass.dval, 0.0);
        assert_eq!(f.x_momentum.val, 0.0);
        assert_eq!(f.x_momentum.dval, 0.0);
        assert_eq!(f.y_momentum.val, 0.0);
        assert_eq!(f.y_momentum.dval, 0.0);
    }
}
