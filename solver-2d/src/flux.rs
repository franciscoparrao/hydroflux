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

use crate::GRAVITY;
use crate::state::Conserved2D;

/// Flux through a face whose outward normal is in the `+x` direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluxX {
    /// Mass flux `h · u`.
    pub mass: f64,
    /// `x`-momentum flux `h u² + g h²/2`.
    pub x_momentum: f64,
    /// `y`-momentum flux `h u v`.
    pub y_momentum: f64,
}

/// Flux through a face whose outward normal is in the `+y` direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluxY {
    /// Mass flux `h · v`.
    pub mass: f64,
    /// `x`-momentum flux `h u v`.
    pub x_momentum: f64,
    /// `y`-momentum flux `h v² + g h²/2`.
    pub y_momentum: f64,
}

impl FluxX {
    /// Zero flux. Used as the trivial value for dry-dry interfaces.
    pub const ZERO: Self = Self {
        mass: 0.0,
        x_momentum: 0.0,
        y_momentum: 0.0,
    };

    /// Compute `F_x(U)` for a conservative state. Returns zero flux
    /// when the cell is dry.
    pub fn from_state(u: Conserved2D) -> Self {
        let h = u.h.max(0.0);
        if h == 0.0 {
            return Self::ZERO;
        }
        Self {
            mass: u.hu,
            x_momentum: u.hu * u.hu / h + 0.5 * GRAVITY * h * h,
            y_momentum: u.hu * u.hv / h,
        }
    }
}

impl FluxY {
    /// Zero flux.
    pub const ZERO: Self = Self {
        mass: 0.0,
        x_momentum: 0.0,
        y_momentum: 0.0,
    };

    /// Compute `F_y(U)` for a conservative state. Returns zero flux
    /// when the cell is dry.
    pub fn from_state(u: Conserved2D) -> Self {
        let h = u.h.max(0.0);
        if h == 0.0 {
            return Self::ZERO;
        }
        Self {
            mass: u.hv,
            x_momentum: u.hv * u.hu / h,
            y_momentum: u.hv * u.hv / h + 0.5 * GRAVITY * h * h,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

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
}
