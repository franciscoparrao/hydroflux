//! 1D channel geometry: bed elevations, cell spacing, and friction.
//!
//! Cross-sections are rectangular and unit-width for the prototype. Variable
//! Manning, trapezoidal sections, and sub-grid topography are follow-up work
//! (Q4 2026 onwards).

use hydroflux_autograd::Real;
use ndarray::Array1;

/// A uniformly-spaced 1D rectangular channel.
///
/// Generic over `T: Real` in the Manning roughness — the calibration
/// parameter of the differentiable workflow. Instantiate with
/// `T = Dual` and `manning = Dual::variable(n)` to thread `∂/∂n`
/// through the solver; the bed stays `f64` (bed calibration is not in
/// scope for the 1D prototype).
#[derive(Debug, Clone)]
pub struct Channel1D<T = f64> {
    /// Bed elevation `z(x)` at each cell center [m].
    pub bed: Array1<f64>,
    /// Uniform cell spacing `Δx` [m].
    pub dx: f64,
    /// Manning roughness coefficient (uniform for now).
    pub manning: T,
}

impl<T: Real> Channel1D<T> {
    /// Build a channel from bed elevations, spacing, and Manning n.
    ///
    /// Panics if `dx <= 0`, `manning < 0`, or `bed` is empty. These are
    /// programming errors; the solver does not silently swallow them.
    pub fn new(bed: Array1<f64>, dx: f64, manning: T) -> Self {
        assert!(dx > 0.0, "dx must be strictly positive (got {dx})");
        assert!(
            manning.value() >= 0.0,
            "Manning n must be non-negative (got {})",
            manning.value()
        );
        assert!(!bed.is_empty(), "channel must have at least one cell");
        Self { bed, dx, manning }
    }

    /// Number of cells in the channel.
    pub fn n_cells(&self) -> usize {
        self.bed.len()
    }

    /// Bed slope `S₀ = (z_i − z_{i+1}) / Δx` at the face between cell `i`
    /// and cell `i+1`. Positive when the bed descends in the downstream
    /// direction.
    pub fn bed_slope_at_face(&self, i: usize) -> f64 {
        debug_assert!(
            i + 1 < self.bed.len(),
            "face index out of range: {i} of {}",
            self.bed.len()
        );
        (self.bed[i] - self.bed[i + 1]) / self.dx
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn slope_of_uniform_bed_is_zero() {
        let ch = Channel1D::new(array![10.0, 10.0, 10.0], 1.0, 0.03);
        assert_relative_eq!(ch.bed_slope_at_face(0), 0.0, epsilon = 1e-12);
        assert_relative_eq!(ch.bed_slope_at_face(1), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn slope_of_descending_bed_is_positive() {
        // bed descends 1 m per 10 m, so S₀ = 0.1.
        let ch = Channel1D::new(array![10.0, 9.0, 8.0, 7.0], 10.0, 0.03);
        assert_relative_eq!(ch.bed_slope_at_face(0), 0.1, epsilon = 1e-12);
        assert_relative_eq!(ch.bed_slope_at_face(2), 0.1, epsilon = 1e-12);
    }

    #[test]
    #[should_panic(expected = "dx must be strictly positive")]
    fn zero_dx_panics() {
        let _ = Channel1D::new(array![1.0], 0.0, 0.03);
    }

    #[test]
    #[should_panic(expected = "Manning n must be non-negative")]
    fn negative_manning_panics() {
        let _ = Channel1D::new(array![1.0], 1.0, -0.01);
    }

    #[test]
    #[should_panic(expected = "channel must have at least one cell")]
    fn empty_bed_panics() {
        let _ = Channel1D::new(Array1::<f64>::zeros(0), 1.0, 0.03);
    }
}
