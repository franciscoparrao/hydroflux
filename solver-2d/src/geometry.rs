//! 2D structured Cartesian mesh: bed elevations, cell spacings,
//! Manning friction.
//!
//! Indexing convention (see also crate-level doc): `bed[(i, j)]` is
//! the bed elevation at row `i` (y direction) and column `j` (x
//! direction). Row centres are at `y = (i + 0.5) · dy`, column centres
//! at `x = (j + 0.5) · dx`. This matches the GeoTIFF/raster convention
//! where the geotransform maps `(col, row) → (x, y)`.

use ndarray::Array2;

/// Structured Cartesian mesh with uniform cell spacings.
///
/// Cross-sections within each cell are flat (single bed elevation, no
/// sub-grid bathymetry yet) and the channel is rectangular per cell.
/// Variable Manning, sub-grid topography, and AMR are deferred to
/// later iterations.
#[derive(Debug, Clone)]
pub struct Mesh2D {
    /// Bed elevation `z(x, y)` at each cell centre [m]. Shape
    /// `(n_rows, n_cols)` with `i` running along `y` and `j` along `x`.
    pub bed: Array2<f64>,
    /// Cell spacing in the `x` direction [m] (between column centres).
    pub dx: f64,
    /// Cell spacing in the `y` direction [m] (between row centres).
    pub dy: f64,
    /// Manning roughness coefficient (uniform for now).
    pub manning: f64,
}

impl Mesh2D {
    /// Build a mesh from bed elevations and spacings.
    ///
    /// Panics if `dx <= 0`, `dy <= 0`, `manning < 0`, or `bed` is
    /// empty. These are programming errors, not runtime conditions.
    pub fn new(bed: Array2<f64>, dx: f64, dy: f64, manning: f64) -> Self {
        assert!(dx > 0.0, "dx must be strictly positive (got {dx})");
        assert!(dy > 0.0, "dy must be strictly positive (got {dy})");
        assert!(
            manning >= 0.0,
            "Manning n must be non-negative (got {manning})"
        );
        assert!(!bed.is_empty(), "mesh must have at least one cell");
        Self {
            bed,
            dx,
            dy,
            manning,
        }
    }

    /// Number of rows (y-direction cells).
    pub fn n_rows(&self) -> usize {
        self.bed.nrows()
    }

    /// Number of columns (x-direction cells).
    pub fn n_cols(&self) -> usize {
        self.bed.ncols()
    }

    /// Total number of cells.
    pub fn n_cells(&self) -> usize {
        self.bed.len()
    }

    /// Bed slope in the `x` direction at the interior face between
    /// cells `(i, j)` and `(i, j+1)`. Positive when the bed descends
    /// in `+x` (column index increases).
    pub fn bed_slope_x(&self, i: usize, j: usize) -> f64 {
        debug_assert!(j + 1 < self.bed.ncols(), "x-face index out of range");
        (self.bed[(i, j)] - self.bed[(i, j + 1)]) / self.dx
    }

    /// Bed slope in the `y` direction at the interior face between
    /// cells `(i, j)` and `(i+1, j)`. Positive when the bed descends
    /// in `+y` (row index increases). With the GeoTIFF convention
    /// `pixel_height < 0`, `+y` in pixel space points south.
    pub fn bed_slope_y(&self, i: usize, j: usize) -> f64 {
        debug_assert!(i + 1 < self.bed.nrows(), "y-face index out of range");
        (self.bed[(i, j)] - self.bed[(i + 1, j)]) / self.dy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn flat_bed_has_zero_slopes() {
        let bed = array![[10.0, 10.0, 10.0], [10.0, 10.0, 10.0]];
        let mesh = Mesh2D::new(bed, 1.0, 1.0, 0.03);
        for i in 0..mesh.n_rows() {
            for j in 0..mesh.n_cols() - 1 {
                assert_relative_eq!(mesh.bed_slope_x(i, j), 0.0, epsilon = 1e-12);
            }
        }
        for i in 0..mesh.n_rows() - 1 {
            for j in 0..mesh.n_cols() {
                assert_relative_eq!(mesh.bed_slope_y(i, j), 0.0, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn x_slope_of_descending_bed() {
        // Bed descends 1 m every 10 m in x. Slope = 0.1.
        let bed = array![[10.0, 9.0, 8.0], [10.0, 9.0, 8.0]];
        let mesh = Mesh2D::new(bed, 10.0, 10.0, 0.03);
        assert_relative_eq!(mesh.bed_slope_x(0, 0), 0.1, epsilon = 1e-12);
        assert_relative_eq!(mesh.bed_slope_x(1, 1), 0.1, epsilon = 1e-12);
        // y slopes are zero on this bed.
        assert_relative_eq!(mesh.bed_slope_y(0, 0), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn y_slope_of_descending_bed() {
        let bed = array![[10.0, 10.0], [9.0, 9.0], [8.0, 8.0]];
        let mesh = Mesh2D::new(bed, 5.0, 10.0, 0.03);
        assert_relative_eq!(mesh.bed_slope_y(0, 0), 0.1, epsilon = 1e-12);
        assert_relative_eq!(mesh.bed_slope_y(1, 1), 0.1, epsilon = 1e-12);
        assert_relative_eq!(mesh.bed_slope_x(0, 0), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn shape_helpers_match_underlying_array() {
        let bed = Array2::<f64>::zeros((5, 7));
        let mesh = Mesh2D::new(bed, 1.0, 1.0, 0.03);
        assert_eq!(mesh.n_rows(), 5);
        assert_eq!(mesh.n_cols(), 7);
        assert_eq!(mesh.n_cells(), 35);
    }

    #[test]
    #[should_panic(expected = "dx must be strictly positive")]
    fn zero_dx_panics() {
        let _ = Mesh2D::new(array![[1.0]], 0.0, 1.0, 0.03);
    }

    #[test]
    #[should_panic(expected = "dy must be strictly positive")]
    fn zero_dy_panics() {
        let _ = Mesh2D::new(array![[1.0]], 1.0, 0.0, 0.03);
    }

    #[test]
    #[should_panic(expected = "Manning n must be non-negative")]
    fn negative_manning_panics() {
        let _ = Mesh2D::new(array![[1.0]], 1.0, 1.0, -0.01);
    }
}
