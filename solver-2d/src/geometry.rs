//! 2D structured Cartesian mesh: bed elevations, cell spacings,
//! Manning friction.
//!
//! Indexing convention (see also crate-level doc): `bed[(i, j)]` is
//! the bed elevation at row `i` (y direction) and column `j` (x
//! direction). Row centres are at `y = (i + 0.5) · dy`, column centres
//! at `x = (j + 0.5) · dx`. This matches the GeoTIFF/raster convention
//! where the geotransform maps `(col, row) → (x, y)`.
//!
//! Generic over a [`Real`] scalar so that bed elevations and Manning
//! coefficients can carry a derivative (`Dual`) for inverse-problem
//! work, while cell spacings stay `f64` — `dx`/`dy` are spatial
//! metadata, not parameters to differentiate over. The legacy `f64`
//! alias `Mesh2D` keeps existing call sites compiling unchanged.

use hydroflux_autograd::Real;
use ndarray::Array2;

/// Structured Cartesian mesh with uniform cell spacings, generic over
/// the scalar type used for bed elevations and Manning roughness.
///
/// Cross-sections within each cell are flat (single bed elevation, no
/// sub-grid bathymetry yet) and the channel is rectangular per cell.
/// Sub-grid topography and AMR are deferred to later iterations.
///
/// Manning friction is stored as a per-cell field; pass a scalar to
/// [`Mesh2DG::new`] for uniform roughness, or build a spatially
/// varying field separately and use [`Mesh2DG::with_manning_field`]
/// for landcover-derived friction maps.
#[derive(Debug, Clone)]
pub struct Mesh2DG<T> {
    /// Bed elevation `z(x, y)` at each cell centre [m]. Shape
    /// `(n_rows, n_cols)` with `i` running along `y` and `j` along `x`.
    pub bed: Array2<T>,
    /// Cell spacing in the `x` direction [m] (between column centres).
    pub dx: f64,
    /// Cell spacing in the `y` direction [m] (between row centres).
    pub dy: f64,
    /// Manning roughness coefficient `n` per cell [s/m^(1/3)]. Same
    /// shape as `bed`. For uniform roughness all entries equal the
    /// scalar passed to [`Mesh2DG::new`]; for landcover-derived maps
    /// use [`Mesh2DG::with_manning_field`].
    pub manning: Array2<T>,
}

/// `f64`-storage mesh — the default for production runs. Existing
/// call sites use this name unchanged.
pub type Mesh2D = Mesh2DG<f64>;

impl<T: Real> Mesh2DG<T> {
    /// Build a mesh from bed elevations and spacings with **uniform**
    /// Manning roughness. The Manning field is filled with `manning`
    /// in every cell — equivalent to the original scalar API.
    ///
    /// Panics if `dx <= 0`, `dy <= 0`, `manning < 0` (compared on the
    /// scalar value), or `bed` is empty. These are programming errors,
    /// not runtime conditions.
    pub fn new(bed: Array2<T>, dx: f64, dy: f64, manning: T) -> Self {
        assert!(dx > 0.0, "dx must be strictly positive (got {dx})");
        assert!(dy > 0.0, "dy must be strictly positive (got {dy})");
        assert!(
            manning.value() >= 0.0,
            "Manning n must be non-negative (got {})",
            manning.value()
        );
        assert!(!bed.is_empty(), "mesh must have at least one cell");
        let manning_field = Array2::from_elem(bed.dim(), manning);
        Self {
            bed,
            dx,
            dy,
            manning: manning_field,
        }
    }

    /// Build a mesh with a **spatially varying** Manning roughness
    /// field. `manning` must have the same shape as `bed`; every cell
    /// value must be non-negative.
    ///
    /// Use this when the friction coefficient is derived from a
    /// landcover raster (e.g. ESA WorldCover via
    /// [`crate::io::mesh_from_geotiff_with_landcover`]) or when the
    /// channel and overbank have distinct calibrated `n` values.
    ///
    /// Panics if `dx <= 0`, `dy <= 0`, `bed` is empty, the shapes of
    /// `bed` and `manning` differ, or any Manning value is negative.
    pub fn with_manning_field(
        bed: Array2<T>,
        dx: f64,
        dy: f64,
        manning: Array2<T>,
    ) -> Self {
        assert!(dx > 0.0, "dx must be strictly positive (got {dx})");
        assert!(dy > 0.0, "dy must be strictly positive (got {dy})");
        assert!(!bed.is_empty(), "mesh must have at least one cell");
        assert_eq!(
            bed.dim(),
            manning.dim(),
            "manning field shape {:?} must match bed shape {:?}",
            manning.dim(),
            bed.dim(),
        );
        assert!(
            manning.iter().all(|&n| n.value() >= 0.0),
            "Manning n must be non-negative at every cell"
        );
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
    pub fn bed_slope_x(&self, i: usize, j: usize) -> T {
        debug_assert!(j + 1 < self.bed.ncols(), "x-face index out of range");
        (self.bed[(i, j)] - self.bed[(i, j + 1)]) / self.dx
    }

    /// Bed slope in the `y` direction at the interior face between
    /// cells `(i, j)` and `(i+1, j)`. Positive when the bed descends
    /// in `+y` (row index increases). With the GeoTIFF convention
    /// `pixel_height < 0`, `+y` in pixel space points south.
    pub fn bed_slope_y(&self, i: usize, j: usize) -> T {
        debug_assert!(i + 1 < self.bed.nrows(), "y-face index out of range");
        (self.bed[(i, j)] - self.bed[(i + 1, j)]) / self.dy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use hydroflux_autograd::Dual;
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

    // ----- Generic-over-Real instantiations. -----

    #[test]
    fn bed_slope_with_dual_carries_derivative_through_subtraction() {
        // bed_slope_x = (z_L − z_R) / dx; with z_L as the seed variable,
        // ∂slope/∂z_L = 1 / dx.
        let z_l = Dual::variable(10.0);
        let z_r = Dual::constant(9.0);
        let bed = Array2::from_shape_vec((1, 2), vec![z_l, z_r]).unwrap();
        let mesh = Mesh2DG::<Dual>::new(bed, 10.0, 10.0, Dual::constant(0.03));
        let slope = mesh.bed_slope_x(0, 0);
        assert_relative_eq!(slope.val, 0.1, epsilon = 1e-12);
        assert_relative_eq!(slope.dval, 1.0 / 10.0, epsilon = 1e-12);
    }

    #[test]
    fn mesh_with_dual_manning_seeds_the_friction_field() {
        // The whole point of generic Mesh: a seeded Dual manning value
        // propagates through downstream Manning friction. Here we only
        // check that the field stores the seed correctly.
        let n_seed = Dual::variable(0.04);
        let bed = Array2::<Dual>::from_elem((2, 3), Dual::constant(0.0));
        let mesh = Mesh2DG::<Dual>::new(bed, 1.0, 1.0, n_seed);
        // Every cell carries the seeded derivative (dval = 1 from variable()).
        for &n in mesh.manning.iter() {
            assert_eq!(n.val, 0.04);
            assert_eq!(n.dval, 1.0);
        }
    }
}
