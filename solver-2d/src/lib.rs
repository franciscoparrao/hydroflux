//! 2D shallow-water solver on structured Cartesian meshes.
//!
//! Discretisation: finite-volume with the HLLC Riemann solver applied
//! per face via the rotational-invariance trick (a 1D normal-direction
//! Riemann problem at each x- and y-face), Audusse well-balanced
//! bed-slope source via hydrostatic reconstruction in 2D, semi-implicit
//! Manning friction as an operator-split fractional step, and forward
//! Euler in time. The design follows `hydroflux-solver-1d` scaled to
//! two spatial dimensions and three conserved variables `(h, hu, hv)`.
//!
//! # Indexing convention
//!
//! Mesh state is stored as `Array2<Conserved2D>` with shape
//! `(n_rows, n_cols)`. By convention:
//!
//! - row index `i` runs along the `y` direction;
//! - column index `j` runs along the `x` direction;
//! - `mesh.bed[(i, j)]` is the bed elevation at the cell with column
//!   centre at `x = (j + 0.5) · dx` and row centre at
//!   `y = (i + 0.5) · dy`.
//!
//! This matches the GeoTIFF / SurtGIS raster convention where the
//! geotransform maps `(col, row) → (x, y)`. The pixel-height of the
//! geotransform is typically negative (north-up); the solver treats
//! the bed array as already row-ordered top-to-bottom in image space.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod boundary;
pub mod flux;
pub mod geometry;
pub mod riemann;
pub mod state;
pub mod update;

/// Standard gravity used throughout the solver, in m/s².
pub const GRAVITY: f64 = 9.81;

pub use boundary::{Boundaries2D, Boundary, Side, ghost_cell};
pub use flux::{FluxX, FluxY};
pub use geometry::Mesh2D;
pub use riemann::{hllc_flux_x, hllc_flux_y};
pub use state::{Conserved2D, Primitive2D};
pub use update::{cfl_time_step, forward_euler_step, max_wave_speeds};
