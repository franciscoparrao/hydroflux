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
pub mod io;
pub mod riemann;
pub mod source;
pub mod state;
pub mod update;

/// Standard gravity used throughout the solver, in m/s².
pub const GRAVITY: f64 = 9.81;

/// Wet/dry threshold [m]. Cells with `h ≤ H_DRY` are treated as dry
/// by every dynamic consumer (CFL accounting, velocity computation,
/// flux evaluation, the dry-dry face short-circuit): their velocity is
/// undefined and their faces carry no flux. Their **mass is kept**,
/// however — a film below the threshold sits inert on the bed until
/// inflow accumulates past `H_DRY`. Destroying it (the previous
/// behaviour) leaked up to `H_DRY` of depth per cell-event along the
/// entire wetting-front perimeter at every step, a systematic volume
/// bias in long inundation runs.
///
/// The choice 10⁻⁶ m is a balance: small enough that a 1 mm puddle is
/// still "wet" (a millimetre of water carries meaningful momentum on
/// metre-scale meshes), but large enough that `u = hu / h` does not
/// blow up under finite-precision arithmetic for cells that are
/// effectively at the wet/dry front. The Riemann solver uses a
/// tighter internal threshold (10⁻¹²) to pick the two-rarefaction
/// wave speed; that is a separate concern from the cell-level
/// definition of "dry".
pub const H_DRY: f64 = 1.0e-6;

/// Velocity cutoff [m], one order of magnitude above [`H_DRY`]. Cells
/// with `h ≤ H_VEL` have their momentum zeroed by the post-update
/// floor and contribute no velocity to the CFL bound: on a film this
/// thin the depth-averaged velocity is not meaningful, and `u = hu/h`
/// with residual momentum just above `H_DRY` produces arbitrarily
/// large wave speeds that collapse `dt` (a stall, not an instability).
/// The dual-threshold pattern follows BASEMENT/SERGHEI practice.
pub const H_VEL: f64 = 1.0e-5;

pub use boundary::{Boundaries2D, Boundary, Side, ghost_cell};
pub use flux::{FluxX, FluxY};
pub use flux::{FluxXG, FluxYG};
pub use geometry::{Mesh2D, Mesh2DG};
pub use io::{
    depth_raster_from_states, esa_worldcover_to_manning, mesh_from_geotiff,
    mesh_from_geotiff_with_landcover, write_depth_geotiff,
};
pub use riemann::{hllc_flux_x, hllc_flux_y};
pub use source::{PointSource, apply_point_sources, apply_rain, manning_friction_step};
pub use state::{Conserved2D, Conserved2DG, Primitive2D, Primitive2DG};
pub use update::{
    StepWorkspace2D, cfl_time_step, cfl_time_step_with_bcs, forward_euler_step,
    forward_euler_step_with, max_wave_speeds, ssprk2_step, ssprk2_step_with,
};
