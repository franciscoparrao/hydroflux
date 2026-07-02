//! Saint-Venant 1D solver for open-channel flow.
//!
//! Implements the conservative shallow water equations in 1D using a
//! finite-volume discretization with an HLL Riemann solver (Toro 2009,
//! §10.5.1; wave-speed estimate after Davis 1988; two-rarefaction
//! dry-front speeds after Toro §10.5.4).
//!
//! The numerical core is generic over [`hydroflux_autograd::Real`]
//! (defaulting to `f64`): instantiated with `Dual` and a seeded Manning
//! roughness, the same production code path yields `∂h/∂n` by
//! forward-mode AD — see `tests/ad_gradient.rs`. Wet/dry and CFL
//! branching decide on `Real::value()` (the primal), the discipline
//! established in `hydroflux-solver-2d`.
//!
//! See `outline.md` § "Plan Año 1 detallado — Fase 2 (2026 Q3)" for the
//! design rationale and the open decisions this prototype closes.

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

/// Wet/dry threshold in metres. Cells with `h ≤ H_DRY` are treated as
/// dry for momentum purposes: their velocity is undefined (`hu/h` would
/// blow up on a vanishing film), so [`update::forward_euler_step`] zeroes
/// `hu` while keeping the residual mass, which preserves conservation.
/// Matches the convention of `hydroflux-solver-2d`.
pub const H_DRY: f64 = 1.0e-6;

/// Velocity cutoff in metres, one order of magnitude above [`H_DRY`].
/// The post-update floor zeroes `hu` for `h ≤ H_VEL`, and the CFL
/// bound ignores the velocity of such films: `hu/h` with residual
/// momentum just above `H_DRY` produces arbitrarily large wave speeds
/// that collapse `dt`. Matches `hydroflux-solver-2d::H_VEL`.
pub const H_VEL: f64 = 1.0e-5;

pub use boundary::{Boundaries, Boundary, Side, ghost_cell};
pub use flux::Flux;
pub use geometry::Channel1D;
pub use io::{IoError, read_channel, write_bed, write_depth, write_discharge};
pub use riemann::hll_flux;
pub use source::manning_friction_step;
pub use state::{Conserved, Primitive};
pub use update::{cfl_time_step, forward_euler_step, max_wave_speed};
