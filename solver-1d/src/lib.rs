//! Saint-Venant 1D solver for open-channel flow.
//!
//! Implements the conservative shallow water equations in 1D using a
//! finite-volume discretization with an HLL Riemann solver (Toro 2009,
//! §10.5.1; wave-speed estimate after Davis 1988).
//!
//! See `outline.md` § "Plan Año 1 detallado — Fase 2 (2026 Q3)" for the
//! design rationale and the open decisions this prototype closes.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod flux;
pub mod geometry;
pub mod riemann;
pub mod state;

/// Standard gravity used throughout the solver, in m/s².
pub const GRAVITY: f64 = 9.81;

pub use flux::Flux;
pub use geometry::Channel1D;
pub use riemann::hll_flux;
pub use state::{Conserved, Primitive};
