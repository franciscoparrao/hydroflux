//! Well-balanced property: lake at rest over piecewise-discontinuous
//! bumpy bathymetry must stay at rest to machine precision.
//!
//! # Test description
//!
//! Greenberg & Le Roux 1996, Audusse et al. 2004: a finite-volume
//! scheme is *well-balanced* if it preserves the steady state
//! `η = h + z = const` and `q = 0` exactly, on arbitrary topography
//! including non-smooth beds. This is the canonical sanity check for
//! the bed-slope source treatment: discretisations that compute the
//! source term independently of the flux divergence introduce a small
//! but persistent imbalance whenever `dz/dx` is large or discontinuous.
//!
//! The bed here is built from random piecewise-constant blocks plus
//! a few sharp ridges, with all elevations below the still water
//! level so every cell remains wet. The Liang & Marche 2009
//! reconstruction (`z_face = midpoint(z_L, z_R)`) combined with the
//! cell-centred algebraic source `S = (g/2)(h_R² − h_L²)/dx` makes
//! the discrete divergence of flux exactly cancel the source for any
//! still state — so this test should pass to round-off.
//!
//! Setup:
//! - Domain: 200 m × 100 m, mesh 100 × 50 (`dx = dy = 2 m`).
//! - Bed: deterministic LCG-style permutation of three height
//!   levels (0.0, 0.3, 0.6 m) on a 10 × 5 macro-block grid, plus
//!   two narrow raised ridges in `x` and `y`.
//! - Initial: `η = 1.0 m` everywhere → `h = η − z ≥ 0.4 m`. `q = 0`.
//! - BC: walls on all sides.
//! - `t_end = 60 s`, CFL 0.4.
//! - Tolerance: `‖η − 1.0‖∞ < 1e-10` and `‖q‖∞ < 1e-10`.
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test lake_at_rest_bumpy
//! ```

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2D, Mesh2D, cfl_time_step, manning_friction_step, ssprk2_step,
};
use ndarray::Array2;

const N_X: usize = 100;
const N_Y: usize = 50;
const DX: f64 = 2.0;
const DY: f64 = 2.0;
const ETA_STILL: f64 = 1.0;
const T_END: f64 = 60.0;
const CFL: f64 = 0.4;
const MANNING: f64 = 0.0; // No friction so the test is purely about the bed-slope source.

/// Deterministic bumpy bed: macro-block heights from a tiny LCG +
/// two raised ridges. Heights are bounded so every cell stays wet
/// at η = 1.0 m.
fn build_bumpy_bed() -> Array2<f64> {
    let levels = [0.0_f64, 0.3, 0.6];
    let macro_x = 10; // 10 columns of macro-blocks
    let macro_y = 5; // 5 rows of macro-blocks
    let block_w = N_X / macro_x;
    let block_h = N_Y / macro_y;

    // Deterministic pseudo-random index for each macro block.
    let mut seed: u32 = 1664525;
    let mut block_choice = vec![0usize; macro_x * macro_y];
    for slot in block_choice.iter_mut() {
        seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
        *slot = (seed as usize) % levels.len();
    }

    let mut bed = Array2::from_elem((N_Y, N_X), 0.0);
    for i in 0..N_Y {
        let bi = i / block_h;
        for j in 0..N_X {
            let bj = j / block_w;
            let idx = bi * macro_x + bj;
            bed[(i, j)] = levels[block_choice[idx]];
        }
    }

    // Sharp ridge running along x at j = 30..32, height 0.5 m.
    for i in 0..N_Y {
        for j in 30..32 {
            bed[(i, j)] = 0.5;
        }
    }
    // Sharp ridge running along y at i = 20..22, height 0.45 m.
    for i in 20..22 {
        for j in 0..N_X {
            bed[(i, j)] = bed[(i, j)].max(0.45);
        }
    }
    bed
}

fn initial_lake(bed: &Array2<f64>) -> Array2<Conserved2D> {
    Array2::from_shape_fn((N_Y, N_X), |(i, j)| {
        let h = (ETA_STILL - bed[(i, j)]).max(0.0);
        Conserved2D::new(h, 0.0, 0.0)
    })
}

fn run_to_t_end(mesh: &Mesh2D, mut states: Array2<Conserved2D>) -> Array2<Conserved2D> {
    let bcs = Boundaries2D::WALLS;
    let mut t = 0.0;
    let mut steps = 0;
    while t < T_END {
        let dt = cfl_time_step(&states, mesh, CFL).min(T_END - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        manning_friction_step(&mut states, MANNING, dt, 1.0e-9);
        t += dt;
        steps += 1;
        if steps > 500_000 {
            panic!("lake_at_rest_bumpy: {steps} steps without reaching t_end");
        }
    }
    states
}

#[test]
fn surface_elevation_stays_flat_to_round_off() {
    let bed = build_bumpy_bed();
    let mesh = Mesh2D::new(bed.clone(), DX, DY, MANNING);
    let initial = initial_lake(&bed);

    // Sanity: all cells start wet at the still water level.
    for ((i, j), s) in initial.indexed_iter() {
        let eta = s.h + bed[(i, j)];
        assert!((eta - ETA_STILL).abs() < 1.0e-15);
        assert_eq!(s.hu, 0.0);
        assert_eq!(s.hv, 0.0);
    }

    let final_states = run_to_t_end(&mesh, initial);

    // η = h + z must stay at ETA_STILL to round-off.
    let mut max_eta_err = 0.0_f64;
    for ((i, j), s) in final_states.indexed_iter() {
        let eta = s.h + bed[(i, j)];
        let err = (eta - ETA_STILL).abs();
        if err > max_eta_err {
            max_eta_err = err;
        }
    }
    assert!(
        max_eta_err < 1.0e-10,
        "‖η − {ETA_STILL}‖∞ = {max_eta_err:.3e} (should be < 1e-10 for well-balanced scheme)"
    );
}

#[test]
fn momentum_stays_zero_to_round_off() {
    let bed = build_bumpy_bed();
    let mesh = Mesh2D::new(bed.clone(), DX, DY, MANNING);
    let initial = initial_lake(&bed);
    let final_states = run_to_t_end(&mesh, initial);

    let mut max_q = 0.0_f64;
    for s in &final_states {
        let q = s.hu.abs().max(s.hv.abs());
        if q > max_q {
            max_q = q;
        }
    }
    assert!(
        max_q < 1.0e-10,
        "‖q‖∞ = {max_q:.3e} (should be < 1e-10 — well-balanced scheme produces no spurious flow over bumpy bed)"
    );
}

#[test]
fn mass_is_conserved_to_machine_precision() {
    let bed = build_bumpy_bed();
    let mesh = Mesh2D::new(bed.clone(), DX, DY, MANNING);
    let initial = initial_lake(&bed);
    let m_initial: f64 = initial.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();

    let final_states = run_to_t_end(&mesh, initial);
    let m_final: f64 = final_states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();

    let rel_err = (m_final - m_initial).abs() / m_initial;
    assert!(
        rel_err < 1.0e-12,
        "mass not conserved over walls: m_initial = {m_initial:.6e}, m_final = {m_final:.6e}, rel_err = {rel_err:.3e}"
    );
}
