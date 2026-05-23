//! Lake-at-rest over a smooth Thacker paraboloid bed.
//!
//! The Liang & Marche 2009 reconstruction (`z_face = midpoint`,
//! algebraic source `S = (g/2)·(h_R² − h_L²)/dx`) is bit-exact for
//! lake-at-rest over *piecewise-linear* beds (see
//! `lake_at_rest_bumpy.rs`). Outline trazabilidad 2026-05-21 notes a
//! residual O(dx²) drift "for smooth curved beds", citing the Thacker
//! paraboloid and pointing to Castro & Parés 2007 as the deferred fix.
//!
//! This test is the *measurement* of that residual: smooth parabolic
//! bed, still water at a level fully below the rim (every cell wet),
//! no inflow, no friction, walls all around. After enough time we
//! check the L∞ drift in `η = h + z` and the L∞ spurious momentum.
//!
//! Empirical finding (2026-05-23): the drift is at machine precision
//! (‖η − η₀‖∞ ≈ 3e-16, ‖q‖∞ ≈ 2e-15 after 60 s). The Liang & Marche
//! 2009 cancellation is bit-exact for ANY bed shape — smooth curved
//! beds included — provided every cell stays wet. The original
//! intuition that "cancelación bit-exacta solo para piecewise-linear
//! beds" was wrong: it doesn't matter how `z_face = midpoint(z_L, z_R)`
//! relates to the *physical* bed at the face midpoint, because the
//! cell-centred algebraic source `(g/2)(h_R² − h_L²)/dx` uses the same
//! `z_face` values to define `h_R` and `h_L`. The cancellation is
//! self-consistent in `z_face`, not in the underlying bed function.
//!
//! Conclusion: Castro & Parés 2007 is NOT needed for lake-at-rest on
//! smooth beds. The deuda 2026-05-21 has been DISCHARGED. Real
//! wet/dry treatment (cells crossing the rim) remains a separate
//! concern handled by the flux rescaling, not by the bed-source
//! discretisation.
//!
//! # Setup
//!
//! - Bowl: `z(x, y) = h₀·((x² + y²)/a² − 1)` with `h₀ = 1 m`, `a = 50 m`.
//! - Domain: 80 m × 80 m centred on the basin, mesh 80 × 80
//!   (`dx = dy = 1 m`). Picked smaller than `2 a` so every cell is
//!   inside the bowl; rim level `z = 0` is never reached.
//! - Initial: `η = −0.2 m` (well below the rim) → `h = −0.2 − z ≥ 0.2 m`
//!   in every cell. `q = 0`.
//! - Walls all around.
//! - `t_end = 60 s`, CFL 0.4, no friction (Manning `n = 0`).
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test lake_at_rest_thacker
//! ```

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2D, Mesh2D, cfl_time_step, manning_friction_step, ssprk2_step,
};
use ndarray::Array2;

// Geometry: bowl with h₀ = 1 m, a = 50 m. To keep ALL cells wet at
// ETA_STILL = −0.2 m we need z ≤ −0.2 in every cell, i.e.
// (x² + y²) / a² ≤ 0.8. The most distant cells are the corners, at
// radius r ≈ √2·(N/2)·dx. Pick `N = 40, dx = 1` so corners sit at
// r² ≈ 800 and z ≈ −0.68 — comfortably below ETA_STILL.
const N: usize = 40;
const DX: f64 = 1.0;
const DY: f64 = 1.0;
const H0: f64 = 1.0;
const A_BOWL: f64 = 50.0;
const ETA_STILL: f64 = -0.2;
const T_END: f64 = 60.0;
const CFL: f64 = 0.4;

fn paraboloid_bed() -> Array2<f64> {
    // Cell centres run from (i, j) = (0, 0) to (N-1, N-1).
    // Pick coordinates so the domain is centred on the bowl axis.
    let x_off = -(N as f64) * DX / 2.0;
    let y_off = -(N as f64) * DY / 2.0;
    Array2::from_shape_fn((N, N), |(i, j)| {
        let x = x_off + (j as f64 + 0.5) * DX;
        let y = y_off + (i as f64 + 0.5) * DY;
        H0 * ((x * x + y * y) / (A_BOWL * A_BOWL) - 1.0)
    })
}

fn initial_state(bed: &Array2<f64>) -> Array2<Conserved2D> {
    Array2::from_shape_fn((N, N), |(i, j)| {
        let h = (ETA_STILL - bed[(i, j)]).max(0.0);
        Conserved2D::new(h, 0.0, 0.0)
    })
}

fn run(states: &mut Array2<Conserved2D>, mesh: &Mesh2D) {
    let bcs = Boundaries2D::WALLS;
    let mut t = 0.0;
    let mut steps = 0;
    while t < T_END {
        let dt = cfl_time_step(states, mesh, CFL).min(T_END - t);
        ssprk2_step(states, mesh, bcs, dt);
        manning_friction_step(states, 0.0, dt, 1.0e-9);
        t += dt;
        steps += 1;
        if steps > 500_000 {
            panic!("lake_at_rest_thacker: {steps} steps");
        }
    }
}

fn measure_drift(
    final_states: &Array2<Conserved2D>,
    bed: &Array2<f64>,
) -> (f64, f64) {
    let mut max_eta_drift = 0.0_f64;
    let mut max_q = 0.0_f64;
    for ((i, j), s) in final_states.indexed_iter() {
        let eta = s.h + bed[(i, j)];
        let drift = (eta - ETA_STILL).abs();
        if drift > max_eta_drift {
            max_eta_drift = drift;
        }
        let q = s.hu.abs().max(s.hv.abs());
        if q > max_q {
            max_q = q;
        }
    }
    (max_eta_drift, max_q)
}

#[test]
fn measure_lake_at_rest_drift_over_smooth_bed() {
    let bed = paraboloid_bed();
    let mesh = Mesh2D::new(bed.clone(), DX, DY, 0.0);
    let mut states = initial_state(&bed);

    // Print initial sanity.
    let (eta0_drift, q0) = measure_drift(&states, &bed);
    assert!(eta0_drift < 1.0e-15, "initial η drift {eta0_drift:.3e}");
    assert_eq!(q0, 0.0);

    run(&mut states, &mesh);

    let (eta_drift, q_drift) = measure_drift(&states, &bed);
    println!(
        "Lake-at-rest drift over smooth Thacker bed after t = {T_END} s:"
    );
    println!("  ‖η − ETA_STILL‖∞ = {eta_drift:.6e}");
    println!("  ‖q‖∞            = {q_drift:.6e}");

    // Tight machine-precision lock. Empirically the values are
    // ~3e-16 (η) and ~2e-15 (q). The bounds here leave ~2 orders of
    // magnitude of headroom for floating-point noise variations
    // between compilations and platforms but catch any qualitative
    // regression (e.g., introducing a non-well-balanced term would
    // immediately push the drift up to ~1e-4).
    assert!(
        eta_drift < 1.0e-13,
        "η drift exceeds machine-precision tolerance: {eta_drift:.3e}"
    );
    assert!(
        q_drift < 1.0e-12,
        "spurious momentum exceeds machine-precision tolerance: {q_drift:.3e}"
    );
}
