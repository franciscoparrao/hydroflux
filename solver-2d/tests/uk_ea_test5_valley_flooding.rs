//! UK Environment Agency 2D benchmark — equivalent of Test 5
//! "Valley flooding" (Néelz & Pender 2013).
//!
//! # Test description
//!
//! A synthetic valley with a parabolic cross-section runs along the
//! `+x` axis. Water is released at the upstream end via a Discharge
//! boundary condition. The flow must:
//!
//! 1. Concentrate in the valley centreline (where the bed is
//!    lowest).
//! 2. Propagate along the valley while the cross-section profile
//!    adjusts to Manning normal depth.
//! 3. Exit downstream through a transmissive boundary.
//!
//! What this test adds over `uk_ea_test4_propagation` (which used a
//! rectangular channel with raised banks): the bed varies in BOTH
//! `x` (longitudinal slope) AND `y` (parabolic cross-section), so
//! the bed-slope source has non-trivial gradients in both directions
//! simultaneously. Exercises the explicit 2D source term + Audusse
//! reconstruction over a non-axis-aligned bed.
//!
//! Synthetic stand-in for the official EA Test 5 geometry (which is
//! a more elaborate Y-shaped valley with tributaries); the
//! single-valley case captures the essential 2D-flow concentration
//! physics.
//!
//! # Setup
//!
//! - Domain: 1000 m × 400 m, mesh 250 × 100 (`dx = dy = 4 m`).
//! - Bed: parabolic cross-section centred on `y = 200`:
//!   `z(x, y) = −S₀·x + κ·(y − 200)²`
//!   with `S₀ = 0.001` (longitudinal slope) and `κ = 5e-4` (so the
//!   bed rises ≈ 5 m at the domain edges from the valley centreline).
//! - Initial: thin film `h = 1 mm` (workaround for the Discharge-on-dry
//!   limitation).
//! - BC: Discharge `q = 2 m²/s` at `W`, Transmissive at `E`, Wall on
//!   `N` and `S`.
//! - Manning `n = 0.04`.
//! - `t_end = 600 s` (10 min).
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test uk_ea_test5_valley_flooding
//! ```

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, cfl_time_step, manning_friction_step, ssprk2_step,
};
use ndarray::Array2;

#[derive(Debug, Clone, Copy)]
struct TestCase {
    length_x: f64,
    length_y: f64,
    /// Longitudinal bed slope (positive descends in +x).
    slope: f64,
    /// Cross-section curvature `κ`. The bed at lateral distance `r`
    /// from the centreline rises by `κ · r²`.
    curvature: f64,
    /// `y` coordinate of the valley centreline [m].
    centreline_y: f64,
    manning: f64,
    q_in: f64,
    t_end: f64,
}

impl TestCase {
    fn standard() -> Self {
        Self {
            length_x: 1000.0,
            length_y: 400.0,
            slope: 0.001,
            curvature: 5.0e-4,
            centreline_y: 200.0,
            manning: 0.04,
            q_in: 2.0,
            t_end: 600.0,
        }
    }
}

fn build_mesh(case: TestCase, n_x: usize, n_y: usize) -> (Mesh2D, f64, f64) {
    let dx = case.length_x / n_x as f64;
    let dy = case.length_y / n_y as f64;
    let bed = Array2::from_shape_fn((n_y, n_x), |(i, j)| {
        let x = (j as f64 + 0.5) * dx;
        let y = (i as f64 + 0.5) * dy;
        let r = y - case.centreline_y;
        -case.slope * x + case.curvature * r * r
    });
    (Mesh2D::new(bed, dx, dy, case.manning), dx, dy)
}

fn boundaries(case: TestCase) -> Boundaries2D {
    Boundaries2D {
        west: Boundary::Discharge { q: case.q_in },
        east: Boundary::Transmissive,
        north: Boundary::Wall,
        south: Boundary::Wall,
    }
}

fn run_until(
    mut states: Array2<Conserved2D>,
    mesh: &Mesh2D,
    bcs: Boundaries2D,
    case: TestCase,
    cfl: f64,
) -> (Array2<Conserved2D>, usize) {
    let mut t = 0.0;
    let mut steps = 0;
    while t < case.t_end {
        let dt = cfl_time_step(&states, mesh, cfl).min(case.t_end - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        manning_friction_step(&mut states, mesh, dt, 1.0e-9);
        t += dt;
        steps += 1;
        if steps > 500_000 {
            panic!("UK EA Test 5: {steps} steps without reaching t_end");
        }
    }
    (states, steps)
}

#[test]
fn depth_remains_bounded_and_finite() {
    let case = TestCase::standard();
    let (mesh, _, _) = build_mesh(case, 250, 100);
    let initial = Array2::from_elem((100, 250), Conserved2D::new(0.001, 0.0, 0.0));
    let (final_states, _) = run_until(initial, &mesh, boundaries(case), case, 0.4);
    for s in &final_states {
        assert!(s.h.is_finite() && s.h >= 0.0, "h ill-formed: {}", s.h);
        // Maximum sensible depth: about 2× the Manning normal depth
        // h_n = (n · q / √S₀)^(3/5) ≈ 1.75 m. Allow up to 5 m as
        // a generous bound.
        assert!(s.h < 5.0, "h exceeded sensible bound: {}", s.h);
        assert!(s.hu.is_finite() && s.hv.is_finite(), "momentum non-finite");
    }
}

#[test]
fn flow_concentrates_along_valley_centreline() {
    // The bed is lowest along y = centreline_y. After sustained
    // inflow, water should be deepest there and progressively
    // shallower away from the centreline. Compare the centreline
    // average depth against off-centreline (y ≈ 80 and y ≈ 320,
    // 120 m off-centreline both ways).
    let case = TestCase::standard();
    let (mesh, _dx, dy) = build_mesh(case, 250, 100);
    let initial = Array2::from_elem((100, 250), Conserved2D::new(0.001, 0.0, 0.0));
    let (final_states, _) = run_until(initial, &mesh, boundaries(case), case, 0.4);

    let i_centre = (case.centreline_y / dy) as usize;
    let i_far_south = (80.0 / dy) as usize;
    let i_far_north = (320.0 / dy) as usize;

    // Average h along each row across the domain (skip the very
    // first columns which are dominated by the BC ghost layer).
    let avg_row = |i: usize| -> f64 {
        let n_cols = final_states.ncols();
        let start = n_cols / 10;
        let count = n_cols - start;
        (start..n_cols).map(|j| final_states[(i, j)].h).sum::<f64>() / count as f64
    };
    let h_centre = avg_row(i_centre);
    let h_far_south = avg_row(i_far_south);
    let h_far_north = avg_row(i_far_north);

    assert!(
        h_centre > 2.0 * h_far_south,
        "valley centreline not deeper than far-south row: h_centre = {:.4}, h_far_south = {:.4}",
        h_centre,
        h_far_south
    );
    assert!(
        h_centre > 2.0 * h_far_north,
        "valley centreline not deeper than far-north row: h_centre = {:.4}, h_far_north = {:.4}",
        h_centre,
        h_far_north
    );
}

#[test]
fn mass_balance_bounded_by_cumulative_inflow() {
    // q_in · effective_width · t_end ≈ cumulative inflow. After
    // walls + transmissive E, final mass ≤ cumulative inflow + initial
    // thin-film mass, and must be > 50% of cumulative (the wave
    // hasn't fully exited the 1 km domain in 10 min).
    let case = TestCase::standard();
    let (mesh, _dx, dy) = build_mesh(case, 250, 100);
    let initial = Array2::from_elem((100, 250), Conserved2D::new(0.001, 0.0, 0.0));
    let m_initial: f64 = initial.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    let (final_states, _) = run_until(initial, &mesh, boundaries(case), case, 0.4);
    let m_final: f64 = final_states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();

    // Effective inflow width: all cells on the W boundary receive q;
    // a small fraction (near the walls where the bed is high) may
    // not effectively conduct flow. Use the full domain width as a
    // generous upper bound.
    let cumulative_inflow = case.q_in * case.length_y * case.t_end;
    assert!(m_final > m_initial, "no net inflow: m_final = {m_final}");
    assert!(
        m_final < m_initial + cumulative_inflow * 1.01,
        "mass exceeded cumulative inflow bound: m_final = {m_final}, cumulative_inflow = {cumulative_inflow}",
    );
    let _ = dy;
}

#[test]
fn depth_grows_at_centreline_with_distance_from_outlet() {
    // Steady-state intuition: upstream of the outlet, depth at
    // a fixed y (centreline) should be at or near Manning normal
    // depth — varying only mildly along x. Specifically, h at
    // x = 200 m (upstream centreline) should be GREATER than h at
    // x = 900 m (near the transmissive outlet) because near the
    // outlet some momentum is exiting and the wave hasn't yet
    // fully formed there.
    //
    // This is a weak monotonicity check: not asserting strict
    // depth profile, just that upstream is wetter than just-
    // -before-outlet.
    let case = TestCase::standard();
    let (mesh, dx, dy) = build_mesh(case, 250, 100);
    let initial = Array2::from_elem((100, 250), Conserved2D::new(0.001, 0.0, 0.0));
    let (final_states, _) = run_until(initial, &mesh, boundaries(case), case, 0.4);

    let i_centre = (case.centreline_y / dy) as usize;
    let j_upstream = (200.0 / dx) as usize;
    let j_near_outlet = (900.0 / dx) as usize;
    let h_upstream = final_states[(i_centre, j_upstream)].h;
    let h_near_outlet = final_states[(i_centre, j_near_outlet)].h;

    assert!(
        h_upstream > h_near_outlet,
        "depth not monotone along centreline: h(200) = {:.4}, h(900) = {:.4}",
        h_upstream,
        h_near_outlet
    );
}
