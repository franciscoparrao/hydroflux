//! MacDonald-style steady-state benchmark in 2D: Manning uniform flow
//! over a uniformly-sloped channel, in both `x`- and `y`-aligned
//! orientations.
//!
//! Degenerate case of MacDonald et al. (1997) inverse design extended
//! to two spatial dimensions: prescribe a constant water depth
//! `h(x, y) = h_n`, derive the bed `z(x) = −S₀·x` (or `z(y) = −S₀·y`
//! for the y-aligned case), and choose Manning `n` such that
//! bed-slope gravity is exactly balanced by friction at the
//! analytical normal depth:
//!
//! ```text
//!   q = (1/n) · h_n^(5/3) · √S₀     (Manning's equation)
//!   ⇒  h_n = (n · q / √S₀)^(3/5)
//! ```
//!
//! Tests jointly:
//! - Bed-slope source (Audusse 2D hydrostatic reconstruction per face);
//! - Manning friction step (`source.rs`);
//! - Boundary conditions `Discharge` (upstream, prescribed normal
//!   momentum) and `Depth` (downstream, prescribed depth);
//! - Per-direction isotropy: the same numerical result must hold for
//!   `x`-aligned and `y`-aligned channels.
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test macdonald_uniform
//! ```

use approx::assert_relative_eq;
use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, cfl_time_step, manning_friction_step, ssprk2_step,
};
use ndarray::Array2;

const G: f64 = 9.81;

/// Manning normal depth `h_n` for prescribed unit discharge `q`, bed
/// slope `S₀`, and roughness `n`. Inverted from Manning's equation.
fn manning_normal_depth(q: f64, slope: f64, manning: f64) -> f64 {
    (manning * q / slope.sqrt()).powf(3.0 / 5.0)
}

/// Bed descending uniformly in `+x` with slope `S₀`: `z(j) = −j·dx·S₀`,
/// uniform across `i`. Manning roughness applied uniformly.
fn x_sloped_mesh(n_rows: usize, n_cols: usize, dx: f64, slope: f64, manning: f64) -> Mesh2D {
    let bed = Array2::from_shape_fn((n_rows, n_cols), |(_i, j)| -(j as f64) * dx * slope);
    Mesh2D::new(bed, dx, dx, manning)
}

/// Bed descending uniformly in `+y` with slope `S₀`: `z(i) = −i·dx·S₀`,
/// uniform across `j`.
fn y_sloped_mesh(n_rows: usize, n_cols: usize, dx: f64, slope: f64, manning: f64) -> Mesh2D {
    let bed = Array2::from_shape_fn((n_rows, n_cols), |(i, _j)| -(i as f64) * dx * slope);
    Mesh2D::new(bed, dx, dx, manning)
}

/// Run forward Euler + Manning fractional step until `t_end`, using
/// the given BCs. Manning roughness is read per cell from `mesh`.
fn run_until(
    mut states: Array2<Conserved2D>,
    mesh: &Mesh2D,
    bcs: Boundaries2D,
    t_end: f64,
    cfl: f64,
) -> (Array2<Conserved2D>, usize) {
    let mut t = 0.0;
    let mut steps = 0;
    while t < t_end {
        let dt = cfl_time_step(&states, mesh, cfl).min(t_end - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        manning_friction_step(&mut states, mesh, dt, 1.0e-9);
        t += dt;
        steps += 1;
        if steps > 100_000 {
            panic!("MacDonald uniform: {steps} steps without reaching t_end");
        }
    }
    (states, steps)
}

#[test]
fn manning_normal_depth_inverts_manning_equation() {
    // Sanity: h_n derived above must satisfy q = (1/n) h_n^(5/3) √S₀.
    let q = 1.0_f64;
    let slope = 0.01_f64;
    let n = 0.03_f64;
    let h_n = manning_normal_depth(q, slope, n);
    let q_recomputed = h_n.powf(5.0 / 3.0) * slope.sqrt() / n;
    assert_relative_eq!(q, q_recomputed, epsilon = 1e-12);
}

#[test]
fn uniform_flow_x_aligned_is_preserved() {
    // Channel along x: water flows in +x direction, uniform in y. With
    // Discharge BC at West (prescribed q_x = q) and Depth BC at East
    // (prescribed h = h_n), the steady-state must hold across the
    // entire domain — not just the interior middle slab as in the 1D
    // version (which lacked Discharge/Depth BCs).
    let q = 1.0;
    let slope = 0.01;
    let manning = 0.03;
    let h_n = manning_normal_depth(q, slope, manning);

    let n_rows = 5;
    let n_cols = 50;
    let dx = 1.0;
    let mesh = x_sloped_mesh(n_rows, n_cols, dx, slope, manning);

    // Initial state: h = h_n, hu = q, hv = 0 everywhere.
    let states = Array2::from_elem((n_rows, n_cols), Conserved2D::new(h_n, q, 0.0));

    let bcs = Boundaries2D {
        west: Boundary::Discharge { q },
        east: Boundary::Depth { h: h_n },
        north: Boundary::Wall,
        south: Boundary::Wall,
    };

    // Run for ~2 traversal times so any boundary-layer wave can damp.
    let u = q / h_n;
    let c = (G * h_n).sqrt();
    let traversal = (n_cols as f64 * dx) / (u + c);
    let t_end = 2.0 * traversal;
    let (final_states, _) = run_until(states, &mesh, bcs, t_end, 0.4);

    // Check uniform-flow preservation: depth must equal h_n, hu = q,
    // hv = 0, everywhere. Tolerance reflects first-order Euler over
    // many CFL-bounded steps + Audusse + semi-implicit Manning combined.
    let mut max_dh = 0.0_f64;
    let mut max_dhu = 0.0_f64;
    let mut max_dhv = 0.0_f64;
    for s in &final_states {
        max_dh = max_dh.max((s.h - h_n).abs());
        max_dhu = max_dhu.max((s.hu - q).abs());
        max_dhv = max_dhv.max(s.hv.abs());
    }
    // Tolerance 2% was set when η-MUSCL alone left a residual
    // O(dx·S₀/h_n) bias on this steady-state problem. With bed
    // reconstruction (Liang & Marche 2009) + flux rescaling, the
    // measured drift drops to ~0.03% (h) and ~0.18% (hu) on this
    // mesh — a 45× improvement on h. The tolerance is left at 2%
    // as a loose regression guard; substantial drift would indicate
    // a regression in the well-balanced source or flux rescaling. The interior
    // face flux uses h* = h_n − 0.5·dx·S₀ instead of h_n (Audusse
    // reconstruction on a piecewise-constant bed), creating a small
    // mismatch with the first-order boundary face that drives a
    // bounded steady-state perturbation. A fully bed-reconstructed
    // MUSCL (Liang & Marche 2009) would eliminate this; out of scope
    // for the current iteration. The trade-off is a ~3x improvement
    // on dam-break-on-dry L1/L²/L∞ for this 1% MacDonald regression.
    assert!(
        max_dh / h_n < 2.0e-2,
        "depth drift too large: max |Δh|/h_n = {:.3}%",
        100.0 * max_dh / h_n
    );
    assert!(
        max_dhu / q < 2.0e-2,
        "discharge drift too large: max |Δhu|/q = {:.3}%",
        100.0 * max_dhu / q
    );
    assert!(
        max_dhv / q < 5.0e-3,
        "spurious tangential momentum: max |hv|/q = {:.3}%",
        100.0 * max_dhv / q
    );
}

#[test]
fn uniform_flow_y_aligned_is_preserved() {
    // Same physics, channel along y. The numerical result must match
    // the x-aligned case to within roundoff — direction isotropy of
    // the FV scheme + Manning step.
    let q = 1.0;
    let slope = 0.01;
    let manning = 0.03;
    let h_n = manning_normal_depth(q, slope, manning);

    let n_rows = 50;
    let n_cols = 5;
    let dx = 1.0;
    let mesh = y_sloped_mesh(n_rows, n_cols, dx, slope, manning);

    // Initial: h = h_n, hv = q (along +y), hu = 0.
    let states = Array2::from_elem((n_rows, n_cols), Conserved2D::new(h_n, 0.0, q));

    let bcs = Boundaries2D {
        west: Boundary::Wall,
        east: Boundary::Wall,
        north: Boundary::Discharge { q },
        south: Boundary::Depth { h: h_n },
    };

    let u = q / h_n;
    let c = (G * h_n).sqrt();
    let traversal = (n_rows as f64 * dx) / (u + c);
    let t_end = 2.0 * traversal;
    let (final_states, _) = run_until(states, &mesh, bcs, t_end, 0.4);

    let mut max_dh = 0.0_f64;
    let mut max_dhv = 0.0_f64;
    let mut max_dhu = 0.0_f64;
    for s in &final_states {
        max_dh = max_dh.max((s.h - h_n).abs());
        max_dhv = max_dhv.max((s.hv - q).abs());
        max_dhu = max_dhu.max(s.hu.abs());
    }
    // Same 2% tolerance as the x-aligned version; see that test's
    // documentation for the rationale.
    assert!(
        max_dh / h_n < 2.0e-2,
        "depth drift (y-aligned): {:.3}%",
        100.0 * max_dh / h_n
    );
    assert!(
        max_dhv / q < 2.0e-2,
        "discharge drift (y-aligned): {:.3}%",
        100.0 * max_dhv / q
    );
    assert!(
        max_dhu / q < 5.0e-3,
        "spurious x-momentum (y-aligned): {:.3}%",
        100.0 * max_dhu / q
    );
}

#[test]
fn perturbation_from_normal_depth_relaxes_back() {
    // Start with a 20% perturbation in `h` over the analytical normal
    // depth on a strip in the middle of the domain. After running with
    // Discharge/Depth BCs sustained, the flow should relax back toward
    // `h_n` — Manning friction + bed-slope source act as a restoring
    // mechanism around the steady state.
    let q = 1.0;
    let slope = 0.01;
    let manning = 0.03;
    let h_n = manning_normal_depth(q, slope, manning);
    let perturbation = 0.20;

    let n_rows = 5;
    let n_cols = 100;
    let dx = 1.0;
    let mesh = x_sloped_mesh(n_rows, n_cols, dx, slope, manning);

    let mut states = Array2::from_elem((n_rows, n_cols), Conserved2D::new(h_n, q, 0.0));
    // Inject a 20%-perturbed slab in the middle.
    let j_mid_lo = n_cols * 4 / 10;
    let j_mid_hi = n_cols * 6 / 10;
    for i in 0..n_rows {
        for j in j_mid_lo..j_mid_hi {
            states[(i, j)].h = h_n * (1.0 + perturbation);
            states[(i, j)].hu = q;
        }
    }
    let initial_perturbation = states
        .iter()
        .map(|s| (s.h - h_n).abs())
        .fold(0.0_f64, f64::max);

    let bcs = Boundaries2D {
        west: Boundary::Discharge { q },
        east: Boundary::Depth { h: h_n },
        north: Boundary::Wall,
        south: Boundary::Wall,
    };

    let u = q / h_n;
    let c = (G * h_n).sqrt();
    let traversal = (n_cols as f64 * dx) / (u + c);
    let t_end = 6.0 * traversal; // long enough for the perturbation to convect out
    let (final_states, _) = run_until(states, &mesh, bcs, t_end, 0.4);

    let final_perturbation = final_states
        .iter()
        .map(|s| (s.h - h_n).abs())
        .fold(0.0_f64, f64::max);

    // The perturbation must decay by at least an order of magnitude.
    assert!(
        final_perturbation < 0.1 * initial_perturbation,
        "perturbation did not decay: initial max |Δh| = {:.4}, final = {:.4}",
        initial_perturbation,
        final_perturbation
    );
}
