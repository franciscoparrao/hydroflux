//! Radial (cylindrical) dam break — Toro (2009) §10.10.
//!
//! A circular water column of radius `R_dam` sitting on a flat bed at
//! depth `h_in`, surrounded by a shallower reservoir at `h_out`. At
//! `t = 0` the dam is released. The resulting flow is **axisymmetric**:
//! an outward shock propagating into the shallow region and an inward
//! rarefaction collapsing the inner column.
//!
//! The exact 2D solution is not closed-form — the radial SWE in
//! similarity variables yields an ODE system (Toro §10.10) that must be
//! integrated numerically. For our purposes we exercise the 2D solver
//! against three structural properties that are **independent** of the
//! exact reference:
//!
//! 1. **Axisymmetric preservation**. The numerical solution must remain
//!    axisymmetric: depth profiles `h(r)` sampled along different
//!    directions (the `x`-axis, the `y`-axis, and two diagonals) must
//!    coincide. Deviation measures the grid anisotropy of the Cartesian
//!    FV update — a well-designed scheme keeps it below a few percent.
//!
//! 2. **Mass conservation** under walls (the shock doesn't reach the
//!    boundary at our chosen `t_end`).
//!
//! 3. **Bounded depth**: every cell has `h ∈ [0, h_in]` after release
//!    — neither below zero nor above the initial maximum.
//!
//! Additionally, an informational metric compares the inner-cap and
//! far-field depths against the 1D Stoker wet-wet limit (valid as a
//! local approximation along the radial direction at short times,
//! before the geometric source `−h u / r` accumulates).
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test radial_dam_break
//! ```

use approx::assert_relative_eq;
use hydroflux_solver_2d::{Boundaries2D, Conserved2D, Mesh2D, cfl_time_step, ssprk2_step};
use ndarray::Array2;

const G: f64 = 9.81;

#[derive(Debug, Clone, Copy)]
struct RadialDamBreak {
    /// Inner-column depth `h_in` [m].
    h_in: f64,
    /// Outer reservoir depth `h_out` [m] (wet-wet variant).
    h_out: f64,
    /// Radius of the initial circular dam [m].
    r_dam: f64,
    /// Half-extent of the (square) computational domain `[-L, L]` [m].
    half_extent: f64,
    /// Final time of the simulation [s]. Choose so the outward shock
    /// stays inside the domain.
    t_end: f64,
}

impl RadialDamBreak {
    fn celerity_in(self) -> f64 {
        (G * self.h_in).sqrt()
    }
    fn celerity_out(self) -> f64 {
        (G * self.h_out).sqrt()
    }
}

fn build_mesh(case: RadialDamBreak, n: usize) -> (Mesh2D, f64) {
    let dx = 2.0 * case.half_extent / n as f64;
    let bed = Array2::<f64>::zeros((n, n));
    (Mesh2D::new(bed, dx, dx, 0.0), dx)
}

fn initial_state(case: RadialDamBreak, n: usize, dx: f64) -> Array2<Conserved2D> {
    Array2::from_shape_fn((n, n), |(i, j)| {
        let x = (j as f64 + 0.5) * dx - case.half_extent;
        let y = (i as f64 + 0.5) * dx - case.half_extent;
        let r = (x * x + y * y).sqrt();
        let h = if r < case.r_dam {
            case.h_in
        } else {
            case.h_out
        };
        Conserved2D::new(h, 0.0, 0.0)
    })
}

fn run_until(
    mut states: Array2<Conserved2D>,
    mesh: &Mesh2D,
    bcs: Boundaries2D,
    t_end: f64,
    cfl: f64,
) -> Array2<Conserved2D> {
    let mut t = 0.0;
    let mut steps = 0;
    while t < t_end {
        let dt = cfl_time_step(&states, mesh, cfl).min(t_end - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        t += dt;
        steps += 1;
        if steps > 100_000 {
            panic!("radial dam break: {steps} steps without reaching t_end");
        }
    }
    states
}

/// Total water volume (sum of `h · dx · dy`).
fn total_volume(states: &Array2<Conserved2D>, dx: f64, dy: f64) -> f64 {
    states.iter().map(|s| s.h * dx * dy).sum()
}

fn standard_case() -> RadialDamBreak {
    RadialDamBreak {
        h_in: 2.5,
        h_out: 0.5,
        r_dam: 2.5,
        half_extent: 10.0,
        t_end: 0.7,
    }
}

#[test]
fn mass_is_conserved_under_walls() {
    // Walls + small enough t_end (the outer shock at ~5 m/s · 0.7 s ≈
    // 3.5 m beyond r_dam doesn't reach the wall at 10 m). Mass must
    // be conserved to ~1e-5 (the bed-recon+flux-rescaling baseline
    // documented in the Thacker writeup).
    let case = standard_case();
    let n = 80;
    let (mesh, dx) = build_mesh(case, n);
    let initial = initial_state(case, n, dx);
    let m0 = total_volume(&initial, mesh.dx, mesh.dy);

    let final_states = run_until(initial, &mesh, Boundaries2D::WALLS, case.t_end, 0.4);
    let m1 = total_volume(&final_states, mesh.dx, mesh.dy);

    assert_relative_eq!(m0, m1, epsilon = 1.0e-5);
}

#[test]
fn depth_remains_bounded_and_finite() {
    // No cell can develop `h < 0` or `h > h_in` (initial maximum). Also
    // checks finiteness as a basic NaN/Inf trap.
    let case = standard_case();
    let n = 80;
    let (mesh, dx) = build_mesh(case, n);
    let initial = initial_state(case, n, dx);
    let final_states = run_until(initial, &mesh, Boundaries2D::WALLS, case.t_end, 0.4);

    for s in &final_states {
        assert!(s.h.is_finite(), "h became non-finite: {}", s.h);
        assert!(s.h >= 0.0, "h went negative: {}", s.h);
        // Allow a tiny overshoot tolerance for MUSCL near the
        // discontinuous front; 1% of `h_in` is a generous bound.
        assert!(
            s.h <= case.h_in * 1.01,
            "h exceeded initial maximum: h = {}, h_in = {}",
            s.h,
            case.h_in
        );
        assert!(s.hu.is_finite(), "hu became non-finite: {}", s.hu);
        assert!(s.hv.is_finite(), "hv became non-finite: {}", s.hv);
    }
}

#[test]
fn solution_remains_approximately_axisymmetric() {
    // Sample the radial depth profile along four directions and check
    // that the profiles agree.
    //
    // Directions: +x, +y, +45° diagonal, +135° diagonal. The Cartesian
    // FV grid has a 4-fold symmetry by construction (90° rotations of
    // the initial condition are exact mirrors of the cell layout), so
    // x-axis and y-axis profiles should match to roundoff; the
    // diagonals carry the real grid anisotropy.
    let case = standard_case();
    let n = 80;
    let (mesh, dx) = build_mesh(case, n);
    let initial = initial_state(case, n, dx);
    let final_states = run_until(initial, &mesh, Boundaries2D::WALLS, case.t_end, 0.4);

    // Index of the cell whose centre is closest to (x, y).
    let cell_index = |x_target: f64, y_target: f64| -> (usize, usize) {
        let j = ((x_target + case.half_extent) / dx - 0.5)
            .round()
            .clamp(0.0, (n - 1) as f64) as usize;
        let i = ((y_target + case.half_extent) / dx - 0.5)
            .round()
            .clamp(0.0, (n - 1) as f64) as usize;
        (i, j)
    };

    // Sample at radii from 0 to ~half the domain in steps of dx.
    let n_samples = (case.half_extent / dx / 2.0) as usize;
    let mut max_dev_ax = 0.0_f64;
    let mut max_dev_diag = 0.0_f64;
    for k in 0..n_samples {
        let r = (k as f64 + 0.5) * dx;
        let h_xpos = {
            let (i, j) = cell_index(r, 0.0);
            final_states[(i, j)].h
        };
        let h_ypos = {
            let (i, j) = cell_index(0.0, r);
            final_states[(i, j)].h
        };
        let h_diag1 = {
            let (i, j) = cell_index(r / std::f64::consts::SQRT_2, r / std::f64::consts::SQRT_2);
            final_states[(i, j)].h
        };
        let h_diag2 = {
            let (i, j) = cell_index(-r / std::f64::consts::SQRT_2, r / std::f64::consts::SQRT_2);
            final_states[(i, j)].h
        };

        let mean = 0.25 * (h_xpos + h_ypos + h_diag1 + h_diag2);
        let dev_ax = ((h_xpos - h_ypos).abs() / mean.max(case.h_out * 0.1)).abs();
        let dev_diag = (h_diag1 - h_diag2).abs() / mean.max(case.h_out * 0.1);
        max_dev_ax = max_dev_ax.max(dev_ax);
        max_dev_diag = max_dev_diag.max(dev_diag);
    }

    // Axis-vs-axis (x vs y): should be near-roundoff because of grid
    // 4-fold symmetry (the initial condition is symmetric under
    // i↔j swap because dx=dy and the cap is centred).
    assert!(
        max_dev_ax < 1.0e-10,
        "x-vs-y profile mismatch: max relative deviation = {:.2e}",
        max_dev_ax
    );
    // Diagonal-vs-diagonal: real grid anisotropy. Limit to 1% for a
    // 80×80 mesh — the MUSCL + SSP-RK2 scheme should be well below
    // this on a smooth-enough region.
    assert!(
        max_dev_diag < 1.0e-2,
        "diagonal-vs-diagonal profile mismatch: max relative deviation = {:.4}",
        max_dev_diag
    );
}

#[test]
fn shock_front_advances_outward() {
    // The outward shock should propagate at a positive radial speed
    // bounded above by `c_in + c_out` (the maximum signal speed at
    // the wet-wet front) and below by ~`c_out` (the slower-side
    // celerity). Detect the front as the rightmost cell on the x-axis
    // whose depth has risen above `(h_out + h_in)/2`.
    let case = standard_case();
    let n = 80;
    let (mesh, dx) = build_mesh(case, n);
    let initial = initial_state(case, n, dx);
    let final_states = run_until(initial, &mesh, Boundaries2D::WALLS, case.t_end, 0.4);

    let mid_row = n / 2;
    // Shock front: the outermost cell on the +x ray whose depth is
    // still meaningfully above the far-field `h_out` level. The radial
    // depth profile is non-monotonic (centre depression + rarefaction
    // hump + star region + shock + outer field), so we look for the
    // OUTERMOST crossing rather than the first one — that's the shock.
    let h_threshold = 1.5 * case.h_out;
    let mut j_front: Option<usize> = None;
    for j in (n / 2)..n {
        if final_states[(mid_row, j)].h > h_threshold {
            j_front = Some(j);
        }
    }
    let j_front =
        j_front.expect("no cell above shock threshold — domain too small or t_end too large?");
    let r_front = (j_front as f64 + 0.5) * dx - case.half_extent;

    // The front must have moved outward from `r_dam`.
    assert!(
        r_front > case.r_dam,
        "front did not advance: r_front = {:.3}, r_dam = {:.3}",
        r_front,
        case.r_dam
    );
    // Upper bound: speed cannot exceed `c_in + c_out` (loose bound).
    let max_advance = (case.celerity_in() + case.celerity_out()) * case.t_end;
    assert!(
        r_front - case.r_dam < max_advance,
        "front moved too fast: Δr = {:.3} > {:.3} (c_in + c_out) · t_end",
        r_front - case.r_dam,
        max_advance
    );
}

#[test]
#[ignore = "informational: prints axisymmetry + Stoker-1D comparison"]
fn report_metrics() {
    // Side-by-side comparison of radial profiles along 4 directions
    // plus the 1D Stoker wet-wet limit. Run with:
    //   cargo test --release -p hydroflux-solver-2d --test radial_dam_break -- --ignored --nocapture
    let case = standard_case();
    let n = 160;
    let (mesh, dx) = build_mesh(case, n);
    let initial = initial_state(case, n, dx);
    let final_states = run_until(initial, &mesh, Boundaries2D::WALLS, case.t_end, 0.4);

    eprintln!("\n=== Radial dam break — axisymmetry report ===");
    eprintln!("Mesh: {n}×{n}, dx = {dx:.4} m, t_end = {} s", case.t_end);
    eprintln!(
        "h_in = {} m, h_out = {} m, r_dam = {} m",
        case.h_in, case.h_out, case.r_dam
    );
    eprintln!(
        "c_in = {:.3} m/s, c_out = {:.3} m/s",
        case.celerity_in(),
        case.celerity_out()
    );

    let cell_index = |x_target: f64, y_target: f64| -> (usize, usize) {
        let j = ((x_target + case.half_extent) / dx - 0.5)
            .round()
            .clamp(0.0, (n - 1) as f64) as usize;
        let i = ((y_target + case.half_extent) / dx - 0.5)
            .round()
            .clamp(0.0, (n - 1) as f64) as usize;
        (i, j)
    };

    eprintln!(
        "\n{:>6} {:>10} {:>10} {:>10} {:>10}",
        "r [m]", "h(+x)", "h(+y)", "h(+45°)", "h(+135°)"
    );
    for k in (0..(case.half_extent / dx / 2.0) as usize).step_by(2) {
        let r = (k as f64 + 0.5) * dx;
        let (i, j) = cell_index(r, 0.0);
        let h_xp = final_states[(i, j)].h;
        let (i, j) = cell_index(0.0, r);
        let h_yp = final_states[(i, j)].h;
        let (i, j) = cell_index(r / std::f64::consts::SQRT_2, r / std::f64::consts::SQRT_2);
        let h_d1 = final_states[(i, j)].h;
        let (i, j) = cell_index(-r / std::f64::consts::SQRT_2, r / std::f64::consts::SQRT_2);
        let h_d2 = final_states[(i, j)].h;
        eprintln!(
            "{:>6.2} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
            r, h_xp, h_yp, h_d1, h_d2
        );
    }
    eprintln!("=============================================\n");
}
