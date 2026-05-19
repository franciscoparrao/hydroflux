//! Planar oscillation in a paraboloidal basin (Thacker 1981) — the
//! canonical 2D analytical benchmark for shallow-water solvers.
//!
//! The water column is a parabolic cap of horizontal radius `a` whose
//! centre orbits the basin axis at radius `B` with angular frequency
//! `ω = √(2 g h₀) / a`. The free-surface profile remains a paraboloid
//! at all times; the centre of mass executes uniform circular motion;
//! the wet region is a closed disk that translates around the basin
//! without changing shape.
//!
//! References:
//! - Thacker (1981), J. Fluid Mech. 107, 499–508. Original derivation.
//! - Sampson, Easton & Singh (2006), ANZIAM J. 47, C373–C387.
//!   Convenient parameterisation used here.
//!
//! Reproducir: `cargo test --test thacker --release`.
//!
//! # Setup
//!
//! Bed: `z_b(x, y) = h₀ · ((x² + y²)/a² − 1)`. Bowl minimum at the
//! origin (`z_b(0, 0) = −h₀`); rim level `z_b = 0` reached at radius
//! `a`.
//!
//! Analytical depth:
//!
//! ```text
//!   h(x, y, t) = (h₀ / a²) · max(0, a² − (x − B cos ωt)² − (y − B sin ωt)²)
//! ```
//!
//! Analytical velocities (spatially uniform throughout the wet region):
//!
//! ```text
//!   u(t) = −B ω sin(ωt)
//!   v(t) =  B ω cos(ωt)
//! ```
//!
//! Total volume `V = π h₀ a² / 2` (paraboloid of revolution), invariant
//! in time.

use approx::assert_relative_eq;
use hydroflux_solver_2d::{Boundaries2D, Conserved2D, Mesh2D, cfl_time_step, forward_euler_step};
use ndarray::Array2;

const G: f64 = 9.81;

/// Parameters of the planar Thacker oscillation.
#[derive(Debug, Clone, Copy)]
struct Thacker {
    /// Bowl depth scale [m]. Bed minimum is at `z = −h0`.
    h0: f64,
    /// Bowl horizontal radius [m]. Rim at `r = a`.
    a: f64,
    /// Orbit radius of the water cap centre [m]. `B < a` for a confined
    /// oscillation that never escapes the basin.
    b: f64,
}

impl Thacker {
    fn omega(self) -> f64 {
        (2.0 * G * self.h0).sqrt() / self.a
    }

    fn period(self) -> f64 {
        2.0 * std::f64::consts::PI / self.omega()
    }

    /// Analytical depth at `(x, y, t)`. Returns 0 outside the wet disk.
    fn depth(self, x: f64, y: f64, t: f64) -> f64 {
        let omega = self.omega();
        let xc = self.b * (omega * t).cos();
        let yc = self.b * (omega * t).sin();
        let r2 = (x - xc).powi(2) + (y - yc).powi(2);
        let raw = self.a * self.a - r2;
        if raw > 0.0 {
            self.h0 * raw / (self.a * self.a)
        } else {
            0.0
        }
    }

    /// Spatially uniform velocity `(u, v)` at time `t`.
    fn velocity(self, t: f64) -> (f64, f64) {
        let omega = self.omega();
        let u = -self.b * omega * (omega * t).sin();
        let v = self.b * omega * (omega * t).cos();
        (u, v)
    }

    /// Bed elevation `z_b(x, y)`.
    fn bed(self, x: f64, y: f64) -> f64 {
        self.h0 * ((x * x + y * y) / (self.a * self.a) - 1.0)
    }

    /// Analytical total volume of the cap. Constant in time.
    fn total_volume(self) -> f64 {
        std::f64::consts::PI * self.h0 * self.a * self.a / 2.0
    }
}

/// Build a square mesh centred on the origin with `n` cells per side
/// over the domain `[-half_extent, half_extent]²`.
fn centred_mesh(thacker: Thacker, n: usize, half_extent: f64) -> (Mesh2D, f64) {
    let dx = 2.0 * half_extent / n as f64;
    let bed = Array2::from_shape_fn((n, n), |(i, j)| {
        let x = -half_extent + (j as f64 + 0.5) * dx;
        let y = -half_extent + (i as f64 + 0.5) * dx;
        thacker.bed(x, y)
    });
    (Mesh2D::new(bed, dx, dx, 0.0), dx)
}

/// Build initial conditions from the analytical solution at `t = 0`.
fn initial_state(thacker: Thacker, n: usize, half_extent: f64, dx: f64) -> Array2<Conserved2D> {
    let (u0, v0) = thacker.velocity(0.0);
    Array2::from_shape_fn((n, n), |(i, j)| {
        let x = -half_extent + (j as f64 + 0.5) * dx;
        let y = -half_extent + (i as f64 + 0.5) * dx;
        let h = thacker.depth(x, y, 0.0);
        Conserved2D::new(h, h * u0, h * v0)
    })
}

/// Total water volume (sum of `h · dx · dy`).
fn total_volume(states: &Array2<Conserved2D>, dx: f64, dy: f64) -> f64 {
    states.iter().map(|s| s.h * dx * dy).sum()
}

/// Centroid of the wet region weighted by depth: `(∑ x·h, ∑ y·h) / ∑ h`.
/// Returns `None` if the total mass is below `min_mass`.
fn depth_weighted_centroid(
    states: &Array2<Conserved2D>,
    half_extent: f64,
    dx: f64,
    min_mass: f64,
) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut m = 0.0;
    for ((i, j), s) in states.indexed_iter() {
        let x = -half_extent + (j as f64 + 0.5) * dx;
        let y = -half_extent + (i as f64 + 0.5) * dx;
        sx += x * s.h;
        sy += y * s.h;
        m += s.h;
    }
    if m > min_mass {
        Some((sx / m, sy / m))
    } else {
        None
    }
}

/// Run the simulation until `t_end`, returning the final state and the
/// number of steps taken. Uses CFL-bounded forward Euler; walls on all
/// four sides (the wet disk never reaches the boundary at the chosen
/// parameters).
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
        let dt_cfl = cfl_time_step(&states, mesh, cfl);
        let dt = dt_cfl.min(t_end - t);
        forward_euler_step(&mut states, mesh, bcs, dt);
        t += dt;
        steps += 1;
        if steps > 100_000 {
            panic!("aborting Thacker run: {steps} steps without reaching t_end");
        }
    }
    (states, steps)
}

#[test]
fn mass_is_conserved_under_wall_boundaries() {
    // The wet disk stays well inside the domain for B = 0.1 a, so wall
    // BCs do not lose mass to outflow. Manning friction is 0; the only
    // source of mass change should be roundoff in the FV update.
    let thacker = Thacker {
        h0: 0.1,
        a: 1.0,
        b: 0.1,
    };
    let n = 50;
    let half_extent = 1.25;
    let (mesh, dx) = centred_mesh(thacker, n, half_extent);
    let states = initial_state(thacker, n, half_extent, dx);

    let v_initial_numerical = total_volume(&states, mesh.dx, mesh.dy);
    let t_end = thacker.period() / 4.0;
    let (final_states, _steps) = run_until(states, &mesh, Boundaries2D::WALLS, t_end, 0.4);
    let v_final_numerical = total_volume(&final_states, mesh.dx, mesh.dy);

    // Wall boundaries: numerical mass must be conserved to roundoff.
    assert_relative_eq!(v_initial_numerical, v_final_numerical, epsilon = 1e-10);
}

#[test]
fn initial_volume_matches_analytical_to_quadrature_error() {
    // The initial cell-centred sum approximates the parabolic-cap
    // volume π·h₀·a²/2. With a 50×50 mesh over [−1.25, 1.25]² the
    // quadrature error is dominated by the discretisation of the
    // circular wet boundary (a few %); we just sanity-check that the
    // sign and order of magnitude are right.
    let thacker = Thacker {
        h0: 0.1,
        a: 1.0,
        b: 0.1,
    };
    let n = 50;
    let half_extent = 1.25;
    let (mesh, dx) = centred_mesh(thacker, n, half_extent);
    let states = initial_state(thacker, n, half_extent, dx);

    let v_numerical = total_volume(&states, mesh.dx, mesh.dy);
    let v_analytical = thacker.total_volume();

    let rel_err = (v_numerical - v_analytical).abs() / v_analytical;
    assert!(
        rel_err < 0.05,
        "initial-volume discretisation error too large: {:.2}% (numerical {:.6}, analytical {:.6})",
        rel_err * 100.0,
        v_numerical,
        v_analytical,
    );
}

#[test]
fn depth_remains_non_negative_and_finite() {
    // Robustness check across one quarter-period of evolution. Catches
    // dry-bed underflow, NaN propagation, or runaway momentum in the
    // wet/dry transition cells at the disk perimeter.
    let thacker = Thacker {
        h0: 0.1,
        a: 1.0,
        b: 0.1,
    };
    let n = 50;
    let half_extent = 1.25;
    let (mesh, dx) = centred_mesh(thacker, n, half_extent);
    let states = initial_state(thacker, n, half_extent, dx);

    let t_end = thacker.period() / 4.0;
    let (final_states, _) = run_until(states, &mesh, Boundaries2D::WALLS, t_end, 0.4);

    for s in &final_states {
        assert!(s.h.is_finite(), "depth became non-finite: {}", s.h);
        assert!(s.h >= 0.0, "depth went negative: {}", s.h);
        assert!(s.hu.is_finite(), "hu became non-finite: {}", s.hu);
        assert!(s.hv.is_finite(), "hv became non-finite: {}", s.hv);
    }
}

#[test]
fn centroid_executes_circular_motion_at_analytical_frequency() {
    // Track the depth-weighted centroid of the cap over one full
    // period. It must trace approximately a circle of radius `B` and
    // return to its initial position at `t = T` (orbit closed).
    let thacker = Thacker {
        h0: 0.1,
        a: 1.0,
        b: 0.1,
    };
    let n = 50;
    let half_extent = 1.25;
    let (mesh, dx) = centred_mesh(thacker, n, half_extent);
    let states = initial_state(thacker, n, half_extent, dx);

    // Initial centroid (analytical): (B, 0).
    let min_mass = 0.1 * total_volume(&states, mesh.dx, mesh.dy);
    let (cx0, cy0) =
        depth_weighted_centroid(&states, half_extent, dx, min_mass).expect("nonzero mass");

    let t_end = thacker.period();
    let (final_states, steps) = run_until(states, &mesh, Boundaries2D::WALLS, t_end, 0.4);
    let (cx1, cy1) =
        depth_weighted_centroid(&final_states, half_extent, dx, min_mass).expect("nonzero mass");

    // After one full period the centroid should return to its starting
    // location. With a 50×50 mesh and first-order forward Euler + HLLC,
    // phase + amplitude error accumulates to ~30% of the orbit radius
    // `B` over one full period (measured: 0.030 m drift for B = 0.10 m).
    // The bound `0.35·B` accepts this and acts as a regression guard;
    // expected to tighten substantially once RK2 / MUSCL slope-limiting
    // are added (roadmap 2027 Q1).
    let drift = ((cx1 - cx0).powi(2) + (cy1 - cy0).powi(2)).sqrt();
    let tol = 0.35 * thacker.b;
    assert!(
        drift < tol,
        "centroid drift {drift:.4} m exceeds tolerance {tol:.4} m after {steps} steps (cx0={cx0:.4}, cy0={cy0:.4}, cx1={cx1:.4}, cy1={cy1:.4})"
    );

    // Sanity-check the initial centroid is close to (B, 0). Cell-size
    // limits the resolution; the analytical position is `(B, 0)` and
    // the numerical centroid should agree within ~half a cell.
    let initial_pos_err = ((cx0 - thacker.b).powi(2) + cy0.powi(2)).sqrt();
    assert!(
        initial_pos_err < dx,
        "initial centroid x={cx0:.4}, y={cy0:.4} not within dx={dx:.4} of analytical (B, 0)=({:.4}, 0.0)",
        thacker.b
    );
}

#[test]
fn lake_at_rest_is_preserved_on_paraboloidal_basin() {
    // Special case `B = 0`: the analytical solution collapses to a
    // perfectly axisymmetric cap at rest. Audusse must hold the
    // free surface exactly flat over the paraboloidal bed at all
    // times — this is the well-balanced property in 2D on a
    // smoothly curved bottom (not just a piecewise plane).
    let thacker = Thacker {
        h0: 0.1,
        a: 1.0,
        b: 0.0,
    };
    let n = 50;
    let half_extent = 1.25;
    let (mesh, dx) = centred_mesh(thacker, n, half_extent);
    let initial = initial_state(thacker, n, half_extent, dx);
    let mut states = initial.clone();

    // Identify "interior wet" cells: at least one cell away from the
    // wet/dry interface, so the assertion isolates the well-balanced
    // property from any dry-front noise.
    let h_threshold = 0.5 * thacker.h0;
    let interior_wet: Vec<(usize, usize)> = initial
        .indexed_iter()
        .filter_map(|((i, j), s)| (s.h > h_threshold).then_some((i, j)))
        .collect();
    assert!(!interior_wet.is_empty(), "no interior wet cells found");

    let t_end = thacker.period() / 4.0;
    for _ in 0..2000 {
        let dt_cfl = cfl_time_step(&states, &mesh, 0.4);
        let dt = dt_cfl.min(t_end);
        forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        if dt < 1e-15 {
            break;
        }
    }

    // For B = 0, all "well inside" cells should retain their initial
    // depth and stay at rest. Tolerance reflects numerical roundoff
    // accumulated over O(100) timesteps.
    for &(i, j) in &interior_wet {
        let diff = (states[(i, j)].h - initial[(i, j)].h).abs();
        assert!(
            diff < 1e-9,
            "interior wet cell ({i},{j}) drifted: h_initial={}, h_final={}, |Δh|={}",
            initial[(i, j)].h,
            states[(i, j)].h,
            diff
        );
        assert!(
            states[(i, j)].hu.abs() < 1e-9 && states[(i, j)].hv.abs() < 1e-9,
            "interior wet cell ({i},{j}) developed momentum: ({}, {})",
            states[(i, j)].hu,
            states[(i, j)].hv
        );
    }
}

#[test]
fn velocity_field_remains_approximately_uniform_through_quarter_period() {
    // The Thacker solution has spatially uniform `(u, v)` at all
    // times. After T/4 the analytical velocity is `(−Bω, 0)`. Average
    // the numerical velocity over deep (well-interior) wet cells and
    // compare. We do NOT expect a tight bound — the wet/dry front
    // contaminates cells near the perimeter and there is dispersive
    // error at the cap apex — but the depth-weighted average should
    // be within order-of-magnitude of the analytical.
    let thacker = Thacker {
        h0: 0.1,
        a: 1.0,
        b: 0.1,
    };
    let n = 50;
    let half_extent = 1.25;
    let (mesh, dx) = centred_mesh(thacker, n, half_extent);
    let states = initial_state(thacker, n, half_extent, dx);

    let t_end = thacker.period() / 4.0;
    let (final_states, _) = run_until(states, &mesh, Boundaries2D::WALLS, t_end, 0.4);

    // Restrict to deep cells (h > 0.5·h_max) to exclude wet/dry front noise.
    let h_max = final_states.iter().fold(0.0_f64, |a, s| a.max(s.h));
    let h_thresh = 0.5 * h_max;
    let mut sum_u = 0.0;
    let mut sum_v = 0.0;
    let mut sum_h = 0.0;
    for s in &final_states {
        if s.h > h_thresh {
            sum_u += s.hu;
            sum_v += s.hv;
            sum_h += s.h;
        }
    }
    let u_avg = sum_u / sum_h;
    let v_avg = sum_v / sum_h;

    // Analytical at t = T/4: u = -B·ω, v = 0.
    let u_an = -thacker.b * thacker.omega();
    let v_an = 0.0_f64;

    // Loose 30% tolerance on `u`, and an absolute bound of 0.2·|u_an|
    // on `v` (which is zero analytically).
    let u_err = (u_avg - u_an).abs() / u_an.abs();
    assert!(
        u_err < 0.30,
        "u (deep-cell depth-weighted avg) error {:.1}% too large: numerical {:.4}, analytical {:.4}",
        u_err * 100.0,
        u_avg,
        u_an
    );
    let v_err = (v_avg - v_an).abs() / u_an.abs();
    assert!(
        v_err < 0.20,
        "v error {:.1}% of |u_an| too large: numerical {:.4}, analytical {:.4}",
        v_err * 100.0,
        v_avg,
        v_an
    );
}

#[test]
#[ignore = "informational benchmark; run with --ignored to print error metrics"]
fn report_error_metrics_for_documentation() {
    // Not a pass/fail test — runs Thacker for T/2 on a moderately fine
    // mesh and prints metrics that the markdown benchmark in
    // `benchmarks/thacker-results.md` can quote. Invoke with:
    //   cargo test --release --test thacker -- --ignored --nocapture
    let thacker = Thacker {
        h0: 0.1,
        a: 1.0,
        b: 0.1,
    };
    let n = 80;
    let half_extent = 1.25;
    let (mesh, dx) = centred_mesh(thacker, n, half_extent);
    let initial = initial_state(thacker, n, half_extent, dx);
    let states = initial.clone();
    let v0 = total_volume(&states, mesh.dx, mesh.dy);

    let t_end = thacker.period() / 2.0;
    let (states, steps) = run_until(states, &mesh, Boundaries2D::WALLS, t_end, 0.4);

    // L2 error on h restricted to cells where the analytical depth
    // exceeds 0.1·h₀ (i.e. the well-interior of the cap, where the
    // numerical scheme is not contaminated by wet/dry).
    let mut l2_num = 0.0;
    let mut l2_den = 0.0;
    let mut linf = 0.0_f64;
    let mut wet_cells = 0usize;
    for ((i, j), s) in states.indexed_iter() {
        let x = -half_extent + (j as f64 + 0.5) * dx;
        let y = -half_extent + (i as f64 + 0.5) * dx;
        let h_an = thacker.depth(x, y, t_end);
        if h_an > 0.1 * thacker.h0 {
            let e = (s.h - h_an).abs();
            l2_num += e * e;
            l2_den += h_an * h_an;
            linf = linf.max(e);
            wet_cells += 1;
        }
    }
    let l2_rel = (l2_num / l2_den).sqrt();
    let v1 = total_volume(&states, mesh.dx, mesh.dy);
    let mass_err = (v1 - v0).abs() / v0;

    eprintln!("\n=== Thacker benchmark report ===");
    eprintln!("Mesh: {n}×{n}, dx = {dx:.4} m");
    eprintln!(
        "h₀ = {} m, a = {} m, B = {} m",
        thacker.h0, thacker.a, thacker.b
    );
    eprintln!(
        "Period T = {:.4} s, integration time t_end = T/2 = {:.4} s, steps = {}",
        thacker.period(),
        t_end,
        steps,
    );
    eprintln!("Interior wet cells (h_an > 0.1·h₀): {}", wet_cells);
    eprintln!("L² relative error on h: {:.4}%", l2_rel * 100.0);
    eprintln!(
        "L∞ error on h: {:.4} m ({:.2}% of h₀)",
        linf,
        100.0 * linf / thacker.h0
    );
    eprintln!("Mass conservation error: {:.2e}", mass_err);
    eprintln!("=================================\n");
}
