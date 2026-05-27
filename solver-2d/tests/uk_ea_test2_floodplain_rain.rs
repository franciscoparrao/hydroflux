//! UK Environment Agency 2D benchmark — equivalent of Test 2
//! "Filling of floodplain depressions" (Néelz & Pender 2013).
//!
//! # Test description
//!
//! A flat (or gently sloping) floodplain receives uniform rainfall
//! over its entire surface. The bed contains several local
//! depressions that must collect water from the surrounding plain.
//! The test exercises:
//!
//! - [`apply_rain`] uniform mass injection over the whole domain.
//! - Bed-slope source driving water laterally into depressions
//!   (wet/dry treatment + Audusse hydrostatic reconstruction).
//! - Steady-state behaviour once rainfall continues at a constant
//!   rate and the system equilibrates.
//!
//! Synthetic stand-in for the official EA Test 2 geometry: 2
//! depressions instead of 3, simple flat bed (the official spec has
//! a mild slope but the essential filling physics is the same).
//!
//! # Setup
//!
//! - Domain: 400 m × 200 m, mesh 100 × 50 (`dx = dy = 4 m`).
//! - Bed: flat (`z = 0`) except 2 quadratic depressions:
//!   - centre (150, 80), radius 25 m, depth 0.3 m
//!   - centre (280, 130), radius 30 m, depth 0.4 m
//! - Initial: dry (`h = 0`).
//! - Rainfall: 50 mm/hour = `1.389e-5 m/s` for `t_end = 3600 s` (1 h).
//! - BC: Wall on all four sides (closed catchment).
//! - Manning `n = 0.06` (rough grassland).
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test uk_ea_test2_floodplain_rain
//! ```

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2D, Mesh2D, apply_rain, cfl_time_step, manning_friction_step,
    ssprk2_step,
};
use ndarray::Array2;

#[derive(Debug, Clone, Copy)]
struct Depression {
    cx: f64,
    cy: f64,
    radius: f64,
    depth: f64,
}

const DEPRESSIONS: &[Depression] = &[
    Depression {
        cx: 150.0,
        cy: 80.0,
        radius: 25.0,
        depth: 0.3,
    },
    Depression {
        cx: 280.0,
        cy: 130.0,
        radius: 30.0,
        depth: 0.4,
    },
];

#[derive(Debug, Clone, Copy)]
struct TestCase {
    length_x: f64,
    length_y: f64,
    /// Rainfall rate [m/s] (50 mm/hour = 1.389e-5 m/s).
    rain_rate: f64,
    manning: f64,
    t_end: f64,
}

impl TestCase {
    fn standard() -> Self {
        Self {
            length_x: 400.0,
            length_y: 200.0,
            rain_rate: 50.0e-3 / 3600.0, // 50 mm / h
            manning: 0.06,
            t_end: 3600.0,
        }
    }
}

fn bed_at(x: f64, y: f64) -> f64 {
    let mut z = 0.0;
    for d in DEPRESSIONS {
        let r2 = (x - d.cx).powi(2) + (y - d.cy).powi(2);
        if r2 < d.radius.powi(2) {
            z -= d.depth * (1.0 - r2 / d.radius.powi(2));
        }
    }
    z
}

fn build_mesh(case: TestCase, n_x: usize, n_y: usize) -> (Mesh2D, f64, f64) {
    let dx = case.length_x / n_x as f64;
    let dy = case.length_y / n_y as f64;
    let bed = Array2::from_shape_fn((n_y, n_x), |(i, j)| {
        let x = (j as f64 + 0.5) * dx;
        let y = (i as f64 + 0.5) * dy;
        bed_at(x, y)
    });
    (Mesh2D::new(bed, dx, dy, case.manning), dx, dy)
}

fn boundaries() -> Boundaries2D {
    Boundaries2D::WALLS
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
        let dt_cfl = cfl_time_step(&states, mesh, cfl);
        // Cap dt by t_end and also by a problem-specific maximum
        // (60 s) so the simulation doesn't take huge steps when
        // the wave speeds are tiny early in the rainfall.
        let dt = dt_cfl.min(60.0).min(case.t_end - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        manning_friction_step(&mut states, mesh, dt, 1.0e-9);
        apply_rain(&mut states, case.rain_rate, dt);
        t += dt;
        steps += 1;
        if steps > 500_000 {
            panic!("UK EA Test 2: {steps} steps without reaching t_end");
        }
    }
    (states, steps)
}

#[test]
fn depth_remains_bounded_and_finite() {
    let case = TestCase::standard();
    let (mesh, _, _) = build_mesh(case, 100, 50);
    let initial = Array2::from_elem((50, 100), Conserved2D::DRY);
    let (final_states, _) = run_until(initial, &mesh, boundaries(), case, 0.4);
    for s in &final_states {
        assert!(s.h.is_finite() && s.h >= 0.0, "h ill-formed: {}", s.h);
        // Maximum possible depth: rate · t + max depression depth
        // ≈ 50 mm + 0.4 m = 0.45 m. Allow 0.6 m as a generous bound.
        assert!(
            s.h < 0.6,
            "h exceeded sensible bound: {} (rain · t = {:.3} m)",
            s.h,
            case.rain_rate * case.t_end
        );
    }
}

#[test]
fn rain_increases_total_mass_at_expected_rate() {
    // Walls all around → no outflow. Total mass at t_end must equal
    // the cumulative rainfall: rate · area · t.
    let case = TestCase::standard();
    let (mesh, _, _) = build_mesh(case, 100, 50);
    let initial = Array2::from_elem((50, 100), Conserved2D::DRY);
    let (final_states, _) = run_until(initial, &mesh, boundaries(), case, 0.4);

    let m_final: f64 = final_states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    let area = case.length_x * case.length_y;
    let cumulative = case.rain_rate * area * case.t_end;
    // Walls + closed catchment: mass conservation should be exact to
    // numerical roundoff over the rain accumulation.
    assert_relative_eq!(m_final, cumulative, epsilon = 1.0e-6 * cumulative);
}

#[test]
fn depressions_accumulate_more_water_than_surroundings() {
    // Cells inside depressions have h_inside = rain_depth + depression_depth
    // (water fills the bowl + the uniform rainfall layer). Cells well
    // away from depressions only carry the rainfall layer. So the
    // average depth in depressions must be substantially higher than
    // the average in non-depression areas.
    let case = TestCase::standard();
    let (mesh, dx, dy) = build_mesh(case, 100, 50);
    let initial = Array2::from_elem((50, 100), Conserved2D::DRY);
    let (final_states, _) = run_until(initial, &mesh, boundaries(), case, 0.4);

    let mut dep_sum = 0.0;
    let mut dep_n = 0_usize;
    let mut plain_sum = 0.0;
    let mut plain_n = 0_usize;
    for ((i, j), s) in final_states.indexed_iter() {
        let x = (j as f64 + 0.5) * dx;
        let y = (i as f64 + 0.5) * dy;
        let inside_dep = DEPRESSIONS.iter().any(|d| {
            let r2 = (x - d.cx).powi(2) + (y - d.cy).powi(2);
            r2 < (d.radius * 0.5).powi(2) // inner third of each depression
        });
        // "Plain" = cells > 60 m from any depression centre.
        let far_from_dep = DEPRESSIONS.iter().all(|d| {
            let r2 = (x - d.cx).powi(2) + (y - d.cy).powi(2);
            r2 > 60.0_f64.powi(2)
        });
        if inside_dep {
            dep_sum += s.h;
            dep_n += 1;
        } else if far_from_dep {
            plain_sum += s.h;
            plain_n += 1;
        }
    }
    assert!(dep_n > 0 && plain_n > 0, "zones empty — geometry issue");
    let dep_avg = dep_sum / dep_n as f64;
    let plain_avg = plain_sum / plain_n as f64;
    assert!(
        dep_avg > plain_avg * 1.5,
        "depressions not noticeably deeper: dep_avg = {:.4} m, plain_avg = {:.4} m",
        dep_avg,
        plain_avg
    );
}

#[test]
fn plain_depth_approximates_uniform_rain_layer() {
    // Far from any depression, the steady state under uniform rain
    // is a uniform layer of depth `rate · t`. The actual cells will
    // deviate slightly because of bed-slope drainage into depressions,
    // but the order of magnitude should match.
    let case = TestCase::standard();
    let (mesh, dx, dy) = build_mesh(case, 100, 50);
    let initial = Array2::from_elem((50, 100), Conserved2D::DRY);
    let (final_states, _) = run_until(initial, &mesh, boundaries(), case, 0.4);

    let expected = case.rain_rate * case.t_end;
    let mut plain_sum = 0.0;
    let mut plain_n = 0_usize;
    for ((i, j), s) in final_states.indexed_iter() {
        let x = (j as f64 + 0.5) * dx;
        let y = (i as f64 + 0.5) * dy;
        let far_from_dep = DEPRESSIONS.iter().all(|d| {
            let r2 = (x - d.cx).powi(2) + (y - d.cy).powi(2);
            r2 > 80.0_f64.powi(2)
        });
        if far_from_dep {
            plain_sum += s.h;
            plain_n += 1;
        }
    }
    let plain_avg = plain_sum / plain_n as f64;
    // Plain depth within a factor of 3 of the analytical (uniform
    // rainfall layer would give EXACTLY `expected`; cell-by-cell
    // drainage into depressions skews this).
    assert!(
        plain_avg > expected * 0.3 && plain_avg < expected * 3.0,
        "plain depth {:.4} m is far from expected rain layer {:.4} m",
        plain_avg,
        expected
    );
}

use approx::assert_relative_eq;
