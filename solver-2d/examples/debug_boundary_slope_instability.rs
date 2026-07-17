//! Reproducer for
//! `docs/bug-report-2026-07-boundary-slope-instability.md`: a thin
//! rain film on a steep ramp blows up under the default SSP-RK2
//! integrator, even though the domain has no boundary-adjacent
//! asymmetry (every column on the ramp carries the same bed jump per
//! cell). Kept as a standalone example (not yet a `#[test]`) pending a
//! decision on the fix — see the bug report for the root-cause
//! write-up and the risk trade-off against the frozen Paper 01 numerics.

use hydroflux_solver_2d::boundary::Boundaries2D;
use hydroflux_solver_2d::geometry::Mesh2D;
use hydroflux_solver_2d::sim::{Simulation, SimulationConfig};
use hydroflux_solver_2d::source::apply_rain;
use hydroflux_solver_2d::state::Conserved2D;
use ndarray::Array2;

fn ramp_mesh(steep: bool) -> Mesh2D {
    let (nr, nc) = (20, 6);
    let slope_per_cell = if steep { 20.0 } else { 0.3 }; // m per 30 m cell
    // Bed is highest at column 0 (West, Transmissive) and descends to
    // zero at the last column (East, also Transmissive).
    let bed = Array2::from_shape_fn((nr, nc), |(_i, j)| (nc - 1 - j) as f64 * slope_per_cell);
    Mesh2D::new(bed, 30.0, 30.0, 0.035)
}

/// Parameterised ramp: `dx` cell size, `dz_per_cell` bed jump per cell
/// (so `dz_per_cell / dx` is the local slope magnitude). Probes the
/// bed-slope source with terrain far steeper than the original 67%
/// reproducer — see
/// docs/bug-report-2026-07-boundary-slope-instability.md §9 (nowcast's
/// Santa Bárbara re-test: 90 m DEM pixels, real precordillera relief)
/// and §10 (the reconstructed-source fix these probes validate).
fn ramp_mesh_scaled(dx: f64, dz_per_cell: f64) -> Mesh2D {
    let (nr, nc) = (20, 6);
    let bed = Array2::from_shape_fn((nr, nc), |(_i, j)| (nc - 1 - j) as f64 * dz_per_cell);
    Mesh2D::new(bed, dx, dx, 0.035)
}

fn run_scaled(label: &str, dx: f64, dz_per_cell: f64, rain_mm_per_day: f64, duration_s: f64, cfl: f64, max_steps: usize) {
    let mesh = ramp_mesh_scaled(dx, dz_per_cell);
    let slope_pct = 100.0 * dz_per_cell / dx;
    let rain_rate_m_s = rain_mm_per_day * 1.0e-3 / 86_400.0;
    let init = Array2::from_elem((mesh.n_rows(), mesh.n_cols()), Conserved2D::DRY);
    let config = SimulationConfig {
        cfl,
        boundaries: Boundaries2D::TRANSMISSIVE,
        ..Default::default()
    };
    let mut sim = Simulation::new(mesh, init, config).unwrap();
    let mut truncated = false;
    while sim.time() < duration_s {
        if sim.steps() >= max_steps {
            truncated = true;
            break;
        }
        let dt = sim.step(duration_s - sim.time()).unwrap();
        apply_rain(sim.states_mut(), rain_rate_m_s, dt);
    }
    let max_depth = sim.states().iter().map(|s| s.h).fold(0.0_f64, f64::max);
    println!(
        "{label}: dx={dx}m slope={slope_pct:.0}% rain={rain_mm_per_day}mm/d dur={duration_s}s steps={} t={:.3} max_depth={:.3} m truncated={}",
        sim.steps(),
        sim.time(),
        max_depth,
        truncated
    );
}

fn run(label: &str, steep: bool, cfl: f64) {
    let mesh = ramp_mesh(steep);
    let rain_rate_m_s = 68.0 * 1.0e-3 / 86_400.0; // 68 mm/day
    let init = Array2::from_elem((mesh.n_rows(), mesh.n_cols()), Conserved2D::DRY);
    let config = SimulationConfig {
        cfl,
        boundaries: Boundaries2D::TRANSMISSIVE,
        ..Default::default()
    };
    let mut sim = Simulation::new(mesh, init, config).unwrap();
    let t_target = 360.0;
    let max_steps = 200_000;
    let mut truncated = false;
    while sim.time() < t_target {
        if sim.steps() >= max_steps {
            truncated = true;
            break;
        }
        let dt = sim.step(t_target - sim.time()).unwrap();
        apply_rain(sim.states_mut(), rain_rate_m_s, dt);
    }
    let max_depth = sim.states().iter().map(|s| s.h).fold(0.0_f64, f64::max);
    println!(
        "{label}: steep={steep} cfl={cfl} steps={} t={:.3} max_depth={:.3} m truncated={}",
        sim.steps(),
        sim.time(),
        max_depth,
        truncated
    );
}

fn run_long_trend(label: &str, steep: bool, cfl: f64, t_target: f64, max_steps: usize) {
    let mesh = ramp_mesh(steep);
    let rain_rate_m_s = 68.0 * 1.0e-3 / 86_400.0;
    let init = Array2::from_elem((mesh.n_rows(), mesh.n_cols()), Conserved2D::DRY);
    let config = SimulationConfig {
        cfl,
        boundaries: Boundaries2D::TRANSMISSIVE,
        ..Default::default()
    };
    let mut sim = Simulation::new(mesh, init, config).unwrap();
    let mut truncated = false;
    let mut next_report = 360.0;
    while sim.time() < t_target {
        if sim.steps() >= max_steps {
            truncated = true;
            break;
        }
        let dt = sim.step(t_target - sim.time()).unwrap();
        apply_rain(sim.states_mut(), rain_rate_m_s, dt);
        if sim.time() >= next_report {
            let max_depth = sim.states().iter().map(|s| s.h).fold(0.0_f64, f64::max);
            println!(
                "  {label}: steps={:7} t={:9.1} max_depth={:12.3} m",
                sim.steps(),
                sim.time(),
                max_depth
            );
            next_report += 360.0;
        }
    }
    let max_depth = sim.states().iter().map(|s| s.h).fold(0.0_f64, f64::max);
    println!(
        "{label} FINAL: steep={steep} cfl={cfl} steps={} t={:.3} max_depth={:.3} m truncated={}",
        sim.steps(),
        sim.time(),
        max_depth,
        truncated
    );
}

fn main() {
    run("A", false, 0.4); // control: gentle slope (~1%)
    run("B", true, 0.4); // steep slope (~70%), normal CFL
    run("C", true, 0.1); // same slope, MORE conservative CFL
    run("D", true, 0.05); // even more conservative CFL

    println!("\n--- long-duration trend probe (case B, 100x the original window) ---");
    run_long_trend("B-long", true, 0.4, 36_000.0, 2_000_000);

    println!("\n--- Santa Barbara-style stress probe (90 m pixels, extreme local relief) ---");
    // Mirrors nowcast's exact reproducer parameters (rain 32.6 mm/day,
    // duration 900 s, cfl 0.4, max-steps 20000) at increasing slope
    // magnitude, up to a deeply pathological near-vertical ramp.
    run_scaled("SB-100pct", 90.0, 90.0, 32.6, 900.0, 0.4, 20_000);
    run_scaled("SB-200pct", 90.0, 180.0, 32.6, 900.0, 0.4, 20_000);
    run_scaled("SB-300pct", 90.0, 270.0, 32.6, 900.0, 0.4, 20_000);
    run_scaled("SB-1000pct", 90.0, 900.0, 32.6, 900.0, 0.4, 20_000);
}
