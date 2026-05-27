//! UK Environment Agency 2D benchmark — equivalent of Test 1
//! "Flooding a disconnected water body" (Néelz & Pender 2013).
//!
//! # Test description
//!
//! A point source upstream releases water onto a sloping plane. The
//! plane contains a topographic depression off the natural flow line;
//! the question is whether the solver routes water down the slope
//! **and** fills the depression once the flow path overtops the
//! depression's rim.
//!
//! What the test exercises:
//!
//! - [`apply_point_sources`] sustained over many timesteps.
//! - Bed-slope source under non-trivial topography (sloping plane +
//!   localised depression).
//! - Wet/dry propagation along the flow path.
//! - "Filling a hole": the depression is dry initially, gets wet via
//!   overflow once the channel reaches it.
//!
//! Synthetic stand-in for the official EA Test 1 geometry (which has
//! three depressions); the essential physics — point inflow + slope
//! + filling — are preserved with one depression for clarity.
//!
//! # Setup
//!
//! - Domain: 200 m × 100 m, mesh 100 × 50 (`dx = dy = 2 m`).
//! - Bed: longitudinal slope `S₀ = 0.005` descending in `+x`, plus a
//!   circular depression of radius 15 m centred at (100, 50) with
//!   bottom 0.5 m below the local slope baseline.
//! - Initial state: thin film `h = 1 mm` (workaround for the
//!   Discharge-on-dry limitation, even though we don't use
//!   Discharge here — the thin film also helps the FV update start
//!   stable when a point source first fires).
//! - Point source at cell (25, 5) — coordinates ≈ (10, 50) — with
//!   `Q = 1 m³/s` sustained.
//! - BC: Wall on `W`, `N`, `S`; Transmissive on `E`.
//! - Manning `n = 0.035`.
//! - `t_end = 400 s` (6:40 — long enough for the channel to reach
//!   the depression and overflow into it).
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test uk_ea_test1_disconnected_pond
//! ```

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, PointSource, apply_point_sources, cfl_time_step,
    manning_friction_step, ssprk2_step,
};
use ndarray::Array2;

#[derive(Debug, Clone, Copy)]
struct TestCase {
    length_x: f64,
    length_y: f64,
    slope: f64,
    /// Depression centre (x, y) [m].
    depression_centre: (f64, f64),
    /// Depression radius [m].
    depression_radius: f64,
    /// Depression depth below the local slope baseline [m].
    depression_depth: f64,
    /// Point source location (x, y) [m] and rate [m³/s].
    source_x: f64,
    source_y: f64,
    source_q: f64,
    manning: f64,
    t_end: f64,
}

impl TestCase {
    fn standard() -> Self {
        Self {
            length_x: 200.0,
            length_y: 100.0,
            slope: 0.005,
            depression_centre: (100.0, 50.0),
            depression_radius: 15.0,
            depression_depth: 0.5,
            source_x: 10.0,
            source_y: 50.0,
            source_q: 1.0,
            manning: 0.035,
            t_end: 400.0,
        }
    }
}

fn build_mesh(case: TestCase, n_x: usize, n_y: usize) -> (Mesh2D, f64, f64) {
    let dx = case.length_x / n_x as f64;
    let dy = case.length_y / n_y as f64;
    let bed = Array2::from_shape_fn((n_y, n_x), |(i, j)| {
        let x = (j as f64 + 0.5) * dx;
        let y = (i as f64 + 0.5) * dy;
        let z_slope = -case.slope * x;
        let (cx, cy) = case.depression_centre;
        let r2 = (x - cx).powi(2) + (y - cy).powi(2);
        // Smooth depression: a quadratic bowl that reaches its
        // minimum at the centre. Avoids a discontinuous step that
        // could degrade well-balancedness.
        let z_depression = if r2 < case.depression_radius.powi(2) {
            -case.depression_depth * (1.0 - r2 / case.depression_radius.powi(2))
        } else {
            0.0
        };
        z_slope + z_depression
    });
    (Mesh2D::new(bed, dx, dy, case.manning), dx, dy)
}

fn boundaries() -> Boundaries2D {
    Boundaries2D {
        west: Boundary::Wall,
        east: Boundary::Transmissive,
        north: Boundary::Wall,
        south: Boundary::Wall,
    }
}

fn source_cell(case: TestCase, dx: f64, dy: f64) -> PointSource {
    PointSource {
        row: (case.source_y / dy) as usize,
        col: (case.source_x / dx) as usize,
        q_mass: case.source_q,
    }
}

fn run_until(
    mut states: Array2<Conserved2D>,
    mesh: &Mesh2D,
    bcs: Boundaries2D,
    sources: &[PointSource],
    case: TestCase,
    cfl: f64,
) -> (Array2<Conserved2D>, usize) {
    let mut t = 0.0;
    let mut steps = 0;
    while t < case.t_end {
        let dt = cfl_time_step(&states, mesh, cfl).min(case.t_end - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        manning_friction_step(&mut states, mesh, dt, 1.0e-9);
        apply_point_sources(&mut states, sources, dt, mesh.dx, mesh.dy);
        t += dt;
        steps += 1;
        if steps > 500_000 {
            panic!("UK EA Test 1: {steps} steps without reaching t_end");
        }
    }
    (states, steps)
}

#[test]
fn depth_remains_bounded_and_finite() {
    let case = TestCase::standard();
    let (mesh, dx, dy) = build_mesh(case, 100, 50);
    let initial = Array2::from_elem((50, 100), Conserved2D::new(0.001, 0.0, 0.0));
    let source = source_cell(case, dx, dy);
    let (final_states, _) = run_until(initial, &mesh, boundaries(), &[source], case, 0.4);
    for s in &final_states {
        assert!(s.h.is_finite() && s.h >= 0.0, "h ill-formed: {}", s.h);
        assert!(s.hu.is_finite() && s.hv.is_finite(), "momentum non-finite");
    }
}

#[test]
fn point_source_injects_mass_into_domain() {
    // After `t_end` with sustained Q = 1 m³/s, the cumulative
    // injected mass is Q · t_end = 400 m³. The domain must hold a
    // non-trivial fraction of this (some leaves through the East
    // transmissive boundary once the channel reaches it).
    let case = TestCase::standard();
    let (mesh, dx, dy) = build_mesh(case, 100, 50);
    let initial = Array2::from_elem((50, 100), Conserved2D::new(0.001, 0.0, 0.0));
    let source = source_cell(case, dx, dy);
    let m_initial: f64 = initial.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    let (final_states, _) = run_until(initial.clone(), &mesh, boundaries(), &[source], case, 0.4);
    let m_final: f64 = final_states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    let cumulative_inflow = case.source_q * case.t_end;
    // Mass in the domain at t_end must be strictly greater than the
    // initial (sustained injection), and below cumulative_inflow +
    // initial (outflow only loses mass).
    assert!(
        m_final > m_initial,
        "mass did not grow: m_initial = {m_initial}, m_final = {m_final}"
    );
    assert!(
        m_final <= m_initial + cumulative_inflow * 1.001,
        "mass exceeded cumulative inflow bound: m_final = {m_final}, bound = {}",
        m_initial + cumulative_inflow * 1.001
    );
    let _ = dx;
    let _ = dy;
}

#[test]
fn depression_fills_via_overflow() {
    // The depression centre is initially below the slope baseline.
    // After sustained inflow + slope-driven flow, the depression
    // must contain water deeper than the depression rim height
    // would imply on a flat bed (i.e. it actually fills, not just
    // gets wet at the bottom).
    let case = TestCase::standard();
    let (mesh, dx, dy) = build_mesh(case, 100, 50);
    let initial = Array2::from_elem((50, 100), Conserved2D::new(0.001, 0.0, 0.0));
    let source = source_cell(case, dx, dy);
    let (final_states, _) = run_until(initial, &mesh, boundaries(), &[source], case, 0.4);

    let (cx, cy) = case.depression_centre;
    let i_centre = (cy / dy) as usize;
    let j_centre = (cx / dx) as usize;
    let h_centre = final_states[(i_centre, j_centre)].h;

    // The depression centre must have depth ≥ half the depression
    // depth (after sustained inflow it should be ≥ the full depth,
    // but allow margin for transient drain).
    assert!(
        h_centre > case.depression_depth * 0.5,
        "depression centre did not fill: h = {:.4}, expected ≥ {:.4}",
        h_centre,
        case.depression_depth * 0.5
    );
}

#[test]
fn off_channel_cells_stay_drier_than_channel() {
    // The flow channel runs along y ≈ 50 (the source row). Cells
    // far from the channel (e.g. y < 20 or y > 80, where the bed
    // slope alone shouldn't direct water) must stay substantially
    // shallower than the channel itself.
    let case = TestCase::standard();
    let (mesh, dx, dy) = build_mesh(case, 100, 50);
    let initial = Array2::from_elem((50, 100), Conserved2D::new(0.001, 0.0, 0.0));
    let source = source_cell(case, dx, dy);
    let (final_states, _) = run_until(initial, &mesh, boundaries(), &[source], case, 0.4);

    let i_channel = (case.source_y / dy) as usize;
    // Sample the channel at three downstream x-positions (outside
    // the depression).
    let xs_channel = [30.0_f64, 60.0_f64, 150.0_f64];
    let h_channel: Vec<f64> = xs_channel
        .iter()
        .map(|x| {
            let j = (x / dx) as usize;
            final_states[(i_channel, j)].h
        })
        .collect();
    let h_channel_avg: f64 = h_channel.iter().sum::<f64>() / h_channel.len() as f64;

    // Off-channel rows.
    let i_far_south = (10.0 / dy) as usize;
    let i_far_north = (90.0 / dy) as usize;
    let mut h_off_max: f64 = 0.0;
    for j in 0..final_states.ncols() {
        h_off_max = h_off_max.max(final_states[(i_far_south, j)].h);
        h_off_max = h_off_max.max(final_states[(i_far_north, j)].h);
    }

    assert!(
        h_off_max < h_channel_avg,
        "off-channel cells wetter than channel: h_off_max = {:.4}, h_channel_avg = {:.4}",
        h_off_max,
        h_channel_avg
    );
    let _ = dx;
    let _ = dy;
}

#[test]
#[ignore = "informational: prints depth snapshot at key positions"]
fn report_depth_snapshot() {
    let case = TestCase::standard();
    let (mesh, dx, dy) = build_mesh(case, 200, 100);
    let initial = Array2::from_elem((100, 200), Conserved2D::new(0.001, 0.0, 0.0));
    let source = source_cell(case, dx, dy);
    let (final_states, steps) = run_until(initial, &mesh, boundaries(), &[source], case, 0.4);

    eprintln!("\n=== UK EA Test 1 (synthetic): disconnected pond ===");
    eprintln!(
        "Mesh: 200×100, dx = {dx:.2} m, t_end = {} s, steps = {}",
        case.t_end, steps
    );
    eprintln!(
        "Slope = {}, Manning n = {}, Q_source = {} m³/s",
        case.slope, case.manning, case.source_q
    );
    eprintln!(
        "Depression: centre ({:.0}, {:.0}), radius {:.0} m, depth {:.2} m",
        case.depression_centre.0,
        case.depression_centre.1,
        case.depression_radius,
        case.depression_depth
    );

    let i_channel = (case.source_y / dy) as usize;
    eprintln!(
        "\nChannel depth profile (y = {} m) at t_end:",
        case.source_y
    );
    eprintln!("{:>6} {:>10}", "x [m]", "h [m]");
    for j in (0..final_states.ncols()).step_by(20) {
        let x = (j as f64 + 0.5) * dx;
        let h = final_states[(i_channel, j)].h;
        eprintln!("{:>6.1} {:>10.4}", x, h);
    }

    let m_final: f64 = final_states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();
    let cumulative_inflow = case.source_q * case.t_end;
    eprintln!(
        "\nMass balance: final = {:.1} m³, cumulative inflow = {:.1} m³, ratio = {:.3}",
        m_final,
        cumulative_inflow,
        m_final / cumulative_inflow
    );
    eprintln!("===================================================\n");
}
