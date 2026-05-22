//! UK Environment Agency 2D benchmark — equivalent of Test 4
//! "Speed of propagation of a flood wave" (Néelz & Pender 2013).
//!
//! # Test description
//!
//! The original UK EA Test 4 measures how fast a flood wave released
//! from upstream propagates over an initially-dry, sloping floodplain,
//! and how the wave deforms under Manning friction. The official
//! geometry is a 1000 m × 2000 m floodplain with a 0.001 longitudinal
//! slope, channel banks raised by 0.5 m, and a hydrograph released
//! from one upstream cell. Gauges measure depth at a sequence of
//! downstream points.
//!
//! Our implementation here is a **synthetic stand-in** that captures
//! the essential physics without the exact UK EA geometry file:
//!
//! - Domain: 2000 m × 500 m, mesh 200 × 50 cells (`dx = dy = 10` m).
//! - Bed: longitudinal slope `S₀ = 0.001` in `+x`, raised banks on
//!   the `N` and `S` edges (`bank_height = 2.5 m`, `bank_strip =
//!   30 m`) to confine the flow.
//! - Manning `n = 0.04` (vegetated floodplain).
//! - Initial state: thin film `h = 1 mm` (workaround for the
//!   Discharge-BC-on-dry limitation, see the note below).
//! - BC: Discharge `q = 2 m²/s` at `W`, Transmissive at `E`, Wall on
//!   `N` and `S`. With these parameters the Manning normal depth in
//!   the central channel is `h_n = (n·q/√S₀)^(3/5) ≈ 1.75 m`, well
//!   below the 2.5 m bank, so the steady state confines flow.
//! - `t_end = 1500 s` (25 min). Wave celerity `c = √(g·h_n) ≈ 4.1
//!   m/s` crosses the 2 km reach in ~ 490 s, so `t_end` gives 3×
//!   the traversal time for the steady state to develop.
//!
//! Substituting the actual EA Test 4 geometry (when the official
//! ASCII grid is downloaded) is a drop-in replacement: only the
//! `build_mesh` function changes.
//!
//! # Limitations surfaced by this test
//!
//! 1. `Boundary::Discharge` with a fully-dry inner cell injects zero
//!    flux — the HLLC sees a dry-dry interface (ghost.h = inner.h =
//!    0 by zero-gradient) and returns 0. Workaround: thin-film init
//!    (1 mm). A critical-depth ghost override was attempted but
//!    interacts badly with raised-bed terrain (places water above
//!    the bank elevation, producing unphysical head differences and
//!    runaway depths). A robust fix requires deriving `h_ghost` from
//!    the local bed slope (Manning normal depth) AND respecting the
//!    inner cell's `η`. Deferred.
//!
//! 2. `Boundary::Discharge` applies `q` UNIFORMLY to every cell on
//!    the boundary face. The raised banks therefore receive some
//!    inflow, just less than the channel. The
//!    `channel_carries_more_water_than_banks` test asserts the
//!    relative magnitude rather than perfect bank confinement.
//!
//! # What it exercises
//!
//! - Discharge BC sustained over many timesteps.
//! - Wet/dry front propagation across a dry floodplain
//!   (two-rarefaction wave speeds, flux rescaling, positivity).
//! - Manning friction over a kilometre-scale reach.
//! - Mass conservation under a sustained inflow with a transmissive
//!   outlet.
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test uk_ea_test4_propagation
//! ```

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, cfl_time_step, manning_friction_step, ssprk2_step,
};
use ndarray::Array2;

const G: f64 = 9.81;

#[derive(Debug, Clone, Copy)]
struct TestCase {
    /// Domain length in `x` [m].
    length_x: f64,
    /// Domain width in `y` [m].
    length_y: f64,
    /// Longitudinal bed slope `S₀` (positive descends in `+x`).
    slope: f64,
    /// Bank height [m] on the `N` and `S` edges (raise the bed by
    /// this amount in a strip of width `bank_strip`).
    bank_height: f64,
    /// Width of the raised-bank strip on each `y` edge [m].
    bank_strip: f64,
    /// Manning roughness `n` [s/m^(1/3)].
    manning: f64,
    /// Prescribed upstream unit discharge `q = hu` [m²/s].
    q_in: f64,
    /// Final simulation time [s].
    t_end: f64,
}

impl TestCase {
    fn standard() -> Self {
        // Manning normal depth in the central channel:
        //   h_n = (n · q / √S₀)^(3/5)
        //       = (0.04 · 2 / √0.001)^(3/5)
        //       = (2.53)^(3/5) ≈ 1.75 m
        // Bank height 2.5 m comfortably contains this. The wave at
        // c = √(g · h_n) ≈ 4.1 m/s crosses 2000 m in ~ 490 s; t_end
        // 1500 s gives 3× the traversal time for the steady state
        // to develop.
        Self {
            length_x: 2000.0,
            length_y: 500.0,
            slope: 0.001,
            bank_height: 2.5,
            bank_strip: 30.0,
            manning: 0.04,
            q_in: 2.0,
            t_end: 1500.0,
        }
    }
}

fn build_mesh(case: TestCase, n_x: usize, n_y: usize) -> (Mesh2D, f64) {
    let dx = case.length_x / n_x as f64;
    let dy = case.length_y / n_y as f64;
    let bed = Array2::from_shape_fn((n_y, n_x), |(i, j)| {
        let x = (j as f64 + 0.5) * dx;
        let y = (i as f64 + 0.5) * dy;
        // Longitudinal slope: bed descends in +x.
        let z_main = -case.slope * x;
        // Banks: raise the bed near y = 0 and y = length_y to
        // confine the flow to the central channel.
        let near_south = y < case.bank_strip;
        let near_north = y > case.length_y - case.bank_strip;
        let z_bank = if near_south || near_north {
            case.bank_height
        } else {
            0.0
        };
        z_main + z_bank
    });
    (Mesh2D::new(bed, dx, dy, case.manning), dx)
}

fn boundaries(case: TestCase) -> Boundaries2D {
    Boundaries2D {
        west: Boundary::Discharge { q: case.q_in },
        east: Boundary::Transmissive,
        north: Boundary::Wall,
        south: Boundary::Wall,
    }
}

/// Run the simulation, recording the depth at a set of gauge cells at
/// every output time. Returns the final state plus the gauge time
/// series.
fn run_with_gauges(
    case: TestCase,
    mesh: &Mesh2D,
    bcs: Boundaries2D,
    gauges: &[(usize, usize)],
    cfl: f64,
) -> (Array2<Conserved2D>, Vec<(f64, Vec<f64>)>) {
    let n_rows = mesh.n_rows();
    let n_cols = mesh.n_cols();
    // Thin-film initialisation (1 mm). Needed because the current
    // Discharge BC injects zero flux when inner is fully dry — see
    // module docstring.
    let mut states = Array2::from_elem((n_rows, n_cols), Conserved2D::new(0.001, 0.0, 0.0));
    let mut t = 0.0;
    let mut steps = 0;
    let mut series: Vec<(f64, Vec<f64>)> = Vec::new();
    series.push((0.0, gauges.iter().map(|&(i, j)| states[(i, j)].h).collect()));
    while t < case.t_end {
        let dt = cfl_time_step(&states, mesh, cfl).min(case.t_end - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        manning_friction_step(&mut states, case.manning, dt, 1.0e-9);
        t += dt;
        steps += 1;
        series.push((t, gauges.iter().map(|&(i, j)| states[(i, j)].h).collect()));
        if steps > 200_000 {
            panic!("UK EA Test 4: {steps} steps without reaching t_end");
        }
    }
    (states, series)
}

fn gauge_locations(case: TestCase, mesh: &Mesh2D) -> Vec<(usize, usize)> {
    let n_cols = mesh.n_cols();
    let mid_row = mesh.n_rows() / 2;
    // Gauges at 200, 500, 1000, 1500, 1900 m along the centreline.
    let xs = [200.0, 500.0, 1000.0, 1500.0, 1900.0];
    xs.iter()
        .map(|x| {
            let j = ((x / case.length_x) * n_cols as f64) as usize;
            (mid_row, j.min(n_cols - 1))
        })
        .collect()
}

fn total_mass(states: &Array2<Conserved2D>, dx: f64, dy: f64) -> f64 {
    states.iter().map(|s| s.h * dx * dy).sum()
}

#[test]
fn depth_remains_bounded_and_finite() {
    // Robustness check. With a sustained Discharge inflow into a
    // sloping dry floodplain we expect a wet wave to fill the domain.
    // No NaN, no negative depth.
    let case = TestCase::standard();
    let n_x = 200;
    let n_y = 50;
    let (mesh, _dx) = build_mesh(case, n_x, n_y);
    let gauges = gauge_locations(case, &mesh);
    let (final_states, _series) = run_with_gauges(case, &mesh, boundaries(case), &gauges, 0.4);
    for s in &final_states {
        assert!(s.h.is_finite(), "h non-finite: {}", s.h);
        assert!(s.h >= 0.0, "h negative: {}", s.h);
        assert!(s.hu.is_finite(), "hu non-finite: {}", s.hu);
        assert!(s.hv.is_finite(), "hv non-finite: {}", s.hv);
    }
}

#[test]
fn wave_arrives_at_progressively_later_times_downstream() {
    // The most fundamental property of a propagating flood wave: it
    // reaches downstream gauges later than upstream ones. Detect
    // "arrival" as the first time a gauge's depth exceeds 1 cm.
    let case = TestCase::standard();
    let n_x = 200;
    let n_y = 50;
    let (mesh, _dx) = build_mesh(case, n_x, n_y);
    let gauges = gauge_locations(case, &mesh);
    let (_, series) = run_with_gauges(case, &mesh, boundaries(case), &gauges, 0.4);

    let n_gauges = gauges.len();
    let mut arrivals: Vec<Option<f64>> = vec![None; n_gauges];
    let threshold = 0.01_f64;
    for (t, depths) in &series {
        for k in 0..n_gauges {
            if arrivals[k].is_none() && depths[k] > threshold {
                arrivals[k] = Some(*t);
            }
        }
    }

    // All gauges must have been reached by the final time.
    for (k, arrival) in arrivals.iter().enumerate() {
        assert!(
            arrival.is_some(),
            "gauge {k} never reached the {} cm depth threshold",
            threshold * 100.0
        );
    }

    // Arrivals must be strictly monotone in the downstream direction.
    let times: Vec<f64> = arrivals.iter().map(|a| a.unwrap()).collect();
    for k in 1..n_gauges {
        assert!(
            times[k] > times[k - 1],
            "non-monotone arrival times: gauge {k} arrived at {:.2} s, gauge {} at {:.2} s",
            times[k],
            k - 1,
            times[k - 1]
        );
    }
}

#[test]
fn mass_balance_consistent_with_inflow_and_outflow() {
    // For a sustained Discharge inflow and transmissive outflow, the
    // mass balance is:
    //   M(t) ≈ Q · width · t − (mass that left through E + clamp loss)
    // We verify that the mass IN the domain at t_end is positive and
    // bounded above by the cumulative inflow (mass can leave, not be
    // created).
    let case = TestCase::standard();
    let n_x = 200;
    let n_y = 50;
    let (mesh, _dx) = build_mesh(case, n_x, n_y);
    let gauges = gauge_locations(case, &mesh);
    let (final_states, _) = run_with_gauges(case, &mesh, boundaries(case), &gauges, 0.4);

    let m_final = total_mass(&final_states, mesh.dx, mesh.dy);
    // Effective width of the inflow boundary: the central channel
    // (length_y − 2 · bank_strip). The Discharge BC is applied to
    // every cell on the W face; cells whose ghost h is below the
    // bank-raised bed do not actually pass flow.
    let effective_width = case.length_y - 2.0 * case.bank_strip;
    let cumulative_inflow = case.q_in * effective_width * case.t_end;

    assert!(m_final > 0.0, "no water in the domain: {m_final}");
    assert!(
        m_final <= cumulative_inflow * 1.2,
        "domain mass {} exceeds bound {} (cumulative inflow × 1.2)",
        m_final,
        cumulative_inflow * 1.2
    );
}

#[test]
fn channel_carries_more_water_than_banks() {
    // The central channel must develop substantially more depth than
    // the raised-bank cells: the bed-slope source pushes water to
    // the lower channel and Manning friction works harder on the
    // shallower bank flow. The current `Boundary::Discharge` applies
    // `q` uniformly across all cells on the W face (no per-cell BC
    // yet), so the banks DO receive some inflow — but the steady
    // state should still favour the channel by a wide margin.
    //
    // Concretely: average depth in cells deep in the bank (within
    // bank_strip/3 of the wall) must be at most HALF the average
    // depth in the central channel (well away from banks).
    let case = TestCase::standard();
    let n_x = 200;
    let n_y = 50;
    let (mesh, _dx) = build_mesh(case, n_x, n_y);
    let dy = mesh.dy;
    let gauges = gauge_locations(case, &mesh);
    let (final_states, _) = run_with_gauges(case, &mesh, boundaries(case), &gauges, 0.4);

    let inner_bank_threshold = case.bank_strip / 3.0;
    let mid_channel_zone = case.bank_strip + 50.0;

    let mut bank_sum = 0.0;
    let mut bank_n = 0_usize;
    let mut channel_sum = 0.0;
    let mut channel_n = 0_usize;
    for ((i, _j), s) in final_states.indexed_iter() {
        let y = (i as f64 + 0.5) * dy;
        if y < inner_bank_threshold || y > case.length_y - inner_bank_threshold {
            bank_sum += s.h;
            bank_n += 1;
        } else if y > mid_channel_zone && y < case.length_y - mid_channel_zone {
            channel_sum += s.h;
            channel_n += 1;
        }
    }
    assert!(bank_n > 0 && channel_n > 0, "zones empty — geometry issue");
    let bank_avg = bank_sum / bank_n as f64;
    let channel_avg = channel_sum / channel_n as f64;
    assert!(
        bank_avg < 0.5 * channel_avg,
        "banks not noticeably drier than channel: bank avg h = {:.3} m, channel avg h = {:.3} m",
        bank_avg,
        channel_avg
    );
}

#[test]
#[ignore = "informational: prints gauge time series + summary"]
fn report_gauge_time_series() {
    // Not a pass/fail test — prints depth at each gauge at logarithmic
    // time intervals so the gauge hydrograph can be inspected.
    // Run with:
    //   cargo test --release -p hydroflux-solver-2d --test uk_ea_test4_propagation \
    //       -- --ignored --nocapture
    let case = TestCase::standard();
    let n_x = 400;
    let n_y = 50;
    let (mesh, dx) = build_mesh(case, n_x, n_y);
    let gauges = gauge_locations(case, &mesh);
    let (final_states, series) = run_with_gauges(case, &mesh, boundaries(case), &gauges, 0.4);

    eprintln!("\n=== UK EA Test 4 (synthetic): flood wave propagation ===");
    eprintln!(
        "Mesh: {n_x}×{n_y}, dx = {dx:.2} m, t_end = {} s",
        case.t_end
    );
    eprintln!(
        "Slope = {}, Manning n = {}, q_in = {} m²/s, bank height = {} m",
        case.slope, case.manning, case.q_in, case.bank_height
    );
    let c_max = (G * case.q_in.powf(2.0 / 3.0)).sqrt();
    eprintln!("Characteristic celerity ~ {:.2} m/s", c_max);

    let times_to_print = [10.0, 60.0, 300.0, 600.0, 900.0, 1200.0, case.t_end];
    eprintln!(
        "\n{:>8} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "t [s]", "h@200m", "h@500m", "h@1000m", "h@1500m", "h@1900m"
    );
    for &t_target in &times_to_print {
        if let Some((t, depths)) = series.iter().find(|(t, _)| *t >= t_target) {
            eprintln!(
                "{:>8.1} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
                t, depths[0], depths[1], depths[2], depths[3], depths[4]
            );
        }
    }
    let m_final = total_mass(&final_states, mesh.dx, mesh.dy);
    let effective_width = case.length_y - 2.0 * case.bank_strip;
    let cumulative_inflow = case.q_in * effective_width * case.t_end;
    eprintln!(
        "\nMass balance: final = {:.1} m³, cumulative inflow ≈ {:.1} m³, ratio = {:.3}",
        m_final,
        cumulative_inflow,
        m_final / cumulative_inflow
    );
    eprintln!("=========================================================\n");
}
