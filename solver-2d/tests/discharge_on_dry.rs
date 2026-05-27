//! Discharge-BC-on-dry inflow over a sloping channel.
//!
//! Before 2026-05-23, [`Boundary::Discharge`] on a fully dry boundary
//! cell did nothing: the ghost depth equalled `inner.h = 0`, so the
//! HLLC face saw a dry-dry state and returned zero flux. The
//! prescribed `q` never entered the domain unless the caller
//! pre-initialised a thin film (see, e.g., UK EA Test 5).
//!
//! The fix (commit [next]) sets `ghost.h` to local Manning normal
//! depth when the inner cell is dry AND the streamwise bed slope is
//! above the safety threshold. This test exercises that path:
//!
//! - Fully dry initial state (`h ≡ 0`).
//! - Sloping channel (`S₀ = 0.005`, comfortably above the 5e-4 threshold).
//! - Discharge BC on the West side, transmissive on the East, walls
//!   on N/S.
//!
//! Expected behaviour:
//! 1. After enough time, the interior is wet — the inflow actually
//!    propagated.
//! 2. The interior depth at the gauge converges toward Manning
//!    normal depth (within Lax-Friedrichs-like wet-front tolerance).
//! 3. Mass balance is bounded by cumulative inflow.
//!
//! Reproducir:
//! ```text
//! cargo test --release -p hydroflux-solver-2d --test discharge_on_dry
//! ```

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, cfl_time_step_with_bcs, manning_friction_step,
    ssprk2_step,
};
use ndarray::Array2;

const N_X: usize = 80;
const N_Y: usize = 20;
const DX: f64 = 2.0;
const DY: f64 = 2.0;
const SLOPE: f64 = 0.005;
const MANNING: f64 = 0.04;
const Q_IN: f64 = 1.5; // m²/s
const T_END: f64 = 300.0;
const CFL: f64 = 0.4;

fn sloped_bed() -> Array2<f64> {
    Array2::from_shape_fn((N_Y, N_X), |(_i, j)| -SLOPE * (j as f64 + 0.5) * DX)
}

fn build_mesh() -> Mesh2D {
    Mesh2D::new(sloped_bed(), DX, DY, MANNING)
}

fn boundaries() -> Boundaries2D {
    Boundaries2D {
        west: Boundary::Discharge { q: Q_IN },
        east: Boundary::Transmissive,
        north: Boundary::Wall,
        south: Boundary::Wall,
    }
}

fn run_to_t_end(mesh: &Mesh2D, mut states: Array2<Conserved2D>) -> Array2<Conserved2D> {
    let bcs = boundaries();
    let mut t = 0.0;
    let mut steps = 0;
    while t < T_END {
        // `cfl_time_step_with_bcs` folds ghost-cell wave speeds into the
        // CFL bound, so the first step from a fully dry interior sees
        // the Manning-normal-depth ghost on the West face and picks a
        // sane dt instead of `INFINITY`.
        let dt = cfl_time_step_with_bcs(&states, mesh, bcs, CFL).min(T_END - t);
        ssprk2_step(&mut states, mesh, bcs, dt);
        manning_friction_step(&mut states, mesh, dt, 1.0e-9);
        t += dt;
        steps += 1;
        if steps > 500_000 {
            panic!("discharge_on_dry: {steps} steps");
        }
    }
    states
}

fn manning_normal_depth() -> f64 {
    (MANNING * Q_IN.abs() / SLOPE.sqrt()).powf(3.0 / 5.0)
}

#[test]
fn inflow_actually_enters_domain_from_fully_dry_state() {
    // Without the Discharge-on-dry fix, this test fails: every cell
    // stays at h = 0 forever because the HLLC face sees dry-dry on
    // the boundary and emits zero flux.
    let mesh = build_mesh();
    let initial = Array2::from_elem((N_Y, N_X), Conserved2D::DRY);
    let final_states = run_to_t_end(&mesh, initial);

    // Sanity: no NaN, no negative depth, no runaway.
    for s in &final_states {
        assert!(s.h.is_finite(), "h non-finite");
        assert!(s.h >= 0.0, "h negative: {}", s.h);
        assert!(s.h < 5.0, "h exceeded sanity bound: {}", s.h);
    }

    // After T_END = 300 s of sustained inflow the wetting front
    // should have moved well past the midpoint (gravity-wave celerity
    // on Manning normal depth ≈ √(g·h_n) ≈ 3 m/s; front would travel
    // ~900 m, the domain is 160 m). Both the inflow column and the
    // midpoint should carry positive depth.
    let mid = N_X / 2;
    let mid_row = N_Y / 2;
    let h_inflow = final_states[(mid_row, 0)].h;
    let h_mid = final_states[(mid_row, mid)].h;
    let h_n = manning_normal_depth();
    assert!(
        h_inflow > 0.5 * h_n,
        "inflow column dry: h = {h_inflow:.4} (expected ≥ {:.4})",
        0.5 * h_n
    );
    assert!(
        h_mid > 0.5 * h_n,
        "midpoint dry: h = {h_mid:.4} (expected ≥ {:.4})",
        0.5 * h_n
    );
}

#[test]
fn interior_depth_approaches_manning_normal_depth() {
    // Tolerance 30 %: the cold-start dam-break-on-dry transient
    // deposits more mass at the leading edge than the asymptotic
    // Manning normal depth, and that excess takes several domain
    // traversals to drain through the transmissive East outlet.
    // Tightening this requires either a much longer `T_END` or
    // initialising at `h_n` (defeats the point of the test).
    let mesh = build_mesh();
    let initial = Array2::from_elem((N_Y, N_X), Conserved2D::DRY);
    let final_states = run_to_t_end(&mesh, initial);

    let mid = N_X / 2;
    let mid_row = N_Y / 2;
    let h_mid = final_states[(mid_row, mid)].h;
    let h_n = manning_normal_depth();
    let ratio = h_mid / h_n;
    assert!(
        (ratio - 1.0).abs() < 0.30,
        "h_mid/h_n = {ratio:.3} (expected within 30 %); h_mid = {h_mid:.4}, h_n = {h_n:.4}"
    );
}

#[test]
fn mass_balance_is_bounded_by_cumulative_inflow() {
    // Walls on N/S, transmissive E. Cumulative inflow over the West
    // face equals Q_IN · L_y · T_END. Final mass must be at most
    // that (some has exited East) and at least a positive fraction
    // (inflow is steady, domain hasn't drained yet).
    let mesh = build_mesh();
    let initial = Array2::from_elem((N_Y, N_X), Conserved2D::DRY);
    let m_initial: f64 = initial.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();

    let final_states = run_to_t_end(&mesh, initial);
    let m_final: f64 = final_states.iter().map(|s| s.h * mesh.dx * mesh.dy).sum();

    let cumulative_inflow = Q_IN * (N_Y as f64 * DY) * T_END;
    assert_eq!(m_initial, 0.0);
    assert!(m_final > 0.0, "final mass is zero; inflow never engaged");
    assert!(
        m_final < cumulative_inflow * 1.05,
        "mass exceeded cumulative inflow: m_final = {m_final:.3e}, inflow = {cumulative_inflow:.3e}"
    );
}
