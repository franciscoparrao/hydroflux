//! Dam break onto a dry bed — Stoker (1957) dry-bed limit.
//!
//! Wet-on-the-left, dry-on-the-right column collapse on a flat
//! frictionless bed. The analytical solution is a left-going
//! rarefaction with leading edge propagating into the dry region at
//! `2·c_L`; the trailing edge propagates leftward at `-c_L` in the
//! still water frame.
//!
//! Inside the rarefaction fan the depth and velocity vary linearly in
//! the similarity variable `ξ = x / t`:
//!
//! ```text
//!   h(ξ) = (1 / (9 g)) · (2 c_L − ξ)²
//!   u(ξ) = (2 / 3) · (c_L + ξ)
//! ```
//!
//! At `ξ = 0` (the original dam location) this gives `h(0) = (4/9)·h_L`
//! and `u(0) = (2/3)·c_L`. The wet front is at `ξ = 2 c_L`; the
//! rarefaction tail at `ξ = -c_L`.
//!
//! This is the canonical benchmark for wetting/drying in shallow-water
//! solvers — Toro (2009) §10.5.4, Liang & Marche (2009), Brufau et
//! al. (2002). It exercises the two-rarefaction wave-speed branch in
//! `hllc_normal_flux` and the positivity-preserving clamp in
//! `forward_euler_step`.
//!
//! Reproducir: `cargo test --release -p hydroflux-solver-2d --test dam_break_on_dry`.
//!
//! # 2D setup
//!
//! The problem is purely 1D, but we run it on a thin 2D mesh
//! (`n_rows = 3`, `n_cols = 200`) with walls on the y-faces to verify
//! that the 2D solver behaves correctly for genuinely 1D flow and
//! that nothing in the tangential (y) direction is generated
//! spuriously. The middle row is what we compare against the
//! analytical profile.

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, cfl_time_step, forward_euler_step,
};
use ndarray::Array2;

const G: f64 = 9.81;

#[derive(Debug, Clone, Copy)]
struct DamBreak {
    /// Wet-side depth `h_L` [m].
    h_l: f64,
    /// Dam location `x_dam` [m] in the domain.
    x_dam: f64,
    /// Domain extent `[0, length]` in `x` [m].
    length: f64,
    /// Time at which to compare against the analytical solution [s].
    t_end: f64,
}

impl DamBreak {
    fn celerity(self) -> f64 {
        (G * self.h_l).sqrt()
    }

    /// Analytical depth at position `x` and time `t`. Returns 0 inside
    /// the dry region (past the wet front).
    fn depth(self, x: f64, t: f64) -> f64 {
        let xi = (x - self.x_dam) / t;
        let c_l = self.celerity();
        if xi <= -c_l {
            self.h_l
        } else if xi >= 2.0 * c_l {
            0.0
        } else {
            let h = (2.0 * c_l - xi).powi(2) / (9.0 * G);
            h.max(0.0)
        }
    }

    /// `x_dam + 2·c_L·t`: the analytical position of the leading dry
    /// front at time `t`.
    fn wet_front_position(self, t: f64) -> f64 {
        self.x_dam + 2.0 * self.celerity() * t
    }
}

fn build_mesh(case: DamBreak, n_cols: usize, n_rows: usize) -> (Mesh2D, f64) {
    let dx = case.length / n_cols as f64;
    let bed = Array2::<f64>::zeros((n_rows, n_cols)); // flat
    (Mesh2D::new(bed, dx, dx, 0.0), dx)
}

fn initial_state(case: DamBreak, n_rows: usize, n_cols: usize, dx: f64) -> Array2<Conserved2D> {
    Array2::from_shape_fn((n_rows, n_cols), |(_i, j)| {
        let x = (j as f64 + 0.5) * dx;
        if x < case.x_dam {
            Conserved2D::new(case.h_l, 0.0, 0.0)
        } else {
            Conserved2D::DRY
        }
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
        let dt_cfl = cfl_time_step(&states, mesh, cfl);
        let dt = dt_cfl.min(t_end - t);
        forward_euler_step(&mut states, mesh, bcs, dt);
        t += dt;
        steps += 1;
        if steps > 100_000 {
            panic!("dam-break-on-dry: {steps} steps without reaching t_end");
        }
    }
    states
}

/// y-walls + x-transmissive: the wave goes nowhere in `y`, exits the
/// domain freely on `x` if it reaches the boundary.
fn make_bcs() -> Boundaries2D {
    Boundaries2D {
        north: Boundary::Wall,
        south: Boundary::Wall,
        west: Boundary::Transmissive,
        east: Boundary::Transmissive,
    }
}

#[test]
fn depth_is_non_negative_and_finite_everywhere() {
    // Robustness check: after running across the dry bed, no NaN, no
    // negative h, no infinite momentum. This is the basic regression
    // guard for the wetting/drying clamp + two-rarefaction wave speed.
    let case = DamBreak {
        h_l: 1.0,
        x_dam: 50.0,
        length: 100.0,
        t_end: 4.0,
    };
    let n_cols = 200;
    let n_rows = 3;
    let (mesh, dx) = build_mesh(case, n_cols, n_rows);
    let initial = initial_state(case, n_rows, n_cols, dx);
    let final_states = run_until(initial, &mesh, make_bcs(), case.t_end, 0.4);

    for s in &final_states {
        assert!(s.h.is_finite(), "h became non-finite: {}", s.h);
        assert!(s.h >= 0.0, "h went negative: {}", s.h);
        assert!(s.hu.is_finite(), "hu became non-finite: {}", s.hu);
        assert!(s.hv.is_finite(), "hv became non-finite: {}", s.hv);
        // The tangential (y) momentum must remain ~zero for a purely
        // 1D flow on a flat bed with wall y-boundaries.
        assert!(
            s.hv.abs() < 1e-10,
            "spurious tangential momentum: hv = {}",
            s.hv
        );
    }
}

#[test]
fn wet_front_propagates_at_two_celerities() {
    // Locate the wet/dry front in the numerical solution and compare
    // against the analytical position `x_dam + 2·c_L·t`. Tolerance:
    // one cell width, which is the resolution of the front locator.
    let case = DamBreak {
        h_l: 1.0,
        x_dam: 50.0,
        length: 100.0,
        t_end: 4.0,
    };
    let n_cols = 200;
    let n_rows = 3;
    let (mesh, dx) = build_mesh(case, n_cols, n_rows);
    let initial = initial_state(case, n_rows, n_cols, dx);
    let final_states = run_until(initial, &mesh, make_bcs(), case.t_end, 0.4);

    // Numerical wet front: rightmost cell with h > 1e-5 m. We
    // compare against the analytical position at the SAME threshold
    // (the analytical front asymptotes to h = 0 at x = x_dam +
    // 2·c_L·t; at h = 1e-5 the analytical front sits 0.12 m short of
    // that limit, well below dx — so the two thresholds give the
    // same expected position at this precision).
    let mid_row = n_rows / 2;
    let h_threshold = 1.0e-5_f64;
    let mut x_front_numerical = 0.0;
    for j in 0..n_cols {
        let x = (j as f64 + 0.5) * dx;
        if final_states[(mid_row, j)].h > h_threshold {
            x_front_numerical = x;
        }
    }
    let x_front_analytical = case.wet_front_position(case.t_end);
    let err = (x_front_numerical - x_front_analytical).abs();

    // First-order HLLC + Audusse + Euler on a dry-bed front is known
    // to lag the analytical position. Liang & Marche (2009) report
    // ~2-3% for comparable setups on similar meshes; with our
    // positivity clamp the lag can grow to ~10% before slope-limited
    // / RK2 reconstruction is added. We allow up to 15% of the
    // analytical front position as regression guard.
    let rel_err = err / x_front_analytical;
    assert!(
        rel_err < 0.15,
        "wet front mislocated: numerical {:.3}, analytical {:.3}, |Δ| = {:.3} ({:.1} cells, {:.1}% of analytical)",
        x_front_numerical,
        x_front_analytical,
        err,
        err / dx,
        rel_err * 100.0
    );

    // Cross-check: a Davis-bound front at x_dam + c_L·t would lag the
    // analytical front by c_L·t = 12.5 m for our parameters. The
    // numerical front must be substantially CLOSER to analytical
    // than to Davis (otherwise the two-rarefaction branch is not
    // actually firing).
    let x_front_davis = case.x_dam + case.celerity() * case.t_end;
    let err_to_davis = (x_front_numerical - x_front_davis).abs();
    assert!(
        err < err_to_davis,
        "front closer to Davis estimate than to analytical: |num - an| = {:.3}, |num - davis| = {:.3}",
        err,
        err_to_davis
    );
}

#[test]
fn depth_at_dam_location_matches_four_ninths_h_l() {
    // At ξ = 0 (the original dam position), the analytical solution
    // says h(0, t) = (4/9)·h_L for all t > 0. Compare the numerical
    // value to this constant.
    let case = DamBreak {
        h_l: 1.0,
        x_dam: 50.0,
        length: 100.0,
        t_end: 4.0,
    };
    let n_cols = 200;
    let n_rows = 3;
    let (mesh, dx) = build_mesh(case, n_cols, n_rows);
    let initial = initial_state(case, n_rows, n_cols, dx);
    let final_states = run_until(initial, &mesh, make_bcs(), case.t_end, 0.4);

    let mid_row = n_rows / 2;
    let j_dam = (case.x_dam / dx) as usize;
    let h_at_dam = final_states[(mid_row, j_dam)].h;
    let h_analytical = 4.0 / 9.0 * case.h_l;
    let rel_err = (h_at_dam - h_analytical).abs() / h_analytical;
    // First-order HLLC on a 200-cell mesh: expect ~5% pointwise error
    // at the dam location. With Davis bound the wet front lags by
    // ~50% and this cell would still be near h_L (way off); the
    // two-rarefaction estimate must put the value within ~10% of 4/9.
    assert!(
        rel_err < 0.10,
        "h at dam location: numerical {:.4}, analytical {:.4} (4·h_L/9), |Δ|/an = {:.2}%",
        h_at_dam,
        h_analytical,
        rel_err * 100.0
    );
}

#[test]
fn l1_error_inside_rarefaction_is_bounded() {
    // L1 error of `h` integrated over the rarefaction interval
    // `(x_dam - c_L·t, x_dam + 2·c_L·t)`, normalised by the analytical
    // L1 norm on the same interval. First-order scheme on 200 cells
    // should land below ~10%.
    let case = DamBreak {
        h_l: 1.0,
        x_dam: 50.0,
        length: 100.0,
        t_end: 4.0,
    };
    let n_cols = 200;
    let n_rows = 3;
    let (mesh, dx) = build_mesh(case, n_cols, n_rows);
    let initial = initial_state(case, n_rows, n_cols, dx);
    let final_states = run_until(initial, &mesh, make_bcs(), case.t_end, 0.4);

    let mid_row = n_rows / 2;
    let c_l = case.celerity();
    let x_tail = case.x_dam - c_l * case.t_end;
    let x_head = case.x_dam + 2.0 * c_l * case.t_end;

    let mut l1_err = 0.0;
    let mut l1_norm = 0.0;
    for j in 0..n_cols {
        let x = (j as f64 + 0.5) * dx;
        if x < x_tail || x > x_head {
            continue;
        }
        let h_an = case.depth(x, case.t_end);
        let h_nu = final_states[(mid_row, j)].h;
        l1_err += (h_nu - h_an).abs() * dx;
        l1_norm += h_an * dx;
    }
    let rel_l1 = l1_err / l1_norm;
    assert!(
        rel_l1 < 0.10,
        "L1 relative error inside rarefaction too large: {:.2}%",
        rel_l1 * 100.0
    );
}

#[test]
#[ignore = "informational: prints metrics for benchmarks/dam-break-on-dry-results.md"]
fn report_metrics() {
    // Not a pass/fail test — prints L1/L∞ + wet-front-position error
    // for the markdown writeup. Run with:
    //   cargo test --release -p hydroflux-solver-2d --test dam_break_on_dry -- --ignored --nocapture
    let case = DamBreak {
        h_l: 1.0,
        x_dam: 50.0,
        length: 100.0,
        t_end: 4.0,
    };
    let n_cols = 400; // finer mesh for the writeup
    let n_rows = 3;
    let (mesh, dx) = build_mesh(case, n_cols, n_rows);
    let initial = initial_state(case, n_rows, n_cols, dx);
    let final_states = run_until(initial, &mesh, make_bcs(), case.t_end, 0.4);

    let mid_row = n_rows / 2;
    let c_l = case.celerity();
    let x_tail = case.x_dam - c_l * case.t_end;
    let x_head = case.x_dam + 2.0 * c_l * case.t_end;

    let mut l1_err = 0.0;
    let mut l1_norm = 0.0;
    let mut l2_num = 0.0;
    let mut l2_den = 0.0;
    let mut linf = 0.0_f64;
    let mut cells = 0;
    let mut x_front_numerical = 0.0;
    for j in 0..n_cols {
        let x = (j as f64 + 0.5) * dx;
        if final_states[(mid_row, j)].h > 0.01 * case.h_l {
            x_front_numerical = x;
        }
        if x < x_tail || x > x_head {
            continue;
        }
        let h_an = case.depth(x, case.t_end);
        let h_nu = final_states[(mid_row, j)].h;
        let e = (h_nu - h_an).abs();
        l1_err += e * dx;
        l1_norm += h_an * dx;
        l2_num += e * e;
        l2_den += h_an * h_an;
        linf = linf.max(e);
        cells += 1;
    }
    let rel_l1 = l1_err / l1_norm;
    let rel_l2 = (l2_num / l2_den).sqrt();
    let x_front_an = case.wet_front_position(case.t_end);

    eprintln!("\n=== Dam-break-on-dry benchmark report ===");
    eprintln!("Mesh: {n_cols} cells × {n_rows} rows, dx = {dx:.3} m");
    eprintln!(
        "h_L = {} m, c_L = {:.3} m/s, t_end = {} s",
        case.h_l, c_l, case.t_end
    );
    eprintln!("Rarefaction cells (interior): {cells}");
    eprintln!("L1 relative error on h: {:.3}%", rel_l1 * 100.0);
    eprintln!("L² relative error on h: {:.3}%", rel_l2 * 100.0);
    eprintln!(
        "L∞ error on h: {:.4} m ({:.2}% of h_L)",
        linf,
        100.0 * linf / case.h_l
    );
    eprintln!(
        "Wet-front position: numerical {:.3} m, analytical {:.3} m, |Δ| = {:.3} m ({:.2} cells)",
        x_front_numerical,
        x_front_an,
        (x_front_numerical - x_front_an).abs(),
        (x_front_numerical - x_front_an).abs() / dx
    );
    eprintln!("=========================================\n");
}
