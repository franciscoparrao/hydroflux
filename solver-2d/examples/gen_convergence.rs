//! Mesh-refinement convergence study on the Thacker oscillating
//! paraboloid — the canonical 2D analytical benchmark with a moving
//! wet/dry shoreline (cf. Liang & Marche 2009). Runs the solver at a
//! sequence of resolutions, measures the relative L1 and L2 error in
//! depth `h` at t = T/2 against the analytical solution, and reports
//! the observed order of accuracy between consecutive grids.
//!
//! Output: papers/01_review/figures/data/convergence_thacker.csv
//!   columns: n, dx, rel_L1, rel_L2, steps, walltime_s
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example gen_convergence

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use ndarray::Array2;

use hydroflux_solver_2d::{Boundaries2D, Conserved2D, Mesh2D, cfl_time_step, ssprk2_step};

const G: f64 = 9.81;
const OUT: &str = "papers/01_review/figures/data/convergence_thacker.csv";

// Thacker parameters (same as the verification test).
const H0: f64 = 0.1;
const A: f64 = 1.0;
const B: f64 = 0.1;
const HALF_EXTENT: f64 = 1.25;
const CFL: f64 = 0.4;

fn thacker_depth(x: f64, y: f64, t: f64) -> f64 {
    let omega = (2.0 * G * H0).sqrt() / A;
    let xc = B * (omega * t).cos();
    let yc = B * (omega * t).sin();
    let r2 = (x - xc).powi(2) + (y - yc).powi(2);
    (H0 / (A * A) * (A * A - r2)).max(0.0)
}

/// Run Thacker on an n×n grid to t = T/2, return (rel_L1, rel_L2, steps).
fn run_resolution(n: usize) -> (f64, f64, usize) {
    let dx = 2.0 * HALF_EXTENT / n as f64;
    let omega = (2.0 * G * H0).sqrt() / A;
    let t_end = std::f64::consts::PI / omega; // T/2

    let bed = Array2::from_shape_fn((n, n), |(i, j)| {
        let x = -HALF_EXTENT + (j as f64 + 0.5) * dx;
        let y = -HALF_EXTENT + (i as f64 + 0.5) * dx;
        H0 * ((x * x + y * y) / (A * A) - 1.0)
    });
    let mesh = Mesh2D::new(bed, dx, dx, 0.0);

    // Initial condition at t = 0 (cap centred at (B, 0), velocity (0, Bω)).
    let v0 = B * omega;
    let mut states = Array2::from_shape_fn((n, n), |(i, j)| {
        let x = -HALF_EXTENT + (j as f64 + 0.5) * dx;
        let y = -HALF_EXTENT + (i as f64 + 0.5) * dx;
        let h = thacker_depth(x, y, 0.0);
        Conserved2D::new(h, 0.0, h * v0)
    });

    let mut t = 0.0;
    let mut steps = 0usize;
    while t < t_end {
        let dt = cfl_time_step(&states, &mesh, CFL).min(t_end - t);
        ssprk2_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        t += dt;
        steps += 1;
    }

    // Relative L1 and L2 over the whole grid (dry cells contribute 0 to
    // both numerator and denominator). Cell area dx² cancels in the
    // relative norms but is kept for clarity.
    let area = dx * dx;
    let mut l1_err = 0.0;
    let mut l1_norm = 0.0;
    let mut l2_err = 0.0;
    let mut l2_norm = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = -HALF_EXTENT + (j as f64 + 0.5) * dx;
            let y = -HALF_EXTENT + (i as f64 + 0.5) * dx;
            let h_an = thacker_depth(x, y, t_end);
            let e = (states[(i, j)].h - h_an).abs();
            l1_err += e * area;
            l1_norm += h_an * area;
            l2_err += e * e * area;
            l2_norm += h_an * h_an * area;
        }
    }
    ((l1_err / l1_norm), (l2_err / l2_norm).sqrt(), steps)
}

fn main() {
    fs::create_dir_all(Path::new(OUT).parent().unwrap()).expect("data dir");
    let resolutions = [32usize, 64, 128, 256];

    let mut rows: Vec<(usize, f64, f64, f64, usize, f64)> = Vec::new();
    println!(
        "{:>5} {:>9} {:>12} {:>12} {:>7} {:>8}",
        "n", "dx [m]", "rel_L1", "rel_L2", "steps", "t [s]"
    );
    for &n in &resolutions {
        let t0 = Instant::now();
        let (l1, l2, steps) = run_resolution(n);
        let wt = t0.elapsed().as_secs_f64();
        let dx = 2.0 * HALF_EXTENT / n as f64;
        println!("{n:>5} {dx:>9.5} {l1:>12.5e} {l2:>12.5e} {steps:>7} {wt:>8.1}");
        rows.push((n, dx, l1, l2, steps, wt));
    }

    // Observed order p between consecutive grids: p = log(e_c/e_f)/log(2)
    // for a 2× refinement.
    println!("\nObserved order (consecutive 2× refinements):");
    println!("{:>12} {:>10} {:>10}", "pair", "order_L1", "order_L2");
    for w in rows.windows(2) {
        let (nc, _, l1c, l2c, _, _) = w[0];
        let (nf, _, l1f, l2f, _, _) = w[1];
        let p_l1 = (l1c / l1f).log2();
        let p_l2 = (l2c / l2f).log2();
        println!("{:>5}→{:<6} {p_l1:>10.3} {p_l2:>10.3}", nc, nf);
    }

    let mut f = File::create(OUT).expect("create csv");
    writeln!(f, "n,dx,rel_L1,rel_L2,steps,walltime_s").unwrap();
    for (n, dx, l1, l2, steps, wt) in &rows {
        writeln!(f, "{n},{dx:.6},{l1:.6e},{l2:.6e},{steps},{wt:.2}").unwrap();
    }
    println!("\nWrote {OUT}");
}
