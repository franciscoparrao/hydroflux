//! M2a: serial performance on large synthetic grids.
//!
//! Runs the full SSP-RK2 + Manning friction pipeline on 256², 512²,
//! and 1024² grids with a Gaussian-bump initial condition that
//! exercises every face (no dry interior) so the cell-mask early-skip
//! does not artificially flatter the timing. Reports cell-step
//! throughput in ns per cell-step and Mcell-steps per second; the
//! latter is the unit other 2D shallow-water solver papers (Rak 2024,
//! Saleem & Norman 2024, Caviedes-Voullième circle) use to publish
//! their headline numbers.
//!
//! Output:
//!   papers/01_review/figures/data/m2_perf_large.csv
//!     columns: n_cells, n_steps, wall_clock_s, ns_per_cell_step,
//!              mcell_steps_per_s
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example m2_perf_large_grid

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use ndarray::Array2;

use hydroflux_solver_2d::{
    Boundaries2D, Conserved2D, Mesh2D, manning_friction_step, ssprk2_step,
};

const OUT_CSV: &str = "papers/01_review/figures/data/m2_perf_large.csv";

const SIZES: &[usize] = &[256, 512, 1024];
const N_STEPS: usize = 50;
const DT_FIXED: f64 = 0.01; // CFL-stable for the Gaussian-bump regime

fn build_state(n: usize) -> (Array2<Conserved2D>, Mesh2D) {
    let bed = Array2::<f64>::zeros((n, n));
    let mesh = Mesh2D::new(bed, 1.0, 1.0, 0.03);
    let cx = (n as f64 - 1.0) / 2.0;
    let cy = (n as f64 - 1.0) / 2.0;
    let w_sq = (n as f64).powi(2) / 50.0;
    let states = Array2::from_shape_fn((n, n), |(i, j)| {
        let dxx = j as f64 - cx;
        let dyy = i as f64 - cy;
        let h = 1.0 + 0.5 * (-(dxx * dxx + dyy * dyy) / w_sq).exp();
        Conserved2D::new(h, 0.0, 0.0)
    });
    (states, mesh)
}

fn run_bench(n: usize) -> (f64, f64, f64) {
    let (mut states, mesh) = build_state(n);
    let bcs = Boundaries2D::TRANSMISSIVE;

    // Warm-up to populate the allocator's hot cache + page-fault any
    // freshly-mapped memory.
    for _ in 0..3 {
        ssprk2_step(&mut states, &mesh, bcs, DT_FIXED);
        manning_friction_step(&mut states, &mesh, DT_FIXED, 1.0e-6);
    }

    // Measure.
    let start = Instant::now();
    for _ in 0..N_STEPS {
        ssprk2_step(&mut states, &mesh, bcs, DT_FIXED);
        manning_friction_step(&mut states, &mesh, DT_FIXED, 1.0e-6);
    }
    let elapsed = start.elapsed().as_secs_f64();
    let cells = (n * n) as f64;
    let cell_steps = cells * N_STEPS as f64;
    let ns_per_cs = elapsed * 1.0e9 / cell_steps;
    let mcs_per_s = cell_steps / elapsed / 1.0e6;
    (elapsed, ns_per_cs, mcs_per_s)
}

fn main() {
    fs::create_dir_all(Path::new(OUT_CSV).parent().unwrap()).expect("create data dir");

    println!(
        "Serial perf bench: {N_STEPS} SSP-RK2 + Manning friction steps per grid\n"
    );
    println!(
        "  {:>5}   {:>11}   {:>9}   {:>14}   {:>17}",
        "n", "cells", "wall (s)", "ns/cell-step", "Mcell-steps/s"
    );

    let mut rows: Vec<(usize, f64, f64, f64)> = Vec::new();
    for &n in SIZES {
        let (elapsed, ns_per_cs, mcs_per_s) = run_bench(n);
        println!(
            "  {:>5}   {:>11}   {:>9.3}   {:>14.1}   {:>17.2}",
            n,
            n * n,
            elapsed,
            ns_per_cs,
            mcs_per_s
        );
        rows.push((n, elapsed, ns_per_cs, mcs_per_s));
    }

    let mut f = File::create(OUT_CSV).expect("create CSV");
    writeln!(f, "n,n_cells,n_steps,wall_clock_s,ns_per_cell_step,mcell_steps_per_s").unwrap();
    for (n, e, ns, mcs) in &rows {
        let cells = n * n;
        writeln!(
            f,
            "{n},{cells},{N_STEPS},{e:.6},{ns:.3},{mcs:.6}"
        )
        .unwrap();
    }
    println!("\nwrote {OUT_CSV}");
}
