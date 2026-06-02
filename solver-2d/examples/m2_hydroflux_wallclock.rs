//! M2d: hydroflux wall-clock benchmark matched to `m2_anuga_wallclock.py`.
//!
//! Same physical problem (Stoker dam-break on a flat-bed channel,
//! 200 m × 5 m, `h_L = 1` m, `t_end = 8` s) at the matched resolution
//! (`Δx = Δy = 0.5` m, giving `400 × 10 = 4000` cells). Wall-clock is
//! measured for the SSP-RK2 time-step loop only and appended to the
//! shared CSV alongside the ANUGA row so the manuscript can read
//! them as a single dataset.
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example m2_hydroflux_wallclock

use std::fs;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use ndarray::Array2;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, cfl_time_step, ssprk2_step,
};

const OUT: &str = "papers/01_review/figures/data/m2_anuga_wallclock.csv";
const LENGTH_X: f64 = 200.0;
const LENGTH_Y: f64 = 5.0;
const DX: f64 = 0.5;
const DY: f64 = 0.5;
const H_L: f64 = 1.0;
const X_DAM: f64 = 100.0;
const T_END: f64 = 8.0;

fn main() {
    fs::create_dir_all(Path::new(OUT).parent().unwrap()).expect("data dir");

    let n_cols = (LENGTH_X / DX) as usize; // 400
    let n_rows = (LENGTH_Y / DY) as usize; // 10
    let n_cells = n_rows * n_cols;
    println!(
        "hydroflux mesh: {n_rows} × {n_cols} = {n_cells} cells, Δx = {DX} m\nPhysical time: {T_END} s"
    );

    let bed = Array2::<f64>::zeros((n_rows, n_cols));
    let mesh = Mesh2D::new(bed, DX, DY, 0.0);
    let mut states = Array2::from_shape_fn((n_rows, n_cols), |(_i, j)| {
        let x = (j as f64 + 0.5) * DX;
        if x < X_DAM {
            Conserved2D::new(H_L, 0.0, 0.0)
        } else {
            Conserved2D::DRY
        }
    });
    let bcs = Boundaries2D {
        north: Boundary::Wall,
        south: Boundary::Wall,
        west: Boundary::Transmissive,
        east: Boundary::Transmissive,
    };

    let mut t = 0.0_f64;
    let mut n_steps = 0_usize;
    let start = Instant::now();
    while t < T_END {
        let dt = cfl_time_step(&states, &mesh, 0.4).min(T_END - t);
        if !dt.is_finite() || dt <= 0.0 {
            break;
        }
        ssprk2_step(&mut states, &mesh, bcs, dt);
        t += dt;
        n_steps += 1;
    }
    let wall = start.elapsed().as_secs_f64();
    let wall_per_sim = wall / T_END;
    let cell_steps = (n_cells * n_steps) as f64;
    let mcs_per_s = cell_steps / wall / 1.0e6;

    println!("Wall clock      : {wall:.3} s");
    println!("SSP-RK2 steps   : {n_steps}");
    println!("Wall / sim-s    : {wall_per_sim:.6}");
    println!("Mcell-steps/s   : {mcs_per_s:.3}");

    // Append a hydroflux row to the shared CSV. If the file does not
    // exist yet (ANUGA not run first), write the header.
    let new_file = !Path::new(OUT).exists();
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(OUT)
        .expect("open CSV");
    if new_file {
        writeln!(
            f,
            "solver,n_cells_or_triangles,n_steps,t_sim_s,wall_clock_s,wall_per_sim_s,mcell_steps_per_s"
        )
        .unwrap();
    }
    writeln!(
        f,
        "hydroflux,{n_cells},{n_steps},{T_END},{wall:.6},{wall_per_sim:.6},{mcs_per_s:.6}"
    )
    .unwrap();
    println!("wrote {OUT}");
}
