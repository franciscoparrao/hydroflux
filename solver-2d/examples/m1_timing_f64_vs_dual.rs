//! M1c: wall-clock comparison between `T = f64` and `T = Dual<f64>`
//! for the 2D solver.
//!
//! Forward-mode automatic differentiation carries a `(value, dval)`
//! pair through every operation. The expected overhead is constant
//! per arithmetic operation — between 2× and 3× the f64 cost,
//! depending on how many of the per-step FLOPs cancel under
//! optimisation (the dval-side of derivatives sometimes ends up
//! multiplied by zero constants).
//!
//! This bench runs a fixed dam-break problem for a fixed number of
//! SSP-RK2 steps and reports wall-clock seconds for the two backends
//! plus the ratio, on the same problem size, after a warm-up. The
//! number that lands in the manuscript is `t_dual / t_f64`.
//!
//! Output:
//!   papers/01_review/figures/data/m1_timing.csv
//!     columns: backend, n_cells, n_steps, wall_clock_s, per_cell_step_ns
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example m1_timing_f64_vs_dual

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use hydroflux_autograd::{Dual, Real};
use ndarray::Array2;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2DG, Mesh2DG, manning_friction_step, ssprk2_step,
};

const OUT_CSV: &str = "papers/01_review/figures/data/m1_timing.csv";

const N_ROWS: usize = 64;
const N_COLS: usize = 64;
const DX: f64 = 1.0;
const DY: f64 = 1.0;
const N_STEPS: usize = 200;
const DT_FIXED: f64 = 0.01; // CFL-stable for the chosen problem; keeps
                             // both runs taking the same number of steps,
                             // so the comparison is fair without coupling
                             // to per-step CFL-time-step calls.

fn build_state<T: Real>() -> (Array2<Conserved2DG<T>>, Mesh2DG<T>) {
    let bed = Array2::<T>::from_elem((N_ROWS, N_COLS), T::zero());
    let mesh = Mesh2DG::<T>::new(bed, DX, DY, T::from_f64(0.03));
    // Gaussian bump initial condition — gives a non-trivial wave that
    // exercises the well-balanced bed-slope source, MUSCL slopes, and
    // the Audusse reconstruction at every face.
    let cx = (N_COLS as f64 - 1.0) / 2.0;
    let cy = (N_ROWS as f64 - 1.0) / 2.0;
    let w_sq = (N_ROWS as f64).powi(2) / 50.0;
    let states = Array2::from_shape_fn((N_ROWS, N_COLS), |(i, j)| {
        let dxx = j as f64 - cx;
        let dyy = i as f64 - cy;
        let h = 1.0 + 0.5 * (-(dxx * dxx + dyy * dyy) / w_sq).exp();
        Conserved2DG::<T>::new_generic(T::from_f64(h), T::zero(), T::zero())
    });
    (states, mesh)
}

fn time_run<T: Real>(label: &str, n_seed: T) -> (f64, f64) {
    let (mut states, mut mesh) = build_state::<T>();
    // Seed the Manning n if requested — for the Dual run this turns
    // every cell into a derivative-carrying value during the friction
    // step (the dominant arithmetic in the Dual cost).
    let manning_seeded = Array2::<T>::from_elem((N_ROWS, N_COLS), n_seed);
    mesh = Mesh2DG::<T>::with_manning_field(
        Array2::<T>::from_elem((N_ROWS, N_COLS), T::zero()),
        DX,
        DY,
        manning_seeded,
    );

    let bcs = Boundaries2D::TRANSMISSIVE;
    let start = Instant::now();
    for _ in 0..N_STEPS {
        ssprk2_step(&mut states, &mesh, bcs, DT_FIXED);
        manning_friction_step(&mut states, &mesh, DT_FIXED, 1.0e-6);
    }
    let elapsed = start.elapsed().as_secs_f64();
    let per_cell_step_ns =
        elapsed * 1.0e9 / (N_ROWS as f64 * N_COLS as f64 * N_STEPS as f64);
    println!(
        "  {label:>14}  elapsed = {elapsed:>7.3} s   per cell-step = {per_cell_step_ns:>7.1} ns"
    );
    (elapsed, per_cell_step_ns)
}

fn main() {
    fs::create_dir_all(Path::new(OUT_CSV).parent().unwrap()).expect("create data dir");

    let cells = N_ROWS * N_COLS;
    println!(
        "Bench: {N_ROWS} × {N_COLS} = {cells} cells, {N_STEPS} SSP-RK2 + friction steps\n"
    );

    // Warm-up — JIT/caches don't apply here (Rust is AOT), but the
    // OS page-fault path and allocator hot-cache do. One warm-up run
    // per backend keeps the timing stable run-to-run.
    println!("Warm-up:");
    let _ = time_run::<f64>("f64 (warm)", 0.03);
    let _ = time_run::<Dual>("Dual (warm)", Dual::variable(0.03));

    println!("\nMeasurement:");
    let (t_f64, pcs_f64) = time_run::<f64>("f64", 0.03);
    let (t_dual, pcs_dual) = time_run::<Dual>("Dual<f64>", Dual::variable(0.03));
    let ratio = t_dual / t_f64;
    println!(
        "\n  ratio (Dual / f64) = {ratio:.2}×  (forward-mode AD overhead)"
    );

    let mut f = File::create(OUT_CSV).expect("create CSV");
    writeln!(f, "backend,n_cells,n_steps,wall_clock_s,per_cell_step_ns").unwrap();
    writeln!(f, "f64,{cells},{N_STEPS},{t_f64:.6},{pcs_f64:.3}").unwrap();
    writeln!(
        f,
        "Dual<f64>,{cells},{N_STEPS},{t_dual:.6},{pcs_dual:.3}"
    )
    .unwrap();
    println!("\n  wrote {OUT_CSV}");
}
