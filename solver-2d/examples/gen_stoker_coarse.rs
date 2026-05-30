//! Hydroflux Stoker dam-break at a coarse resolution matched to the
//! ANUGA reference run (`anuga_stoker_compare.py`) for the head-to-head
//! comparison in §3.8 of the methods paper. Mesh 100×3, dx = 1 m.
//!
//! Output: papers/01_review/figures/data/verif_stoker_coarse.csv
//!   columns: x, h_analytical, h_sim
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example gen_stoker_coarse

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use ndarray::Array2;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, cfl_time_step, ssprk2_step,
};

const G: f64 = 9.81;
const OUT: &str = "papers/01_review/figures/data/verif_stoker_coarse.csv";

fn stoker_depth(x: f64, t: f64, h_l: f64, x_dam: f64) -> f64 {
    let c_l = (G * h_l).sqrt();
    let xi = (x - x_dam) / t;
    if xi <= -c_l {
        h_l
    } else if xi >= 2.0 * c_l {
        0.0
    } else {
        ((2.0 * c_l - xi).powi(2) / (9.0 * G)).max(0.0)
    }
}

fn main() {
    fs::create_dir_all(Path::new(OUT).parent().unwrap()).expect("data dir");

    let (h_l, x_dam, length, t_end) = (1.0_f64, 50.0_f64, 100.0_f64, 4.0_f64);
    let (n_cols, n_rows) = (100usize, 3usize);
    let dx = length / n_cols as f64;
    let bed = Array2::<f64>::zeros((n_rows, n_cols));
    let mesh = Mesh2D::new(bed, dx, dx, 0.0);
    let mut states = Array2::from_shape_fn((n_rows, n_cols), |(_i, j)| {
        let x = (j as f64 + 0.5) * dx;
        if x < x_dam {
            Conserved2D::new(h_l, 0.0, 0.0)
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
    let mut t = 0.0;
    while t < t_end {
        let dt = cfl_time_step(&states, &mesh, 0.4).min(t_end - t);
        ssprk2_step(&mut states, &mesh, bcs, dt);
        t += dt;
    }

    let mid = n_rows / 2;
    let mut f = File::create(OUT).expect("create");
    writeln!(f, "x,h_analytical,h_sim").unwrap();
    for j in 0..n_cols {
        let x = (j as f64 + 0.5) * dx;
        let h_an = stoker_depth(x, t_end, h_l, x_dam);
        let h_sim = states[(mid, j)].h;
        writeln!(f, "{x:.4},{h_an:.6},{h_sim:.6}").unwrap();
    }
    println!("wrote {OUT}");
}
