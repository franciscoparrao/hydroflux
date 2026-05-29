//! Generate verification-figure data: runs the Stoker dam-break,
//! MacDonald uniform-flow, and Thacker oscillating benchmarks and
//! writes (position, analytical, simulated) profiles as CSV for the
//! methods-paper verification panel (Figure 2).
//!
//! Output CSVs go to `papers/01_review/figures/data/`:
//!   - `verif_stoker.csv`    (x, h_analytical, h_sim) at t = 4 s
//!   - `verif_macdonald.csv` (x, h_analytical, h_sim) at steady state
//!   - `verif_thacker.csv`   (x, h_analytical, h_sim) on the centre row
//!
//! The analytical solutions and setups mirror the integration tests
//! (`dam_break_on_dry.rs`, `macdonald_uniform.rs`, `thacker.rs`); this
//! example exists only to export the profiles that the asserting tests
//! consume internally.
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example gen_verification_data

use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use ndarray::Array2;

use hydroflux_solver_2d::{
    Boundaries2D, Boundary, Conserved2D, Mesh2D, cfl_time_step, manning_friction_step, ssprk2_step,
};

const G: f64 = 9.81;
const OUT_DIR: &str = "papers/01_review/figures/data";

// --- Stoker dam-break on dry bed -----------------------------------

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

fn gen_stoker(out: &Path) -> std::io::Result<()> {
    let (h_l, x_dam, length, t_end) = (1.0_f64, 50.0_f64, 100.0_f64, 4.0_f64);
    let (n_cols, n_rows) = (400usize, 3usize);
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
    let mut f = File::create(out)?;
    writeln!(f, "x,h_analytical,h_sim")?;
    for j in 0..n_cols {
        let x = (j as f64 + 0.5) * dx;
        let h_an = stoker_depth(x, t_end, h_l, x_dam);
        let h_sim = states[(mid, j)].h;
        writeln!(f, "{x:.4},{h_an:.6},{h_sim:.6}")?;
    }
    Ok(())
}

// --- MacDonald uniform flow (well-balanced steady state) -----------

fn gen_macdonald(out: &Path) -> std::io::Result<()> {
    let (q, slope, manning) = (1.0_f64, 0.01_f64, 0.03_f64);
    let h_n = (manning * q / slope.sqrt()).powf(3.0 / 5.0);
    let (n_rows, n_cols, dx) = (5usize, 50usize, 1.0_f64);
    let bed = Array2::from_shape_fn((n_rows, n_cols), |(_i, j)| -(j as f64) * dx * slope);
    let mesh = Mesh2D::new(bed, dx, dx, manning);
    let mut states = Array2::from_elem((n_rows, n_cols), Conserved2D::new(h_n, q, 0.0));
    let bcs = Boundaries2D {
        west: Boundary::Discharge { q },
        east: Boundary::Depth { h: h_n },
        north: Boundary::Wall,
        south: Boundary::Wall,
    };
    let u = q / h_n;
    let c = (G * h_n).sqrt();
    let t_end = 2.0 * (n_cols as f64 * dx) / (u + c);
    let mut t = 0.0;
    while t < t_end {
        let dt = cfl_time_step(&states, &mesh, 0.4).min(t_end - t);
        ssprk2_step(&mut states, &mesh, bcs, dt);
        manning_friction_step(&mut states, &mesh, dt, 1.0e-9);
        t += dt;
    }
    let mid = n_rows / 2;
    let mut f = File::create(out)?;
    writeln!(f, "x,h_analytical,h_sim")?;
    for j in 0..n_cols {
        let x = (j as f64 + 0.5) * dx;
        let h_sim = states[(mid, j)].h;
        writeln!(f, "{x:.4},{h_n:.6},{h_sim:.6}")?;
    }
    Ok(())
}

// --- Thacker planar oscillation ------------------------------------

fn thacker_depth(x: f64, y: f64, t: f64, h0: f64, a: f64, b: f64) -> f64 {
    let omega = (2.0 * G * h0).sqrt() / a;
    let xc = b * (omega * t).cos();
    let yc = b * (omega * t).sin();
    let r2 = (x - xc).powi(2) + (y - yc).powi(2);
    (h0 / (a * a) * (a * a - r2)).max(0.0)
}

fn gen_thacker(out: &Path) -> std::io::Result<()> {
    let (h0, a, b) = (0.1_f64, 1.0_f64, 0.1_f64);
    let n = 80usize;
    let half_extent = 1.25_f64;
    let dx = 2.0 * half_extent / n as f64;
    let omega = (2.0 * G * h0).sqrt() / a;
    let period = 2.0 * std::f64::consts::PI / omega;
    let t_end = period / 2.0;

    let bed = Array2::from_shape_fn((n, n), |(i, j)| {
        let x = -half_extent + (j as f64 + 0.5) * dx;
        let y = -half_extent + (i as f64 + 0.5) * dx;
        h0 * ((x * x + y * y) / (a * a) - 1.0)
    });
    let mesh = Mesh2D::new(bed, dx, dx, 0.0);

    // Initial state at t = 0: cap centred at (B, 0), uniform velocity.
    let u0 = -b * omega * 0.0_f64.sin();
    let v0 = b * omega * 0.0_f64.cos();
    let mut states = Array2::from_shape_fn((n, n), |(i, j)| {
        let x = -half_extent + (j as f64 + 0.5) * dx;
        let y = -half_extent + (i as f64 + 0.5) * dx;
        let h = thacker_depth(x, y, 0.0, h0, a, b);
        Conserved2D::new(h, h * u0, h * v0)
    });

    let mut t = 0.0;
    while t < t_end {
        let dt = cfl_time_step(&states, &mesh, 0.4).min(t_end - t);
        ssprk2_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        t += dt;
    }

    // Centre row (y ≈ 0): the cap centre at t = T/2 sits at x = −B.
    let mid = n / 2;
    let y_mid = -half_extent + (mid as f64 + 0.5) * dx;
    let mut f = File::create(out)?;
    writeln!(f, "x,h_analytical,h_sim")?;
    for j in 0..n {
        let x = -half_extent + (j as f64 + 0.5) * dx;
        let h_an = thacker_depth(x, y_mid, t_end, h0, a, b);
        let h_sim = states[(mid, j)].h;
        writeln!(f, "{x:.4},{h_an:.6},{h_sim:.6}")?;
    }
    Ok(())
}

fn main() {
    fs::create_dir_all(OUT_DIR).expect("create data dir");
    let dir = Path::new(OUT_DIR);
    gen_stoker(&dir.join("verif_stoker.csv")).expect("stoker");
    println!("wrote verif_stoker.csv");
    gen_macdonald(&dir.join("verif_macdonald.csv")).expect("macdonald");
    println!("wrote verif_macdonald.csv");
    gen_thacker(&dir.join("verif_thacker.csv")).expect("thacker");
    println!("wrote verif_thacker.csv");
    println!("Verification profiles written to {OUT_DIR}/");
}
