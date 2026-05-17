//! `cargo run --release --example run_reach -- <bed.tif> [manning] [q] [t_end]`
//!
//! Run hydroflux-solver-1d on a 1D bed profile extracted from a real DEM.
//! Reads the bed GeoTIFF, integrates to steady state, and writes
//! `depth.tif` and `discharge.tif` next to the input.
//!
//! Defaults (Manning 0.04, `q = 3 m²/s`, `t_end = 5000 s`) target an Andean
//! foothill reach like the Río Maule tributary in `examples/maule_reach_demo`.
//! Override any of the three positional args for different scenarios — e.g.
//! Huasco-style semiarid reaches use higher Manning and lower `q`.

use std::path::PathBuf;

use hydroflux_solver_1d::{
    Boundaries, Boundary, Conserved, cfl_time_step, forward_euler_step, manning_friction_step,
    read_channel, write_depth, write_discharge,
};

const DEFAULT_MANNING: f64 = 0.04; // s/m^(1/3) — natural rocky channel
const DEFAULT_Q: f64 = 3.0; // m²/s unit discharge — moderate flood scenario
const DEFAULT_T_END: f64 = 5_000.0; // s — ≈ 1.5 wave transits over 10 km at c≈3 m/s
const CFL: f64 = 0.4;
const G: f64 = 9.81;

fn main() {
    let mut args = std::env::args().skip(1);
    let bed_path: PathBuf = args
        .next()
        .map(PathBuf::from)
        .expect("usage: run_reach <bed.tif> [manning] [q] [t_end]");
    let manning: f64 = args
        .next()
        .map(|s| s.parse().expect("manning must be a number"))
        .unwrap_or(DEFAULT_MANNING);
    let q: f64 = args
        .next()
        .map(|s| s.parse().expect("q must be a number"))
        .unwrap_or(DEFAULT_Q);
    let t_end: f64 = args
        .next()
        .map(|s| s.parse().expect("t_end must be a number"))
        .unwrap_or(DEFAULT_T_END);
    println!("Params: manning = {manning}, q = {q} m²/s, t_end = {t_end} s");
    let out_dir = bed_path
        .parent()
        .expect("bed path has no parent directory")
        .to_path_buf();

    println!("Loading bed from {}", bed_path.display());
    let channel = read_channel(&bed_path, manning).expect("failed to read bed");
    let n = channel.n_cells();
    let length_m = (n as f64) * channel.dx;
    println!(
        "Channel: {n} cells, dx = {:.2} m, length = {:.0} m",
        channel.dx, length_m
    );

    // Mean and downstream slopes to pick initial state and outflow BC.
    let mean_slope = (channel.bed[0] - channel.bed[n - 1]) / length_m;
    let h_n_mean = (manning * q / mean_slope.sqrt()).powf(3.0 / 5.0);
    let local_slope_ds = ((channel.bed[n - 2] - channel.bed[n - 1]) / channel.dx).max(1e-6);
    let h_n_ds = (manning * q / local_slope_ds.sqrt()).powf(3.0 / 5.0);
    println!(
        "Mean slope {:.4}, h_n_mean = {:.3} m, u_n = {:.2} m/s",
        mean_slope,
        h_n_mean,
        q / h_n_mean
    );
    println!(
        "Downstream slope {:.4}, h_n_ds = {:.3} m (used for outflow BC)",
        local_slope_ds, h_n_ds
    );

    let mut states: Vec<Conserved> = vec![Conserved::new(h_n_mean, q); n];
    let bcs = Boundaries {
        left: Boundary::Discharge { q },
        right: Boundary::Depth { h: h_n_ds },
    };

    let mut t = 0.0;
    let mut step = 0u64;
    while t < t_end {
        let dt = cfl_time_step(&states, channel.dx, CFL).min(t_end - t);
        forward_euler_step(&mut states, &channel, bcs, dt);
        manning_friction_step(&mut states, manning, dt, 1e-9);
        t += dt;
        step += 1;
        if step % 200 == 0 {
            let max_h = states.iter().map(|s| s.h).fold(0.0_f64, f64::max);
            let max_fr = states
                .iter()
                .map(|s| (s.hu / s.h).abs() / (G * s.h).sqrt())
                .fold(0.0_f64, f64::max);
            println!("  step {step:5}, t = {t:6.0} s, max h = {max_h:.3} m, max Fr = {max_fr:.3}");
        }
    }
    println!(
        "Solver finished after {} steps, simulated t = {:.0} s",
        step, t
    );

    let depth_path = out_dir.join("depth.tif");
    let discharge_path = out_dir.join("discharge.tif");
    write_depth(&states, &channel, &depth_path).expect("write_depth failed");
    write_discharge(&states, &channel, &discharge_path).expect("write_discharge failed");
    println!("Wrote {}", depth_path.display());
    println!("Wrote {}", discharge_path.display());
}
