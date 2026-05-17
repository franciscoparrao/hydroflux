//! End-to-end smoke test of the SurtGIS I/O layer composed with the solver.
//!
//! Writes a synthetic bed profile as a GeoTIFF, reads it back as a
//! [`Channel1D`], runs the solver on the loaded channel, and writes the
//! outputs as GeoTIFFs. Verifies that the I/O layer interoperates
//! correctly with the solver across an actual file round-trip — the
//! typical user-facing workflow.

use hydroflux_solver_1d::{
    Boundaries, Boundary, Channel1D, Conserved, cfl_time_step, forward_euler_step,
    manning_friction_step, read_channel, write_bed, write_depth, write_discharge,
};
use ndarray::Array1;
use std::fs;
use tempfile::tempdir;

#[test]
fn write_bed_run_solver_write_outputs() {
    let dir = tempdir().unwrap();
    let bed_path = dir.path().join("bed.tif");
    let depth_path = dir.path().join("depth.tif");
    let discharge_path = dir.path().join("discharge.tif");

    // ---- 1. Synthesise a uniformly-sloped channel and write its bed.
    let n = 80;
    let dx = 1.0;
    let slope = 0.005;
    let manning = 0.03;
    let bed = Array1::from_iter((0..n).map(|i| -(i as f64) * dx * slope));
    let channel = Channel1D::new(bed, dx, manning);
    write_bed(&channel, &bed_path).unwrap();
    assert!(fs::metadata(&bed_path).unwrap().len() > 0);

    // ---- 2. Read it back; check the round-trip preserved the geometry.
    let channel_loaded = read_channel(&bed_path, manning).unwrap();
    assert_eq!(channel_loaded.n_cells(), n);
    assert!((channel_loaded.dx - dx).abs() < 1e-9);

    // ---- 3. Run the solver on the loaded channel at the Manning
    //         normal-depth equilibrium for a few wave transits.
    let q = 1.0;
    let h_n = (manning * q / slope.sqrt()).powf(3.0 / 5.0);
    let mut states: Vec<Conserved> = vec![Conserved::new(h_n, q); n];
    let bcs = Boundaries {
        left: Boundary::Discharge { q },
        right: Boundary::Depth { h: h_n },
    };
    let mut t = 0.0;
    let t_end = 10.0;
    while t < t_end {
        let dt = cfl_time_step(&states, dx, 0.4).min(t_end - t);
        forward_euler_step(&mut states, &channel_loaded, bcs, dt);
        manning_friction_step(&mut states, manning, dt, 1e-9);
        t += dt;
    }

    // ---- 4. Write outputs and verify the files exist and are non-empty.
    write_depth(&states, &channel_loaded, &depth_path).unwrap();
    write_discharge(&states, &channel_loaded, &discharge_path).unwrap();
    assert!(fs::metadata(&depth_path).unwrap().len() > 0);
    assert!(fs::metadata(&discharge_path).unwrap().len() > 0);

    // ---- 5. Round-trip the depth GeoTIFF back through `read_channel`
    //         (we abuse it to load any 1×N raster) and check the values
    //         landed close to the Manning equilibrium.
    let depth_loaded = read_channel(&depth_path, 0.0).unwrap();
    let max_dev: f64 = depth_loaded
        .bed
        .iter()
        .map(|&h| (h - h_n).abs() / h_n)
        .fold(0.0, f64::max);
    // Tolerance includes the f32 storage error and the solver's drift
    // over t_end; both are bounded well below 1 %.
    assert!(
        max_dev < 0.01,
        "depth round-trip max relative deviation = {max_dev:.4}"
    );
}
