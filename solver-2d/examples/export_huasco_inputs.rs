//! Export the Huasco setup as plain GeoTIFFs so an independent solver
//! can be driven from *exactly* the same inputs.
//!
//! This is step 0 of the SynxFlow/HiPIMS cross-validation
//! (`docs/xval-synxflow-huasco.md`). The point of that comparison is
//! that any difference between the two depth fields is attributable to
//! the numerical scheme, not to the setup — which only holds if both
//! solvers see the same friction. hydroflux derives its Manning field
//! from the ESA WorldCover raster through `esa_worldcover_to_manning`;
//! letting the other solver do its own land-cover mapping would
//! introduce a difference we could not separate from the scheme.
//!
//! So we write the *resolved* `n(x, y)` that hydroflux actually
//! integrates with, cell for cell, on the DEM grid:
//!
//!   examples/huasco_2d_phase2/output/huasco_manning.tif
//!
//! The DEM itself needs no export (both solvers read the same file),
//! but we re-emit it alongside so the hand-off directory is
//! self-contained.
//!
//! Run:
//!   cargo run --release -p hydroflux-solver-2d --example export_huasco_inputs

use std::path::PathBuf;

use surtgis_core::io::{read_geotiff, write_geotiff};
use surtgis_core::raster::Raster;

use hydroflux_solver_2d::{esa_worldcover_to_manning, mesh_from_geotiff_with_landcover};

const SUBSET_DEM: &str = "examples/huasco_2d_phase2/data/huasco_subset_dem.tif";
const SUBSET_LC: &str = "examples/huasco_2d_phase2/data/huasco_subset_landcover.tif";
const OUTPUT_DIR: &str = "examples/huasco_2d_phase2/output";

fn main() {
    let (mesh, transform) =
        mesh_from_geotiff_with_landcover(SUBSET_DEM, SUBSET_LC, esa_worldcover_to_manning)
            .expect("failed to load DEM + landcover");

    let n_rows = mesh.n_rows();
    let n_cols = mesh.n_cols();
    println!("Mesh: {n_rows}×{n_cols}, dx = {} m, dy = {} m", mesh.dx, mesh.dy);

    // The Manning field as hydroflux resolved it — not the land-cover
    // codes, and not a re-derivation. This is the array the solver
    // multiplies into the friction term.
    let mut manning = Raster::from_array(mesh.manning.clone());
    manning.set_transform(transform);
    let out_manning = PathBuf::from(OUTPUT_DIR).join("huasco_manning.tif");
    write_geotiff(&manning, &out_manning, None).expect("failed to write manning raster");

    let n_min = mesh.manning.iter().copied().fold(f64::INFINITY, f64::min);
    let n_max = mesh.manning.iter().copied().fold(0.0_f64, f64::max);
    let n_mean = mesh.manning.iter().sum::<f64>() / mesh.manning.len() as f64;
    println!(
        "Manning field written: n_min = {n_min:.4}, n_mean = {n_mean:.4}, n_max = {n_max:.4}"
    );
    println!("  {}", out_manning.display());

    // Re-emit the bed so the hand-off directory stands alone. Read back
    // from the source rather than from `mesh.bed` so any pit-filling or
    // no-data handling the loader applied is visible in the file the
    // other solver consumes.
    let dem: Raster<f64> = read_geotiff(SUBSET_DEM, None).expect("failed to re-read DEM");
    let out_dem = PathBuf::from(OUTPUT_DIR).join("huasco_dem_for_xval.tif");
    write_geotiff(&dem, &out_dem, None).expect("failed to write DEM copy");
    println!("  {}", out_dem.display());

    // The inflow cell in map coordinates, so the other solver can place
    // its source at the same physical location rather than trusting
    // that its own row/col indexing matches ours.
    let (inflow_row, inflow_col) = (135usize, 66usize);
    let x = transform.origin_x + (inflow_col as f64 + 0.5) * transform.pixel_width;
    let y = transform.origin_y + (inflow_row as f64 + 0.5) * transform.pixel_height;
    println!(
        "\nInflow point source: cell (row {inflow_row}, col {inflow_col}) \
         -> map ({x:.2}, {y:.2})"
    );
    println!("Boundaries: west Transmissive (outflow), N/S/E Wall.");
}
