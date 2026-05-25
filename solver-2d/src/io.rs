//! GeoTIFF I/O bridging `surtgis-core::Raster<f64>` ↔ solver `Mesh2D` + state.
//!
//! Two directions:
//!
//! - [`mesh_from_geotiff`]: read a DEM GeoTIFF (any GDAL-readable raster),
//!   wrap it as a [`Mesh2D`] with the user-supplied uniform Manning `n`,
//!   and return the original [`GeoTransform`] so output rasters land on
//!   the same grid.
//!
//! - [`write_depth_geotiff`]: take the solver's state array and write
//!   the depth band as a single-channel float32 GeoTIFF, with optional
//!   no-data masking for cells below the wet/dry threshold.
//!
//! Manning is uniform in this first iteration (matches the current
//! [`Mesh2D`] API). Variable Manning fields from a land-use raster are
//! Phase 2 (`mesh_from_geotiff_with_manning_raster`).
//!
//! # GeoTransform conventions
//!
//! GeoTIFFs (and the SurtGIS [`Raster`] wrapper) use the GDAL-style
//! affine transform: `(col, row) → (x, y)` with `pixel_height` typically
//! **negative** for north-up images. The solver internally uses
//! row-major indexing where row `i` increases top-to-bottom in image
//! space; this means the solver's "`+y` direction" maps to geographic
//! **south**. We take `dy = pixel_height.abs()` to give the solver a
//! positive spacing and rely on the symmetry of the SWE operator: the
//! sign of `y` is irrelevant to the physics, only the per-cell distance
//! `dy` matters. Output rasters written through [`write_depth_geotiff`]
//! carry the *same* [`GeoTransform`] as the input DEM, so the geographic
//! correspondence is preserved without any per-cell flipping.

use std::path::Path;

use ndarray::Array2;
use surtgis_core::Result;
use surtgis_core::io::{read_geotiff, write_geotiff};
use surtgis_core::raster::{GeoTransform, Raster};

use crate::geometry::Mesh2D;
use crate::state::Conserved2D;
use crate::H_DRY;

/// Read a DEM GeoTIFF as a [`Mesh2D`] with uniform Manning `n`.
///
/// Returns the mesh and the source [`GeoTransform`] so subsequent
/// output rasters from the same simulation can be georeferenced
/// consistently.
///
/// The Manning argument is uniform; variable Manning is Phase 2.
pub fn mesh_from_geotiff<P: AsRef<Path>>(
    path: P,
    manning: f64,
) -> Result<(Mesh2D, GeoTransform)> {
    let dem: Raster<f64> = read_geotiff(path, None)?;
    let transform = *dem.transform();
    let dx = transform.pixel_width.abs();
    let dy = transform.pixel_height.abs();
    let bed = dem.into_array();
    let mesh = Mesh2D::new(bed, dx, dy, manning);
    Ok((mesh, transform))
}

/// Build a `Raster<f64>` of cell depths from the solver state, with the
/// given [`GeoTransform`].
///
/// If `nodata` is `Some(v)`, cells with `h ≤ H_DRY` are written as `v`
/// (and `nodata` is also recorded as the raster's no-data sentinel for
/// downstream GIS tools). If `nodata` is `None`, dry cells retain their
/// numeric `h` value (typically 0 after the solver's dry-clamp step).
pub fn depth_raster_from_states(
    states: &Array2<Conserved2D>,
    transform: GeoTransform,
    nodata: Option<f64>,
) -> Raster<f64> {
    let depth = Array2::from_shape_fn(states.dim(), |(i, j)| {
        let h = states[(i, j)].h;
        match nodata {
            Some(nd) if h <= H_DRY => nd,
            _ => h,
        }
    });
    let mut raster = Raster::from_array(depth);
    raster.set_transform(transform);
    if let Some(nd) = nodata {
        raster.set_nodata(Some(nd));
    }
    raster
}

/// Write the depth field as a single-band GeoTIFF.
///
/// Convenience wrapper over [`depth_raster_from_states`] + the SurtGIS
/// `write_geotiff` writer. Uses the default `GeoTiffOptions`
/// (uncompressed float32).
pub fn write_depth_geotiff<P: AsRef<Path>>(
    path: P,
    states: &Array2<Conserved2D>,
    transform: GeoTransform,
    nodata: Option<f64>,
) -> Result<()> {
    let raster = depth_raster_from_states(states, transform, nodata);
    write_geotiff(&raster, path, None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn synthetic_dem() -> Raster<f64> {
        // 4×5 sloping bed: descends 1 m/cell in +col direction.
        let data = Array2::from_shape_fn((4, 5), |(_i, j)| 100.0 - j as f64);
        let mut r = Raster::from_array(data);
        r.set_transform(GeoTransform::new(
            500_000.0, // origin_x
            6_800_000.0, // origin_y
            30.0,        // pixel_width
            -30.0,       // pixel_height (negative, north-up)
        ));
        r
    }

    #[test]
    fn write_then_read_geotiff_preserves_data_and_transform() {
        let dem_in = synthetic_dem();
        let tmp = NamedTempFile::new().unwrap();
        write_geotiff(&dem_in, tmp.path(), None).unwrap();
        let dem_out: Raster<f64> = read_geotiff(tmp.path(), None).unwrap();

        assert_eq!(dem_in.shape(), dem_out.shape());
        let t_in = dem_in.transform();
        let t_out = dem_out.transform();
        assert!((t_in.origin_x - t_out.origin_x).abs() < 1e-9);
        assert!((t_in.origin_y - t_out.origin_y).abs() < 1e-9);
        assert!((t_in.pixel_width - t_out.pixel_width).abs() < 1e-9);
        assert!((t_in.pixel_height - t_out.pixel_height).abs() < 1e-9);
        for ((i, j), &v_in) in dem_in.data().indexed_iter() {
            let v_out = dem_out.data()[(i, j)];
            assert!(
                (v_in - v_out).abs() < 1e-5,
                "cell ({i},{j}): in={v_in}, out={v_out}"
            );
        }
    }

    #[test]
    fn mesh_from_geotiff_recovers_bed_shape_and_dx_dy() {
        let dem_in = synthetic_dem();
        let tmp = NamedTempFile::new().unwrap();
        write_geotiff(&dem_in, tmp.path(), None).unwrap();

        let (mesh, transform) = mesh_from_geotiff(tmp.path(), 0.035).unwrap();
        assert_eq!(mesh.n_rows(), 4);
        assert_eq!(mesh.n_cols(), 5);
        assert!((mesh.dx - 30.0).abs() < 1e-9);
        assert!((mesh.dy - 30.0).abs() < 1e-9);
        assert!((mesh.manning - 0.035).abs() < 1e-9);
        // Bed values match within float32 GeoTIFF precision.
        for ((i, j), &z_in) in dem_in.data().indexed_iter() {
            let z_mesh = mesh.bed[(i, j)];
            assert!(
                (z_in - z_mesh).abs() < 1e-5,
                "cell ({i},{j}): dem={z_in}, mesh={z_mesh}"
            );
        }
        // Transform passes through.
        assert!((transform.origin_x - 500_000.0).abs() < 1e-9);
        assert!((transform.pixel_height.abs() - 30.0).abs() < 1e-9);
    }

    #[test]
    fn depth_raster_writes_dry_cells_as_nodata_when_requested() {
        let mut states = Array2::from_elem((3, 3), Conserved2D::new(0.5, 0.0, 0.0));
        // Mark one cell as dry.
        states[(1, 1)] = Conserved2D::DRY;
        let transform = GeoTransform::new(0.0, 0.0, 1.0, -1.0);

        let r_with_nd = depth_raster_from_states(&states, transform, Some(-9999.0));
        assert_eq!(r_with_nd.data()[(0, 0)], 0.5);
        assert_eq!(r_with_nd.data()[(1, 1)], -9999.0);
        assert_eq!(r_with_nd.nodata(), Some(-9999.0));

        let r_no_nd = depth_raster_from_states(&states, transform, None);
        assert_eq!(r_no_nd.data()[(0, 0)], 0.5);
        assert_eq!(r_no_nd.data()[(1, 1)], 0.0);
        assert_eq!(r_no_nd.nodata(), None);
    }

    #[test]
    fn write_depth_geotiff_roundtrip_recovers_h() {
        let mut states = Array2::from_elem((3, 4), Conserved2D::new(1.25, 0.0, 0.0));
        states[(2, 3)] = Conserved2D::DRY;
        let transform = GeoTransform::new(100.0, 200.0, 5.0, -5.0);
        let tmp = NamedTempFile::new().unwrap();
        write_depth_geotiff(tmp.path(), &states, transform, Some(-9999.0)).unwrap();

        let read_back: Raster<f64> = read_geotiff(tmp.path(), None).unwrap();
        assert_eq!(read_back.shape(), (3, 4));
        for ((i, j), &v) in read_back.data().indexed_iter() {
            let expected = if (i, j) == (2, 3) { -9999.0 } else { 1.25 };
            assert!(
                (v - expected).abs() < 1e-5,
                "cell ({i},{j}): expected {expected}, got {v}"
            );
        }
    }
}
