//! GeoTIFF I/O bridging `surtgis-core::Raster<f64>` ↔ solver `Mesh2D` + state.
//!
//! Three directions:
//!
//! - [`mesh_from_geotiff`]: read a DEM GeoTIFF (any GDAL-readable raster),
//!   wrap it as a [`Mesh2D`] with a user-supplied **uniform** Manning `n`,
//!   and return the original [`GeoTransform`] so output rasters land on
//!   the same grid.
//!
//! - [`mesh_from_geotiff_with_landcover`]: same as above but reads a
//!   second raster of landcover class codes (`u8`) and maps each cell
//!   to its Manning `n` via a user-supplied lookup closure. See
//!   [`esa_worldcover_to_manning`] for a ready-made mapping for ESA
//!   WorldCover (10 m global landcover, 11 classes).
//!
//! - [`write_depth_geotiff`]: take the solver's state array and write
//!   the depth band as a single-channel float32 GeoTIFF, with optional
//!   no-data masking for cells below the wet/dry threshold.
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
/// consistently. For a landcover-derived Manning field instead of a
/// uniform value, use [`mesh_from_geotiff_with_landcover`].
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

/// Read a DEM GeoTIFF + a landcover GeoTIFF and build a [`Mesh2D`]
/// with a **spatially varying** Manning roughness field. The
/// `landcover_to_manning` closure maps each landcover class code
/// (`u8`) to its Manning `n` value [s/m^(1/3)].
///
/// The two rasters must have the same shape (number of rows × cols);
/// the function panics if they don't. The DEM's [`GeoTransform`] is
/// taken as the authoritative grid and returned to the caller.
///
/// For ESA WorldCover (10 m global landcover, 11 classes), pass
/// [`esa_worldcover_to_manning`] as the closure. For a custom mapping
/// (e.g. CORINE Land Cover, a calibrated channel-vs-overbank split),
/// supply your own `Fn(u8) -> f64`.
///
/// # Example
///
/// ```ignore
/// use hydroflux_solver_2d::io::{mesh_from_geotiff_with_landcover, esa_worldcover_to_manning};
///
/// let (mesh, transform) = mesh_from_geotiff_with_landcover(
///     "dem.tif",
///     "esa_worldcover.tif",
///     esa_worldcover_to_manning,
/// )?;
/// ```
pub fn mesh_from_geotiff_with_landcover<P, Q, F>(
    dem_path: P,
    landcover_path: Q,
    landcover_to_manning: F,
) -> Result<(Mesh2D, GeoTransform)>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    F: Fn(u8) -> f64,
{
    let dem: Raster<f64> = read_geotiff(dem_path, None)?;
    let landcover: Raster<u8> = read_geotiff(landcover_path, None)?;
    assert_eq!(
        dem.shape(),
        landcover.shape(),
        "DEM shape {:?} and landcover shape {:?} must match",
        dem.shape(),
        landcover.shape(),
    );
    let transform = *dem.transform();
    let dx = transform.pixel_width.abs();
    let dy = transform.pixel_height.abs();
    let bed = dem.into_array();
    let landcover_data = landcover.into_array();
    let manning_field =
        Array2::from_shape_fn(bed.dim(), |idx| landcover_to_manning(landcover_data[idx]));
    let mesh = Mesh2D::with_manning_field(bed, dx, dy, manning_field);
    Ok((mesh, transform))
}

/// ESA WorldCover (10 m global landcover, 11 classes) → Manning `n`
/// lookup. Values come from a compilation of Chow (1959), Arcement &
/// Schneider (1989), and Liu et al. (2019) for the most commonly
/// applied roughness coefficients in distributed hydrodynamic
/// simulations.
///
/// | Code | Class                    | n [s·m^(−1/3)] |
/// |------|--------------------------|----------------|
/// |   10 | Tree cover (dense forest)| 0.100          |
/// |   20 | Shrubland                | 0.060          |
/// |   30 | Grassland                | 0.040          |
/// |   40 | Cropland                 | 0.035          |
/// |   50 | Built-up (impervious)    | 0.015          |
/// |   60 | Bare / sparse vegetation | 0.025          |
/// |   70 | Snow and ice             | 0.030          |
/// |   80 | Permanent water bodies   | 0.030          |
/// |   90 | Herbaceous wetland       | 0.050          |
/// |   95 | Mangroves                | 0.100          |
/// |  100 | Moss and lichen          | 0.045          |
/// |other | (fallback)               | 0.040          |
///
/// The fallback `n = 0.040` matches the default for unclassified
/// natural channels (gravel-bed Andean rivers, our Huasco baseline).
/// Tune the table to the watershed's measured roughness if a
/// calibration is available.
pub fn esa_worldcover_to_manning(code: u8) -> f64 {
    match code {
        10 => 0.100,
        20 => 0.060,
        30 => 0.040,
        40 => 0.035,
        50 => 0.015,
        60 => 0.025,
        70 => 0.030,
        80 => 0.030,
        90 => 0.050,
        95 => 0.100,
        100 => 0.045,
        _ => 0.040,
    }
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
        for &n in mesh.manning.iter() {
            assert!((n - 0.035).abs() < 1e-9);
        }
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
    fn esa_worldcover_lookup_returns_expected_values() {
        // Spot-check the table against the documented mapping. The
        // exact numbers come from Chow (1959), Arcement & Schneider
        // (1989), and Liu et al. (2019).
        assert_eq!(esa_worldcover_to_manning(10), 0.100); // tree cover
        assert_eq!(esa_worldcover_to_manning(30), 0.040); // grassland
        assert_eq!(esa_worldcover_to_manning(50), 0.015); // built-up
        assert_eq!(esa_worldcover_to_manning(60), 0.025); // bare
        assert_eq!(esa_worldcover_to_manning(80), 0.030); // water
        assert_eq!(esa_worldcover_to_manning(95), 0.100); // mangroves
        // Unknown code falls back to the gravel-bed-channel default.
        assert_eq!(esa_worldcover_to_manning(7), 0.040);
        assert_eq!(esa_worldcover_to_manning(255), 0.040);
    }

    #[test]
    fn mesh_from_geotiff_with_landcover_maps_codes_to_manning() {
        // Synthetic 4×5 DEM + landcover raster with mixed codes:
        // mostly grassland (30 → 0.040) with one built-up cell
        // (50 → 0.015) and one tree-cover cell (10 → 0.100).
        let dem_in = synthetic_dem();
        let mut lc_data = Array2::<u8>::from_elem((4, 5), 30);
        lc_data[(0, 0)] = 50;
        lc_data[(3, 4)] = 10;
        let mut lc_raster = Raster::from_array(lc_data.clone());
        lc_raster.set_transform(*dem_in.transform());

        let dem_path = NamedTempFile::new().unwrap();
        let lc_path = NamedTempFile::new().unwrap();
        write_geotiff(&dem_in, dem_path.path(), None).unwrap();
        write_geotiff(&lc_raster, lc_path.path(), None).unwrap();

        let (mesh, _transform) = mesh_from_geotiff_with_landcover(
            dem_path.path(),
            lc_path.path(),
            esa_worldcover_to_manning,
        )
        .unwrap();
        assert_eq!(mesh.n_rows(), 4);
        assert_eq!(mesh.n_cols(), 5);
        assert!((mesh.manning[(0, 0)] - 0.015).abs() < 1e-12);
        assert!((mesh.manning[(3, 4)] - 0.100).abs() < 1e-12);
        for ((i, j), &n) in mesh.manning.indexed_iter() {
            if (i, j) == (0, 0) || (i, j) == (3, 4) {
                continue;
            }
            assert!((n - 0.040).abs() < 1e-12, "cell ({i},{j}) n = {n}");
        }
    }

    #[test]
    #[should_panic(expected = "DEM shape")]
    fn mesh_from_geotiff_with_landcover_panics_on_shape_mismatch() {
        // DEM 4×5 vs landcover 3×3 — must abort cleanly rather than
        // silently producing a zero-padded mesh.
        let dem_in = synthetic_dem();
        let mut lc_raster = Raster::from_array(Array2::<u8>::from_elem((3, 3), 30));
        lc_raster.set_transform(*dem_in.transform());
        let dem_path = NamedTempFile::new().unwrap();
        let lc_path = NamedTempFile::new().unwrap();
        write_geotiff(&dem_in, dem_path.path(), None).unwrap();
        write_geotiff(&lc_raster, lc_path.path(), None).unwrap();
        let _ = mesh_from_geotiff_with_landcover(
            dem_path.path(),
            lc_path.path(),
            esa_worldcover_to_manning,
        );
    }

    #[test]
    fn write_depth_geotiff_roundtrip_recovers_h() {
        let mut states = Array2::from_elem((3, 4), Conserved2D::new(1.25, 0.0, 0.0));
        states[(2, 3)] = Conserved2D::DRY;
        let transform = GeoTransform::new(100.0, 200.0, 5.0, -5.0);
        let tmp = NamedTempFile::new().unwrap();
        write_depth_geotiff(tmp.path(), &states, transform, Some(-9999.0)).unwrap();

        // SurtGIS ≥ Sprint-1 normalizes the nodata sentinel to NaN on
        // read (float rasters), so the dry cell written as -9999 must
        // come back as NaN, not as the on-disk sentinel.
        let read_back: Raster<f64> = read_geotiff(tmp.path(), None).unwrap();
        assert_eq!(read_back.shape(), (3, 4));
        for ((i, j), &v) in read_back.data().indexed_iter() {
            if (i, j) == (2, 3) {
                assert!(v.is_nan(), "dry cell ({i},{j}): expected NaN, got {v}");
            } else {
                assert!((v - 1.25).abs() < 1e-5, "cell ({i},{j}): expected 1.25, got {v}");
            }
        }
    }
}
