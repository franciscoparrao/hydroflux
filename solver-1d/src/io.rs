//! GeoTIFF I/O for 1D channels and solver state via SurtGIS.
//!
//! A 1D channel is stored as a single-band 1×N GeoTIFF whose `pixel_width`
//! encodes `dx`. Cell `i` of the FV grid maps to column `i` of the raster.
//! Outputs (depth, discharge) are written as 1×N rasters with the same
//! geotransform as the input bed, so QGIS / `gdalinfo` align them
//! pixel-by-pixel for visualisation and post-processing.
//!
//! Manning roughness is not encoded in the raster — it is a scalar input
//! to [`read_channel`].
//!
//! # Precision note
//!
//! SurtGIS's **native** GeoTIFF writer stores data as `Float32` (single
//! precision, ~7 significant digits). The solver runs in `f64` and the
//! round-trip lose-recover error is bounded by `~1e-7` relative. For
//! visualisation in QGIS / processing in Python this is more than enough;
//! for bit-exact `f64` round-trips, enable SurtGIS's `gdal` feature.

use std::path::Path;

use ndarray::Array1;
use surtgis_core::Error as SurtgisError;
use surtgis_core::io::{read_geotiff, write_geotiff};
use surtgis_core::raster::{GeoTransform, Raster};

use crate::geometry::Channel1D;
use crate::state::Conserved;

/// Errors returned by the I/O module.
#[derive(Debug, thiserror::Error)]
pub enum IoError {
    /// Underlying SurtGIS error (file format, libtiff, geotransform parsing, …).
    #[error("SurtGIS error: {0}")]
    Surtgis(#[from] SurtgisError),
    /// The raster is not 1 row tall, so it cannot represent a 1D channel.
    #[error("raster must have exactly 1 row, got {0}")]
    NotOneDimensional(usize),
    /// `pixel_width` is not strictly positive.
    #[error("pixel_width must be strictly positive, got {0}")]
    InvalidDx(f64),
}

/// Convenience alias for `Result<T, IoError>`.
pub type Result<T> = std::result::Result<T, IoError>;

/// Read a 1D bed profile from a single-band GeoTIFF and build a
/// [`Channel1D`]. The raster must have exactly one row; `dx` is taken
/// from `pixel_width` of the geotransform. Manning roughness is passed
/// in by the caller (the raster does not carry it).
pub fn read_channel<P: AsRef<Path>>(path: P, manning: f64) -> Result<Channel1D> {
    let raster: Raster<f64> = read_geotiff(path, Some(1))?;
    if raster.rows() != 1 {
        return Err(IoError::NotOneDimensional(raster.rows()));
    }
    let dx = raster.transform().pixel_width;
    if dx <= 0.0 {
        return Err(IoError::InvalidDx(dx));
    }
    let row = raster
        .row(0)
        .expect("rows() >= 1 already verified, row(0) cannot fail");
    let bed = Array1::from_iter(row.iter().copied());
    Ok(Channel1D::new(bed, dx, manning))
}

/// Write the bed profile of `channel` as a single-band 1×N GeoTIFF.
/// Geotransform: `origin = (0, 0)`, `pixel_width = dx`, `pixel_height = -dx`
/// (north-up convention).
pub fn write_bed<P: AsRef<Path>>(channel: &Channel1D, path: P) -> Result<()> {
    let bed_vec: Vec<f64> = channel.bed.iter().copied().collect();
    write_row_geotiff(&bed_vec, channel.dx, path)
}

/// Write the depth `h` of each cell as a single-band 1×N GeoTIFF aligned
/// with the channel's geotransform.
pub fn write_depth<P: AsRef<Path>>(
    states: &[Conserved],
    channel: &Channel1D,
    path: P,
) -> Result<()> {
    assert_eq!(
        states.len(),
        channel.n_cells(),
        "states.len() ({}) must match channel.n_cells() ({})",
        states.len(),
        channel.n_cells()
    );
    let data: Vec<f64> = states.iter().map(|s| s.h).collect();
    write_row_geotiff(&data, channel.dx, path)
}

/// Write the unit discharge `hu` of each cell as a single-band 1×N
/// GeoTIFF aligned with the channel's geotransform.
pub fn write_discharge<P: AsRef<Path>>(
    states: &[Conserved],
    channel: &Channel1D,
    path: P,
) -> Result<()> {
    assert_eq!(
        states.len(),
        channel.n_cells(),
        "states.len() ({}) must match channel.n_cells() ({})",
        states.len(),
        channel.n_cells()
    );
    let data: Vec<f64> = states.iter().map(|s| s.hu).collect();
    write_row_geotiff(&data, channel.dx, path)
}

/// Internal helper: wrap a `Vec<f64>` of length `n` as a 1×N raster with
/// `pixel_width = dx`, `pixel_height = -dx`, and write it.
fn write_row_geotiff<P: AsRef<Path>>(data: &[f64], dx: f64, path: P) -> Result<()> {
    let n = data.len();
    let mut raster = Raster::<f64>::from_vec(data.to_vec(), 1, n)?;
    raster.set_transform(GeoTransform::new(0.0, 0.0, dx, -dx));
    write_geotiff(&raster, path, None)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;
    use tempfile::tempdir;

    #[test]
    fn round_trip_bed() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("bed.tif");

        let bed = array![10.0, 9.5, 9.0, 8.5, 8.0];
        let dx = 0.25;
        let manning = 0.03;
        let channel = Channel1D::new(bed.clone(), dx, manning);
        write_bed(&channel, &path).unwrap();

        let loaded = read_channel(&path, manning).unwrap();
        assert_eq!(loaded.n_cells(), 5);
        assert_relative_eq!(loaded.dx, dx, epsilon = 1e-6);
        assert_eq!(loaded.manning, manning);
        for i in 0..5 {
            assert_relative_eq!(loaded.bed[i], bed[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn round_trip_depth_and_discharge() {
        let dir = tempdir().unwrap();
        let path_h = dir.path().join("depth.tif");
        let path_hu = dir.path().join("discharge.tif");

        let channel = Channel1D::new(array![10.0, 9.5, 9.0, 8.5], 1.0, 0.03);
        let states = vec![
            Conserved::new(1.5, 0.7),
            Conserved::new(1.6, 0.8),
            Conserved::new(1.7, 0.9),
            Conserved::new(1.8, 1.0),
        ];

        write_depth(&states, &channel, &path_h).unwrap();
        write_discharge(&states, &channel, &path_hu).unwrap();

        // Read back as raw rasters and check the values + geotransform.
        let r_h: Raster<f64> = read_geotiff(&path_h, Some(1)).unwrap();
        let r_hu: Raster<f64> = read_geotiff(&path_hu, Some(1)).unwrap();
        for r in [&r_h, &r_hu] {
            assert_eq!(r.rows(), 1);
            assert_eq!(r.cols(), 4);
            assert_relative_eq!(r.transform().pixel_width, 1.0, epsilon = 1e-6);
        }
        let row_h = r_h.row(0).unwrap();
        let row_hu = r_hu.row(0).unwrap();
        for (i, s) in states.iter().enumerate() {
            assert_relative_eq!(row_h[i], s.h, epsilon = 1e-6);
            assert_relative_eq!(row_hu[i], s.hu, epsilon = 1e-6);
        }
    }

    #[test]
    fn rejects_multi_row_raster() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("2d.tif");
        // Build a 3×4 raster by hand and write it directly with SurtGIS.
        let raster = Raster::<f64>::from_vec(vec![0.0; 12], 3, 4).unwrap();
        write_geotiff(&raster, &path, None).unwrap();
        let err = read_channel(&path, 0.03).unwrap_err();
        match err {
            IoError::NotOneDimensional(rows) => assert_eq!(rows, 3),
            other => panic!("expected NotOneDimensional, got {other:?}"),
        }
    }
}
