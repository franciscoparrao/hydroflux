//! ESRI ASCII grid (`.asc`) reader, optionally gzip-compressed.
//!
//! A minimal, dependency-light reader for the ESRI ASCII grid format
//! used by LISFLOOD-FP and widely across the flood-modelling
//! community. Added to reproduce the official UK Environment Agency
//! 2D benchmark inputs (Néelz & Pender 2013, EA report SC120002)
//! directly from their distributed DEMs, without a GDAL round-trip.
//! GeoTIFF via [`crate::io`] stays the primary I/O path for
//! production use; this is a narrow, format-specific companion for
//! benchmark reproduction.
//!
//! Row 0 of the returned array is the **first data row after the
//! header** — the northernmost row for a standard north-up grid, same
//! as [`crate::io::mesh_from_geotiff`] and the `Side::North = i = 0`
//! convention documented in [`crate::boundary`]. No row flip is
//! needed to build a [`crate::geometry::Mesh2D`] from the result.

use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use flate2::read::GzDecoder;
use ndarray::Array2;
use thiserror::Error;

/// Header fields of an ESRI ASCII grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AsciiGridHeader {
    /// Number of columns (`x` direction).
    pub ncols: usize,
    /// Number of rows (`y` direction).
    pub nrows: usize,
    /// Grid origin `x`, south-west corner [m].
    pub xllcorner: f64,
    /// Grid origin `y`, south-west corner [m].
    pub yllcorner: f64,
    /// Uniform cell size [m].
    pub cellsize: f64,
    /// Sentinel value marking missing data.
    pub nodata_value: f64,
}

impl AsciiGridHeader {
    /// Map a geographic point `(x, y)` to the `(row, col)` of the
    /// solver-mesh cell containing it (row 0 = north, matching
    /// [`crate::geometry::Mesh2DG`]). Returns `None` if the point
    /// falls outside the grid extent.
    pub fn cell_at(&self, x: f64, y: f64) -> Option<(usize, usize)> {
        let col = (x - self.xllcorner) / self.cellsize;
        let row_from_south = (y - self.yllcorner) / self.cellsize;
        if col < 0.0 || row_from_south < 0.0 {
            return None;
        }
        let col = col.floor() as usize;
        let row_from_south = row_from_south.floor() as usize;
        if col >= self.ncols || row_from_south >= self.nrows {
            return None;
        }
        Some((self.nrows - 1 - row_from_south, col))
    }

    /// Rows (mesh convention, row 0 = north) whose cell extent
    /// overlaps the geographic `y`-range `[y_min, y_max)` — the rows a
    /// boundary segment restricted to that range should touch.
    /// Returns an empty vector if the range falls outside the grid.
    pub fn rows_overlapping_y_range(&self, y_min: f64, y_max: f64) -> Vec<usize> {
        let lo = ((y_min - self.yllcorner) / self.cellsize)
            .floor()
            .max(0.0) as usize;
        let hi_exclusive = (((y_max - self.yllcorner) / self.cellsize).ceil() as usize).min(self.nrows);
        (lo..hi_exclusive)
            .map(|row_from_south| self.nrows - 1 - row_from_south)
            .collect()
    }
}

/// Errors reading an ASCII grid.
#[derive(Debug, Error)]
pub enum AsciiGridError {
    /// Underlying I/O failure (file not found, gzip corruption, …).
    #[error("I/O error reading ASCII grid: {0}")]
    Io(#[from] io::Error),
    /// A header line was missing, out of order, or malformed.
    #[error("malformed ASCII grid header: {0}")]
    Header(String),
    /// The data section did not contain exactly `ncols * nrows` values.
    #[error("expected {expected} data values (ncols * nrows), found {found}")]
    ValueCount {
        /// Values required by the header.
        expected: usize,
        /// Values actually present.
        found: usize,
    },
    /// A token in the data section did not parse as `f64`.
    #[error("value {0:?} is not a valid number")]
    Parse(String),
}

/// Read an ESRI ASCII grid into a row-major `Array2<f64>` plus its
/// header. Transparently gunzips when `path` ends in `.gz`.
pub fn read_ascii_grid<P: AsRef<Path>>(
    path: P,
) -> Result<(Array2<f64>, AsciiGridHeader), AsciiGridError> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let mut text = String::new();
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        GzDecoder::new(file).read_to_string(&mut text)?;
    } else {
        std::io::BufReader::new(file).read_to_string(&mut text)?;
    }
    parse_ascii_grid(&text)
}

fn parse_ascii_grid(text: &str) -> Result<(Array2<f64>, AsciiGridHeader), AsciiGridError> {
    let mut lines = text.lines();
    let ncols = header_value(&mut lines, "ncols")? as usize;
    let nrows = header_value(&mut lines, "nrows")? as usize;
    let xllcorner = header_value(&mut lines, "xllcorner")?;
    let yllcorner = header_value(&mut lines, "yllcorner")?;
    let cellsize = header_value(&mut lines, "cellsize")?;
    let nodata_value = header_value(&mut lines, "NODATA_value")?;

    // The remaining text is whitespace-separated data values; ESRI
    // ASCII grid conventionally writes one row per text line, but
    // some writers wrap long rows — parsing as one flat token stream
    // is correct either way.
    let rest: String = lines.collect::<Vec<_>>().join(" ");
    let mut data = Vec::with_capacity(ncols * nrows);
    for tok in rest.split_whitespace() {
        data.push(
            tok.parse::<f64>()
                .map_err(|_| AsciiGridError::Parse(tok.to_string()))?,
        );
    }
    if data.len() != ncols * nrows {
        return Err(AsciiGridError::ValueCount {
            expected: ncols * nrows,
            found: data.len(),
        });
    }
    let array = Array2::from_shape_vec((nrows, ncols), data)
        .expect("length checked above to equal nrows * ncols");

    Ok((
        array,
        AsciiGridHeader {
            ncols,
            nrows,
            xllcorner,
            yllcorner,
            cellsize,
            nodata_value,
        },
    ))
}

fn header_value<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    expected_key: &str,
) -> Result<f64, AsciiGridError> {
    let line = lines
        .next()
        .ok_or_else(|| AsciiGridError::Header(format!("missing {expected_key}")))?;
    let mut parts = line.split_whitespace();
    let key = parts
        .next()
        .ok_or_else(|| AsciiGridError::Header(format!("empty line, expected {expected_key}")))?;
    if !key.eq_ignore_ascii_case(expected_key) {
        return Err(AsciiGridError::Header(format!(
            "expected key {expected_key:?}, got {key:?} (line: {line:?})"
        )));
    }
    let value = parts
        .next()
        .ok_or_else(|| AsciiGridError::Header(format!("missing value for {expected_key}")))?;
    value
        .parse::<f64>()
        .map_err(|_| AsciiGridError::Parse(value.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
ncols         4
nrows         3
xllcorner     0.0
yllcorner     0.0
cellsize      10
NODATA_value  -9999
1 2 3 4
5 6 7 8
9 10 11 12
";

    #[test]
    fn parses_header_fields() {
        let (_, h) = parse_ascii_grid(SAMPLE).unwrap();
        assert_eq!(h.ncols, 4);
        assert_eq!(h.nrows, 3);
        assert_eq!(h.cellsize, 10.0);
        assert_eq!(h.nodata_value, -9999.0);
    }

    #[test]
    fn row_zero_is_the_first_data_row_after_header() {
        // First data row after the header ("1 2 3 4") must land at
        // array row 0 — the north row under the mesh convention.
        let (data, _) = parse_ascii_grid(SAMPLE).unwrap();
        assert_eq!(data.row(0).to_vec(), vec![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(data.row(2).to_vec(), vec![9.0, 10.0, 11.0, 12.0]);
    }

    #[test]
    fn rejects_wrong_value_count() {
        let truncated = SAMPLE.replace("9 10 11 12\n", "");
        let err = parse_ascii_grid(&truncated).unwrap_err();
        assert!(matches!(err, AsciiGridError::ValueCount { .. }));
    }

    #[test]
    fn rejects_out_of_order_header() {
        let swapped = SAMPLE.replacen("ncols         4\n", "", 1);
        let err = parse_ascii_grid(&swapped).unwrap_err();
        assert!(matches!(err, AsciiGridError::Header(_)));
    }

    #[test]
    fn cell_at_maps_south_west_corner_to_bottom_row() {
        let (_, h) = parse_ascii_grid(SAMPLE).unwrap();
        // (5, 5) is inside the SW-most cell -> row_from_south = 0 ->
        // mesh row = nrows - 1 = 2 (the LAST array row, since array
        // row 0 is north).
        assert_eq!(h.cell_at(5.0, 5.0), Some((2, 0)));
        // (35, 25): row_from_south = floor(25/10) = 2 -> mesh row 0.
        assert_eq!(h.cell_at(35.0, 25.0), Some((0, 3)));
        assert_eq!(h.cell_at(-1.0, 5.0), None);
        assert_eq!(h.cell_at(5.0, 1000.0), None);
    }

    #[test]
    fn rows_overlapping_y_range_matches_expected_band() {
        // cellsize = 10, nrows = 3 -> south rows at y in [0,10),
        // [10,20), [20,30). Range [5, 25) overlaps south-rows 0,1,2
        // -> ALL rows (small grid); check on a case that excludes one.
        let (_, h) = parse_ascii_grid(SAMPLE).unwrap();
        let rows = h.rows_overlapping_y_range(12.0, 22.0);
        // south rows 1 (10-20) and 2 (20-30) overlap -> mesh rows
        // nrows-1-1=1 and nrows-1-2=0.
        let mut sorted = rows.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1]);
    }
}
