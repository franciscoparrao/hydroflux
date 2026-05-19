//! Boundary conditions for the 2D solver: 4 sides, 4 kinds, ghost cells.
//!
//! Boundary conditions are imposed via a single ghost cell adjacent to
//! each boundary face. The numerical flux at the boundary face is then
//! the standard Audusse well-balanced HLLC flux between the ghost cell
//! and the first inner cell.
//!
//! # Side convention
//!
//! With the GeoTIFF row-major layout (`bed[(i, j)]`, row `i` ↑ along
//! `+y`, column `j` ↑ along `+x`), the four sides are:
//!
//! | Side     | Boundary row/col            | Face is normal to | Outward normal |
//! |----------|-----------------------------|-------------------|----------------|
//! | `North`  | `i = 0`                     | `y`               | `−y`           |
//! | `South`  | `i = n_rows − 1`            | `y`               | `+y`           |
//! | `West`   | `j = 0`                     | `x`               | `−x`           |
//! | `East`   | `j = n_cols − 1`            | `x`               | `+x`           |
//!
//! With the typical north-up GeoTIFF orientation (`pixel_height < 0`),
//! `+y` in matrix space points to image south — so `Side::North` is
//! image-top, `Side::South` is image-bottom. The solver does not care
//! about the geographic interpretation; it only needs to know which
//! row/column to read.
//!
//! # Sign convention for [`Boundary::Discharge`]
//!
//! `q` is prescribed in the *coordinate-axis* direction, **not** as
//! "inflow positive". For `West`/`East` faces, `q` sets `hu` (so `q > 0`
//! means flow in `+x`, which is inflow at `West` and outflow at `East`).
//! For `North`/`South` faces, `q` sets `hv` (so `q > 0` means flow in
//! `+y`, inflow at `North`, outflow at `South`). This matches the 1D
//! convention.
//!
//! # Tangential momentum at physical BCs
//!
//! At an `x`-normal face the tangential momentum is `hv`; at a
//! `y`-normal face it is `hu`. `Discharge` and `Depth` extrapolate the
//! tangential component as zero-gradient (i.e. `tangential_ghost =
//! tangential_inner`). `Wall` reverses only the normal component. This
//! is consistent with the HLLC contact-wave property: the tangential
//! velocity is preserved through the contact wave anyway, so leaving
//! it zero-gradient at the ghost does not contaminate the interior.
//!
//! # Bed elevation at the ghost
//!
//! Computational BCs (`Transmissive`, `Wall`) extend the bed elevation
//! as zero-gradient. Physical BCs (`Discharge`, `Depth`) extend it
//! linearly across the boundary so that the boundary face carries the
//! same bed jump as an interior face. This is the same trick that
//! restored uniform-flow preservation in the 1D solver and is essential
//! for analytical steady-state benchmarks on a slope.

use crate::geometry::Mesh2D;
use crate::state::Conserved2D;

/// Boundary kind at one side of the 2D domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Boundary {
    /// Zero-gradient outflow: ghost equals the adjacent inner cell.
    /// Waves leave the domain (approximately) without reflecting.
    Transmissive,
    /// Reflective wall: the ghost mirrors the inner cell with reversed
    /// *normal* momentum, preserving the *tangential* momentum and the
    /// depth. Mass flux through the boundary is exactly zero by
    /// symmetry of the HLLC flux at a state of equal depth with
    /// opposite normal velocities.
    Wall,
    /// Prescribed unit discharge along the coordinate axis normal to the
    /// face: `hu = q` for `West`/`East` faces, `hv = q` for `North`/
    /// `South` faces. The tangential momentum and the depth are
    /// extrapolated as zero-gradient. Bed is extended linearly.
    Discharge {
        /// Prescribed normal-direction discharge `q` [m²/s] in the
        /// coordinate-axis convention (see module docs).
        q: f64,
    },
    /// Prescribed depth at the boundary. Both momentum components are
    /// extrapolated as zero-gradient. Bed is extended linearly.
    /// Intended for sub-critical outflow with a tailwater depth.
    Depth {
        /// Prescribed water depth [m] imposed in the ghost cell.
        h: f64,
    },
}

/// Which side of the rectangular domain a boundary lives on. Used by
/// [`ghost_cell`] to pick the correct inner / next-inner cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    /// Top of the matrix (`i = 0`); face is normal to `y`.
    North,
    /// Bottom of the matrix (`i = n_rows − 1`); face is normal to `y`.
    South,
    /// First column (`j = 0`); face is normal to `x`.
    West,
    /// Last column (`j = n_cols − 1`); face is normal to `x`.
    East,
}

impl Side {
    /// Is the face normal to the `x` axis (West/East)?
    pub const fn is_x_face(self) -> bool {
        matches!(self, Side::West | Side::East)
    }
}

/// Set of boundary conditions for all four sides of the 2D domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boundaries2D {
    /// Boundary on the `North` side (`i = 0`).
    pub north: Boundary,
    /// Boundary on the `South` side (`i = n_rows − 1`).
    pub south: Boundary,
    /// Boundary on the `West` side (`j = 0`).
    pub west: Boundary,
    /// Boundary on the `East` side (`j = n_cols − 1`).
    pub east: Boundary,
}

impl Boundaries2D {
    /// All four sides transmissive (open box). Default for prototype runs.
    pub const TRANSMISSIVE: Self = Self {
        north: Boundary::Transmissive,
        south: Boundary::Transmissive,
        west: Boundary::Transmissive,
        east: Boundary::Transmissive,
    };

    /// All four sides walls (closed box). Useful for mass-conservation
    /// tests: the total water volume must be preserved to machine
    /// precision modulo the time-integration error.
    pub const WALLS: Self = Self {
        north: Boundary::Wall,
        south: Boundary::Wall,
        west: Boundary::Wall,
        east: Boundary::Wall,
    };
}

/// Build the ghost cell adjacent to a given boundary face. Returns the
/// conservative state `(h, hu, hv)` and the bed elevation `z` of the
/// ghost. Both are needed by the Audusse hydrostatic reconstruction at
/// the boundary face.
///
/// `inner` is the adjacent inner-cell state. `idx` is the index *along*
/// the boundary: row `i` for `West`/`East` faces, column `j` for
/// `North`/`South` faces.
pub fn ghost_cell(
    mesh: &Mesh2D,
    inner: Conserved2D,
    kind: Boundary,
    side: Side,
    idx: usize,
) -> (Conserved2D, f64) {
    let state = ghost_state(inner, kind, side);
    let bed = ghost_bed(mesh, kind, side, idx);
    (state, bed)
}

fn ghost_state(inner: Conserved2D, kind: Boundary, side: Side) -> Conserved2D {
    match kind {
        Boundary::Transmissive => inner,
        Boundary::Wall => {
            // Reverse the normal component, preserve the tangential.
            if side.is_x_face() {
                Conserved2D {
                    h: inner.h,
                    hu: -inner.hu,
                    hv: inner.hv,
                }
            } else {
                Conserved2D {
                    h: inner.h,
                    hu: inner.hu,
                    hv: -inner.hv,
                }
            }
        }
        Boundary::Discharge { q } => {
            // Set the normal-direction momentum to q; keep depth and
            // tangential momentum zero-gradient.
            if side.is_x_face() {
                Conserved2D {
                    h: inner.h,
                    hu: q,
                    hv: inner.hv,
                }
            } else {
                Conserved2D {
                    h: inner.h,
                    hu: inner.hu,
                    hv: q,
                }
            }
        }
        Boundary::Depth { h } => Conserved2D {
            h,
            hu: inner.hu,
            hv: inner.hv,
        },
    }
}

fn ghost_bed(mesh: &Mesh2D, kind: Boundary, side: Side, idx: usize) -> f64 {
    let n_rows = mesh.n_rows();
    let n_cols = mesh.n_cols();
    match kind {
        Boundary::Transmissive | Boundary::Wall => match side {
            Side::North => mesh.bed[(0, idx)],
            Side::South => mesh.bed[(n_rows - 1, idx)],
            Side::West => mesh.bed[(idx, 0)],
            Side::East => mesh.bed[(idx, n_cols - 1)],
        },
        Boundary::Discharge { .. } | Boundary::Depth { .. } => match side {
            // Linear extrapolation across the boundary face. Requires at
            // least 2 cells in the normal direction; the assertion fires
            // only on programming errors (single-cell domains).
            Side::North => {
                debug_assert!(n_rows >= 2, "linear bed extrapolation requires n_rows ≥ 2");
                2.0 * mesh.bed[(0, idx)] - mesh.bed[(1, idx)]
            }
            Side::South => {
                debug_assert!(n_rows >= 2, "linear bed extrapolation requires n_rows ≥ 2");
                2.0 * mesh.bed[(n_rows - 1, idx)] - mesh.bed[(n_rows - 2, idx)]
            }
            Side::West => {
                debug_assert!(n_cols >= 2, "linear bed extrapolation requires n_cols ≥ 2");
                2.0 * mesh.bed[(idx, 0)] - mesh.bed[(idx, 1)]
            }
            Side::East => {
                debug_assert!(n_cols >= 2, "linear bed extrapolation requires n_cols ≥ 2");
                2.0 * mesh.bed[(idx, n_cols - 1)] - mesh.bed[(idx, n_cols - 2)]
            }
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    fn flat_mesh() -> Mesh2D {
        Mesh2D::new(
            array![[10.0, 10.0, 10.0], [10.0, 10.0, 10.0], [10.0, 10.0, 10.0],],
            1.0,
            1.0,
            0.03,
        )
    }

    /// Bed descends 0.01 m per metre in `+y` (matrix-row increases
    /// southward) and is uniform in `x`. Slope along `y` = 0.01.
    fn y_sloped_mesh() -> Mesh2D {
        Mesh2D::new(
            array![
                [0.0, 0.0, 0.0],
                [-0.01, -0.01, -0.01],
                [-0.02, -0.02, -0.02],
                [-0.03, -0.03, -0.03],
            ],
            1.0,
            1.0,
            0.03,
        )
    }

    /// Bed descends 0.01 m per metre in `+x`. Slope along `x` = 0.01.
    fn x_sloped_mesh() -> Mesh2D {
        Mesh2D::new(
            array![[0.0, -0.01, -0.02, -0.03], [0.0, -0.01, -0.02, -0.03],],
            1.0,
            1.0,
            0.03,
        )
    }

    #[test]
    fn transmissive_ghost_equals_inner_on_all_sides() {
        let mesh = flat_mesh();
        let u = Conserved2D::new(1.5, 0.6, -0.3);
        for side in [Side::North, Side::South, Side::West, Side::East] {
            let (g, _) = ghost_cell(&mesh, u, Boundary::Transmissive, side, 0);
            assert_eq!(g, u, "transmissive must be identity at side {side:?}");
        }
    }

    #[test]
    fn wall_on_x_face_reverses_hu_preserves_hv() {
        let mesh = flat_mesh();
        let u = Conserved2D::new(1.5, 0.6, -0.3);
        let (g, _) = ghost_cell(&mesh, u, Boundary::Wall, Side::West, 1);
        assert_eq!(g.h, u.h);
        assert_eq!(g.hu, -u.hu);
        assert_eq!(g.hv, u.hv);
    }

    #[test]
    fn wall_on_y_face_reverses_hv_preserves_hu() {
        let mesh = flat_mesh();
        let u = Conserved2D::new(1.5, 0.6, -0.3);
        let (g, _) = ghost_cell(&mesh, u, Boundary::Wall, Side::South, 1);
        assert_eq!(g.h, u.h);
        assert_eq!(g.hu, u.hu);
        assert_eq!(g.hv, -u.hv);
    }

    #[test]
    fn discharge_on_x_face_sets_hu_keeps_hv_zero_gradient() {
        let mesh = x_sloped_mesh();
        let inner = Conserved2D::new(0.6, 1.2, 0.4);
        let q_in = 0.8;
        let (g, _) = ghost_cell(&mesh, inner, Boundary::Discharge { q: q_in }, Side::West, 0);
        assert_relative_eq!(g.h, inner.h, epsilon = 1e-12);
        assert_relative_eq!(g.hu, q_in, epsilon = 1e-12);
        assert_relative_eq!(g.hv, inner.hv, epsilon = 1e-12);
    }

    #[test]
    fn discharge_on_y_face_sets_hv_keeps_hu_zero_gradient() {
        let mesh = y_sloped_mesh();
        let inner = Conserved2D::new(0.6, 1.2, 0.4);
        let q_in = 0.7;
        let (g, _) = ghost_cell(
            &mesh,
            inner,
            Boundary::Discharge { q: q_in },
            Side::North,
            0,
        );
        assert_relative_eq!(g.h, inner.h, epsilon = 1e-12);
        assert_relative_eq!(g.hu, inner.hu, epsilon = 1e-12);
        assert_relative_eq!(g.hv, q_in, epsilon = 1e-12);
    }

    #[test]
    fn depth_prescribes_h_keeps_both_momenta_zero_gradient() {
        let mesh = y_sloped_mesh();
        let inner = Conserved2D::new(0.6, 1.2, 0.4);
        let h_out = 0.45;
        let (g, _) = ghost_cell(&mesh, inner, Boundary::Depth { h: h_out }, Side::South, 2);
        assert_relative_eq!(g.h, h_out, epsilon = 1e-12);
        assert_relative_eq!(g.hu, inner.hu, epsilon = 1e-12);
        assert_relative_eq!(g.hv, inner.hv, epsilon = 1e-12);
    }

    #[test]
    fn computational_bcs_use_zero_gradient_bed_on_all_sides() {
        let mesh = y_sloped_mesh(); // varies in y, uniform in x
        let inner = Conserved2D::new(0.5, 0.0, 0.0);
        let n_rows = mesh.n_rows();
        let n_cols = mesh.n_cols();

        let (_, z_n) = ghost_cell(&mesh, inner, Boundary::Transmissive, Side::North, 1);
        assert_eq!(z_n, mesh.bed[(0, 1)]);

        let (_, z_s) = ghost_cell(&mesh, inner, Boundary::Wall, Side::South, 2);
        assert_eq!(z_s, mesh.bed[(n_rows - 1, 2)]);

        let (_, z_w) = ghost_cell(&mesh, inner, Boundary::Transmissive, Side::West, 0);
        assert_eq!(z_w, mesh.bed[(0, 0)]);

        let (_, z_e) = ghost_cell(&mesh, inner, Boundary::Wall, Side::East, 1);
        assert_eq!(z_e, mesh.bed[(1, n_cols - 1)]);
    }

    #[test]
    fn physical_bcs_extend_bed_linearly_in_y() {
        let mesh = y_sloped_mesh(); // bed descends 0.01 per row in +y
        let inner = Conserved2D::new(0.5, 0.0, 0.0);

        // Upstream (North): 2·bed[0,j] − bed[1,j] = 2·0 − (−0.01) = +0.01.
        let (_, z_north) = ghost_cell(&mesh, inner, Boundary::Discharge { q: 1.0 }, Side::North, 1);
        assert_relative_eq!(z_north, 0.01, epsilon = 1e-12);

        // Downstream (South): 2·bed[3,j] − bed[2,j] = 2·(−0.03) − (−0.02) = −0.04.
        let (_, z_south) = ghost_cell(&mesh, inner, Boundary::Depth { h: 0.5 }, Side::South, 1);
        assert_relative_eq!(z_south, -0.04, epsilon = 1e-12);
    }

    #[test]
    fn physical_bcs_extend_bed_linearly_in_x() {
        let mesh = x_sloped_mesh(); // bed descends 0.01 per column in +x
        let inner = Conserved2D::new(0.5, 0.0, 0.0);

        // West: 2·bed[i,0] − bed[i,1] = 2·0 − (−0.01) = +0.01.
        let (_, z_west) = ghost_cell(&mesh, inner, Boundary::Discharge { q: 1.0 }, Side::West, 0);
        assert_relative_eq!(z_west, 0.01, epsilon = 1e-12);

        // East: 2·bed[i,3] − bed[i,2] = 2·(−0.03) − (−0.02) = −0.04.
        let (_, z_east) = ghost_cell(&mesh, inner, Boundary::Depth { h: 0.5 }, Side::East, 1);
        assert_relative_eq!(z_east, -0.04, epsilon = 1e-12);
    }

    #[test]
    fn wall_at_rest_coincides_with_transmissive_on_all_sides() {
        // h > 0, both velocities zero: wall reflection is a no-op so the
        // ghost equals the inner cell, same as transmissive.
        let mesh = flat_mesh();
        let u = Conserved2D::new(2.0, 0.0, 0.0);
        for side in [Side::North, Side::South, Side::West, Side::East] {
            let (g_w, _) = ghost_cell(&mesh, u, Boundary::Wall, side, 1);
            let (g_t, _) = ghost_cell(&mesh, u, Boundary::Transmissive, side, 1);
            assert_eq!(g_w, g_t, "wall ≠ transmissive at rest on {side:?}");
        }
    }

    #[test]
    fn boundaries2d_constants_apply_uniform_kind() {
        assert_eq!(Boundaries2D::TRANSMISSIVE.north, Boundary::Transmissive);
        assert_eq!(Boundaries2D::TRANSMISSIVE.east, Boundary::Transmissive);
        assert_eq!(Boundaries2D::WALLS.south, Boundary::Wall);
        assert_eq!(Boundaries2D::WALLS.west, Boundary::Wall);
    }

    #[test]
    fn side_is_x_face_partitions_sides_correctly() {
        assert!(Side::West.is_x_face());
        assert!(Side::East.is_x_face());
        assert!(!Side::North.is_x_face());
        assert!(!Side::South.is_x_face());
    }
}
