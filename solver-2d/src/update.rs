//! Time-stepping for the 2D Saint-Venant solver.
//!
//! Explicit forward Euler in time, finite-volume in space, with the
//! HLLC interface flux and the hydrostatic reconstruction of Audusse
//! et al. (2004) extended per face direction. Manning friction is
//! applied separately as an operator-split fractional step (see
//! `crate::source`, pending).
//!
//! # Audusse in 2D
//!
//! At each face the reconstruction is one-dimensional in the
//! face-normal direction (the SWE flux only couples bed-slope source
//! to the normal-direction momentum). On an `x`-face between cells
//! `L = (i, j−1)` and `R = (i, j)`, we set:
//!
//! ```text
//!   z_max = max(z_L, z_R)
//!   h*_L  = max(h_L + z_L − z_max, 0)
//!   h*_R  = max(h_R + z_R − z_max, 0)
//! ```
//!
//! The HLLC flux is then evaluated on the reconstructed states
//! `(h*, h*·u, h*·v)` (tangential velocity is unchanged by the
//! reconstruction). The hydrostatic-pressure correction
//! `(g/2)(h² − h*²)` is added to the `x`-momentum component on each
//! side; **mass and `y`-momentum are not corrected** because the
//! bed-slope source enters only the `x`-momentum equation at an
//! `x`-face. Symmetric statement holds for `y`-faces (correction goes
//! to `y`-momentum).
//!
//! # CFL condition
//!
//! For an unsplit explicit FV scheme in 2D the time step is bounded by
//!
//! ```text
//!   dt · (s_x_max / dx + s_y_max / dy) ≤ cfl
//! ```
//!
//! with `s_x_max = max(|u| + c)`, `s_y_max = max(|v| + c)`. We return
//! `f64::INFINITY` for an entirely dry domain.

use crate::boundary::{Boundaries2D, Side, ghost_cell};
use crate::flux::{FluxX, FluxY};
use crate::geometry::Mesh2D;
use crate::riemann::{hllc_flux_x, hllc_flux_y};
use crate::state::Conserved2D;
use crate::{GRAVITY, H_DRY};
use ndarray::Array2;

/// Maximum signal speeds `(s_x, s_y)` across the state field, where
/// `s_x = max(|u| + c)` and `s_y = max(|v| + c)`. Returns `(0, 0)` for
/// an empty or all-dry state.
///
/// Dry cells (`h ≤ H_DRY`) contribute zero to both maxima — the CFL
/// must not be tightened by spurious `hu / h` blow-ups in essentially
/// dry cells.
pub fn max_wave_speeds(states: &Array2<Conserved2D>) -> (f64, f64) {
    let mut s_x = 0.0_f64;
    let mut s_y = 0.0_f64;
    for s in states {
        if s.h <= H_DRY {
            continue;
        }
        let c = (GRAVITY * s.h).sqrt();
        let u = s.hu / s.h;
        let v = s.hv / s.h;
        s_x = s_x.max(u.abs() + c);
        s_y = s_y.max(v.abs() + c);
    }
    (s_x, s_y)
}

/// CFL-bounded time step `dt = cfl / (s_x/dx + s_y/dy)`. Returns
/// `f64::INFINITY` when the domain is entirely dry (no signal can
/// propagate); callers should clamp against a problem-specific maximum.
///
/// `cfl` is typically 0.4–0.5 for an explicit FV solver with HLLC in 2D.
pub fn cfl_time_step(states: &Array2<Conserved2D>, mesh: &Mesh2D, cfl: f64) -> f64 {
    let (s_x, s_y) = max_wave_speeds(states);
    let denom = s_x / mesh.dx + s_y / mesh.dy;
    if denom > 0.0 {
        cfl / denom
    } else {
        f64::INFINITY
    }
}

/// Numerical fluxes at a single `x`-face in the well-balanced
/// two-sided formulation. `minus` is consumed by the cell on the LEFT
/// of the face (lower column index); `plus` by the cell on the RIGHT.
/// They share `mass` and `y_momentum` but differ in `x_momentum` by the
/// hydrostatic-pressure correction.
#[derive(Debug, Clone, Copy)]
struct FaceFluxX {
    minus: FluxX,
    plus: FluxX,
}

/// Numerical fluxes at a single `y`-face. `minus` is consumed by the
/// cell ABOVE the face (lower row index); `plus` by the cell BELOW.
/// They share `mass` and `x_momentum` but differ in `y_momentum`.
#[derive(Debug, Clone, Copy)]
struct FaceFluxY {
    minus: FluxY,
    plus: FluxY,
}

/// Audusse well-balanced HLLC flux on an `x`-face. See module docs for
/// the reconstruction. On a flat bed (`z_left == z_right`) the
/// correction vanishes and the flux is the plain HLLC flux of the
/// original states.
fn well_balanced_x_face(
    left: Conserved2D,
    z_left: f64,
    right: Conserved2D,
    z_right: f64,
) -> FaceFluxX {
    let z_max = z_left.max(z_right);
    let h_star_left = (left.h + z_left - z_max).max(0.0);
    let h_star_right = (right.h + z_right - z_max).max(0.0);

    let (u_left, v_left) = if left.h > 0.0 {
        (left.hu / left.h, left.hv / left.h)
    } else {
        (0.0, 0.0)
    };
    let (u_right, v_right) = if right.h > 0.0 {
        (right.hu / right.h, right.hv / right.h)
    } else {
        (0.0, 0.0)
    };

    let recon_left = Conserved2D::new(h_star_left, h_star_left * u_left, h_star_left * v_left);
    let recon_right =
        Conserved2D::new(h_star_right, h_star_right * u_right, h_star_right * v_right);

    let f = hllc_flux_x(recon_left, recon_right);
    let half_g = 0.5 * GRAVITY;
    let corr_left = half_g * (left.h * left.h - h_star_left * h_star_left);
    let corr_right = half_g * (right.h * right.h - h_star_right * h_star_right);
    FaceFluxX {
        minus: FluxX {
            mass: f.mass,
            x_momentum: f.x_momentum + corr_left,
            y_momentum: f.y_momentum,
        },
        plus: FluxX {
            mass: f.mass,
            x_momentum: f.x_momentum + corr_right,
            y_momentum: f.y_momentum,
        },
    }
}

/// Audusse well-balanced HLLC flux on a `y`-face. `left` (alias top,
/// lower row index) is at `z_left`; `right` (bottom, higher row index)
/// is at `z_right`.
fn well_balanced_y_face(
    left: Conserved2D,
    z_left: f64,
    right: Conserved2D,
    z_right: f64,
) -> FaceFluxY {
    let z_max = z_left.max(z_right);
    let h_star_left = (left.h + z_left - z_max).max(0.0);
    let h_star_right = (right.h + z_right - z_max).max(0.0);

    let (u_left, v_left) = if left.h > 0.0 {
        (left.hu / left.h, left.hv / left.h)
    } else {
        (0.0, 0.0)
    };
    let (u_right, v_right) = if right.h > 0.0 {
        (right.hu / right.h, right.hv / right.h)
    } else {
        (0.0, 0.0)
    };

    let recon_left = Conserved2D::new(h_star_left, h_star_left * u_left, h_star_left * v_left);
    let recon_right =
        Conserved2D::new(h_star_right, h_star_right * u_right, h_star_right * v_right);

    let f = hllc_flux_y(recon_left, recon_right);
    let half_g = 0.5 * GRAVITY;
    let corr_left = half_g * (left.h * left.h - h_star_left * h_star_left);
    let corr_right = half_g * (right.h * right.h - h_star_right * h_star_right);
    FaceFluxY {
        minus: FluxY {
            mass: f.mass,
            x_momentum: f.x_momentum,
            y_momentum: f.y_momentum + corr_left,
        },
        plus: FluxY {
            mass: f.mass,
            x_momentum: f.x_momentum,
            y_momentum: f.y_momentum + corr_right,
        },
    }
}

/// One forward-Euler update of the 2D FV solution. Modifies `states` in
/// place.
///
/// Includes the well-balanced bed-slope source via per-face hydrostatic
/// reconstruction. Friction is **not** included — call the Manning
/// friction step (forthcoming) separately.
///
/// Panics if the shape of `states` does not match `(mesh.n_rows(),
/// mesh.n_cols())`. The caller is responsible for keeping `dt` below
/// the CFL bound (see [`cfl_time_step`]).
pub fn forward_euler_step(
    states: &mut Array2<Conserved2D>,
    mesh: &Mesh2D,
    bcs: Boundaries2D,
    dt: f64,
) {
    let n_rows = mesh.n_rows();
    let n_cols = mesh.n_cols();
    assert_eq!(
        states.shape(),
        [n_rows, n_cols],
        "states shape {:?} must match mesh ({}, {})",
        states.shape(),
        n_rows,
        n_cols,
    );
    if n_rows == 0 || n_cols == 0 {
        return;
    }

    // Precompute all x-faces and y-faces. x-faces have shape
    // (n_rows, n_cols + 1); y-faces have shape (n_rows + 1, n_cols).
    // Boundary faces use ghost cells from `bcs`.
    let faces_x = Array2::<FaceFluxX>::from_shape_fn((n_rows, n_cols + 1), |(i, j)| {
        if j == 0 {
            let (g, z_g) = ghost_cell(mesh, states[(i, 0)], bcs.west, Side::West, i);
            well_balanced_x_face(g, z_g, states[(i, 0)], mesh.bed[(i, 0)])
        } else if j == n_cols {
            let (g, z_g) = ghost_cell(mesh, states[(i, n_cols - 1)], bcs.east, Side::East, i);
            well_balanced_x_face(states[(i, n_cols - 1)], mesh.bed[(i, n_cols - 1)], g, z_g)
        } else {
            well_balanced_x_face(
                states[(i, j - 1)],
                mesh.bed[(i, j - 1)],
                states[(i, j)],
                mesh.bed[(i, j)],
            )
        }
    });

    let faces_y = Array2::<FaceFluxY>::from_shape_fn((n_rows + 1, n_cols), |(i, j)| {
        if i == 0 {
            let (g, z_g) = ghost_cell(mesh, states[(0, j)], bcs.north, Side::North, j);
            well_balanced_y_face(g, z_g, states[(0, j)], mesh.bed[(0, j)])
        } else if i == n_rows {
            let (g, z_g) = ghost_cell(mesh, states[(n_rows - 1, j)], bcs.south, Side::South, j);
            well_balanced_y_face(states[(n_rows - 1, j)], mesh.bed[(n_rows - 1, j)], g, z_g)
        } else {
            well_balanced_y_face(
                states[(i - 1, j)],
                mesh.bed[(i - 1, j)],
                states[(i, j)],
                mesh.bed[(i, j)],
            )
        }
    });

    // FV update. For cell (i, j):
    //   right x-face is faces_x[(i, j+1)] — cell is on its LEFT side → .minus
    //   left  x-face is faces_x[(i, j)]   — cell is on its RIGHT side → .plus
    //   bottom y-face is faces_y[(i+1, j)] — cell is on its TOP side → .minus
    //   top    y-face is faces_y[(i, j)]   — cell is on its BOTTOM side → .plus
    //
    // Positivity preservation: cells whose updated depth would fall
    // at or below H_DRY are clamped to DRY (depth = 0, both momentum
    // components zeroed). This is the simplest wetting/drying
    // treatment — see Liang & Marche (2009) for the flux-rescaling
    // alternative that is strictly mass-conservative. The clamp can
    // lose a small amount of mass at the wet/dry front (bounded
    // pointwise by H_DRY per cell, globally by H_DRY · dx · dy per
    // step times the number of cells that crossed the threshold);
    // this is acceptable for first-iteration robustness and will be
    // tightened with flux-rescaling when needed.
    let dt_dx = dt / mesh.dx;
    let dt_dy = dt / mesh.dy;
    for i in 0..n_rows {
        for j in 0..n_cols {
            let fx_right = faces_x[(i, j + 1)].minus;
            let fx_left = faces_x[(i, j)].plus;
            let fy_bottom = faces_y[(i + 1, j)].minus;
            let fy_top = faces_y[(i, j)].plus;

            let dh =
                dt_dx * (fx_right.mass - fx_left.mass) + dt_dy * (fy_bottom.mass - fy_top.mass);
            let new_h = states[(i, j)].h - dh;
            if new_h <= H_DRY {
                states[(i, j)] = Conserved2D::DRY;
            } else {
                let dhu = dt_dx * (fx_right.x_momentum - fx_left.x_momentum)
                    + dt_dy * (fy_bottom.x_momentum - fy_top.x_momentum);
                let dhv = dt_dx * (fx_right.y_momentum - fx_left.y_momentum)
                    + dt_dy * (fy_bottom.y_momentum - fy_top.y_momentum);
                states[(i, j)].h = new_h;
                states[(i, j)].hu -= dhu;
                states[(i, j)].hv -= dhv;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::assert_relative_eq;
    use ndarray::Array2;

    fn flat_mesh(n_rows: usize, n_cols: usize, dx: f64, dy: f64) -> Mesh2D {
        Mesh2D::new(Array2::<f64>::zeros((n_rows, n_cols)), dx, dy, 0.0)
    }

    fn x_sloped_mesh(n_rows: usize, n_cols: usize, dx: f64, dy: f64, slope: f64) -> Mesh2D {
        let bed = Array2::from_shape_fn((n_rows, n_cols), |(_i, j)| -(j as f64) * dx * slope);
        Mesh2D::new(bed, dx, dy, 0.0)
    }

    fn y_sloped_mesh(n_rows: usize, n_cols: usize, dx: f64, dy: f64, slope: f64) -> Mesh2D {
        let bed = Array2::from_shape_fn((n_rows, n_cols), |(i, _j)| -(i as f64) * dy * slope);
        Mesh2D::new(bed, dx, dy, 0.0)
    }

    fn diagonal_sloped_mesh(
        n_rows: usize,
        n_cols: usize,
        dx: f64,
        dy: f64,
        sx: f64,
        sy: f64,
    ) -> Mesh2D {
        let bed = Array2::from_shape_fn((n_rows, n_cols), |(i, j)| {
            -(j as f64) * dx * sx - (i as f64) * dy * sy
        });
        Mesh2D::new(bed, dx, dy, 0.0)
    }

    fn gaussian_bed(n_rows: usize, n_cols: usize, dx: f64, dy: f64, amp: f64) -> Mesh2D {
        // Centre on cell ((n_rows-1)/2, (n_cols-1)/2): puts the peak
        // exactly on a cell for odd `n`, and keeps the bump symmetric
        // under index reflection for any `n`.
        let cx = (n_cols as f64 - 1.0) / 2.0;
        let cy = (n_rows as f64 - 1.0) / 2.0;
        let w_sq = ((n_rows.min(n_cols)) as f64).powi(2) / 25.0;
        let bed = Array2::from_shape_fn((n_rows, n_cols), |(i, j)| {
            let dxx = j as f64 - cx;
            let dyy = i as f64 - cy;
            amp * (-(dxx * dxx + dyy * dyy) / w_sq).exp()
        });
        Mesh2D::new(bed, dx, dy, 0.0)
    }

    fn lake_at_rest_on(mesh: &Mesh2D, eta: f64) -> Array2<Conserved2D> {
        Array2::from_shape_fn((mesh.n_rows(), mesh.n_cols()), |(i, j)| {
            Conserved2D::new(eta - mesh.bed[(i, j)], 0.0, 0.0)
        })
    }

    fn total_mass(states: &Array2<Conserved2D>, dx: f64, dy: f64) -> f64 {
        states.iter().map(|s| s.h * dx * dy).sum()
    }

    fn gaussian_bump_2d(
        n_rows: usize,
        n_cols: usize,
        h_base: f64,
        h_amp: f64,
    ) -> Array2<Conserved2D> {
        // Centre on cell ((n_rows-1)/2, (n_cols-1)/2): symmetric under
        // index reflection so isotropy tests are well-posed.
        let cx = (n_cols as f64 - 1.0) / 2.0;
        let cy = (n_rows as f64 - 1.0) / 2.0;
        let w_sq = ((n_rows.min(n_cols)) as f64).powi(2) / 25.0;
        Array2::from_shape_fn((n_rows, n_cols), |(i, j)| {
            let dxx = j as f64 - cx;
            let dyy = i as f64 - cy;
            let h = h_base + h_amp * (-(dxx * dxx + dyy * dyy) / w_sq).exp();
            Conserved2D::new(h, 0.0, 0.0)
        })
    }

    #[test]
    fn max_wave_speeds_zero_for_dry_domain() {
        let states = Array2::from_elem((4, 5), Conserved2D::DRY);
        assert_eq!(max_wave_speeds(&states), (0.0, 0.0));
    }

    #[test]
    fn max_wave_speeds_match_textbook_formula() {
        let mut states = Array2::from_elem((2, 2), Conserved2D::DRY);
        states[(0, 0)] = Conserved2D::new(1.0, 2.0, 3.0); // u=2, v=3
        let (sx, sy) = max_wave_speeds(&states);
        let c = (GRAVITY * 1.0).sqrt();
        assert_relative_eq!(sx, 2.0 + c, epsilon = 1e-12);
        assert_relative_eq!(sy, 3.0 + c, epsilon = 1e-12);
    }

    #[test]
    fn cfl_time_step_dry_domain_returns_infinity() {
        let states = Array2::from_elem((4, 5), Conserved2D::DRY);
        let mesh = flat_mesh(4, 5, 1.0, 1.0);
        assert_eq!(cfl_time_step(&states, &mesh, 0.5), f64::INFINITY);
    }

    #[test]
    fn lake_at_rest_on_flat_bed_is_preserved_exactly() {
        let n_rows = 10;
        let n_cols = 12;
        let mesh = flat_mesh(n_rows, n_cols, 1.0, 1.0);
        let mut states = Array2::from_elem((n_rows, n_cols), Conserved2D::new(2.0, 0.0, 0.0));
        let dt = 0.01;
        for _ in 0..100 {
            forward_euler_step(&mut states, &mesh, Boundaries2D::TRANSMISSIVE, dt);
        }
        for s in &states {
            assert_relative_eq!(s.h, 2.0, epsilon = 1e-12);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-12);
            assert_relative_eq!(s.hv, 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn lake_at_rest_on_x_sloped_bed_is_preserved() {
        // η = h + z = const, u = v = 0, bed varies in x only.
        let n_rows = 6;
        let n_cols = 30;
        let dx = 1.0;
        let dy = 1.0;
        let mesh = x_sloped_mesh(n_rows, n_cols, dx, dy, 0.05);
        let eta = 5.0;
        let initial = lake_at_rest_on(&mesh, eta);
        let mut states = initial.clone();

        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        for ((i, j), s) in states.indexed_iter() {
            assert_relative_eq!(s.h, initial[(i, j)].h, epsilon = 1e-10);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-10);
            assert_relative_eq!(s.hv, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn lake_at_rest_on_y_sloped_bed_is_preserved() {
        // Same as above but the slope is along y. Catches bugs where
        // the y-face hydrostatic correction was wired to the wrong
        // momentum component.
        let n_rows = 30;
        let n_cols = 6;
        let mesh = y_sloped_mesh(n_rows, n_cols, 1.0, 1.0, 0.05);
        let eta = 5.0;
        let initial = lake_at_rest_on(&mesh, eta);
        let mut states = initial.clone();

        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        for ((i, j), s) in states.indexed_iter() {
            assert_relative_eq!(s.h, initial[(i, j)].h, epsilon = 1e-10);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-10);
            assert_relative_eq!(s.hv, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn lake_at_rest_on_diagonal_bed_is_preserved() {
        // η = const on a bed sloping in BOTH x and y. This is the
        // canonical 2D well-balanced test — both face directions must
        // cancel their hydrostatic-pressure term against the bed jump
        // simultaneously.
        let n_rows = 15;
        let n_cols = 20;
        let mesh = diagonal_sloped_mesh(n_rows, n_cols, 1.0, 1.0, 0.03, 0.04);
        let eta = 5.0;
        let initial = lake_at_rest_on(&mesh, eta);
        let mut states = initial.clone();

        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        for ((i, j), s) in states.indexed_iter() {
            assert_relative_eq!(s.h, initial[(i, j)].h, epsilon = 1e-10);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-10);
            assert_relative_eq!(s.hv, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn lake_at_rest_on_bumpy_bed_is_preserved() {
        // η constant above a submerged 2D Gaussian hill. Non-monotonic
        // bed in both directions; catches asymmetric reconstruction bugs.
        let n_rows = 20;
        let n_cols = 20;
        let mesh = gaussian_bed(n_rows, n_cols, 1.0, 1.0, 1.5);
        let eta = 3.0; // safely above max(bed) = 1.5
        let initial = lake_at_rest_on(&mesh, eta);
        let mut states = initial.clone();

        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        for ((i, j), s) in states.indexed_iter() {
            assert_relative_eq!(s.h, initial[(i, j)].h, epsilon = 1e-10);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-10);
            assert_relative_eq!(s.hv, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn mass_is_conserved_exactly_with_wall_boundaries() {
        // 2D Gaussian bump on a flat bed in a closed box. Total volume
        // h·dx·dy summed over the domain must be preserved to roundoff.
        let n_rows = 25;
        let n_cols = 25;
        let dx = 1.0;
        let dy = 1.0;
        let mesh = flat_mesh(n_rows, n_cols, dx, dy);
        let mut states = gaussian_bump_2d(n_rows, n_cols, 1.0, 0.5);
        let m0 = total_mass(&states, dx, dy);

        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        let m1 = total_mass(&states, dx, dy);
        assert_relative_eq!(m0, m1, epsilon = 1e-10);
    }

    #[test]
    fn bump_remains_bounded_and_finite_with_transmissive_bc() {
        let n_rows = 20;
        let n_cols = 20;
        let mesh = flat_mesh(n_rows, n_cols, 1.0, 1.0);
        let mut states = gaussian_bump_2d(n_rows, n_cols, 1.0, 0.5);

        for _ in 0..300 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            if !dt.is_finite() {
                break;
            }
            forward_euler_step(&mut states, &mesh, Boundaries2D::TRANSMISSIVE, dt);
        }
        for s in &states {
            assert!(s.h.is_finite(), "h became non-finite: {}", s.h);
            assert!(s.h >= 0.0, "h went negative: {}", s.h);
            assert!(s.hu.is_finite(), "hu became non-finite: {}", s.hu);
            assert!(s.hv.is_finite(), "hv became non-finite: {}", s.hv);
        }
    }

    #[test]
    fn bump_propagates_outward_and_central_depth_decays() {
        let n_rows = 21;
        let n_cols = 21;
        let mesh = flat_mesh(n_rows, n_cols, 1.0, 1.0);
        let h_base = 1.0;
        let h_amp = 0.5;
        let mut states = gaussian_bump_2d(n_rows, n_cols, h_base, h_amp);
        let center = (n_rows / 2, n_cols / 2);
        let h_center_initial = states[center].h;

        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::TRANSMISSIVE, dt);
        }
        assert!(
            states[center].h < h_center_initial,
            "central depth did not decay: started {}, ended {}",
            h_center_initial,
            states[center].h
        );
        assert!(
            states[center].h > 0.0,
            "central depth went non-positive: {}",
            states[center].h
        );
    }

    #[test]
    fn bump_propagates_isotropically_on_uniform_mesh() {
        // Round Gaussian bump on a square flat mesh with dx=dy. Depth
        // profile must remain symmetric under (i ↔ n-1-i) and
        // (j ↔ n-1-j) up to roundoff. Catches asymmetric x/y wiring
        // bugs in the FV update or in the face-flux assembly.
        let n = 21;
        let mesh = flat_mesh(n, n, 1.0, 1.0);
        let mut states = gaussian_bump_2d(n, n, 1.0, 0.5);
        for _ in 0..50 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::TRANSMISSIVE, dt);
        }
        for i in 0..n {
            for j in 0..n {
                assert_relative_eq!(states[(i, j)].h, states[(n - 1 - i, j)].h, epsilon = 1e-10);
                assert_relative_eq!(states[(i, j)].h, states[(i, n - 1 - j)].h, epsilon = 1e-10);
                assert_relative_eq!(states[(i, j)].h, states[(j, i)].h, epsilon = 1e-10);
            }
        }
    }

    #[test]
    #[should_panic(expected = "must match mesh")]
    fn mismatched_shape_panics() {
        let mesh = flat_mesh(4, 5, 1.0, 1.0);
        let mut states = Array2::from_elem((3, 3), Conserved2D::new(1.0, 0.0, 0.0));
        forward_euler_step(&mut states, &mesh, Boundaries2D::TRANSMISSIVE, 0.01);
    }
}
