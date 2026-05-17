//! Time-stepping for the 1D Saint-Venant solver.
//!
//! Explicit forward Euler in time, finite-volume in space, with the HLL
//! interface flux and the hydrostatic reconstruction of Audusse et al.
//! (2004) for well-balanced treatment of the bed-slope source. Manning
//! friction is applied separately as an operator-split fractional step
//! (see [`crate::source`]).

use crate::GRAVITY;
use crate::boundary::{Boundaries, Side, ghost_cell};
use crate::flux::Flux;
use crate::geometry::Channel1D;
use crate::riemann::hll_flux;
use crate::state::Conserved;

/// Maximum signal speed `|u| + c` across the state vector. Used to bound
/// the time step under the CFL condition. Returns 0 for an empty or all-dry
/// state.
pub fn max_wave_speed(states: &[Conserved]) -> f64 {
    states
        .iter()
        .map(|s| {
            let c = (GRAVITY * s.h.max(0.0)).sqrt();
            let u = if s.h > 0.0 { s.hu / s.h } else { 0.0 };
            u.abs() + c
        })
        .fold(0.0, f64::max)
}

/// CFL-bounded time step: `dt = cfl · dx / max_wave_speed`. Returns
/// `f64::INFINITY` when the domain is entirely dry, since no signal can
/// propagate; callers should clamp this against a problem-specific maximum.
///
/// `cfl` is typically 0.5 for an explicit FV solver with HLL.
pub fn cfl_time_step(states: &[Conserved], dx: f64, cfl: f64) -> f64 {
    let smax = max_wave_speed(states);
    if smax > 0.0 {
        cfl * dx / smax
    } else {
        f64::INFINITY
    }
}

/// Numerical fluxes at a single interface, in the well-balanced "two-sided"
/// formulation. `minus` is consumed by the cell on the LEFT of the face;
/// `plus` is consumed by the cell on the RIGHT. They share `mass` but
/// differ in `momentum` by a hydrostatic-pressure correction term.
#[derive(Debug, Clone, Copy)]
struct FaceFluxes {
    minus: Flux,
    plus: Flux,
}

/// Audusse (2004) well-balanced face flux. Reconstructs the depths at the
/// face using the higher of the two bed elevations as a virtual interface,
/// computes the HLL flux on those reconstructed states, and adds the
/// (`g/2`)(`h² − h*²`) pressure correction on each side.
///
/// On a flat bed (`z_left == z_right`) `h*_L = h_left` and `h*_R = h_right`,
/// the correction vanishes, and the flux degenerates to the plain HLL flux
/// of the original states — i.e. this generalises the flat-bed update.
fn well_balanced_face(left: Conserved, z_left: f64, right: Conserved, z_right: f64) -> FaceFluxes {
    let z_max = z_left.max(z_right);
    let h_star_left = (left.h + z_left - z_max).max(0.0);
    let h_star_right = (right.h + z_right - z_max).max(0.0);

    // Reconstructed states carry the original velocities of each side.
    let u_left = if left.h > 0.0 { left.hu / left.h } else { 0.0 };
    let u_right = if right.h > 0.0 {
        right.hu / right.h
    } else {
        0.0
    };

    let f_hll = hll_flux(
        Conserved::new(h_star_left, h_star_left * u_left),
        Conserved::new(h_star_right, h_star_right * u_right),
    );

    let half_g = 0.5 * GRAVITY;
    FaceFluxes {
        minus: Flux {
            mass: f_hll.mass,
            momentum: f_hll.momentum + half_g * (left.h * left.h - h_star_left * h_star_left),
        },
        plus: Flux {
            mass: f_hll.mass,
            momentum: f_hll.momentum + half_g * (right.h * right.h - h_star_right * h_star_right),
        },
    }
}

/// One forward-Euler update of the FV solution. Modifies `states` in place.
///
/// Includes the well-balanced bed-slope source via hydrostatic
/// reconstruction. Friction is **not** included here — call
/// [`crate::source::manning_friction_step`] before or after this step.
///
/// Panics if `states.len() != channel.n_cells()`. The caller is responsible
/// for keeping `dt` below the CFL bound (see [`cfl_time_step`]).
pub fn forward_euler_step(states: &mut [Conserved], channel: &Channel1D, bcs: Boundaries, dt: f64) {
    assert_eq!(
        states.len(),
        channel.n_cells(),
        "states.len() ({}) must match channel.n_cells() ({})",
        states.len(),
        channel.n_cells()
    );

    let n = states.len();
    if n == 0 {
        return;
    }

    // n+1 face fluxes: face 0 is the left boundary, face n is the right.
    // Ghost cell state and bed come from the boundary kind via `ghost_cell`:
    // computational BCs (Transmissive, Wall) leave the bed flat across the
    // boundary; physical BCs (Discharge, Depth) extend it linearly so the
    // boundary face carries the same bed jump as interior faces.
    let mut faces: Vec<FaceFluxes> = Vec::with_capacity(n + 1);

    let (ghost_left, z_ghost_left) = ghost_cell(channel, states[0], bcs.left, Side::Left);
    faces.push(well_balanced_face(
        ghost_left,
        z_ghost_left,
        states[0],
        channel.bed[0],
    ));

    for i in 0..n.saturating_sub(1) {
        faces.push(well_balanced_face(
            states[i],
            channel.bed[i],
            states[i + 1],
            channel.bed[i + 1],
        ));
    }

    let (ghost_right, z_ghost_right) = ghost_cell(channel, states[n - 1], bcs.right, Side::Right);
    faces.push(well_balanced_face(
        states[n - 1],
        channel.bed[n - 1],
        ghost_right,
        z_ghost_right,
    ));

    // FV update: U_i^{n+1} = U_i^n - (dt/dx)(F^-_{i+1/2} - F^+_{i-1/2}).
    // Cell `i` consumes the `minus` of its right face and the `plus` of its
    // left face. Mass is identical in both (no bed correction in continuity);
    // momentum carries the hydrostatic balance term.
    let dt_dx = dt / channel.dx;
    for i in 0..n {
        let f_right = faces[i + 1].minus;
        let f_left = faces[i].plus;
        states[i].h -= dt_dx * (f_right.mass - f_left.mass);
        states[i].hu -= dt_dx * (f_right.momentum - f_left.momentum);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array1;

    fn flat_channel(n: usize, dx: f64) -> Channel1D {
        Channel1D::new(Array1::zeros(n), dx, 0.0)
    }

    fn sloped_channel(n: usize, dx: f64, slope: f64) -> Channel1D {
        let bed = Array1::from_iter((0..n).map(|i| -(i as f64) * dx * slope));
        Channel1D::new(bed, dx, 0.0)
    }

    fn bumpy_channel(n: usize, dx: f64, bump_amp: f64) -> Channel1D {
        let center = (n as f64) / 2.0;
        let width_sq = (n as f64).powi(2) / 25.0;
        let bed = Array1::from_iter((0..n).map(|i| {
            let x = i as f64;
            bump_amp * (-((x - center).powi(2) / width_sq)).exp()
        }));
        Channel1D::new(bed, dx, 0.0)
    }

    fn lake_at_rest_on(channel: &Channel1D, eta: f64) -> Vec<Conserved> {
        channel
            .bed
            .iter()
            .map(|&z| Conserved::new(eta - z, 0.0))
            .collect()
    }

    fn total_mass(states: &[Conserved], dx: f64) -> f64 {
        states.iter().map(|s| s.h * dx).sum()
    }

    fn gaussian_bump(n: usize, h_base: f64, h_amp: f64) -> Vec<Conserved> {
        let center = (n as f64) / 2.0;
        let width_sq = (n as f64).powi(2) / 25.0;
        (0..n)
            .map(|i| {
                let x = i as f64;
                let h = h_base + h_amp * (-((x - center).powi(2) / width_sq)).exp();
                Conserved::new(h, 0.0)
            })
            .collect()
    }

    #[test]
    fn max_wave_speed_is_zero_for_dry_domain() {
        let states = [Conserved::DRY; 5];
        assert_eq!(max_wave_speed(&states), 0.0);
    }

    #[test]
    fn max_wave_speed_matches_textbook_formula() {
        let states = [Conserved::new(1.0, 2.0)];
        let expected = 2.0 + (GRAVITY * 1.0).sqrt();
        assert_relative_eq!(max_wave_speed(&states), expected, epsilon = 1e-12);
    }

    #[test]
    fn cfl_time_step_dry_domain_returns_infinity() {
        let states = [Conserved::DRY; 5];
        assert_eq!(cfl_time_step(&states, 1.0, 0.5), f64::INFINITY);
    }

    #[test]
    fn lake_at_rest_on_flat_bed_is_preserved_exactly() {
        let n = 20;
        let dx = 1.0;
        let h0 = 2.0;
        let channel = flat_channel(n, dx);
        let mut states: Vec<Conserved> = vec![Conserved::new(h0, 0.0); n];

        let dt = 0.01;
        for _ in 0..100 {
            forward_euler_step(&mut states, &channel, Boundaries::TRANSMISSIVE, dt);
        }
        for s in &states {
            assert_relative_eq!(s.h, h0, epsilon = 1e-12);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn lake_at_rest_on_sloped_bed_is_preserved() {
        // η = h + z = const, u = 0 on a linearly descending bed. With
        // hydrostatic reconstruction the bed-slope source exactly cancels
        // the hydrostatic-pressure flux, so the lake stays at rest.
        let n = 30;
        let dx = 1.0;
        let slope = 0.05;
        let eta = 5.0;
        let channel = sloped_channel(n, dx, slope);
        let initial = lake_at_rest_on(&channel, eta);
        let mut states = initial.clone();

        for _ in 0..200 {
            let dt = cfl_time_step(&states, dx, 0.4);
            forward_euler_step(&mut states, &channel, Boundaries::WALLS, dt);
        }
        for (i, s) in states.iter().enumerate() {
            assert_relative_eq!(s.h, initial[i].h, epsilon = 1e-10);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn lake_at_rest_on_bumpy_bed_is_preserved() {
        // η = const above a fully-submerged Gaussian hill. Same test as
        // above with a non-monotonic bed; catches asymmetric bugs in the
        // reconstruction.
        let n = 50;
        let dx = 1.0;
        let bump_amp = 1.5;
        let eta = 3.0; // safely above max(bed)
        let channel = bumpy_channel(n, dx, bump_amp);
        let initial = lake_at_rest_on(&channel, eta);
        let mut states = initial.clone();

        for _ in 0..200 {
            let dt = cfl_time_step(&states, dx, 0.4);
            forward_euler_step(&mut states, &channel, Boundaries::WALLS, dt);
        }
        for (i, s) in states.iter().enumerate() {
            assert_relative_eq!(s.h, initial[i].h, epsilon = 1e-10);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn mass_is_conserved_exactly_with_wall_boundaries() {
        // Re-asserted under the new well-balanced flux: the wall HLL
        // mass flux is still 0 by anti-symmetry, so global mass is
        // conserved to roundoff.
        let n = 50;
        let dx = 1.0;
        let cfl = 0.4;
        let channel = flat_channel(n, dx);
        let mut states = gaussian_bump(n, 1.0, 0.5);
        let m0 = total_mass(&states, dx);

        for _ in 0..200 {
            let dt = cfl_time_step(&states, dx, cfl);
            forward_euler_step(&mut states, &channel, Boundaries::WALLS, dt);
        }
        let m1 = total_mass(&states, dx);
        assert_relative_eq!(m0, m1, epsilon = 1e-10);
    }

    #[test]
    fn bump_remains_bounded_with_transmissive_bc() {
        let n = 50;
        let dx = 1.0;
        let cfl = 0.4;
        let channel = flat_channel(n, dx);
        let mut states = gaussian_bump(n, 1.0, 0.5);

        for _ in 0..300 {
            let dt = cfl_time_step(&states, dx, cfl);
            if !dt.is_finite() {
                break;
            }
            forward_euler_step(&mut states, &channel, Boundaries::TRANSMISSIVE, dt);
        }
        for s in &states {
            assert!(s.h.is_finite(), "h became non-finite: {}", s.h);
            assert!(s.h >= 0.0, "h went negative: {}", s.h);
            assert!(s.hu.is_finite(), "hu became non-finite: {}", s.hu);
        }
    }

    #[test]
    fn bump_propagates_outward_and_decays() {
        let n = 50;
        let dx = 1.0;
        let cfl = 0.4;
        let h_base = 1.0;
        let h_amp = 0.5;
        let channel = flat_channel(n, dx);
        let mut states = gaussian_bump(n, h_base, h_amp);
        let center = n / 2;
        let h_center_initial = states[center].h;

        for _ in 0..200 {
            let dt = cfl_time_step(&states, dx, cfl);
            forward_euler_step(&mut states, &channel, Boundaries::TRANSMISSIVE, dt);
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
    #[should_panic(expected = "must match channel.n_cells()")]
    fn mismatched_lengths_panic() {
        let channel = flat_channel(10, 1.0);
        let mut states: Vec<Conserved> = vec![Conserved::new(1.0, 0.0); 5];
        forward_euler_step(&mut states, &channel, Boundaries::TRANSMISSIVE, 0.01);
    }
}
