//! Time-stepping for the 1D Saint-Venant solver.
//!
//! Explicit forward Euler in time, finite-volume in space, with the HLL
//! interface flux. Source terms (bed slope, Manning friction) are deferred
//! to a follow-up commit; this loop assumes a flat, frictionless prismatic
//! channel.

use crate::GRAVITY;
use crate::boundary::{Boundaries, ghost_state};
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

/// One forward-Euler update of the FV solution. Modifies `states` in place.
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

    // n+1 interface fluxes: 1 left-boundary face + (n-1) internal faces +
    // 1 right-boundary face.
    let mut fluxes: Vec<Flux> = Vec::with_capacity(n + 1);

    let ghost_left = ghost_state(states[0], bcs.left);
    fluxes.push(hll_flux(ghost_left, states[0]));

    for i in 0..n.saturating_sub(1) {
        fluxes.push(hll_flux(states[i], states[i + 1]));
    }

    let ghost_right = ghost_state(states[n - 1], bcs.right);
    fluxes.push(hll_flux(states[n - 1], ghost_right));

    // FV update: U_i^{n+1} = U_i^n - (dt/dx)(F_{i+1/2} - F_{i-1/2}).
    let dt_dx = dt / channel.dx;
    for i in 0..n {
        let f_left = fluxes[i];
        let f_right = fluxes[i + 1];
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
        // h = 1, u = 2 → max signal speed = |u| + sqrt(g h) = 2 + sqrt(9.81).
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
        // Uniform depth, zero velocity, flat bed, transmissive BC: each
        // interface flux equals F(U) (consistency of HLL), so all flux
        // differences are exactly zero. Mass and momentum updates subtract
        // 0.0, an exact IEEE operation. After N steps the state is bit-exact.
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
    fn mass_is_conserved_exactly_with_wall_boundaries() {
        // Gaussian bump + walls: the wall HLL flux has F*.mass = 0 by
        // anti-symmetry, so the telescoping sum over cells gives ΔM = 0
        // per step (to roundoff).
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
        // Sanity: gaussian bump on flat bed, transmissive BC, many steps.
        // Depths must remain finite and non-negative; momentum must stay
        // finite. Mass leaves the domain so total mass is NOT conserved.
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
        // After enough time with transmissive BC, the central depth
        // perturbation must decrease (waves carry mass away from center).
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
        // Central depth must have dropped from its initial peak.
        assert!(
            states[center].h < h_center_initial,
            "central depth did not decay: started {}, ended {}",
            h_center_initial,
            states[center].h
        );
        // And must still be above the base (not yet fully drained).
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
