//! High-level simulation driver.
//!
//! Every consumer of the solver used to hand-roll the same ~50-line
//! time loop (CFL step → FV update → friction → bookkeeping), each
//! with its own subtle choices of integrator ordering and stopping
//! logic. `Simulation` owns that loop once: construct it from a mesh,
//! an initial state and a [`SimulationConfig`], then call
//! [`Simulation::run_until`] (or [`Simulation::step`] for manual
//! control, e.g. to mutate boundary conditions per step for a
//! hydrograph).
//!
//! The driver is `f64`-only by design: it is the operational
//! entry point. Calibration workflows that thread `Dual` through the
//! solver drive the loop manually (see `solver-1d/tests/ad_gradient.rs`
//! for the pattern) — a differentiable driver would freeze design
//! decisions (what is a parameter, what is a constant) that belong to
//! the calibration harness.
//!
//! Runtime failures (a degenerate time step, a step-budget blowout)
//! are typed [`SimError`]s, not panics: a long-running simulation
//! driven by an optimiser or a service must be able to observe the
//! failure and keep the process alive.

use ndarray::Array2;
use thiserror::Error;

use crate::boundary::Boundaries2D;
use crate::geometry::Mesh2D;
use crate::source::manning_friction_step;
use crate::state::Conserved2D;
use crate::update::{
    StepWorkspace2D, cfl_time_step_with_bcs, forward_euler_step_with, ssprk2_step_with,
};
use crate::H_DRY;

/// Time integrator for the hyperbolic part of the step.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Integrator {
    /// First-order forward Euler.
    ForwardEuler,
    /// Second-order SSP-RK2 (Heun). The default: it is what the MUSCL
    /// spatial reconstruction is designed to pair with.
    #[default]
    SspRk2,
}

/// Configuration of a [`Simulation`].
#[derive(Debug, Clone, Copy)]
pub struct SimulationConfig {
    /// CFL number for the adaptive time step. Default 0.4.
    pub cfl: f64,
    /// Time integrator. Default SSP-RK2.
    pub integrator: Integrator,
    /// Boundary conditions on the four sides. Default: walls.
    pub boundaries: Boundaries2D,
    /// Upper bound on the number of steps a single `run_until` may
    /// take before erroring out (guards against a collapsed dt turning
    /// a batch job into an infinite loop). Default 10 million.
    pub max_steps: usize,
    /// Cap on any single dt [s]. A fully dry domain has no interior
    /// signal and an unbounded CFL step; the cap keeps forcing terms
    /// (boundary inflow, rain applied by the caller between steps)
    /// well-resolved. Default 60 s.
    pub max_dt: f64,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            cfl: 0.4,
            integrator: Integrator::default(),
            boundaries: Boundaries2D::WALLS,
            max_steps: 10_000_000,
            max_dt: 60.0,
        }
    }
}

/// Errors surfaced by the simulation driver.
#[derive(Debug, Error)]
pub enum SimError {
    /// States shape does not match the mesh.
    #[error("states shape {states:?} must match mesh ({rows}, {cols})")]
    ShapeMismatch {
        /// Shape of the provided state array.
        states: (usize, usize),
        /// Mesh rows.
        rows: usize,
        /// Mesh cols.
        cols: usize,
    },
    /// Configuration value out of range.
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),
    /// The CFL machinery produced a non-positive or non-finite dt.
    #[error("degenerate time step dt = {dt} at t = {t} (step {step})")]
    DegenerateDt {
        /// Offending dt [s].
        dt: f64,
        /// Simulation time [s].
        t: f64,
        /// Step counter.
        step: usize,
    },
    /// `max_steps` exhausted before reaching the target time.
    #[error("step budget exhausted: {steps} steps, t = {t} of {target} s")]
    StepBudgetExhausted {
        /// Steps taken.
        steps: usize,
        /// Simulation time reached [s].
        t: f64,
        /// Requested target time [s].
        target: f64,
    },
}

/// Owns the mesh, the evolving state, the scratch workspace and the
/// simulation clock. See the module docs for the design rationale.
pub struct Simulation {
    mesh: Mesh2D,
    config: SimulationConfig,
    states: Array2<Conserved2D>,
    workspace: StepWorkspace2D<f64>,
    t: f64,
    steps: usize,
}

impl Simulation {
    /// Build a simulation from a mesh, an initial state and a config.
    pub fn new(
        mesh: Mesh2D,
        states: Array2<Conserved2D>,
        config: SimulationConfig,
    ) -> Result<Self, SimError> {
        let (rows, cols) = (mesh.n_rows(), mesh.n_cols());
        if states.dim() != (rows, cols) {
            return Err(SimError::ShapeMismatch {
                states: states.dim(),
                rows,
                cols,
            });
        }
        if !(config.cfl > 0.0 && config.cfl <= 1.0) {
            return Err(SimError::InvalidConfig(format!(
                "cfl must be in (0, 1], got {}",
                config.cfl
            )));
        }
        if !(config.max_dt > 0.0) {
            return Err(SimError::InvalidConfig(format!(
                "max_dt must be positive, got {}",
                config.max_dt
            )));
        }
        let workspace = StepWorkspace2D::for_mesh(&mesh);
        Ok(Self {
            mesh,
            config,
            states,
            workspace,
            t: 0.0,
            steps: 0,
        })
    }

    /// Current simulation time [s].
    pub fn time(&self) -> f64 {
        self.t
    }

    /// Steps taken since construction.
    pub fn steps(&self) -> usize {
        self.steps
    }

    /// Read access to the evolving state.
    pub fn states(&self) -> &Array2<Conserved2D> {
        &self.states
    }

    /// Mutable access to the state (initial-condition surgery, point
    /// sources or rain applied between steps by the caller).
    pub fn states_mut(&mut self) -> &mut Array2<Conserved2D> {
        &mut self.states
    }

    /// The mesh this simulation runs on.
    pub fn mesh(&self) -> &Mesh2D {
        &self.mesh
    }

    /// Replace the boundary conditions (hydrograph forcing: update the
    /// `Discharge`/`Depth` values as the event evolves).
    pub fn set_boundaries(&mut self, bcs: Boundaries2D) {
        self.config.boundaries = bcs;
    }

    /// Advance one step, bounded by `dt_cap` (in addition to the CFL
    /// bound and `config.max_dt`). Returns the dt actually taken.
    ///
    /// One step = hyperbolic update (Euler or SSP-RK2) followed by the
    /// operator-split semi-implicit Manning friction, the ordering
    /// used across the validation suite.
    pub fn step(&mut self, dt_cap: f64) -> Result<f64, SimError> {
        let bcs = self.config.boundaries;
        let dt = cfl_time_step_with_bcs(&self.states, &self.mesh, bcs, self.config.cfl)
            .min(self.config.max_dt)
            .min(dt_cap);
        if !(dt.is_finite() && dt > 0.0) {
            return Err(SimError::DegenerateDt {
                dt,
                t: self.t,
                step: self.steps,
            });
        }
        match self.config.integrator {
            Integrator::ForwardEuler => {
                forward_euler_step_with(&mut self.states, &self.mesh, bcs, dt, &mut self.workspace);
            }
            Integrator::SspRk2 => {
                ssprk2_step_with(&mut self.states, &self.mesh, bcs, dt, &mut self.workspace);
            }
        }
        manning_friction_step(&mut self.states, &self.mesh, dt, H_DRY);
        self.t += dt;
        self.steps += 1;
        Ok(dt)
    }

    /// Advance until the simulation clock reaches `t_target`.
    pub fn run_until(&mut self, t_target: f64) -> Result<(), SimError> {
        let mut budget = self.config.max_steps;
        while self.t < t_target {
            if budget == 0 {
                return Err(SimError::StepBudgetExhausted {
                    steps: self.steps,
                    t: self.t,
                    target: t_target,
                });
            }
            self.step(t_target - self.t)?;
            budget -= 1;
        }
        Ok(())
    }
}

/// Cheap plausibility check for a simulated depth field, independent
/// of any `NaN`/`inf` guard. A numerical blow-up on steep terrain can
/// be entirely finite (see
/// `docs/bug-report-2026-07-boundary-slope-instability.md` §4 in the
/// repo root) — depths of thousands of metres that are still normal
/// floating-point numbers, so a caller checking only for `NaN`/`inf`
/// sees a "valid" result.
///
/// Returns `true` when the maximum depth in `states` exceeds `factor`
/// times the mesh's own bed relief (`bed.max() - bed.min()`, floored
/// at 1 m so a near-flat mesh does not make the check degenerate): no
/// physically sane flood holds more water than some modest multiple of
/// the terrain's own vertical range. Meant for real-DEM windows with
/// genuine relief, not flat synthetic test meshes. `factor` is
/// caller-chosen; 2-5 is a reasonable default for flood applications —
/// generous enough to allow a closed depression to pond well above its
/// own rim without flagging every legitimately deep lake.
pub fn max_depth_exceeds_relief(states: &Array2<Conserved2D>, mesh: &Mesh2D, factor: f64) -> bool {
    let max_depth = states.iter().map(|s| s.h).fold(0.0_f64, f64::max);
    let bed_max = mesh.bed.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let bed_min = mesh.bed.iter().cloned().fold(f64::INFINITY, f64::min);
    let relief = (bed_max - bed_min).max(1.0);
    max_depth > factor * relief
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::update::cfl_time_step;
    use approx::assert_relative_eq;

    fn bump_setup(n: usize) -> (Mesh2D, Array2<Conserved2D>) {
        let mesh = Mesh2D::new(Array2::<f64>::zeros((n, n)), 1.0, 1.0, 0.03);
        let states = Array2::from_shape_fn((n, n), |(i, j)| {
            let x = j as f64 - n as f64 / 2.0;
            let y = i as f64 - n as f64 / 2.0;
            Conserved2D::new(1.0 + 0.3 * (-(x * x + y * y) / 8.0).exp(), 0.0, 0.0)
        });
        (mesh, states)
    }

    #[test]
    fn driver_matches_hand_rolled_loop_bitwise() {
        // The driver must be a refactoring of the canonical loop, not
        // a reinterpretation: same CFL, same integrator, same friction
        // ordering → bit-identical states.
        let (mesh, init) = bump_setup(16);
        let t_end = 2.0;

        let mut sim = Simulation::new(mesh.clone(), init.clone(), SimulationConfig::default())
            .expect("valid setup");
        sim.run_until(t_end).expect("run");

        let mut states = init;
        let mut ws = StepWorkspace2D::for_mesh(&mesh);
        let mut t = 0.0;
        while t < t_end {
            let dt = cfl_time_step_with_bcs(&states, &mesh, Boundaries2D::WALLS, 0.4)
                .min(60.0)
                .min(t_end - t);
            ssprk2_step_with(&mut states, &mesh, Boundaries2D::WALLS, dt, &mut ws);
            manning_friction_step(&mut states, &mesh, dt, H_DRY);
            t += dt;
        }

        assert_eq!(sim.time(), t);
        for ((i, j), a) in sim.states().indexed_iter() {
            let b = states[(i, j)];
            assert!(
                a.h == b.h && a.hu == b.hu && a.hv == b.hv,
                "diverged at ({i},{j}): {a:?} vs {b:?}"
            );
        }
    }

    #[test]
    fn mass_is_conserved_through_the_driver() {
        let (mesh, init) = bump_setup(20);
        let m0: f64 = init.iter().map(|s| s.h).sum();
        let mut sim =
            Simulation::new(mesh, init, SimulationConfig::default()).expect("valid setup");
        sim.run_until(5.0).expect("run");
        let m1: f64 = sim.states().iter().map(|s| s.h).sum();
        assert_relative_eq!(m0, m1, epsilon = 1e-10);
        assert!(sim.steps() > 0);
    }

    #[test]
    fn shape_mismatch_is_a_typed_error() {
        let (mesh, _) = bump_setup(8);
        let bad = Array2::from_elem((4, 4), Conserved2D::DRY);
        match Simulation::new(mesh, bad, SimulationConfig::default()) {
            Err(SimError::ShapeMismatch { states, rows, cols }) => {
                assert_eq!(states, (4, 4));
                assert_eq!((rows, cols), (8, 8));
            }
            other => panic!("expected ShapeMismatch, got {:?}", other.map(|_| ())),
        }
    }

    #[test]
    fn invalid_cfl_is_rejected() {
        let (mesh, init) = bump_setup(8);
        let config = SimulationConfig {
            cfl: 1.5,
            ..Default::default()
        };
        assert!(matches!(
            Simulation::new(mesh, init, config),
            Err(SimError::InvalidConfig(_))
        ));
    }

    #[test]
    fn step_budget_exhaustion_is_a_typed_error() {
        let (mesh, init) = bump_setup(8);
        let config = SimulationConfig {
            max_steps: 3,
            ..Default::default()
        };
        let mut sim = Simulation::new(mesh, init, config).expect("valid setup");
        match sim.run_until(1.0e9) {
            Err(SimError::StepBudgetExhausted { steps, .. }) => assert_eq!(steps, 3),
            other => panic!("expected StepBudgetExhausted, got {:?}", other.map(|_| ())),
        }
    }

    #[test]
    fn max_depth_exceeds_relief_flags_a_finite_blow_up() {
        // 97 m of relief (matches the Curacautín reproducer in the bug
        // report), a plausible 5 m flood, then an implausible one.
        let bed = Array2::from_shape_fn((4, 4), |(i, _j)| i as f64 * 32.0); // 0..96 m
        let mesh = Mesh2D::new(bed, 30.0, 30.0, 0.035);
        let plausible = Array2::from_elem((4, 4), Conserved2D::new(5.0, 0.0, 0.0));
        let implausible = Array2::from_elem((4, 4), Conserved2D::new(3000.0, 0.0, 0.0));
        assert!(!max_depth_exceeds_relief(&plausible, &mesh, 3.0));
        assert!(max_depth_exceeds_relief(&implausible, &mesh, 3.0));
    }

    #[test]
    fn max_depth_exceeds_relief_floors_a_near_flat_mesh() {
        // Flat bed: relief floors at 1 m instead of dividing by ~0.
        let mesh = Mesh2D::new(Array2::<f64>::zeros((3, 3)), 1.0, 1.0, 0.03);
        let modest = Array2::from_elem((3, 3), Conserved2D::new(2.0, 0.0, 0.0));
        let extreme = Array2::from_elem((3, 3), Conserved2D::new(50.0, 0.0, 0.0));
        assert!(!max_depth_exceeds_relief(&modest, &mesh, 3.0));
        assert!(max_depth_exceeds_relief(&extreme, &mesh, 3.0));
    }
}
