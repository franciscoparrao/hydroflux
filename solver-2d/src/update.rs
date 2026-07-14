//! Time-stepping for the 2D Saint-Venant solver.
//!
//! Explicit forward Euler in time, finite-volume in space, with the
//! HLLC interface flux and the hydrostatic reconstruction of Audusse
//! et al. (2004) extended per face direction. Manning friction is
//! applied separately as an operator-split fractional step (see
//! `crate::source`).
//!
//! The kernel is generic over a [`Real`] scalar so the same code
//! evaluates `f64` for production runs and `Dual` for forward-mode
//! AD: `Conserved2DG<T>` cells flow through a `Mesh2DG<T>` mesh and
//! produce gradient-carrying outputs without a parallel adjoint
//! implementation. Branching on wet/dry uses `.value()` to keep
//! control flow real-valued; arithmetic on each branch propagates
//! derivatives normally.
//!
//! # MUSCL slope-limited reconstruction
//!
//! At each interior face the cell states are linearly extrapolated to
//! the face midpoint using minmod-limited slopes (van Leer 1979). The
//! reconstruction is performed on the **primitive** vector `(η, u, v)`
//! where `η = h + z` is the water-surface elevation and `u, v` are
//! velocity components — *not* on the conserved `(h, hu, hv)` directly,
//! and *not* on `(η, hu, hv)` either. Reconstructing velocities rather
//! than momenta is what gives consistent face states on a non-flat
//! bed: when the analytical solution has `u` uniform but `h` varying
//! (e.g. Manning normal flow), `(η, u, v)`-MUSCL gives equal velocity
//! on both sides of every face, and the HLLC sees a consistent state.
//! Reconstructing `(η, hu, hv)` would shift the velocity across the
//! face because `h` differs between the two sides of an Audusse
//! reconstruction, producing a 1% steady-state drift in MacDonald
//! uniform flow. See Liang & Marche (2009), Kurganov & Petrova (2007),
//! and Bouchut (2004) for the same choice in the SWE literature.
//!
//! The bed is reconstructed linearly at each face as the midpoint of
//! the two adjacent cell-centered bed elevations (Liang & Marche
//! 2009). Both sides of the face see the SAME `z_face`, which makes
//! the Audusse hydrostatic correction vanish (the Audusse term
//! `(g/2)(h² − h*²)` is exactly zero when `z_L = z_R = z_face`). The
//! bed-slope source then moves out of the face flux and into an
//! explicit cell-centered term `S = −g · h · ∇z` evaluated with
//! central differences on the bed (one-sided at the boundaries). The
//! flux divergence and the explicit source cancel exactly for the
//! lake-at-rest configuration (C-property), and the net contribution
//! reduces to the analytical bed-slope force for non-trivial flows
//! over a smooth bed — eliminating the `O(dx · S₀ / h_n)` steady-state
//! bias that earlier η-MUSCL iterations had without bed
//! reconstruction.
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

use hydroflux_autograd::Real;

use crate::boundary::{Boundaries2D, Side, ghost_cell};
use crate::flux::{FluxXG, FluxYG};
use crate::geometry::Mesh2DG;
use crate::parallel::{MaybeSendSync, zip_for_each};
use crate::riemann::{hllc_flux_x, hllc_flux_y};
use crate::state::Conserved2DG;
use crate::{GRAVITY, H_DRY, H_VEL};
use ndarray::{Array2, Zip};

/// Ratio threshold for the explicit bed-slope source's bounded-slope
/// safety valve (see the comment at its use site in
/// [`forward_euler_step_with`], and
/// docs/bug-report-2026-07-boundary-slope-instability.md §7). An
/// ordinary (non-shoreline) face's `h_face` exceeding `h_old` by more
/// than this multiple only happens when the neighbouring bed sits far
/// enough below this cell's own bed that the quadratic source form
/// scales with the bed jump rather than the depth.
///
/// A ratio of 5 looked generous against the flat/uniformly-sloped
/// fixtures (worst case ~1.03), but `lake_at_rest_with_emerged_island_
/// is_preserved` broke that intuition: on a smooth Gaussian shoreline,
/// `h -> 0` continuously approaching the shore while the local bed
/// gradient does NOT vanish, so an ordinary (non-shoreline-rule) wet
/// cell one ring outside the actual wet/dry boundary can legitimately
/// see a ratio around 6 -- "thin water near a real shoreline", not
/// "thin film on steep terrain far from any shoreline", and the two
/// are indistinguishable by ratio alone at that scale. Calibrated
/// instead well above the worst ratio empirically observed across the
/// full lake-at-rest suite (~6) and still ~100x below the pathological
/// ratio in the bug report's reproducer (> 10^4): the gap between
/// "smooth benchmark bed, however close to a shoreline" and "thin film
/// on steep terrain" is still several orders of magnitude, just not as
/// many as first assumed.
const STEEP_SOURCE_RATIO: f64 = 500.0;

/// Maximum signal speeds `(s_x, s_y)` across the state field, where
/// `s_x = max(|u| + c)` and `s_y = max(|v| + c)`. Returns `(0, 0)` for
/// an empty or all-dry state. Wave speeds are reduced to `f64` (via
/// `Real::value()`) because they feed CFL bookkeeping, which is a
/// scalar timing decision rather than a differentiated quantity.
///
/// Dry cells (`h ≤ H_DRY`) contribute zero to both maxima — the CFL
/// must not be tightened by spurious `hu / h` blow-ups in essentially
/// dry cells.
pub fn max_wave_speeds<T: Real>(states: &Array2<Conserved2DG<T>>) -> (f64, f64) {
    let mut s_x = 0.0_f64;
    let mut s_y = 0.0_f64;
    for s in states {
        if s.h.value() <= H_DRY {
            continue;
        }
        let h = s.h.value();
        let c = (GRAVITY * h).sqrt();
        // Velocity cutoff: below H_VEL the post-update floor keeps the
        // momentum at zero, so hu/h is either zero or a stale residual
        // (e.g. a user-supplied initial condition) that would collapse
        // dt; either way the film carries no meaningful signal beyond
        // its gravity-wave celerity.
        let (u, v) = if h > H_VEL {
            (s.hu.value() / h, s.hv.value() / h)
        } else {
            (0.0, 0.0)
        };
        s_x = s_x.max(u.abs() + c);
        s_y = s_y.max(v.abs() + c);
    }
    (s_x, s_y)
}

/// CFL-bounded time step `dt = cfl / (s_x/dx + s_y/dy)` from the
/// interior state. Returns `f64::INFINITY` when the domain is entirely
/// dry (no signal can propagate); callers should clamp against a
/// problem-specific maximum.
///
/// Use [`cfl_time_step_with_bcs`] when starting from a fully dry
/// domain with a wet inflow boundary — that variant peeks at the
/// boundary ghosts so the first time step does not run away.
///
/// `cfl` is typically 0.4–0.5 for an explicit FV solver with HLLC in 2D.
pub fn cfl_time_step<T: Real>(states: &Array2<Conserved2DG<T>>, mesh: &Mesh2DG<T>, cfl: f64) -> f64 {
    let (s_x, s_y) = max_wave_speeds(states);
    let denom = s_x / mesh.dx + s_y / mesh.dy;
    if denom > 0.0 {
        cfl / denom
    } else {
        f64::INFINITY
    }
}

/// CFL-bounded time step that also considers wave activity carried by
/// the boundary ghost cells. Use this when the interior may be dry
/// while a [`Boundary::Discharge`] BC is injecting flow — in that
/// case [`cfl_time_step`] returns `INFINITY` (no interior signal) and
/// the solver would attempt a step long enough to inject unphysical
/// mass from the wet ghost.
///
/// The implementation samples the four boundary ghost cells along each
/// side and folds their wave speeds into the interior maxima. When the
/// interior is already active, the ghost contribution is dominated by
/// it and this function returns essentially the same value as
/// [`cfl_time_step`].
pub fn cfl_time_step_with_bcs<T: Real>(
    states: &Array2<Conserved2DG<T>>,
    mesh: &Mesh2DG<T>,
    bcs: Boundaries2D,
    cfl: f64,
) -> f64 {
    let (mut s_x, mut s_y) = max_wave_speeds(states);
    let n_rows = mesh.n_rows();
    let n_cols = mesh.n_cols();
    let mut update_from_ghost = |ghost: Conserved2DG<T>| {
        let h = ghost.h.value();
        if h <= H_DRY {
            return;
        }
        let c = (GRAVITY * h).sqrt();
        let u = ghost.hu.value() / h;
        let v = ghost.hv.value() / h;
        s_x = s_x.max(u.abs() + c);
        s_y = s_y.max(v.abs() + c);
    };
    for i in 0..n_rows {
        let (g_w, _) = ghost_cell(mesh, states[(i, 0)], bcs.west, Side::West, i);
        update_from_ghost(g_w);
        let (g_e, _) = ghost_cell(mesh, states[(i, n_cols - 1)], bcs.east, Side::East, i);
        update_from_ghost(g_e);
    }
    for j in 0..n_cols {
        let (g_n, _) = ghost_cell(mesh, states[(0, j)], bcs.north, Side::North, j);
        update_from_ghost(g_n);
        let (g_s, _) = ghost_cell(mesh, states[(n_rows - 1, j)], bcs.south, Side::South, j);
        update_from_ghost(g_s);
    }
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
struct FaceFluxXG<T> {
    minus: FluxXG<T>,
    plus: FluxXG<T>,
}

impl<T: Real> FaceFluxXG<T> {
    /// Zero face flux on both sides. Used as the cheap return path for
    /// dry-dry interior faces, where the Audusse HLLC flux is exactly
    /// zero (both `h*_L = h*_R = 0` ⇒ HLLC mass = 0, and both
    /// `(g/2)(h² − h*²) = 0` ⇒ corrections vanish identically).
    fn zero() -> Self {
        Self {
            minus: FluxXG::zero(),
            plus: FluxXG::zero(),
        }
    }
}

/// Numerical fluxes at a single `y`-face. `minus` is consumed by the
/// cell ABOVE the face (lower row index); `plus` by the cell BELOW.
/// They share `mass` and `x_momentum` but differ in `y_momentum`.
#[derive(Debug, Clone, Copy)]
struct FaceFluxYG<T> {
    minus: FluxYG<T>,
    plus: FluxYG<T>,
}

impl<T: Real> FaceFluxYG<T> {
    fn zero() -> Self {
        Self {
            minus: FluxYG::zero(),
            plus: FluxYG::zero(),
        }
    }
}

/// Audusse well-balanced HLLC flux on an `x`-face. See module docs for
/// the reconstruction. On a flat bed (`z_left == z_right`) the
/// correction vanishes and the flux is the plain HLLC flux of the
/// original states.
fn well_balanced_x_face<T: Real>(
    left: Conserved2DG<T>,
    z_left: T,
    right: Conserved2DG<T>,
    z_right: T,
) -> FaceFluxXG<T> {
    let z_max = z_left.max(z_right);
    let h_star_left = (left.h + z_left - z_max).max(T::zero());
    let h_star_right = (right.h + z_right - z_max).max(T::zero());

    // Velocity extraction under the H_VEL cutoff: a film's hu/h is
    // either exactly zero (floor invariant) or a stale residual that
    // must not enter the Riemann problem. Also removes the former
    // `> 0.0` guard, whose window (0, H_DRY] could overflow on
    // caller-supplied states with tiny h and non-zero momentum.
    let (u_left, v_left) = if left.h.value() > H_VEL {
        (left.hu / left.h, left.hv / left.h)
    } else {
        (T::zero(), T::zero())
    };
    let (u_right, v_right) = if right.h.value() > H_VEL {
        (right.hu / right.h, right.hv / right.h)
    } else {
        (T::zero(), T::zero())
    };

    let recon_left = Conserved2DG::new_generic(h_star_left, h_star_left * u_left, h_star_left * v_left);
    let recon_right = Conserved2DG::new_generic(
        h_star_right,
        h_star_right * u_right,
        h_star_right * v_right,
    );

    let f = hllc_flux_x(recon_left, recon_right);
    let corr_left = (left.h * left.h - h_star_left * h_star_left) * (0.5 * GRAVITY);
    let corr_right = (right.h * right.h - h_star_right * h_star_right) * (0.5 * GRAVITY);
    FaceFluxXG {
        minus: FluxXG {
            mass: f.mass,
            x_momentum: f.x_momentum + corr_left,
            y_momentum: f.y_momentum,
        },
        plus: FluxXG {
            mass: f.mass,
            x_momentum: f.x_momentum + corr_right,
            y_momentum: f.y_momentum,
        },
    }
}

/// Audusse well-balanced HLLC flux on a `y`-face. `left` (alias top,
/// lower row index) is at `z_left`; `right` (bottom, higher row index)
/// is at `z_right`.
fn well_balanced_y_face<T: Real>(
    left: Conserved2DG<T>,
    z_left: T,
    right: Conserved2DG<T>,
    z_right: T,
) -> FaceFluxYG<T> {
    let z_max = z_left.max(z_right);
    let h_star_left = (left.h + z_left - z_max).max(T::zero());
    let h_star_right = (right.h + z_right - z_max).max(T::zero());

    // See well_balanced_x_face for the H_VEL rationale.
    let (u_left, v_left) = if left.h.value() > H_VEL {
        (left.hu / left.h, left.hv / left.h)
    } else {
        (T::zero(), T::zero())
    };
    let (u_right, v_right) = if right.h.value() > H_VEL {
        (right.hu / right.h, right.hv / right.h)
    } else {
        (T::zero(), T::zero())
    };

    let recon_left = Conserved2DG::new_generic(h_star_left, h_star_left * u_left, h_star_left * v_left);
    let recon_right = Conserved2DG::new_generic(
        h_star_right,
        h_star_right * u_right,
        h_star_right * v_right,
    );

    let f = hllc_flux_y(recon_left, recon_right);
    let corr_left = (left.h * left.h - h_star_left * h_star_left) * (0.5 * GRAVITY);
    let corr_right = (right.h * right.h - h_star_right * h_star_right) * (0.5 * GRAVITY);
    FaceFluxYG {
        minus: FluxYG {
            mass: f.mass,
            x_momentum: f.x_momentum,
            y_momentum: f.y_momentum + corr_left,
        },
        plus: FluxYG {
            mass: f.mass,
            x_momentum: f.x_momentum,
            y_momentum: f.y_momentum + corr_right,
        },
    }
}

/// Minmod slope limiter (van Leer 1974, Roe 1986). Returns the
/// smaller-magnitude argument when both have the same sign; zero
/// otherwise. The simplest TVD limiter — most dissipative of the
/// common choices but unconditionally stable and trivial to verify.
///
/// The branch decision uses `.value()`; only one of the two operands
/// is selected as the result, so its derivative carries through.
fn minmod<T: Real>(a: T, b: T) -> T {
    if a.value() * b.value() <= 0.0 {
        T::zero()
    } else if a.value().abs() < b.value().abs() {
        a
    } else {
        b
    }
}

/// Slopes per cell for the **primitive** reconstruction vector
/// `(η, u, v)` in a single coordinate direction.
///
/// Reconstructing velocities `u, v` (rather than momenta `hu, hv`) is
/// essential for well-balancedness on flows with non-uniform depth
/// over a sloped bed: when the analytical solution has constant `u`
/// but `h` varies (e.g. MacDonald uniform Manning flow with `η` linear
/// in the bed slope direction), velocity-reconstruction gives the
/// same `u` on both sides of every face, so the HLLC sees consistent
/// face states. Momentum-reconstruction would produce a steady-state
/// drift of order `dx·S₀/h_n`.
#[derive(Debug, Clone, Copy, Default)]
struct CellSlopesG<T> {
    eta: T,
    u: T,
    v: T,
}

/// Cell-centered primitive values `(η, u, v)`, materialised once per
/// step. Films carry `u = v = 0`: velocity is not meaningful below the
/// `H_VEL` cutoff, and dividing residual momentum by a near-`H_DRY`
/// depth produces unphysical velocities (see [`crate::H_VEL`]).
#[derive(Debug, Clone, Copy, Default)]
struct PrimG<T> {
    eta: T,
    u: T,
    v: T,
}

/// Fill the per-cell primitive buffer. One pass replaces the ~10
/// redundant `hu/h`, `hv/h` divisions per cell per step that the
/// slope and reconstruction stencils used to pay by recomputing
/// primitives at every stencil visit.
fn fill_primitives<T: Real + MaybeSendSync>(
    states: &Array2<Conserved2DG<T>>,
    bed: &Array2<T>,
    out: &mut Array2<PrimG<T>>,
) {
    zip_for_each!(Zip::indexed(out), |(i, j), p| {
        let s = states[(i, j)];
        *p = if s.h.value() > H_VEL {
            PrimG {
                eta: s.h + bed[(i, j)],
                u: s.hu / s.h,
                v: s.hv / s.h,
            }
        } else {
            PrimG {
                eta: s.h + bed[(i, j)],
                u: T::zero(),
                v: T::zero(),
            }
        };
    });
}

/// Returns `true` if any cell within distance 2 of `(i, j)` along the
/// `x` axis is dry (depth ≤ `H_DRY`, per the pre-step snapshot). Used
/// to drop slopes to zero in a 2-cell buffer around wet/dry fronts so
/// that MUSCL reconstruction from a wet cell cannot extrapolate over a
/// nearby dry cell and produce spurious overshoots at the front.
fn any_neighbor_dry_x(dry: &Array2<bool>, i: usize, j: usize, n_cols: usize) -> bool {
    let j_lo = j.saturating_sub(2);
    let j_hi = (j + 2).min(n_cols - 1);
    (j_lo..=j_hi).any(|jj| dry[(i, jj)])
}

/// Same as [`any_neighbor_dry_x`] but along the `y` axis.
fn any_neighbor_dry_y(dry: &Array2<bool>, i: usize, j: usize, n_rows: usize) -> bool {
    let i_lo = i.saturating_sub(2);
    let i_hi = (i + 2).min(n_rows - 1);
    (i_lo..=i_hi).any(|ii| dry[(ii, j)])
}

/// Compute minmod-limited primitive slopes per cell in the `x`
/// direction.
///
/// Interior cells (`1 ≤ j ≤ n_cols − 2`) use central minmod between
/// forward and backward differences. Boundary cells use the
/// available one-sided difference (forward for `j = 0`, backward for
/// `j = n_cols − 1`) — this lets a smooth steady-state flow on a
/// sloped bed carry its bed-slope through the first/last cells
/// consistently with the interior. If any cell within a 2-cell
/// buffer of the slope's stencil is dry the slope drops to zero
/// (first-order at wet/dry fronts) — this prevents a MUSCL
/// reconstruction from a wet cell from extrapolating over a nearby
/// dry cell and producing a spurious overshoot.
fn fill_slopes_x<T: Real + MaybeSendSync>(
    prim: &Array2<PrimG<T>>,
    dry: &Array2<bool>,
    mesh: &Mesh2DG<T>,
    out: &mut Array2<CellSlopesG<T>>,
) {
    let n_cols = mesh.n_cols();
    zip_for_each!(Zip::indexed(out), |(i, j), slope| {
        *slope = if n_cols < 2 || any_neighbor_dry_x(dry, i, j, n_cols) {
            CellSlopesG::default()
        } else if j == 0 {
            // Forward one-sided difference: slope = (right − centre) / dx.
            let c = prim[(i, 0)];
            let r = prim[(i, 1)];
            CellSlopesG {
                eta: (r.eta - c.eta) / mesh.dx,
                u: (r.u - c.u) / mesh.dx,
                v: (r.v - c.v) / mesh.dx,
            }
        } else if j + 1 == n_cols {
            // Backward one-sided difference.
            let l = prim[(i, j - 1)];
            let c = prim[(i, j)];
            CellSlopesG {
                eta: (c.eta - l.eta) / mesh.dx,
                u: (c.u - l.u) / mesh.dx,
                v: (c.v - l.v) / mesh.dx,
            }
        } else {
            // Interior: central minmod.
            let l = prim[(i, j - 1)];
            let c = prim[(i, j)];
            let r = prim[(i, j + 1)];
            CellSlopesG {
                eta: minmod((c.eta - l.eta) / mesh.dx, (r.eta - c.eta) / mesh.dx),
                u: minmod((c.u - l.u) / mesh.dx, (r.u - c.u) / mesh.dx),
                v: minmod((c.v - l.v) / mesh.dx, (r.v - c.v) / mesh.dx),
            }
        };
    });
}

/// Fill minmod-limited primitive slopes per cell in the `y`
/// direction. Boundary cells use the available one-sided difference.
fn fill_slopes_y<T: Real + MaybeSendSync>(
    prim: &Array2<PrimG<T>>,
    dry: &Array2<bool>,
    mesh: &Mesh2DG<T>,
    out: &mut Array2<CellSlopesG<T>>,
) {
    let n_rows = mesh.n_rows();
    zip_for_each!(Zip::indexed(out), |(i, j), slope| {
        *slope = if n_rows < 2 || any_neighbor_dry_y(dry, i, j, n_rows) {
            CellSlopesG::default()
        } else if i == 0 {
            let c = prim[(0, j)];
            let b = prim[(1, j)];
            CellSlopesG {
                eta: (b.eta - c.eta) / mesh.dy,
                u: (b.u - c.u) / mesh.dy,
                v: (b.v - c.v) / mesh.dy,
            }
        } else if i + 1 == n_rows {
            let t = prim[(i - 1, j)];
            let c = prim[(i, j)];
            CellSlopesG {
                eta: (c.eta - t.eta) / mesh.dy,
                u: (c.u - t.u) / mesh.dy,
                v: (c.v - t.v) / mesh.dy,
            }
        } else {
            let t = prim[(i - 1, j)];
            let c = prim[(i, j)];
            let b = prim[(i + 1, j)];
            CellSlopesG {
                eta: minmod((c.eta - t.eta) / mesh.dy, (b.eta - c.eta) / mesh.dy),
                u: minmod((c.u - t.u) / mesh.dy, (b.u - c.u) / mesh.dy),
                v: minmod((c.v - t.v) / mesh.dy, (b.v - c.v) / mesh.dy),
            }
        };
    });
}

/// Reconstruct the left/right cell states at an interior `x`-face
/// between cells `(i, j_left)` and `(i, j_right)` using MUSCL slopes,
/// **with shared bed-reconstruction at the face**.
///
/// The reconstruction is on primitives `(η, u, v)`. The bed at the
/// face is the linear-interpolation midpoint `z_face = ½(z_L + z_R)`,
/// SHARED between both sides. The depth at the face is therefore
/// `h = max(η_face − z_face, 0)` from both sides — no asymmetry due
/// to a piecewise-constant bed jump. This is the Liang & Marche
/// (2009) construction that eliminates the steady-state bias of
/// η-MUSCL on a sloped bed (the bias that motivated this iteration).
///
/// Returns `(recon_l, recon_r, z_face)`. The caller passes `z_face`
/// to BOTH `z_left` and `z_right` of [`well_balanced_x_face`], which
/// makes the Audusse correction vanish (the bed-slope source is
/// captured by the explicit cell-centered source term in
/// [`forward_euler_step`] instead).
#[allow(clippy::too_many_arguments)]
fn reconstruct_x_face_states<T: Real>(
    prim: &Array2<PrimG<T>>,
    slopes_x: &Array2<CellSlopesG<T>>,
    z_face: T,
    i: usize,
    j_left: usize,
    j_right: usize,
    dx: f64,
) -> (Conserved2DG<T>, Conserved2DG<T>) {
    let half_dx = 0.5 * dx;

    let l = prim[(i, j_left)];
    let eta_minus = l.eta + slopes_x[(i, j_left)].eta * half_dx;
    let u_minus = l.u + slopes_x[(i, j_left)].u * half_dx;
    let v_minus = l.v + slopes_x[(i, j_left)].v * half_dx;

    let r = prim[(i, j_right)];
    let eta_plus = r.eta - slopes_x[(i, j_right)].eta * half_dx;
    let u_plus = r.u - slopes_x[(i, j_right)].u * half_dx;
    let v_plus = r.v - slopes_x[(i, j_right)].v * half_dx;

    let h_minus = (eta_minus - z_face).max(T::zero());
    let h_plus = (eta_plus - z_face).max(T::zero());

    (
        Conserved2DG::new_generic(h_minus, h_minus * u_minus, h_minus * v_minus),
        Conserved2DG::new_generic(h_plus, h_plus * u_plus, h_plus * v_plus),
    )
}

/// Reconstruct the top/bottom cell states at an interior `y`-face
/// between cells `(i_top, j)` and `(i_bottom, j)` using MUSCL slopes,
/// with shared bed-reconstruction at the face. See
/// [`reconstruct_x_face_states`] for the rationale.
#[allow(clippy::too_many_arguments)]
fn reconstruct_y_face_states<T: Real>(
    prim: &Array2<PrimG<T>>,
    slopes_y: &Array2<CellSlopesG<T>>,
    z_face: T,
    i_top: usize,
    i_bottom: usize,
    j: usize,
    dy: f64,
) -> (Conserved2DG<T>, Conserved2DG<T>) {
    let half_dy = 0.5 * dy;

    let t = prim[(i_top, j)];
    let eta_minus = t.eta + slopes_y[(i_top, j)].eta * half_dy;
    let u_minus = t.u + slopes_y[(i_top, j)].u * half_dy;
    let v_minus = t.v + slopes_y[(i_top, j)].v * half_dy;

    let b = prim[(i_bottom, j)];
    let eta_plus = b.eta - slopes_y[(i_bottom, j)].eta * half_dy;
    let u_plus = b.u - slopes_y[(i_bottom, j)].u * half_dy;
    let v_plus = b.v - slopes_y[(i_bottom, j)].v * half_dy;

    let h_minus = (eta_minus - z_face).max(T::zero());
    let h_plus = (eta_plus - z_face).max(T::zero());

    (
        Conserved2DG::new_generic(h_minus, h_minus * u_minus, h_minus * v_minus),
        Conserved2DG::new_generic(h_plus, h_plus * u_plus, h_plus * v_plus),
    )
}

/// Face bed elevation at an interior face given the two adjacent cell
/// states and beds.
///
/// - **Both sides wet (or both dry)**: the linear-interpolation
///   midpoint `½(z_L + z_R)` — the Liang & Marche (2009) construction
///   that eliminates the steady-state bias of η-MUSCL on a sloped bed.
/// - **Exactly one side wet** (shoreline / wetting front): the Audusse
///   `max(z_L, z_R)`. The midpoint is wrong here in two ways. If the
///   dry cell's bed stands above the wet surface (`z_dry ≥ η_wet`) the
///   face is a physical wall, but the midpoint can sit below `η_wet`:
///   the dry side then reconstructs a spurious column of half the bed
///   jump, the α-rescaling zeroes the resulting flux, and the wet
///   cell's source keeps an unbalanced `g·h_face²/(2Δ)` term — a
///   permanent spurious acceleration toward the shore (broken
///   C-property at shorelines). With `max`, both reconstructed depths
///   vanish and flux and source balance exactly. If instead
///   `z_dry < η_wet` (legitimate wetting), `max` removes the same
///   spurious half-jump column from the dry side and lets the wet side
///   carry its full head into the front. Slopes are already zeroed in
///   a 2-cell buffer around wet/dry fronts, so the wet-wet midpoint
///   rationale (MUSCL steady bias) does not apply at these faces.
///
/// State-dependent by necessity — the C-property with wetting/drying
/// cannot be enforced by a purely static face bed. `l_dry`/`r_dry`
/// come from the pre-step wet/dry snapshot.
fn interior_z_face<T: Real>(l_dry: bool, z_l: T, r_dry: bool, z_r: T) -> T {
    if l_dry != r_dry {
        z_l.max(z_r)
    } else {
        (z_l + z_r) * 0.5
    }
}

/// Fill the array of face bed elevations `z_face_x[i, j]` for every
/// `x`-face (interior and boundary).
///
/// Interior `j ∈ 1..n_cols`: [`interior_z_face`] (midpoint on wet-wet,
/// `max` on wet/dry).
/// West boundary `j = 0`: `z_face = ½(z_ghost_west + z[i, 0])`.
/// East boundary `j = n_cols`: `z_face = ½(z[i, n_cols-1] + z_ghost_east)`.
/// Boundary faces keep the midpoint unconditionally: the ghost mirrors
/// the interior wetness for Wall/Transmissive (no wet/dry face arises),
/// and the Discharge-on-dry override manages its own ghost depth.
///
/// The cell-centered bed gradient used by the explicit source is then
/// `(z_face[i, j+1] − z_face[i, j]) / dx`, which is **consistent with
/// the flux divergence**: the same `z_face` values feed both the
/// pressure-flux term `g·h_face²/2` and the source. This consistency
/// is what makes lake-at-rest exact on a general bed — the cell-by-
/// cell cancellation between flux divergence and source only holds
/// when both use the same face beds.
fn fill_z_face_x<T: Real + MaybeSendSync>(
    states: &Array2<Conserved2DG<T>>,
    dry: &Array2<bool>,
    mesh: &Mesh2DG<T>,
    bcs: Boundaries2D,
    out: &mut Array2<T>,
) {
    let n_cols = mesh.n_cols();
    zip_for_each!(Zip::indexed(out), |(i, j), z| {
        *z = if j == 0 {
            let (_, z_g) = ghost_cell(mesh, states[(i, 0)], bcs.west, Side::West, i);
            (z_g + mesh.bed[(i, 0)]) * 0.5
        } else if j == n_cols {
            let (_, z_g) = ghost_cell(mesh, states[(i, n_cols - 1)], bcs.east, Side::East, i);
            (mesh.bed[(i, n_cols - 1)] + z_g) * 0.5
        } else {
            interior_z_face(
                dry[(i, j - 1)],
                mesh.bed[(i, j - 1)],
                dry[(i, j)],
                mesh.bed[(i, j)],
            )
        };
    });
}

/// Fill the array of face bed elevations `z_face_y[i, j]` for every
/// `y`-face. See [`fill_z_face_x`] for the rationale.
fn fill_z_face_y<T: Real + MaybeSendSync>(
    states: &Array2<Conserved2DG<T>>,
    dry: &Array2<bool>,
    mesh: &Mesh2DG<T>,
    bcs: Boundaries2D,
    out: &mut Array2<T>,
) {
    let n_rows = mesh.n_rows();
    zip_for_each!(Zip::indexed(out), |(i, j), z| {
        *z = if i == 0 {
            let (_, z_g) = ghost_cell(mesh, states[(0, j)], bcs.north, Side::North, j);
            (z_g + mesh.bed[(0, j)]) * 0.5
        } else if i == n_rows {
            let (_, z_g) = ghost_cell(mesh, states[(n_rows - 1, j)], bcs.south, Side::South, j);
            (mesh.bed[(n_rows - 1, j)] + z_g) * 0.5
        } else {
            interior_z_face(
                dry[(i - 1, j)],
                mesh.bed[(i - 1, j)],
                dry[(i, j)],
                mesh.bed[(i, j)],
            )
        };
    });
}

/// Reusable scratch buffers for the 2D time step.
///
/// [`forward_euler_step`] allocated seven fresh arrays per Euler step
/// (slopes ×2, face beds ×2, face fluxes ×2, α) — ≈ 170 MB of
/// alloc+write+read traffic per step on a 1024² mesh, which made the
/// step memory-bandwidth-bound. Constructing one `StepWorkspace2D` per
/// simulation and calling [`forward_euler_step_with`] /
/// [`ssprk2_step_with`] reuses the buffers across steps (and holds the
/// SSP-RK2 state snapshot, removing that clone too).
///
/// The buffers are an implementation detail: their contents are
/// overwritten at the start of every step and carry no state between
/// steps beyond capacity.
pub struct StepWorkspace2D<T> {
    prim: Array2<PrimG<T>>,
    slopes_x: Array2<CellSlopesG<T>>,
    slopes_y: Array2<CellSlopesG<T>>,
    z_face_x: Array2<T>,
    z_face_y: Array2<T>,
    faces_x: Array2<FaceFluxXG<T>>,
    faces_y: Array2<FaceFluxYG<T>>,
    alpha: Array2<T>,
    was_dry: Array2<bool>,
    u_n: Array2<Conserved2DG<T>>,
}

impl<T: Real> StepWorkspace2D<T> {
    /// Workspace for an `n_rows × n_cols` mesh.
    pub fn new(n_rows: usize, n_cols: usize) -> Self {
        Self {
            prim: Array2::from_elem((n_rows, n_cols), PrimG::default()),
            slopes_x: Array2::from_elem((n_rows, n_cols), CellSlopesG::default()),
            slopes_y: Array2::from_elem((n_rows, n_cols), CellSlopesG::default()),
            z_face_x: Array2::from_elem((n_rows, n_cols + 1), T::zero()),
            z_face_y: Array2::from_elem((n_rows + 1, n_cols), T::zero()),
            faces_x: Array2::from_elem((n_rows, n_cols + 1), FaceFluxXG::zero()),
            faces_y: Array2::from_elem((n_rows + 1, n_cols), FaceFluxYG::zero()),
            alpha: Array2::from_elem((n_rows, n_cols), T::one()),
            was_dry: Array2::from_elem((n_rows, n_cols), false),
            u_n: Array2::from_elem((n_rows, n_cols), Conserved2DG::dry()),
        }
    }

    /// Workspace sized for `mesh`.
    pub fn for_mesh(mesh: &Mesh2DG<T>) -> Self {
        Self::new(mesh.n_rows(), mesh.n_cols())
    }
}

/// One forward-Euler update of the 2D FV solution. Modifies `states` in
/// place.
///
/// Includes the well-balanced bed-slope source via per-face hydrostatic
/// reconstruction. Friction is **not** included — call
/// [`crate::source::manning_friction_step`] separately.
///
/// Allocates a fresh [`StepWorkspace2D`] on every call; loops that
/// step many times should construct one workspace and call
/// [`forward_euler_step_with`] instead.
///
/// Panics if the shape of `states` does not match `(mesh.n_rows(),
/// mesh.n_cols())`. The caller is responsible for keeping `dt` below
/// the CFL bound (see [`cfl_time_step`]).
pub fn forward_euler_step<T: Real + MaybeSendSync>(
    states: &mut Array2<Conserved2DG<T>>,
    mesh: &Mesh2DG<T>,
    bcs: Boundaries2D,
    dt: f64,
) {
    let mut ws = StepWorkspace2D::for_mesh(mesh);
    forward_euler_step_with(states, mesh, bcs, dt, &mut ws);
}

/// [`forward_euler_step`] over caller-owned scratch buffers — the
/// allocation-free hot path.
pub fn forward_euler_step_with<T: Real + MaybeSendSync>(
    states: &mut Array2<Conserved2DG<T>>,
    mesh: &Mesh2DG<T>,
    bcs: Boundaries2D,
    dt: f64,
    ws: &mut StepWorkspace2D<T>,
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
    // Split borrows: every pass below reads some buffers and writes
    // others; destructuring hands the borrow checker disjoint fields.
    let StepWorkspace2D {
        prim,
        slopes_x,
        slopes_y,
        z_face_x,
        z_face_y,
        faces_x,
        faces_y,
        alpha,
        was_dry,
        u_n: _,
    } = ws;
    assert_eq!(
        prim.shape(),
        [n_rows, n_cols],
        "workspace shape {:?} must match mesh ({}, {})",
        prim.shape(),
        n_rows,
        n_cols,
    );

    // Pre-step wet/dry snapshot. Feeds the slope front-buffer checks,
    // the shoreline face-bed rule, the dry-dry face short-circuit, the
    // face-scaling skip and the isolated-dry fast path — all of which
    // must see the SAME (pre-step) pattern for the step to be
    // consistent (the update loop mutates `states` in place).
    zip_for_each!(Zip::indexed(&mut *was_dry), |(i, j), d| {
        *d = states[(i, j)].h.value() <= H_DRY;
    });

    // Cell primitives (η, u, v), then MUSCL slopes per cell.
    fill_primitives(states, &mesh.bed, prim);
    fill_slopes_x(prim, was_dry, mesh, slopes_x);
    fill_slopes_y(prim, was_dry, mesh, slopes_y);

    // Face bed elevations: `z_face_x[i, j]` for every `x`-face and
    // `z_face_y[i, j]` for every `y`-face. These arrays feed BOTH the
    // face-flux pass (via Audusse-HLLC with `z_left = z_right =
    // z_face`) and the explicit cell-centered source term that
    // follows. The well-balancedness of the scheme hinges on this
    // single source of truth for the face beds. State-dependent at
    // shoreline faces (see `interior_z_face`), so filled per step.
    fill_z_face_x(states, was_dry, mesh, bcs, z_face_x);
    fill_z_face_y(states, was_dry, mesh, bcs, z_face_y);

    // Precompute all x-faces and y-faces. x-faces have shape
    // (n_rows, n_cols + 1); y-faces have shape (n_rows + 1, n_cols).
    // Interior dry-dry faces have identically zero flux: with both
    // adjacent cells at `h ≤ H_DRY` ⇒ `h*_L = h*_R = 0` ⇒ HLLC mass
    // & momentum vanish, and the Audusse correction `(g/2)(h² − h*²)`
    // is also exactly zero on each side. Short-circuiting these faces
    // skips the HLLC + MUSCL reconstruction work without changing
    // the numerics by even a roundoff (the path is mathematically
    // equivalent, not just "close to zero"). For Huasco-style
    // single-reach problems where ~97 % of cells are dry, the
    // skipped face count dominates the face-flux pass cost.
    //
    // Boundary faces (`j == 0` and `j == n_cols`) keep the full path:
    // the ghost cell may carry mass via `Discharge` or `Depth` BCs
    // independently of the interior — checking the ghost AND the
    // interior cell for dryness is cheap but adds a branch that
    // never fires for the Wall / Transmissive cases that dominate
    // our event simulations, so the win there is marginal.
    zip_for_each!(Zip::indexed(&mut *faces_x), |(i, j), face| {
        let z_face = z_face_x[(i, j)];
        *face = if j == 0 {
            let (g, _) = ghost_cell(mesh, states[(i, 0)], bcs.west, Side::West, i);
            well_balanced_x_face(g, z_face, states[(i, 0)], z_face)
        } else if j == n_cols {
            let (g, _) = ghost_cell(mesh, states[(i, n_cols - 1)], bcs.east, Side::East, i);
            well_balanced_x_face(states[(i, n_cols - 1)], z_face, g, z_face)
        } else if was_dry[(i, j - 1)] && was_dry[(i, j)] {
            FaceFluxXG::zero()
        } else {
            let (recon_l, recon_r) =
                reconstruct_x_face_states(prim, slopes_x, z_face, i, j - 1, j, mesh.dx);
            well_balanced_x_face(recon_l, z_face, recon_r, z_face)
        };
    });

    zip_for_each!(Zip::indexed(&mut *faces_y), |(i, j), face| {
        let z_face = z_face_y[(i, j)];
        *face = if i == 0 {
            let (g, _) = ghost_cell(mesh, states[(0, j)], bcs.north, Side::North, j);
            well_balanced_y_face(g, z_face, states[(0, j)], z_face)
        } else if i == n_rows {
            let (g, _) = ghost_cell(mesh, states[(n_rows - 1, j)], bcs.south, Side::South, j);
            well_balanced_y_face(states[(n_rows - 1, j)], z_face, g, z_face)
        } else if was_dry[(i - 1, j)] && was_dry[(i, j)] {
            FaceFluxYG::zero()
        } else {
            let (recon_t, recon_b) =
                reconstruct_y_face_states(prim, slopes_y, z_face, i - 1, i, j, mesh.dy);
            well_balanced_y_face(recon_t, z_face, recon_b, z_face)
        };
    });

    // FV update. For cell (i, j):
    //   right x-face is faces_x[(i, j+1)] — cell is on its LEFT side → .minus
    //   left  x-face is faces_x[(i, j)]   — cell is on its RIGHT side → .plus
    //   bottom y-face is faces_y[(i+1, j)] — cell is on its TOP side → .minus
    //   top    y-face is faces_y[(i, j)]   — cell is on its BOTTOM side → .plus
    //
    let dt_dx = dt / mesh.dx;
    let dt_dy = dt / mesh.dy;

    // Liang & Marche (2009) FLUX RESCALING for strictly mass-
    // conservative wetting/drying. Pass 1: for each cell, compute
    // the total OUTGOING mass over this timestep assuming the raw
    // face fluxes were applied. If the outflow would drain more
    // mass than the cell holds, compute a rescaling factor
    // `α ∈ [0, 1]` that caps the outflow at the available mass.
    // The same `α` is then applied (in pass 2) to ALL three
    // components of the face flux (mass, x-mom, y-mom) at each
    // outgoing face — momentum advects with the rescaled mass.
    //
    // For each face, the scaling factor used is that of the
    // UPSTREAM cell — the cell that is LOSING mass through this
    // face. This guarantees: (a) no cell drains below H_DRY → 0,
    // (b) the same scaled flux is seen by both cells sharing the
    // face → mass is conserved exactly across the face, (c) on
    // wet-wet flows where no cell hits the drain limit, α ≡ 1
    // and the scheme reduces to the unrescaled FV update.
    //
    // α is T-typed so its dependence on h (the available mass)
    // carries the derivative through the rescaling.
    zip_for_each!(Zip::indexed(&mut *alpha), |(i, j), a| {
        let fx_right = faces_x[(i, j + 1)].minus.mass;
        let fx_left = faces_x[(i, j)].plus.mass;
        let fy_bottom = faces_y[(i + 1, j)].minus.mass;
        let fy_top = faces_y[(i, j)].plus.mass;

        // Cell (i,j) outflow: positive mass flux on the right side
        // means mass flows L→R (cell loses); negative on the left
        // means mass flows R→L (cell loses); analogously in y.
        let mut out_mass = T::zero();
        if fx_right.value() > 0.0 {
            out_mass = out_mass + fx_right * dt_dx;
        }
        if fx_left.value() < 0.0 {
            out_mass = out_mass + (-fx_left) * dt_dx;
        }
        if fy_bottom.value() > 0.0 {
            out_mass = out_mass + fy_bottom * dt_dy;
        }
        if fy_top.value() < 0.0 {
            out_mass = out_mass + (-fy_top) * dt_dy;
        }

        let available = (states[(i, j)].h - H_DRY).max(T::zero());
        *a = if out_mass.value() > available.value() && out_mass.value() > 0.0 {
            available / out_mass
        } else {
            T::one()
        };
    });

    // Pass 2: for each face, determine the upstream cell (the side
    // losing mass) and scale by that side's α, IN PLACE and once per
    // face. Same α applies to .minus and .plus (they differ only in
    // Audusse's hydrostatic-pressure correction, which is also
    // rescaled). The former formulation evaluated this scaling lazily
    // per cell, re-scaling every interior face twice (once from each
    // adjacent cell). Interior dry-dry faces are identically zero and
    // α of a dry cell is one, so they are skipped.
    zip_for_each!(Zip::indexed(&mut *faces_x), |(i, j), face| 'cell: {
        if j > 0 && j < n_cols && was_dry[(i, j - 1)] && was_dry[(i, j)] {
            break 'cell;
        }
        let alpha_up = if j == 0 {
            // West boundary face. Outflow from cell (i, 0) when face.mass < 0
            // (mass moves R→L, i.e. cell 0 loses mass into the ghost).
            if face.minus.mass.value() < 0.0 {
                alpha[(i, 0)]
            } else {
                T::one()
            }
        } else if j == n_cols {
            // East boundary face. Outflow from cell (i, n_cols-1) when
            // face.mass > 0.
            if face.minus.mass.value() > 0.0 {
                alpha[(i, n_cols - 1)]
            } else {
                T::one()
            }
        } else if face.minus.mass.value() > 0.0 {
            // Interior, L→R: upstream is left cell (j-1).
            alpha[(i, j - 1)]
        } else {
            // Interior, R→L (or zero): upstream is right cell (j).
            alpha[(i, j)]
        };
        face.minus.mass = face.minus.mass * alpha_up;
        face.minus.x_momentum = face.minus.x_momentum * alpha_up;
        face.minus.y_momentum = face.minus.y_momentum * alpha_up;
        face.plus.mass = face.plus.mass * alpha_up;
        face.plus.x_momentum = face.plus.x_momentum * alpha_up;
        face.plus.y_momentum = face.plus.y_momentum * alpha_up;
    });
    zip_for_each!(Zip::indexed(&mut *faces_y), |(i, j), face| 'cell: {
        if i > 0 && i < n_rows && was_dry[(i - 1, j)] && was_dry[(i, j)] {
            break 'cell;
        }
        let alpha_up = if i == 0 {
            if face.minus.mass.value() < 0.0 {
                alpha[(0, j)]
            } else {
                T::one()
            }
        } else if i == n_rows {
            if face.minus.mass.value() > 0.0 {
                alpha[(n_rows - 1, j)]
            } else {
                T::one()
            }
        } else if face.minus.mass.value() > 0.0 {
            alpha[(i - 1, j)]
        } else {
            alpha[(i, j)]
        };
        face.minus.mass = face.minus.mass * alpha_up;
        face.minus.x_momentum = face.minus.x_momentum * alpha_up;
        face.minus.y_momentum = face.minus.y_momentum * alpha_up;
        face.plus.mass = face.plus.mass * alpha_up;
        face.plus.x_momentum = face.plus.x_momentum * alpha_up;
        face.plus.y_momentum = face.plus.y_momentum * alpha_up;
    });

    // FV update with scaled face fluxes + explicit bed-slope source.
    //
    // Isolated-dry-cell fast path: if the cell and all four immediate
    // neighbours were dry (`h ≤ H_DRY`) at the START of the step, the
    // four face fluxes around the cell are all zero by the dry-dry
    // short-circuit above. The standard path would compute
    // `dh = dhu = dhv = 0` and only the bed-slope source `s_hu, s_hv`
    // could be non-zero — but those would feed momentum into an
    // essentially dry cell, which the `new_h ≤ H_VEL` floor below
    // zeroes anyway. Keeping `h` and zeroing just the momentum is
    // mathematically identical to running the full body (dh = 0, so
    // the floor keeps the same mass) and lets us skip ~12 FLOPs per
    // dry interior cell. We restrict the shortcut to STRICT interior
    // cells (`1 ≤ i ≤ n_rows − 2`, `1 ≤ j ≤ n_cols − 2`) so that
    // boundary cells — which could receive injected mass from a
    // `Discharge` / `Depth` ghost — keep the full path.
    zip_for_each!(Zip::indexed(&mut *states), |(i, j), cell| 'cell: {
        if was_dry[(i, j)]
            && i > 0
            && i + 1 < n_rows
            && j > 0
            && j + 1 < n_cols
            && was_dry[(i - 1, j)]
            && was_dry[(i + 1, j)]
            && was_dry[(i, j - 1)]
            && was_dry[(i, j + 1)]
        {
            // Keep the moisture film (see the floor below); only
            // the momentum is reset, matching what the full path
            // does for h ≤ H_VEL.
            cell.hu = T::zero();
            cell.hv = T::zero();
            break 'cell;
        }
        let fx_right = faces_x[(i, j + 1)].minus;
        let fx_left = faces_x[(i, j)].plus;
        let fy_bottom = faces_y[(i + 1, j)].minus;
        let fy_top = faces_y[(i, j)].plus;

        let dh = (fx_right.mass - fx_left.mass) * dt_dx + (fy_bottom.mass - fy_top.mass) * dt_dy;
        let new_h = cell.h - dh;

        // After flux rescaling new_h ≥ H_DRY by construction
        // for cells that started wet. A residual floor catches
        // floating-point roundoff that could nudge a barely-
        // wet cell below the threshold.
        let dhu = (fx_right.x_momentum - fx_left.x_momentum) * dt_dx
            + (fy_bottom.x_momentum - fy_top.x_momentum) * dt_dy;
        let dhv = (fx_right.y_momentum - fx_left.y_momentum) * dt_dx
            + (fy_bottom.y_momentum - fy_top.y_momentum) * dt_dy;

        // Explicit bed-slope source — see module doc. Algebraic
        // form `S = (g/2) · (h_R_at_face² − h_L_at_face²)/Δx`
        // that cancels the pressure-flux divergence exactly for
        // lake-at-rest on any bed. `h_face` from cell's own `η`.
        let h_old = cell.h;
        let eta_cell = h_old + mesh.bed[(i, j)];

        let h_face_left = (eta_cell - z_face_x[(i, j)]).max(T::zero());
        let h_face_right = (eta_cell - z_face_x[(i, j + 1)]).max(T::zero());

        // Bounded-slope safety valve — see
        // docs/bug-report-2026-07-boundary-slope-instability.md §7.
        // An EARLIER version of this guard fired whenever `h_face` on
        // either side exceeded `h_old`, on the (wrong) assumption that
        // this could only happen at a bed jump much larger than the
        // depth. That is false in general: on an ORDINARY (non-
        // shoreline) sloped face, `h_face_downhill = h_old + Δz/2` for
        // ANY `Δz > 0` — lake-at-rest on a sloped bed legitimately
        // exceeds `h_old` on the downhill side by design (that
        // asymmetry is exactly what the well-balanced cancellation
        // needs). A SECOND earlier version required BOTH faces of the
        // cell to be interior/ordinary before applying the cap, which
        // excluded the very cells the bug report's reproducer blows up
        // first (column adjacent to an open boundary, where one face
        // is the boundary itself but the OTHER face is an ordinary
        // interior face with the same pathology).
        //
        // The corrected guard caps each SIDE of the quadratic term
        // independently: `h_face_side` is soft-clamped to at most
        // `STEEP_SOURCE_RATIO · h_old` whenever that specific side (a)
        // is an ordinary face — the neighbour across it shares this
        // cell's pre-step wet/dry state, i.e. `z_face` came from the
        // midpoint rule, never the shoreline `max(z_L, z_R)` rule —
        // and (b) actually exceeds the cap. `max(z_L, z_R) ≥` this
        // cell's own bed by construction, so a genuine shoreline
        // (buildings, coastlines) always has `h_face ≤ h_old` on the
        // wet side and can never trip the cap; boundary faces
        // (Transmissive/Wall) carry zero jump by construction and
        // never trip it either — those C-property tests, and both
        // domain edges, are untouched by construction, not by tuning
        // the ratio. The quadratic FORM is preserved unchanged (still
        // a difference of squares, still exact in the smooth-bed
        // limit); only a face whose implied depth has run away to many
        // times the cell's own depth — the signature of a bed jump
        // dominating a thin film, not a resolved free surface — gets
        // pulled back before squaring, so its contribution collapses
        // toward zero (∝ h_old²) instead of growing with the bed jump
        // squared, independent of `h`.
        let cap_x = h_old * STEEP_SOURCE_RATIO;
        let ordinary_left_x = j > 0 && was_dry[(i, j - 1)] == was_dry[(i, j)];
        let ordinary_right_x = j + 1 < n_cols && was_dry[(i, j)] == was_dry[(i, j + 1)];
        let h_face_left_capped = if ordinary_left_x && h_face_left.value() > cap_x.value() {
            cap_x
        } else {
            h_face_left
        };
        let h_face_right_capped = if ordinary_right_x && h_face_right.value() > cap_x.value() {
            cap_x
        } else {
            h_face_right
        };
        let s_hu =
            (h_face_right_capped.powi(2) - h_face_left_capped.powi(2)) * (0.5 * GRAVITY) / mesh.dx;

        let h_face_top = (eta_cell - z_face_y[(i, j)]).max(T::zero());
        let h_face_bottom = (eta_cell - z_face_y[(i + 1, j)]).max(T::zero());
        let cap_y = h_old * STEEP_SOURCE_RATIO;
        let ordinary_top_y = i > 0 && was_dry[(i - 1, j)] == was_dry[(i, j)];
        let ordinary_bottom_y = i + 1 < n_rows && was_dry[(i, j)] == was_dry[(i + 1, j)];
        let h_face_top_capped = if ordinary_top_y && h_face_top.value() > cap_y.value() {
            cap_y
        } else {
            h_face_top
        };
        let h_face_bottom_capped = if ordinary_bottom_y && h_face_bottom.value() > cap_y.value() {
            cap_y
        } else {
            h_face_bottom
        };
        let s_hv = (h_face_bottom_capped.powi(2) - h_face_top_capped.powi(2)) * (0.5 * GRAVITY)
            / mesh.dy;

        if new_h.value() <= H_VEL {
            // Moisture floor. The mass is KEPT: a wetting-front
            // cell that received δ ≤ H_DRY of inflow had that mass
            // already leave its neighbour through the shared face,
            // and zeroing it here destroyed volume all along the
            // front perimeter at every step. The film stays inert
            // on the bed (its faces are dry-dry for h ≤ H_DRY, and
            // the α available-mass term is zero) until inflow
            // accumulates past the threshold. Momentum is dropped
            // for anything up to H_VEL: velocity on such a film is
            // meaningless, and residual hu with h barely above
            // H_DRY collapses dt through hu/h (see H_VEL doc).
            *cell = Conserved2DG::new_generic(new_h.max(T::zero()), T::zero(), T::zero());
        } else {
            // A cell that was a film inherits no momentum: the
            // floor invariant keeps evolved films at hu = hv = 0,
            // and a caller-supplied initial state violating it
            // must not leak its stale momentum into the wetted
            // cell (hu/h on the old film depth is unphysical).
            let (hu_old, hv_old) = if h_old.value() > H_VEL {
                (cell.hu, cell.hv)
            } else {
                (T::zero(), T::zero())
            };
            cell.h = new_h;
            cell.hu = hu_old - dhu + s_hu * dt;
            cell.hv = hv_old - dhv + s_hv * dt;
        }
    });
}

/// Strong-stability-preserving Runge-Kutta second-order step
/// (Shu & Osher 1988):
///
/// ```text
///   U^(1)    = U^n + dt · L(U^n)              (forward-Euler predictor)
///   U^(2)    = U^(1) + dt · L(U^(1))          (forward-Euler corrector)
///   U^{n+1}  = ½ · U^n + ½ · U^(2)            (convex combination)
/// ```
///
/// Equivalent to the Heun method, but written so it is manifestly a
/// convex combination of two forward-Euler updates. This is what makes
/// it **strong-stability-preserving**: every property that the
/// forward-Euler step preserves under a CFL bound — non-negative
/// depth, finite momentum, mass conservation under wall BCs — is
/// inherited by the SSP-RK2 update without any extra CFL constraint
/// (`cfl_SSPRK2 ≤ cfl_FE`). The positivity clamp inside
/// [`forward_euler_step`] composes cleanly: the predictor produces a
/// valid wet/dry state, the corrector produces another, the average
/// is the convex combination of two non-negative states (thus
/// non-negative). See Gottlieb, Ketcheson & Shu (2009) for the
/// general theory of SSP methods.
///
/// In combination with MUSCL slope-limited spatial reconstruction
/// (Audusse + minmod), the scheme is second-order in space AND time
/// — the natural pairing that gives the full benefit of MUSCL.
/// Without SSP-RK2 the time-error of forward Euler dominates and the
/// asymptotic order drops to ~1.5 in practice.
///
/// `dt` is bounded by the same `cfl_time_step` as forward Euler.
/// Panics on shape mismatch (same as [`forward_euler_step`]).
pub fn ssprk2_step<T: Real + MaybeSendSync>(
    states: &mut Array2<Conserved2DG<T>>,
    mesh: &Mesh2DG<T>,
    bcs: Boundaries2D,
    dt: f64,
) {
    let mut ws = StepWorkspace2D::for_mesh(mesh);
    ssprk2_step_with(states, mesh, bcs, dt, &mut ws);
}

/// [`ssprk2_step`] over caller-owned scratch buffers — the
/// allocation-free hot path (the U^n snapshot lives in the workspace).
pub fn ssprk2_step_with<T: Real + MaybeSendSync>(
    states: &mut Array2<Conserved2DG<T>>,
    mesh: &Mesh2DG<T>,
    bcs: Boundaries2D,
    dt: f64,
    ws: &mut StepWorkspace2D<T>,
) {
    // Snapshot U^n into the workspace buffer (copy, no allocation).
    ws.u_n.assign(states);

    // Predictor: states becomes U^(1) = U^n + dt L(U^n).
    forward_euler_step_with(states, mesh, bcs, dt, ws);

    // Corrector: states becomes U^(2) = U^(1) + dt L(U^(1)).
    forward_euler_step_with(states, mesh, bcs, dt, ws);

    // Convex combination: U^{n+1} = ½(U^n + U^(2)).
    let u_n = &ws.u_n;
    zip_for_each!(Zip::indexed(&mut *states), |(i, j), s| {
        let prev = u_n[(i, j)];
        s.h = (prev.h + s.h) * 0.5;
        s.hu = (prev.hu + s.hu) * 0.5;
        s.hv = (prev.hv + s.hv) * 0.5;
        // Re-apply the moisture floor to the averaged state (one of
        // the two halves may have been wet but very shallow): keep the
        // mass, drop the momentum — same rule as the Euler floor.
        if s.h.value() <= H_VEL {
            s.h = s.h.max(T::zero());
            s.hu = T::zero();
            s.hv = T::zero();
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Mesh2D;
    use crate::state::Conserved2D;

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
    fn lake_at_rest_with_emerged_island_is_preserved() {
        // C-property at shorelines: a lake at rest around an island
        // whose crest stands ABOVE the free surface must stay exactly
        // at rest. With the interpolated face bed ½(z_L + z_R) at
        // wet/dry faces, the face bed can sit below the wet-side
        // surface while the dry cell stands above it: the α-rescaling
        // zeroes the (spurious) face flux but the cell-centered source
        // keeps a g·h_face²/(2Δ) pressure term, accelerating the
        // shoreline cells toward the island every step. The wet/dry
        // face rule z_face = max(z_L, z_R) closes both sides to zero
        // depth and restores the exact balance.
        let n = 21;
        let eta = 1.0;
        let mesh = gaussian_bed(n, n, 1.0, 1.0, 2.5); // crest 2.5 m > η
        let initial = Array2::from_shape_fn((n, n), |(i, j)| {
            let h = (eta - mesh.bed[(i, j)]).max(0.0);
            Conserved2D::new(h, 0.0, 0.0)
        });
        let mut states = initial.clone();
        let m0 = total_mass(&states, 1.0, 1.0);

        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        for ((i, j), s) in states.indexed_iter() {
            assert!(
                (s.h - initial[(i, j)].h).abs() < 1e-10,
                "h drifted at ({i},{j}): {} vs {}",
                s.h,
                initial[(i, j)].h
            );
            assert!(
                s.hu.abs() < 1e-10 && s.hv.abs() < 1e-10,
                "spurious shoreline momentum at ({i},{j}): hu = {:e}, hv = {:e}",
                s.hu,
                s.hv
            );
        }
        let m1 = total_mass(&states, 1.0, 1.0);
        assert_relative_eq!(m0, m1, epsilon = 1e-12);
    }

    #[test]
    fn lake_at_rest_against_emerged_bank_is_preserved() {
        // Same C-property, 1D-like geometry: a linear bank rising
        // through the free surface (the "thin layer against a slope"
        // configuration reviewers ask for). Exercises the x-axis
        // shoreline correction in isolation.
        let n_rows = 4;
        let n_cols = 30;
        let eta = 0.5;
        // Bed rises from -0.5 to +1.5 across the columns; the
        // shoreline sits mid-domain.
        let bed = Array2::from_shape_fn((n_rows, n_cols), |(_i, j)| {
            -0.5 + 2.0 * j as f64 / (n_cols - 1) as f64
        });
        let mesh = Mesh2D::new(bed, 1.0, 1.0, 0.0);
        let initial = Array2::from_shape_fn((n_rows, n_cols), |(i, j)| {
            let h = (eta - mesh.bed[(i, j)]).max(0.0);
            Conserved2D::new(h, 0.0, 0.0)
        });
        let mut states = initial.clone();

        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        for ((i, j), s) in states.indexed_iter() {
            assert!(
                (s.h - initial[(i, j)].h).abs() < 1e-10,
                "h drifted at ({i},{j}): {} vs {}",
                s.h,
                initial[(i, j)].h
            );
            assert!(
                s.hu.abs() < 1e-10 && s.hv.abs() < 1e-10,
                "spurious shoreline momentum at ({i},{j}): hu = {:e}, hv = {:e}",
                s.hu,
                s.hv
            );
        }
    }

    #[test]
    fn draining_cell_inflow_to_dry_neighbour_is_not_discarded() {
        // Regression test for the isolated-dry-cell fast path: it must
        // consult the PRE-step wet/dry pattern, not the in-place
        // mutated `states`. Cell (0, 1) starts wet with strong +y
        // momentum; the α-rescaling caps its outflow at `h − H_DRY`,
        // so it drains to the floor and is reset to dry IN THE SAME
        // STEP in which it exports its mass into the dry interior
        // cell (1, 1). A fast path reading post-update states sees
        // (0, 1) already dry, classifies (1, 1) as isolated-dry, and
        // discards the inflow — destroying most of the domain's mass.
        let mesh = flat_mesh(3, 3, 1.0, 1.0);
        let mut states = Array2::from_elem((3, 3), Conserved2D::DRY);
        let h0 = 0.05;
        states[(0, 1)] = Conserved2D::new(h0, 0.0, h0 * 5.0);
        let m0 = total_mass(&states, 1.0, 1.0);

        // dt chosen so the raw outflow through the (0,1)/(1,1) face
        // exceeds the available mass and the rescaling engages. CFL
        // is irrelevant here: a single step, and the α cap bounds the
        // update regardless of dt.
        forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, 0.25);

        assert!(
            states[(0, 1)].h <= H_DRY,
            "source cell should have drained, h = {}",
            states[(0, 1)].h
        );
        assert!(
            states[(1, 1)].h > 0.5 * h0,
            "inflow into the dry interior cell was discarded: h(1,1) = {}",
            states[(1, 1)].h
        );
        // With the moisture floor the drained source cell keeps its
        // residual film, so the step conserves mass to roundoff.
        let m1 = total_mass(&states, 1.0, 1.0);
        assert_relative_eq!(m0, m1, epsilon = 1e-12);
    }

    #[test]
    fn reused_workspace_is_bit_identical_to_allocating_step() {
        // The workspace buffers are overwritten in full at every step;
        // nothing may leak from one step into the next. Run a
        // wetting/drying-heavy scenario (dam break over an emerged
        // bump, walls) twice — once through the allocating wrapper,
        // once through a single reused workspace — and require
        // bit-identical states after every step, for both integrators.
        let n = 20;
        let mesh = gaussian_bed(n, n, 1.0, 1.0, 1.5);
        let init = Array2::from_shape_fn((n, n), |(_i, j)| {
            if j < n / 3 {
                Conserved2D::new(2.0, 0.0, 0.0)
            } else {
                Conserved2D::DRY
            }
        });

        for use_ssprk2 in [false, true] {
            let mut a = init.clone();
            let mut b = init.clone();
            let mut ws = StepWorkspace2D::for_mesh(&mesh);
            for step in 0..50 {
                let dt = cfl_time_step(&a, &mesh, 0.4);
                if use_ssprk2 {
                    ssprk2_step(&mut a, &mesh, Boundaries2D::WALLS, dt);
                    ssprk2_step_with(&mut b, &mesh, Boundaries2D::WALLS, dt, &mut ws);
                } else {
                    forward_euler_step(&mut a, &mesh, Boundaries2D::WALLS, dt);
                    forward_euler_step_with(&mut b, &mesh, Boundaries2D::WALLS, dt, &mut ws);
                }
                for ((i, j), sa) in a.indexed_iter() {
                    let sb = b[(i, j)];
                    assert!(
                        sa.h == sb.h && sa.hu == sb.hu && sa.hv == sb.hv,
                        "step {step} (ssprk2 = {use_ssprk2}) diverged at ({i},{j}): \
                         {sa:?} vs {sb:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn wetting_front_conserves_mass_in_closed_box() {
        // Dam break onto a dry half of a closed box. Every cell the
        // front touches first receives δ ≤ H_DRY of mass; the old
        // floor reset those cells to exact dry, destroying δ along
        // the whole front perimeter at every step — a systematic
        // volume leak this test bounds at roundoff now that the
        // moisture floor keeps sub-threshold films.
        let n = 30;
        let mesh = flat_mesh(n, n, 1.0, 1.0);
        let mut states = Array2::from_shape_fn((n, n), |(_i, j)| {
            if j < n / 2 {
                Conserved2D::new(1.0, 0.0, 0.0)
            } else {
                Conserved2D::DRY
            }
        });
        let m0 = total_mass(&states, 1.0, 1.0);

        for _ in 0..300 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        for s in &states {
            assert!(s.h.is_finite() && s.h >= 0.0, "bad depth: {}", s.h);
        }
        let m1 = total_mass(&states, 1.0, 1.0);
        assert_relative_eq!(m0, m1, epsilon = 1e-12);
    }

    #[test]
    fn residual_momentum_on_thin_film_does_not_collapse_dt() {
        // A film barely above H_DRY carrying stale momentum (as a
        // user-supplied IC can produce) must not drive u = hu/h into
        // the CFL bound: with the H_VEL cutoff the film contributes
        // only its celerity, and the first Euler step zeroes its
        // momentum.
        let mesh = flat_mesh(4, 4, 1.0, 1.0);
        let mut states = Array2::from_elem((4, 4), Conserved2D::new(1.0, 0.0, 0.0));
        states[(1, 1)] = Conserved2D::new(1.1 * H_DRY, 1.0, 1.0); // u ~ 1e6 m/s if honoured
        let dt = cfl_time_step(&states, &mesh, 0.4);
        // Wet-cell celerity dominates: dt is the ordinary bound, not
        // a hu/h-collapsed one.
        let c = (GRAVITY * 1.0_f64).sqrt();
        let dt_expected = 0.4 / (2.0 * c); // dx = dy = 1
        assert!(
            dt > 0.9 * dt_expected,
            "dt collapsed by thin-film residual momentum: dt = {dt:e}, expected ~{dt_expected:e}"
        );
        forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        // The film may legitimately wet up past H_VEL in this step
        // (deep neighbours dump mass into it); what must NOT survive
        // is the stale 1e6 m/s velocity. Bound |u| by the physical
        // dam-break scale (2c ≈ 6.3 m/s for h = 1 m).
        let s = states[(1, 1)];
        if s.h > H_VEL {
            let u = (s.hu / s.h).abs().max((s.hv / s.h).abs());
            assert!(
                u < 10.0,
                "stale film momentum leaked into the wetted cell: |u| = {u:e}"
            );
        } else {
            assert_eq!(s.hu, 0.0);
            assert_eq!(s.hv, 0.0);
        }
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

    // -----------------------------------------------------------------
    // SSP-RK2 tests (the predictor-corrector + convex-combination time
    // integrator). The forward-Euler tests above cover the FV operator
    // L(U); these focus on the time-integration shell.
    // -----------------------------------------------------------------

    #[test]
    fn ssprk2_lake_at_rest_on_flat_bed_is_preserved_exactly() {
        // The convex combination of two forward-Euler updates of an
        // exactly-preserved state must itself preserve the state.
        let n_rows = 10;
        let n_cols = 12;
        let mesh = flat_mesh(n_rows, n_cols, 1.0, 1.0);
        let mut states = Array2::from_elem((n_rows, n_cols), Conserved2D::new(2.0, 0.0, 0.0));
        let dt = 0.01;
        for _ in 0..100 {
            ssprk2_step(&mut states, &mesh, Boundaries2D::TRANSMISSIVE, dt);
        }
        for s in &states {
            assert_relative_eq!(s.h, 2.0, epsilon = 1e-12);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-12);
            assert_relative_eq!(s.hv, 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn ssprk2_lake_at_rest_on_diagonal_bed_is_preserved() {
        // The well-balanced 2D test under SSP-RK2 — the C-property
        // must hold under the full predictor-corrector scheme, not
        // just the single Euler step.
        let n_rows = 15;
        let n_cols = 20;
        let mesh = diagonal_sloped_mesh(n_rows, n_cols, 1.0, 1.0, 0.03, 0.04);
        let eta = 5.0;
        let initial = lake_at_rest_on(&mesh, eta);
        let mut states = initial.clone();
        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            ssprk2_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        for ((i, j), s) in states.indexed_iter() {
            assert_relative_eq!(s.h, initial[(i, j)].h, epsilon = 1e-10);
            assert_relative_eq!(s.hu, 0.0, epsilon = 1e-10);
            assert_relative_eq!(s.hv, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn ssprk2_mass_is_conserved_with_walls() {
        // SSP-RK2 inherits the telescoping property of forward Euler:
        // each sub-step is conservative under walls, the convex
        // combination of conservative states is conservative.
        let n_rows = 25;
        let n_cols = 25;
        let dx = 1.0;
        let dy = 1.0;
        let mesh = flat_mesh(n_rows, n_cols, dx, dy);
        let mut states = gaussian_bump_2d(n_rows, n_cols, 1.0, 0.5);
        let m0 = total_mass(&states, dx, dy);
        for _ in 0..200 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            ssprk2_step(&mut states, &mesh, Boundaries2D::WALLS, dt);
        }
        let m1 = total_mass(&states, dx, dy);
        assert_relative_eq!(m0, m1, epsilon = 1e-10);
    }

    #[test]
    fn ssprk2_remains_isotropic_on_uniform_mesh() {
        // The averaging step is variable-by-variable and direction-
        // agnostic; isotropy must survive.
        let n = 21;
        let mesh = flat_mesh(n, n, 1.0, 1.0);
        let mut states = gaussian_bump_2d(n, n, 1.0, 0.5);
        for _ in 0..50 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            ssprk2_step(&mut states, &mesh, Boundaries2D::TRANSMISSIVE, dt);
        }
        for i in 0..n {
            for j in 0..n {
                assert_relative_eq!(states[(i, j)].h, states[(n - 1 - i, j)].h, epsilon = 1e-10);
                assert_relative_eq!(states[(i, j)].h, states[(j, i)].h, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn ssprk2_preserves_positivity_at_dry_bed_interface() {
        // A wet-on-left, dry-on-right configuration must not produce
        // negative depths under SSP-RK2 (each Euler sub-step clamps
        // to DRY when needed; the convex combination of two clamped
        // states is also non-negative).
        let n_rows = 3;
        let n_cols = 100;
        let mesh = flat_mesh(n_rows, n_cols, 1.0, 1.0);
        let mut states = Array2::from_elem((n_rows, n_cols), Conserved2D::DRY);
        for i in 0..n_rows {
            for j in 0..n_cols / 2 {
                states[(i, j)] = Conserved2D::new(1.0, 0.0, 0.0);
            }
        }
        for _ in 0..300 {
            let dt = cfl_time_step(&states, &mesh, 0.4);
            if !dt.is_finite() {
                break;
            }
            ssprk2_step(&mut states, &mesh, Boundaries2D::TRANSMISSIVE, dt);
        }
        for s in &states {
            assert!(s.h.is_finite() && s.h >= 0.0, "h ill-formed: {}", s.h);
            assert!(s.hu.is_finite(), "hu non-finite: {}", s.hu);
            assert!(s.hv.is_finite(), "hv non-finite: {}", s.hv);
        }
    }

    // ----- Generic-over-Real: M1 evidence on the full FV update. -----

    #[test]
    fn forward_euler_ad_seed_on_uniform_depth_matches_finite_difference() {
        // Seed the initial uniform lake depth h_init as the variable.
        // After one forward-Euler step on a flat-bed walled domain,
        // the total mass remains h_init · (n_rows · n_cols · dx · dy).
        // d(total_mass)/d(h_init) = n_rows · n_cols · dx · dy.
        // AD must recover this analytic value; FD provides the
        // independent cross-check.
        use hydroflux_autograd::{Dual, Real as _};

        fn one_step_mass<T: hydroflux_autograd::Real + MaybeSendSync>(h_init: T) -> T {
            let n_rows = 4usize;
            let n_cols = 5usize;
            let dx = 1.0;
            let dy = 1.0;
            let bed = Array2::<T>::from_elem((n_rows, n_cols), T::zero());
            let mesh = Mesh2DG::<T>::new(bed, dx, dy, T::zero());
            let mut states = Array2::from_elem(
                (n_rows, n_cols),
                Conserved2DG::<T>::new_generic(h_init, T::zero(), T::zero()),
            );
            forward_euler_step(&mut states, &mesh, Boundaries2D::WALLS, 0.01);
            let mut m = T::zero();
            for s in &states {
                m = m + s.h * (dx * dy);
            }
            m
        }
        let h_val = 2.0_f64;
        let eps = 1e-6_f64;
        let fd =
            (one_step_mass::<f64>(h_val + eps) - one_step_mass::<f64>(h_val - eps)) / (2.0 * eps);
        let ad = one_step_mass::<Dual>(Dual::variable(h_val));
        // FD against the analytical value (the cell count times area).
        // The AD result must match the FD result to ~1e-7 relative.
        assert_relative_eq!(ad.val, 4.0 * 5.0 * 1.0 * 1.0 * 2.0, epsilon = 1e-10);
        assert_relative_eq!(ad.dval, fd, epsilon = 1e-6);
        // And independently to the analytic value.
        assert_relative_eq!(ad.dval, 20.0, epsilon = 1e-10);
    }

    #[test]
    fn ssprk2_ad_seed_on_bed_elevation_matches_finite_difference() {
        // Seed a SINGLE bed-elevation value through the whole bed
        // (uniform bed shift z₀) on a wet domain at rest. After one
        // SSP-RK2 step, η = h + z₀ stays uniform so the lake-at-rest
        // property holds; mass is preserved. The interesting test:
        // d(total_mass) / d(z₀) = 0 because mass is conserved
        // regardless of the bed reference level. AD must return
        // approximately zero (modulo floating-point noise from the
        // many arithmetic ops in the SSP-RK2 update). This is the
        // strongest kind of M1 evidence: a known *invariance* that
        // the gradient must respect.
        use hydroflux_autograd::{Dual, Real as _};

        fn ssprk_mass<T: hydroflux_autograd::Real + MaybeSendSync>(z_seed: T) -> T {
            let n_rows = 6usize;
            let n_cols = 6usize;
            let dx = 1.0;
            let dy = 1.0;
            let bed = Array2::<T>::from_elem((n_rows, n_cols), z_seed);
            let mesh = Mesh2DG::<T>::new(bed, dx, dy, T::zero());
            // η = 3 constant ⇒ h = 3 − z_seed
            let h = T::from_f64(3.0) - z_seed;
            let mut states = Array2::from_elem(
                (n_rows, n_cols),
                Conserved2DG::<T>::new_generic(h, T::zero(), T::zero()),
            );
            ssprk2_step(&mut states, &mesh, Boundaries2D::WALLS, 0.005);
            let mut m = T::zero();
            for s in &states {
                m = m + s.h * (dx * dy);
            }
            m
        }
        let z_val = 0.5_f64;
        let eps = 1e-6_f64;
        let fd = (ssprk_mass::<f64>(z_val + eps) - ssprk_mass::<f64>(z_val - eps)) / (2.0 * eps);
        let ad = ssprk_mass::<Dual>(Dual::variable(z_val));
        // Analytic: d(total mass)/dz = −(n_rows · n_cols · dx · dy)
        // because h = 3 − z applied uniformly. Mass is preserved by
        // the FV step, so it depends linearly on the initial uniform h.
        let analytic = -(6.0 * 6.0 * 1.0 * 1.0);
        assert_relative_eq!(ad.dval, fd, epsilon = 1e-6);
        assert_relative_eq!(ad.dval, analytic, epsilon = 1e-8);
    }
}
