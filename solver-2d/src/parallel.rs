//! Feature-gated CPU parallelism for the 2D step's per-cell passes
//! (WP4, 2026-07 pre-submission revision).
//!
//! The step's hot-loop passes (primitives, MUSCL slopes, face beds,
//! HLLC face fluxes, the Liang–Marche α rescaling, and the final FV
//! update) are pure per-cell maps over the reusable
//! [`crate::update::StepWorkspace2D`] buffers — see
//! `docs/auditoria-motor-2026-07.md` §3.2-3.3: after the dry-mask
//! snapshot fix, every pass reads only pre-step state, so cells can be
//! visited in any order. This module supplies the two pieces every
//! such pass needs, gated behind the `parallel` Cargo feature so the
//! default (serial) build carries no rayon dependency and no
//! `Send + Sync` requirement on [`crate::Real`] itself:
//!
//! - [`MaybeSendSync`]: a conditional bound — `Send + Sync` when
//!   `parallel` is on, a no-op otherwise. Applied to the specific
//!   functions that need it (`T: Real + MaybeSendSync`), never to the
//!   `Real` trait definition. A future reverse-mode tape type is not
//!   expected to be `Sync` (tapes are typically `Rc<RefCell<..>>`-
//!   backed); it should still compile and run on the serial path.
//! - [`zip_for_each!`]: dispatches `Zip::par_for_each` (rayon) or
//!   `Zip::for_each` (serial) for the exact same closure. ndarray's
//!   rayon `Zip` splits by recursively bisecting the outer (row) axis
//!   until a chunk is small enough for the work-stealing scheduler —
//!   i.e. row-chunked in effect, not one rayon task per cell/face.
//!   That distinction is the point: a per-face-granularity rayon
//!   experiment (pre-workspace, 2026-05) failed to scale because task
//!   overhead dominated ~50-100 FLOP cell bodies; row-sized chunks
//!   amortise that overhead over hundreds of cells.

/// `Send + Sync` when built with `feature = "parallel"`, a no-op bound
/// otherwise. See the module docs for why this lives on individual
/// functions rather than on `hydroflux_autograd::Real` itself.
#[cfg(feature = "parallel")]
pub trait MaybeSendSync: Send + Sync {}
#[cfg(feature = "parallel")]
impl<T: Send + Sync> MaybeSendSync for T {}

/// `Send + Sync` when built with `feature = "parallel"`, a no-op bound
/// otherwise. See the module docs for why this lives on individual
/// functions rather than on `hydroflux_autograd::Real` itself.
#[cfg(not(feature = "parallel"))]
pub trait MaybeSendSync {}
#[cfg(not(feature = "parallel"))]
impl<T> MaybeSendSync for T {}

/// `$zip.par_for_each($f)` under `feature = "parallel"`,
/// `$zip.for_each($f)` otherwise. See the module docs for why this
/// is row-chunked in effect despite calling a per-cell closure.
macro_rules! zip_for_each {
    ($zip:expr, $f:expr) => {{
        #[cfg(feature = "parallel")]
        {
            $zip.par_for_each($f);
        }
        #[cfg(not(feature = "parallel"))]
        {
            $zip.for_each($f);
        }
    }};
}
pub(crate) use zip_for_each;
