//! Criterion regression benches for the 2D time-step hot loop.
//!
//! Two régimes, matching how the solver is actually used:
//!
//! - `all_wet_256`: fully wet 256×256 flat-bed domain — the m2-paper
//!   configuration. Measures the raw face-flux + rescaling + update
//!   cost with no dry-cell short-circuits firing.
//! - `mostly_dry_256`: a 256×256 valley with a wet channel strip
//!   (~6 % wet), the Huasco-like régime where the dry-cell skips
//!   dominate. Guards the skip machinery against regressions.
//!
//! Baseline at the time of writing (2026-07-02, serial, quiet
//! machine): `euler_all_wet` ≈ 16 ms/step (≈ 246 ns/cell),
//! `ssprk2_all_wet` ≈ 31 ms/step, `euler_mostly_dry` ≈ 6.3 ms/step
//! (the dry-cell short-circuits buy ~2.6× on the ~94 %-dry valley).
//! Run with:
//!
//! ```text
//! cargo bench -p hydroflux-solver-2d
//! ```

use criterion::{Criterion, criterion_group, criterion_main};
use hydroflux_solver_2d::{
    Boundaries2D, Conserved2D, Mesh2D, cfl_time_step, forward_euler_step, ssprk2_step,
};
use ndarray::Array2;
use std::hint::black_box;

const N: usize = 256;

fn all_wet() -> (Mesh2D, Array2<Conserved2D>) {
    let mesh = Mesh2D::new(Array2::<f64>::zeros((N, N)), 1.0, 1.0, 0.03);
    // Gentle depth gradient so the fluxes are non-trivial (a perfectly
    // uniform lake would exercise only the well-balanced fast paths).
    let states = Array2::from_shape_fn((N, N), |(i, j)| {
        let h = 1.0 + 0.2 * ((i as f64) / N as f64) + 0.1 * ((j as f64) / N as f64);
        Conserved2D::new(h, 0.1 * h, 0.05 * h)
    });
    (mesh, states)
}

fn mostly_dry() -> (Mesh2D, Array2<Conserved2D>) {
    // V-shaped valley in x with a wet channel strip in the middle
    // ~6 % of columns; the rest of the domain is dry hillslope.
    let bed = Array2::from_shape_fn((N, N), |(_i, j)| {
        let x = (j as f64 - N as f64 / 2.0).abs();
        0.05 * x
    });
    let mesh = Mesh2D::new(bed, 1.0, 1.0, 0.03);
    let states = Array2::from_shape_fn((N, N), |(i, j)| {
        let eta = 0.4;
        let h = (eta - mesh.bed[(i, j)]).max(0.0);
        if h > 0.0 {
            let _ = i;
            Conserved2D::new(h, 0.0, 0.1 * h)
        } else {
            Conserved2D::DRY
        }
    });
    (mesh, states)
}

fn bench_steps(c: &mut Criterion) {
    let mut group = c.benchmark_group("step_256");
    group.sample_size(20);

    let (mesh, states0) = all_wet();
    let dt = cfl_time_step(&states0, &mesh, 0.4);
    group.bench_function("euler_all_wet", |b| {
        b.iter_batched(
            || states0.clone(),
            |mut s| {
                forward_euler_step(&mut s, &mesh, Boundaries2D::WALLS, black_box(dt));
                s
            },
            criterion::BatchSize::LargeInput,
        )
    });
    group.bench_function("ssprk2_all_wet", |b| {
        b.iter_batched(
            || states0.clone(),
            |mut s| {
                ssprk2_step(&mut s, &mesh, Boundaries2D::WALLS, black_box(dt));
                s
            },
            criterion::BatchSize::LargeInput,
        )
    });

    let (mesh_d, states_d) = mostly_dry();
    let dt_d = cfl_time_step(&states_d, &mesh_d, 0.4);
    group.bench_function("euler_mostly_dry", |b| {
        b.iter_batched(
            || states_d.clone(),
            |mut s| {
                forward_euler_step(&mut s, &mesh_d, Boundaries2D::WALLS, black_box(dt_d));
                s
            },
            criterion::BatchSize::LargeInput,
        )
    });

    group.finish();
}

criterion_group!(benches, bench_steps);
criterion_main!(benches);
