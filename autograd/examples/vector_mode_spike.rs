//! Spike: what would a vector-mode dual number actually buy?
//!
//! The scalar `Dual` carries one derivative, so a gradient over `P`
//! parameters costs `P` forward passes at 2.0x the primal each, which
//! puts the break-even against reverse mode at about two parameters
//! (`m1_forward_scaling`). A vector-mode dual carrying `N` derivatives
//! in one pass would amortise everything the `P` passes currently
//! repeat: control flow, mesh traversal, memory traffic, and the value
//! arithmetic itself.
//!
//! How much that is worth depends on a quantity we have not measured —
//! the cost of one additional derivative component once the shared work
//! is paid. Assuming it equals the scalar figure of 1.04x the primal
//! gives a break-even near four. Assuming perfect SIMD-like sharing
//! gives something much larger. The honest range is wide enough that
//! implementing the real thing on a guess would be careless.
//!
//! So this measures it on a kernel with the solver's arithmetic mix
//! (sqrt, division, powf, fused multiply-add over a state array) before
//! committing to the full implementation in `dual.rs`.
//!
//! Run:
//!   cargo run --release -p hydroflux-autograd --example vector_mode_spike

use std::time::Instant;

const CELLS: usize = 4096;
const STEPS: usize = 400;
const G: f64 = 9.81;

/// Vector-mode dual: one value, `N` derivative components.
#[derive(Clone, Copy)]
struct DualN<const N: usize> {
    val: f64,
    dval: [f64; N],
}

impl<const N: usize> DualN<N> {
    fn constant(v: f64) -> Self {
        Self { val: v, dval: [0.0; N] }
    }
    fn seeded(v: f64, k: usize) -> Self {
        let mut d = [0.0; N];
        d[k] = 1.0;
        Self { val: v, dval: d }
    }
    fn mul(self, o: Self) -> Self {
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = self.dval[i] * o.val + self.val * o.dval[i];
        }
        Self { val: self.val * o.val, dval: d }
    }
    fn add(self, o: Self) -> Self {
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = self.dval[i] + o.dval[i];
        }
        Self { val: self.val + o.val, dval: d }
    }
    fn div(self, o: Self) -> Self {
        let inv = 1.0 / o.val;
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = (self.dval[i] * o.val - self.val * o.dval[i]) * inv * inv;
        }
        Self { val: self.val * inv, dval: d }
    }
    fn sqrt(self) -> Self {
        let s = self.val.sqrt();
        let f = if s == 0.0 { 0.0 } else { 0.5 / s };
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = self.dval[i] * f;
        }
        Self { val: s, dval: d }
    }
    fn powf(self, p: f64) -> Self {
        let v = self.val.powf(p);
        let f = p * self.val.powf(p - 1.0);
        let mut d = [0.0; N];
        for i in 0..N {
            d[i] = self.dval[i] * f;
        }
        Self { val: v, dval: d }
    }
}

/// One step of arithmetic shaped like the solver's inner loop: a wave
/// speed (sqrt), a Manning-like friction factor (powf and division),
/// and the update itself. The control flow and the array traversal are
/// the part a vector mode gets to share.
fn kernel<const N: usize>(h: &mut [DualN<N>], n_man: DualN<N>) {
    let dt = DualN::<N>::constant(0.01);
    let g = DualN::<N>::constant(G);
    for c in h.iter_mut() {
        let c_wave = c.mul(g).sqrt();
        let fric = n_man.mul(n_man).div(c.powf(4.0 / 3.0));
        let upd = c_wave.mul(fric).mul(dt);
        *c = c.add(upd);
    }
}

fn time_scalar(passes: usize) -> f64 {
    let t = Instant::now();
    for k in 0..passes {
        let mut h: Vec<DualN<1>> = (0..CELLS)
            .map(|i| DualN::<1>::constant(1.0 + (i % 7) as f64 * 0.1))
            .collect();
        let n = DualN::<1>::seeded(0.03 + k as f64 * 1e-9, 0);
        for _ in 0..STEPS {
            kernel(&mut h, n);
        }
        std::hint::black_box(&h);
    }
    t.elapsed().as_secs_f64()
}

fn time_vector<const N: usize>() -> f64 {
    let t = Instant::now();
    let mut h: Vec<DualN<N>> = (0..CELLS)
        .map(|i| DualN::<N>::constant(1.0 + (i % 7) as f64 * 0.1))
        .collect();
    let n = DualN::<N>::seeded(0.03, 0);
    for _ in 0..STEPS {
        kernel(&mut h, n);
    }
    std::hint::black_box(&h);
    t.elapsed().as_secs_f64()
}

fn main() {
    println!("Vector-mode spike: {CELLS} cells x {STEPS} steps\n");

    // Warm-up.
    let _ = time_scalar(1);
    let _ = time_vector::<4>();

    let base = time_scalar(1);
    println!("  one scalar pass (N=1): {base:.4} s  — the unit below\n");
    println!("{:>4}  {:>12}  {:>12}  {:>9}  {:>16}", "N", "scalar xN", "vector", "saving", "cost/component");

    for n in [2usize, 4, 8, 16] {
        let sc = time_scalar(n);
        let ve = match n {
            2 => time_vector::<2>(),
            4 => time_vector::<4>(),
            8 => time_vector::<8>(),
            _ => time_vector::<16>(),
        };
        // Cost of one further derivative component, in units of the
        // scalar pass: (vector - shared) / N, with shared taken as the
        // value-only part of one pass.
        let per = (ve / base - 1.0) / n as f64;
        println!(
            "{n:>4}  {:>11.2}x  {:>11.2}x  {:>8.0}%  {per:>15.2}x",
            sc / base,
            ve / base,
            100.0 * (1.0 - ve / sc)
        );
    }
    println!(
        "\n  break-even vs reverse mode (3-5x primal) follows from cost/component:\n  \
         N* = (k - 1) / cost_per_component"
    );
}
