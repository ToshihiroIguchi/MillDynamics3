//! Benchmarks for `Simulation::step`.
//!
//! Through M0 this only measures drum-kinematics-only stepping. M1 replaces/extends this with the
//! 500/1000/2000-ball DEM sweep described in docs/PLAN.md ss5 (M1) and ss6.

use criterion::{criterion_group, criterion_main, Criterion};
use mill_core::{Params, Simulation};

fn bench_drum_only_step(c: &mut Criterion) {
    let mut sim = Simulation::new(Params::default()).expect("default params are valid");
    c.bench_function("drum_only_step", |b| {
        b.iter(|| sim.step(mill_core::FIXED_DT));
    });
}

criterion_group!(benches, bench_drum_only_step);
criterion_main!(benches);
