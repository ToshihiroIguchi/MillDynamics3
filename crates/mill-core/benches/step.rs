//! Benchmarks for `Simulation::step` and the underlying ball solver, at the ball-count sweep
//! called for in docs/PLAN.md ss5 (M1) and ss6.

use criterion::{criterion_group, criterion_main, Criterion};
use mill_core::dem::DemState;
use mill_core::geometry::Drum;
use mill_core::params::{LiftersParams, MediaParams};
use mill_core::{EffectiveMedia, Params, Simulation, FIXED_DT};

const DRUM_RADIUS_M: f32 = 0.5;
const FILL_FRACTION: f32 = 0.30;

/// Builds an effective media population with exactly `ball_count` balls at `FILL_FRACTION` of the
/// benchmark drum's area, bypassing `Params::effective_media`'s coarse-graining derivation so the
/// sweep hits exact, round ball counts (500/1000/2000) regardless of `max_balls`.
fn effective_media_with_count(ball_count: u32) -> EffectiveMedia {
    let drum_area = std::f32::consts::PI * DRUM_RADIUS_M * DRUM_RADIUS_M;
    let r = (FILL_FRACTION * drum_area / (std::f32::consts::PI * ball_count as f32)).sqrt();
    EffectiveMedia {
        true_diameter_m: 2.0 * r,
        diameter_m: 2.0 * r,
        density_kg_m3: 7800.0,
        ball_count,
        scale_factor: 1.0,
    }
}

fn bench_ball_solver_step(c: &mut Criterion) {
    let media = MediaParams::default();
    let drum = Drum::new(
        DRUM_RADIUS_M,
        3.0,
        LiftersParams {
            count: 0,
            ..LiftersParams::default()
        },
    );

    let mut group = c.benchmark_group("dem_step");
    for &ball_count in &[500u32, 1000, 2000] {
        let effective = effective_media_with_count(ball_count);
        let mut state = DemState::new(&effective, DRUM_RADIUS_M, 1);
        // Warm up so the timed iterations measure a representative settled/moving arrangement,
        // not the initial perfectly-ordered lattice.
        for _ in 0..60 {
            state.step(&drum, 0.0, &media, 4, FIXED_DT);
        }

        group.bench_function(format!("{ball_count}_balls"), |b| {
            let mut drum_angle = 0.0f32;
            b.iter(|| {
                state.step(&drum, drum_angle, &media, 4, FIXED_DT);
                drum_angle += drum.omega * FIXED_DT;
            });
        });
    }
    group.finish();
}

fn bench_drum_only_step(c: &mut Criterion) {
    let mut sim = Simulation::new(Params::default()).expect("default params are valid");
    c.bench_function("drum_only_step", |b| {
        b.iter(|| sim.step(FIXED_DT));
    });
}

criterion_group!(benches, bench_ball_solver_step, bench_drum_only_step);
criterion_main!(benches);
