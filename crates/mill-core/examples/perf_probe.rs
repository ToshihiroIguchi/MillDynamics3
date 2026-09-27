//! Performance-regression probe for the M6 solver pass (docs/PERF.md).
//!
//! Not part of the test suite: a standalone diagnostic run, e.g.
//!
//! ```text
//! cargo run -p mill-core --release --example perf_probe -- --config accuracy --frames 300
//! cargo run -p mill-core --release --example perf_probe -- --config all
//! ```
//!
//! For each named configuration this settles the charge for `--settle-s` seconds, then times
//! `--frames` calls to `Simulation::step(1/60)`, reporting:
//! - mean wall-clock ms/frame,
//! - mean fluid solver diagnostics (viscosity CG iterations, mean shear rate) over the measured
//!   window,
//! - an FNV-1a hash over the raw bits of every ball/fluid `x`/`v`/`omega` component after the
//!   measured window.
//!
//! The hash is this project's bit-exactness oracle while restructuring the solver's internals
//! (`grid.rs`/`pbf.rs`/`dem.rs`) for speed without changing any floating-point summation order:
//! an optimization pass that is meant to be purely mechanical must reproduce the exact same hash,
//! for every configuration below, both before and after the change (alongside `cargo test`
//! passing) -- see docs/PERF.md's 2026-09-27 "M6 solver pass" section for how this was used.

use std::env;
use std::time::Instant;

use mill_core::params::{Direction, LiftersParams, MediaParams, SlurryParams, SpeedMode};
use mill_core::{Params, Simulation};

struct Config {
    name: &'static str,
    max_balls: u32,
    resolution: u32,
    slurry_enabled: bool,
    lifters_count: u32,
}

const CONFIGS: [Config; 5] = [
    Config {
        name: "realtime",
        max_balls: 150,
        resolution: 15,
        slurry_enabled: true,
        lifters_count: 0,
    },
    Config {
        name: "balanced",
        max_balls: 300,
        resolution: 25,
        slurry_enabled: true,
        lifters_count: 0,
    },
    Config {
        name: "accuracy",
        max_balls: 600,
        resolution: 40,
        slurry_enabled: true,
        lifters_count: 0,
    },
    Config {
        name: "dry",
        max_balls: 600,
        resolution: 40,
        slurry_enabled: false,
        lifters_count: 0,
    },
    Config {
        name: "lifters",
        max_balls: 600,
        resolution: 40,
        slurry_enabled: true,
        lifters_count: 4,
    },
];

struct Args {
    config: String,
    settle_s: f32,
    frames: u32,
    seed: u64,
}

impl Args {
    fn parse() -> Self {
        let mut a = Args {
            config: "all".to_string(),
            settle_s: 2.0,
            frames: 300,
            seed: 1,
        };
        let mut args = env::args().skip(1);
        while let Some(flag) = args.next() {
            match flag.as_str() {
                "--config" => a.config = args.next().expect("--config needs a value"),
                "--settle-s" => {
                    a.settle_s = args
                        .next()
                        .expect("--settle-s needs a value")
                        .parse()
                        .expect("expected a number")
                }
                "--frames" => {
                    a.frames = args
                        .next()
                        .expect("--frames needs a value")
                        .parse()
                        .expect("expected an integer")
                }
                "--seed" => {
                    a.seed = args
                        .next()
                        .expect("--seed needs a value")
                        .parse()
                        .expect("expected an integer")
                }
                "--help" | "-h" => {
                    println!(
                        "Flags: --config <realtime|balanced|accuracy|dry|lifters|all> \
                         --settle-s <f32> --frames <u32> --seed <u64>"
                    );
                    std::process::exit(0);
                }
                other => panic!("unknown flag: {other}"),
            }
        }
        a
    }
}

/// FNV-1a over the raw bits of every ball/fluid state component, in the same fixed field/index
/// order every run -- deterministic run-to-run for byte-identical solver behaviour (this crate's
/// own "reproducible from Params + seed" convention, see `grid.rs`'s module doc comment), and
/// sensitive to any change in floating-point summation order, not just to a change in final
/// physical outcome.
fn state_hash(sim: &Simulation) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;
    let mut h = FNV_OFFSET;
    let mut mix = |v: f32| {
        for b in v.to_bits().to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(FNV_PRIME);
        }
    };
    let balls = sim.balls();
    for i in 0..balls.len() {
        mix(balls.x[i].x);
        mix(balls.x[i].y);
        mix(balls.v[i].x);
        mix(balls.v[i].y);
        mix(balls.omega[i]);
    }
    let fluid = sim.fluid();
    for i in 0..fluid.x.len() {
        mix(fluid.x[i].x);
        mix(fluid.x[i].y);
        mix(fluid.v[i].x);
        mix(fluid.v[i].y);
    }
    h
}

fn run_config(cfg: &Config, args: &Args) {
    let mut params = Params {
        mill: mill_core::params::MillParams {
            diameter_m: 1.0,
            speed_mode: SpeedMode::PercentCritical,
            speed_value: 70.0,
            direction: Direction::CounterClockwise,
        },
        media: MediaParams::default(),
        lifters: LiftersParams {
            count: cfg.lifters_count,
            ..LiftersParams::default()
        },
        slurry: SlurryParams {
            enabled: cfg.slurry_enabled,
            ..SlurryParams::default()
        },
        ..Params::default()
    };
    params.simulation.max_balls = cfg.max_balls;
    params.simulation.resolution = cfg.resolution;
    params.simulation.seed = args.seed;
    params.validate().expect("invalid params");

    let mut sim = Simulation::new(params).expect("failed to build simulation");

    let settle_frames = (args.settle_s * 60.0).round() as u32;
    for _ in 0..settle_frames {
        sim.step(1.0 / 60.0);
    }

    let mut viscosity_iters_sum = 0u64;
    let mut shear_rate_sum = 0.0f64;
    let start = Instant::now();
    for _ in 0..args.frames {
        sim.step(1.0 / 60.0);
        viscosity_iters_sum += sim.viscosity_iterations() as u64;
        shear_rate_sum += sim.mean_shear_rate_per_s() as f64;
    }
    let elapsed = start.elapsed();
    let ms_per_frame = elapsed.as_secs_f64() * 1000.0 / args.frames as f64;
    let mean_viscosity_iters = viscosity_iters_sum as f64 / args.frames as f64;
    let mean_shear_rate = shear_rate_sum / args.frames as f64;

    let hash = state_hash(&sim);

    println!(
        "config={:<9} balls={:<4} fluid={:<5} ms/frame={:8.4} mean_cg_iters={:6.2} \
         mean_shear_rate={:9.3} hash=0x{:016x}",
        cfg.name,
        sim.balls().len(),
        sim.fluid().len(),
        ms_per_frame,
        mean_viscosity_iters,
        mean_shear_rate,
        hash,
    );
}

fn main() {
    let args = Args::parse();
    println!(
        "perf_probe: settle_s={} frames={} seed={}",
        args.settle_s, args.frames, args.seed
    );
    for cfg in &CONFIGS {
        if args.config == "all" || args.config == cfg.name {
            run_config(cfg, &args);
        }
    }
}
