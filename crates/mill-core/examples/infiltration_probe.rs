//! Diagnostic for a user report: after lowering `slurry.viscosity_pa_s` and stopping the drum
//! (`mill.speed_value = 0`), the slurry sitting on top of the settled ball bed does not seep down
//! into the bed's interstitial gaps, even given a long time to do so. Not part of the test suite:
//! a standalone measurement of the existing model, run e.g.
//!
//! ```text
//! cargo run -p mill-core --release --example infiltration_probe -- --max-balls 150 --resolution 15 --viscosity 1
//! ```
//!
//! Spins the charge up briefly to build a realistic packed/mixed bed, then stops the drum and
//! tracks how many fluid particles sit meaningfully below the bed's own top surface over a long
//! idle period. Also reports the fluid lattice spacing `dx` against the (possibly coarse-grained)
//! effective ball diameter `d_eff`, since `dx` vs. `d_eff` is the quantity that determines whether
//! the fluid discretisation can even geometrically represent flow through gaps between balls.
//!
//! This binary makes no solver changes; it only measures the existing model. Run with `--help`
//! for every flag.

use std::env;

use mill_core::params::{Direction, SpeedMode};
use mill_core::{Params, Simulation};

struct Args {
    max_balls: u32,
    resolution: u32,
    viscosity_pa_s: f32,
    fill_fraction: f32,
    spin_up_percent_critical: f32,
    spin_up_s: f32,
    idle_s: f32,
    sample_every_s: f32,
    seed: u64,
}

impl Args {
    fn parse() -> Self {
        let mut a = Args {
            max_balls: 150,
            resolution: 15,
            viscosity_pa_s: 1.0,
            fill_fraction: 0.35,
            spin_up_percent_critical: 30.0,
            spin_up_s: 8.0,
            idle_s: 20.0,
            sample_every_s: 2.0,
            seed: 1,
        };
        let mut args = env::args().skip(1);
        while let Some(flag) = args.next() {
            let mut next_f32 = || -> f32 {
                args.next()
                    .expect("flag needs a value")
                    .parse()
                    .expect("expected a number")
            };
            match flag.as_str() {
                "--max-balls" => a.max_balls = next_f32() as u32,
                "--resolution" => a.resolution = next_f32() as u32,
                "--viscosity" => a.viscosity_pa_s = next_f32(),
                "--fill" => a.fill_fraction = next_f32(),
                "--spin-up-percent-critical" => a.spin_up_percent_critical = next_f32(),
                "--spin-up-s" => a.spin_up_s = next_f32(),
                "--idle-s" => a.idle_s = next_f32(),
                "--sample-every-s" => a.sample_every_s = next_f32(),
                "--seed" => a.seed = next_f32() as u64,
                "--help" | "-h" => {
                    println!(
                        "infiltration_probe: measures whether slurry seeps into a settled ball \
                         bed once the drum stops.\n\n\
                         Flags: --max-balls --resolution --viscosity --fill \
                         --spin-up-percent-critical --spin-up-s --idle-s --sample-every-s --seed"
                    );
                    std::process::exit(0);
                }
                other => panic!("unknown flag: {other}"),
            }
        }
        a
    }
}

fn main() {
    let args = Args::parse();

    let mut params = Params::default();
    params.simulation.max_balls = args.max_balls;
    params.simulation.resolution = args.resolution;
    params.simulation.seed = args.seed;
    params.slurry.enabled = true;
    params.slurry.viscosity_pa_s = args.viscosity_pa_s;
    params.slurry.fill_fraction = args.fill_fraction;
    params.mill.speed_mode = SpeedMode::PercentCritical;
    params.mill.speed_value = args.spin_up_percent_critical;
    params.mill.direction = Direction::CounterClockwise;
    params.lifters.count = 0;
    params.validate().expect("invalid params");

    let effective = params.effective_media();
    let dx = params.mill.radius_m() / params.effective_fluid_resolution() as f32;
    println!(
        "config: max_balls={} resolution={} (effective {}) viscosity={} Pa*s, fill={}",
        args.max_balls,
        args.resolution,
        params.effective_fluid_resolution(),
        args.viscosity_pa_s,
        args.fill_fraction,
    );
    println!(
        "effective ball diameter d_eff={:.4} m (true {:.4} m, scale x{:.2}), fluid lattice spacing dx={:.4} m, dx/d_eff={:.3}",
        effective.diameter_m,
        effective.true_diameter_m,
        effective.scale_factor,
        dx,
        dx / effective.diameter_m,
    );

    let mut sim = Simulation::new(params).expect("failed to build simulation");
    let dt = 1.0 / 60.0;
    let ball_radius = sim.balls().radius;

    let report = |sim: &Simulation, label: &str| {
        let balls = sim.balls();
        let fluid = sim.fluid();
        if balls.is_empty() {
            println!("{label}: no balls");
            return;
        }
        let bed_top_y = balls
            .x
            .iter()
            .fold(f32::NEG_INFINITY, |m, p| m.max(p.y + ball_radius));
        let bed_bottom_y = balls
            .x
            .iter()
            .fold(f32::INFINITY, |m, p| m.min(p.y - ball_radius));
        // "Embedded": at least 1.5 ball diameters below the bed's own topmost point -- a coarse
        // bulk-position check (naive: the pack tilts/slopes, so this alone conflates "beside the
        // slope" with "actually under a ball").
        let embed_threshold = bed_top_y - 3.0 * ball_radius;
        let total = fluid.len();
        let embedded = fluid.x.iter().filter(|p| p.y < embed_threshold).count();
        // "Covered": there exists a ball horizontally overlapping this particle (|dx| < radius)
        // whose own bottom surface sits above the particle -- i.e. this particle has a ball
        // directly over it, meaning it is physically underneath/between balls rather than sitting
        // in the open pool with a clear line to the surface above. This is the metric that
        // actually matches "has slurry gotten into the packed bed", independent of the bed's slope.
        let covered = fluid
            .x
            .iter()
            .filter(|p| {
                balls
                    .x
                    .iter()
                    .any(|b| (b.x - p.x).abs() < ball_radius && (b.y - ball_radius) > p.y)
            })
            .count();
        println!(
            "{label}: bed y in [{bed_bottom_y:.3}, {bed_top_y:.3}] m | fluid particles: {total} total, \
             {embedded} ({:.1}%) >1.5 diam below bed top, {covered} ({:.1}%) covered (a ball sits directly over them)",
            100.0 * embedded as f32 / total.max(1) as f32,
            100.0 * covered as f32 / total.max(1) as f32,
        );
    };

    report(&sim, "t=0.0s (seeded)");

    let mut t = 0.0f32;
    while t < args.spin_up_s {
        sim.step(dt);
        t += dt;
    }
    report(&sim, &format!("t={t:.1}s (end of spin-up)"));

    // Stop the drum and hold it there for the rest of the run.
    let mut stopped_params = params;
    stopped_params.mill.speed_value = 0.0;
    sim.set_params(stopped_params).expect("valid params");

    let mut next_sample = args.sample_every_s;
    let idle_end = t + args.idle_s;
    while t < idle_end {
        sim.step(dt);
        t += dt;
        if t >= next_sample {
            report(&sim, &format!("t={t:.1}s (idle, drum stopped)"));
            next_sample += args.sample_every_s;
        }
    }
    report(&sim, &format!("t={t:.1}s (final)"));
}
