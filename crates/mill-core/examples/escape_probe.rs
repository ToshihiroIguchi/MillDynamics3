//! Containment probe: runs the default configuration (optionally with overrides) and reports, per
//! sampled second, how far the outermost ball pokes through the drum wall and how fast balls get.
//! `cargo run --release -p mill-core --example escape_probe -- [--t 6] [--no-slurry] [--trace]`

use mill_core::{Params, Simulation};

#[allow(clippy::needless_range_loop)]
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let flag = |name: &str| args.iter().any(|a| a == name);
    let value = |name: &str, default: f32| -> f32 {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    };
    let mut params = Params::default();
    if flag("--no-slurry") {
        params.slurry.enabled = false;
    }
    let t_end = value("--t", 6.0) as f64;
    let mut sim = Simulation::new(params).expect("valid params");
    let radius = sim.params().mill.radius_m();
    let mut next = 0.0;
    println!(
        "{:>6} {:>8} {:>10} {:>10} {:>8}",
        "t", "balls", "max out", "max speed", "outside"
    );
    while sim.sim_time() < t_end {
        let before: Vec<(glam::Vec2, glam::Vec2)> = sim
            .balls()
            .x
            .iter()
            .zip(sim.balls().v.iter())
            .map(|(a, b)| (*a, *b))
            .collect();
        sim.step_fixed();
        if flag("--trace") {
            let b = sim.balls();
            for i in 0..b.x.len() {
                let o = b.x[i].length() + b.radius - radius;
                if o > 2e-3 {
                    println!(
                        "t {:.4} ball {i} out {o:.2e} | x {:?} -> {:?} | v {:?} -> {:?} | hits {}",
                        sim.sim_time(),
                        before[i].0,
                        b.x[i],
                        before[i].1,
                        b.v[i],
                        sim.coupling_clamp_hits()
                    );
                }
            }
        }
        if sim.sim_time() >= next {
            next += 0.25;
            let b = sim.balls();
            let mut worst = f32::NEG_INFINITY;
            let mut vmax = 0.0f32;
            let mut outside = 0;
            for i in 0..b.x.len() {
                let o = b.x[i].length() + b.radius - radius;
                worst = worst.max(o);
                vmax = vmax.max(b.v[i].length());
                if b.x[i].length() > radius {
                    outside += 1;
                }
            }
            println!(
                "{:>6.2} {:>8} {:>10.2e} {:>10.3} {:>8}",
                sim.sim_time(),
                b.x.len(),
                worst,
                vmax,
                outside
            );
        }
    }
}
