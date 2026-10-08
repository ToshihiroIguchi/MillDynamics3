//! Still-pool probe: a slurry pool in a stopped drum (no balls) should come to rest. Prints the
//! fluid's RMS and maximum speed over time, so a decaying slosh (initial lattice collapse) can be
//! told apart from persistent particle noise.
//! `cargo run --release -p mill-core --example still_pool_probe -- [--mu 0.001] [--res 50]
//!  [--substeps 16] [--t 10] [--st 0.072] [--wet 0.6] [--iters 3]`

use mill_core::params::SpeedMode;
use mill_core::{Params, Simulation};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |key: &str, default: f32| -> f32 {
        args.iter()
            .position(|a| a == key)
            .and_then(|i| args.get(i + 1))
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    };
    let mut p = Params::default();
    p.mill.speed_mode = SpeedMode::Rpm;
    p.mill.speed_value = 0.0;
    p.media.fill_fraction = 0.0;
    p.slurry.viscosity_pa_s = get("--mu", 0.001);
    p.slurry.surface_tension_n_m = get("--st", 0.072);
    p.slurry.wettability = get("--wet", 0.6);
    p.simulation.resolution = get("--res", 50.0) as u32;
    p.simulation.substeps = get("--substeps", 16.0) as u32;
    p.simulation.pbf_iterations = get("--iters", p.simulation.pbf_iterations as f32) as u32;
    let t_end = get("--t", 10.0) as f64;
    let mut sim = Simulation::new(p).expect("valid params");
    let n = sim.fluid().len();
    println!(
        "still_pool_probe mu {} res {} substeps {} particles {n}",
        get("--mu", 0.001),
        get("--res", 50.0),
        get("--substeps", 16.0)
    );
    println!("{:>6} {:>10} {:>10}", "t", "v_rms", "v_max");
    let mut next = 0.0;
    while sim.sim_time() < t_end {
        sim.step_fixed();
        if sim.sim_time() >= next {
            next += 0.5;
            let v = &sim.fluid().v;
            let rms = (v.iter().map(|v| v.length_squared()).sum::<f32>() / n as f32).sqrt();
            let vmax = v.iter().map(|v| v.length()).fold(0.0f32, f32::max);
            println!("{:>6.2} {:>10.5} {:>10.5}", sim.sim_time(), rms, vmax);
        }
    }
}
