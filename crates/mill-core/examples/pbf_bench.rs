//! PBF side of the grid-vs-PBF benchmark (docs/VERIFICATION.md): a small heavy disc charge in a
//! slurry-filled drum, reporting the net (buoyancy-corrected) gravity moment of the discs about
//! the drum axis -- the same metric `grid_probe --exp E8a` prints. Dimensionless grid units map
//! to SI as drum radius 0.5 m, fluid density 1000 kg/m^3, mu = nu * 1000.
//! `cargo run --release -p mill-core --example pbf_bench -- [--nu 0.1] [--rho 8] [--omega 1]
//!  [--a 0.04] [--fill 0.077] [--t 4] [--res 40] [--substeps 8] [--k 1] [--tau 0]
//!  [--match 1 --m 10]` (`--match`: grid-identical start)

use mill_core::params::{CoarseGrainingMode, Rheology, SpeedMode};
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
    let (nu, rho, omega, a) = (
        get("--nu", 0.1),
        get("--rho", 8.0),
        get("--omega", 1.0),
        get("--a", 0.04),
    );
    let tau = get("--tau", 0.0);
    let rho_f = 1000.0f32;
    let mut p = Params::default();
    p.mill.diameter_m = 1.0;
    p.mill.speed_mode = SpeedMode::Rpm;
    p.mill.speed_value = omega * 60.0 / std::f32::consts::TAU;
    p.media.ball_diameter_m = 2.0 * a;
    p.media.density_kg_m3 = rho * rho_f;
    p.media.fill_fraction = get("--fill", 0.077);
    p.slurry.enabled = true;
    p.slurry.fill_fraction = get("--sfill", 0.9);
    p.slurry.density_kg_m3 = rho_f;
    p.slurry.viscosity_pa_s = nu * rho_f;
    if tau > 0.0 {
        p.slurry.rheology = Rheology::Bingham;
        p.slurry.yield_stress_pa = tau * rho_f;
    }
    p.simulation.resolution = get("--res", 40.0) as u32;
    p.simulation.substeps = get("--substeps", 8.0) as u32;
    p.simulation.max_balls = 10_000;
    let k = get("--k", 1.0);
    if k > 1.0 {
        p.simulation.coarse_graining_mode = CoarseGrainingMode::Manual;
        p.simulation.coarse_graining_k = k;
    }
    let t_end = get("--t", 4.0) as f64;
    let mut sim = if get("--match", 0.0) > 0.0 {
        // Same start as `grid_probe --exp E8a`: hex block of `--m` discs, solid-body rotation.
        let count = get("--m", 10.0) as usize;
        let pitch = 2.4 * a;
        let mut sites: Vec<(f32, f32)> = Vec::new();
        for j in -5i32..=5 {
            for i in -5i32..=5 {
                sites.push((
                    pitch * (i as f32 + 0.5 * (j & 1) as f32) + 0.0137,
                    pitch * 0.866_025_4 * j as f32 - 0.0091 - 0.25,
                ));
            }
        }
        sites.sort_by(|p, q| {
            (p.0 * p.0 + (p.1 + 0.25).powi(2)).total_cmp(&(q.0 * q.0 + (q.1 + 0.25).powi(2)))
        });
        sites.truncate(count);
        let pos: Vec<glam::Vec2> = sites.iter().map(|&(x, y)| glam::Vec2::new(x, y)).collect();
        Simulation::with_initial_balls(p, &pos).expect("valid params")
    } else {
        Simulation::new(p).expect("valid params")
    };
    let b = sim.balls();
    let (n, r, m) = (b.x.len(), b.radius, b.mass);
    let force = (m - rho_f * std::f32::consts::PI * r * r) * 9.81;
    println!(
        "pbf_bench nu {nu} rho_s {rho} omega {omega} balls {n} radius {r:.4} k {k} res {} substeps {}",
        sim.params().simulation.resolution,
        sim.params().simulation.substeps
    );
    println!("{:>7} {:>12} {:>10}", "t", "gravity mom", "P draw");
    let (mut next, mut peak, mut sum, mut cnt) = (0.0, 0.0f64, 0.0f64, 0usize);
    while sim.sim_time() < t_end {
        sim.step_fixed();
        let tg: f64 = sim.balls().x.iter().map(|x| -(force * x.x) as f64).sum();
        if tg.abs() > peak.abs() {
            peak = tg;
        }
        if sim.sim_time() > 0.5 * t_end {
            sum += tg;
            cnt += 1;
        }
        if sim.sim_time() >= next {
            next += 0.25;
            println!(
                "{:>7.2} {:>12.5e} {:>10.2}",
                sim.sim_time(),
                tg,
                sim.power_draw_w()
            );
        }
    }
    println!(
        "peak gravity moment {peak:.5e}, mean over 2nd half {:.5e}",
        sum / cnt.max(1) as f64
    );
}
