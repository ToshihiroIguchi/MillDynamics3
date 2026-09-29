//! Prints the per-frame distribution of the fluid compression-error metrics at `Params::default()`.
//! `cargo run --release -p mill-core --example compression_probe -- [warmup_s] [measure_s]`

use mill_core::{Params, Simulation};

fn pct(sorted: &[f32], p: f32) -> f32 {
    let i = ((sorted.len() - 1) as f32 * p).round() as usize;
    sorted[i]
}

fn frac_above(v: &[f32], t: f32) -> f32 {
    v.iter().filter(|x| **x > t).count() as f32 / v.len() as f32
}

fn main() {
    let args: Vec<f32> = std::env::args()
        .skip(1)
        .filter_map(|a| a.parse().ok())
        .collect();
    let warmup_s = args.first().copied().unwrap_or(2.0);
    let measure_s = args.get(1).copied().unwrap_or(12.0);
    let mut params = Params::default();
    if let Ok(v) = std::env::var("PROBE_DIAMETER_M") {
        params.mill.diameter_m = v.parse().expect("PROBE_DIAMETER_M");
    }
    if let Ok(v) = std::env::var("PROBE_RPM") {
        params.mill.speed_value = v.parse().expect("PROBE_RPM");
    }
    if let Ok(v) = std::env::var("PROBE_MEDIA_FILL") {
        params.media.fill_fraction = v.parse().expect("PROBE_MEDIA_FILL");
    }
    if let Ok(v) = std::env::var("PROBE_SLURRY_FILL") {
        params.slurry.fill_fraction = v.parse().expect("PROBE_SLURRY_FILL");
    }
    if let Ok(v) = std::env::var("PROBE_SURFACE_TENSION") {
        params.slurry.surface_tension_n_m = v.parse().expect("PROBE_SURFACE_TENSION");
    }
    if let Ok(v) = std::env::var("PROBE_WETTABILITY") {
        params.slurry.wettability = v.parse().expect("PROBE_WETTABILITY");
    }
    if let Ok(v) = std::env::var("PROBE_MAX_BALLS") {
        params.simulation.max_balls = v.parse().expect("PROBE_MAX_BALLS");
    }
    println!(
        "max_balls={} slurry={} rpm-mode value={}",
        params.simulation.max_balls, params.slurry.enabled, params.mill.speed_value
    );
    let mut sim = Simulation::new(params).expect("valid params");
    let dt = 1.0 / 60.0;
    let m0 = sim.metrics();
    println!(
        "t=0: max={:?} mean={:?}",
        m0.max_fluid_compression_error_fraction, m0.mean_fluid_compression_error_fraction
    );
    for _ in 0..(warmup_s * 60.0) as usize {
        sim.step(dt);
    }
    let (mut maxv, mut meanv) = (Vec::new(), Vec::new());
    for _ in 0..(measure_s * 60.0) as usize {
        sim.step(dt);
        let m = sim.metrics();
        maxv.push(m.max_fluid_compression_error_fraction.unwrap_or(0.0));
        meanv.push(m.mean_fluid_compression_error_fraction.unwrap_or(0.0));
    }
    // 1 s rolling median (60 frames) of max-error.
    let mut roll = Vec::new();
    for w in maxv.windows(60) {
        let mut s = w.to_vec();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roll.push(s[30]);
    }
    let mut ms = maxv.clone();
    ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut ns = meanv.clone();
    ns.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!("frames={}", maxv.len());
    println!(
        "max-error p50={:.4} p95={:.4} p99={:.4} max={:.4}",
        pct(&ms, 0.5),
        pct(&ms, 0.95),
        pct(&ms, 0.99),
        ms[ms.len() - 1]
    );
    println!(
        "mean-error mean={:.4} p95={:.4}",
        meanv.iter().sum::<f32>() / meanv.len() as f32,
        pct(&ns, 0.95)
    );
    for t in [0.05, 0.10, 0.20] {
        println!(
            "frac frames max>{:.0}% = {:.3}",
            t * 100.0,
            frac_above(&maxv, t)
        );
    }
    println!(
        "frac rolling-median(max)>20% = {:.3}",
        frac_above(&roll, 0.20)
    );
}
