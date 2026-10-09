//! Convergence probe for the 3D-converted impact statistics of the main simulation
//! (`mill_core::stress`): how they depend on sub-steps, contact iterations and the
//! impact-episode debounce, dry 63 mm drum / 2 mm ZrO2 at the shipped defaults.
//!
//! `cargo run --release -p mill-core --example stress_probe -- [--seeds 2] [--cases dt,iters]`
//! Prints markdown tables (mean over the measurement window, then over seeds).

use mill_core::{Params, Simulation};

struct Row {
    impacts: f64,
    se_drop: f64,
    se_shear: f64,
    intensity: f64,
    spec_power: f64,
}

fn run(substeps: u32, iters: u32, debounce: Option<f64>, seed: u64) -> Row {
    let mut p = Params::default();
    p.slurry.enabled = false;
    p.simulation.substeps = substeps;
    p.simulation.dem_iterations = iters;
    p.simulation.seed = seed;
    let rev_s = 60.0 / p.mill.rpm() as f64;
    let mut sim = Simulation::new(p).unwrap();
    if let Some(f) = debounce {
        sim.set_stress_debounce_factor(f);
    }
    let frame = 1.0f32 / 60.0;
    let warm = (4.0 * rev_s * 60.0).round() as usize;
    let meas = (4.0 * rev_s * 60.0).round() as usize;
    for _ in 0..warm {
        sim.step(frame);
    }
    let mut acc = [0.0f64; 5];
    for _ in 0..meas {
        sim.step(frame);
        let s = sim.metrics().stress;
        acc[0] += s.impacts_total_per_ball_per_s as f64;
        acc[1] += s.mean_se_scaled_drop as f64;
        acc[2] += s.mean_se_scaled_shear as f64;
        acc[3] += s.stress_intensity_scaled as f64;
        acc[4] += s.specific_power_w_per_kg as f64;
    }
    let n = meas as f64;
    Row {
        impacts: acc[0] / n,
        se_drop: acc[1] / n,
        se_shear: acc[2] / n,
        intensity: acc[3] / n,
        spec_power: acc[4] / n,
    }
}

fn mean_rows(rows: &[Row]) -> [f64; 5] {
    let n = rows.len() as f64;
    let f = |g: &dyn Fn(&Row) -> f64| rows.iter().map(g).sum::<f64>() / n;
    [
        f(&|r| r.impacts),
        f(&|r| r.se_drop),
        f(&|r| r.se_shear),
        f(&|r| r.intensity),
        f(&|r| r.spec_power),
    ]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .and_then(|i| args.get(i + 1).cloned())
    };
    let n_seeds: u64 = arg("--seeds").and_then(|s| s.parse().ok()).unwrap_or(2);
    let cases = arg("--cases").unwrap_or_else(|| "dt,iters".into());
    let hdr = "| impacts/ball/s | SE drop | SE shear | intensity | W/kg |\n|---|---|---|---|---|";
    let line = |label: String, rows: Vec<Row>| {
        let m = mean_rows(&rows);
        println!(
            "| {label} | {:.1} | {:.4} | {:.4} | {:.4} | {:.4} |",
            m[0], m[1], m[2], m[3], m[4]
        );
    };
    if cases.contains("dt") {
        println!("## dt: sub-steps per 1/60 s frame (12 iterations)\n");
        println!("| substeps, debounce {hdr}");
        for debounce in [Some(0.0), None] {
            for sub in [8u32, 16, 32, 64] {
                let rows = (1..=n_seeds).map(|s| run(sub, 12, debounce, s)).collect();
                let tag = if debounce.is_some() { "off" } else { "on" };
                line(format!("{sub}, {tag}"), rows);
            }
        }
        println!();
    }
    if cases.contains("deb") {
        println!(
            "## deb: debounce window factor x sub-steps (12 iterations)
"
        );
        println!("| factor, substeps {hdr}");
        for f in [0.25f64, 0.5, 1.0, 2.0] {
            for sub in [16u32, 32, 64] {
                let rows = (1..=n_seeds).map(|s| run(sub, 12, Some(f), s)).collect();
                line(format!("{f}, {sub}"), rows);
            }
        }
        println!();
    }
    if cases.contains("iters") {
        println!("## iters: contact iterations (32 sub-steps, debounce on)\n");
        println!("| iterations {hdr}");
        for it in [2u32, 6, 12, 20] {
            let rows = (1..=n_seeds).map(|s| run(32, it, None, s)).collect();
            line(format!("{it}"), rows);
        }
    }
}
