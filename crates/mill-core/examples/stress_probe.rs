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
    balls: usize,
    k: f32,
}

fn run(substeps: u32, iters: u32, debounce: Option<f64>, seed: u64) -> Row {
    run_with(substeps, iters, debounce, seed, |_| {})
}

fn run_with(
    substeps: u32,
    iters: u32,
    debounce: Option<f64>,
    seed: u64,
    tweak: impl Fn(&mut Params),
) -> Row {
    let mut p = Params::default();
    p.slurry.enabled = false;
    tweak(&mut p);
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
    let (balls, k) = {
        let m = sim.metrics();
        (m.simulated_ball_count as usize, m.coarse_graining_factor)
    };
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
        balls,
        k,
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
        let n = rows.len() as f64;
        let col = |g: &dyn Fn(&Row) -> f64| rows.iter().map(g).collect::<Vec<f64>>();
        let cols = [
            col(&|r| r.impacts),
            col(&|r| r.se_drop),
            col(&|r| r.se_shear),
            col(&|r| r.intensity),
            col(&|r| r.spec_power),
        ];
        let se = |k: usize| {
            if rows.len() < 2 {
                0.0
            } else {
                let v = cols[k].iter().map(|x| (x - m[k]).powi(2)).sum::<f64>() / (n - 1.0);
                (v / n).sqrt()
            }
        };
        println!(
            "| {label} (balls {}, k {:.2}) | {:.1} ± {:.1} | {:.4} ± {:.4} | {:.4} ± {:.4} | {:.4} ± {:.4} | {:.4} ± {:.4} |",
            rows[0].balls, rows[0].k, m[0], se(0), m[1], se(1), m[2], se(2), m[3], se(3), m[4], se(4)
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
    if cases.contains("dd") {
        use mill_core::params::SpeedMode;
        println!(
            "## dd: D/d sweep (2 mm ball, 16 sub-steps, 12 iterations, seed 1)
"
        );
        println!("| %Nc, D/d {hdr}");
        let pcs: Vec<f32> = arg("--pcs")
            .map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect())
            .unwrap_or_else(|| vec![59.0, 85.0]);
        let ratios: Vec<f32> = arg("--ratios")
            .map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect())
            .unwrap_or_else(|| vec![20.0, 30.0, 45.0, 60.0]);
        for &pc in &pcs {
            for &ratio in &ratios {
                let rows: Vec<Row> = (1..=n_seeds)
                    .map(|seed| {
                        run_with(16, 12, None, seed, |p| {
                            p.simulation.coarse_graining_mode =
                                mill_core::params::CoarseGrainingMode::Off;
                            p.mill.diameter_m = ratio * p.media.ball_diameter_m;
                            p.mill.speed_mode = SpeedMode::PercentCritical;
                            p.mill.speed_value = pc;
                        })
                    })
                    .collect();
                line(format!("{pc}, {ratio}"), rows);
            }
        }
        println!();
    }
    if cases.contains("cg") {
        use mill_core::params::{CoarseGrainingMode, SpeedMode};
        println!(
            "## cg: coarse-graining factor k at D/d = 45 (59 % Nc, 16 sub-steps, 12 iterations)
"
        );
        println!("| k | balls | power W/m | dissipated W/m | KE J/m | toe deg | shoulder deg | centroid y mm | intensity | shear SE | drop SE |");
        println!("|---|---|---|---|---|---|---|---|---|---|---|");
        for k in [1.0f32, 1.5, 2.0, 3.0] {
            let mut acc = [0.0f64; 9];
            let mut balls = 0usize;
            for seed in 1..=n_seeds {
                let mut p = Params::default();
                p.slurry.enabled = false;
                p.mill.diameter_m = 45.0 * p.media.ball_diameter_m;
                p.mill.speed_mode = SpeedMode::PercentCritical;
                p.mill.speed_value = 59.0;
                p.simulation.substeps = 16;
                p.simulation.dem_iterations = 12;
                p.simulation.seed = seed;
                p.simulation.coarse_graining_mode = CoarseGrainingMode::Manual;
                p.simulation.coarse_graining_k = k;
                let rev_s = 60.0 / p.mill.rpm() as f64;
                let omega = p.mill.omega();
                let mut sim = Simulation::new(p).unwrap();
                for _ in 0..(4.0 * rev_s * 60.0) as usize {
                    sim.step(1.0 / 60.0);
                }
                let frames = (4.0 * rev_s * 60.0) as usize;
                let mut a = [0.0f64; 9];
                let mut cnt_toe: f64 = 0.0;
                for _ in 0..frames {
                    sim.step(1.0 / 60.0);
                    let m = sim.metrics();
                    balls = m.simulated_ball_count as usize;
                    a[0] += m.power_draw_w as f64;
                    a[1] += m.dissipated_power_w as f64;
                    a[2] += m.total_kinetic_energy_j as f64;
                    if let (Some(t), Some(sh)) = (m.toe_angle_rad, m.shoulder_angle_rad) {
                        a[3] += mill_core::metrics::to_vertical_degrees(t, omega) as f64;
                        a[4] += mill_core::metrics::to_vertical_degrees(sh, omega) as f64;
                        cnt_toe += 1.0;
                    }
                    a[5] += m
                        .charge_centroid_m
                        .map(|c| c.1 as f64 * 1000.0)
                        .unwrap_or(0.0);
                    a[6] += m.stress.stress_intensity_scaled as f64;
                    a[7] += m.stress.mean_se_scaled_shear as f64;
                    a[8] += m.stress.mean_se_scaled_drop as f64;
                }
                let n = frames as f64;
                for (i, v) in a.iter().enumerate() {
                    let d = if i == 3 || i == 4 {
                        cnt_toe.max(1.0)
                    } else {
                        n
                    };
                    acc[i] += v / d;
                }
            }
            let ns = n_seeds as f64;
            println!(
                "| {k} | {balls} | {:.3} | {:.3} | {:.4} | {:.1} | {:.1} | {:.2} | {:.3} | {:.3} | {:.3} |",
                acc[0] / ns, acc[1] / ns, acc[2] / ns, acc[3] / ns, acc[4] / ns, acc[5] / ns,
                acc[6] / ns, acc[7] / ns, acc[8] / ns
            );
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
