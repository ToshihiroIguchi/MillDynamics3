//! A/B tolerance gate for physical metrics.
//!
//! Not part of the test suite: a standalone diagnostic. Record a baseline from one build, change
//! the solver, then compare:
//!
//! ```text
//! cargo run --release -p mill-core --example tolerance_probe -- --baseline base.csv
//! cargo run --release -p mill-core --example tolerance_probe -- --compare base.csv
//! ```
//!
//! Configurations are the `ui_fast` (300 balls / res 25) and `ui_accuracy` (600 balls / res 50)
//! configs of `perf_probe` (1 m drum, 70 % Nc, CCW, wet, coarse-graining `Auto`), seeds 1..=3,
//! 2 settle revolutions then 3 measured revolutions. Runs are sequential (so ms/frame is honest)
//! and deterministic. Metrics are window means per seed, reported as mean +- standard error over
//! seeds (n - 1). Gate rule per metric (`--compare`):
//!
//! - power, torque, specific power, impacts/ball/s: PASS if
//!   `|delta| <= max(2 sqrt(se_a^2 + se_b^2), 2 % of |baseline mean|)`;
//! - toe / shoulder angle (degrees from vertical): same with an absolute 2 degree floor, delta
//!   taken modulo 360;
//! - unattributed / shaft energy fraction: `|new| <= |baseline| * 1.10 + 0.005` (the solver's
//!   energy closure is already off by 0.17 / -0.04 at the reference commit, so the gate guards
//!   against it getting worse rather than against an absolute limit);
//! - compression errors (mean and worst-5 %): `new <= baseline * 1.10 + 1e-4`;
//! - mean shear rate and ms/frame: informational only.
//!
//! Exit code 1 if any gated metric fails.

use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::time::Instant;

use mill_core::metrics::to_vertical_degrees;
use mill_core::params::{
    CoarseGrainingMode, Direction, LiftersParams, MediaParams, MillParams, SlurryParams, SpeedMode,
};
use mill_core::{Params, Simulation};

const FPS: f32 = 60.0;
/// Heavier metrics (toe/shoulder, stress, compression) are sampled every this many frames.
const SAMPLE_EVERY_FRAMES: u32 = 6;

struct Config {
    name: &'static str,
    max_balls: u32,
    resolution: u32,
    /// `true`: the shipped default mill (63 mm drum, 2 mm balls, coarse-graining Off), which is what
    /// the web app runs; `false`: the legacy 1 m drum with coarse-grained media.
    shipped: bool,
    /// Slurry viscosity (Pa s); `None` keeps the default (water, 0.001).
    viscosity_pa_s: Option<f32>,
}

/// Same Params as `perf_probe`'s configs of the same name.
const CONFIGS: [Config; 8] = [
    Config {
        name: "ui_fast",
        max_balls: 300,
        resolution: 25,
        shipped: false,
        viscosity_pa_s: None,
    },
    Config {
        name: "ui_accuracy",
        max_balls: 600,
        resolution: 50,
        shipped: false,
        viscosity_pa_s: None,
    },
    Config {
        name: "ship_fast",
        max_balls: 300,
        resolution: 25,
        shipped: true,
        viscosity_pa_s: None,
    },
    Config {
        name: "ship_acc",
        max_balls: 600,
        resolution: 50,
        shipped: true,
        viscosity_pa_s: None,
    },
    // Viscosity sweep on the shipped mill (Fast tier): the approximate speedups must hold for a
    // thin (0.01), a syrup-like (1) and a very viscous (50 Pa s) slurry, not only for water.
    Config {
        name: "ship_v0.01",
        max_balls: 300,
        resolution: 25,
        shipped: true,
        viscosity_pa_s: Some(0.01),
    },
    Config {
        name: "ship_v1",
        max_balls: 300,
        resolution: 25,
        shipped: true,
        viscosity_pa_s: Some(1.0),
    },
    Config {
        name: "ship_v50",
        max_balls: 300,
        resolution: 25,
        shipped: true,
        viscosity_pa_s: Some(50.0),
    },
    Config {
        name: "ship_acc_v50",
        max_balls: 600,
        resolution: 50,
        shipped: true,
        viscosity_pa_s: Some(50.0),
    },
];

struct Args {
    seeds: u64,
    settle_revs: f32,
    measure_revs: f32,
    baseline: Option<String>,
    compare: Option<String>,
    /// Comma-separated config names to run (default: all).
    only: Option<Vec<String>>,
}

impl Args {
    fn parse() -> Self {
        let mut a = Args {
            seeds: 3,
            settle_revs: 2.0,
            measure_revs: 3.0,
            baseline: None,
            compare: None,
            only: None,
        };
        let mut args = env::args().skip(1);
        while let Some(flag) = args.next() {
            let mut val = || args.next().expect("flag needs a value");
            match flag.as_str() {
                "--seeds" => a.seeds = val().parse().expect("integer"),
                "--settle-revs" => a.settle_revs = val().parse().expect("number"),
                "--measure-revs" => a.measure_revs = val().parse().expect("number"),
                "--baseline" => a.baseline = Some(val()),
                "--compare" => a.compare = Some(val()),
                "--only" => a.only = Some(val().split(',').map(str::to_string).collect()),
                "--help" | "-h" => {
                    println!(
                        "Flags: --seeds <u64, default 3> --settle-revs <f32, default 2> \
                         --measure-revs <f32, default 3> --baseline <csv to write> \
                         --compare <csv to compare against> --only <comma-separated config names>"
                    );
                    std::process::exit(0);
                }
                other => panic!("unknown flag: {other}"),
            }
        }
        a
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Gate {
    /// Relative 2 % floor.
    Rel,
    /// Absolute 2 degree floor, circular delta.
    Angle,
    /// `|new| <= |baseline| * 1.10 + 0.005`.
    EnergyFraction,
    /// `new <= baseline * 1.10 + 1e-4`.
    Compression,
    Info,
}

/// Metric names in report order, with their gate.
const METRICS: [(&str, Gate); 11] = [
    ("power_w", Gate::Rel),
    ("torque_nm", Gate::Rel),
    ("specific_power_w_per_kg", Gate::Rel),
    ("impacts_per_ball_per_s", Gate::Rel),
    ("toe_angle_deg", Gate::Angle),
    ("shoulder_angle_deg", Gate::Angle),
    ("unattributed_over_shaft", Gate::EnergyFraction),
    ("compression_mean", Gate::Compression),
    ("compression_worst5pct", Gate::Compression),
    ("mean_shear_rate_per_s", Gate::Info),
    ("ms_per_frame", Gate::Info),
];

fn gate_of(name: &str) -> Gate {
    METRICS
        .iter()
        .find(|(n, _)| *n == name)
        .map_or(Gate::Info, |(_, g)| *g)
}

/// Wrap-safe circular mean, in degrees `[0, 360)`; NaN if empty.
fn circular_mean_deg(x: &[f64]) -> f64 {
    if x.is_empty() {
        return f64::NAN;
    }
    let (s, c) = x.iter().fold((0.0, 0.0), |(s, c), d| {
        (s + d.to_radians().sin(), c + d.to_radians().cos())
    });
    s.atan2(c).to_degrees().rem_euclid(360.0)
}

fn mean_or_nan(x: &[f64]) -> f64 {
    if x.is_empty() {
        f64::NAN
    } else {
        x.iter().sum::<f64>() / x.len() as f64
    }
}

fn params_for(cfg: &Config, seed: u64) -> Params {
    let mut params = Params {
        mill: if cfg.shipped {
            MillParams::default()
        } else {
            MillParams {
                diameter_m: 1.0,
                speed_mode: SpeedMode::PercentCritical,
                speed_value: 70.0,
                direction: Direction::CounterClockwise,
            }
        },
        media: MediaParams::default(),
        lifters: LiftersParams {
            count: 0,
            ..LiftersParams::default()
        },
        slurry: SlurryParams {
            enabled: true,
            viscosity_pa_s: cfg
                .viscosity_pa_s
                .unwrap_or(SlurryParams::default().viscosity_pa_s),
            ..SlurryParams::default()
        },
        ..Params::default()
    };
    if !cfg.shipped {
        params.simulation.coarse_graining_mode = CoarseGrainingMode::Auto;
    }
    params.simulation.max_balls = cfg.max_balls;
    params.simulation.resolution = cfg.resolution;
    params.simulation.seed = seed;
    params.validate().expect("invalid params");
    params
}

/// One (config, seed) run: window means, keyed by metric name.
fn run_one(
    cfg: &Config,
    seed: u64,
    settle_revs: f32,
    measure_revs: f32,
) -> Vec<(&'static str, f64)> {
    let params = params_for(cfg, seed);
    let omega = params.mill.omega();
    let rev_s = 60.0 / params.mill.rpm();
    let settle_frames = (settle_revs * rev_s * FPS).round() as u32;
    let measure_frames = (measure_revs * rev_s * FPS).round().max(1.0) as u32;
    let mut sim = Simulation::new(params).expect("failed to build simulation");

    for _ in 0..settle_frames {
        sim.step(1.0 / FPS);
    }
    sim.reset_energy_budget();
    let spec_e0 = sim.metrics().stress.cumulative_specific_energy_j_per_kg;

    let mut step_ns = 0u128;
    let mut impacts = Vec::new();
    let mut toe = Vec::new();
    let mut shoulder = Vec::new();
    let mut comp_mean = Vec::new();
    let mut comp_worst = Vec::new();
    let mut shear_sum = 0.0f64;
    for f in 0..measure_frames {
        let t = Instant::now();
        sim.step(1.0 / FPS);
        step_ns += t.elapsed().as_nanos();
        shear_sum += sim.mean_shear_rate_per_s() as f64;
        if f % SAMPLE_EVERY_FRAMES == 0 {
            let m = sim.metrics();
            impacts.push(m.stress.impacts_total_per_ball_per_s as f64);
            if let Some(a) = m.toe_angle_rad {
                toe.push(to_vertical_degrees(a, omega) as f64);
            }
            if let Some(a) = m.shoulder_angle_rad {
                shoulder.push(to_vertical_degrees(a, omega) as f64);
            }
            let fluid = sim.fluid();
            if !fluid.is_empty() && fluid.rest_density > 0.0 {
                let rho0 = fluid.rest_density;
                let mut c: Vec<f64> = fluid
                    .densities()
                    .iter()
                    .map(|r| (r / rho0 - 1.0).max(0.0) as f64)
                    .collect();
                comp_mean.push(c.iter().sum::<f64>() / c.len() as f64);
                c.sort_by(|a, b| b.partial_cmp(a).unwrap());
                let top = &c[..(c.len() / 20).max(1)];
                comp_worst.push(top.iter().sum::<f64>() / top.len() as f64);
            }
        }
    }

    let eb = sim.energy_budget();
    let el = eb.elapsed_s.max(1e-12);
    let power = eb.shaft_work_j / el;
    let spec_e1 = sim.metrics().stress.cumulative_specific_energy_j_per_kg;
    vec![
        ("power_w", power),
        ("torque_nm", power / (omega as f64).abs().max(1e-12)),
        ("specific_power_w_per_kg", (spec_e1 - spec_e0) / el),
        ("impacts_per_ball_per_s", mean_or_nan(&impacts)),
        ("toe_angle_deg", circular_mean_deg(&toe)),
        ("shoulder_angle_deg", circular_mean_deg(&shoulder)),
        (
            "unattributed_over_shaft",
            eb.unattributed_j / eb.shaft_work_j,
        ),
        ("compression_mean", mean_or_nan(&comp_mean)),
        ("compression_worst5pct", mean_or_nan(&comp_worst)),
        ("mean_shear_rate_per_s", shear_sum / measure_frames as f64),
        ("ms_per_frame", step_ns as f64 / 1e6 / measure_frames as f64),
    ]
}

/// Mean and standard error (n - 1) of the finite values; NaN if none.
fn mean_se(v: &[f64]) -> (f64, f64) {
    let v: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let se = if v.len() > 1 {
        ((v.iter().map(|a| (a - m).powi(2)).sum::<f64>() / (n - 1.0)) / n).sqrt()
    } else {
        0.0
    };
    (m, se)
}

type Table = BTreeMap<(String, String), (f64, f64)>;

fn run_all(args: &Args) -> Table {
    let mut table = Table::new();
    for cfg in CONFIGS.iter().filter(|c| {
        args.only
            .as_ref()
            .is_none_or(|o| o.iter().any(|n| n == c.name))
    }) {
        let mut per_metric: BTreeMap<&'static str, Vec<f64>> = BTreeMap::new();
        for s in 0..args.seeds {
            let t = Instant::now();
            let r = run_one(cfg, 1 + s, args.settle_revs, args.measure_revs);
            eprintln!(
                "[{}] seed {} done in {:.0}s",
                cfg.name,
                1 + s,
                t.elapsed().as_secs_f32()
            );
            for (name, v) in r {
                per_metric.entry(name).or_default().push(v);
            }
        }
        for (name, vals) in per_metric {
            table.insert((cfg.name.to_string(), name.to_string()), mean_se(&vals));
        }
    }
    table
}

fn write_csv(path: &str, table: &Table) {
    let mut out = String::from("config,metric,mean,se\n");
    for cfg in &CONFIGS {
        for (name, _) in &METRICS {
            if let Some((m, se)) = table.get(&(cfg.name.to_string(), name.to_string())) {
                out.push_str(&format!("{},{},{:e},{:e}\n", cfg.name, name, m, se));
            }
        }
    }
    fs::write(path, out).expect("cannot write baseline CSV");
}

fn read_csv(path: &str) -> Table {
    let text = fs::read_to_string(path).expect("cannot read baseline CSV");
    let mut table = Table::new();
    for line in text.lines().skip(1).filter(|l| !l.trim().is_empty()) {
        let p: Vec<&str> = line.split(',').collect();
        assert!(p.len() == 4, "bad CSV line: {line}");
        table.insert(
            (p[0].to_string(), p[1].to_string()),
            (
                p[2].parse().expect("mean"),
                p[3].parse().expect("standard error"),
            ),
        );
    }
    table
}

fn print_table(table: &Table) {
    println!(
        "{:<12} {:<26} {:>14} {:>12}",
        "config", "metric", "mean", "se"
    );
    for cfg in &CONFIGS {
        for (name, _) in &METRICS {
            if let Some((m, se)) = table.get(&(cfg.name.to_string(), name.to_string())) {
                println!("{:<12} {:<26} {:>14.6e} {:>12.3e}", cfg.name, name, m, se);
            }
        }
    }
}

/// Returns `(delta, pass)`; `pass` is `None` for informational metrics.
fn judge(gate: Gate, a: (f64, f64), b: (f64, f64)) -> (f64, Option<bool>) {
    let (am, ase) = a;
    let (bm, bse) = b;
    let delta = if gate == Gate::Angle && am.is_finite() && bm.is_finite() {
        let d = (bm - am).rem_euclid(360.0);
        if d > 180.0 {
            d - 360.0
        } else {
            d
        }
    } else {
        bm - am
    };
    if am.is_nan() || bm.is_nan() {
        let both = am.is_nan() && bm.is_nan();
        return (
            if both { 0.0 } else { f64::NAN },
            if gate == Gate::Info { None } else { Some(both) },
        );
    }
    let stat = 2.0 * (ase * ase + bse * bse).sqrt();
    let pass = match gate {
        Gate::Rel => Some(delta.abs() <= stat.max(0.02 * am.abs())),
        Gate::Angle => Some(delta.abs() <= stat.max(2.0)),
        Gate::EnergyFraction => Some(bm.abs() <= am.abs() * 1.10 + 0.005),
        Gate::Compression => Some(bm <= am * 1.10 + 1e-4),
        Gate::Info => None,
    };
    (delta, pass)
}

fn main() {
    let args = Args::parse();
    let wall = Instant::now();
    println!(
        "tolerance_probe: seeds={} settle_revs={} measure_revs={}",
        args.seeds, args.settle_revs, args.measure_revs
    );
    let table = run_all(&args);
    println!();
    print_table(&table);
    if let Some(path) = &args.baseline {
        write_csv(path, &table);
        println!("\nbaseline written to {path}");
    }
    let mut failed = 0;
    if let Some(path) = &args.compare {
        let base = read_csv(path);
        println!("\nCOMPARE against {path}");
        println!(
            "{:<12} {:<26} {:>13} {:>13} {:>10}  verdict",
            "config", "metric", "baseline", "new", "delta %"
        );
        for cfg in &CONFIGS {
            if !table.keys().any(|k| k.0 == cfg.name) {
                continue; // not run this time (--only)
            }
            for (name, _) in &METRICS {
                let key = (cfg.name.to_string(), name.to_string());
                let (Some(&a), Some(&b)) = (base.get(&key), table.get(&key)) else {
                    println!(
                        "{:<12} {:<26} missing in baseline or new run",
                        cfg.name, name
                    );
                    continue;
                };
                let gate = gate_of(name);
                let (delta, pass) = judge(gate, a, b);
                let pct = if a.0.abs() > 1e-300 {
                    delta / a.0.abs() * 100.0
                } else {
                    f64::NAN
                };
                let verdict = match pass {
                    Some(true) => "PASS",
                    Some(false) => {
                        failed += 1;
                        "FAIL"
                    }
                    None => "info",
                };
                println!(
                    "{:<12} {:<26} {:>13.6e} {:>13.6e} {:>+10.3}  {}",
                    cfg.name, name, a.0, b.0, pct, verdict
                );
            }
        }
        println!("\n{failed} gated metric(s) FAILED");
    }
    println!("total wall time: {:.0}s", wall.elapsed().as_secs_f32());
    if failed > 0 {
        std::process::exit(1);
    }
}
