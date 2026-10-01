//! Coarse-graining invariance probe: which output metrics are truly independent of the
//! coarse-graining factor `k` (docs/PHYSICS.md ss3), and which drift with it.
//!
//! Not part of the test suite: a standalone diagnostic, e.g.
//!
//! ```text
//! cargo run -p mill-core --release --example coarse_graining_probe
//! cargo run -p mill-core --release --example coarse_graining_probe -- --drum-mm 500 --ball-mm 5 \
//!     --ks 1,1.5,2,3,4 --speeds 60,75 --fills 0.30 --slurry both --seeds 3
//! ```
//!
//! Fixes one physical case (drum diameter, true ball diameter, fill), then forces the coarse-
//! graining factor with `CoarseGrainingMode::Manual` (equivalent to the `Auto` mode's
//! `k = sqrt(N_true / max_balls)`, but hits exact `k` values). Every (condition, k, seed) run is
//! settled for `--settle-revs` drum revolutions and then sampled for `--measure-revs`; per-run
//! time averages are aggregated over seeds and printed as mean +/- seed standard deviation, with
//! the relative deviation from the first (lowest) `k` and a `*` when that deviation exceeds twice
//! the combined standard error of the two means (i.e. is unlikely to be seed noise alone).
//!
//! `collision_rate_per_s` and the impact-energy histogram are deliberately not covered: they are
//! known to scale as `1/k^2` and `k^2` respectively.
//!
//! Runs are executed in parallel on `--threads` worker threads (default: available cores - 1);
//! results are deterministic per (Params, seed) regardless of thread count.

use std::collections::BTreeMap;
use std::env;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use glam::Vec2;
use mill_core::metrics::to_vertical_degrees;
use mill_core::params::{CoarseGrainingMode, Direction, SpeedMode};
use mill_core::{Params, Simulation};

const FPS: f32 = 60.0;
/// Metrics are sampled every this many frames (10 Hz).
const SAMPLE_EVERY_FRAMES: u32 = 6;
/// The (more expensive) fluid/ball proximity statistics are sampled every this many frames (1 Hz).
const PROXIMITY_EVERY_FRAMES: u32 = 60;

struct Args {
    drum_mm: f32,
    ball_mm: f32,
    ks: Vec<f32>,
    speeds: Vec<f32>,
    fills: Vec<f32>,
    slurry: Vec<bool>,
    seeds: u32,
    settle_revs: f32,
    measure_revs: f32,
    resolution: u32,
    threads: usize,
    only_metrics: Option<Vec<String>>,
    /// Explicit `(k, resolution)` pairs (`--combos 1:50,2:25`); overrides `--ks`/`--resolution`.
    combos: Vec<(f32, u32)>,
    csv: bool,
    /// Shipped browser defaults (63 mm drum, 2 mm ball, 0.30 fill, 50/60/70 %Nc, slurry on/off,
    /// presets Realtime/Balanced/Accuracy = res 15/25/40).
    browser_defaults: bool,
    /// Reference mode: add a k=1 fine-lattice reference run per condition and print an error map.
    reference: bool,
    ref_res: u32,
    ref_substeps: u32,
    /// Optional `simulation.substeps` override for the non-reference runs.
    substeps: Option<u32>,
    /// Optional `simulation.pbf_iterations` override (the adaptive floor; validated in 1..=20).
    pbf_iters: Option<u32>,
    ks_given: bool,
}

fn parse_list(s: &str) -> Vec<f32> {
    s.split(',')
        .map(|v| {
            v.trim()
                .parse()
                .expect("expected a comma-separated list of numbers")
        })
        .collect()
}

impl Args {
    fn parse() -> Self {
        let mut a = Args {
            drum_mm: 500.0,
            ball_mm: 5.0,
            ks: vec![1.0, 1.5, 2.0, 3.0, 4.0],
            speeds: vec![60.0, 75.0],
            fills: vec![0.30],
            slurry: vec![false, true],
            seeds: 3,
            settle_revs: 6.0,
            measure_revs: 8.0,
            resolution: 50,
            threads: thread::available_parallelism()
                .map(|n| n.get().saturating_sub(1).max(1))
                .unwrap_or(1),
            only_metrics: None,
            combos: Vec::new(),
            csv: false,
            browser_defaults: false,
            reference: false,
            ref_res: 100,
            ref_substeps: 16,
            substeps: None,
            pbf_iters: None,
            ks_given: false,
        };
        if env::args().any(|f| f == "--browser-defaults") {
            let d = Params::default();
            a.browser_defaults = true;
            a.drum_mm = d.mill.diameter_m * 1000.0;
            a.ball_mm = d.media.ball_diameter_m * 1000.0;
            a.fills = vec![d.media.fill_fraction];
            a.speeds = vec![50.0, 60.0, 70.0];
            a.slurry = vec![false, true];
        }
        let mut args = env::args().skip(1);
        while let Some(flag) = args.next() {
            let mut val = || args.next().expect("flag needs a value");
            match flag.as_str() {
                "--drum-mm" => a.drum_mm = val().parse().expect("number"),
                "--ball-mm" => a.ball_mm = val().parse().expect("number"),
                "--ks" => {
                    a.ks = parse_list(&val());
                    a.ks_given = true;
                }
                "--browser-defaults" => {}
                "--reference" => a.reference = true,
                "--ref-res" => a.ref_res = val().parse().expect("number"),
                "--ref-substeps" => a.ref_substeps = val().parse().expect("number"),
                "--substeps" => a.substeps = Some(val().parse().expect("number")),
                "--pbf-iters" => a.pbf_iters = Some(val().parse().expect("number")),
                "--speeds" => a.speeds = parse_list(&val()),
                "--fills" => a.fills = parse_list(&val()),
                "--slurry" => {
                    a.slurry = match val().as_str() {
                        "on" => vec![true],
                        "off" => vec![false],
                        "both" => vec![false, true],
                        other => panic!("--slurry expects on|off|both, got {other}"),
                    }
                }
                "--seeds" => a.seeds = val().parse().expect("number"),
                "--settle-revs" => a.settle_revs = val().parse().expect("number"),
                "--measure-revs" => a.measure_revs = val().parse().expect("number"),
                "--resolution" => a.resolution = val().parse().expect("number"),
                "--threads" => a.threads = val().parse().expect("number"),
                "--metrics" => {
                    a.only_metrics = Some(val().split(',').map(|s| s.trim().to_string()).collect())
                }
                "--combos" => {
                    a.combos = val()
                        .split(',')
                        .map(|c| {
                            let (k, r) = c.split_once(':').expect("--combos expects k:res,k:res");
                            (k.trim().parse().expect("k"), r.trim().parse().expect("res"))
                        })
                        .collect()
                }
                "--csv" => a.csv = true,
                "--help" | "-h" => {
                    println!(
                        "Flags: --drum-mm <f32> --ball-mm <f32> --ks <list> --speeds <list %Nc> \
                         --fills <list> --slurry <on|off|both> --seeds <u32, >=3 recommended> \
                         --settle-revs <f32> --measure-revs <f32> --resolution <u32> \
                         --threads <usize> --metrics <comma list of names to print>                          --browser-defaults --reference --ref-res <u32, default 100>                          --ref-substeps <u32, default 16> --substeps <u32> --pbf-iters <u32, 1..=20> --combos k:res,..."
                    );
                    std::process::exit(0);
                }
                other => panic!("unknown flag: {other}"),
            }
        }
        a
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Condition {
    slurry: bool,
    /// Percent of critical speed, times 10 (integer key for ordering).
    speed_x10: u32,
    /// Fill fraction, times 1000.
    fill_x1000: u32,
}

impl Condition {
    fn label(&self) -> String {
        format!(
            "{} speed={:.0}%Nc fill={:.2}",
            if self.slurry { "WET" } else { "DRY" },
            self.speed_x10 as f32 / 10.0,
            self.fill_x1000 as f32 / 1000.0
        )
    }
}

#[derive(Clone, Copy)]
struct Job {
    cond: Condition,
    k: f32,
    res: u32,
    seed: u64,
    /// Reference run (k=1, fine lattice): raised sub-steps and a 2x longer measure window.
    is_ref: bool,
}

/// Running mean/count of an `Option`-valued or plain sample stream.
#[derive(Default)]
struct Acc {
    sum: f64,
    n: u32,
}

impl Acc {
    fn push(&mut self, v: f64) {
        if v.is_finite() {
            self.sum += v;
            self.n += 1;
        }
    }
    fn mean(&self) -> f64 {
        if self.n == 0 {
            f64::NAN
        } else {
            self.sum / self.n as f64
        }
    }
}

/// Mean and standard deviation of a raw sample series.
fn mean_std_of(x: &[f64]) -> (f64, f64) {
    let x: Vec<f64> = x.iter().copied().filter(|v| v.is_finite()).collect();
    if x.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let m = x.iter().sum::<f64>() / x.len() as f64;
    let v = x.iter().map(|a| (a - m) * (a - m)).sum::<f64>() / x.len() as f64;
    (m, v.sqrt())
}

/// Wrap-safe circular mean, in degrees, of angles already in `[0, 360)`.
fn circular_mean_deg(x: &[f64]) -> f64 {
    let x: Vec<f64> = x.iter().copied().filter(|v| v.is_finite()).collect();
    if x.is_empty() {
        return f64::NAN;
    }
    let (s, c) = x.iter().fold((0.0, 0.0), |(s, c), d| {
        (s + d.to_radians().sin(), c + d.to_radians().cos())
    });
    s.atan2(c).to_degrees().rem_euclid(360.0)
}

fn percentile(sorted: &[f32], q: f32) -> f32 {
    if sorted.is_empty() {
        return f32::NAN;
    }
    let idx = ((sorted.len() - 1) as f32 * q).round() as usize;
    sorted[idx]
}

/// Fluid/ball proximity statistics: fraction of balls with at least one fluid particle within
/// `r_ball + dx` of their centre ("wetted"), and the fluid-particle count in the annulus
/// `[r_ball, r_ball + 2 dx]` around each ball relative to the count a full hexagonal lattice would
/// have there ("ring occupancy", 1.0 = ball fully surrounded by slurry).
fn ball_slurry_proximity(
    balls: &mill_core::dem::Balls,
    fluid: &mill_core::pbf::FluidParticles,
) -> (f64, f64) {
    if balls.is_empty() || fluid.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let dx = fluid.h * 0.5;
    let r = balls.radius;
    let cell = (r + 2.0 * dx).max(dx);
    let key = |p: Vec2| ((p.x / cell).floor() as i32, (p.y / cell).floor() as i32);
    let mut grid: BTreeMap<(i32, i32), Vec<Vec2>> = BTreeMap::new();
    for &p in &fluid.x {
        grid.entry(key(p)).or_default().push(p);
    }
    let lattice_density = 2.0 / (3.0f32.sqrt() * dx * dx);
    let ring_area = std::f32::consts::PI * ((r + 2.0 * dx).powi(2) - r * r);
    let expected = lattice_density * ring_area;
    let mut wetted = 0u32;
    let mut ring_sum = 0.0f64;
    for &b in &balls.x {
        let (cx, cy) = key(b);
        let mut ring = 0u32;
        let mut any = false;
        for gx in cx - 1..=cx + 1 {
            for gy in cy - 1..=cy + 1 {
                if let Some(cellv) = grid.get(&(gx, gy)) {
                    for &p in cellv {
                        let d = (p - b).length();
                        if d <= r + dx {
                            any = true;
                        }
                        if d >= r && d <= r + 2.0 * dx {
                            ring += 1;
                        }
                    }
                }
            }
        }
        if any {
            wetted += 1;
        }
        ring_sum += (ring as f32 / expected) as f64;
    }
    let n = balls.len() as f64;
    (wetted as f64 / n, ring_sum / n)
}

/// One run: returns `(metric name, time-averaged value)` pairs in a fixed order.
fn run_one(args: &Args, job: Job) -> Vec<(&'static str, f64)> {
    let mut p = Params::default();
    p.mill.diameter_m = args.drum_mm / 1000.0;
    p.mill.speed_mode = SpeedMode::PercentCritical;
    p.mill.speed_value = job.cond.speed_x10 as f32 / 10.0;
    p.mill.direction = Direction::CounterClockwise;
    p.lifters.count = 0;
    p.media.ball_diameter_m = args.ball_mm / 1000.0;
    p.media.fill_fraction = job.cond.fill_x1000 as f32 / 1000.0;
    p.slurry.enabled = job.cond.slurry;
    p.simulation.resolution = job.res;
    p.simulation.seed = job.seed;
    p.simulation.coarse_graining_mode = CoarseGrainingMode::Manual;
    p.simulation.coarse_graining_k = job.k;
    // Manual mode ignores max_balls for k, but keep validation happy with a generous value.
    p.simulation.max_balls = 50_000;
    if job.is_ref {
        p.simulation.substeps = args.ref_substeps;
    } else if let Some(n) = args.substeps {
        p.simulation.substeps = n;
    }
    if let Some(n) = args.pbf_iters {
        p.simulation.pbf_iterations = n;
    }
    p.validate().expect("invalid params");

    let mut sim = Simulation::new(p).expect("failed to build simulation");
    let omega = p.mill.omega();
    let rev_s = 60.0 / p.mill.rpm();
    let settle_frames = (args.settle_revs * rev_s * FPS).round() as u32;
    let measure_revs = if job.is_ref {
        2.0 * args.measure_revs
    } else {
        args.measure_revs
    };
    let measure_frames = (measure_revs * rev_s * FPS).round() as u32;

    for _ in 0..settle_frames {
        sim.step(1.0 / FPS);
    }

    sim.reset_energy_budget();
    let mut acc: BTreeMap<&'static str, Acc> = BTreeMap::new();
    let mut centroid_deg: Vec<f64> = Vec::new();
    let mut ke_series: Vec<f64> = Vec::new();
    let mut toe_deg: Vec<f64> = Vec::new();
    let mut shoulder_deg: Vec<f64> = Vec::new();
    let mut pool_min_deg: Vec<f64> = Vec::new();
    let mut pool_max_deg: Vec<f64> = Vec::new();
    let mut fs_angle_deg: Vec<f64> = Vec::new();
    let mut fluid_centroid_deg: Vec<f64> = Vec::new();
    let mut overlap_peak = 0.0f64;
    let mut push = |name: &'static str, v: f64| acc.entry(name).or_default().push(v);

    for f in 0..measure_frames {
        sim.step(1.0 / FPS);
        if f % SAMPLE_EVERY_FRAMES == 0 {
            let m = sim.metrics();
            push("power_draw_w", m.power_draw_w as f64);
            push("torque_nm", m.torque_nm as f64);
            push("dissipated_power_w", m.dissipated_power_w as f64);
            push("total_kinetic_energy_j", m.total_kinetic_energy_j as f64);
            ke_series.push(m.total_kinetic_energy_j as f64);
            push("mean_shear_rate_per_s", m.mean_shear_rate_per_s as f64);
            push("max_ball_overlap_frac", m.max_ball_overlap_fraction as f64);
            overlap_peak = overlap_peak.max(m.max_ball_overlap_fraction as f64);
            push(
                "max_ball_wall_overlap_frac",
                m.max_ball_wall_overlap_fraction as f64,
            );
            push("coupling_clamp_hits", m.coupling_clamp_hits as f64);
            push("viscosity_iters", m.viscosity_solver_iterations as f64);
            push(
                "max_substep_disp_over_d",
                m.max_substep_displacement_over_diameter as f64,
            );
            if let Some(v) = m.mixing_index {
                push("mixing_index", v as f64);
            }
            if let Some(v) = m.pool_depth_m {
                push("pool_depth_m", v as f64);
            }
            if let Some(v) = m.free_surface_offset_m {
                push("free_surface_offset_m", v as f64);
            }
            if let Some(v) = m.max_fluid_compression_error_fraction {
                push("fluid_max_compression_err", v as f64);
            }
            if let Some(v) = m.mean_fluid_compression_error_fraction {
                push("fluid_mean_compression_err", v as f64);
            }
            if let Some(a) = m.toe_angle_rad {
                toe_deg.push(to_vertical_degrees(a, omega) as f64);
            }
            if let Some(a) = m.shoulder_angle_rad {
                shoulder_deg.push(to_vertical_degrees(a, omega) as f64);
            }
            if let (Some(lo), Some(hi)) = (m.pool_angle_min_rad, m.pool_angle_max_rad) {
                pool_min_deg.push(to_vertical_degrees(lo, omega) as f64);
                pool_max_deg.push(to_vertical_degrees(hi, omega) as f64);
            }
            if let Some(a) = m.free_surface_angle_rad {
                // Line direction, mod 180 deg; fold to a tilt from horizontal in (-90, 90].
                let d = (a.to_degrees()).rem_euclid(180.0);
                fs_angle_deg.push(if d > 90.0 { d as f64 - 180.0 } else { d as f64 });
            }
            if let Some((cx, cy)) = m.charge_centroid_m {
                let rr = (cx * cx + cy * cy).sqrt();
                push("charge_centroid_radius_m", rr as f64);
                let deg = to_vertical_degrees(cy.atan2(cx), omega) as f64;
                centroid_deg.push(deg);
            }

            // Ball speed distribution (linear speed, and spin normalised by wall speed).
            let balls = sim.balls();
            if !balls.is_empty() {
                let mut speeds: Vec<f32> = balls.v.iter().map(|v| v.length()).collect();
                let n = speeds.len() as f64;
                let mean = speeds.iter().map(|&s| s as f64).sum::<f64>() / n;
                let rms = (speeds.iter().map(|&s| (s as f64) * (s as f64)).sum::<f64>() / n).sqrt();
                speeds.sort_by(|a, b| a.partial_cmp(b).unwrap());
                push("ball_speed_mean_m_s", mean);
                push("ball_speed_rms_m_s", rms);
                push("ball_speed_p50_m_s", percentile(&speeds, 0.5) as f64);
                push("ball_speed_p90_m_s", percentile(&speeds, 0.9) as f64);
                push("ball_speed_p99_m_s", percentile(&speeds, 0.99) as f64);
                push(
                    "ball_speed_max_m_s",
                    *speeds.last().unwrap_or(&f32::NAN) as f64,
                );
                let spin = balls.omega.iter().map(|w| w.abs() as f64).sum::<f64>() / n;
                push("ball_spin_mean_rad_s", spin);
                // Spin in units of the geometric roll-without-slip rate |v| / r.
                let roll_ratio = balls
                    .v
                    .iter()
                    .zip(&balls.omega)
                    .map(|(v, w)| (w.abs() * balls.radius) as f64 / (v.length() as f64 + 1e-3))
                    .sum::<f64>()
                    / n;
                push("ball_spin_over_roll_rate", roll_ratio);
            }

            let fluid = sim.fluid();
            if !fluid.is_empty() {
                let n = fluid.len() as f64;
                let fmean = fluid.v.iter().map(|v| v.length() as f64).sum::<f64>() / n;
                push("fluid_speed_mean_m_s", fmean);
                let c = fluid.x.iter().fold(Vec2::ZERO, |a, &p| a + p) / fluid.len() as f32;
                push(
                    "fluid_centroid_radius_m",
                    (c.x * c.x + c.y * c.y).sqrt() as f64,
                );
                fluid_centroid_deg.push(to_vertical_degrees(c.y.atan2(c.x), omega) as f64);
            }
        }
        if f % PROXIMITY_EVERY_FRAMES == 0 {
            let (wetted, ring) = ball_slurry_proximity(sim.balls(), sim.fluid());
            push("balls_wetted_fraction", wetted);
            push("ball_ring_slurry_occupancy", ring);
        }
    }

    let mut out: Vec<(&'static str, f64)> = acc.iter().map(|(k, v)| (*k, v.mean())).collect();
    // Window-mean energy budget (unsmoothed), in W = J/s per metre depth.
    let eb = sim.energy_budget();
    let el = eb.elapsed_s.max(1e-12);
    let power_mean_w = eb.shaft_work_j / el;
    out.push(("power_mean_w", power_mean_w));
    out.push((
        "torque_mean_nm",
        power_mean_w / (omega as f64).abs().max(1e-12),
    ));
    out.push(("budget_shaft_work_w", power_mean_w));
    out.push((
        "budget_ball_contact_dissipation_w",
        eb.ball_contact_dissipation_j / el,
    ));
    out.push(("budget_fluid_wall_slip_w", eb.fluid_wall_slip_j / el));
    out.push(("budget_fluid_viscous_w", eb.fluid_viscous_j / el));
    out.push((
        "budget_fluid_clamp_removed_w",
        eb.fluid_clamp_removed_j / el,
    ));
    out.push(("budget_interface_created_w", eb.interface_created_j / el));
    out.push(("budget_unattributed_w", eb.unattributed_j / el));
    out.push((
        "unattributed_pct_of_shaft",
        100.0 * eb.unattributed_j / eb.shaft_work_j,
    ));
    out.push((
        "torque_from_impulse_nm",
        -eb.fluid_wall_angular_impulse / el,
    ));
    out.push(("max_ball_overlap_peak", overlap_peak));
    out.push(("toe_angle_deg", circular_mean_deg(&toe_deg)));
    out.push(("shoulder_angle_deg", circular_mean_deg(&shoulder_deg)));
    out.push(("dynamic_repose_span_deg", {
        let t = circular_mean_deg(&toe_deg);
        let s = circular_mean_deg(&shoulder_deg);
        if t.is_finite() && s.is_finite() {
            (s - t).rem_euclid(360.0).min((t - s).rem_euclid(360.0))
        } else {
            f64::NAN
        }
    }));
    out.push((
        "toe_defined_fraction",
        toe_deg.len() as f64 / ke_series.len().max(1) as f64,
    ));
    out.push((
        "charge_centroid_angle_deg",
        circular_mean_deg(&centroid_deg),
    ));
    out.push(("pool_angle_min_deg", circular_mean_deg(&pool_min_deg)));
    out.push(("pool_angle_max_deg", circular_mean_deg(&pool_max_deg)));
    out.push(("free_surface_tilt_deg", mean_std_of(&fs_angle_deg).0));
    out.push((
        "fluid_centroid_angle_deg",
        circular_mean_deg(&fluid_centroid_deg),
    ));
    // Sloshing / oscillation amplitude: std of the (circularly unwrapped about its mean) centroid
    // angle over the window, and relative fluctuation of total KE.
    let cmean = circular_mean_deg(&centroid_deg);
    let dev: Vec<f64> = centroid_deg
        .iter()
        .map(|d| {
            let raw = d - cmean;
            raw - 360.0 * (raw / 360.0).round()
        })
        .collect();
    out.push(("centroid_angle_std_deg", mean_std_of(&dev).1));
    let (ke_m, ke_s) = mean_std_of(&ke_series);
    out.push(("ke_relative_fluctuation", ke_s / ke_m.max(1e-12)));
    out
}

fn main() {
    let args = Arc::new(Args::parse());
    if args.seeds < 3 {
        eprintln!("warning: fewer than 3 seeds; seed spread is not meaningful");
    }

    let mut args = args;
    let wall_start = std::time::Instant::now();
    {
        let a = Arc::get_mut(&mut args).expect("sole owner");
        if a.combos.is_empty() {
            a.combos = if a.browser_defaults && !a.ks_given {
                // Shipped presets Realtime / Balanced / Accuracy (web/src/params/presets.ts).
                vec![(1.0, 15), (1.0, 25), (1.0, 40)]
            } else {
                a.ks.iter().map(|&k| (k, a.resolution)).collect()
            };
        }
        if a.reference {
            let ref_res = a.ref_res;
            a.combos.retain(|&(k, r)| !(k == 1.0 && r == ref_res));
            // The reference goes first: the generic table then reports deviations against it.
            a.combos.insert(0, (1.0, ref_res));
        }
    }
    let ref_combo = (1.0f32, args.ref_res);
    let mut jobs: Vec<Job> = Vec::new();
    for &slurry in &args.slurry {
        for &speed in &args.speeds {
            for &fill in &args.fills {
                for &(k, res) in &args.combos {
                    for s in 0..args.seeds {
                        jobs.push(Job {
                            cond: Condition {
                                slurry,
                                speed_x10: (speed * 10.0).round() as u32,
                                fill_x1000: (fill * 1000.0).round() as u32,
                            },
                            k,
                            res,
                            seed: 1 + s as u64,
                            is_ref: args.reference && (k, res) == ref_combo,
                        });
                    }
                }
            }
        }
    }

    // Most expensive first (finest fluid), so the worker pool's tail is short.
    jobs.sort_by_key(|j| std::cmp::Reverse(j.res));

    // Header: report true ball count and the ball count at each k.
    {
        let mut p = Params::default();
        p.mill.diameter_m = args.drum_mm / 1000.0;
        p.media.ball_diameter_m = args.ball_mm / 1000.0;
        for &fill in &args.fills {
            p.media.fill_fraction = fill;
            println!(
                "case: drum {:.0} mm, true ball {:.2} mm, fill {:.2}, N_true = {:.0}, fluid res base {}; \
                 balls per k: {}",
                args.drum_mm,
                args.ball_mm,
                fill,
                p.true_ball_count(),
                args.resolution,
                args.combos
                    .iter()
                    .map(|&(k, _)| format!("k={k}:{}", (p.true_ball_count() / (k * k)).round()))
                    .collect::<Vec<_>>()
                    .join(" ")
            );
        }
        if args.reference {
            println!(
                "REFERENCE mode: k=1 res={} substeps={} measure window {} revs (2x); \
                 non-reference substeps: {}",
                args.ref_res,
                args.ref_substeps,
                2.0 * args.measure_revs,
                args.substeps
                    .map_or("Params default".to_string(), |n| n.to_string())
            );
        }
        println!(
            "settle {} revs, measure {} revs, {} seeds, {} jobs, {} threads",
            args.settle_revs,
            args.measure_revs,
            args.seeds,
            jobs.len(),
            args.threads
        );
    }

    let next = Arc::new(AtomicUsize::new(0));
    let jobs = Arc::new(jobs);
    type Results = Vec<(Job, Vec<(&'static str, f64)>)>;
    let results: Arc<Mutex<Results>> = Arc::new(Mutex::new(Vec::new()));
    let mut handles = Vec::new();
    for _ in 0..args.threads {
        let (next, jobs, results, args) =
            (next.clone(), jobs.clone(), results.clone(), args.clone());
        handles.push(thread::spawn(move || loop {
            let i = next.fetch_add(1, Ordering::SeqCst);
            if i >= jobs.len() {
                break;
            }
            let t0 = std::time::Instant::now();
            let r = run_one(&args, jobs[i]);
            eprintln!(
                "[{}/{}] {} k={} res={} seed={} done in {:.0}s",
                i + 1,
                jobs.len(),
                jobs[i].cond.label(),
                jobs[i].k,
                jobs[i].res,
                jobs[i].seed,
                t0.elapsed().as_secs_f32()
            );
            let mut r = r;
            r.push(("wall_s", t0.elapsed().as_secs_f64()));
            results.lock().unwrap().push((jobs[i], r));
        }));
    }
    for h in handles {
        h.join().unwrap();
    }
    let results = results.lock().unwrap();

    // Group: condition -> metric -> k(bits) -> per-seed values.
    let mut conds: Vec<Condition> = results.iter().map(|(j, _)| j.cond).collect();
    conds.sort();
    conds.dedup();
    for cond in conds {
        println!("\n=== {} ===", cond.label());
        let mut table: BTreeMap<&'static str, BTreeMap<(u32, u32), Vec<f64>>> = BTreeMap::new();
        for (j, r) in results.iter().filter(|(j, _)| j.cond == cond) {
            for (name, v) in r {
                table
                    .entry(name)
                    .or_default()
                    .entry((j.k.to_bits(), j.res))
                    .or_default()
                    .push(*v);
            }
        }
        let ks: Vec<(f32, u32)> = args.combos.clone();
        print!("{:<30}", "metric (mean +/- seed sd)");
        for (k, r) in &ks {
            print!("{:>26}", format!("k={k},res={r}"));
        }
        println!();
        for (name, per_k) in &table {
            if let Some(only) = &args.only_metrics {
                if !only.iter().any(|o| o == name) {
                    continue;
                }
            }
            let stat = |k: (f32, u32)| -> Option<(f64, f64, f64)> {
                let v: Vec<f64> = per_k
                    .get(&(k.0.to_bits(), k.1))?
                    .iter()
                    .copied()
                    .filter(|x| x.is_finite())
                    .collect();
                if v.is_empty() {
                    return None;
                }
                let n = v.len() as f64;
                let m = v.iter().sum::<f64>() / n;
                let sd = if v.len() > 1 {
                    (v.iter().map(|a| (a - m) * (a - m)).sum::<f64>() / (n - 1.0)).sqrt()
                } else {
                    f64::NAN
                };
                Some((m, sd, sd / n.sqrt()))
            };
            let reference = stat(ks[0]);
            print!("{name:<30}");
            for (i, &k) in ks.iter().enumerate() {
                if args.csv {
                    if let Some((m, sd, _)) = stat(k) {
                        eprintln!("CSV,{name},{},{},{m:e},{sd:e}", k.0, k.1);
                    }
                }
                match (stat(k), reference) {
                    (Some((m, sd, se)), Some((rm, _, rse))) => {
                        if i == 0 {
                            print!("{:>26}", format!("{m:.4e}+/-{sd:.1e}"));
                        } else {
                            let dev = if rm.abs() > 1e-12 {
                                (m - rm) / rm.abs() * 100.0
                            } else {
                                f64::NAN
                            };
                            let sig = (m - rm).abs() > 2.0 * (se * se + rse * rse).sqrt();
                            print!(
                                "{:>26}",
                                format!(
                                    "{m:.3e}+/-{sd:.0e} {dev:+.0}%{}",
                                    if sig { "*" } else { "" }
                                )
                            );
                        }
                    }
                    _ => print!("{:>26}", "n/a"),
                }
            }
            println!();
        }
    }
    println!(
        "\n('*' = deviation from the first k exceeds 2x combined standard error of the means)"
    );
    if args.reference {
        error_map(&args, &results, ref_combo);
    }
    println!(
        "\ntotal wall time: {:.0}s",
        wall_start.elapsed().as_secs_f32()
    );
}

/// Solver diagnostics rather than physical outputs: excluded from the pass / <=10% / fail tally
/// (still printed in the generic table).
const DIAGNOSTIC_METRICS: &[&str] = &[
    "wall_s",
    "coupling_clamp_hits",
    "viscosity_iters",
    "max_substep_disp_over_d",
    "max_ball_overlap_frac",
    "max_ball_overlap_peak",
    "max_ball_wall_overlap_frac",
    "fluid_max_compression_err",
    "fluid_mean_compression_err",
    "toe_defined_fraction",
];

/// Mean and standard error of the mean of the finite values.
fn mean_se(v: &[f64]) -> Option<(f64, f64)> {
    let v: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return None;
    }
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let se = if v.len() > 1 {
        ((v.iter().map(|a| (a - m) * (a - m)).sum::<f64>() / (n - 1.0)) / n).sqrt()
    } else {
        f64::NAN
    };
    Some((m, se))
}

type RunResults = [(Job, Vec<(&'static str, f64)>)];

/// Acceptance error map: each non-reference (k, res) against the k=1 fine-lattice reference of the
/// same condition. Power and torque PASS only if |dev| <= 3% and the deviation is not significant
/// at 2x combined SE; other physical metrics are classed pass / <=10% / fail.
fn error_map(args: &Args, results: &RunResults, ref_combo: (f32, u32)) {
    let values = |cond: Condition, combo: (f32, u32), name: &str| -> Vec<f64> {
        results
            .iter()
            .filter(|(j, _)| j.cond == cond && (j.k, j.res) == combo)
            .filter_map(|(_, r)| r.iter().find(|(n, _)| *n == name).map(|(_, v)| *v))
            .collect()
    };
    let mut conds: Vec<Condition> = results.iter().map(|(j, _)| j.cond).collect();
    conds.sort();
    conds.dedup();
    println!(
        "\n================ ERROR MAP vs reference (k=1, res={}) ================",
        ref_combo.1
    );
    println!("PASS = |dev| <= 3% AND deviation <= 2x combined SE.");
    let mut summary: Vec<String> = Vec::new();
    for cond in conds {
        println!("\n--- {} ---", cond.label());
        if let Some((m, _)) = mean_se(&values(cond, ref_combo, "max_ball_overlap_peak")) {
            println!(
                "reference peak ball overlap = {:.1}% of radius (target < 10%): {}",
                m * 100.0,
                if m < 0.10 { "OK" } else { "EXCEEDS" }
            );
        }
        let mut names: Vec<&str> = Vec::new();
        for (_, r) in results.iter().filter(|(j, _)| j.cond == cond) {
            for (n, _) in r {
                if !names.contains(n) {
                    names.push(n);
                }
            }
        }
        for &combo in args.combos.iter().filter(|&&c| c != ref_combo) {
            println!("[k={}, res={}]", combo.0, combo.1);
            let mut head = format!("[k={},res={}] {} ", combo.0, combo.1, cond.label());
            let (mut n_pass, mut n_le10, mut n_fail, mut n_na) = (0, 0, 0, 0);
            let mut fails: Vec<String> = Vec::new();
            for name in &names {
                if DIAGNOSTIC_METRICS.contains(name) {
                    continue;
                }
                let (Some((rm, rse)), Some((m, se))) = (
                    mean_se(&values(cond, ref_combo, name)),
                    mean_se(&values(cond, combo, name)),
                ) else {
                    n_na += 1;
                    continue;
                };
                if rm.abs() < 1e-12 {
                    n_na += 1;
                    continue;
                }
                let dev = (m - rm) / rm.abs() * 100.0;
                let sig = (m - rm).abs() > 2.0 * (se * se + rse * rse).sqrt();
                let pass = dev.abs() <= 3.0 && !sig;
                if *name == "power_draw_w" || *name == "torque_nm" {
                    let verdict = if pass { "PASS" } else { "FAIL" };
                    println!(
                        "  {name:<14} {m:.4e} +/- {se:.1e} (SE)  ref {rm:.4e} +/- {rse:.1e}  \
                         dev {dev:+.1}%  sig={}  {verdict}",
                        if sig { "yes" } else { "no" },
                    );
                    let short = if *name == "power_draw_w" {
                        "power"
                    } else {
                        "torque"
                    };
                    head.push_str(&format!("{short}={dev:+.1}% {verdict} "));
                } else if pass {
                    n_pass += 1;
                } else if dev.abs() <= 10.0 {
                    n_le10 += 1;
                } else {
                    n_fail += 1;
                    fails.push(format!("{name}{dev:+.0}%"));
                }
            }
            println!("  other metrics: {n_pass} pass, {n_le10} <=10%, {n_fail} fail, {n_na} n/a");
            if !fails.is_empty() {
                println!("  failing: {}", fails.join(" "));
            }
            head.push_str(&format!(
                "| other: {n_pass} pass / {n_le10} <=10% / {n_fail} fail"
            ));
            summary.push(head);
        }
    }
    println!("\n================ ACCEPTANCE SUMMARY ================");
    for l in summary {
        println!("{l}");
    }
}
