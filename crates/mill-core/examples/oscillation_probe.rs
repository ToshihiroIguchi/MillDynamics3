//! Reproduction probe for the periodic media-charge oscillation ("surging" or "slumping")
//! reported at low mill speed with a smooth (lifter-less) wall -- see docs/PHYSICS.md ss9's
//! "Periodic media-charge oscillation" entry, which this binary's measurements are cited from.
//! Not part of the test suite: a standalone diagnostic run, e.g.
//!
//! ```text
//! cargo run -p mill-core --release --example oscillation_probe -- --percent-critical 20 --seed 1
//! ```
//!
//! Runs a fixed geometry (default `Params`, no lifters) at the requested speed/slurry config, then
//! samples the charge's motion at 60 Hz after a settle period. It reports whether a periodic
//! signal is present and, if so, classifies it as "surging" (whole charge sticks to the wall then
//! releases as a block -- wall slip and centroid angle oscillate in phase) or "slumping" (the bed
//! surface periodically avalanches while the wall-adjacent layer keeps rolling/sliding -- wall
//! slip stays roughly flat while the shoulder/centroid angle saw-tooths).
//!
//! This binary makes no solver changes; it only measures the existing model. Besides the geometry/
//! speed/friction flags above, it also accepts `--rolling-friction`, `--fill` (ball fill fraction),
//! `--slurry-fill` (slurry fill fraction), and `--restitution-wall` (ball-wall restitution) for
//! sweeping those parameters too. Run with `--help` for every flag.

use std::env;
use std::fs::File;
use std::io::Write as _;

use mill_core::geometry::Drum;
use mill_core::metrics::{charge_centroid, charge_toe_shoulder, to_vertical_degrees};
use mill_core::params::{Direction, SpeedMode};
use mill_core::{Params, Simulation};

/// Distance from the drum wall, as a multiple of the ball radius, within which a ball is counted
/// as "wall-adjacent" for the slip-ratio measurement. Slightly looser than metrics.rs's own
/// `WALL_MARGIN_BALL_RADII` (2.5) since this probe wants the outermost layer specifically, not the
/// whole toe/shoulder cluster.
const WALL_ADJACENT_BALL_RADII: f32 = 1.5;

struct Args {
    percent_critical: f32,
    slurry_enabled: bool,
    viscosity_pa_s: f32,
    friction_ball_wall: f32,
    friction_ball_ball: f32,
    rolling_friction: f32,
    fill_fraction: f32,
    slurry_fill_fraction: f32,
    restitution_wall: f32,
    lifters_count: u32,
    max_balls: u32,
    resolution: u32,
    seed: u64,
    settle_s: f32,
    measure_s: f32,
    csv_path: Option<String>,
}

impl Args {
    fn parse() -> Self {
        let mut a = Args {
            percent_critical: 20.0,
            slurry_enabled: true,
            viscosity_pa_s: 50.0,
            friction_ball_wall: 0.35,
            friction_ball_ball: 0.25,
            rolling_friction: 0.01,
            fill_fraction: 0.30,
            slurry_fill_fraction: 0.35,
            restitution_wall: 0.5,
            lifters_count: 0,
            max_balls: 600,
            resolution: 40,
            seed: 1,
            settle_s: 6.0,
            measure_s: 8.0,
            csv_path: None,
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
                "--percent-critical" => a.percent_critical = next_f32(),
                "--slurry" => {
                    let v = args.next().expect("--slurry needs on/off");
                    a.slurry_enabled = v == "on";
                }
                "--viscosity" => a.viscosity_pa_s = next_f32(),
                "--friction-ball-wall" => a.friction_ball_wall = next_f32(),
                "--friction-ball-ball" => a.friction_ball_ball = next_f32(),
                "--rolling-friction" => a.rolling_friction = next_f32(),
                "--fill" => a.fill_fraction = next_f32(),
                "--slurry-fill" => a.slurry_fill_fraction = next_f32(),
                "--restitution-wall" => a.restitution_wall = next_f32(),
                "--lifters" => a.lifters_count = next_f32() as u32,
                "--max-balls" => a.max_balls = next_f32() as u32,
                "--resolution" => a.resolution = next_f32() as u32,
                "--seed" => a.seed = next_f32() as u64,
                "--settle-s" => a.settle_s = next_f32(),
                "--measure-s" => a.measure_s = next_f32(),
                "--csv" => a.csv_path = Some(args.next().expect("--csv needs a path")),
                "--help" | "-h" => {
                    println!(
                        "Flags: --percent-critical <f32> --slurry <on|off> --viscosity <f32> \
                         --friction-ball-wall <f32> --friction-ball-ball <f32> \
                         --rolling-friction <f32> --fill <f32> --slurry-fill <f32> \
                         --restitution-wall <f32> --lifters <u32> --max-balls <u32> \
                         --resolution <u32> --seed <u64> --settle-s <f32> --measure-s <f32> \
                         --csv <path>"
                    );
                    std::process::exit(0);
                }
                other => panic!("unknown flag: {other}"),
            }
        }
        a
    }
}

/// One 60 Hz sample of the charge's motion, taken during the measurement window.
struct Sample {
    t_s: f32,
    centroid_deg: f32,
    /// Charge centroid's distance from the drum axis (m) -- the effective pendulum arm length for
    /// the sloshing-mode sanity check in `report`.
    centroid_radius_m: f32,
    toe_deg: Option<f32>,
    shoulder_deg: Option<f32>,
    wall_slip_ratio: f32,
    total_ke_j: f32,
}

fn main() {
    let args = Args::parse();

    let mut params = Params::default();
    params.mill.speed_mode = SpeedMode::PercentCritical;
    params.mill.speed_value = args.percent_critical;
    params.mill.direction = Direction::CounterClockwise;
    params.lifters.count = args.lifters_count;
    params.slurry.enabled = args.slurry_enabled;
    params.slurry.viscosity_pa_s = args.viscosity_pa_s;
    params.media.friction_ball_wall = args.friction_ball_wall;
    params.media.friction_ball_ball = args.friction_ball_ball;
    params.media.rolling_friction = args.rolling_friction;
    params.media.fill_fraction = args.fill_fraction;
    params.media.restitution_ball_wall = args.restitution_wall;
    params.slurry.fill_fraction = args.slurry_fill_fraction;
    params.simulation.max_balls = args.max_balls;
    params.simulation.resolution = args.resolution;
    params.simulation.seed = args.seed;
    params.validate().expect("invalid params");

    let mut sim = Simulation::new(params).expect("failed to build simulation");
    let drum = Drum::new(params.mill.radius_m(), params.mill.omega(), params.lifters);

    let settle_frames = (args.settle_s * 60.0).round() as u32;
    for _ in 0..settle_frames {
        sim.step(1.0 / 60.0);
    }

    let measure_frames = (args.measure_s * 60.0).round() as u32;
    let mut samples = Vec::with_capacity(measure_frames as usize);
    for f in 0..measure_frames {
        sim.step(1.0 / 60.0);
        let balls = sim.balls();
        let drum_angle = sim.drum_angle();

        let centroid = charge_centroid(balls);
        let centroid_deg = centroid
            .map(|(x, y)| to_vertical_degrees(y.atan2(x), drum.omega))
            .unwrap_or(f32::NAN);
        let centroid_radius_m = centroid
            .map(|(x, y)| (x * x + y * y).sqrt())
            .unwrap_or(f32::NAN);
        let (toe, shoulder) = charge_toe_shoulder(balls, &drum, drum_angle);
        let toe_deg = toe.map(|a| to_vertical_degrees(a, drum.omega));
        let shoulder_deg = shoulder.map(|a| to_vertical_degrees(a, drum.omega));

        let wall_slip_ratio = mean_wall_slip_ratio(balls, &drum, drum_angle);
        let total_ke_j = mill_core::metrics::total_kinetic_energy_j(balls);

        samples.push(Sample {
            t_s: f as f32 / 60.0,
            centroid_deg,
            centroid_radius_m,
            toe_deg,
            shoulder_deg,
            wall_slip_ratio,
            total_ke_j,
        });
    }

    if let Some(path) = &args.csv_path {
        write_csv(path, &samples);
    }

    report(&args, &samples);
}

/// Mean, over balls within `WALL_ADJACENT_BALL_RADII * radius` of the wall, of the contact-point
/// tangential slip speed relative to the wall's own contact-point speed -- the same `v_t` formula
/// `dem.rs` step 5 uses for its wall-friction pass, so this measures exactly the quantity that
/// pass drives toward zero (0 = perfectly co-rotating with the wall, ~1 = sliding at roughly the
/// wall's own speed, since it is normalized by `|drum.omega| * radius_m`).
fn mean_wall_slip_ratio(balls: &mill_core::dem::Balls, drum: &Drum, drum_angle: f32) -> f32 {
    if balls.is_empty() || drum.omega.abs() < 1e-6 {
        return f32::NAN;
    }
    let r = balls.radius;
    let margin = WALL_ADJACENT_BALL_RADII * r;
    let wall_speed_scale = drum.omega.abs() * drum.radius_m;

    let mut sum = 0.0f32;
    let mut count = 0u32;
    for i in 0..balls.len() {
        let (dist, n_hat) = drum.sdf_world(balls.x[i], drum_angle);
        if dist >= margin {
            continue;
        }
        let t_hat = glam::Vec2::new(-n_hat.y, n_hat.x);
        let v_wall = drum.wall_velocity(balls.x[i] - r * n_hat);
        let v_t = (balls.v[i] - v_wall).dot(t_hat) - r * balls.omega[i];
        sum += v_t.abs();
        count += 1;
    }
    if count == 0 {
        return f32::NAN;
    }
    (sum / count as f32) / wall_speed_scale
}

fn write_csv(path: &str, samples: &[Sample]) {
    let mut f = File::create(path).expect("failed to create csv file");
    writeln!(
        f,
        "t_s,centroid_deg,centroid_radius_m,toe_deg,shoulder_deg,wall_slip_ratio,total_ke_j"
    )
    .unwrap();
    for s in samples {
        writeln!(
            f,
            "{},{},{},{},{},{},{}",
            s.t_s,
            s.centroid_deg,
            s.centroid_radius_m,
            s.toe_deg.map(|v| v.to_string()).unwrap_or_default(),
            s.shoulder_deg.map(|v| v.to_string()).unwrap_or_default(),
            s.wall_slip_ratio,
            s.total_ke_j,
        )
        .unwrap();
    }
}

/// Unwraps a sequence of angles (degrees, each in `[0, 360)`) into a continuous series by
/// accumulating the shortest signed step between consecutive samples -- needed because the raw
/// centroid angle can cross the `0`/`360` seam even during small, ordinary back-and-forth motion,
/// which would otherwise look like a huge spurious jump to the statistics below.
fn unwrap_degrees(angles: &[f32]) -> Vec<f32> {
    if angles.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(angles.len());
    out.push(angles[0]);
    for w in angles.windows(2) {
        let prev_unwrapped = *out.last().unwrap();
        let raw_delta = w[1] - w[0];
        // Shortest signed equivalent of `raw_delta` modulo 360.
        let wrapped_delta = raw_delta - 360.0 * (raw_delta / 360.0).round();
        out.push(prev_unwrapped + wrapped_delta);
    }
    out
}

/// First difference of `x` (length `x.len() - 1`), scaled to a per-second rate at the given
/// sample rate. A slow monotonic drift in the raw signal (this project's charge does not start at
/// its final dynamic angle -- it keeps rotating up to it during the settle period, and even
/// "settled" cascading has a slowly evolving toe/shoulder) shows up in the *raw* signal's
/// autocorrelation as a spurious near-1.0 correlation at every short lag, since nearby samples of
/// a smooth trend are highly correlated regardless of whether anything is actually periodic --
/// differencing removes that trend almost entirely (a linear drift in `x` becomes a constant in
/// the derivative), leaving genuine oscillation as the dominant remaining structure.
fn finite_difference_rate(x: &[f32], sample_rate_hz: f32) -> Vec<f32> {
    x.windows(2)
        .map(|w| (w[1] - w[0]) * sample_rate_hz)
        .collect()
}

/// Removes the best-fit linear trend from `x` (least squares against sample index), guarding
/// against whatever slow drift `finite_difference_rate` didn't fully remove (e.g. a drift whose
/// *rate* is itself still slowly changing, i.e. a curved, not linear, original trend).
fn linear_detrend(x: &[f32]) -> Vec<f32> {
    let n = x.len() as f32;
    if x.len() < 2 {
        return x.to_vec();
    }
    let t_mean = (n - 1.0) / 2.0;
    let y_mean = x.iter().sum::<f32>() / n;
    let mut num = 0.0f32;
    let mut den = 0.0f32;
    for (i, &y) in x.iter().enumerate() {
        let t = i as f32 - t_mean;
        num += t * (y - y_mean);
        den += t * t;
    }
    let slope = if den > 1e-9 { num / den } else { 0.0 };
    x.iter()
        .enumerate()
        .map(|(i, &y)| y - (y_mean + slope * (i as f32 - t_mean)))
        .collect()
}

/// Biased-normalized autocorrelation of `x` (already mean-centered by this function) at lags
/// `0..=max_lag`, each in `[-1, 1]` (1.0 at lag 0 by construction unless the series is constant).
fn autocorrelation(x: &[f32], max_lag: usize) -> Vec<f32> {
    let n = x.len();
    let mean = x.iter().sum::<f32>() / n as f32;
    let centered: Vec<f32> = x.iter().map(|v| v - mean).collect();
    let var: f32 = centered.iter().map(|v| v * v).sum();
    (0..=max_lag.min(n.saturating_sub(1)))
        .map(|lag| {
            if var < 1e-9 {
                return 0.0;
            }
            let s: f32 = (0..n - lag).map(|i| centered[i] * centered[i + lag]).sum();
            s / var
        })
        .collect()
}

fn mean_std(x: &[f32]) -> (f32, f32) {
    let n = x.len() as f32;
    let mean = x.iter().sum::<f32>() / n;
    let var = x.iter().map(|v| (v - mean) * (v - mean)).sum::<f32>() / n;
    (mean, var.sqrt())
}

fn report(args: &Args, samples: &[Sample]) {
    println!(
        "config: {:.1}% Nc, slurry={} ({:.1} Pa*s), mu_wall={:.2}, mu_ball={:.2}, lifters={}, seed={}",
        args.percent_critical,
        if args.slurry_enabled { "on" } else { "off" },
        args.viscosity_pa_s,
        args.friction_ball_wall,
        args.friction_ball_ball,
        args.lifters_count,
        args.seed,
    );

    let centroid_raw: Vec<f32> = samples.iter().map(|s| s.centroid_deg).collect();
    let centroid = unwrap_degrees(&centroid_raw);
    let (centroid_mean, centroid_std) = mean_std(&centroid);
    let centroid_p2p = centroid.iter().cloned().fold(f32::MIN, f32::max)
        - centroid.iter().cloned().fold(f32::MAX, f32::min);
    println!(
        "centroid angle (deg from vertical): mean={centroid_mean:.2} std={centroid_std:.3} p2p={centroid_p2p:.3}"
    );

    let slip: Vec<f32> = samples
        .iter()
        .map(|s| s.wall_slip_ratio)
        .filter(|v| v.is_finite())
        .collect();
    if !slip.is_empty() {
        let (slip_mean, slip_std) = mean_std(&slip);
        println!("wall slip ratio: mean={slip_mean:.4} std={slip_std:.4}");
    } else {
        println!("wall slip ratio: no wall-adjacent balls sampled");
    }

    // Work on the *rate* of centroid motion, linearly detrended, not the raw angle: the raw
    // signal's slow settle-in drift otherwise dominates the autocorrelation at short lags and
    // reads as spurious "periodicity" (see `finite_difference_rate`'s doc comment). A genuine
    // relaxation oscillation (stick-slip surging, or slumping) shows up here as the rate itself
    // swinging between two regimes once per cycle.
    let rate = linear_detrend(&finite_difference_rate(&centroid, 60.0));
    let min_lag = 12usize; // ignore lags under 0.2 s: not a meaningful bulk-charge timescale here.
    let max_lag = (rate.len() / 2).max(min_lag + 1);
    let ac = autocorrelation(&rate, max_lag);

    // A genuine oscillation's autocorrelation dips to a negative trough at roughly half its
    // period, then rises back to a secondary positive peak at roughly the full period -- this
    // trough-then-rebound shape is what distinguishes real periodicity from a one-off transient or
    // leftover smooth drift, which just decays toward zero without a clear rebound.
    let trough = ac
        .iter()
        .enumerate()
        .skip(min_lag)
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap());
    let rebound = trough.and_then(|(trough_lag, _)| {
        ac.iter()
            .enumerate()
            .skip(trough_lag + 1)
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
    });

    let detected_period_s: Option<f32> = match (trough, rebound) {
        (Some((trough_lag, &trough_val)), Some((rebound_lag, &rebound_val)))
            if trough_val <= -0.3 && rebound_val >= 0.3 =>
        {
            let period_s = rebound_lag as f32 / 60.0;
            println!(
                "detrended-rate autocorr: trough {trough_val:.3} at lag {trough_lag} ({:.3} s), \
                 rebound {rebound_val:.3} at lag {rebound_lag} ({period_s:.3} s)",
                trough_lag as f32 / 60.0
            );
            println!(
                "=> periodic signal DETECTED, dominant period ~= {period_s:.3} s (trough<=-0.3 and rebound>=0.3)"
            );
            Some(period_s)
        }
        (Some((trough_lag, &trough_val)), rebound) => {
            let rebound_desc = rebound
                .map(|(lag, &v)| format!("{v:.3} at lag {lag}"))
                .unwrap_or_else(|| "none".to_string());
            println!(
                "detrended-rate autocorr: trough {trough_val:.3} at lag {trough_lag}, best rebound {rebound_desc}"
            );
            println!("=> no clear periodic signal (trough/rebound below the +-0.3 threshold)");
            None
        }
        _ => {
            println!("=> no clear periodic signal (autocorrelation too flat to assess)");
            None
        }
    };

    // Amplitude: median peak-to-peak swing of the (linearly detrended) *raw* centroid angle,
    // measured per detected oscillation period -- distinct from `centroid_p2p` above, which is a
    // single global peak-to-peak over the whole measurement window and can be inflated by settle-in
    // drift or a single outlier swing. Chopping into windows of the detected period and taking the
    // median is more robust to both of those, and gives a number that maps directly to "how big is
    // one typical surge swing", which is what a large-amplitude-surging search cares about.
    match detected_period_s {
        Some(period_s) if period_s > 0.0 => {
            let detrended_centroid = linear_detrend(&centroid);
            let window_len = (period_s * 60.0).round() as usize;
            if let Some(n_windows) = detrended_centroid.len().checked_div(window_len) {
                if n_windows < 2 {
                    println!("amplitude: n/a (no periodic signal / too few periods)");
                } else {
                    let mut swings: Vec<f32> = (0..n_windows)
                        .map(|w| {
                            let chunk = &detrended_centroid[w * window_len..(w + 1) * window_len];
                            chunk.iter().cloned().fold(f32::MIN, f32::max)
                                - chunk.iter().cloned().fold(f32::MAX, f32::min)
                        })
                        .collect();
                    swings.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let median = if swings.len() % 2 == 1 {
                        swings[swings.len() / 2]
                    } else {
                        (swings[swings.len() / 2 - 1] + swings[swings.len() / 2]) / 2.0
                    };
                    println!("amplitude (median swing per period): {median:.2} deg");
                }
            } else {
                println!("amplitude: n/a (no periodic signal / too few periods)");
            }
        }
        _ => println!("amplitude: n/a (no periodic signal / too few periods)"),
    }

    // Sanity check against a physical-pendulum estimate: if the observed period tracks
    // `2*pi*sqrt(R/g)` (R = charge centroid's mean distance from the drum axis) rather than the
    // drum's own rotation period, the oscillation is better explained as a gravity-driven sloshing
    // mode of the charge's centroid (a pendulum bob at radius R) than as a wall stick-slip cycle
    // tied to `omega`, regardless of what the wall-slip-ratio phase relationship above suggests.
    let radii: Vec<f32> = samples
        .iter()
        .map(|s| s.centroid_radius_m)
        .filter(|v| v.is_finite())
        .collect();
    if !radii.is_empty() {
        let (radius_mean, _) = mean_std(&radii);
        let pendulum_period_s = std::f32::consts::TAU * (radius_mean / 9.81).sqrt();
        println!(
            "pendulum sanity check: mean centroid radius={radius_mean:.3} m => 2*pi*sqrt(R/g) = {pendulum_period_s:.3} s"
        );
    }

    let toe_count = samples.iter().filter(|s| s.toe_deg.is_some()).count();
    println!(
        "toe/shoulder defined in {toe_count}/{} samples (centrifuged or too-few-balls frames report None)",
        samples.len()
    );

    println!(
        "hint: in-phase centroid + wall-slip oscillation => surging (stick-slip against the wall); \
         flat wall slip with saw-toothing centroid/shoulder => slumping (bed-surface avalanching)"
    );
}
