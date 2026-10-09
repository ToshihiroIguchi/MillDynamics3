//! Relative grinding-proxy runner for scale-up studies (CPU only, WASM-safe, deterministic).
//!
//! Grinding itself is not simulated. What is measured instead are the quantities mechanistic
//! grinding models are driven by: the distribution of **stress energies** (`SE`, the relative
//! kinetic energy of an impact computed with the *3D sphere* mass, so it does not depend on the 2D
//! unit-depth disc convention), their frequency per ball, and the power drawn. Everything is also
//! reported dimensionless (`SE / (m g d)`, `m` the 3D ball mass, `d` the ball diameter), which is
//! what makes conditions of different size comparable. Results are meant for **relative**
//! comparison between conditions run with this same code; absolute values carry the 2D-model bias.
//!
//! Conventions:
//! - Impacts are split into **drop** type (normal approach speed >= tangential slip speed) and
//!   **shear** type (slip dominates), for ball-ball and ball-wall separately, because their
//!   scaling with drum/ball size differs (drop: `SE/(m g d) ~ D/d`; shear: `~ d/D`).
//! - "Per ball" powers use a one-ball-thick slab (depth = ball diameter): the 2D per-metre-depth
//!   power is multiplied by `d` and divided by the ball count.
//! - Numerical settings are ball-relative: sub-step `dt = dt_ratio * sqrt(d / g)`.
//! - The slurry is not resolved. Wet impact damping is the lubrication-cushion law of
//!   [`crate::dem::wet_cushion_restitution`] (`media.wet_cushion_viscosity_pa_s`); particle
//!   slurry, drag and buoyancy are *not* included.

use crate::dem::{DemState, ImpactKind, ImpactRecord};
use crate::geometry::Drum;
use crate::params::{CoarseGrainingMode, Params, SpeedMode};
use serde::{Deserialize, Serialize};

const G: f64 = 9.81;

/// Number of log-spaced bins of the dimensionless stress-energy histogram.
pub const SE_BINS: usize = 28;
const SE_BINS_PER_DECADE: f64 = 4.0;
const SE_LOG10_MIN: f64 = -3.0;
/// Impact classes, in the order used by every per-class array.
pub const CLASS_NAMES: [&str; 4] = [
    "ball_ball_drop",
    "ball_ball_shear",
    "ball_wall_drop",
    "ball_wall_shear",
];

/// Edges of the dimensionless stress-energy bins (`SE / (m g d)`), `SE_BINS + 1` values.
pub fn se_bin_edges() -> Vec<f64> {
    (0..=SE_BINS)
        .map(|i| 10f64.powf(SE_LOG10_MIN + i as f64 / SE_BINS_PER_DECADE))
        .collect()
}

fn se_bin(e_hat: f64) -> usize {
    if e_hat.is_nan() || e_hat <= 0.0 {
        return 0;
    }
    let t = (e_hat.log10() - SE_LOG10_MIN) * SE_BINS_PER_DECADE;
    (t.max(0.0) as usize).min(SE_BINS - 1)
}

/// Run configuration. All fields have defaults (`#[serde(default)]`), so a partial JSON works.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ScaleupConfig {
    /// Revolutions simulated before recording starts (charge settling).
    pub warmup_revolutions: f32,
    /// Revolutions recorded; use a whole number (the charge oscillates once per revolution).
    pub measure_revolutions: f32,
    /// Recording intervals per revolution (time-series resolution).
    pub samples_per_revolution: u32,
    /// Sub-step as a multiple of `sqrt(d / g)`.
    pub dt_ratio: f32,
    /// XPBD contact iterations per sub-step.
    pub dem_iterations: u32,
    /// Lubrication-cushion viscosity (Pa*s); `0` keeps `params.media.wet_cushion_viscosity_pa_s`.
    pub wet_viscosity_pa_s: f32,
    /// Absolute stress-energy thresholds (J) for the "effective" impact power. Empty selects
    /// `[0.03, 0.1, 0.3, 1, 3, 10] * m g d`.
    pub thresholds_j: Vec<f64>,
    /// Seeds of an ensemble; empty uses `params.simulation.seed` only.
    pub seeds: Vec<u64>,
}

impl Default for ScaleupConfig {
    fn default() -> Self {
        Self {
            warmup_revolutions: 4.0,
            measure_revolutions: 6.0,
            samples_per_revolution: 8,
            dt_ratio: 0.025,
            dem_iterations: 12,
            wet_viscosity_pa_s: 0.0,
            thresholds_j: Vec::new(),
            seeds: Vec::new(),
        }
    }
}

/// One recording interval (a fraction of a revolution). Powers are W per ball (one-ball-thick slab).
#[derive(Debug, Clone, Serialize)]
pub struct Interval {
    pub t_s: f64,
    pub revolutions: f64,
    pub wall_power_w: f64,
    pub dissipated_power_w: f64,
    /// Impacts per ball per second, per class ([`CLASS_NAMES`]).
    pub impacts_per_ball_per_s: [f64; 4],
    /// Mean dimensionless stress energy per class (0 when the class had no impacts).
    pub mean_se_scaled: [f64; 4],
    /// Rate (W per ball) of stress energy carried by impacts above each threshold.
    pub effective_power_w: Vec<f64>,
    /// Cumulative wall energy per ball since recording started (J).
    pub cum_wall_energy_j: f64,
    /// Cumulative stress energy per ball above each threshold since recording started (J).
    pub cum_effective_energy_j: Vec<f64>,
    /// Dimensionless stress-energy histogram counts per class (`SE_BINS` bins, see [`se_bin_edges`]).
    pub histogram: Vec<Vec<u32>>,
}

/// Scalar summary of one run (means over the recorded window).
#[derive(Debug, Clone, Serialize)]
pub struct Summary {
    pub drum_diameter_m: f64,
    pub ball_diameter_m: f64,
    pub dd_ratio: f64,
    pub n_balls: usize,
    pub percent_critical: f64,
    pub rpm: f64,
    pub ball_mass_3d_kg: f64,
    /// `m g d` (J): the energy unit used by every `*_scaled` quantity.
    pub scale_energy_j: f64,
    pub dt_s: f64,
    pub max_substep_displacement_over_diameter: f64,
    /// `rho_ball * sqrt(g D) * d / (9 mu)`; `None` for a dry run.
    pub reference_stokes_number: Option<f64>,
    pub wall_power_per_ball_w: f64,
    pub dissipated_power_per_ball_w: f64,
    pub specific_power_w_per_kg: f64,
    pub impacts_per_ball_per_s: [f64; 4],
    pub impacts_total_per_ball_per_s: f64,
    pub mean_se_scaled: [f64; 4],
    pub mean_se_scaled_drop: f64,
    pub mean_se_scaled_shear: f64,
    /// Mean dimensionless stress energy over all impacts (Kwade-type stress intensity).
    pub stress_intensity_scaled: f64,
    pub thresholds_j: Vec<f64>,
    pub effective_power_per_ball_w: Vec<f64>,
    /// Fraction of wall power carried by impacts above each threshold.
    pub effective_power_fraction: Vec<f64>,
    /// Least-squares slope of wall power over the window, relative to its mean, per revolution.
    pub power_drift_rel_per_rev: f64,
    /// Standard deviation of interval wall power relative to its mean.
    pub power_rel_std: f64,
}

impl Summary {
    /// Named scalars used for ensembles, sweeps, sensitivity and the UI's ratio table.
    pub fn scalars(&self) -> Vec<(String, f64)> {
        let mut v: Vec<(String, f64)> = vec![
            ("wall_power_per_ball_w".into(), self.wall_power_per_ball_w),
            (
                "dissipated_power_per_ball_w".into(),
                self.dissipated_power_per_ball_w,
            ),
            (
                "specific_power_w_per_kg".into(),
                self.specific_power_w_per_kg,
            ),
            (
                "impacts_total_per_ball_per_s".into(),
                self.impacts_total_per_ball_per_s,
            ),
            ("mean_se_scaled_drop".into(), self.mean_se_scaled_drop),
            ("mean_se_scaled_shear".into(), self.mean_se_scaled_shear),
            (
                "stress_intensity_scaled".into(),
                self.stress_intensity_scaled,
            ),
        ];
        for (i, name) in CLASS_NAMES.iter().enumerate() {
            v.push((
                format!("impacts_per_ball_per_s.{name}"),
                self.impacts_per_ball_per_s[i],
            ));
        }
        for (i, thr) in self.thresholds_j.iter().enumerate() {
            v.push((
                format!("effective_power_per_ball_w.thr{i}({thr:.3e}J)"),
                self.effective_power_per_ball_w[i],
            ));
            v.push((
                format!("effective_power_fraction.thr{i}"),
                self.effective_power_fraction[i],
            ));
        }
        v
    }
}

/// Everything one run produces.
#[derive(Debug, Clone, Serialize)]
pub struct ScaleupResult {
    pub seed: u64,
    pub summary: Summary,
    pub intervals: Vec<Interval>,
    /// Pooled histogram counts over the recorded window, per class.
    pub histogram_total: Vec<Vec<u64>>,
    /// Dimensionless bin edges shared by every histogram.
    pub se_bin_edges: Vec<f64>,
    pub class_names: Vec<String>,
}

impl ScaleupResult {
    /// The interval time series as CSV (scalars only; histograms are in the JSON).
    pub fn intervals_csv(&self) -> String {
        let mut s = String::from("t_s,revolutions,wall_power_w,dissipated_power_w");
        for n in CLASS_NAMES {
            s.push_str(&format!(",impacts_per_ball_per_s.{n}"));
        }
        for n in CLASS_NAMES {
            s.push_str(&format!(",mean_se_scaled.{n}"));
        }
        let nthr = self.summary.thresholds_j.len();
        for i in 0..nthr {
            s.push_str(&format!(",effective_power_w.thr{i}"));
        }
        s.push_str(",cum_wall_energy_j");
        for i in 0..nthr {
            s.push_str(&format!(",cum_effective_energy_j.thr{i}"));
        }
        s.push('\n');
        for iv in &self.intervals {
            s.push_str(&format!(
                "{},{},{},{}",
                iv.t_s, iv.revolutions, iv.wall_power_w, iv.dissipated_power_w
            ));
            for x in iv.impacts_per_ball_per_s {
                s.push_str(&format!(",{x}"));
            }
            for x in iv.mean_se_scaled {
                s.push_str(&format!(",{x}"));
            }
            for x in &iv.effective_power_w {
                s.push_str(&format!(",{x}"));
            }
            s.push_str(&format!(",{}", iv.cum_wall_energy_j));
            for x in &iv.cum_effective_energy_j {
                s.push_str(&format!(",{x}"));
            }
            s.push('\n');
        }
        s
    }

    /// Time (s since recording started) at which the cumulative effective energy above threshold
    /// `thr_index` reaches `target_j_per_ball`, linearly interpolated; `None` if never reached.
    /// The ratio of two conditions' times at the same target is the proxy for relative grinding time.
    pub fn time_to_cumulative_effective_energy(
        &self,
        thr_index: usize,
        target_j_per_ball: f64,
    ) -> Option<f64> {
        let t0 = self.intervals.first().map(|i| i.t_s)?;
        let mut prev_t = t0 - interval_len(&self.intervals);
        let mut prev_e = 0.0;
        for iv in &self.intervals {
            let e = *iv.cum_effective_energy_j.get(thr_index)?;
            if e >= target_j_per_ball {
                let f = if e > prev_e {
                    (target_j_per_ball - prev_e) / (e - prev_e)
                } else {
                    1.0
                };
                return Some(prev_t + f * (iv.t_s - prev_t) - (t0 - interval_len(&self.intervals)));
            }
            prev_t = iv.t_s;
            prev_e = e;
        }
        None
    }
}

fn interval_len(intervals: &[Interval]) -> f64 {
    if intervals.len() >= 2 {
        intervals[1].t_s - intervals[0].t_s
    } else {
        intervals.first().map(|i| i.t_s).unwrap_or(0.0)
    }
}

/// Prepares the parameters shared by every run (no coarse-graining, validated).
fn prepare(params: &Params, cfg: &ScaleupConfig) -> Result<Params, String> {
    let mut p = *params;
    p.simulation.coarse_graining_mode = CoarseGrainingMode::Off;
    p.simulation.dem_iterations = cfg.dem_iterations.clamp(1, 20);
    if cfg.wet_viscosity_pa_s > 0.0 {
        p.media.wet_cushion_viscosity_pa_s = cfg.wet_viscosity_pa_s;
    }
    // Ball-relative versions of the solver's absolute-speed constants, so that geometrically
    // similar mills of different size run identically (see `scaleup_probe --cases sim`).
    const FROUDE: f32 = 0.1428;
    p.media.restitution_threshold_froude = FROUDE;
    p.media.friction_velocity_scale_m_s = FROUDE * (9.81 * p.media.ball_diameter_m).sqrt();
    p.validate()?;
    let n = p.true_ball_count();
    if n < 2.0 {
        return Err("fewer than 2 balls: raise fill fraction or lower ball diameter".into());
    }
    if n > 50_000.0 {
        return Err(format!(
            "{n:.0} balls is too many for an in-browser run (limit 50000): use a smaller D/d"
        ));
    }
    if p.mill.omega_magnitude() <= 1e-6 {
        return Err("mill speed must be > 0".into());
    }
    if !(cfg.measure_revolutions >= 1.0 && cfg.warmup_revolutions >= 0.0) {
        return Err("measure_revolutions must be >= 1 and warmup_revolutions >= 0".into());
    }
    if cfg.samples_per_revolution == 0 || cfg.dt_ratio.is_nan() || cfg.dt_ratio <= 0.0 {
        return Err("samples_per_revolution and dt_ratio must be > 0".into());
    }
    Ok(p)
}

#[derive(Default, Clone)]
struct Accum {
    count: [u64; 4],
    sum_se_scaled: [f64; 4],
    hist: [[u32; SE_BINS]; 4],
    se_above: Vec<f64>,
}

impl Accum {
    fn new(nthr: usize) -> Self {
        Self {
            se_above: vec![0.0; nthr],
            ..Default::default()
        }
    }
}

/// Runs one deterministic simulation and returns its summary and time series.
pub fn run(params: &Params, cfg: &ScaleupConfig) -> Result<ScaleupResult, String> {
    run_seeded(params, cfg, params.simulation.seed)
}

fn run_seeded(params: &Params, cfg: &ScaleupConfig, seed: u64) -> Result<ScaleupResult, String> {
    let mut p = prepare(params, cfg)?;
    p.simulation.seed = seed;
    let effective = p.effective_media();
    let radius_m = p.mill.radius_m();
    let omega = p.mill.omega();
    let drum = Drum::new(radius_m, omega, p.lifters);
    let mut state = DemState::new(&effective, radius_m, seed);
    let n_balls = state.balls.len();

    let d = effective.diameter_m as f64;
    let r = d * 0.5;
    let rho = effective.density_kg_m3 as f64;
    let m3 = rho * 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
    let scale_j = m3 * G * d;
    let thresholds_j: Vec<f64> = if cfg.thresholds_j.is_empty() {
        [0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
            .iter()
            .map(|k| k * scale_j)
            .collect()
    } else {
        cfg.thresholds_j.clone()
    };
    let nthr = thresholds_j.len();

    let dt = cfg.dt_ratio as f64 * (d / G).sqrt();
    let period = std::f64::consts::TAU / omega.abs() as f64;
    let interval_s = period / cfg.samples_per_revolution as f64;
    let warmup_s = cfg.warmup_revolutions as f64 * period;
    let n_intervals =
        (cfg.measure_revolutions * cfg.samples_per_revolution as f32).round() as usize;
    let slab_depth_m = d;
    let per_ball = slab_depth_m / n_balls as f64;

    let mut angle = 0.0f32;
    let mut t = 0.0f64;
    let dt32 = dt as f32;
    // Warm-up: no recording.
    while t < warmup_s {
        state.step_recorded(
            &drum,
            angle,
            &p.media,
            p.simulation.dem_iterations,
            dt32,
            None,
            None,
        );
        angle = (angle + omega * dt32).rem_euclid(std::f32::consts::TAU);
        t += dt;
    }

    let mut log: Vec<ImpactRecord> = Vec::new();
    let mut intervals: Vec<Interval> = Vec::with_capacity(n_intervals);
    let mut histogram_total = vec![vec![0u64; SE_BINS]; 4];
    let mut acc = Accum::new(nthr);
    let (mut wall_j, mut diss_j) = (0.0f64, 0.0f64);
    let mut interval_t = 0.0f64;
    let mut rec_t = 0.0f64;
    let (mut cum_wall, mut cum_eff) = (0.0f64, vec![0.0f64; nthr]);
    let mut max_disp = 0.0f64;
    let masses = [m3 * 0.5, m3 * 0.5, m3, m3]; // reduced mass per class (bb, bb, bw, bw)

    while intervals.len() < n_intervals {
        log.clear();
        let stats = state.step_recorded(
            &drum,
            angle,
            &p.media,
            p.simulation.dem_iterations,
            dt32,
            None,
            Some(&mut log),
        );
        angle = (angle + omega * dt32).rem_euclid(std::f32::consts::TAU);
        wall_j += stats.wall_work_j as f64;
        diss_j += stats.dissipated_energy_j as f64;
        max_disp = max_disp.max(stats.max_substep_displacement_over_diameter as f64);
        for rec in &log {
            let bb = rec.kind == ImpactKind::BallBall;
            let (vn, vs) = (rec.normal_speed as f64, rec.slip_speed as f64);
            let shear = vs > vn;
            let class = match (bb, shear) {
                (true, false) => 0,
                (true, true) => 1,
                (false, false) => 2,
                (false, true) => 3,
            };
            let m_red = masses[class];
            let se = 0.5 * m_red * (vn * vn + vs * vs);
            let e_hat = se / scale_j;
            acc.count[class] += 1;
            acc.sum_se_scaled[class] += e_hat;
            let b = se_bin(e_hat);
            acc.hist[class][b] += 1;
            histogram_total[class][b] += 1;
            for (k, thr) in thresholds_j.iter().enumerate() {
                if se >= *thr {
                    acc.se_above[k] += se;
                }
            }
        }
        interval_t += dt;
        rec_t += dt;
        if rec_t >= (intervals.len() + 1) as f64 * interval_s {
            let len = interval_t;
            let wall_power = wall_j * per_ball / len;
            let mut impacts = [0.0; 4];
            let mut mean_se = [0.0; 4];
            for c in 0..4 {
                impacts[c] = acc.count[c] as f64 / n_balls as f64 / len;
                if acc.count[c] > 0 {
                    mean_se[c] = acc.sum_se_scaled[c] / acc.count[c] as f64;
                }
            }
            let eff_power: Vec<f64> = acc
                .se_above
                .iter()
                .map(|e| e / n_balls as f64 / len)
                .collect();
            cum_wall += wall_j * per_ball;
            for (c, e) in cum_eff.iter_mut().zip(&acc.se_above) {
                *c += e / n_balls as f64;
            }
            intervals.push(Interval {
                t_s: rec_t,
                revolutions: rec_t / period,
                wall_power_w: wall_power,
                dissipated_power_w: diss_j * per_ball / len,
                impacts_per_ball_per_s: impacts,
                mean_se_scaled: mean_se,
                effective_power_w: eff_power,
                cum_wall_energy_j: cum_wall,
                cum_effective_energy_j: cum_eff.clone(),
                histogram: acc.hist.iter().map(|h| h.to_vec()).collect(),
            });
            acc = Accum::new(nthr);
            wall_j = 0.0;
            diss_j = 0.0;
            interval_t = 0.0;
        }
    }
    if intervals.is_empty() {
        return Err("no recording interval completed".into());
    }

    let summary = summarise(
        &p,
        &intervals,
        &histogram_total,
        n_balls,
        m3,
        scale_j,
        dt,
        max_disp,
        thresholds_j,
        period,
    );
    Ok(ScaleupResult {
        seed,
        summary,
        intervals,
        histogram_total,
        se_bin_edges: se_bin_edges(),
        class_names: CLASS_NAMES.iter().map(|s| s.to_string()).collect(),
    })
}

#[allow(clippy::too_many_arguments)]
fn summarise(
    p: &Params,
    intervals: &[Interval],
    hist_total: &[Vec<u64>],
    n_balls: usize,
    m3: f64,
    scale_j: f64,
    dt: f64,
    max_disp: f64,
    thresholds_j: Vec<f64>,
    period: f64,
) -> Summary {
    let n = intervals.len() as f64;
    let mean = |f: &dyn Fn(&Interval) -> f64| intervals.iter().map(f).sum::<f64>() / n;
    let wall = mean(&|i| i.wall_power_w);
    let nthr = thresholds_j.len();
    let eff: Vec<f64> = (0..nthr)
        .map(|k| mean(&|i| i.effective_power_w[k]))
        .collect();
    let eff_frac: Vec<f64> = eff
        .iter()
        .map(|e| if wall > 0.0 { e / wall } else { 0.0 })
        .collect();
    let mut impacts = [0.0; 4];
    for (c, imp) in impacts.iter_mut().enumerate() {
        *imp = mean(&|i| i.impacts_per_ball_per_s[c]);
    }
    // Pooled mean stress energy per class from the (count-weighted) interval means.
    let mut mean_se = [0.0; 4];
    let mut total_count = [0.0f64; 4];
    for (c, ts) in total_count.iter_mut().enumerate() {
        *ts = hist_total[c].iter().sum::<u64>() as f64;
        let w: f64 = intervals.iter().map(|i| i.impacts_per_ball_per_s[c]).sum();
        if w > 0.0 {
            mean_se[c] = intervals
                .iter()
                .map(|i| i.mean_se_scaled[c] * i.impacts_per_ball_per_s[c])
                .sum::<f64>()
                / w;
        }
    }
    let pooled = |cs: [usize; 2]| {
        let w: f64 = cs.iter().map(|&c| impacts[c]).sum();
        if w > 0.0 {
            cs.iter().map(|&c| mean_se[c] * impacts[c]).sum::<f64>() / w
        } else {
            0.0
        }
    };
    let total_rate: f64 = impacts.iter().sum();
    let stress_intensity = if total_rate > 0.0 {
        (0..4).map(|c| mean_se[c] * impacts[c]).sum::<f64>() / total_rate
    } else {
        0.0
    };

    // Drift and scatter of the wall power over the window.
    let revs: Vec<f64> = intervals.iter().map(|i| i.t_s / period).collect();
    let rev_mean = revs.iter().sum::<f64>() / n;
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (x, i) in revs.iter().zip(intervals) {
        sxy += (x - rev_mean) * (i.wall_power_w - wall);
        sxx += (x - rev_mean) * (x - rev_mean);
    }
    let slope = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    let var = intervals
        .iter()
        .map(|i| (i.wall_power_w - wall).powi(2))
        .sum::<f64>()
        / n;
    let rho = p.media.density_kg_m3 as f64;
    let d = p.media.ball_diameter_m as f64;
    let big_d = p.mill.diameter_m as f64;
    let mu = p.media.wet_cushion_viscosity_pa_s as f64;
    Summary {
        drum_diameter_m: big_d,
        ball_diameter_m: d,
        dd_ratio: big_d / d,
        n_balls,
        percent_critical: p.mill.percent_critical() as f64,
        rpm: p.mill.rpm() as f64,
        ball_mass_3d_kg: m3,
        scale_energy_j: scale_j,
        dt_s: dt,
        max_substep_displacement_over_diameter: max_disp,
        reference_stokes_number: if mu > 0.0 {
            Some(rho * (G * big_d).sqrt() * d / (9.0 * mu))
        } else {
            None
        },
        wall_power_per_ball_w: wall,
        dissipated_power_per_ball_w: mean(&|i| i.dissipated_power_w),
        specific_power_w_per_kg: wall / m3,
        impacts_per_ball_per_s: impacts,
        impacts_total_per_ball_per_s: total_rate,
        mean_se_scaled: mean_se,
        mean_se_scaled_drop: pooled([0, 2]),
        mean_se_scaled_shear: pooled([1, 3]),
        stress_intensity_scaled: stress_intensity,
        thresholds_j,
        effective_power_per_ball_w: eff,
        effective_power_fraction: eff_frac,
        power_drift_rel_per_rev: if wall.abs() > 0.0 { slope / wall } else { 0.0 },
        power_rel_std: if wall.abs() > 0.0 {
            var.sqrt() / wall.abs()
        } else {
            0.0
        },
    }
}

/// Mean and standard error of every [`Summary::scalars`] entry over an ensemble.
#[derive(Debug, Clone, Serialize)]
pub struct EnsembleResult {
    pub seeds: Vec<u64>,
    pub scalar_names: Vec<String>,
    pub mean: Vec<f64>,
    pub std_err: Vec<f64>,
    /// The first member's full result (time series, histogram); the others are summary-only.
    pub first: ScaleupResult,
    pub summaries: Vec<Summary>,
}

/// Runs the seeds in `cfg.seeds` (or the single `params.simulation.seed`).
pub fn run_ensemble(params: &Params, cfg: &ScaleupConfig) -> Result<EnsembleResult, String> {
    let seeds = if cfg.seeds.is_empty() {
        vec![params.simulation.seed]
    } else {
        cfg.seeds.clone()
    };
    let mut results = Vec::with_capacity(seeds.len());
    for &s in &seeds {
        results.push(run_seeded(params, cfg, s)?);
    }
    let summaries: Vec<Summary> = results.iter().map(|r| r.summary.clone()).collect();
    let names: Vec<String> = summaries[0].scalars().into_iter().map(|(n, _)| n).collect();
    let cols: Vec<Vec<f64>> = (0..names.len())
        .map(|k| summaries.iter().map(|s| s.scalars()[k].1).collect())
        .collect();
    let n = seeds.len() as f64;
    let mean: Vec<f64> = cols.iter().map(|c| c.iter().sum::<f64>() / n).collect();
    let std_err: Vec<f64> = cols
        .iter()
        .zip(&mean)
        .map(|(c, m)| {
            if seeds.len() < 2 {
                0.0
            } else {
                let var = c.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0);
                (var / n).sqrt()
            }
        })
        .collect();
    let first = results.swap_remove(0);
    Ok(EnsembleResult {
        seeds,
        scalar_names: names,
        mean,
        std_err,
        first,
        summaries,
    })
}

/// A power law `y = prefactor * x^exponent` fitted in log-log space.
#[derive(Debug, Clone, Serialize)]
pub struct PowerLawFit {
    pub name: String,
    pub exponent: f64,
    pub prefactor: f64,
    pub r2: f64,
    /// Exponent expected from the drop/shear velocity scaling, where one exists.
    pub expected_exponent: Option<f64>,
}

impl PowerLawFit {
    pub fn predict(&self, x: f64) -> f64 {
        self.prefactor * x.powf(self.exponent)
    }
}

fn fit_power_law(
    name: &str,
    xs: &[f64],
    ys: &[f64],
    expected_exponent: Option<f64>,
) -> Option<PowerLawFit> {
    let pts: Vec<(f64, f64)> = xs
        .iter()
        .zip(ys)
        .filter(|(x, y)| **x > 0.0 && **y > 0.0)
        .map(|(x, y)| (x.ln(), y.ln()))
        .collect();
    if pts.len() < 2 {
        return None;
    }
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = pts.iter().map(|p| (p.0 - mx).powi(2)).sum();
    if sxx <= 0.0 {
        return None;
    }
    let sxy: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let syy: f64 = pts.iter().map(|p| (p.1 - my).powi(2)).sum();
    let slope = sxy / sxx;
    let intercept = my - slope * mx;
    let r2 = if syy > 0.0 {
        sxy * sxy / (sxx * syy)
    } else {
        1.0
    };
    Some(PowerLawFit {
        name: name.to_string(),
        exponent: slope,
        prefactor: intercept.exp(),
        r2,
        expected_exponent,
    })
}

/// One drum-to-ball-diameter ratio of a sweep.
#[derive(Debug, Clone, Serialize)]
pub struct SweepPoint {
    pub dd_ratio: f64,
    pub n_balls: usize,
    pub ensemble: EnsembleResult,
}

/// A `D/d` sweep (ball diameter and every dimensionless setting fixed, drum diameter varied) and
/// the power-law fits used to extrapolate to a larger `D/d` without simulating its ball count.
#[derive(Debug, Clone, Serialize)]
pub struct SweepResult {
    pub points: Vec<SweepPoint>,
    pub fits: Vec<PowerLawFit>,
}

/// Runs `base` at each drum-to-ball ratio. The drum speed is held at the base `% critical`
/// (so the Froude number is constant), everything else unchanged.
pub fn sweep_dd(base: &Params, cfg: &ScaleupConfig, ratios: &[f64]) -> Result<SweepResult, String> {
    let pc = base.mill.percent_critical();
    let mut points = Vec::new();
    for &ratio in ratios {
        let mut p = *base;
        p.mill.diameter_m = (ratio * p.media.ball_diameter_m as f64) as f32;
        p.mill.speed_mode = SpeedMode::PercentCritical;
        p.mill.speed_value = pc;
        let ensemble = run_ensemble(&p, cfg)?;
        points.push(SweepPoint {
            dd_ratio: ratio,
            n_balls: ensemble.first.summary.n_balls,
            ensemble,
        });
    }
    let xs: Vec<f64> = points.iter().map(|p| p.dd_ratio).collect();
    let column = |name: &str| -> Vec<f64> {
        points
            .iter()
            .map(|p| {
                let k = p
                    .ensemble
                    .scalar_names
                    .iter()
                    .position(|n| n == name)
                    .unwrap_or(0);
                p.ensemble.mean[k]
            })
            .collect()
    };
    let mut fits = Vec::new();
    for (name, expected) in [
        ("mean_se_scaled_drop", Some(1.0)),
        ("mean_se_scaled_shear", Some(-1.0)),
        ("stress_intensity_scaled", None),
        ("impacts_total_per_ball_per_s", None),
        ("specific_power_w_per_kg", None),
        ("wall_power_per_ball_w", None),
    ] {
        if let Some(f) = fit_power_law(name, &xs, &column(name), expected) {
            fits.push(f);
        }
    }
    Ok(SweepResult { points, fits })
}

/// One sensitivity case: a labelled perturbation of the contact parameters.
#[derive(Debug, Clone, Serialize)]
pub struct SensitivityCase {
    pub label: String,
    pub scalar_names: Vec<String>,
    pub mean: Vec<f64>,
}

/// Re-runs `params` with each contact coefficient scaled by `1 +- rel` (restitution and friction,
/// ball-ball and ball-wall), so the caller can check that condition-to-condition *ratios* are
/// insensitive to the assumed values. The first case is the unperturbed base.
pub fn sensitivity(
    params: &Params,
    cfg: &ScaleupConfig,
    rel: f32,
) -> Result<Vec<SensitivityCase>, String> {
    type Setter = fn(&mut Params, f32);
    let knobs: [(&str, Setter); 4] = [
        ("restitution_ball_ball", |p, f| {
            p.media.restitution_ball_ball = (p.media.restitution_ball_ball * f).min(1.0)
        }),
        ("restitution_ball_wall", |p, f| {
            p.media.restitution_ball_wall = (p.media.restitution_ball_wall * f).min(1.0)
        }),
        ("friction_ball_ball", |p, f| {
            p.media.friction_ball_ball *= f;
            p.media.friction_ball_ball_static *= f;
        }),
        ("friction_ball_wall", |p, f| {
            p.media.friction_ball_wall *= f;
            p.media.friction_ball_wall_static *= f;
        }),
    ];
    let mut cases = Vec::new();
    let base = run_ensemble(params, cfg)?;
    cases.push(SensitivityCase {
        label: "base".into(),
        scalar_names: base.scalar_names.clone(),
        mean: base.mean.clone(),
    });
    for (name, set) in knobs {
        for (sign, f) in [("-", 1.0 - rel), ("+", 1.0 + rel)] {
            let mut p = *params;
            set(&mut p, f);
            let e = run_ensemble(&p, cfg)?;
            cases.push(SensitivityCase {
                label: format!("{name} {sign}{:.0}%", rel * 100.0),
                scalar_names: e.scalar_names,
                mean: e.mean,
            });
        }
    }
    Ok(cases)
}

/// Ball diameter that keeps the drop-type stress energy `SE ~ rho d^3 g D` of a lab condition
/// (`d_lab`, `D_lab`) when the drum grows to `D_target` at the same `% critical`:
/// `d = d_lab * (D_lab / D_target)^(1/3)`. A starting point for a candidate search, not a result.
pub fn matching_ball_diameter(d_lab_m: f64, drum_lab_m: f64, drum_target_m: f64) -> f64 {
    d_lab_m * (drum_lab_m / drum_target_m).powf(1.0 / 3.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quick_cfg() -> ScaleupConfig {
        ScaleupConfig {
            warmup_revolutions: 1.0,
            measure_revolutions: 2.0,
            ..Default::default()
        }
    }

    #[test]
    fn bins_are_log_spaced_and_cover_the_range() {
        let e = se_bin_edges();
        assert_eq!(e.len(), SE_BINS + 1);
        assert!((e[0] - 1e-3).abs() < 1e-9);
        assert!((e[SE_BINS] - 1e4).abs() / 1e4 < 1e-9);
        assert_eq!(se_bin(1e-9), 0);
        assert_eq!(se_bin(1e9), SE_BINS - 1);
        assert_eq!(se_bin(1.0), 12);
    }

    #[test]
    fn matching_diameter_keeps_drop_stress_energy() {
        // SE ~ d^3 D: tenfold drum -> d shrinks by 10^(1/3).
        let d = matching_ball_diameter(0.002, 0.063, 0.63);
        let ratio = 0.002f64.powi(3) * 0.063 / (d.powi(3) * 0.63);
        assert!((ratio - 1.0).abs() < 1e-9);
    }

    #[test]
    fn power_law_fit_recovers_an_exact_exponent() {
        let xs = [30.0, 60.0, 120.0];
        let ys: Vec<f64> = xs.iter().map(|x: &f64| 0.5 * x.powf(-1.0)).collect();
        let f = fit_power_law("t", &xs, &ys, Some(-1.0)).unwrap();
        assert!((f.exponent + 1.0).abs() < 1e-9);
        assert!((f.prefactor - 0.5).abs() < 1e-9);
        assert!((f.r2 - 1.0).abs() < 1e-9);
        assert!((f.predict(240.0) - 0.5 / 240.0).abs() < 1e-12);
    }

    #[test]
    fn default_scene_produces_a_sane_deterministic_summary() {
        let mut p = Params::default();
        p.slurry.enabled = false;
        let cfg = quick_cfg();
        let a = run(&p, &cfg).unwrap();
        let b = run(&p, &cfg).unwrap();
        assert_eq!(
            a.summary.wall_power_per_ball_w,
            b.summary.wall_power_per_ball_w
        );
        assert_eq!(a.intervals.len(), 16);
        let s = &a.summary;
        assert!(s.wall_power_per_ball_w > 0.0, "{s:?}");
        assert!(s.impacts_total_per_ball_per_s > 0.0, "{s:?}");
        assert!(s.stress_intensity_scaled > 0.0);
        // Histogram totals equal the per-interval sums.
        let from_intervals: u64 = a
            .intervals
            .iter()
            .flat_map(|i| i.histogram.iter().flatten())
            .map(|&c| c as u64)
            .sum();
        let pooled: u64 = a.histogram_total.iter().flatten().sum();
        assert_eq!(from_intervals, pooled);
        // Cumulative wall energy equals the sum of interval powers times the interval length.
        let dt = a.intervals[1].t_s - a.intervals[0].t_s;
        let e: f64 = a.intervals.iter().map(|i| i.wall_power_w * dt).sum();
        let last = a.intervals.last().unwrap().cum_wall_energy_j;
        assert!(
            (e - last).abs() / last.abs().max(1e-12) < 0.02,
            "{e} vs {last}"
        );
        assert!(a.intervals_csv().lines().count() == 17);
    }

    #[test]
    fn wet_cushion_removes_bouncing_energy() {
        let mut p = Params::default();
        p.slurry.enabled = false;
        let cfg = quick_cfg();
        let dry = run(&p, &cfg).unwrap();
        let wet_cfg = ScaleupConfig {
            wet_viscosity_pa_s: 1.0,
            ..quick_cfg()
        };
        let wet = run(&p, &wet_cfg).unwrap();
        assert!(wet.summary.reference_stokes_number.unwrap() < 20.0);
        eprintln!(
            "DRY {:?}
WET {:?}",
            dry.summary, wet.summary
        );
        assert!(
            wet.summary.impacts_total_per_ball_per_s < dry.summary.impacts_total_per_ball_per_s,
            "wet {} dry {}",
            wet.summary.impacts_total_per_ball_per_s,
            dry.summary.impacts_total_per_ball_per_s
        );
    }

    #[test]
    fn geometrically_similar_mills_give_the_same_dimensionless_stress_energy() {
        // Drum and ball scaled together at constant % critical: with ball-relative numerical
        // settings the rigid-contact model is exactly Froude-similar, so SE / (m g d) must agree.
        let make = |scale: f32| {
            let mut p = Params::default();
            p.slurry.enabled = false;
            let pc = p.mill.percent_critical();
            p.mill.speed_mode = SpeedMode::PercentCritical;
            p.mill.speed_value = pc;
            p.mill.diameter_m *= scale;
            p.media.ball_diameter_m *= scale;
            run(&p, &quick_cfg()).unwrap().summary
        };
        let (a, b) = (make(1.0), make(2.0));
        let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(y.abs());
        assert!(
            rel(a.mean_se_scaled_drop, b.mean_se_scaled_drop) < 0.1,
            "{} vs {}",
            a.mean_se_scaled_drop,
            b.mean_se_scaled_drop
        );
        assert!(rel(a.mean_se_scaled_shear, b.mean_se_scaled_shear) < 0.1);
        // Specific power scales as sqrt(length).
        let ratio = b.specific_power_w_per_kg / a.specific_power_w_per_kg;
        assert!((ratio - 2f64.sqrt()).abs() / 2f64.sqrt() < 0.1, "{ratio}");
    }

    #[test]
    fn prepare_rejects_absurd_ball_counts() {
        let mut p = Params::default();
        p.mill.diameter_m = 5.0;
        assert!(run(&p, &quick_cfg()).is_err());
    }
}
