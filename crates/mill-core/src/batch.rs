//! Headless repeat runs of the **same** [`Simulation`] the main screen runs, for averaging the
//! grinding-proxy statistics ([`crate::stress`]) over a whole number of revolutions and over seeds.
//!
//! Nothing here changes the physics or the numerical settings: `params` is used as given (slurry,
//! sub-steps, iterations, coarse-graining included), only `simulation.seed` is replaced. One call
//! runs one seed; the caller (the web worker) combines seeds with [`combine`], so it can report
//! progress and cancel between seeds.

use crate::{stress, Params, Simulation};
use serde::Serialize;

/// Window means of one seed's run.
#[derive(Debug, Clone, Serialize)]
pub struct SeedResult {
    pub seed: u64,
    /// Simulated seconds in the averaging window.
    pub window_s: f64,
    /// Coarse-graining factor `k` of the run; the statistics are meaningless when `k > 1`.
    pub coarse_graining_factor: f32,
    pub names: Vec<String>,
    pub values: Vec<f64>,
    /// Least-squares slope of the wall power over the window relative to its mean, per revolution
    /// (steady-state indicator).
    pub power_drift_rel_per_rev: f64,
}

/// Mean and standard error (n - 1) of one scalar over seeds.
#[derive(Debug, Clone, Serialize)]
pub struct Combined {
    pub names: Vec<String>,
    pub mean: Vec<f64>,
    pub std_err: Vec<f64>,
    pub seeds: usize,
}

const FRAME_S: f32 = 1.0 / 60.0;

/// Runs `params` with `seed` for `warmup_revolutions`, then averages the stress metrics over
/// `measure_revolutions` (use a whole number: the charge oscillates once per revolution).
pub fn run_seed(
    params: &Params,
    seed: u64,
    warmup_revolutions: f32,
    measure_revolutions: f32,
) -> Result<SeedResult, String> {
    let mut p = *params;
    p.simulation.seed = seed;
    let omega = p.mill.omega_magnitude() as f64;
    if omega <= 1e-6 {
        return Err("mill speed must be > 0".into());
    }
    if !(measure_revolutions >= 1.0 && warmup_revolutions >= 0.0) {
        return Err("measure_revolutions must be >= 1 and warmup_revolutions >= 0".into());
    }
    let rev_s = std::f64::consts::TAU / omega;
    let mut sim = Simulation::new(p)?;
    let warm_frames = (warmup_revolutions as f64 * rev_s * 60.0).round() as usize;
    let meas_frames = (measure_revolutions as f64 * rev_s * 60.0).round().max(1.0) as usize;
    for _ in 0..warm_frames {
        sim.step(FRAME_S);
    }

    let start = sim.metrics().stress;
    let names = names(&start);
    let mut sum = vec![0.0f64; names.len()];
    let mut powers = Vec::with_capacity(meas_frames);
    for _ in 0..meas_frames {
        sim.step(FRAME_S);
        let m = sim.metrics();
        for (s, v) in sum.iter_mut().zip(values(&m.stress)) {
            *s += v;
        }
        powers.push(m.stress.wall_power_per_ball_w as f64);
    }
    let end = sim.metrics();
    let n = meas_frames as f64;
    let window_s = n * FRAME_S as f64;
    let mut vals: Vec<f64> = sum.iter().map(|s| s / n).collect();
    // Exact window means from the cumulative (not smoothed) energies.
    let dm = &end.stress;
    vals.push(
        (dm.cumulative_specific_energy_j_per_kg - start.cumulative_specific_energy_j_per_kg)
            / window_s,
    );
    for (e1, e0) in dm
        .cumulative_effective_energy_per_ball_j
        .iter()
        .zip(&start.cumulative_effective_energy_per_ball_j)
    {
        vals.push((e1 - e0) / window_s);
    }
    let mut names = names;
    names.push("window_specific_power_w_per_kg".into());
    for i in 0..start.thresholds_scaled.len() {
        names.push(format!("window_effective_power_per_ball_w.thr{i}"));
    }

    Ok(SeedResult {
        seed,
        window_s,
        coarse_graining_factor: end.coarse_graining_factor,
        names,
        values: vals,
        power_drift_rel_per_rev: drift_per_rev(&powers, rev_s * 60.0),
    })
}

fn names(m: &stress::StressMetrics) -> Vec<String> {
    let mut v: Vec<String> = vec![
        "specific_power_w_per_kg".into(),
        "wall_power_per_ball_w".into(),
        "impacts_total_per_ball_per_s".into(),
        "mean_se_scaled_drop".into(),
        "mean_se_scaled_shear".into(),
        "stress_intensity_scaled".into(),
    ];
    for c in stress::CLASS_NAMES {
        v.push(format!("impacts_per_ball_per_s.{c}"));
    }
    for i in 0..m.thresholds_scaled.len() {
        v.push(format!("effective_power_per_ball_w.thr{i}"));
    }
    v
}

fn values(m: &stress::StressMetrics) -> Vec<f64> {
    let mut v = vec![
        m.specific_power_w_per_kg as f64,
        m.wall_power_per_ball_w as f64,
        m.impacts_total_per_ball_per_s as f64,
        m.mean_se_scaled_drop as f64,
        m.mean_se_scaled_shear as f64,
        m.stress_intensity_scaled as f64,
    ];
    v.extend(m.impacts_per_ball_per_s.iter().map(|&x| x as f64));
    v.extend(m.effective_power_per_ball_w.iter().map(|&x| x as f64));
    v
}

/// Slope of `y` against frame index, relative to the mean of `y`, per revolution
/// (`frames_per_rev` frames).
fn drift_per_rev(y: &[f64], frames_per_rev: f64) -> f64 {
    let n = y.len() as f64;
    if y.len() < 2 {
        return 0.0;
    }
    let mean = y.iter().sum::<f64>() / n;
    let xm = (n - 1.0) / 2.0;
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (i, v) in y.iter().enumerate() {
        sxy += (i as f64 - xm) * (v - mean);
        sxx += (i as f64 - xm) * (i as f64 - xm);
    }
    if mean.abs() > 0.0 && sxx > 0.0 {
        sxy / sxx * frames_per_rev / mean.abs()
    } else {
        0.0
    }
}

/// Mean and standard error over seeds. All results must share the same names.
pub fn combine(results: &[SeedResult]) -> Result<Combined, String> {
    let first = results.first().ok_or("no results to combine")?;
    if results.iter().any(|r| r.names != first.names) {
        return Err("results have different scalar names".into());
    }
    let n = results.len() as f64;
    let mut mean = Vec::new();
    let mut std_err = Vec::new();
    for k in 0..first.names.len() {
        let m = results.iter().map(|r| r.values[k]).sum::<f64>() / n;
        let se = if results.len() < 2 {
            0.0
        } else {
            let var = results
                .iter()
                .map(|r| (r.values[k] - m).powi(2))
                .sum::<f64>()
                / (n - 1.0);
            (var / n).sqrt()
        };
        mean.push(m);
        std_err.push(se);
    }
    Ok(Combined {
        names: first.names.clone(),
        mean,
        std_err,
        seeds: results.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dry() -> Params {
        let mut p = Params::default();
        p.slurry.enabled = false;
        p
    }

    #[test]
    fn run_seed_is_deterministic_and_seeds_differ() {
        let p = dry();
        let a = run_seed(&p, 1, 1.0, 1.0).unwrap();
        let b = run_seed(&p, 1, 1.0, 1.0).unwrap();
        let c = run_seed(&p, 2, 1.0, 1.0).unwrap();
        assert_eq!(a.values, b.values);
        assert_ne!(a.values, c.values);
        assert_eq!(a.names.len(), a.values.len());
        let k = a
            .names
            .iter()
            .position(|n| n == "specific_power_w_per_kg")
            .unwrap();
        assert!(a.values[k] > 0.0, "{:?}", a.values);
        assert!((a.window_s - 0.6).abs() < 0.02, "{}", a.window_s);
    }

    #[test]
    fn exact_window_power_agrees_with_the_smoothed_mean() {
        let r = run_seed(&dry(), 1, 2.0, 2.0).unwrap();
        let get = |n: &str| r.values[r.names.iter().position(|x| x == n).unwrap()];
        let smoothed = get("specific_power_w_per_kg");
        let exact = get("window_specific_power_w_per_kg");
        assert!(
            (smoothed - exact).abs() / exact < 0.2,
            "{smoothed} vs {exact}"
        );
    }

    #[test]
    fn combine_gives_mean_and_standard_error() {
        let mk = |v: f64| SeedResult {
            seed: 0,
            window_s: 1.0,
            coarse_graining_factor: 1.0,
            names: vec!["x".into()],
            values: vec![v],
            power_drift_rel_per_rev: 0.0,
        };
        let c = combine(&[mk(1.0), mk(2.0), mk(3.0)]).unwrap();
        assert!((c.mean[0] - 2.0).abs() < 1e-12);
        // sample sd 1, n = 3 -> se = 1/sqrt(3)
        assert!((c.std_err[0] - 1.0 / 3f64.sqrt()).abs() < 1e-12);
        assert!(combine(&[]).is_err());
    }

    #[test]
    fn drift_is_zero_for_flat_and_detects_a_ramp() {
        assert_eq!(drift_per_rev(&[2.0; 10], 10.0), 0.0);
        // y = 1 + 0.1 * i over 10 frames per revolution: slope 0.1/frame, mean 1.45 -> 1 / 1.45.
        let y: Vec<f64> = (0..10).map(|i| 1.0 + 0.1 * i as f64).collect();
        assert!((drift_per_rev(&y, 10.0) - 1.0 / 1.45).abs() < 1e-9);
    }
}
