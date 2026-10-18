//! Grinding-proxy statistics of the running simulation: the impacts the ball solver records
//! ([`crate::dem::ImpactRecord`]) converted to **3D** quantities.
//!
//! Grinding itself is not simulated. What is reported are the quantities mechanistic grinding
//! models are driven by: the distribution of **stress energies** (`SE`, the relative kinetic energy
//! of an impact), their frequency per ball, and the power drawn. `SE` uses the mass of the
//! *3D sphere* a disc stands for (`m = rho * 4/3 * pi * r^3`), not the 2D unit-depth disc mass
//! ([`crate::dem::ball_mass`]): the disc mass overstates `SE` by a ball-size-dependent factor
//! (about 750x for a 2 mm ball), which would corrupt comparisons between ball sizes.
//!
//! Everything is also given dimensionless (`SE / (m g d)`, `d` the ball diameter) so conditions of
//! different size can be compared. Impacts are split into **drop** type (normal approach speed >=
//! tangential slip speed) and **shear** type, for ball-ball and ball-wall separately
//! ([`CLASS_NAMES`]); the two scale differently with drum and ball size.
//!
//! Conventions: "per ball" powers use a one-ball-thick slab (depth = ball diameter), i.e. the 2D
//! per-metre-depth power times `d` divided by the ball count. This is equivalent to an implied 3D
//! solid fraction `phi_3D = (2/3) * packing_fraction_2d` (0.547 at the default 0.82); ratios between
//! conditions do not depend on this depth basis. Results are meant for **relative**
//! comparison between conditions run with the same code and the same numerical settings.
//! With coarse-graining (`k > 1`) the simulated ball is not the real ball: do not use them.

use crate::dem::{ImpactKind, ImpactRecord};
use serde::Serialize;
use std::collections::BTreeMap;

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

/// Stress-energy thresholds (multiples of `m g d`) above which an impact counts as "effective" in
/// [`StressMetrics::effective_power_per_ball_w`].
pub const THRESHOLDS_REL: [f64; 6] = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0];

/// Re-contacts of the same pair (or of a ball with the wall) closer in time than
/// `DEBOUNCE_FACTOR * sqrt(d / g)` count as one impact episode. Without it the impact count keeps
/// rising as the sub-step shrinks (finer steps resolve more micro-bounces of a resting contact).
/// The window is ~10x the sub-step of the validated accuracy settings, so it is sub-step independent.
pub const DEBOUNCE_FACTOR: f64 = 0.25;

/// Edges of the dimensionless stress-energy bins (`SE / (m g d)`), `SE_BINS + 1` values.
pub fn se_bin_edges() -> Vec<f64> {
    (0..=SE_BINS)
        .map(|i| 10f64.powf(SE_LOG10_MIN + i as f64 / SE_BINS_PER_DECADE))
        .collect()
}

/// Bin index of a dimensionless stress energy (values outside the range go to the end bins).
pub fn se_bin(e_hat: f64) -> usize {
    if e_hat.is_nan() || e_hat <= 0.0 {
        return 0;
    }
    let t = (e_hat.log10() - SE_LOG10_MIN) * SE_BINS_PER_DECADE;
    (t.max(0.0) as usize).min(SE_BINS - 1)
}

/// Class index (see [`CLASS_NAMES`]) and stress energy (J) of one impact for a 3D ball mass `m3`
/// (kg). Ball-ball uses the reduced mass `m3 / 2`, ball-wall `m3`.
pub fn classify(rec: &ImpactRecord, m3: f64) -> (usize, f64) {
    let bb = rec.kind == ImpactKind::BallBall;
    let (vn, vs) = (rec.normal_speed as f64, rec.slip_speed as f64);
    let class = match (bb, vs > vn) {
        (true, false) => 0,
        (true, true) => 1,
        (false, false) => 2,
        (false, true) => 3,
    };
    let m_red = if bb { 0.5 * m3 } else { m3 };
    (class, 0.5 * m_red * (vn * vn + vs * vs))
}

/// EMA-smoothed accumulators (rates per second for the whole simulated ball population).
#[derive(Debug, Clone, Default)]
pub struct StressAccum {
    count_rate: [f64; 4],
    /// Rate of summed dimensionless stress energy per class.
    se_hat_rate: [f64; 4],
    hist_rate: [f64; SE_BINS],
    /// Rate (J/s) of stress energy carried by impacts above each [`THRESHOLDS_REL`] entry.
    above_rate_j: [f64; THRESHOLDS_REL.len()],
    /// Cumulative (not smoothed) wall work since construction, per metre of mill depth (J/m).
    cum_wall_j: f64,
    /// Cumulative stress energy (J, whole population) carried by impacts above each threshold.
    cum_above_j: [f64; THRESHOLDS_REL.len()],
    /// Time of the last counted impact per ball pair / ball-wall (see [`DEBOUNCE_FACTOR`]).
    last_counted: BTreeMap<(u32, u32), f64>,
    /// Debounce window as a multiple of `sqrt(d / g)`; [`DEBOUNCE_FACTOR`] unless overridden.
    debounce_factor: Option<f64>,
    updates: u64,
}

impl StressAccum {
    /// Overrides the debounce window factor (`0` disables it). Diagnostic use.
    pub fn set_debounce_factor(&mut self, factor: f64) {
        self.debounce_factor = Some(factor);
    }

    /// Folds one sub-step's impacts into the EMAs. `m3` is the 3D ball mass (kg), `scale_j = m3 g d`,
    /// `alpha` the EMA weight of this sub-step of length `sub_dt`.
    #[allow(clippy::too_many_arguments)]
    pub fn update(
        &mut self,
        records: &[ImpactRecord],
        m3: f64,
        scale_j: f64,
        d: f64,
        now_s: f64,
        sub_dt: f64,
        alpha: f64,
        wall_work_j: f64,
    ) {
        let window = self.debounce_factor.unwrap_or(DEBOUNCE_FACTOR) * (d / 9.81).sqrt();
        self.updates += 1;
        if self.updates.is_multiple_of(4096) {
            self.last_counted.retain(|_, t| now_s - *t < 4.0 * window);
        }
        let mut count = [0.0f64; 4];
        let mut se_hat = [0.0f64; 4];
        let mut hist = [0.0f64; SE_BINS];
        let mut above = [0.0f64; THRESHOLDS_REL.len()];
        for rec in records {
            if window > 0.0 {
                let key = (rec.ball.min(rec.partner), rec.ball.max(rec.partner));
                let key = if rec.partner == ImpactRecord::WALL {
                    (rec.ball, ImpactRecord::WALL)
                } else {
                    key
                };
                match self.last_counted.get(&key) {
                    Some(&t0) if now_s - t0 < window => continue,
                    _ => {
                        self.last_counted.insert(key, now_s);
                    }
                }
            }
            let (class, se) = classify(rec, m3);
            let e_hat = se / scale_j;
            count[class] += 1.0;
            se_hat[class] += e_hat;
            hist[se_bin(e_hat)] += 1.0;
            for (a, thr) in above.iter_mut().zip(THRESHOLDS_REL) {
                if e_hat >= thr {
                    *a += se;
                }
            }
        }
        self.cum_wall_j += wall_work_j;
        for (c, a) in self.cum_above_j.iter_mut().zip(above) {
            *c += a;
        }
        let ema = |state: &mut f64, instant: f64| *state += alpha * (instant / sub_dt - *state);
        for c in 0..4 {
            ema(&mut self.count_rate[c], count[c]);
            ema(&mut self.se_hat_rate[c], se_hat[c]);
        }
        for (s, h) in self.hist_rate.iter_mut().zip(hist) {
            ema(s, h);
        }
        for (s, a) in self.above_rate_j.iter_mut().zip(above) {
            ema(s, a);
        }
    }

    /// Total impact rate (impacts/s, whole simulated population).
    pub fn total_rate(&self) -> f64 {
        self.count_rate.iter().sum()
    }

    /// The smoothed pooled stress-energy histogram: rates (impacts/s) per bin.
    pub fn histogram_rate(&self) -> &[f64; SE_BINS] {
        &self.hist_rate
    }

    /// Builds the reported metrics. `n_balls` is the simulated ball count, `m3` the 3D ball mass
    /// (kg), `scale_j = m3 g d`, `d` the ball diameter (m), `power_draw_w` the EMA wall power per
    /// metre of mill depth.
    pub fn metrics(
        &self,
        n_balls: usize,
        m3: f64,
        scale_j: f64,
        d: f64,
        power_draw_w: f64,
    ) -> StressMetrics {
        let n = n_balls.max(1) as f64;
        let mut impacts = [0.0f32; 4];
        let mut mean_se = [0.0f32; 4];
        for c in 0..4 {
            impacts[c] = (self.count_rate[c] / n) as f32;
            if self.count_rate[c] > 0.0 {
                mean_se[c] = (self.se_hat_rate[c] / self.count_rate[c]) as f32;
            }
        }
        let pooled = |a: usize, b: usize| {
            let w = self.count_rate[a] + self.count_rate[b];
            if w > 0.0 {
                ((self.se_hat_rate[a] + self.se_hat_rate[b]) / w) as f32
            } else {
                0.0
            }
        };
        let total = self.total_rate();
        let power_per_ball = power_draw_w * d / n;
        StressMetrics {
            ball_mass_3d_kg: m3 as f32,
            scale_energy_j: scale_j as f32,
            impacts_per_ball_per_s: impacts,
            impacts_total_per_ball_per_s: (total / n) as f32,
            mean_se_scaled: mean_se,
            mean_se_scaled_drop: pooled(0, 2),
            mean_se_scaled_shear: pooled(1, 3),
            stress_intensity_scaled: if total > 0.0 {
                (self.se_hat_rate.iter().sum::<f64>() / total) as f32
            } else {
                0.0
            },
            wall_power_per_ball_w: power_per_ball as f32,
            specific_power_w_per_kg: (power_per_ball / m3) as f32,
            thresholds_scaled: THRESHOLDS_REL.iter().map(|&t| t as f32).collect(),
            effective_power_per_ball_w: self.above_rate_j.iter().map(|&r| (r / n) as f32).collect(),
            cumulative_wall_energy_per_ball_j: self.cum_wall_j * d / n,
            cumulative_specific_energy_j_per_kg: self.cum_wall_j * d / n / m3,
            cumulative_effective_energy_per_ball_j: self
                .cum_above_j
                .iter()
                .map(|&e| e / n)
                .collect(),
        }
    }
}

/// 3D-converted grinding-proxy readings (see the module doc comment for conventions).
#[derive(Debug, Clone, Default, Serialize)]
pub struct StressMetrics {
    /// Mass of the 3D sphere each simulated ball stands for (kg).
    pub ball_mass_3d_kg: f32,
    /// `m g d` (J): the energy unit of every `*_scaled` value.
    pub scale_energy_j: f32,
    /// Impacts per ball per second, per class ([`CLASS_NAMES`]).
    pub impacts_per_ball_per_s: [f32; 4],
    pub impacts_total_per_ball_per_s: f32,
    /// Mean `SE / (m g d)` per class (0 when the class has no impacts).
    pub mean_se_scaled: [f32; 4],
    pub mean_se_scaled_drop: f32,
    pub mean_se_scaled_shear: f32,
    /// Mean `SE / (m g d)` over all impacts (Kwade-type stress intensity).
    pub stress_intensity_scaled: f32,
    /// Wall power (including slurry drag on the wall) per ball, one-ball-thick slab (W).
    pub wall_power_per_ball_w: f32,
    /// `wall_power_per_ball_w / ball_mass_3d_kg` (W/kg). Scales as `sqrt(L)` between similar mills.
    pub specific_power_w_per_kg: f32,
    /// Thresholds, in multiples of `m g d`, used for `effective_power_per_ball_w`.
    pub thresholds_scaled: Vec<f32>,
    /// Rate of stress energy carried by impacts above each threshold (W per ball).
    pub effective_power_per_ball_w: Vec<f32>,
    /// Wall work per ball since the simulation was created (J, one-ball-thick slab). Monotone:
    /// take differences between two times for a window. `f64`: it grows without bound.
    pub cumulative_wall_energy_per_ball_j: f64,
    /// `cumulative_wall_energy_per_ball_j / ball_mass_3d_kg` (J/kg).
    pub cumulative_specific_energy_j_per_kg: f64,
    /// Cumulative stress energy per ball carried by impacts above each threshold (J). The time two
    /// conditions need to reach the same value is the relative-grinding-time proxy.
    pub cumulative_effective_energy_per_ball_j: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(kind: ImpactKind, vn: f32, vs: f32) -> ImpactRecord {
        ImpactRecord {
            kind,
            ball: 0,
            partner: 1,
            normal_speed: vn,
            slip_speed: vs,
            restitution: 0.7,
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
    fn stress_energy_uses_the_reduced_3d_mass() {
        let m3 = 2.0e-5;
        // Ball-ball head-on at 2 m/s: 1/2 * (m/2) * 4.
        let (c, se) = classify(&rec(ImpactKind::BallBall, 2.0, 0.0), m3);
        assert_eq!(c, 0);
        assert!((se - 0.5 * 0.5 * m3 * 4.0).abs() < 1e-12);
        // Ball-wall, slip-dominated: shear class, full mass, both components.
        let (c, se) = classify(&rec(ImpactKind::BallWall, 1.0, 3.0), m3);
        assert_eq!(c, 3);
        assert!((se - 0.5 * m3 * 10.0).abs() < 1e-12);
    }

    #[test]
    fn ema_converges_to_the_per_second_rate_and_normalises_per_ball() {
        let (m3, d) = (2.0e-5, 0.002);
        let scale = m3 * 9.81 * d;
        let mut acc = StressAccum::default();
        acc.set_debounce_factor(0.0);
        let sub_dt = 1e-3;
        // Two identical ball-ball impacts per sub-step -> 2000 impacts/s for the population.
        let r = [rec(ImpactKind::BallBall, 0.5, 0.0); 2];
        for _ in 0..20_000 {
            acc.update(
                &r,
                m3,
                scale,
                d,
                0.0,
                sub_dt,
                1.0 - (-sub_dt / 1.0f64).exp(),
                0.0,
            );
        }
        let m = acc.metrics(100, m3, scale, d, 0.0);
        assert!((m.impacts_per_ball_per_s[0] - 20.0).abs() < 0.05, "{m:?}");
        let e_hat = 0.5 * 0.5 * m3 * 0.25 / scale;
        assert!((m.mean_se_scaled[0] as f64 - e_hat).abs() / e_hat < 1e-4);
        assert_eq!(m.mean_se_scaled[1], 0.0);
        // Cumulative energy above the lowest threshold: all 2 * 20000 impacts carry the same SE.
        let n_imp = 2.0 * 20_000.0;
        let se_each = 0.5 * 0.5 * m3 * 0.25;
        assert!(
            ((m.cumulative_effective_energy_per_ball_j[0] * 100.0) - n_imp * se_each).abs()
                / (n_imp * se_each)
                < 1e-9,
            "{m:?}"
        );
        let hist_total: f64 = acc.histogram_rate().iter().sum();
        assert!((hist_total - 2000.0).abs() < 5.0);
    }

    #[test]
    fn debounce_counts_one_impact_per_episode() {
        let (m3, d) = (2.0e-5, 0.002);
        let scale = m3 * 9.81 * d;
        let window = DEBOUNCE_FACTOR * (d / 9.81f64).sqrt(); // ~3.6 ms
        let r = [rec(ImpactKind::BallBall, 0.5, 0.0)];
        let mut acc = StressAccum::default();
        let sub_dt = 1e-4;
        // Re-contacts every 0.1 ms for 2 windows -> 2 counted episodes, not 72.
        let mut t = 0.0;
        while t < 2.0 * window - sub_dt {
            acc.update(&r, m3, scale, d, t, sub_dt, 1.0, 0.0);
            t += sub_dt;
        }
        let counted: f64 = acc.histogram_rate().iter().sum::<f64>() * sub_dt;
        // With alpha = 1 the EMA holds only the last sub-step; use the cumulative energy instead.
        assert!(counted <= 1.0 + 1e-9);
        let m = acc.metrics(1, m3, scale, d, 0.0);
        let se_each = 0.5 * 0.5 * m3 * 0.25;
        let n_counted = m.cumulative_effective_energy_per_ball_j[0] / se_each;
        assert!((n_counted - 2.0).abs() < 1e-6, "{n_counted}");
    }

    #[test]
    fn specific_power_is_per_ball_slab_over_3d_mass() {
        let (m3, d) = (2.0e-5, 0.002);
        let m = StressAccum::default().metrics(200, m3, m3 * 9.81 * d, d, 4.0);
        // 4 W/m * 0.002 m / 200 balls = 4e-5 W per ball; / 2e-5 kg = 2 W/kg.
        assert!((m.specific_power_w_per_kg - 2.0).abs() < 1e-5, "{m:?}");
    }
}
