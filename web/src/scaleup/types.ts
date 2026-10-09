// TypeScript mirrors of the JSON produced by crates/mill-core/src/scaleup.rs (serde of
// ScaleupConfig, Interval, Summary, ScaleupResult, EnsembleResult, PowerLawFit, SweepResult and
// SensitivityCase). Option<T> serialises as `T | null`; fixed-size [f64; 4] arrays as number[4].

/** Impact class order used by every per-class array. */
export const CLASS_NAMES = ["ball_ball_drop", "ball_ball_shear", "ball_wall_drop", "ball_wall_shear"] as const;
export type ClassName = (typeof CLASS_NAMES)[number];

/** `ScaleupConfig` (all fields optional on the wasm side; this client always sends all of them). */
export interface ScaleupConfigJson {
  warmup_revolutions: number;
  measure_revolutions: number;
  samples_per_revolution: number;
  dt_ratio: number;
  dem_iterations: number;
  wet_viscosity_pa_s: number;
  thresholds_j: number[];
  seeds: number[];
}

export interface Interval {
  t_s: number;
  revolutions: number;
  wall_power_w: number;
  dissipated_power_w: number;
  impacts_per_ball_per_s: number[];
  mean_se_scaled: number[];
  effective_power_w: number[];
  cum_wall_energy_j: number;
  cum_effective_energy_j: number[];
  histogram: number[][];
}

export interface Summary {
  drum_diameter_m: number;
  ball_diameter_m: number;
  dd_ratio: number;
  n_balls: number;
  percent_critical: number;
  rpm: number;
  ball_mass_3d_kg: number;
  scale_energy_j: number;
  dt_s: number;
  max_substep_displacement_over_diameter: number;
  reference_stokes_number: number | null;
  wall_power_per_ball_w: number;
  dissipated_power_per_ball_w: number;
  specific_power_w_per_kg: number;
  impacts_per_ball_per_s: number[];
  impacts_total_per_ball_per_s: number;
  mean_se_scaled: number[];
  mean_se_scaled_drop: number;
  mean_se_scaled_shear: number;
  stress_intensity_scaled: number;
  thresholds_j: number[];
  effective_power_per_ball_w: number[];
  effective_power_fraction: number[];
  power_drift_rel_per_rev: number;
  power_rel_std: number;
}

export interface ScaleupResult {
  seed: number;
  summary: Summary;
  intervals: Interval[];
  histogram_total: number[][];
  se_bin_edges: number[];
  class_names: string[];
}

export interface EnsembleResult {
  seeds: number[];
  scalar_names: string[];
  mean: number[];
  std_err: number[];
  first: ScaleupResult;
  summaries: Summary[];
}

export interface PowerLawFit {
  name: string;
  exponent: number;
  prefactor: number;
  r2: number;
  expected_exponent: number | null;
}

export interface SweepPoint {
  dd_ratio: number;
  n_balls: number;
  ensemble: EnsembleResult;
}

export interface SweepResult {
  points: SweepPoint[];
  fits: PowerLawFit[];
}

export interface SensitivityCase {
  label: string;
  scalar_names: string[];
  mean: number[];
}
