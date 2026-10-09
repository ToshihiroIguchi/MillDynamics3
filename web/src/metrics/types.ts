// Typed mirror of `mill_core::metrics::Metrics` (crates/mill-core/src/metrics.rs), as serialized
// by `Simulation::metrics_json()` / `serde_json`. Field names are the exact Rust serde names
// (snake_case, no renames on either side) -- see that struct's own doc comment for which fields
// come from a pure function of the current ball/fluid state versus which come from the
// `Simulation`'s own accumulated diagnostics (grinding/solver stats, coarse-graining).

/** Mirrors `mill_core::metrics::ImpactEnergyHistogram`. */
export interface ImpactEnergyHistogram {
  bin_edges_j: number[];
  counts_per_s: number[];
}

/** Names of the four impact classes, in the order of every per-class array in `StressMetrics`. */
export const STRESS_CLASS_NAMES = ["ball_ball_drop", "ball_ball_shear", "ball_wall_drop", "ball_wall_shear"] as const;

/** Mirrors `mill_core::stress::StressMetrics` (crates/mill-core/src/stress.rs): 3D-converted
 * grinding-proxy statistics of the running simulation. Per-class arrays have 4 entries in the
 * order of `STRESS_CLASS_NAMES`; per-threshold arrays are index-aligned with `thresholds_scaled`. */
export interface StressMetrics {
  ball_mass_3d_kg: number;
  scale_energy_j: number;
  impacts_per_ball_per_s: number[];
  impacts_total_per_ball_per_s: number;
  mean_se_scaled: number[];
  mean_se_scaled_drop: number;
  mean_se_scaled_shear: number;
  stress_intensity_scaled: number;
  wall_power_per_ball_w: number;
  specific_power_w_per_kg: number;
  thresholds_scaled: number[];
  effective_power_per_ball_w: number[];
  cumulative_wall_energy_per_ball_j: number;
  cumulative_specific_energy_j_per_kg: number;
  cumulative_effective_energy_per_ball_j: number[];
}

/** Mirrors `mill_core::metrics::Metrics`. */
export interface Metrics {
  // Charge geometry (angle convention: atan2(y, x) in [0, 2*pi), CCW from +x -- see
  // metrics/angles.ts's toVerticalDegrees for the "from vertical" conversion these need for
  // display).
  toe_angle_rad: number | null;
  shoulder_angle_rad: number | null;
  charge_centroid_m: [number, number] | null;

  // Slurry pool geometry.
  pool_angle_min_rad: number | null;
  pool_angle_max_rad: number | null;
  free_surface_angle_rad: number | null;
  free_surface_offset_m: number | null;
  pool_depth_m: number | null;
  mixing_index: number | null;

  // Debug/health checks.
  total_kinetic_energy_j: number;
  max_ball_overlap_fraction: number;
  max_ball_wall_overlap_fraction: number;
  max_fluid_density_error_fraction: number | null;
  mean_fluid_density_error_fraction: number | null;
  max_fluid_compression_error_fraction: number | null;
  mean_fluid_compression_error_fraction: number | null;

  drum_omega_rad_s: number;
  fluid_particle_count: number;

  // Grinding/solver diagnostics (EMA-smoothed by Simulation, see its doc comments).
  power_draw_w: number;
  torque_nm: number;
  collision_rate_per_s: number;
  dissipated_power_w: number;
  impact_energy_histogram: ImpactEnergyHistogram;
  coupling_clamp_hits: number;

  /** 3D grinding-proxy statistics; absent (or null) when the build/state does not provide them. */
  stress?: StressMetrics | null;

  // Coarse-graining / media.
  effective_ball_diameter_m: number;
  simulated_ball_count: number;
  true_ball_count: number;
  coarse_graining_factor: number;

  // Solver health.
  max_substep_displacement_over_diameter: number;
  mean_shear_rate_per_s: number;
  viscosity_solver_iterations: number;
}
