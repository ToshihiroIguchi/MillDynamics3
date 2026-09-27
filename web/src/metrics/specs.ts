// Declarative table driving the metrics panel (ui/metricsPanel.ts): one entry per displayed
// value, used to build the panel's rows, pick which values get a sparkline, and generate the CSV
// export's columns (metrics/history.ts). Keeping this as data (not scattered DOM-building code)
// is what lets those three consumers stay in sync without duplicating the metric list three times.

import { toVerticalDegrees } from "./angles";
import type { Metrics } from "./types";

/** Everything a spec's `value()` function might need, beyond `Metrics` itself (things derived
 * from `Params`/`AppState` that mill-core's `Metrics` doesn't carry, e.g. sim time or fps). */
export interface MetricContext {
  metrics: Metrics | null;
  simTime: number;
  fps: number;
  achievedTimeScale: number;
  /** See protocol.ts's `FrameMessage.subStepsPerSecondAchieved`/`subStepsPerSecondRequired`. */
  subStepsPerSecondAchieved: number;
  subStepsPerSecondRequired: number;
  /** Current rotation speed (rpm), or `null` before the first frame/params arrive. */
  rpm: number | null;
  percentCritical: number | null;
  criticalSpeedRpm: number | null;
  /** `params.slurry.enabled`, or `null` before the first params message arrives. Used to hide
   * slurry-only key metrics (e.g. mixing_index) rather than show a permanent "-". */
  slurryEnabled: boolean | null;
}

export type MetricGroupName = "Drum" | "Media" | "Grinding" | "Slurry" | "Solver";

/** %Nc at/above which the charge is considered centrifuging (pinned to the wall, no free-fall/
 * cascading region to measure a toe/shoulder edge from) -- see `unavailable` on `toe_angle_deg`/
 * `shoulder_angle_deg` below. */
const CENTRIFUGING_PERCENT_NC_THRESHOLD = 100;

export interface MetricSpec {
  id: string;
  group: MetricGroupName;
  label: string;
  unit?: string;
  /** Extracts this metric's raw numeric value from the current context, or `null` if
   * unavailable/undefined for the current state (e.g. toe angle while centrifuged). */
  value(ctx: MetricContext): number | null;
  /** Formats a raw value (already extracted by `value()`) for display; receives `null` when
   * `value()` returned `null` and should render the usual "-" placeholder unless overridden. */
  format?(v: number): string;
  /** Include this metric's history in a sparkline (ui/sparkline.ts) in the panel. */
  sparkline?: boolean;
  /** Include this metric as a column in the CSV export (metrics/history.ts). */
  csv?: boolean;
  /** Short reason (e.g. "Centrifuging", "Coarse-grained") this metric is not meaningful right
   * now, or null when it is valid. Receives the value value() already extracted. */
  unavailable?(value: number | null, ctx: MetricContext): string | null;
  /** When unavailable() returns non-null, hide the whole row instead of showing the reason text
   * in the value cell. */
  hideWhenUnavailable?: boolean;
  /** Render this metric in the panel's top "Key results" block (ui/metricsPanel.ts) instead of
   * only inside its group's collapsible `<details>`. Display order within the key block comes
   * from `KEY_METRIC_IDS`, not from this table's order. */
  key?: boolean;
}

/** The key-results block's metric ids, in display order. Kept as an explicit list (rather than
 * `METRIC_SPECS.filter((s) => s.key)`) so the block's order is independent of where each spec
 * happens to sit in the full table -- reordering METRIC_SPECS (e.g. to regroup Media rows) must
 * not silently reorder the key block too. `specs.test.ts` asserts these two views never drift. */
export const KEY_METRIC_IDS = [
  "power_draw",
  "dissipated_power",
  "toe_angle_deg",
  "shoulder_angle_deg",
  "total_kinetic_energy",
  "mixing_index",
] as const;

/** `unavailable` reason for any metric that scales directly with the coarse-graining factor `k`
 * (impact counts as `1/k^2`, per-impact energies as `k^2`, docs/PHYSICS.md §9) and is therefore
 * not comparable to a real mill whenever coarse-graining is active. `k` is exactly `1.0` when
 * inactive, so `> 1` needs no epsilon. */
export function coarseGrainingReason(ctx: MetricContext): string | null {
  const k = ctx.metrics?.coarse_graining_factor;
  return k !== undefined && k > 1 ? "Coarse-grained" : null;
}

function fixed(decimals: number): (v: number) => string {
  return (v) => v.toFixed(decimals);
}

function percentOf1(decimals: number): (v: number) => string {
  return (v) => `${(v * 100).toFixed(decimals)}%`;
}

const m = <K extends keyof Metrics>(key: K) =>
  (ctx: MetricContext): number | null => {
    const value = ctx.metrics?.[key];
    return typeof value === "number" ? value : null;
  };

export const METRIC_SPECS: MetricSpec[] = [
  // --- Drum -------------------------------------------------------------------------------
  { id: "sim_time", group: "Drum", label: "Sim time", unit: "s", value: (ctx) => ctx.simTime, format: fixed(1), csv: true },
  { id: "fps", group: "Drum", label: "FPS", value: (ctx) => ctx.fps, format: fixed(0) },
  { id: "achieved_speed", group: "Drum", label: "Achieved speed", unit: "x", value: (ctx) => ctx.achievedTimeScale, format: fixed(2), csv: true },
  { id: "substeps_achieved_per_s", group: "Drum", label: "Sub-steps/s (achieved)", value: (ctx) => ctx.subStepsPerSecondAchieved, format: fixed(0), csv: true },
  { id: "substeps_required_per_s", group: "Drum", label: "Sub-steps/s (required for 1x)", value: (ctx) => ctx.subStepsPerSecondRequired, format: fixed(0), csv: true },
  { id: "rpm", group: "Drum", label: "Speed", unit: "rpm", value: (ctx) => ctx.rpm, format: fixed(1), csv: true },
  { id: "percent_critical", group: "Drum", label: "Speed", unit: "% Nc", value: (ctx) => ctx.percentCritical, format: fixed(0), csv: true },
  { id: "critical_speed_rpm", group: "Drum", label: "Critical speed", unit: "rpm", value: (ctx) => ctx.criticalSpeedRpm, format: fixed(1) },
  { id: "drum_omega", group: "Drum", label: "Angular velocity", unit: "rad/s", value: m("drum_omega_rad_s"), format: fixed(2) },

  // --- Media (charge geometry) --------------------------------------------------------------
  { id: "true_ball_count", group: "Media", label: "True ball count", value: m("true_ball_count"), format: fixed(0), csv: true },
  { id: "simulated_ball_count", group: "Media", label: "Simulated balls", value: m("simulated_ball_count"), format: fixed(0), csv: true },
  { id: "coarse_graining_factor", group: "Media", label: "Coarse-graining (k)", value: m("coarse_graining_factor"), format: fixed(2), csv: true },
  { id: "effective_ball_diameter_mm", group: "Media", label: "Effective diameter", unit: "mm", value: (ctx) => (ctx.metrics ? ctx.metrics.effective_ball_diameter_m * 1000 : null), format: fixed(2) },
  {
    id: "toe_angle_deg",
    group: "Media",
    label: "Toe angle",
    unit: "° from vertical",
    value: (ctx) =>
      ctx.metrics?.toe_angle_rad != null ? toVerticalDegrees(ctx.metrics.toe_angle_rad, ctx.metrics.drum_omega_rad_s) : null,
    format: fixed(1),
    sparkline: true,
    csv: true,
    key: true,
    // Physically correct (metrics.rs's charge_toe_shoulder, CENTRIFUGE_GAP_THRESHOLD_RAD): once
    // %Nc is at or above 100 the charge is pinned to the wall with no leading/trailing edge to
    // measure, and this is by far the common cause of a null reading at that speed -- so a bare
    // "-" reads as a broken readout rather than an expected state without this explanation.
    unavailable: (value, ctx) =>
      value === null && ctx.percentCritical !== null && ctx.percentCritical >= CENTRIFUGING_PERCENT_NC_THRESHOLD
        ? "Centrifuging"
        : null,
  },
  {
    id: "shoulder_angle_deg",
    group: "Media",
    label: "Shoulder angle",
    unit: "° from vertical",
    value: (ctx) =>
      ctx.metrics?.shoulder_angle_rad != null
        ? toVerticalDegrees(ctx.metrics.shoulder_angle_rad, ctx.metrics.drum_omega_rad_s)
        : null,
    format: fixed(1),
    sparkline: true,
    csv: true,
    key: true,
    unavailable: (value, ctx) =>
      value === null && ctx.percentCritical !== null && ctx.percentCritical >= CENTRIFUGING_PERCENT_NC_THRESHOLD
        ? "Centrifuging"
        : null,
  },
  { id: "centroid_x", group: "Media", label: "Charge centroid X", unit: "m", value: (ctx) => ctx.metrics?.charge_centroid_m?.[0] ?? null, format: fixed(3), csv: true },
  { id: "centroid_y", group: "Media", label: "Charge centroid Y", unit: "m", value: (ctx) => ctx.metrics?.charge_centroid_m?.[1] ?? null, format: fixed(3), csv: true },
  { id: "total_kinetic_energy", group: "Media", label: "Total kinetic energy", unit: "J/m", value: m("total_kinetic_energy_j"), format: fixed(2), sparkline: true, csv: true, key: true },
  { id: "max_ball_overlap", group: "Media", label: "Max ball overlap (of radius)", value: m("max_ball_overlap_fraction"), format: percentOf1(2), csv: true },
  { id: "max_ball_wall_overlap", group: "Media", label: "Max ball-wall overlap (of radius)", value: m("max_ball_wall_overlap_fraction"), format: percentOf1(2), csv: true },

  // --- Grinding (power draw / collisions / dissipation) -------------------------------------
  { id: "power_draw", group: "Grinding", label: "Power draw", unit: "W/m", value: m("power_draw_w"), format: fixed(2), sparkline: true, csv: true, key: true },
  { id: "torque", group: "Grinding", label: "Torque", unit: "N·m/m", value: m("torque_nm"), format: fixed(3), csv: true },
  {
    id: "collision_rate",
    group: "Grinding",
    label: "Collision rate",
    unit: "1/s/m",
    value: m("collision_rate_per_s"),
    format: fixed(1),
    csv: true,
    // Impact count scales as 1/k^2 relative to the true (uncoarsened) population (docs/PHYSICS.md
    // §9), so this reading cannot be compared to a real mill while coarse-graining is active --
    // hide the row entirely rather than show a number that would be misread as physical.
    unavailable: (_value, ctx) => coarseGrainingReason(ctx),
    hideWhenUnavailable: true,
  },
  { id: "dissipated_power", group: "Grinding", label: "Dissipated power", unit: "W/m", value: m("dissipated_power_w"), format: fixed(2), sparkline: true, csv: true, key: true },

  // --- Slurry -------------------------------------------------------------------------------
  { id: "fluid_particle_count", group: "Slurry", label: "Fluid particles", value: m("fluid_particle_count"), format: fixed(0), csv: true },
  {
    id: "pool_angle_range",
    group: "Slurry",
    label: "Pool angular extent",
    unit: "°",
    value: (ctx) =>
      ctx.metrics?.pool_angle_min_rad != null && ctx.metrics.pool_angle_max_rad != null
        ? ((ctx.metrics.pool_angle_max_rad - ctx.metrics.pool_angle_min_rad + Math.PI * 2) % (Math.PI * 2)) * (180 / Math.PI)
        : null,
    format: fixed(0),
  },
  { id: "free_surface_angle", group: "Slurry", label: "Free-surface angle", unit: "°", value: (ctx) => (ctx.metrics?.free_surface_angle_rad != null ? (ctx.metrics.free_surface_angle_rad * 180) / Math.PI : null), format: fixed(1) },
  { id: "free_surface_offset", group: "Slurry", label: "Free-surface offset", unit: "m", value: m("free_surface_offset_m"), format: fixed(3) },
  { id: "pool_depth", group: "Slurry", label: "Pool depth (bottom)", unit: "m", value: m("pool_depth_m"), format: fixed(3), csv: true },
  {
    id: "mixing_index",
    group: "Slurry",
    label: "Mixing index",
    value: m("mixing_index"),
    format: fixed(2),
    sparkline: true,
    csv: true,
    key: true,
    // Promoted to the key-results block (it's the only Slurry-group metric with direct process
    // meaning), so it needs its own reason rather than a bare "-" when there is no slurry to mix.
    unavailable: (_value, ctx) => (ctx.slurryEnabled === false ? "Slurry off" : null),
    hideWhenUnavailable: true,
  },
  { id: "max_compression_error", group: "Slurry", label: "Max compression error", value: m("max_fluid_compression_error_fraction"), format: percentOf1(2), sparkline: true, csv: true },
  { id: "mean_compression_error", group: "Slurry", label: "Mean compression error", value: m("mean_fluid_compression_error_fraction"), format: percentOf1(2), csv: true },
  { id: "max_density_error", group: "Slurry", label: "Density spread, max (incl. free surface)", value: m("max_fluid_density_error_fraction"), format: percentOf1(1), csv: true },
  { id: "mean_density_error", group: "Slurry", label: "Density spread, mean (incl. free surface)", value: m("mean_fluid_density_error_fraction"), format: percentOf1(1), csv: true },
  { id: "mean_shear_rate", group: "Slurry", label: "Mean shear rate", unit: "1/s", value: m("mean_shear_rate_per_s"), format: fixed(2), csv: true },
  { id: "viscosity_iterations", group: "Slurry", label: "Viscosity CG iterations", value: m("viscosity_solver_iterations"), format: fixed(0) },

  // --- Solver health --------------------------------------------------------------------------
  { id: "coupling_clamp_hits", group: "Solver", label: "Coupling clamp hits", value: m("coupling_clamp_hits"), format: fixed(0), sparkline: true, csv: true },
  { id: "substep_displacement", group: "Solver", label: "Substep displacement / diameter", value: m("max_substep_displacement_over_diameter"), format: fixed(3), csv: true },
];

export const METRIC_GROUPS: MetricGroupName[] = ["Drum", "Media", "Grinding", "Slurry", "Solver"];
