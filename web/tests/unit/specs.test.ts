import { describe, expect, it } from "vitest";
import { METRIC_GROUPS, METRIC_SPECS, type MetricContext, type MetricGroupName } from "../../src/metrics/specs";
import type { Metrics } from "../../src/metrics/types";

const CONTEXT_ONLY_IDS = new Set([
  "sim_time",
  "fps",
  "achieved_speed",
  "substeps_achieved_per_s",
  "substeps_required_per_s",
  "rpm",
  "percent_critical",
  "critical_speed_rpm",
]);

const CONTEXT_ONLY_EXPECTED: Record<string, number> = {
  sim_time: 1.5,
  fps: 60,
  achieved_speed: 1,
  substeps_achieved_per_s: 480,
  substeps_required_per_s: 480,
  rpm: 30,
  percent_critical: 71,
  critical_speed_rpm: 42.3,
};

const baseCtx: MetricContext = {
  metrics: null,
  simTime: 1.5,
  fps: 60,
  achievedTimeScale: 1,
  subStepsPerSecondAchieved: 480,
  subStepsPerSecondRequired: 480,
  rpm: 30,
  percentCritical: 71,
  criticalSpeedRpm: 42.3,
};

describe("METRIC_SPECS with a metrics-less context", () => {
  it("returns the raw context value for context-only specs, and null for every metrics-derived spec", () => {
    for (const spec of METRIC_SPECS) {
      const v = spec.value(baseCtx);
      if (CONTEXT_ONLY_IDS.has(spec.id)) {
        expect(v).toBe(CONTEXT_ONLY_EXPECTED[spec.id]);
      } else {
        expect(v).toBeNull();
      }
    }
  });
});

describe("METRIC_GROUPS", () => {
  it("is exactly Drum, Media, Grinding, Slurry, Solver in that order", () => {
    expect(METRIC_GROUPS).toEqual(["Drum", "Media", "Grinding", "Slurry", "Solver"]);
  });

  it("covers every spec's group", () => {
    const groups: MetricGroupName[] = METRIC_GROUPS;
    for (const spec of METRIC_SPECS) {
      expect(groups).toContain(spec.group);
    }
  });
});

describe("format", () => {
  it("formats effective_ball_diameter_mm with fixed(2)", () => {
    const spec = METRIC_SPECS.find((s) => s.id === "effective_ball_diameter_mm");
    expect(spec).toBeDefined();
    expect(spec?.format).toBeDefined();
    expect(spec?.format?.(11.09)).toContain("11.09");
  });

  it("formats max_ball_overlap as a percentage", () => {
    const spec = METRIC_SPECS.find((s) => s.id === "max_ball_overlap");
    expect(spec).toBeDefined();
    expect(spec?.format?.(0.9386)).toBe("93.86%");
  });
});

describe("sparkline specs", () => {
  it("matches the exact set of ids flagged sparkline: true in specs.ts", () => {
    const sparklineIds = METRIC_SPECS.filter((s) => s.sparkline).map((s) => s.id);
    expect(sparklineIds).toEqual([
      "toe_angle_deg",
      "shoulder_angle_deg",
      "total_kinetic_energy",
      "power_draw",
      "dissipated_power",
      "mixing_index",
      "max_compression_error",
      "coupling_clamp_hits",
    ]);
  });
});

// Minimal but complete Metrics fixture (every field of the interface needs a value) for
// `unavailable` tests below that need a non-null `ctx.metrics`.
const fakeMetrics: Metrics = {
  toe_angle_rad: null,
  shoulder_angle_rad: null,
  charge_centroid_m: null,
  pool_angle_min_rad: null,
  pool_angle_max_rad: null,
  free_surface_angle_rad: null,
  free_surface_offset_m: null,
  pool_depth_m: null,
  mixing_index: null,
  total_kinetic_energy_j: 0,
  max_ball_overlap_fraction: 0,
  max_ball_wall_overlap_fraction: 0,
  max_fluid_density_error_fraction: null,
  mean_fluid_density_error_fraction: null,
  max_fluid_compression_error_fraction: null,
  mean_fluid_compression_error_fraction: null,
  drum_omega_rad_s: 0,
  fluid_particle_count: 0,
  power_draw_w: 0,
  torque_nm: 0,
  collision_rate_per_s: 0,
  dissipated_power_w: 0,
  impact_energy_histogram: { bin_edges_j: [], counts_per_s: [] },
  coupling_clamp_hits: 0,
  effective_ball_diameter_m: 0,
  simulated_ball_count: 0,
  true_ball_count: 0,
  coarse_graining_factor: 1,
  max_substep_displacement_over_diameter: 0,
  mean_shear_rate_per_s: 0,
  viscosity_solver_iterations: 0,
};

describe("unavailable / hideWhenUnavailable specs", () => {
  it("matches the exact set of ids flagged unavailable/hideWhenUnavailable in specs.ts", () => {
    const unavailableIds = METRIC_SPECS.filter((s) => s.unavailable).map((s) => s.id);
    expect(unavailableIds).toEqual(["toe_angle_deg", "shoulder_angle_deg", "collision_rate"]);

    const hideWhenUnavailableIds = METRIC_SPECS.filter((s) => s.hideWhenUnavailable).map((s) => s.id);
    expect(hideWhenUnavailableIds).toEqual(["collision_rate"]);
  });

  it("collision_rate: 'Coarse-grained' when k > 1, null at k = 1, null with no metrics", () => {
    const spec = METRIC_SPECS.find((s) => s.id === "collision_rate");
    expect(spec?.unavailable).toBeDefined();

    const ctxK2: MetricContext = { ...baseCtx, metrics: { ...fakeMetrics, coarse_graining_factor: 2 } };
    expect(spec?.unavailable?.(null, ctxK2)).toBe("Coarse-grained");

    const ctxK1: MetricContext = { ...baseCtx, metrics: { ...fakeMetrics, coarse_graining_factor: 1 } };
    expect(spec?.unavailable?.(null, ctxK1)).toBeNull();

    expect(spec?.unavailable?.(null, baseCtx)).toBeNull();
  });

  it("toe_angle_deg/shoulder_angle_deg: 'Centrifuging' only when the value is null and percentCritical >= 100", () => {
    for (const id of ["toe_angle_deg", "shoulder_angle_deg"]) {
      const spec = METRIC_SPECS.find((s) => s.id === id);
      expect(spec?.unavailable).toBeDefined();

      expect(spec?.unavailable?.(null, { ...baseCtx, percentCritical: 100 })).toBe("Centrifuging");
      expect(spec?.unavailable?.(null, { ...baseCtx, percentCritical: 150 })).toBe("Centrifuging");
      expect(spec?.unavailable?.(null, { ...baseCtx, percentCritical: 99 })).toBeNull();
      expect(spec?.unavailable?.(null, { ...baseCtx, percentCritical: null })).toBeNull();

      // Regression guard: a genuinely computed angle at >= 100% Nc must never be overwritten.
      expect(spec?.unavailable?.(12.3, { ...baseCtx, percentCritical: 120 })).toBeNull();
    }
  });
});
