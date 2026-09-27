import { describe, expect, it } from "vitest";
import { KEY_METRIC_IDS, METRIC_GROUPS, METRIC_SPECS, type MetricContext, type MetricGroupName } from "../../src/metrics/specs";
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
  slurryEnabled: true,
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

describe("KEY_METRIC_IDS / key-flagged specs (ui/metricsPanel.ts's key-results block)", () => {
  it("every KEY_METRIC_IDS entry has a matching spec flagged key: true", () => {
    for (const id of KEY_METRIC_IDS) {
      const spec = METRIC_SPECS.find((s) => s.id === id);
      expect(spec, `no MetricSpec with id "${id}"`).toBeDefined();
      expect(spec?.key, `spec "${id}" is in KEY_METRIC_IDS but not flagged key: true`).toBe(true);
    }
  });

  it("no spec is flagged key: true outside of KEY_METRIC_IDS (the two views never drift)", () => {
    const keyFlaggedIds = METRIC_SPECS.filter((s) => s.key).map((s) => s.id);
    expect(new Set(keyFlaggedIds)).toEqual(new Set(KEY_METRIC_IDS));
  });
});

describe("CSV column order (metrics/history.ts's toCsv, driven by METRIC_SPECS's own order)", () => {
  it("matches the current, intentional column order -- a reorder of METRIC_SPECS that changes this\n" +
    "silently changes every exported CSV's header row, so this guard must be updated deliberately", () => {
    const csvIds = METRIC_SPECS.filter((s) => s.csv).map((s) => s.id);
    expect(csvIds).toEqual([
      "sim_time",
      "achieved_speed",
      "substeps_achieved_per_s",
      "substeps_required_per_s",
      "rpm",
      "percent_critical",
      "true_ball_count",
      "simulated_ball_count",
      "coarse_graining_factor",
      "toe_angle_deg",
      "shoulder_angle_deg",
      "centroid_x",
      "centroid_y",
      "total_kinetic_energy",
      "max_ball_overlap",
      "max_ball_wall_overlap",
      "power_draw",
      "torque",
      "collision_rate",
      "dissipated_power",
      "fluid_particle_count",
      "pool_depth",
      "mixing_index",
      "max_compression_error",
      "mean_compression_error",
      "max_density_error",
      "mean_density_error",
      "mean_shear_rate",
      "coupling_clamp_hits",
      "substep_displacement",
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
    expect(unavailableIds).toEqual(["toe_angle_deg", "shoulder_angle_deg", "collision_rate", "mixing_index"]);

    const hideWhenUnavailableIds = METRIC_SPECS.filter((s) => s.hideWhenUnavailable).map((s) => s.id);
    expect(hideWhenUnavailableIds).toEqual(["collision_rate", "mixing_index"]);
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

  it("mixing_index: 'Slurry off' only when slurryEnabled === false, null when true or unknown (null)", () => {
    const spec = METRIC_SPECS.find((s) => s.id === "mixing_index");
    expect(spec?.unavailable).toBeDefined();

    expect(spec?.unavailable?.(null, { ...baseCtx, slurryEnabled: false })).toBe("Slurry off");
    expect(spec?.unavailable?.(0.62, { ...baseCtx, slurryEnabled: false })).toBe("Slurry off");
    expect(spec?.unavailable?.(0.62, { ...baseCtx, slurryEnabled: true })).toBeNull();
    // Before the first params message arrives, slurryEnabled is null -- don't hide the row over
    // an unknown state, only over a confirmed "off".
    expect(spec?.unavailable?.(null, { ...baseCtx, slurryEnabled: null })).toBeNull();
  });
});
