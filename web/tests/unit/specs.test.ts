import { describe, expect, it } from "vitest";
import { METRIC_GROUPS, METRIC_SPECS, type MetricContext, type MetricGroupName } from "../../src/metrics/specs";

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
