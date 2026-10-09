import { describe, expect, it } from "vitest";
import {
  ballCount,
  buildParams,
  buildScaleupConfig,
  combineEnsembles,
  csvCell,
  defaultCondition,
  DEFAULT_STUDY_CONFIG,
  ensembleWarnings,
  extrapolate,
  getMetric,
  intervalsCsv,
  isSensitive,
  metricDefs,
  metricKey,
  metricsCsv,
  percentChange,
  primaryMetricKeys,
  ratioWithError,
  sensitivityChange,
  validateCondition,
  validateStudyConfig,
} from "../../src/scaleup/logic";
import type { EnsembleResult, Interval, PowerLawFit, SensitivityCase, Summary } from "../../src/scaleup/types";

const baseParams = {
  mill: { diameter_m: 0.063, speed_mode: "rpm", speed_value: 100, direction: "counter_clockwise" },
  media: {
    ball_diameter_m: 0.002,
    fill_fraction: 0.3,
    packing_fraction_2d: 0.82,
    density_kg_m3: 6000,
    restitution_ball_ball: 0.7,
    restitution_ball_wall: 0.5,
    friction_ball_ball: 0.25,
    friction_ball_wall: 0.35,
    friction_ball_ball_static: 0.25,
    friction_ball_wall_static: 0.35,
    wet_cushion_viscosity_pa_s: 0,
  },
  slurry: { enabled: true },
  lifters: { count: 0, height_m: 0.02, base_width_m: 0.03, top_width_m: 0.02, phase_deg: 0 },
  simulation: { coarse_graining_mode: "auto", max_balls: 600 },
};

function summary(over: Partial<Summary> = {}): Summary {
  return {
    drum_diameter_m: 0.063,
    ball_diameter_m: 0.002,
    dd_ratio: 31.5,
    n_balls: 244,
    percent_critical: 60,
    rpm: 100,
    ball_mass_3d_kg: 2.5e-5,
    scale_energy_j: 0.5e-3,
    dt_s: 1e-3,
    max_substep_displacement_over_diameter: 0.1,
    reference_stokes_number: null,
    wall_power_per_ball_w: 1,
    dissipated_power_per_ball_w: 1,
    specific_power_w_per_kg: 1,
    impacts_per_ball_per_s: [1, 1, 1, 1],
    impacts_total_per_ball_per_s: 4,
    mean_se_scaled: [1, 1, 1, 1],
    mean_se_scaled_drop: 1,
    mean_se_scaled_shear: 1,
    stress_intensity_scaled: 1,
    thresholds_j: [0.1e-3, 1e-3],
    effective_power_per_ball_w: [0.5, 0.2],
    effective_power_fraction: [0.5, 0.2],
    power_drift_rel_per_rev: 0.01,
    power_rel_std: 0.1,
    ...over,
  };
}

function interval(i: number): Interval {
  return {
    t_s: (i + 1) * 0.1,
    revolutions: (i + 1) / 8,
    wall_power_w: 1 + i,
    dissipated_power_w: 0.5,
    impacts_per_ball_per_s: [1, 2, 3, 4],
    mean_se_scaled: [1, 2, 3, 4],
    effective_power_w: [0.1, 0.2],
    cum_wall_energy_j: i,
    cum_effective_energy_j: [i, 2 * i],
    histogram: [[1], [1], [1], [1]],
  };
}

function ensemble(
  means: Record<string, number>,
  se: Record<string, number> = {},
  sOver: Partial<Summary> = {},
): EnsembleResult {
  const names = Object.keys(means);
  const s = summary(sOver);
  return {
    seeds: [1, 2, 3],
    scalar_names: names,
    mean: names.map((n) => means[n] ?? 0),
    std_err: names.map((n) => se[n] ?? 0),
    first: {
      seed: 1,
      summary: s,
      intervals: [interval(0), interval(1)],
      histogram_total: [[1], [1], [1], [1]],
      se_bin_edges: [1, 10],
      class_names: [],
    },
    summaries: [s],
  };
}

describe("buildParams", () => {
  it("applies every condition field through the schema paths and forces a dry run", () => {
    const c = { ...defaultCondition("big"), drumDiameterMm: 400, ballDiameterMm: 8, percentCritical: 70, lifterCount: 6, viscosityPaS: 2 };
    const p = buildParams(baseParams, c) as typeof baseParams;
    expect(p.mill.diameter_m).toBeCloseTo(0.4);
    expect(p.mill.speed_mode).toBe("percent_critical");
    expect(p.mill.speed_value).toBe(70);
    expect(p.media.ball_diameter_m).toBeCloseTo(0.008);
    expect(p.media.wet_cushion_viscosity_pa_s).toBe(2);
    expect(p.slurry.enabled).toBe(false);
    expect(p.lifters.count).toBe(6);
    expect(p.lifters.height_m).toBeCloseTo(0.008);
    expect(p.simulation.coarse_graining_mode).toBe("off");
    // The base object is not mutated.
    expect(baseParams.slurry.enabled).toBe(true);
    expect(baseParams.mill.diameter_m).toBe(0.063);
  });

  it("throws on a missing field instead of silently ignoring it", () => {
    const broken = { ...baseParams, media: { ...baseParams.media } } as Record<string, unknown>;
    delete (broken.media as Record<string, unknown>).friction_ball_ball_static;
    expect(() => buildParams(broken, defaultCondition("x"))).toThrow();
  });
});

describe("validation", () => {
  it("computes the ball count like Params::true_ball_count", () => {
    expect(ballCount(defaultCondition("lab"), 0.82)).toBeCloseTo(0.3 * 0.82 * 31.5 * 31.5, 3);
  });

  it("refuses more than 50000 balls and fewer than 2", () => {
    const big = { ...defaultCondition("big"), drumDiameterMm: 2000, ballDiameterMm: 2 };
    expect(validateCondition(big).join(" ")).toMatch(/exceeds/);
    const tiny = { ...defaultCondition("tiny"), drumDiameterMm: 4, ballDiameterMm: 2 };
    expect(validateCondition(tiny).join(" ")).toMatch(/below the minimum/);
    expect(validateCondition(defaultCondition("ok"))).toEqual([]);
  });

  it("rejects bad numbers and bad run settings", () => {
    expect(validateCondition({ ...defaultCondition("x"), restitutionBallBall: 1.5 })).not.toEqual([]);
    expect(validateCondition({ ...defaultCondition("x"), ballDiameterMm: Number.NaN })).not.toEqual([]);
    expect(validateStudyConfig(DEFAULT_STUDY_CONFIG)).toEqual([]);
    expect(validateStudyConfig({ ...DEFAULT_STUDY_CONFIG, measureRevolutions: 2.5 })).not.toEqual([]);
  });

  it("builds a wasm config with the requested seeds", () => {
    const c = buildScaleupConfig(DEFAULT_STUDY_CONFIG, [7]);
    expect(c.seeds).toEqual([7]);
    expect(c.measure_revolutions).toBe(6);
    expect(c.thresholds_j).toEqual([]);
  });
});

describe("combineEnsembles", () => {
  it("matches the mean and sample standard error of the Rust ensemble", () => {
    const parts = [1, 2, 3].map((x) => ensemble({ a: x }));
    const e = combineEnsembles(parts);
    expect(e.mean[0]).toBeCloseTo(2);
    // sample std = 1, standard error = 1 / sqrt(3)
    expect(e.std_err[0]).toBeCloseTo(1 / Math.sqrt(3));
    expect(e.summaries).toHaveLength(3);
  });

  it("gives zero error for a single seed", () => {
    expect(combineEnsembles([ensemble({ a: 5 })]).std_err[0]).toBe(0);
  });
});

describe("ratioWithError", () => {
  it("propagates both standard errors", () => {
    const r = ratioWithError({ mean: 6, se: 0.3 }, { mean: 3, se: 0.15 });
    expect(r.ratio).toBeCloseTo(2);
    // relative errors 5 % and 5 % -> sqrt(2) * 5 %
    expect(r.relErr).toBeCloseTo(Math.sqrt(2) * 0.05, 6);
    expect(r.err).toBeCloseTo(2 * Math.sqrt(2) * 0.05, 6);
    expect(r.noisy).toBe(false);
  });

  it("flags a ratio whose error exceeds 10 %", () => {
    expect(ratioWithError({ mean: 1, se: 0.2 }, { mean: 1, se: 0 }).noisy).toBe(true);
  });

  it("returns NaN for a zero reference", () => {
    expect(ratioWithError({ mean: 1, se: 0 }, { mean: 0, se: 0 }).ratio).toBeNaN();
  });
});

describe("metrics", () => {
  it("strips the Joule suffix from threshold names", () => {
    expect(metricKey("effective_power_per_ball_w.thr0(1.000e-4J)")).toBe("effective_power_per_ball_w.thr0");
    expect(metricKey("specific_power_w_per_kg")).toBe("specific_power_w_per_kg");
  });

  it("looks metrics up across conditions with different Joule suffixes", () => {
    const e = ensemble({ "effective_power_per_ball_w.thr0(2.0e-3J)": 4 }, { "effective_power_per_ball_w.thr0(2.0e-3J)": 0.1 });
    expect(getMetric(e, "effective_power_per_ball_w.thr0")).toEqual({ mean: 4, se: 0.1 });
    expect(getMetric(e, "nope")).toBeNull();
  });

  it("lists fixed metrics plus a pair per threshold", () => {
    const defs = metricDefs(summary());
    expect(defs).toHaveLength(5 + 4);
    expect(defs[5]?.label).toContain("m*g*d");
    expect(primaryMetricKeys(["wall_power_per_ball_w", "specific_power_w_per_kg", "effective_power_fraction.thr1(1e-3J)"])).toEqual([
      "specific_power_w_per_kg",
      "effective_power_fraction.thr1",
    ]);
  });
});

describe("steady-state warnings", () => {
  it("warns on drift above 5 %/rev and displacement above 0.5", () => {
    expect(ensembleWarnings("A", ensemble({ a: 1 }))).toEqual([]);
    const w = ensembleWarnings("A", ensemble({ a: 1 }, {}, { power_drift_rel_per_rev: -0.08, max_substep_displacement_over_diameter: 0.7 }));
    expect(w).toHaveLength(2);
  });
});

describe("extrapolate", () => {
  const fit: PowerLawFit = { name: "x", exponent: 1, prefactor: 2, r2: 0.99, expected_exponent: 1 };
  const swept = [20, 30, 45, 60, 90];

  it("predicts prefactor * x^exponent without warnings inside the range", () => {
    const e = extrapolate(fit, 100, swept);
    expect(e.value).toBeCloseTo(200);
    expect(e.warnings).toEqual([]);
  });

  it("warns beyond 1.5x outside the swept range, on either side", () => {
    expect(extrapolate(fit, 136, swept).warnings).toHaveLength(1);
    expect(extrapolate(fit, 135, swept).warnings).toHaveLength(0);
    expect(extrapolate(fit, 10, swept).warnings).toHaveLength(1);
  });

  it("warns on a poor fit", () => {
    expect(extrapolate({ ...fit, r2: 0.9 }, 50, swept).warnings).toHaveLength(1);
  });

  it("rejects a non-positive target", () => {
    expect(extrapolate(fit, 0, swept).value).toBeNaN();
  });
});

describe("sensitivity", () => {
  const base: SensitivityCase = { label: "base", scalar_names: ["a", "b(1J)"], mean: [10, 4] };
  const pert: SensitivityCase = { label: "e +20%", scalar_names: ["a", "b(1J)"], mean: [11.5, 4.2] };

  it("computes the percentage change relative to base and highlights above 10 %", () => {
    expect(percentChange(10, 11.5)).toBeCloseTo(15);
    expect(sensitivityChange(base, pert, "a")).toBeCloseTo(15);
    expect(sensitivityChange(base, pert, "b")).toBeCloseTo(5);
    expect(isSensitive(15)).toBe(true);
    expect(isSensitive(-9.9)).toBe(false);
    expect(isSensitive(-10.1)).toBe(true);
    expect(percentChange(0, 1)).toBeNaN();
  });
});

describe("CSV", () => {
  it("escapes commas, quotes and newlines and blanks non-finite numbers", () => {
    expect(csvCell('a,"b"')).toBe('"a,""b"""');
    expect(csvCell(Number.NaN)).toBe("");
    expect(csvCell(1.5)).toBe("1.5");
  });

  const names = {
    specific_power_w_per_kg: 2,
    impacts_total_per_ball_per_s: 4,
    mean_se_scaled_drop: 1,
    mean_se_scaled_shear: 1,
    stress_intensity_scaled: 1,
  };

  it("builds the metric table with ratio columns for candidates only", () => {
    const ref = ensemble(names);
    const cand = ensemble({ ...names, specific_power_w_per_kg: 4 });
    const lines = metricsCsv([
      { name: "Lab", ens: ref },
      { name: "Big", ens: cand },
    ])
      .trim()
      .split("\n");
    expect(lines[0]).toBe("metric,key,Lab mean,Lab std_err,Big mean,Big std_err,Big ratio_to_reference,Big ratio_err");
    expect(lines[1]).toContain("specific_power_w_per_kg,2,0,4,0,2,0");
    expect(lines.length).toBeGreaterThan(5);
  });

  it("builds the interval time series with a condition column", () => {
    const csv = intervalsCsv([{ name: "Lab", ens: ensemble(names) }]).trim().split("\n");
    expect(csv[0]?.startsWith("condition,t_s,revolutions,wall_power_w,dissipated_power_w,impacts_per_ball_per_s.ball_ball_drop")).toBe(true);
    expect(csv[0]?.endsWith("cum_effective_energy_j.thr1")).toBe(true);
    expect(csv).toHaveLength(3);
    expect(csv[1]?.split(",").length).toBe(csv[0]?.split(",").length);
  });
});
