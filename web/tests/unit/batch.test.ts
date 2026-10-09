import { describe, expect, it } from "vitest";
import { buildBatchCsv, buildBatchJson, type BatchRecord } from "../../src/batch/exporters";
import { estimateSeconds, refusalReason, seedList, summarizeParams, validateInputs } from "../../src/batch/plan";
import { metricIdOf, referenceKey, rowInfo } from "../../src/batch/rows";
import {
  buildSensitivityConditions,
  buildSensitivityRows,
  exceedsSensitivity,
  relativeChange,
} from "../../src/batch/sensitivity";
import {
  combineSeeds,
  parseSeedResult,
  ratioWithError,
  steadyState,
  type SeedResult,
} from "../../src/batch/stats";
import { buildTable, referenceFromBatch } from "../../src/batch/table";
import { METRIC_SPECS } from "../../src/metrics/specs";
import { parseReference } from "../../src/metrics/reference";

function result(seed: number, values: number[], drift = 0): SeedResult {
  return {
    seed,
    window_s: 2,
    coarse_graining_factor: 1,
    names: ["specific_power_w_per_kg", "stress_intensity_scaled"].slice(0, values.length),
    values,
    power_drift_rel_per_rev: drift,
  };
}

describe("combineSeeds", () => {
  it("gives mean and sample standard error (n - 1)", () => {
    const c = combineSeeds([result(1, [1]), result(2, [2]), result(3, [6])]);
    expect(c.seeds).toBe(3);
    expect(c.mean[0]).toBeCloseTo(3, 12);
    // sample variance = (4 + 1 + 9) / 2 = 7 ; SE = sqrt(7 / 3)
    expect(c.stdErr[0]).toBeCloseTo(Math.sqrt(7 / 3), 12);
  });

  it("has zero standard error for a single seed", () => {
    const c = combineSeeds([result(1, [5, 7])]);
    expect(c.mean).toEqual([5, 7]);
    expect(c.stdErr).toEqual([0, 0]);
  });

  it("rejects empty input and mismatched names", () => {
    expect(() => combineSeeds([])).toThrow();
    const a = result(1, [1, 2]);
    const b = { ...result(2, [1, 2]), names: ["x", "y"] };
    expect(() => combineSeeds([a, b])).toThrow(/different/);
  });
});

describe("parseSeedResult", () => {
  it("round-trips a Rust-shaped result and rejects malformed ones", () => {
    const r = result(4, [1, 2], 0.01);
    expect(parseSeedResult(JSON.stringify(r))).toEqual(r);
    expect(() => parseSeedResult("{}")).toThrow();
    expect(() => parseSeedResult(JSON.stringify({ ...r, values: [1] }))).toThrow();
  });
});

describe("ratioWithError", () => {
  it("treats the reference as exact when refSe is null", () => {
    const r = ratioWithError(1.2, 0.06, 1, null)!;
    expect(r.ratio).toBeCloseTo(1.2, 12);
    expect(r.relErr).toBeCloseTo(0.05, 12);
  });

  it("propagates both relative errors in quadrature", () => {
    const r = ratioWithError(2, 0.06, 1, 0.04)!; // 3 % and 4 %
    expect(r.relErr).toBeCloseTo(0.05, 12);
  });

  it("returns null for a missing or zero reference", () => {
    expect(ratioWithError(1, 0.1, undefined, null)).toBeNull();
    expect(ratioWithError(1, 0.1, 0, null)).toBeNull();
  });

  it("has unknown error when the current mean is 0 but its SE is not", () => {
    expect(ratioWithError(0, 0.1, 1, null)).toEqual({ ratio: 0, relErr: null });
  });
});

describe("steadyState", () => {
  it("averages drift over seeds and warns above 5 %/rev", () => {
    const ok = steadyState([result(1, [1], 0.01), result(2, [1], -0.02)]);
    expect(ok.mean).toBeCloseTo(-0.005, 12);
    expect(ok.meanAbs).toBeCloseTo(0.015, 12);
    expect(ok.warn).toBe(false);
    expect(steadyState([result(1, [1], 0.08), result(2, [1], 0.06)]).warn).toBe(true);
  });
});

describe("rows", () => {
  it("maps batch names to live spec ids that exist", () => {
    const names = [
      "specific_power_w_per_kg",
      "wall_power_per_ball_w",
      "impacts_total_per_ball_per_s",
      "mean_se_scaled_drop",
      "mean_se_scaled_shear",
      "stress_intensity_scaled",
      "impacts_per_ball_per_s.ball_ball_drop",
      "impacts_per_ball_per_s.ball_ball_shear",
      "impacts_per_ball_per_s.ball_wall_drop",
      "impacts_per_ball_per_s.ball_wall_shear",
      "effective_power_per_ball_w.thr3",
    ];
    for (const n of names) {
      const id = metricIdOf(n);
      expect(id, n).not.toBeNull();
      expect(METRIC_SPECS.some((s) => s.id === id), n).toBe(true);
    }
    expect(metricIdOf("window_specific_power_w_per_kg")).toBeNull();
    expect(referenceKey("window_specific_power_w_per_kg")).toBe("window_specific_power_w_per_kg");
  });

  it("labels threshold rows from the live thresholds", () => {
    expect(rowInfo("effective_power_per_ball_w.thr1", [0.03, 0.1]).label).toBe("Eff. power, SE > 0.1 m·g·d");
    expect(rowInfo("window_effective_power_per_ball_w.thr0", [0.03]).label).toContain("SE > 0.03");
    expect(rowInfo("specific_power_w_per_kg", []).unit).toBe("W/kg");
  });
});

describe("table and reference", () => {
  const combined = combineSeeds([result(1, [10, 1]), result(2, [12, 1.2]), result(3, [11, 1.1])]);

  it("computes ratios against an exact (non-batch) reference", () => {
    const ref = parseReference({ values: { proxy_specific_power: 11 }, settings: {} })!;
    const t = buildTable(combined, [], ref);
    expect(t[0]!.ratio).toBeCloseTo(1, 12);
    expect(t[0]!.ratioRelErr).toBeCloseTo(t[0]!.relSe!, 12);
    expect(t[1]!.ratio).toBeNull();
  });

  it("propagates both errors when the reference came from a batch, and flags > 10 %", () => {
    const ref = referenceFromBatch(combined, { mill: {} }, 1);
    expect(ref.fromBatch).toBe(true);
    expect(ref.seeds).toBe(3);
    const t = buildTable(combined, [], ref);
    expect(t[0]!.ratio).toBeCloseTo(1, 12);
    expect(t[0]!.ratioRelErr).toBeCloseTo(Math.SQRT2 * t[0]!.relSe!, 12);
    const noisy = combineSeeds([result(1, [1, 1]), result(2, [3, 1])]);
    const t2 = buildTable(noisy, [], ref);
    expect(t2[0]!.ratioUncertain).toBe(true);
    expect(t[0]!.ratioUncertain).toBe(false);
  });

  it("round-trips a batch reference through parseReference and stays backward compatible", () => {
    const ref = referenceFromBatch(combined, { simulation: { substeps: 32, dem_iterations: 12 } }, 1);
    const parsed = parseReference(JSON.parse(JSON.stringify(ref)))!;
    expect(parsed.fromBatch).toBe(true);
    expect(parsed.errors?.proxy_specific_power).toBeCloseTo(ref.errors!.proxy_specific_power!, 12);
    expect(parsed.settings.substeps).toBe(32);

    const old = parseReference({ values: { a: 1 }, settings: {}, params: null, simTime: 3, windowed: true })!;
    expect(old.fromBatch).toBeUndefined();
    expect(old.errors).toBeUndefined();
    expect(old.values).toEqual({ a: 1 });
    // Malformed errors are dropped, not fatal.
    expect(parseReference({ values: { a: 1 }, errors: { a: -1, b: 2 }, fromBatch: "yes" })!.errors).toBeUndefined();
  });
});

describe("exporters", () => {
  const results = [result(1, [10, 1]), result(2, [12, 1.2])];
  const rec: BatchRecord = {
    params: { simulation: { seed: 7 }, note: 'has "quote", comma' },
    warmupRevolutions: 4,
    measureRevolutions: 4,
    results,
    combined: combineSeeds(results),
    thresholdsScaled: [0.03],
    reference: null,
    createdAt: "2026-01-01T00:00:00.000Z",
  };

  it("builds a CSV with comment header, per-seed columns and combined stats", () => {
    const csv = buildBatchCsv(rec);
    const lines = csv.trimEnd().split("\n");
    expect(lines[0]).toMatch(/^# /);
    expect(csv).toContain("# seeds: 1 2");
    expect(csv).toContain("# warmup_revolutions: 4");
    const header = lines.find((l) => l.startsWith("name,"))!;
    expect(header).toBe(
      "name,label,unit,seed1,seed2,mean,std_err,rel_se_percent,reference,ratio,ratio_rel_err_percent",
    );
    const row = lines.find((l) => l.startsWith("specific_power_w_per_kg,"))!;
    expect(row.split(",").slice(3, 6)).toEqual(["10", "12", "11"]);
  });

  it("quotes CSV fields containing commas and quotes (labels with units)", () => {
    const csv = buildBatchCsv({ ...rec, results: [{ ...results[0]!, names: ["a,b", "stress_intensity_scaled"] }], combined: combineSeeds([{ ...results[0]!, names: ["a,b", "stress_intensity_scaled"] }]) });
    expect(csv).toContain('"a,b"');
  });

  it("builds JSON with params, seeds, per-seed values and the combined table", () => {
    const j = JSON.parse(buildBatchJson(rec));
    expect(j.seeds).toEqual([1, 2]);
    expect(j.params.simulation.seed).toBe(7);
    expect(j.perSeed).toHaveLength(2);
    expect(j.table[0].mean).toBeCloseTo(11, 12);
    expect(j.combined.seeds).toBe(2);
  });
});

describe("plan", () => {
  it("validates inputs", () => {
    expect(validateInputs("3", "4", "4")).toEqual({ seeds: 3, warmup: 4, measure: 4 });
    expect(typeof validateInputs("0", "4", "4")).toBe("string");
    expect(typeof validateInputs("3", "-1", "4")).toBe("string");
    expect(typeof validateInputs("3", "4", "0")).toBe("string");
    expect(typeof validateInputs("3", "4", "1.5")).toBe("string");
    expect(typeof validateInputs("", "4", "4")).toBe("string");
  });

  it("uses the literal seeds 1..N", () => {
    expect(seedList(3)).toEqual([1, 2, 3]);
  });

  const params = {
    mill: { diameter_m: 0.063, speed_mode: "rpm", speed_value: 100 },
    media: { ball_diameter_m: 0.002 },
    simulation: { substeps: 32, dem_iterations: 12 },
    slurry: { enabled: true },
  };

  it("summarises the applied params", () => {
    expect(summarizeParams(params)).toBe("63 mm drum, 2 mm balls, 100 rpm, sub-steps 32, iterations 12, slurry on");
    expect(summarizeParams(null)).toBe("no parameters yet");
  });

  it("refuses when k > 1, speed 0, or not ready", () => {
    const live = { params, coarseGrainingFactor: 1, thresholdsScaled: [], subStepsPerSecondAchieved: 0 };
    expect(refusalReason(live)).toBeNull();
    expect(refusalReason({ ...live, coarseGrainingFactor: 1.5 })).toMatch(/Coarse-graining/);
    expect(refusalReason({ ...live, coarseGrainingFactor: null })).not.toBeNull();
    expect(refusalReason({ ...live, params: { ...params, mill: { ...params.mill, speed_value: 0 } } })).toMatch(/speed/);
    expect(refusalReason(null)).not.toBeNull();
  });

  it("estimates duration from sub-step rate, or null when unknown", () => {
    // 100 rpm -> 0.6 s/rev -> 36 frames/rev; 8 revolutions x 3 runs x 36 x 32 sub-steps
    const est = estimateSeconds(params, { warmup: 4, measure: 4 }, 3, 1000)!;
    expect(est).toBeCloseTo((3 * 8 * 36 * 32) / 1000, 9);
    expect(estimateSeconds(params, { warmup: 4, measure: 4 }, 3, 0)).toBeNull();
  });
});

describe("sensitivity", () => {
  const params = {
    media: {
      friction_ball_ball: 0.3,
      friction_ball_ball_static: 0.35,
      friction_ball_wall: 0.4,
      friction_ball_wall_static: 0.4,
      restitution_ball_ball: 0.9,
      restitution_ball_wall: 0.5,
    },
    simulation: { seed: 9 },
  };

  it("builds the base plus 8 perturbed conditions without mutating the input", () => {
    const snapshot = JSON.stringify(params);
    const conds = buildSensitivityConditions(params);
    expect(conds).toHaveLength(9);
    expect(conds[0]!.label).toBe("base");
    expect(JSON.stringify(params)).toBe(snapshot);
    expect(JSON.stringify(conds[0]!.params)).toBe(snapshot);
  });

  it("scales one coefficient per condition, clamps restitution, keeps static >= kinetic", () => {
    const conds = buildSensitivityConditions(params);
    const get = (label: string) => (conds.find((c) => c.label === label)!.params.media as Record<string, number>);
    expect(get("friction ball-ball x1.2").friction_ball_ball).toBeCloseTo(0.36, 12);
    expect(get("friction ball-ball x1.2").friction_ball_ball_static).toBeCloseTo(0.42, 12);
    expect(get("friction ball-ball x1.2").friction_ball_wall).toBe(0.4); // untouched
    for (const c of conds) {
      const m = c.params.media as Record<string, number>;
      expect(m.friction_ball_ball_static!).toBeGreaterThanOrEqual(m.friction_ball_ball! - 1e-12);
      expect(m.friction_ball_wall_static!).toBeGreaterThanOrEqual(m.friction_ball_wall! - 1e-12);
      expect(m.restitution_ball_ball!).toBeLessThanOrEqual(1);
      expect(m.restitution_ball_wall!).toBeLessThanOrEqual(1);
    }
    expect(get("restitution ball-ball x1.2").restitution_ball_ball).toBe(1); // 1.08 clamped
    expect(get("restitution ball-ball x0.8").restitution_ball_ball).toBeCloseTo(0.72, 12);
    expect((conds[1]!.params.simulation as Record<string, number>).seed).toBe(9);
  });

  it("raises static friction up to a scaled-up kinetic value when it would fall below", () => {
    const p = { media: { ...params.media, friction_ball_ball: 0.5, friction_ball_ball_static: 0.5 } };
    const up = buildSensitivityConditions(p).find((c) => c.label === "friction ball-ball x1.2")!;
    expect((up.params.media as Record<string, number>).friction_ball_ball_static).toBeCloseTo(0.6, 12);
  });

  it("throws a clear error when a coefficient is missing", () => {
    expect(() => buildSensitivityConditions({ media: {} })).toThrow(/missing/);
  });

  it("computes relative changes and highlights above 10 %", () => {
    expect(relativeChange(10, 11)).toBeCloseTo(0.1, 12);
    expect(relativeChange(0, 1)).toBeNull();
    expect(exceedsSensitivity(0.11)).toBe(true);
    expect(exceedsSensitivity(-0.05)).toBe(false);
    expect(exceedsSensitivity(null)).toBe(false);
    const rows = buildSensitivityRows([result(1, [10, 1]), result(1, [12, 0.9])]);
    expect(rows[0]!.changes[0]).toBeCloseTo(0.2, 12);
    expect(rows[1]!.changes[0]).toBeCloseTo(-0.1, 12);
  });
});
