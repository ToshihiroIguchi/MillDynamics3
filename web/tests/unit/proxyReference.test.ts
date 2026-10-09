import { describe, expect, it } from "vitest";
import { parseMetrics, parseStress } from "../../src/metrics/parse";
import { BINS_PER_DECADE, formatEdgeJ, tickEdgeIndices } from "../../src/metrics/histogramAxis";
import {
  buildReference,
  formatRatio,
  parseReference,
  ratio,
  settingsMismatch,
  settingsOf,
} from "../../src/metrics/reference";

const params = { simulation: { substeps: 4, dem_iterations: 8 } };

describe("ratio / formatRatio", () => {
  it("divides current by reference and rejects bad inputs", () => {
    expect(ratio(3, 2)).toBe(1.5);
    expect(ratio(null, 2)).toBeNull();
    expect(ratio(1, undefined)).toBeNull();
    expect(ratio(1, 0)).toBeNull();
    expect(ratio(Number.NaN, 1)).toBeNull();
  });

  it("formats with 2 significant digits", () => {
    expect(formatRatio(0.8734)).toBe("0.87");
    expect(formatRatio(12.34)).toBe("12");
    expect(formatRatio(1.04)).toBe("1");
    expect(formatRatio(1234)).toBe("1200");
    expect(formatRatio(1e-5)).toBe("1.0e-5");
    expect(formatRatio(null)).toBe("-");
  });
});

describe("settings mismatch", () => {
  it("reports nothing for identical settings", () => {
    const a = settingsOf(params, 1);
    expect(settingsMismatch(a, settingsOf(params, 1))).toEqual([]);
  });

  it("reports each differing numerical setting", () => {
    const a = settingsOf(params, 1);
    const b = settingsOf({ simulation: { substeps: 8, dem_iterations: 12 } }, 2);
    expect(settingsMismatch(a, b)).toEqual(["sub-steps", "ball solver iterations", "coarse-graining factor"]);
  });

  it("ignores unknown settings", () => {
    expect(settingsMismatch(settingsOf(null, null), settingsOf(params, 1))).toEqual([]);
  });
});

describe("buildReference / parseReference", () => {
  it("prefers window means over current values and omits ids with neither", () => {
    const ref = buildReference(["a", "b", "c"], { a: 1, b: 2, c: null }, (id) => (id === "a" ? 10 : null), params, 1, 5);
    expect(ref.values).toEqual({ a: 10, b: 2 });
    expect(ref.windowed).toBe(true);
    expect(ref.settings.substeps).toBe(4);
    expect(ref.simTime).toBe(5);
  });

  it("marks a reference built only from current values as not windowed", () => {
    const ref = buildReference(["a"], { a: 1 }, () => null, params, 1, 0);
    expect(ref.windowed).toBe(false);
  });

  it("round-trips through JSON and rejects malformed input", () => {
    const ref = buildReference(["a"], { a: 1.5 }, () => null, params, 1, 3);
    expect(parseReference(JSON.parse(JSON.stringify(ref)))).toEqual(ref);
    expect(parseReference(null)).toBeNull();
    expect(parseReference({ values: "x" })).toBeNull();
    expect(parseReference({ values: {} })).toBeNull();
    expect(parseReference({ values: { a: "nope" } })).toBeNull();
  });
});

describe("tolerant stress parsing", () => {
  const stress = {
    ball_mass_3d_kg: 1,
    scale_energy_j: 1,
    impacts_per_ball_per_s: [1, 2, 3, 4],
    impacts_total_per_ball_per_s: 10,
    mean_se_scaled: [1, 2, 3, 4],
    mean_se_scaled_drop: 1,
    mean_se_scaled_shear: 1,
    stress_intensity_scaled: 1,
    wall_power_per_ball_w: 1,
    specific_power_w_per_kg: 1,
    thresholds_scaled: [0.1, 1],
    effective_power_per_ball_w: [1, 2],
    cumulative_wall_energy_per_ball_j: 1,
    cumulative_specific_energy_j_per_kg: 1,
    cumulative_effective_energy_per_ball_j: [1, 2],
  };

  it("accepts a well-formed block", () => {
    expect(parseStress(stress)).toBeDefined();
  });

  it("returns undefined for missing or malformed blocks", () => {
    expect(parseStress(undefined)).toBeUndefined();
    expect(parseStress(null)).toBeUndefined();
    expect(parseStress({ ...stress, impacts_per_ball_per_s: [1, 2] })).toBeUndefined();
    expect(parseStress({ ...stress, effective_power_per_ball_w: [1] })).toBeUndefined();
    expect(parseStress({ ...stress, specific_power_w_per_kg: null })).toBeUndefined();
  });

  it("parseMetrics keeps the rest of Metrics when stress is missing", () => {
    const m = parseMetrics(JSON.stringify({ power_draw_w: 3, coarse_graining_factor: 1 }));
    expect(m.power_draw_w).toBe(3);
    expect(m.stress).toBeUndefined();
    expect(parseMetrics(JSON.stringify({ stress })).stress).toBeDefined();
  });
});

describe("histogram axis", () => {
  it("labels one edge per decade of the 29-edge stress histogram", () => {
    const edges = Array.from({ length: 29 }, (_, i) => 1e-9 * 10 ** (i / BINS_PER_DECADE));
    expect(tickEdgeIndices(edges)).toEqual([0, 4, 8, 12, 16, 20, 24, 28]);
    expect(formatEdgeJ(edges[0]!)).toBe("1.0e-9");
  });

  it("skips non-positive edges", () => {
    expect(tickEdgeIndices([0, 1, 2, 3, 4])).toEqual([4]);
  });
});
