import { describe, expect, it } from "vitest";
import { buildBatchCsv, buildBatchJson, type BatchRecord } from "../../src/batch/exporters";
import { combineSeeds, type SeedResult } from "../../src/batch/stats";
import { MetricsHistory } from "../../src/metrics/history";
import { METRIC_SPECS, THRESHOLDS_REL } from "../../src/metrics/specs";
import { modelBasisLines, numericsLines } from "../../src/report/assumptions";

const DEFAULT_LIKE_PARAMS = {
  media: { packing_fraction_2d: 0.82 },
  simulation: { substeps: 8, dem_iterations: 2, seed: 1 },
};

describe("modelBasisLines / numericsLines", () => {
  it("states phi_3D = 0.547 at the default packing fraction", () => {
    const text = modelBasisLines(DEFAULT_LIKE_PARAMS).join("\n");
    expect(text).toContain("0.547");
    expect(text).toContain("phi_3D");
  });

  it("falls back to defaults for missing params", () => {
    expect(modelBasisLines(null).join("\n")).toContain("0.547");
    const n = numericsLines(null).join("\n");
    expect(n).toContain("Sub-steps per frame: 8");
    expect(n).toContain("Random seed: 1");
  });
});

describe("threshold metric labels", () => {
  it("are distinct for the effective-power and cumulative-energy rows", () => {
    const power = METRIC_SPECS.filter((s) => s.id.startsWith("proxy_eff_power_")).map((s) => s.label);
    const cum = METRIC_SPECS.filter((s) => s.id.startsWith("proxy_cum_eff_energy_")).map((s) => s.label);
    expect(power).toHaveLength(THRESHOLDS_REL.length);
    expect(new Set(power).size).toBe(6);
    expect(new Set(cum).size).toBe(6);
  });

  it("every csv metric label is unique", () => {
    const labels = METRIC_SPECS.filter((s) => s.csv).map((s) => (s.unit ? `${s.label} (${s.unit})` : s.label));
    expect(new Set(labels).size).toBe(labels.length);
  });
});

describe("toCsv header comments", () => {
  it("writes # lines before the header and keeps the body unchanged", () => {
    const h = new MetricsHistory(10, 0);
    h.push(0, { a: 1 });
    const csv = h.toCsv([{ id: "a", label: "A" }], ["first", "multi\nline"]);
    expect(csv.split("\n")).toEqual(["# first", "# multi", "# line", "sim_time_s,A", "0.000,1"]);
    expect(h.toCsv([{ id: "a", label: "A" }]).startsWith("sim_time_s,A")).toBe(true);
  });
});

describe("batch exports carry the basis", () => {
  const results: SeedResult[] = [1, 2].map((seed) => ({
    seed,
    window_s: 2,
    coarse_graining_factor: 1,
    names: ["specific_power_w_per_kg"],
    values: [10 + seed],
    power_drift_rel_per_rev: 0,
  }));
  const rec: BatchRecord = {
    params: DEFAULT_LIKE_PARAMS,
    warmupRevolutions: 4,
    measureRevolutions: 4,
    results,
    combined: combineSeeds(results),
    thresholdsScaled: [0.03],
    reference: null,
    createdAt: "2026-01-01T00:00:00.000Z",
  };

  it("CSV contains software and basis comment lines", () => {
    const csv = buildBatchCsv(rec);
    expect(csv).toContain("# software: ");
    expect(csv).toContain("# basis: ");
    expect(csv).toContain("0.547");
  });

  it("JSON has model_basis and software_version", () => {
    const j = JSON.parse(buildBatchJson(rec));
    expect(Array.isArray(j.model_basis)).toBe(true);
    expect(typeof j.software_version).toBe("string");
  });
});
