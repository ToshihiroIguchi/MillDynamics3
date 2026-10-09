import { describe, expect, it } from "vitest";
import { QUALITY_PRESETS, presetFieldValues, presetMatchesFields } from "../../src/params/presets";
import { SCHEMA } from "../../src/params/schema";

function preset(id: string) {
  const p = QUALITY_PRESETS.find((q) => q.id === id);
  if (!p) throw new Error(`missing preset ${id}`);
  return p;
}

describe("quality presets", () => {
  it("existing presets set exactly max_balls and resolution", () => {
    for (const id of ["realtime", "balanced", "accuracy"]) {
      expect(presetFieldValues(preset(id)).map((f) => f.path)).toEqual(["simulation.max_balls", "simulation.resolution"]);
    }
  });

  it("impact_stats sets the converged-statistics settings", () => {
    const values = Object.fromEntries(presetFieldValues(preset("impact_stats")).map((f) => [f.path, f.value]));
    expect(values).toEqual({
      "simulation.max_balls": "600",
      "simulation.resolution": "25",
      "simulation.substeps": "32",
      "simulation.dem_iterations": "12",
      "simulation.coarse_graining_mode": "off",
    });
    expect(preset("impact_stats").label).toBe("Impact statistics (slow)");
  });

  it("every preset field path exists in the schema and respects its range", () => {
    for (const p of QUALITY_PRESETS) {
      for (const { path, value } of presetFieldValues(p)) {
        const field = SCHEMA.find((f) => f.path === path);
        expect(field, path).toBeDefined();
        if (field!.type === "select") {
          expect(field!.options!.map((o) => o.value)).toContain(value);
        } else {
          expect(Number(value)).toBeGreaterThanOrEqual(field!.min!);
          expect(Number(value)).toBeLessThanOrEqual(field!.max!);
        }
      }
    }
  });

  it("matching ignores fields the preset does not define", () => {
    const form: Record<string, string> = {
      "simulation.max_balls": "300",
      "simulation.resolution": "25",
      "simulation.substeps": "32",
      "simulation.coarse_graining_mode": "off",
    };
    expect(presetMatchesFields(preset("realtime"), (p) => form[p])).toBe(true);
  });

  it("impact_stats stops matching when any of its fields is edited", () => {
    const base: Record<string, string> = {
      "simulation.max_balls": "600",
      "simulation.resolution": "25",
      "simulation.substeps": "32",
      "simulation.dem_iterations": "12",
      "simulation.coarse_graining_mode": "off",
    };
    expect(presetMatchesFields(preset("impact_stats"), (p) => base[p])).toBe(true);
    const edits: [string, string][] = [
      ["simulation.substeps", "16"],
      ["simulation.dem_iterations", "2"],
      ["simulation.coarse_graining_mode", "auto"],
      ["simulation.max_balls", "500"],
    ];
    for (const [path, other] of edits) {
      const edited = { ...base, [path]: other };
      expect(presetMatchesFields(preset("impact_stats"), (p) => edited[p]), path).toBe(false);
    }
    expect(presetMatchesFields(preset("impact_stats"), () => undefined)).toBe(false);
  });
});
