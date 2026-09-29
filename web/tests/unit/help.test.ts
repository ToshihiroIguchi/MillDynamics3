import { describe, expect, it } from "vitest";
import { METRIC_SPECS } from "../../src/metrics/specs";
import { GROUPS, SCHEMA } from "../../src/params/schema";

describe("help text", () => {
  it("every SCHEMA field has non-empty help", () => {
    for (const field of SCHEMA) {
      expect(field.help?.trim(), `SCHEMA field ${field.path}`).toBeTruthy();
    }
  });

  it("every METRIC_SPECS entry has non-empty help", () => {
    for (const spec of METRIC_SPECS) {
      expect(spec.help?.trim(), `metric ${spec.id}`).toBeTruthy();
    }
  });

  it("uses the renamed Numerical accuracy group and no legacy Simulation group", () => {
    expect(GROUPS).toContain("Numerical accuracy");
    expect(GROUPS as string[]).not.toContain("Simulation");
    expect(SCHEMA.every((f) => (GROUPS as string[]).includes(f.group))).toBe(true);
  });
});
