import { describe, expect, it } from "vitest";
import { getPath, paramsChangeRequiresReset, withPath, type FieldSchema, SCHEMA } from "../../src/params/schema";
import type { ParamsJson } from "../../src/protocol";

// A minimal but structurally-real params tree: one representative field per group SCHEMA
// actually touches, enough for withPath/getPath round-trips on every SCHEMA path below.
function baseParams(): ParamsJson {
  return {
    mill: { diameter_m: 1.0, speed_mode: "rpm", speed_value: 30, direction: "counter_clockwise" },
    media: {
      ball_diameter_m: 0.02,
      fill_fraction: 0.3,
      packing_fraction_2d: 0.82,
      density_kg_m3: 6000,
      restitution_ball_ball: 0.7,
      restitution_ball_wall: 0.7,
      friction_ball_ball: 0.3,
      friction_ball_wall: 0.3,
      rolling_friction: 0.01,
    },
    lifters: { count: 0, height_m: 0.02, base_width_m: 0.03, top_width_m: 0.02, phase_deg: 0 },
    slurry: {
      enabled: true,
      fill_fraction: 0.15,
      density_kg_m3: 1800,
      viscosity_pa_s: 50,
      wall_no_slip: 1.0,
      wettability: 0.6,
      surface_tension_n_m: 0.072,
      dye_pattern: "left_right",
    },
    simulation: { substeps: 8, dem_iterations: 2, max_balls: 600, resolution: 40, time_scale: 1.0, seed: 1 },
  };
}

describe("paramsChangeRequiresReset", () => {
  it("returns null when nothing changed", () => {
    const p = baseParams();
    expect(paramsChangeRequiresReset(p, p)).toBeNull();
  });

  it("returns null when only hot (resetRequired: false) fields changed", () => {
    const before = baseParams();
    let after = withPath(before, "mill.speed_value", 45);
    after = withPath(after, "media.friction_ball_ball", 0.5);
    after = withPath(after, "lifters.count", 8);
    after = withPath(after, "slurry.viscosity_pa_s", 100);
    after = withPath(after, "simulation.substeps", 4);
    expect(paramsChangeRequiresReset(before, after)).toBeNull();
  });

  it("names the field's label when a resetRequired field changed", () => {
    const before = baseParams();
    const after = withPath(before, "mill.diameter_m", 1.5);
    expect(paramsChangeRequiresReset(before, after)).toBe("Drum diameter");
  });

  it("still flags a reset field even when hot fields also changed alongside it", () => {
    const before = baseParams();
    let after = withPath(before, "mill.speed_value", 45); // hot
    after = withPath(after, "simulation.max_balls", 300); // resetRequired
    expect(paramsChangeRequiresReset(before, after)).toBe("Max balls (coarse-graining target)");
  });

  it("every SCHEMA path round-trips through getPath/withPath on the base fixture", () => {
    const params = baseParams();
    for (const field of SCHEMA) {
      const value = getPath(params, field.path);
      expect(value, `missing fixture value for ${field.path}`).not.toBeUndefined();
      const updated = withPath(params, field.path, value);
      expect(getPath(updated, field.path)).toEqual(value);
    }
  });

  it("every SCHEMA entry declares resetRequired explicitly", () => {
    for (const field of SCHEMA as FieldSchema[]) {
      expect(typeof field.resetRequired, `${field.path} is missing resetRequired`).toBe("boolean");
    }
  });
});
