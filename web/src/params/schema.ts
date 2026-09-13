// Parameter field definitions driving the parameters modal (src/ui/paramsModal.ts). Mirrors
// crates/mill-core/src/params.rs field-for-field for the groups implemented so far (Mill, Media,
// Simulation); Slurry/Lifters/Display land with M2/M3/M5 (see docs/PLAN.md ss4.3).
//
// v1 (M1) simplification: every field is treated as requiring a full simulation reset on Apply
// (no "hot" live-apply distinction yet -- see docs/PLAN.md ss4.1's hot/reset split, completed in
// M5). This is safe by construction: mill-core's Simulation::set_params does not reseed the ball
// population even for changes that would need it, so defaulting to "always reset" avoids a class
// of silent-desync bugs while that distinction isn't implemented on the UI side yet.

import type { ParamsJson } from "../protocol";

export type FieldType = "number" | "select";

export interface SelectOption {
  value: string;
  label: string;
}

export interface FieldSchema {
  /** Dot path into the Params JSON tree, e.g. "mill.diameter_m". */
  path: string;
  group: "Mill" | "Media" | "Simulation";
  label: string;
  unit?: string;
  type: FieldType;
  min?: number;
  max?: number;
  step?: number;
  options?: SelectOption[];
}

export const GROUPS: FieldSchema["group"][] = ["Mill", "Media", "Simulation"];

export const SCHEMA: FieldSchema[] = [
  // Mill
  { path: "mill.diameter_m", group: "Mill", label: "Drum diameter", unit: "m", type: "number", min: 0.1, max: 5, step: 0.05 },
  {
    path: "mill.speed_mode",
    group: "Mill",
    label: "Speed mode",
    type: "select",
    options: [
      { value: "percent_critical", label: "% of critical speed" },
      { value: "rpm", label: "rpm" },
    ],
  },
  { path: "mill.speed_value", group: "Mill", label: "Speed", type: "number", min: 0, max: 300, step: 1 },
  {
    path: "mill.direction",
    group: "Mill",
    label: "Direction",
    type: "select",
    options: [
      { value: "counter_clockwise", label: "Counter-clockwise" },
      { value: "clockwise", label: "Clockwise" },
    ],
  },
  // Media
  { path: "media.ball_diameter_m", group: "Media", label: "Ball diameter", unit: "m", type: "number", min: 0.0005, max: 0.2, step: 0.0005 },
  { path: "media.fill_fraction", group: "Media", label: "Fill fraction (J)", type: "number", min: 0, max: 0.9, step: 0.01 },
  { path: "media.density_kg_m3", group: "Media", label: "Media density", unit: "kg/m3", type: "number", min: 100, max: 20000, step: 100 },
  { path: "media.restitution_ball_ball", group: "Media", label: "Restitution (ball-ball)", type: "number", min: 0, max: 1, step: 0.01 },
  { path: "media.restitution_ball_wall", group: "Media", label: "Restitution (ball-wall)", type: "number", min: 0, max: 1, step: 0.01 },
  { path: "media.friction_ball_ball", group: "Media", label: "Friction (ball-ball)", type: "number", min: 0, max: 2, step: 0.01 },
  { path: "media.friction_ball_wall", group: "Media", label: "Friction (ball-wall)", type: "number", min: 0, max: 2, step: 0.01 },
  { path: "media.rolling_friction", group: "Media", label: "Rolling friction", type: "number", min: 0, max: 1, step: 0.001 },
  // Simulation
  { path: "simulation.substeps", group: "Simulation", label: "Sub-steps / frame", type: "number", min: 1, max: 16, step: 1 },
  { path: "simulation.dem_iterations", group: "Simulation", label: "Ball solver iterations", type: "number", min: 1, max: 20, step: 1 },
  { path: "simulation.max_balls", group: "Simulation", label: "Max balls (coarse-graining target)", type: "number", min: 10, max: 50000, step: 10 },
  { path: "simulation.time_scale", group: "Simulation", label: "Time scale", type: "number", min: 0.1, max: 5, step: 0.1 },
  { path: "simulation.seed", group: "Simulation", label: "Random seed", type: "number", min: 0, max: 1_000_000_000, step: 1 },
];

/** Reads a dot-path (e.g. "mill.diameter_m") out of a nested params object. */
export function getPath(params: ParamsJson, path: string): unknown {
  return path.split(".").reduce<unknown>((node, key) => {
    if (node && typeof node === "object") {
      return (node as Record<string, unknown>)[key];
    }
    return undefined;
  }, params);
}

/** Returns a deep-cloned copy of `params` with the value at `path` replaced by `value`. */
export function withPath(params: ParamsJson, path: string, value: unknown): ParamsJson {
  const clone = structuredClone(params);
  const keys = path.split(".");
  const last = keys.pop();
  if (last === undefined) return clone;
  let node = clone as Record<string, unknown>;
  for (const key of keys) {
    const next = node[key];
    if (!next || typeof next !== "object") {
      throw new Error(`withPath: missing intermediate object at "${key}" for path "${path}"`);
    }
    node = next as Record<string, unknown>;
  }
  node[last] = value;
  return clone;
}
