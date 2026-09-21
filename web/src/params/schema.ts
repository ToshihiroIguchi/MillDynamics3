// Parameter field definitions driving the parameters left panel (src/ui/paramsPanel.ts). Mirrors
// crates/mill-core/src/params.rs field-for-field for the groups implemented so far (Mill, Media,
// Lifters, Slurry, Simulation); a Display group lands with M5 (see docs/PLAN.md ss4.3).
// `slurry.rheology`/`yield_stress_pa` are intentionally not exposed: pbf.rs's v1 solver always
// runs an implicit (conjugate-gradient) Newtonian viscosity solve regardless of `rheology`
// (Bingham is unimplemented, M7), so surfacing those fields would imply behaviour that doesn't
// exist yet.
//
// **`resetRequired`**: whether a change to this field takes effect only through a full simulation
// reset (`InitMessage`, reconstructing `Simulation`) or can be applied live
// (`SetParamsMessage`/`Simulation::set_params`, which does not reseed/resize the ball or fluid
// population or the fluid lattice's baked-in spacing -- see that method's own doc comment in
// lib.rs for the exact, authoritative list this mirrors). Getting this wrong in the "should be
// `true`" direction produces a silent no-op (the UI shows the new value, nothing actually
// changes); getting it wrong in the "should be `false`" direction forces an unnecessary reset.
// `paramsChangeRequiresReset` below is what main.ts actually calls -- these annotations are its
// data, not consulted directly.

import type { ParamsJson } from "../protocol";

export type FieldType = "number" | "select" | "boolean";

export interface SelectOption {
  value: string;
  label: string;
}

export interface FieldSchema {
  /** Dot path into the Params JSON tree, e.g. "mill.diameter_m". */
  path: string;
  group: "Mill" | "Media" | "Lifters" | "Slurry" | "Simulation";
  label: string;
  unit?: string;
  type: FieldType;
  min?: number;
  max?: number;
  step?: number;
  options?: SelectOption[];
  /**
   * Multiplier from the SI value stored in Params to the value shown/edited in the UI (e.g. 1000
   * for a meters field displayed in mm: `display = raw * displayScale`). `min`/`max`/`step` are
   * already expressed in display units. Omitted (or 1) means no conversion.
   */
  displayScale?: number;
  /** See this module's doc comment. */
  resetRequired: boolean;
}

/** Converts a raw (SI) field value to the value shown in the UI input. */
export function toDisplayValue(field: FieldSchema, raw: unknown): unknown {
  if (field.type === "number" && typeof raw === "number" && field.displayScale) {
    return raw * field.displayScale;
  }
  return raw;
}

/** Converts a UI input's displayed number back to the raw (SI) value stored in Params. */
export function fromDisplayValue(field: FieldSchema, display: number): number {
  return field.displayScale ? display / field.displayScale : display;
}

export const GROUPS: FieldSchema["group"][] = ["Mill", "Media", "Lifters", "Slurry", "Simulation"];

export const SCHEMA: FieldSchema[] = [
  // Mill
  { path: "mill.diameter_m", group: "Mill", label: "Drum diameter", unit: "m", type: "number", min: 0.03, max: 5, step: 0.01, resetRequired: true },
  {
    path: "mill.speed_mode",
    group: "Mill",
    label: "Speed mode",
    type: "select",
    options: [
      { value: "percent_critical", label: "% of critical speed" },
      { value: "rpm", label: "rpm" },
    ],
    resetRequired: false,
  },
  { path: "mill.speed_value", group: "Mill", label: "Speed", type: "number", min: 0, max: 300, step: 1, resetRequired: false },
  {
    path: "mill.direction",
    group: "Mill",
    label: "Direction",
    type: "select",
    options: [
      { value: "counter_clockwise", label: "Counter-clockwise" },
      { value: "clockwise", label: "Clockwise" },
    ],
    resetRequired: false,
  },
  // Media: geometry/mass fields (baked into the ball population at seed time, params.rs's
  // `Params::effective_media`) need a reset; the physics coefficients below them are read fresh
  // every sub-step (`dem.rs`'s `step_with_external_forces`) and hot-apply correctly.
  { path: "media.ball_diameter_m", group: "Media", label: "Ball diameter", unit: "mm", type: "number", min: 0.5, max: 200, step: 0.1, displayScale: 1000, resetRequired: true },
  { path: "media.fill_fraction", group: "Media", label: "Fill fraction (J)", type: "number", min: 0, max: 0.9, step: 0.01, resetRequired: true },
  { path: "media.packing_fraction_2d", group: "Media", label: "2D packing fraction", type: "number", min: 0.5, max: 0.907, step: 0.001, resetRequired: true },
  { path: "media.density_kg_m3", group: "Media", label: "Media density", unit: "kg/m3", type: "number", min: 100, max: 20000, step: 100, resetRequired: true },
  { path: "media.restitution_ball_ball", group: "Media", label: "Restitution (ball-ball)", type: "number", min: 0, max: 1, step: 0.01, resetRequired: false },
  { path: "media.restitution_ball_wall", group: "Media", label: "Restitution (ball-wall)", type: "number", min: 0, max: 1, step: 0.01, resetRequired: false },
  { path: "media.friction_ball_ball", group: "Media", label: "Friction (ball-ball)", type: "number", min: 0, max: 2, step: 0.01, resetRequired: false },
  { path: "media.friction_ball_wall", group: "Media", label: "Friction (ball-wall)", type: "number", min: 0, max: 2, step: 0.01, resetRequired: false },
  { path: "media.rolling_friction", group: "Media", label: "Rolling friction", type: "number", min: 0, max: 1, step: 0.001, resetRequired: false },
  // Lifters (count = 0 is the default: a perfectly smooth wall, see docs/PLAN.md ss3.1). The drum
  // geometry these describe is re-derived from `params.lifters` every sub-step (lib.rs), not baked
  // into any seeded state, so every field here hot-applies -- including `count` itself, letting a
  // user watch the charge's transient response to lifters appearing/disappearing mid-run.
  { path: "lifters.count", group: "Lifters", label: "Lifter count", type: "number", min: 0, max: 64, step: 1, resetRequired: false },
  { path: "lifters.height_m", group: "Lifters", label: "Height", unit: "m", type: "number", min: 0, max: 0.2, step: 0.005, resetRequired: false },
  { path: "lifters.base_width_m", group: "Lifters", label: "Base width", unit: "m", type: "number", min: 0, max: 0.3, step: 0.005, resetRequired: false },
  { path: "lifters.top_width_m", group: "Lifters", label: "Top width", unit: "m", type: "number", min: 0, max: 0.3, step: 0.005, resetRequired: false },
  { path: "lifters.phase_deg", group: "Lifters", label: "Phase offset", unit: "deg", type: "number", min: -180, max: 180, step: 1, resetRequired: false },
  // Slurry: `enabled`/`fill_fraction`/`density_kg_m3` are only read once, at `seed_lattice` time
  // (pbf.rs) -- toggling/changing them live neither seeds, clears, nor re-derives the existing
  // fluid population's mass/rest density (`FluidParticles::rest_density`/`particle_mass` are baked
  // in then, not re-read from `SlurryParams` on every step the way viscosity/no-slip/wettability
  // are). `dye_pattern` is likewise baked per-particle at seed time (no live re-tagging exists).
  { path: "slurry.enabled", group: "Slurry", label: "Enabled", type: "boolean", resetRequired: true },
  { path: "slurry.fill_fraction", group: "Slurry", label: "Fill fraction", type: "number", min: 0, max: 0.9, step: 0.01, resetRequired: true },
  { path: "slurry.density_kg_m3", group: "Slurry", label: "Slurry density", unit: "kg/m3", type: "number", min: 100, max: 5000, step: 50, resetRequired: true },
  // `step` is a UI granularity hint only (spinner increment); it is never enforced as a validity
  // constraint (the form uses novalidate and its own JS range check -- see ui/paramsPanel.ts).
  { path: "slurry.viscosity_pa_s", group: "Slurry", label: "Viscosity", unit: "Pa·s", type: "number", min: 0, max: 200, step: 0.1, resetRequired: false },
  { path: "slurry.wall_no_slip", group: "Slurry", label: "Wall no-slip (β)", type: "number", min: 0, max: 1, step: 0.05, resetRequired: false },
  { path: "slurry.wettability", group: "Slurry", label: "Wettability", type: "number", min: 0, max: 1, step: 0.05, resetRequired: false },
  // Real water is 0.072 N/m; see this field's params.rs doc comment for why the mapping into the
  // solver's internal cohesion strength is a calibrated proportionality, not a unit conversion.
  { path: "slurry.surface_tension_n_m", group: "Slurry", label: "Surface tension", unit: "N/m", type: "number", min: 0, max: 0.2, step: 0.005, resetRequired: false },
  {
    path: "slurry.dye_pattern",
    group: "Slurry",
    label: "Dye pattern",
    type: "select",
    options: [
      { value: "left_right", label: "Left / right" },
      { value: "top_bottom", label: "Top / bottom" },
      { value: "none", label: "None" },
    ],
    resetRequired: true,
  },
  // Simulation: `substeps`/`dem_iterations`/`time_scale` (and `pbf_iterations`, not yet exposed
  // here) are read fresh every sub-step/call; `max_balls`, `seed`, and `resolution` are baked into
  // the ball/fluid population at seed time (same reason as `media.*`/`slurry.*` geometry above --
  // `resolution` sets the fluid lattice's spacing/kernel radius/particle mass in `seed_lattice`).
  { path: "simulation.substeps", group: "Simulation", label: "Sub-steps / frame", type: "number", min: 1, max: 16, step: 1, resetRequired: false },
  { path: "simulation.dem_iterations", group: "Simulation", label: "Ball solver iterations", type: "number", min: 1, max: 20, step: 1, resetRequired: false },
  { path: "simulation.max_balls", group: "Simulation", label: "Max balls (coarse-graining target)", type: "number", min: 10, max: 50000, step: 10, resetRequired: true },
  { path: "simulation.resolution", group: "Simulation", label: "Slurry resolution (particles/radius)", type: "number", min: 4, max: 200, step: 1, resetRequired: true },
  { path: "simulation.time_scale", group: "Simulation", label: "Time scale", type: "number", min: 0.1, max: 5, step: 0.1, resetRequired: false },
  { path: "simulation.seed", group: "Simulation", label: "Random seed", type: "number", min: 0, max: 1_000_000_000, step: 1, resetRequired: true },
];

/**
 * Whether applying `newParams` over `oldParams` needs a full simulation reset (`InitMessage`)
 * rather than a hot apply (`SetParamsMessage`) -- true if any field marked `resetRequired: true`
 * in `SCHEMA` actually changed value. Returns the *label* of the first such field found (for a
 * status message naming the cause), or `null` if every changed field can be hot-applied.
 *
 * Deliberately walks `SCHEMA`, not a generic deep-diff of the two param trees: fields not exposed
 * in the UI (e.g. `simulation.resolution`, `slurry.rheology`) cannot have changed via this panel,
 * and a generic diff would have no `resetRequired` annotation for them to consult anyway.
 */
export function paramsChangeRequiresReset(oldParams: ParamsJson, newParams: ParamsJson): string | null {
  for (const field of SCHEMA) {
    if (!field.resetRequired) continue;
    const before = getPath(oldParams, field.path);
    const after = getPath(newParams, field.path);
    if (before !== after) return field.label;
  }
  return null;
}

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
  if (!(last in node)) {
    throw new Error(`withPath: missing leaf key "${last}" for path "${path}"`);
  }
  node[last] = value;
  return clone;
}
