// Contact-coefficient sensitivity: the base params plus ball-ball / ball-wall friction and
// restitution each scaled by 0.8 and 1.2. Pure; the caller runs one seed per condition.

import type { ParamsJson } from "../protocol";
import type { SeedResult } from "./stats";

export const SENSITIVITY_FACTORS = [0.8, 1.2] as const;
/** A change of a row beyond this fraction is highlighted. */
export const SENSITIVITY_WARN = 0.1;

interface Coefficient {
  id: string;
  label: string;
  kind: "friction" | "restitution";
  /** Kinetic / main field in `media`. */
  field: string;
  /** Static-friction field that must stay >= the kinetic one (friction only). */
  staticField?: string;
}

const COEFFICIENTS: readonly Coefficient[] = [
  { id: "fbb", label: "friction ball-ball", kind: "friction", field: "friction_ball_ball", staticField: "friction_ball_ball_static" },
  { id: "fbw", label: "friction ball-wall", kind: "friction", field: "friction_ball_wall", staticField: "friction_ball_wall_static" },
  { id: "ebb", label: "restitution ball-ball", kind: "restitution", field: "restitution_ball_ball" },
  { id: "ebw", label: "restitution ball-wall", kind: "restitution", field: "restitution_ball_wall" },
];

export interface Condition {
  id: string;
  label: string;
  params: ParamsJson;
}

/** Returns a copy of `params` with one coefficient scaled by `factor`. Restitution is clamped to
 * <= 1; a static friction field present in `media` is scaled too and kept >= the kinetic one. */
export function scaleCoefficient(params: ParamsJson, coeff: Coefficient, factor: number): ParamsJson {
  const clone = structuredClone(params);
  const media = clone.media as Record<string, unknown> | undefined;
  if (!media || typeof media[coeff.field] !== "number") {
    throw new Error(`params.media.${coeff.field} is missing`);
  }
  const base = media[coeff.field] as number;
  if (coeff.kind === "restitution") {
    media[coeff.field] = Math.min(1, base * factor);
    return clone;
  }
  const kinetic = base * factor;
  media[coeff.field] = kinetic;
  if (coeff.staticField && typeof media[coeff.staticField] === "number") {
    media[coeff.staticField] = Math.max((media[coeff.staticField] as number) * factor, kinetic);
  }
  return clone;
}

/** Base condition first, then each coefficient at x0.8 and x1.2. */
export function buildSensitivityConditions(params: ParamsJson): Condition[] {
  const out: Condition[] = [{ id: "base", label: "base", params: structuredClone(params) }];
  for (const c of COEFFICIENTS) {
    for (const f of SENSITIVITY_FACTORS) {
      out.push({ id: `${c.id}_x${f}`, label: `${c.label} x${f}`, params: scaleCoefficient(params, c, f) });
    }
  }
  return out;
}

/** (value - base) / base as a fraction, or null when the base is 0 or a value is non-finite. */
export function relativeChange(base: number, value: number): number | null {
  if (!Number.isFinite(base) || !Number.isFinite(value) || base === 0) return null;
  return (value - base) / Math.abs(base);
}

export interface SensitivityRow {
  name: string;
  base: number;
  /** One entry per non-base condition, in condition order: fractional change or null. */
  changes: (number | null)[];
}

/** `results[0]` is the base; all results must share names. */
export function buildSensitivityRows(results: readonly SeedResult[]): SensitivityRow[] {
  const base = results[0];
  if (!base) return [];
  return base.names.map((name, k) => ({
    name,
    base: base.values[k]!,
    changes: results.slice(1).map((r) => {
      if (r.names[k] !== name) throw new Error("sensitivity results have different scalar names");
      return relativeChange(base.values[k]!, r.values[k]!);
    }),
  }));
}

export function exceedsSensitivity(change: number | null): boolean {
  return change !== null && Math.abs(change) > SENSITIVITY_WARN;
}
