// Tolerant parsing of the `Metrics` JSON produced by `Simulation::metrics_json()`. Only the
// optional `stress` block is validated here (the rest is trusted as before): a missing or
// malformed `stress` becomes `undefined` instead of poisoning the rows that read it.

import type { Metrics, StressMetrics } from "./types";

const NUMBER_FIELDS = [
  "ball_mass_3d_kg",
  "scale_energy_j",
  "impacts_total_per_ball_per_s",
  "mean_se_scaled_drop",
  "mean_se_scaled_shear",
  "stress_intensity_scaled",
  "wall_power_per_ball_w",
  "specific_power_w_per_kg",
  "cumulative_wall_energy_per_ball_j",
  "cumulative_specific_energy_j_per_kg",
] as const;

function isNumberArray(v: unknown, length?: number): v is number[] {
  return Array.isArray(v) && (length === undefined || v.length === length) && v.every((x) => typeof x === "number");
}

/** Returns a validated `StressMetrics`, or `undefined` if `raw` is missing or does not have the
 * expected shape (wrong types, per-class arrays not of length 4, per-threshold arrays whose length
 * differs from `thresholds_scaled`). Null numbers (serde's encoding of NaN/inf) count as invalid. */
export function parseStress(raw: unknown): StressMetrics | undefined {
  if (raw === null || typeof raw !== "object") return undefined;
  const o = raw as Record<string, unknown>;
  for (const f of NUMBER_FIELDS) {
    if (typeof o[f] !== "number") return undefined;
  }
  if (!isNumberArray(o.impacts_per_ball_per_s, 4) || !isNumberArray(o.mean_se_scaled, 4)) return undefined;
  if (!isNumberArray(o.thresholds_scaled)) return undefined;
  const n = o.thresholds_scaled.length;
  if (!isNumberArray(o.effective_power_per_ball_w, n) || !isNumberArray(o.cumulative_effective_energy_per_ball_j, n)) {
    return undefined;
  }
  return o as unknown as StressMetrics;
}

/** Parses `Simulation::metrics_json()` text; `stress` is normalised via `parseStress`. */
export function parseMetrics(json: string): Metrics {
  const m = JSON.parse(json) as Metrics;
  m.stress = parseStress((m as { stress?: unknown }).stress);
  return m;
}
