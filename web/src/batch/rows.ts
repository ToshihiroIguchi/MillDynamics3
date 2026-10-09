// Row labels, units and number formatting for the repeat-runs table. Labels and units are reused
// from the live proxy rows (metrics/specs.ts) where the batch scalar names map onto a spec; the
// mapping is also what lets a batch be stored as the "Set reference" reference (reference.ts keys
// its values by spec id).

import { METRIC_SPECS } from "../metrics/specs";

export interface RowInfo {
  name: string;
  label: string;
  unit: string;
  /** Live metric id the row corresponds to, or null for window_* rows (batch-only). */
  metricId: string | null;
}

const NAME_TO_SPEC_ID: Record<string, string> = {
  specific_power_w_per_kg: "proxy_specific_power",
  wall_power_per_ball_w: "proxy_wall_power",
  impacts_total_per_ball_per_s: "proxy_impacts_total",
  mean_se_scaled_drop: "proxy_mean_se_drop",
  mean_se_scaled_shear: "proxy_mean_se_shear",
  stress_intensity_scaled: "proxy_stress_intensity",
  "impacts_per_ball_per_s.ball_ball_drop": "proxy_impacts_bb_drop",
  "impacts_per_ball_per_s.ball_ball_shear": "proxy_impacts_bb_shear",
  "impacts_per_ball_per_s.ball_wall_drop": "proxy_impacts_bw_drop",
  "impacts_per_ball_per_s.ball_wall_shear": "proxy_impacts_bw_shear",
};

const EFF_POWER = /^effective_power_per_ball_w\.thr(\d+)$/;
const WINDOW_EFF_POWER = /^window_effective_power_per_ball_w\.thr(\d+)$/;

function thresholdText(thresholds: readonly number[], i: number): string {
  const t = thresholds[i];
  return t === undefined ? `threshold ${i}` : `SE > ${t} m·g·d`;
}

/** Live metric id (the key used in `ProxyReference.values`) for a batch scalar name, or null. */
export function metricIdOf(name: string): string | null {
  const direct = NAME_TO_SPEC_ID[name];
  if (direct) return direct;
  const m = EFF_POWER.exec(name);
  return m ? `proxy_eff_power_${m[1]}` : null;
}

/** Key under which a batch row is stored in / looked up from a reference. */
export function referenceKey(name: string): string {
  return metricIdOf(name) ?? name;
}

export function rowInfo(name: string, thresholds: readonly number[]): RowInfo {
  const metricId = metricIdOf(name);
  if (metricId) {
    const spec = METRIC_SPECS.find((s) => s.id === metricId);
    const m = EFF_POWER.exec(name);
    if (m) {
      return {
        name,
        metricId,
        label: `Eff. power, ${thresholdText(thresholds, Number(m[1]))}`,
        unit: spec?.unit ?? "W/ball",
      };
    }
    if (spec) return { name, metricId, label: spec.label, unit: spec.unit ?? "" };
  }
  if (name === "window_specific_power_w_per_kg") {
    return { name, metricId: null, label: "Window specific power (exact)", unit: "W/kg" };
  }
  const w = WINDOW_EFF_POWER.exec(name);
  if (w) {
    return {
      name,
      metricId: null,
      label: `Window eff. power (exact), ${thresholdText(thresholds, Number(w[1]))}`,
      unit: "W/ball",
    };
  }
  return { name, metricId: null, label: name, unit: "" };
}

/** Significant-digit formatter matching specs.ts's `sig`. */
export function formatSig(v: number, digits = 3): string {
  if (!Number.isFinite(v)) return "-";
  if (v === 0) return "0";
  const a = Math.abs(v);
  return a >= 1e5 || a < 1e-3 ? v.toExponential(digits - 1) : v.toPrecision(digits);
}

/** Percentage text of a fraction (0.123 -> "12.3 %"); "-" when null or non-finite. */
export function formatPercent(fraction: number | null, digits = 1): string {
  return fraction === null || !Number.isFinite(fraction) ? "-" : `${(fraction * 100).toFixed(digits)} %`;
}
