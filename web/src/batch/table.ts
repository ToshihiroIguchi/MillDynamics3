// Pure construction of the repeat-runs result table and of the reference a batch can store.

import { settingsOf, type ProxyReference } from "../metrics/reference";
import type { ParamsJson } from "../protocol";
import { referenceKey, rowInfo, type RowInfo } from "./rows";
import { RATIO_ERROR_WARN, ratioWithError, relativeError, type Combined } from "./stats";

export interface TableRow {
  info: RowInfo;
  mean: number;
  se: number;
  /** Relative standard error (fraction), null when the mean is 0. */
  relSe: number | null;
  refValue: number | null;
  ratio: number | null;
  /** Propagated relative error of the ratio (fraction), null when unknown. */
  ratioRelErr: number | null;
  /** True when `ratioRelErr` exceeds `RATIO_ERROR_WARN`. */
  ratioUncertain: boolean;
}

export function buildTable(
  combined: Combined,
  thresholds: readonly number[],
  reference: ProxyReference | null,
): TableRow[] {
  return combined.names.map((name, i) => {
    const mean = combined.mean[i]!;
    const se = combined.stdErr[i]!;
    const key = referenceKey(name);
    const refValue = reference?.values[key] ?? null;
    // The reference carries no error unless it was itself stored from a batch.
    const refSe = reference?.fromBatch ? (reference.errors?.[key] ?? null) : null;
    const r = reference && refValue !== null ? ratioWithError(mean, se, refValue, refSe) : null;
    return {
      info: rowInfo(name, thresholds),
      mean,
      se,
      relSe: relativeError(mean, se),
      refValue,
      ratio: r?.ratio ?? null,
      ratioRelErr: r?.relErr ?? null,
      ratioUncertain: r?.relErr != null && r.relErr > RATIO_ERROR_WARN,
    };
  });
}

/** The reference a batch stores via "Use as reference": means + standard errors keyed like the
 * live rows (reference.ts), with the exact params/settings the batch used. */
export function referenceFromBatch(
  combined: Combined,
  params: ParamsJson,
  coarseGrainingFactor: number | null,
): ProxyReference {
  const values: Record<string, number> = {};
  const errors: Record<string, number> = {};
  combined.names.forEach((name, i) => {
    const mean = combined.mean[i]!;
    if (!Number.isFinite(mean)) return;
    const key = referenceKey(name);
    values[key] = mean;
    errors[key] = combined.stdErr[i]!;
  });
  return {
    values,
    errors,
    fromBatch: true,
    seeds: combined.seeds,
    settings: settingsOf(params, coarseGrainingFactor),
    params,
    simTime: 0,
    windowed: true,
  };
}
