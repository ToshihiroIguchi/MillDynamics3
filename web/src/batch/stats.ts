// Pure statistics for repeat runs: combining per-seed results (mirrors `batch::combine` in
// crates/mill-core/src/batch.rs) and propagating relative errors into ratios. No DOM here.

/** Mirrors `batch::SeedResult` (crates/mill-core/src/batch.rs), as serialised by `stressRunJson`. */
export interface SeedResult {
  seed: number;
  window_s: number;
  coarse_graining_factor: number;
  names: string[];
  values: number[];
  power_drift_rel_per_rev: number;
}

export interface Combined {
  names: string[];
  mean: number[];
  /** Standard error of the mean: sample standard deviation (n - 1) over sqrt(n); 0 for one seed. */
  stdErr: number[];
  seeds: number;
}

/** Parses and validates the JSON returned by `stressRunJson`; throws a descriptive error. */
export function parseSeedResult(json: string): SeedResult {
  const o = JSON.parse(json) as Partial<SeedResult> | null;
  if (
    !o ||
    typeof o.seed !== "number" ||
    !Array.isArray(o.names) ||
    !Array.isArray(o.values) ||
    o.names.length !== o.values.length ||
    typeof o.coarse_graining_factor !== "number" ||
    typeof o.window_s !== "number"
  ) {
    throw new Error("Malformed batch result from the solver");
  }
  return {
    seed: o.seed,
    window_s: o.window_s,
    coarse_graining_factor: o.coarse_graining_factor,
    names: o.names as string[],
    values: o.values as number[],
    power_drift_rel_per_rev: typeof o.power_drift_rel_per_rev === "number" ? o.power_drift_rel_per_rev : 0,
  };
}

/** Mean and standard error over seeds; all results must share the same names. */
export function combineSeeds(results: readonly SeedResult[]): Combined {
  const first = results[0];
  if (!first) throw new Error("no results to combine");
  for (const r of results) {
    if (r.names.length !== first.names.length || r.names.some((n, i) => n !== first.names[i])) {
      throw new Error("results have different scalar names");
    }
    if (r.values.length !== r.names.length) throw new Error("result values do not match names");
  }
  const n = results.length;
  const mean: number[] = [];
  const stdErr: number[] = [];
  for (let k = 0; k < first.names.length; k++) {
    let sum = 0;
    for (const r of results) sum += r.values[k]!;
    const m = sum / n;
    let se = 0;
    if (n >= 2) {
      let ss = 0;
      for (const r of results) ss += (r.values[k]! - m) ** 2;
      se = Math.sqrt(ss / (n - 1) / n);
    }
    mean.push(m);
    stdErr.push(se);
  }
  return { names: [...first.names], mean, stdErr, seeds: n };
}

/** |se / value| as a fraction, or null when the value is 0 or either is non-finite. */
export function relativeError(value: number, se: number): number | null {
  if (!Number.isFinite(value) || !Number.isFinite(se) || value === 0) return null;
  return Math.abs(se / value);
}

export interface RatioWithError {
  ratio: number;
  /** Relative error of the ratio (fraction), or null if it cannot be determined. */
  relErr: number | null;
}

/** current / reference with first-order propagated relative error
 * sqrt((se_c / c)^2 + (se_r / r)^2). `refSe` null/undefined means "the reference has no error".
 * Returns null when the reference is missing or 0, or a value is non-finite. */
export function ratioWithError(
  current: number,
  currentSe: number,
  reference: number | undefined,
  refSe: number | null | undefined,
): RatioWithError | null {
  if (reference === undefined || !Number.isFinite(current) || !Number.isFinite(reference) || reference === 0) {
    return null;
  }
  const ratio = current / reference;
  if (current === 0 && currentSe === 0) return { ratio, relErr: 0 };
  const ec = relativeError(current, currentSe);
  const er = refSe === null || refSe === undefined ? 0 : relativeError(reference, refSe);
  if (ec === null || er === null) return { ratio, relErr: null };
  return { ratio, relErr: Math.sqrt(ec * ec + er * er) };
}

/** A ratio's relative error above this fraction is flagged as too uncertain to read. */
export const RATIO_ERROR_WARN = 0.1;
/** Wall-power drift above this fraction per revolution means the window is not steady state. */
export const DRIFT_WARN_PER_REV = 0.05;

export interface SteadyState {
  /** Mean of the signed per-seed drifts (fraction per revolution). */
  mean: number;
  meanAbs: number;
  maxAbs: number;
  warn: boolean;
}

export function steadyState(results: readonly SeedResult[]): SteadyState {
  const d = results.map((r) => r.power_drift_rel_per_rev);
  if (d.length === 0) return { mean: 0, meanAbs: 0, maxAbs: 0, warn: false };
  const mean = d.reduce((a, b) => a + b, 0) / d.length;
  const meanAbs = d.reduce((a, b) => a + Math.abs(b), 0) / d.length;
  const maxAbs = Math.max(...d.map(Math.abs));
  return { mean, meanAbs, maxAbs, warn: meanAbs > DRIFT_WARN_PER_REV };
}
