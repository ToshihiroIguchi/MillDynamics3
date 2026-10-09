// Pure logic for the "reference" feature of the Grinding proxy (3D) rows: store the window-mean of
// each proxy value of one run, then show current / reference ratios for later runs. The proxy
// values are meant for RELATIVE comparison only, so this is the main way to read them.
// No DOM here; localStorage access is isolated in `loadReference` / `saveReference`.

import type { ParamsJson } from "../protocol";

const STORAGE_KEY = "milldynamics.proxyReference";

/** Numerical settings that change impact statistics by several percent. */
export interface NumericalSettings {
  substeps: number | null;
  demIterations: number | null;
  coarseGrainingFactor: number | null;
}

export interface ProxyReference {
  /** Metric id -> reference value (window mean when a window was active, else the current value). */
  values: Record<string, number>;
  /** Numerical settings of the run that produced the reference. */
  settings: NumericalSettings;
  /** The full params that produced the reference (for the record). */
  params: ParamsJson | null;
  /** Sim time at which the reference was taken (s). */
  simTime: number;
  /** True when at least some values are window means rather than instantaneous values. */
  windowed: boolean;
  /** Optional per-row standard error of the mean (same keys as `values`); present only for
   * references stored from a repeat-runs batch (batch/). Absent in references stored earlier. */
  errors?: Record<string, number>;
  /** True when the reference was stored from a repeat-runs batch (mean over seeds). */
  fromBatch?: boolean;
  /** Number of seeds averaged when `fromBatch`. */
  seeds?: number;
}

function num(v: unknown): number | null {
  return typeof v === "number" && Number.isFinite(v) ? v : null;
}

export function settingsOf(params: ParamsJson | null, coarseGrainingFactor: number | null): NumericalSettings {
  const sim = params?.simulation as { substeps?: unknown; dem_iterations?: unknown } | undefined;
  return {
    substeps: num(sim?.substeps),
    demIterations: num(sim?.dem_iterations),
    coarseGrainingFactor: num(coarseGrainingFactor),
  };
}

/** Names of the settings that differ between two runs (empty when comparable). A setting that is
 * unknown (null) on either side is not reported. Coarse-graining factors are compared to 1e-6. */
export function settingsMismatch(a: NumericalSettings, b: NumericalSettings): string[] {
  const out: string[] = [];
  if (a.substeps !== null && b.substeps !== null && a.substeps !== b.substeps) out.push("sub-steps");
  if (a.demIterations !== null && b.demIterations !== null && a.demIterations !== b.demIterations) {
    out.push("ball solver iterations");
  }
  if (
    a.coarseGrainingFactor !== null &&
    b.coarseGrainingFactor !== null &&
    Math.abs(a.coarseGrainingFactor - b.coarseGrainingFactor) > 1e-6
  ) {
    out.push("coarse-graining factor");
  }
  return out;
}

/** current / reference, or null if either is missing/non-finite or the reference is 0. */
export function ratio(current: number | null | undefined, reference: number | null | undefined): number | null {
  if (current === null || current === undefined || reference === null || reference === undefined) return null;
  if (!Number.isFinite(current) || !Number.isFinite(reference) || reference === 0) return null;
  return current / reference;
}

/** Formats a ratio with 2 significant digits, e.g. 0.8734 -> "0.87", 12.3 -> "12", 1234 -> "1200". */
export function formatRatio(r: number | null): string {
  if (r === null || !Number.isFinite(r)) return "-";
  if (r === 0) return "0";
  const abs = Math.abs(r);
  if (abs >= 1e5 || abs < 1e-3) return r.toExponential(1);
  return String(Number(r.toPrecision(2)));
}

/** Builds a reference from `ids`: window mean where `windowMean(id)` has one, else `current[id]`.
 * Ids with neither are omitted. */
export function buildReference(
  ids: readonly string[],
  current: Record<string, number | null | undefined>,
  windowMean: (id: string) => number | null,
  params: ParamsJson | null,
  coarseGrainingFactor: number | null,
  simTime: number,
): ProxyReference {
  const values: Record<string, number> = {};
  let windowed = false;
  for (const id of ids) {
    const w = windowMean(id);
    if (w !== null && Number.isFinite(w)) {
      values[id] = w;
      windowed = true;
      continue;
    }
    const c = current[id];
    if (c !== null && c !== undefined && Number.isFinite(c)) values[id] = c;
  }
  return { values, settings: settingsOf(params, coarseGrainingFactor), params, simTime, windowed };
}

/** Validates a parsed JSON value as a `ProxyReference`; returns null when malformed. */
export function parseReference(raw: unknown): ProxyReference | null {
  if (raw === null || typeof raw !== "object") return null;
  const o = raw as Record<string, unknown>;
  if (o.values === null || typeof o.values !== "object") return null;
  const values: Record<string, number> = {};
  for (const [k, v] of Object.entries(o.values as Record<string, unknown>)) {
    if (typeof v === "number" && Number.isFinite(v)) values[k] = v;
  }
  if (Object.keys(values).length === 0) return null;
  const st = (o.settings ?? {}) as Record<string, unknown>;
  let errors: Record<string, number> | undefined;
  if (o.errors !== null && typeof o.errors === "object") {
    const e: Record<string, number> = {};
    for (const [k, v] of Object.entries(o.errors as Record<string, unknown>)) {
      if (typeof v === "number" && Number.isFinite(v) && v >= 0 && k in values) e[k] = v;
    }
    if (Object.keys(e).length > 0) errors = e;
  }
  const seeds = num(o.seeds);
  return {
    values,
    ...(errors ? { errors } : {}),
    ...(o.fromBatch === true ? { fromBatch: true } : {}),
    ...(seeds !== null && seeds >= 1 ? { seeds } : {}),
    settings: {
      substeps: num(st.substeps),
      demIterations: num(st.demIterations),
      coarseGrainingFactor: num(st.coarseGrainingFactor),
    },
    params: o.params !== null && typeof o.params === "object" ? (o.params as ParamsJson) : null,
    simTime: num(o.simTime) ?? 0,
    windowed: o.windowed === true,
  };
}

export function loadReference(): ProxyReference | null {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    return raw === null ? null : parseReference(JSON.parse(raw));
  } catch {
    return null;
  }
}

export function saveReference(ref: ProxyReference | null): void {
  try {
    if (ref === null) localStorage.removeItem(STORAGE_KEY);
    else localStorage.setItem(STORAGE_KEY, JSON.stringify(ref));
  } catch {
    // ignore (private browsing / disabled storage)
  }
}
