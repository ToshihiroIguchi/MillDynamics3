// Pure (DOM-free, worker-free) logic of the scale-up study: condition -> Params JSON, validation,
// ensemble combination, ratio tables with propagated errors, extrapolation warnings and CSV
// building. Everything here is unit-tested (tests/unit/scaleup.test.ts).

import { withPath } from "../params/schema";
import type { ParamsJson } from "../protocol";
import {
  CLASS_NAMES,
  type EnsembleResult,
  type Interval,
  type PowerLawFit,
  type ScaleupConfigJson,
  type SensitivityCase,
  type Summary,
} from "./types";

export const MAX_BALLS = 50_000;
export const MIN_BALLS = 2;
/** Ratio errors above this relative size get a note. */
export const RATIO_ERROR_NOTE_REL = 0.1;
export const DRIFT_WARN_REL_PER_REV = 0.05;
export const DISPLACEMENT_WARN = 0.5;
export const SENSITIVITY_HIGHLIGHT_PCT = 10;
export const EXTRAPOLATION_RANGE_FACTOR = 1.5;
export const EXTRAPOLATION_MIN_R2 = 0.95;
/** Fallback 2D packing fraction (MediaParams default) until the wasm defaults arrive. */
export const DEFAULT_PACKING_FRACTION_2D = 0.82;

export interface Condition {
  name: string;
  drumDiameterMm: number;
  ballDiameterMm: number;
  percentCritical: number;
  fillFraction: number;
  ballDensity: number;
  /** 0 = smooth wall. Lifter height/widths follow the ball diameter (see buildParams). */
  lifterCount: number;
  restitutionBallBall: number;
  restitutionBallWall: number;
  frictionBallBall: number;
  frictionBallWall: number;
  /** Lubrication-cushion viscosity in Pa*s; 0 = dry. */
  viscosityPaS: number;
}

export interface StudyConfig {
  seedCount: number;
  warmupRevolutions: number;
  measureRevolutions: number;
  samplesPerRevolution: number;
  dtRatio: number;
  demIterations: number;
}

/** Mirrors `ScaleupConfig::default()` in scaleup.rs plus 3 seeds. */
export const DEFAULT_STUDY_CONFIG: StudyConfig = {
  seedCount: 3,
  warmupRevolutions: 4,
  measureRevolutions: 6,
  samplesPerRevolution: 8,
  dtRatio: 0.05,
  demIterations: 6,
};

/** Lab-scale defaults mirroring `MediaParams`/`MillParams` defaults (63 mm drum, 2 mm ZrO2 balls). */
export function defaultCondition(name: string): Condition {
  return {
    name,
    drumDiameterMm: 63,
    ballDiameterMm: 2,
    percentCritical: 60,
    fillFraction: 0.3,
    ballDensity: 6000,
    lifterCount: 0,
    restitutionBallBall: 0.7,
    restitutionBallWall: 0.5,
    frictionBallBall: 0.25,
    frictionBallWall: 0.35,
    viscosityPaS: 0,
  };
}

/** Ball count of a condition (mirrors `Params::true_ball_count`). */
export function ballCount(c: Condition, packing: number = DEFAULT_PACKING_FRACTION_2D): number {
  if (!(c.ballDiameterMm > 0) || !(c.drumDiameterMm > 0)) return 0;
  const ratio = c.drumDiameterMm / c.ballDiameterMm;
  return c.fillFraction * packing * ratio * ratio;
}

/** Ball count at a given D/d for a fill fraction. */
export function ballCountAtRatio(fillFraction: number, ratio: number, packing: number = DEFAULT_PACKING_FRACTION_2D): number {
  return fillFraction * packing * ratio * ratio;
}

export function ddRatio(c: Condition): number {
  return c.ballDiameterMm > 0 ? c.drumDiameterMm / c.ballDiameterMm : Number.NaN;
}

/** Human-readable problems that make `c` impossible to run (empty = OK). */
export function validateCondition(c: Condition, packing: number = DEFAULT_PACKING_FRACTION_2D): string[] {
  const problems: string[] = [];
  const positive: [string, number][] = [
    ["drum diameter", c.drumDiameterMm],
    ["ball diameter", c.ballDiameterMm],
    ["% critical speed", c.percentCritical],
    ["fill fraction", c.fillFraction],
    ["ball density", c.ballDensity],
  ];
  for (const [label, v] of positive) {
    if (!Number.isFinite(v) || v <= 0) problems.push(`${label} must be > 0`);
  }
  if (c.fillFraction > 0.9) problems.push("fill fraction must be <= 0.9");
  if (!Number.isInteger(c.lifterCount) || c.lifterCount < 0) problems.push("lifter count must be a non-negative integer");
  if (!Number.isFinite(c.viscosityPaS) || c.viscosityPaS < 0) problems.push("viscosity must be >= 0");
  const restitutions: [string, number][] = [
    ["restitution (ball-ball)", c.restitutionBallBall],
    ["restitution (ball-wall)", c.restitutionBallWall],
  ];
  for (const [label, v] of restitutions) {
    if (!Number.isFinite(v) || v < 0 || v > 1) problems.push(`${label} must be in [0, 1]`);
  }
  const frictions: [string, number][] = [
    ["friction (ball-ball)", c.frictionBallBall],
    ["friction (ball-wall)", c.frictionBallWall],
  ];
  for (const [label, v] of frictions) {
    if (!Number.isFinite(v) || v < 0) problems.push(`${label} must be >= 0`);
  }
  if (problems.length === 0) {
    const n = ballCount(c, packing);
    if (n > MAX_BALLS) problems.push(`${Math.round(n)} balls exceeds the ${MAX_BALLS} limit: use a smaller D/d`);
    if (n < MIN_BALLS) problems.push(`${n.toFixed(1)} balls is below the minimum of ${MIN_BALLS}: raise fill fraction or lower the ball diameter`);
  }
  return problems;
}

export function validateStudyConfig(g: StudyConfig): string[] {
  const problems: string[] = [];
  if (!Number.isInteger(g.seedCount) || g.seedCount < 1 || g.seedCount > 20) problems.push("seed count must be an integer in 1..20");
  if (!(g.warmupRevolutions >= 0)) problems.push("warm-up revolutions must be >= 0");
  if (!(g.measureRevolutions >= 1) || !Number.isInteger(g.measureRevolutions)) problems.push("measured revolutions must be a whole number >= 1");
  if (!Number.isInteger(g.samplesPerRevolution) || g.samplesPerRevolution < 1) problems.push("samples per revolution must be an integer >= 1");
  if (!(g.dtRatio > 0)) problems.push("dt ratio must be > 0");
  if (!Number.isInteger(g.demIterations) || g.demIterations < 1 || g.demIterations > 20) problems.push("DEM iterations must be an integer in 1..20");
  return problems;
}

/**
 * Builds the Params JSON of a condition from `base` (the wasm `defaultParamsJson()`), using the
 * schema's `withPath` helper so a renamed field fails loudly instead of being silently ignored.
 * Lifters (when count > 0) are sized from the ball diameter: height d, base width 1.5 d, top width d.
 */
export function buildParams(base: ParamsJson, c: Condition): ParamsJson {
  const d = c.ballDiameterMm / 1000;
  const edits: [string, unknown][] = [
    ["mill.diameter_m", c.drumDiameterMm / 1000],
    ["mill.speed_mode", "percent_critical"],
    ["mill.speed_value", c.percentCritical],
    ["media.ball_diameter_m", d],
    ["media.fill_fraction", c.fillFraction],
    ["media.density_kg_m3", c.ballDensity],
    ["media.restitution_ball_ball", c.restitutionBallBall],
    ["media.restitution_ball_wall", c.restitutionBallWall],
    ["media.friction_ball_ball", c.frictionBallBall],
    ["media.friction_ball_wall", c.frictionBallWall],
    ["media.friction_ball_ball_static", c.frictionBallBall],
    ["media.friction_ball_wall_static", c.frictionBallWall],
    ["media.wet_cushion_viscosity_pa_s", c.viscosityPaS],
    ["slurry.enabled", false],
    ["lifters.count", c.lifterCount],
    ["simulation.coarse_graining_mode", "off"],
    ["simulation.max_balls", MAX_BALLS],
  ];
  if (c.lifterCount > 0) {
    edits.push(["lifters.height_m", d], ["lifters.base_width_m", 1.5 * d], ["lifters.top_width_m", d]);
  }
  let p = base;
  for (const [path, value] of edits) p = withPath(p, path, value);
  return p;
}

export function seedList(count: number): number[] {
  return Array.from({ length: count }, (_, i) => i + 1);
}

/** The wasm `ScaleupConfig` for the given seeds (thresholds empty = the default m g d multiples). */
export function buildScaleupConfig(g: StudyConfig, seeds: number[]): ScaleupConfigJson {
  return {
    warmup_revolutions: g.warmupRevolutions,
    measure_revolutions: g.measureRevolutions,
    samples_per_revolution: g.samplesPerRevolution,
    dt_ratio: g.dtRatio,
    dem_iterations: g.demIterations,
    wet_viscosity_pa_s: 0,
    thresholds_j: [],
    seeds,
  };
}

/**
 * Combines single-seed ensembles (one wasm call each, so the UI can report per-seed progress) into
 * the same statistics `run_ensemble` produces: mean and sample standard error (n-1) per scalar,
 * the first seed's time series, and every seed's summary.
 */
export function combineEnsembles(parts: EnsembleResult[]): EnsembleResult {
  const head = parts[0];
  if (!head) throw new Error("combineEnsembles: no results");
  const n = parts.length;
  const mean: number[] = [];
  const stdErr: number[] = [];
  head.scalar_names.forEach((_, k) => {
    const col = parts.map((p) => p.mean[k] ?? Number.NaN);
    const m = col.reduce((a, b) => a + b, 0) / n;
    mean.push(m);
    if (n < 2) {
      stdErr.push(0);
    } else {
      const v = col.reduce((a, x) => a + (x - m) ** 2, 0) / (n - 1);
      stdErr.push(Math.sqrt(v / n));
    }
  });
  return {
    seeds: parts.flatMap((p) => p.seeds),
    scalar_names: head.scalar_names,
    mean,
    std_err: stdErr,
    first: head.first,
    summaries: parts.flatMap((p) => p.summaries),
  };
}

// ---------------------------------------------------------------------------------------------
// Metrics

export interface MetricDef {
  /** Stable key: the scalar name with the threshold's Joule suffix removed. */
  key: string;
  label: string;
  help: string;
}

/** `effective_power_per_ball_w.thr0(1.234e-3J)` -> `effective_power_per_ball_w.thr0`. */
export function metricKey(scalarName: string): string {
  return scalarName.replace(/\(.*\)$/, "");
}

const BASE_METRICS: MetricDef[] = [
  {
    key: "specific_power_w_per_kg",
    label: "Specific power (W/kg)",
    help: "Wall power per ball divided by the 3D ball mass; per-ball power uses a one-ball-thick slab (depth = ball diameter).",
  },
  {
    key: "impacts_total_per_ball_per_s",
    label: "Impacts per ball per s",
    help: "Total impact frequency (all four classes) per ball.",
  },
  {
    key: "mean_se_scaled_drop",
    label: "Mean SE, drop impacts (SE/mgd)",
    help: "Mean stress energy of drop-type impacts (normal approach speed >= slip speed). SE is the relative kinetic energy of an impact with the 3D ball mass; scaled = divided by m*g*d.",
  },
  {
    key: "mean_se_scaled_shear",
    label: "Mean SE, shear impacts (SE/mgd)",
    help: "Mean stress energy of shear-type impacts (slip speed dominates), scaled by m*g*d.",
  },
  {
    key: "stress_intensity_scaled",
    label: "Stress intensity (SE/mgd)",
    help: "Mean scaled stress energy over all impacts (Kwade-type stress intensity).",
  },
];

/** Metric rows for the table: the fixed ones plus one pair per energy threshold of `summary`. */
export function metricDefs(summary: Summary): MetricDef[] {
  const defs = [...BASE_METRICS];
  summary.thresholds_j.forEach((thr, i) => {
    const rel = summary.scale_energy_j > 0 ? thr / summary.scale_energy_j : Number.NaN;
    const relText = Number.isFinite(rel) ? `${formatNumber(rel, 3)} m*g*d` : `${formatNumber(thr, 3)} J`;
    defs.push({
      key: `effective_power_per_ball_w.thr${i}`,
      label: `Effective power, SE > ${relText} (W/ball)`,
      help: `Power carried by impacts whose stress energy exceeds ${relText} (one-ball-thick slab, per ball).`,
    });
    defs.push({
      key: `effective_power_fraction.thr${i}`,
      label: `Effective power fraction, SE > ${relText}`,
      help: `Fraction of the wall power carried by impacts above ${relText}.`,
    });
  });
  return defs;
}

/** Label of a metric key (falls back to the key itself for unknown ones). */
export function metricLabelForKey(key: string): string {
  const b = BASE_METRICS.find((m) => m.key === key);
  if (b) return b.label;
  const m = /^effective_power_(per_ball_w|fraction)\.thr(\d+)$/.exec(key);
  if (m) return `${m[1] === "fraction" ? "Effective power fraction" : "Effective power (W/ball)"}, threshold ${m[2]}`;
  return key;
}

/** The keys of `scalarNames` that appear in the metric table (fixed metrics plus effective-power pairs). */
export function primaryMetricKeys(scalarNames: string[]): string[] {
  const keys = scalarNames.map(metricKey);
  const fixed = BASE_METRICS.map((m) => m.key).filter((k) => keys.includes(k));
  const eff = keys.filter((k) => /^effective_power_(per_ball_w|fraction)\.thr\d+$/.test(k));
  return [...fixed, ...eff];
}

export interface MeanSe {
  mean: number;
  se: number;
}

export function getMetric(ens: EnsembleResult, key: string): MeanSe | null {
  const k = ens.scalar_names.findIndex((n) => metricKey(n) === key);
  if (k < 0) return null;
  return { mean: ens.mean[k] ?? Number.NaN, se: ens.std_err[k] ?? 0 };
}

export interface RatioResult {
  ratio: number;
  /** Absolute standard error of the ratio (first-order propagation). */
  err: number;
  /** err / |ratio|; Infinity when the ratio is 0 but uncertain. */
  relErr: number;
  /** True when relErr exceeds RATIO_ERROR_NOTE_REL. */
  noisy: boolean;
}

/** candidate / reference with first-order error propagation of both standard errors. */
export function ratioWithError(cand: MeanSe, ref: MeanSe): RatioResult {
  if (!Number.isFinite(cand.mean) || !Number.isFinite(ref.mean) || ref.mean === 0) {
    return { ratio: Number.NaN, err: Number.NaN, relErr: Number.NaN, noisy: false };
  }
  const ratio = cand.mean / ref.mean;
  const err = Math.hypot(cand.se / ref.mean, (cand.mean * ref.se) / (ref.mean * ref.mean));
  const relErr = ratio === 0 ? (err > 0 ? Infinity : 0) : err / Math.abs(ratio);
  return { ratio, err, relErr, noisy: relErr > RATIO_ERROR_NOTE_REL };
}

export interface Warning {
  level: "warn";
  text: string;
}

export interface SteadyIndicators {
  driftMax: number;
  relStdMax: number;
  displacementMax: number;
}

export function steadyIndicators(ens: EnsembleResult): SteadyIndicators {
  return {
    driftMax: Math.max(...ens.summaries.map((s) => Math.abs(s.power_drift_rel_per_rev))),
    relStdMax: Math.max(...ens.summaries.map((s) => s.power_rel_std)),
    displacementMax: Math.max(...ens.summaries.map((s) => s.max_substep_displacement_over_diameter)),
  };
}

/** Steady-state / numerical-validity warnings for one condition's ensemble. */
export function ensembleWarnings(name: string, ens: EnsembleResult): Warning[] {
  const out: Warning[] = [];
  const { driftMax, displacementMax } = steadyIndicators(ens);
  if (driftMax > DRIFT_WARN_REL_PER_REV) {
    out.push({
      level: "warn",
      text: `${name}: power drift ${(driftMax * 100).toFixed(1)}%/rev exceeds ${DRIFT_WARN_REL_PER_REV * 100}%/rev: the charge has not reached a steady state (raise warm-up revolutions).`,
    });
  }
  if (displacementMax > DISPLACEMENT_WARN) {
    out.push({
      level: "warn",
      text: `${name}: max sub-step displacement ${displacementMax.toFixed(2)} d exceeds ${DISPLACEMENT_WARN} d: contacts may tunnel (lower dt ratio).`,
    });
  }
  return out;
}

// ---------------------------------------------------------------------------------------------
// Sweep extrapolation and sensitivity

export interface Extrapolation {
  value: number;
  warnings: string[];
}

/** Predicts `fit` at `target` D/d and warns when that is far outside the swept range or the fit is poor. */
export function extrapolate(fit: PowerLawFit, target: number, swept: number[]): Extrapolation {
  const warnings: string[] = [];
  if (!(target > 0) || !Number.isFinite(target)) {
    return { value: Number.NaN, warnings: ["target D/d must be > 0"] };
  }
  if (swept.length > 0) {
    const lo = Math.min(...swept);
    const hi = Math.max(...swept);
    if (target > hi * EXTRAPOLATION_RANGE_FACTOR || target < lo / EXTRAPOLATION_RANGE_FACTOR) {
      warnings.push(
        `target D/d ${target} is more than ${EXTRAPOLATION_RANGE_FACTOR}x outside the swept range ${lo}-${hi}: extrapolation is unreliable`,
      );
    }
  }
  if (fit.r2 < EXTRAPOLATION_MIN_R2) {
    warnings.push(`power-law fit R^2 = ${fit.r2.toFixed(3)} is below ${EXTRAPOLATION_MIN_R2}: the scaling is not a clean power law`);
  }
  return { value: fit.prefactor * Math.pow(target, fit.exponent), warnings };
}

/** Percentage change of `value` relative to `base`; NaN when base is 0. */
export function percentChange(base: number, value: number): number {
  if (!Number.isFinite(base) || !Number.isFinite(value) || base === 0) return Number.NaN;
  return ((value - base) / Math.abs(base)) * 100;
}

export function isSensitive(pct: number): boolean {
  return Math.abs(pct) > SENSITIVITY_HIGHLIGHT_PCT;
}

/** Base-relative % change of metric `key` in sensitivity case `c`. */
export function sensitivityChange(base: SensitivityCase, c: SensitivityCase, key: string): number {
  const kb = base.scalar_names.findIndex((n) => metricKey(n) === key);
  const kc = c.scalar_names.findIndex((n) => metricKey(n) === key);
  if (kb < 0 || kc < 0) return Number.NaN;
  return percentChange(base.mean[kb] ?? Number.NaN, c.mean[kc] ?? Number.NaN);
}

// ---------------------------------------------------------------------------------------------
// Histogram helpers and CSV

/** Duration (s) of the recorded window of a run (recording time starts at 0). */
export function recordedDuration(intervals: Interval[]): number {
  const last = intervals[intervals.length - 1];
  return last ? last.t_s : 0;
}

/** Histogram of one impact class as impacts per ball per second per bin. */
export function histogramRates(histogramTotal: number[][], classIndex: number, nBalls: number, durationS: number): number[] {
  const row = histogramTotal[classIndex] ?? [];
  const norm = nBalls > 0 && durationS > 0 ? nBalls * durationS : Number.NaN;
  return row.map((c) => c / norm);
}

/** Geometric bin centres for log-spaced `edges`. */
export function binCentres(edges: number[]): number[] {
  const out: number[] = [];
  for (let i = 0; i + 1 < edges.length; i++) out.push(Math.sqrt((edges[i] ?? 0) * (edges[i + 1] ?? 0)));
  return out;
}

export function csvCell(v: string | number): string {
  if (typeof v === "number") return Number.isFinite(v) ? String(v) : "";
  return /[",\n\r]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v;
}

export function csvRow(cells: (string | number)[]): string {
  return cells.map(csvCell).join(",");
}

export interface CsvCondition {
  name: string;
  ens: EnsembleResult;
}

/** The metric table as CSV: per condition mean and standard error, plus ratio and error to condition 0. */
export function metricsCsv(conds: CsvCondition[]): string {
  const ref = conds[0];
  if (!ref) return "";
  const header: string[] = ["metric", "key"];
  conds.forEach((c, i) => {
    header.push(`${c.name} mean`, `${c.name} std_err`);
    if (i > 0) header.push(`${c.name} ratio_to_reference`, `${c.name} ratio_err`);
  });
  const lines = [csvRow(header)];
  for (const def of metricDefs(ref.ens.first.summary)) {
    const refM = getMetric(ref.ens, def.key);
    const row: (string | number)[] = [def.label, def.key];
    conds.forEach((c, i) => {
      const m = getMetric(c.ens, def.key);
      row.push(m ? m.mean : "", m ? m.se : "");
      if (i > 0) {
        const r = m && refM ? ratioWithError(m, refM) : null;
        row.push(r ? r.ratio : "", r ? r.err : "");
      }
    });
    lines.push(csvRow(row));
  }
  return lines.join("\n") + "\n";
}

/** Interval time series of several conditions in one CSV (first seed of each); `condition` column first. */
export function intervalsCsv(conds: CsvCondition[]): string {
  const first = conds[0];
  if (!first) return "";
  const nthr = first.ens.first.summary.thresholds_j.length;
  const header = ["condition", "t_s", "revolutions", "wall_power_w", "dissipated_power_w"];
  for (const n of CLASS_NAMES) header.push(`impacts_per_ball_per_s.${n}`);
  for (const n of CLASS_NAMES) header.push(`mean_se_scaled.${n}`);
  for (let i = 0; i < nthr; i++) header.push(`effective_power_w.thr${i}`);
  header.push("cum_wall_energy_j");
  for (let i = 0; i < nthr; i++) header.push(`cum_effective_energy_j.thr${i}`);
  const lines = [csvRow(header)];
  for (const c of conds) {
    for (const iv of c.ens.first.intervals) {
      lines.push(
        csvRow([
          c.name,
          iv.t_s,
          iv.revolutions,
          iv.wall_power_w,
          iv.dissipated_power_w,
          ...iv.impacts_per_ball_per_s,
          ...iv.mean_se_scaled,
          ...iv.effective_power_w,
          iv.cum_wall_energy_j,
          ...iv.cum_effective_energy_j,
        ]),
      );
    }
  }
  return lines.join("\n") + "\n";
}

/** Compact number formatting for tables: scientific for very small/large magnitudes. */
export function formatNumber(v: number, digits = 3): string {
  if (!Number.isFinite(v)) return "-";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-2) return v.toExponential(digits - 1);
  return Number(v.toPrecision(digits)).toString();
}
