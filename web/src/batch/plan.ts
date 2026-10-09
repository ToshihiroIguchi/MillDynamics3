// Pure helpers for planning a repeat run: input validation, the settings summary, refusal rules,
// and a rough duration estimate. No DOM here.

import { rpmOf } from "../params/derived";
import type { ParamsJson } from "../protocol";

export interface BatchInputs {
  seeds: number;
  warmup: number;
  measure: number;
}

export const MAX_SEEDS = 50;

/** Parses the three form fields; returns inputs or a message naming the first problem. */
export function validateInputs(seeds: string, warmup: string, measure: string): BatchInputs | string {
  const s = Number(seeds);
  const w = Number(warmup);
  const m = Number(measure);
  if (seeds.trim() === "" || !Number.isInteger(s) || s < 1 || s > MAX_SEEDS) {
    return `Seeds must be an integer from 1 to ${MAX_SEEDS}.`;
  }
  if (warmup.trim() === "" || !Number.isFinite(w) || w < 0) return "Warm-up revolutions must be a number >= 0.";
  if (measure.trim() === "" || !Number.isInteger(m) || m < 1) return "Measured revolutions must be an integer >= 1.";
  return { seeds: s, warmup: w, measure: m };
}

/** Seeds are the literal numbers 1..N (never offset by the main simulation's seed). */
export function seedList(n: number): number[] {
  return Array.from({ length: n }, (_, i) => i + 1);
}

type Section = Record<string, unknown>;
const sec = (p: ParamsJson | null, key: string): Section => (p?.[key] as Section | undefined) ?? {};
const numOf = (v: unknown): number | null => (typeof v === "number" && Number.isFinite(v) ? v : null);

/** Rotation speed (rpm) of the params, or null if the mill section is incomplete. */
export function rpmFromParams(params: ParamsJson | null): number | null {
  const mill = sec(params, "mill");
  const d = numOf(mill.diameter_m);
  const v = numOf(mill.speed_value);
  if (d === null || v === null || typeof mill.speed_mode !== "string") return null;
  return rpmOf({ diameter_m: d, speed_mode: mill.speed_mode, speed_value: v });
}

/** One-line summary of the settings a run uses, e.g.
 * "63 mm drum, 2 mm balls, 100 rpm, sub-steps 32, iterations 12, slurry on". */
export function summarizeParams(params: ParamsJson | null): string {
  if (!params) return "no parameters yet";
  const mill = sec(params, "mill");
  const media = sec(params, "media");
  const sim = sec(params, "simulation");
  const slurry = sec(params, "slurry");
  const drum = numOf(mill.diameter_m);
  const ball = numOf(media.ball_diameter_m);
  const rpm = rpmFromParams(params);
  const parts: string[] = [];
  if (drum !== null) parts.push(`${Number((drum * 1000).toPrecision(4))} mm drum`);
  if (ball !== null) parts.push(`${Number((ball * 1000).toPrecision(3))} mm balls`);
  if (rpm !== null) parts.push(`${Number(rpm.toPrecision(4))} rpm`);
  const sub = numOf(sim.substeps);
  if (sub !== null) parts.push(`sub-steps ${sub}`);
  const it = numOf(sim.dem_iterations);
  if (it !== null) parts.push(`iterations ${it}`);
  parts.push(`slurry ${slurry.enabled === true ? "on" : "off"}`);
  return parts.join(", ");
}

/** What the batch section needs to know about the live (main) simulation. */
export interface LiveSnapshot {
  /** The CURRENT APPLIED params of the main simulation (never unapplied form edits). */
  params: ParamsJson | null;
  coarseGrainingFactor: number | null;
  thresholdsScaled: readonly number[];
  subStepsPerSecondAchieved: number;
}

/** Reason a batch must not start, or null. */
export function refusalReason(live: LiveSnapshot | null): string | null {
  if (!live || !live.params) return "The simulation is not ready yet.";
  if (live.coarseGrainingFactor === null) return "No metrics from the main simulation yet; wait for the first frame.";
  if (live.coarseGrainingFactor > 1) {
    return (
      `Coarse-graining is active (k = ${live.coarseGrainingFactor.toFixed(2)}): the simulated ball is not the real ball, ` +
      `so the statistics are meaningless. Raise Max balls so that k = 1.`
    );
  }
  const rpm = rpmFromParams(live.params);
  if (rpm === null || !(rpm > 0)) return "The drum speed is 0; whole-revolution windows are undefined.";
  return null;
}

/** Rough wall time (s) for `runs` headless runs from the main simulation's achieved sub-steps/s, or
 * null if unknown (e.g. the main simulation is paused, so nothing is being achieved). A frame is
 * 1/60 s of simulated time and holds `substeps` sub-steps. */
export function estimateSeconds(
  params: ParamsJson | null,
  inputs: Pick<BatchInputs, "warmup" | "measure">,
  runs: number,
  subStepsPerSecondAchieved: number,
): number | null {
  const rpm = rpmFromParams(params);
  const sub = numOf(sec(params, "simulation").substeps);
  if (rpm === null || rpm <= 0 || sub === null || !(subStepsPerSecondAchieved > 0)) return null;
  const revS = 60 / rpm;
  const frames = (inputs.warmup + inputs.measure) * revS * 60;
  return (runs * frames * sub) / subStepsPerSecondAchieved;
}

export function formatDuration(s: number): string {
  if (s < 90) return `${Math.round(s)} s`;
  if (s < 5400) return `${Math.round(s / 60)} min`;
  return `${(s / 3600).toFixed(1)} h`;
}
