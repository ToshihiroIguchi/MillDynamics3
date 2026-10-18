// Pure CSV / JSON builders for a finished repeat run (reproducible: includes the exact params JSON
// and the seeds used). The DOM download helper lives in ui/batchSection.ts.

import type { ProxyReference } from "../metrics/reference";
import type { ParamsJson } from "../protocol";
import { modelBasisLines } from "../report/assumptions";
import { versionLine } from "../version";
import { buildTable } from "./table";
import type { Combined, SeedResult } from "./stats";

export interface BatchRecord {
  params: ParamsJson;
  warmupRevolutions: number;
  measureRevolutions: number;
  results: SeedResult[];
  combined: Combined;
  thresholdsScaled: readonly number[];
  reference: ProxyReference | null;
  /** Creation time as an ISO string (passed in so the builders stay pure). */
  createdAt: string;
  /** Software version (defaults to this build's commit). */
  softwareVersion?: string;
  /** Model-basis statements (defaults to `modelBasisLines(params)`). */
  modelBasis?: string[];
}

function csvField(v: string): string {
  return /[",\n\r]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v;
}

function num(v: number | null): string {
  return v === null || !Number.isFinite(v) ? "" : String(v);
}

/** CSV with `#` comment lines (settings, seeds, params), then one row per scalar:
 * name,label,unit,seed<i>...,mean,std_err,rel_se_percent,reference,ratio,ratio_rel_err_percent. */
export function buildBatchCsv(rec: BatchRecord): string {
  const table = buildTable(rec.combined, rec.thresholdsScaled, rec.reference);
  const lines: string[] = [
    "# MillDynamics repeat runs (same Simulation and params as the main screen; only the seed differs)",
    `# created: ${rec.createdAt}`,
    `# software: ${rec.softwareVersion ?? versionLine()}`,
    `# seeds: ${rec.results.map((r) => r.seed).join(" ")}`,
    `# warmup_revolutions: ${rec.warmupRevolutions}`,
    `# measure_revolutions: ${rec.measureRevolutions}`,
    `# window_s: ${num(rec.results[0]?.window_s ?? null)}`,
    `# params_json: ${JSON.stringify(rec.params)}`,
    ...(rec.modelBasis ?? modelBasisLines(rec.params)).map((l) => `# basis: ${l}`),
  ];
  const header = [
    "name",
    "label",
    "unit",
    ...rec.results.map((r) => `seed${r.seed}`),
    "mean",
    "std_err",
    "rel_se_percent",
    "reference",
    "ratio",
    "ratio_rel_err_percent",
  ];
  lines.push(header.map(csvField).join(","));
  table.forEach((row, i) => {
    const cells = [
      row.info.name,
      row.info.label,
      row.info.unit,
      ...rec.results.map((r) => num(r.values[i] ?? null)),
      num(row.mean),
      num(row.se),
      num(row.relSe === null ? null : row.relSe * 100),
      num(row.refValue),
      num(row.ratio),
      num(row.ratioRelErr === null ? null : row.ratioRelErr * 100),
    ];
    lines.push(cells.map(csvField).join(","));
  });
  return lines.join("\n") + "\n";
}

/** Full-fidelity JSON: per-seed results, combined table, params and seeds. */
export function buildBatchJson(rec: BatchRecord): string {
  const table = buildTable(rec.combined, rec.thresholdsScaled, rec.reference);
  return JSON.stringify(
    {
      format: "milldynamics-repeat-runs",
      version: 1,
      createdAt: rec.createdAt,
      software_version: rec.softwareVersion ?? versionLine(),
      model_basis: rec.modelBasis ?? modelBasisLines(rec.params),
      seeds: rec.results.map((r) => r.seed),
      warmupRevolutions: rec.warmupRevolutions,
      measureRevolutions: rec.measureRevolutions,
      thresholdsScaled: rec.thresholdsScaled,
      params: rec.params,
      perSeed: rec.results,
      combined: rec.combined,
      table: table.map((r) => ({
        name: r.info.name,
        label: r.info.label,
        unit: r.info.unit,
        mean: r.mean,
        stdErr: r.se,
        relSe: r.relSe,
        reference: r.refValue,
        ratio: r.ratio,
        ratioRelErr: r.ratioRelErr,
      })),
    },
    null,
    2,
  );
}
