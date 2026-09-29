// Builds a client-side PDF "simulation report" (docs on why: a researcher wants a handful of
// representative values -- mean power draw, collision rate, etc. -- summarized in a form suitable
// for citing in a paper's methods/results section) from a `ReportAggregator`'s running statistics
// plus the current run's configuration. No server involved: jsPDF renders entirely in the browser.
//
// The actual download trigger mirrors ui/metricsPanel.ts's CSV export idiom exactly (Blob + a
// synthetic `<a download>` click), just fed jsPDF's own Blob output instead of a hand-built string.

import { jsPDF } from "jspdf";
import { METRIC_GROUPS, METRIC_SPECS, type MetricSpec } from "../metrics/specs";
import type { ImpactEnergyHistogram } from "../metrics/types";
import {
  criticalSpeedRpm,
  effectiveMedia,
  percentCriticalOf,
  rpmOf,
  trueBallCount,
  type CoarseGrainingMode,
  type MediaLike,
  type MillLike,
} from "../params/derived";
import type { ParamsJson } from "../protocol";
import type { ReportAggregator } from "./aggregator";

const PAGE_MARGIN_MM = 14;
const PAGE_WIDTH_MM = 210; // A4 portrait
const PAGE_HEIGHT_MM = 297;
const CONTENT_RIGHT_MM = PAGE_WIDTH_MM - PAGE_MARGIN_MM;
const ACCENT_RGB: [number, number, number] = [58, 130, 168]; // matches ui/metricsPanel.ts's histogram bar color

/** Defensive-cast readout of the parts of `ParamsJson` this report needs, one group at a time --
 * same idiom as ui/metricsPanel.ts's `millOf`/`slurryEnabledOf` and ui/paramsPanel.ts's derived
 * panel. All optional: an older saved params blob predating a field (e.g.
 * `simulation.coarse_graining_mode`, added alongside this feature) falls back to that field's
 * `Params::default()` value so this never throws on a stale params.json. */
interface RunConfig {
  mill: MillLike & { direction: string };
  media: MediaLike;
  maxBalls: number;
  seed: number;
  coarseGrainingMode: CoarseGrainingMode;
  coarseGrainingK: number;
  slurryEnabled: boolean;
  slurryFillFraction: number;
  slurryViscosityPaS: number;
}

function readRunConfig(params: ParamsJson | null): RunConfig | null {
  const mill = params?.mill as
    | { diameter_m?: number; speed_mode?: string; speed_value?: number; direction?: string }
    | undefined;
  const media = params?.media as
    | { ball_diameter_m?: number; fill_fraction?: number; packing_fraction_2d?: number; density_kg_m3?: number }
    | undefined;
  const simulation = params?.simulation as
    | { max_balls?: number; seed?: number; coarse_graining_mode?: string; coarse_graining_k?: number }
    | undefined;
  const slurry = params?.slurry as
    | { enabled?: boolean; fill_fraction?: number; viscosity_pa_s?: number }
    | undefined;

  if (
    !mill ||
    mill.diameter_m === undefined ||
    mill.speed_mode === undefined ||
    mill.speed_value === undefined ||
    !media ||
    media.ball_diameter_m === undefined ||
    media.fill_fraction === undefined ||
    media.packing_fraction_2d === undefined ||
    media.density_kg_m3 === undefined
  ) {
    return null;
  }

  return {
    mill: {
      diameter_m: mill.diameter_m,
      speed_mode: mill.speed_mode,
      speed_value: mill.speed_value,
      direction: mill.direction ?? "counter_clockwise",
    },
    media: {
      ball_diameter_m: media.ball_diameter_m,
      fill_fraction: media.fill_fraction,
      packing_fraction_2d: media.packing_fraction_2d,
      density_kg_m3: media.density_kg_m3,
    },
    maxBalls: simulation?.max_balls ?? 600,
    seed: simulation?.seed ?? 1,
    // Backward-compatible with params saved before `simulation.coarse_graining_mode`/`_k` existed
    // (Params::default(): CoarseGrainingMode::Auto, k = 1.0).
    coarseGrainingMode: (simulation?.coarse_graining_mode as CoarseGrainingMode | undefined) ?? "auto",
    coarseGrainingK: simulation?.coarse_graining_k ?? 1,
    slurryEnabled: slurry?.enabled ?? false,
    slurryFillFraction: slurry?.fill_fraction ?? 0,
    slurryViscosityPaS: slurry?.viscosity_pa_s ?? 0,
  };
}

function directionLabel(direction: string): string {
  return direction === "clockwise" ? "Clockwise" : "Counter-clockwise";
}

/** Draws the impact-energy histogram to an offscreen canvas and returns a PNG data URL, or `null`
 * if there is nothing to draw (no histogram yet, or every bin count is zero) -- callers fall back
 * to a text table in that case. Deliberately larger/plainer than ui/metricsPanel.ts's on-screen
 * sparkline-sized version (that one is tuned for a ~300x80 CSS panel strip); this one is sized for
 * a print figure and drawn on a white background (the PDF page). */
function renderHistogramImage(histogram: ImpactEnergyHistogram | undefined | null): string | null {
  const counts = histogram?.counts_per_s ?? [];
  const edges = histogram?.bin_edges_j ?? [];
  const maxCount = counts.length > 0 ? Math.max(...counts) : 0;
  if (counts.length === 0 || maxCount <= 0) return null;

  const width = 900;
  const height = 340;
  const canvas = document.createElement("canvas");
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext("2d");
  if (!ctx) return null;

  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  const padLeft = 60;
  const padRight = 20;
  const padTop = 40;
  const padBottom = 60;
  const plotW = width - padLeft - padRight;
  const plotH = height - padTop - padBottom;
  const barW = plotW / counts.length;

  ctx.fillStyle = "#3a82a8";
  counts.forEach((count, i) => {
    const h = (count / maxCount) * plotH;
    const x = padLeft + i * barW;
    const y = padTop + (plotH - h);
    ctx.fillRect(x, y, Math.max(1, barW - 2), h);
  });

  ctx.strokeStyle = "#333333";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(padLeft, padTop + plotH);
  ctx.lineTo(padLeft + plotW, padTop + plotH);
  ctx.stroke();

  ctx.fillStyle = "#222222";
  ctx.font = "16px sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(`Peak: ${maxCount.toFixed(2)} impacts/s`, padLeft, 22);

  ctx.font = "12px sans-serif";
  ctx.fillStyle = "#444444";
  for (let i = 0; i < edges.length; i++) {
    const edge = edges[i];
    if (edge === undefined || edge <= 0) continue;
    const log10 = Math.log10(edge);
    if (Math.abs(log10 - Math.round(log10)) > 1e-6) continue; // only label bin edges landing on a decade
    const label = edge.toExponential(0);
    const x = padLeft + i * barW;
    ctx.textAlign = "center";
    ctx.fillText(label, Math.min(Math.max(x, padLeft + 24), width - padRight - 24), height - padBottom + 22);
  }

  ctx.textAlign = "center";
  ctx.font = "13px sans-serif";
  ctx.fillText("Impact energy (J)", padLeft + plotW / 2, height - 8);

  ctx.save();
  ctx.translate(18, padTop + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText("Impacts / s", 0, 0);
  ctx.restore();

  return canvas.toDataURL("image/png");
}

/** Formats one `MetricStats` value with `spec.format` if given (falls back to 2-decimal fixed),
 * or "-" for `null` (no data). */
function fmt(spec: MetricSpec, value: number | null): string {
  if (value === null) return "-";
  return spec.format ? spec.format(value) : value.toFixed(2);
}

interface Cursor {
  y: number;
}

/** Starts a new page and resets `cursor.y` to the top margin if `needed` mm of vertical space
 * would overflow the page -- the only page-break strategy this report needs, since every element
 * (a table row, a heading + its content) is drawn top-to-bottom in one pass. */
function ensureSpace(doc: jsPDF, cursor: Cursor, needed: number): void {
  if (cursor.y + needed > PAGE_HEIGHT_MM - PAGE_MARGIN_MM) {
    doc.addPage();
    cursor.y = PAGE_MARGIN_MM;
  }
}

function drawKeyValueRow(doc: jsPDF, cursor: Cursor, label: string, value: string): void {
  ensureSpace(doc, cursor, 6);
  doc.setFont("helvetica", "normal");
  doc.setFontSize(10);
  doc.setTextColor(90, 90, 90);
  doc.text(label, PAGE_MARGIN_MM, cursor.y);
  doc.setTextColor(20, 20, 20);
  doc.text(value, CONTENT_RIGHT_MM, cursor.y, { align: "right" });
  cursor.y += 5.5;
}

export interface ReportInput {
  /** Current params (may be `null` before the worker's first "ready" message -- the report is
   * still generated, just without a run-configuration section). */
  params: ParamsJson | null;
  /** The report's aggregation window and per-metric running statistics. */
  aggregator: ReportAggregator;
  /** The most recent frame's impact-energy histogram (already EMA-smoothed by mill-core), or
   * `null`/`undefined` if none has arrived yet. */
  impactEnergyHistogram: ImpactEnergyHistogram | undefined | null;
  /** Sim time at generation, used only for the header line (the filename is built separately by
   * the caller, mirroring the CSV export's `lastSimTime.toFixed(1)` idiom). */
  simTimeAtGeneration: number;
}

/** The headline results block's metric ids, in display order -- a subset of
 * metrics/specs.ts's `KEY_METRIC_IDS` (that list also includes charge-geometry angles and total
 * kinetic energy, which read more like debug/health values than "what a paper's abstract leads
 * with"). `mixing_index` is appended conditionally (slurry-only). */
const HEADLINE_METRIC_IDS = ["power_draw", "torque", "collision_rate", "dissipated_power"] as const;

/** Builds the full report as a jsPDF document. Pure with respect to global state except for
 * `document.createElement("canvas")` in `renderHistogramImage` (browser-only; this module is never
 * imported by the aggregator's unit tests). */
export function buildReportPdf(input: ReportInput): jsPDF {
  const { params, aggregator, impactEnergyHistogram, simTimeAtGeneration } = input;
  const doc = new jsPDF({ unit: "mm", format: "a4" });
  const cursor: Cursor = { y: PAGE_MARGIN_MM };

  // --- Header ------------------------------------------------------------------------------
  doc.setFont("helvetica", "bold");
  doc.setFontSize(18);
  doc.setTextColor(20, 20, 20);
  doc.text("MillDynamics3 Simulation Report", PAGE_MARGIN_MM, cursor.y + 4);
  cursor.y += 11;

  doc.setFont("helvetica", "normal");
  doc.setFontSize(9.5);
  doc.setTextColor(110, 110, 110);
  doc.text(`Generated: ${new Date().toISOString()}`, PAGE_MARGIN_MM, cursor.y);
  cursor.y += 5;

  const hasWindow = aggregator.startTime !== null && aggregator.endTime !== null;
  const windowText = hasWindow
    ? `Aggregation window: t = ${aggregator.startTime!.toFixed(1)} s -> ${aggregator.endTime!.toFixed(1)} s ` +
      `(duration ${(aggregator.endTime! - aggregator.startTime!).toFixed(1)} s, ${aggregator.samples} samples)`
    : "Aggregation window: no data aggregated yet";
  doc.text(windowText, PAGE_MARGIN_MM, cursor.y);
  cursor.y += 5;
  doc.text(`Sim time at generation: ${simTimeAtGeneration.toFixed(1)} s`, PAGE_MARGIN_MM, cursor.y);
  cursor.y += 8;

  // --- Run configuration ---------------------------------------------------------------------
  doc.setFont("helvetica", "bold");
  doc.setFontSize(12);
  doc.setTextColor(20, 20, 20);
  doc.text("Run configuration", PAGE_MARGIN_MM, cursor.y);
  cursor.y += 6;

  const config = readRunConfig(params);
  if (!config) {
    doc.setFont("helvetica", "italic");
    doc.setFontSize(10);
    doc.setTextColor(110, 110, 110);
    doc.text("Parameters unavailable (report generated before the simulation finished loading).", PAGE_MARGIN_MM, cursor.y);
    cursor.y += 8;
  } else {
    const eff = effectiveMedia(config.mill.diameter_m, config.media, config.maxBalls, config.coarseGrainingMode, config.coarseGrainingK);
    const coarseGrainingActive = eff.scaleFactor > 1;
    const nc = criticalSpeedRpm(config.mill.diameter_m);

    drawKeyValueRow(doc, cursor, "Drum diameter", `${(config.mill.diameter_m * 1000).toFixed(1)} mm`);
    drawKeyValueRow(
      doc,
      cursor,
      "Rotation speed",
      `${rpmOf(config.mill).toFixed(1)} rpm (${percentCriticalOf(config.mill).toFixed(0)}% Nc, Nc = ${nc.toFixed(1)} rpm)`,
    );
    drawKeyValueRow(doc, cursor, "Direction", directionLabel(config.mill.direction));
    drawKeyValueRow(doc, cursor, "Ball diameter (true)", `${(config.media.ball_diameter_m * 1000).toFixed(2)} mm`);
    drawKeyValueRow(doc, cursor, "Ball diameter (effective)", `${(eff.diameterM * 1000).toFixed(2)} mm`);
    drawKeyValueRow(doc, cursor, "Fill fraction (J)", config.media.fill_fraction.toFixed(3));
    drawKeyValueRow(doc, cursor, "Media density", `${config.media.density_kg_m3.toFixed(0)} kg/m^3`);
    drawKeyValueRow(doc, cursor, "True ball count (N_true)", trueBallCount(config.mill.diameter_m, config.media).toFixed(0));
    drawKeyValueRow(doc, cursor, "Simulated ball count (N_sim)", String(eff.ballCount));
    drawKeyValueRow(
      doc,
      cursor,
      "Coarse-graining",
      coarseGrainingActive ? `ON (k = ${eff.scaleFactor.toFixed(3)})` : "OFF (k = 1, true size)",
    );
    if (config.slurryEnabled) {
      drawKeyValueRow(doc, cursor, "Slurry", "Enabled");
      drawKeyValueRow(doc, cursor, "Slurry fill fraction", config.slurryFillFraction.toFixed(3));
      drawKeyValueRow(doc, cursor, "Slurry viscosity", `${config.slurryViscosityPaS.toFixed(2)} Pa*s`);
    } else {
      drawKeyValueRow(doc, cursor, "Slurry", "Disabled");
    }
    drawKeyValueRow(doc, cursor, "Random seed", String(config.seed));
    cursor.y += 2;
  }

  // --- Headline results (visually distinguished box) ------------------------------------------
  const headlineIds = config?.slurryEnabled ? [...HEADLINE_METRIC_IDS, "mixing_index"] : [...HEADLINE_METRIC_IDS];
  const headlineSpecs = headlineIds
    .map((id) => METRIC_SPECS.find((s) => s.id === id))
    .filter((s): s is MetricSpec => s !== undefined);

  const boxLineHeight = 7;
  const boxHeight = 12 + headlineSpecs.length * boxLineHeight;
  ensureSpace(doc, cursor, boxHeight + 6);

  doc.setFillColor(230, 240, 247);
  doc.setDrawColor(...ACCENT_RGB);
  doc.setLineWidth(0.6);
  doc.roundedRect(PAGE_MARGIN_MM, cursor.y, CONTENT_RIGHT_MM - PAGE_MARGIN_MM, boxHeight, 2.5, 2.5, "FD");

  let boxY = cursor.y + 8;
  doc.setFont("helvetica", "bold");
  doc.setFontSize(13);
  doc.setTextColor(20, 20, 20);
  doc.text("Headline results", PAGE_MARGIN_MM + 5, boxY);
  boxY += boxLineHeight;

  doc.setFontSize(11.5);
  for (const spec of headlineSpecs) {
    const stats = aggregator.statsFor(spec.id);
    const unitSuffix = spec.unit ? ` ${spec.unit}` : "";
    const valueText = stats ? `${fmt(spec, stats.mean)} +/- ${fmt(spec, stats.std)}${unitSuffix}` : "No data";
    doc.setFont("helvetica", "normal");
    doc.text(spec.label, PAGE_MARGIN_MM + 5, boxY);
    doc.setFont("helvetica", "bold");
    doc.text(valueText, CONTENT_RIGHT_MM - 5, boxY, { align: "right" });
    boxY += boxLineHeight;
  }
  cursor.y += boxHeight + 8;

  // --- Full metrics table ----------------------------------------------------------------------
  doc.setFont("helvetica", "bold");
  doc.setFontSize(12);
  doc.setTextColor(20, 20, 20);
  ensureSpace(doc, cursor, 10);
  doc.text("All tracked metrics (mean / std / min / max over the aggregation window)", PAGE_MARGIN_MM, cursor.y);
  cursor.y += 7;

  const colMetricX = PAGE_MARGIN_MM;
  const colMeanX = 110;
  const colStdX = 140;
  const colMinX = 165;
  const colMaxX = 190;

  function drawTableHeader(): void {
    ensureSpace(doc, cursor, 6);
    doc.setFont("courier", "bold");
    doc.setFontSize(8.5);
    doc.setTextColor(90, 90, 90);
    doc.text("Metric", colMetricX, cursor.y);
    doc.text("Mean", colMeanX, cursor.y, { align: "right" });
    doc.text("Std", colStdX, cursor.y, { align: "right" });
    doc.text("Min", colMinX, cursor.y, { align: "right" });
    doc.text("Max", colMaxX, cursor.y, { align: "right" });
    cursor.y += 4.5;
    doc.setDrawColor(180, 180, 180);
    doc.setLineWidth(0.2);
    doc.line(PAGE_MARGIN_MM, cursor.y - 3, CONTENT_RIGHT_MM, cursor.y - 3);
  }

  drawTableHeader();
  for (const group of METRIC_GROUPS) {
    const specsInGroup = METRIC_SPECS.filter((s) => s.group === group);
    if (specsInGroup.length === 0) continue;

    ensureSpace(doc, cursor, 6);
    doc.setFont("helvetica", "bold");
    doc.setFontSize(9);
    doc.setTextColor(58, 130, 168);
    doc.text(group, colMetricX, cursor.y);
    cursor.y += 4.5;

    doc.setFont("courier", "normal");
    doc.setFontSize(8.5);
    doc.setTextColor(30, 30, 30);
    for (const spec of specsInGroup) {
      ensureSpace(doc, cursor, 4.5);
      const stats = aggregator.statsFor(spec.id);
      const label = spec.unit ? `${spec.label} (${spec.unit})` : spec.label;
      doc.text(label.length > 44 ? `${label.slice(0, 43)}...` : label, colMetricX, cursor.y);
      doc.text(fmt(spec, stats?.mean ?? null), colMeanX, cursor.y, { align: "right" });
      doc.text(fmt(spec, stats?.std ?? null), colStdX, cursor.y, { align: "right" });
      doc.text(fmt(spec, stats?.min ?? null), colMinX, cursor.y, { align: "right" });
      doc.text(fmt(spec, stats?.max ?? null), colMaxX, cursor.y, { align: "right" });
      cursor.y += 4.5;
    }
    cursor.y += 2;
  }

  // --- Impact energy histogram -----------------------------------------------------------------
  ensureSpace(doc, cursor, 14);
  doc.setFont("helvetica", "bold");
  doc.setFontSize(12);
  doc.setTextColor(20, 20, 20);
  doc.text("Impact energy distribution", PAGE_MARGIN_MM, cursor.y);
  cursor.y += 7;

  const imageDataUrl = renderHistogramImage(impactEnergyHistogram);
  if (imageDataUrl) {
    const imgWidthMm = CONTENT_RIGHT_MM - PAGE_MARGIN_MM;
    const imgHeightMm = imgWidthMm * (340 / 900);
    ensureSpace(doc, cursor, imgHeightMm + 4);
    doc.addImage(imageDataUrl, "PNG", PAGE_MARGIN_MM, cursor.y, imgWidthMm, imgHeightMm);
    cursor.y += imgHeightMm + 6;
  } else {
    const counts = impactEnergyHistogram?.counts_per_s ?? [];
    const edges = impactEnergyHistogram?.bin_edges_j ?? [];
    if (counts.length === 0) {
      doc.setFont("helvetica", "italic");
      doc.setFontSize(10);
      doc.setTextColor(110, 110, 110);
      doc.text("No impact-energy data available.", PAGE_MARGIN_MM, cursor.y);
      cursor.y += 6;
    } else {
      doc.setFont("courier", "bold");
      doc.setFontSize(8.5);
      doc.setTextColor(90, 90, 90);
      doc.text("Bin range (J)", colMetricX, cursor.y);
      doc.text("Count (1/s)", colMeanX, cursor.y, { align: "right" });
      cursor.y += 4.5;
      doc.setFont("courier", "normal");
      doc.setTextColor(30, 30, 30);
      for (let i = 0; i < counts.length; i++) {
        ensureSpace(doc, cursor, 4.5);
        const lo = edges[i];
        const hi = edges[i + 1];
        const rangeText = lo !== undefined && hi !== undefined ? `${lo.toExponential(2)} - ${hi.toExponential(2)}` : `bin ${i}`;
        doc.text(rangeText, colMetricX, cursor.y);
        doc.text(counts[i]!.toFixed(2), colMeanX, cursor.y, { align: "right" });
        cursor.y += 4.5;
      }
    }
  }

  return doc;
}

/** Downloads `doc` as a PDF file, mirroring ui/metricsPanel.ts's CSV export idiom exactly: a Blob
 * from the document's own output, an object URL, a synthetic `<a download>` click, then revoke. */
export function downloadReportPdf(doc: jsPDF, simTimeAtGeneration: number): void {
  const blob = doc.output("blob");
  const url = URL.createObjectURL(blob);
  try {
    const a = document.createElement("a");
    a.href = url;
    a.download = `milldynamics-report-${simTimeAtGeneration.toFixed(1)}s.pdf`;
    a.click();
  } finally {
    URL.revokeObjectURL(url);
  }
}
