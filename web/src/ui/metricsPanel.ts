// Right-side metrics panel (docs/PLAN.md §4.4): a top "Key results" block (a short, defended
// shortlist -- see the comment on `KEY_METRIC_IDS` in metrics/specs.ts for what's in it and why)
// followed by every other `MetricSpec` grouped into collapsible `<details>` (closed by default,
// mirroring ui/paramsPanel.ts's group pattern), a rolling history feeding a few sparklines plus a
// CSV export, and a small impact-energy histogram chart. All the data plumbing (spec table,
// history ring buffer, sparkline renderer, `Metrics` type) lives elsewhere; this file only
// builds/updates the DOM.

import {
  coarseGrainingReason,
  KEY_METRIC_IDS,
  METRIC_GROUPS,
  METRIC_SPECS,
  type MetricContext,
  type MetricGroupName,
  type MetricSpec,
} from "../metrics/specs";
import { MetricsHistory, type CsvColumn } from "../metrics/history";
import type { ImpactEnergyHistogram } from "../metrics/types";
import { criticalSpeedRpm, percentCriticalOf, rpmOf } from "../params/derived";
import type { ParamsJson } from "../protocol";
import { ReportAggregator } from "../report/aggregator";
import { buildReportPdf, downloadReportPdf, energyCumulativeFraction } from "../report/pdf";
import type { AppState } from "../state";
import { createHelpIcon } from "./helpIcon";
import { drawSparkline } from "./sparkline";

export interface MetricsPanel {
  update(state: AppState, fps: number, nowMs: number): void;
  reset(): void;
}

/** Same defensive-cast pattern as ui/hud.ts's `millOf`: pulls the bits of `params.mill` this
 * panel needs to derive rpm/%Nc/critical speed, or `null` if any are missing (e.g. before the
 * worker's "ready" message arrives). */
function millOf(params: ParamsJson | null): { diameter_m: number; speed_mode: string; speed_value: number } | null {
  const mill = params?.mill as { diameter_m?: number; speed_mode?: string; speed_value?: number } | undefined;
  if (!mill || mill.diameter_m === undefined || mill.speed_mode === undefined || mill.speed_value === undefined) {
    return null;
  }
  return { diameter_m: mill.diameter_m, speed_mode: mill.speed_mode, speed_value: mill.speed_value };
}

/** `params.slurry.enabled`, or `null` before the first params message arrives -- feeds the
 * Mixing index key row's `unavailable` reason (specs.ts) so it disappears cleanly rather than
 * showing a permanent "-" once slurry is switched off. */
function slurryEnabledOf(params: ParamsJson | null): boolean | null {
  const slurry = params?.slurry as { enabled?: boolean } | undefined;
  return slurry?.enabled ?? null;
}

const HIST_CSS_WIDTH = 300;
const HIST_CSS_HEIGHT = 80;

/** (Re)sizes the impact-histogram canvas's backing store for its fixed CSS footprint, mirroring
 * ui/sparkline.ts's `ensureSized` technique (only resize when the DPR-scaled size actually
 * changed, then reset the transform so drawing code can work in CSS pixels). */
function ensureHistSized(canvas: HTMLCanvasElement): CanvasRenderingContext2D | null {
  const dpr = window.devicePixelRatio || 1;
  const wantW = Math.round(HIST_CSS_WIDTH * dpr);
  const wantH = Math.round(HIST_CSS_HEIGHT * dpr);
  if (canvas.width !== wantW || canvas.height !== wantH) {
    canvas.width = wantW;
    canvas.height = wantH;
    canvas.style.width = `${HIST_CSS_WIDTH}px`;
    canvas.style.height = `${HIST_CSS_HEIGHT}px`;
  }
  const ctx = canvas.getContext("2d");
  if (!ctx) return null;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return ctx;
}

/** Bar chart of `Metrics.impact_energy_histogram`: equal-width bars (the underlying bins are
 * already log-spaced, so equal screen width per bin is correct), a peak-count label, and decade
 * tick labels on the x-axis where a bin edge lands on (close to) a power of ten. */
function drawImpactHistogram(canvas: HTMLCanvasElement, histogram: ImpactEnergyHistogram | undefined | null): void {
  const ctx = ensureHistSized(canvas);
  if (!ctx) return;
  ctx.clearRect(0, 0, HIST_CSS_WIDTH, HIST_CSS_HEIGHT);

  const counts = histogram?.counts_per_s ?? [];
  const edges = histogram?.bin_edges_j ?? [];
  const maxCount = counts.length > 0 ? Math.max(...counts) : 0;
  if (counts.length === 0 || maxCount <= 0) return;

  const padLeft = 4;
  const padRight = 4;
  const padTop = 14;
  const padBottom = 12;
  const plotW = HIST_CSS_WIDTH - padLeft - padRight;
  const plotH = HIST_CSS_HEIGHT - padTop - padBottom;
  const barW = plotW / counts.length;

  ctx.fillStyle = "#3a82a8";
  counts.forEach((count, i) => {
    const h = (count / maxCount) * plotH;
    const x = padLeft + i * barW;
    const y = padTop + (plotH - h);
    ctx.fillRect(x, y, Math.max(1, barW - 1), h);
  });

  // Energy-weighted cumulative fraction (0 at the bottom, 1 at the top of the plot area).
  const cum = energyCumulativeFraction(edges, counts);
  if (cum) {
    ctx.strokeStyle = "#e0803a";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    cum.forEach((f, i) => {
      const x = padLeft + (i + 0.5) * barW;
      const y = padTop + plotH * (1 - f);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();
  }

  ctx.fillStyle = "#9aa3b0";
  ctx.font = "9px system-ui, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(`${maxCount.toFixed(1)}/s`, padLeft, 10);
  if (cum) {
    ctx.fillStyle = "#e0803a";
    ctx.textAlign = "right";
    ctx.fillText("cum. energy", HIST_CSS_WIDTH - padRight, 10);
  }

  // Decade tick labels: anchored `center` except within one label-width of either canvas edge,
  // where centering would draw half the text past the edge and clip it (this was most visible on
  // the rightmost label, e.g. "1e+1" losing its last character -- reported by an external
  // review). Switching anchor near the edges keeps every label fully inside the canvas instead.
  for (let i = 0; i < edges.length; i++) {
    const edge = edges[i];
    if (edge === undefined || edge <= 0) continue;
    const log10 = Math.log10(edge);
    if (Math.abs(log10 - Math.round(log10)) > 1e-6) continue;
    const label = edge.toExponential(0);
    const x = padLeft + i * barW;
    const halfWidth = ctx.measureText(label).width / 2;
    if (x + halfWidth > HIST_CSS_WIDTH) {
      ctx.textAlign = "right";
      ctx.fillText(label, HIST_CSS_WIDTH, HIST_CSS_HEIGHT - 2);
    } else if (x - halfWidth < 0) {
      ctx.textAlign = "left";
      ctx.fillText(label, 0, HIST_CSS_HEIGHT - 2);
    } else {
      ctx.textAlign = "center";
      ctx.fillText(label, x, HIST_CSS_HEIGHT - 2);
    }
  }
}

// Warning-highlight thresholds (raw metric values, not formatted strings): a row is flagged
// `is-warn` once its value exceeds the threshold below. All ids below use the same "value is
// greater than threshold" test -- coupling_clamp_hits's "warn on any hit" semantics is just the
// threshold-0 case of that same rule, not a different rule. The Solver-health key-block roll-up
// (below) reads this same map, so there is exactly one place these four diagnostics are defined
// as "alarms" rather than "continuously-read numbers".
const COUPLING_CLAMP_HITS_WARN_THRESHOLD = 0;
const SUBSTEP_DISPLACEMENT_WARN_THRESHOLD = 0.5;
// Mean, not max: the single worst particle spikes to tens of percent at any free surface even on a
// healthy run (measured: max ~20-35% at rest, ~65% at 100 rpm, mean 2% / 6% -- docs/METRICS.md), so a
// max-based alarm fires constantly. The mean tracks whether the whole fluid is actually compressed.
const MEAN_COMPRESSION_ERROR_WARN_THRESHOLD = 0.1;
const MAX_BALL_OVERLAP_WARN_THRESHOLD = 0.5;

const WARN_THRESHOLDS = new Map<string, number>([
  ["coupling_clamp_hits", COUPLING_CLAMP_HITS_WARN_THRESHOLD],
  ["substep_displacement", SUBSTEP_DISPLACEMENT_WARN_THRESHOLD], // confirmed against metrics/specs.ts
  ["mean_compression_error", MEAN_COMPRESSION_ERROR_WARN_THRESHOLD],
  ["max_ball_overlap", MAX_BALL_OVERLAP_WARN_THRESHOLD],
]);

const DEFAULT_REPORT_WARMUP_S = 2;
const DEFAULT_REPORT_END_S = 7;
const GENERATE_LABEL = "Generate report now";

const GROUPS_STORAGE_KEY = "milldynamics.metricsGroups";

// All five groups start collapsed: the key-results block above them already carries the numbers
// that matter at a glance (see specs.ts's `KEY_METRIC_IDS`), so everything else is opt-in detail.
// Mirrors ui/paramsPanel.ts's per-group open/closed persistence pattern one-for-one, under its own
// storage key (paramsPanel.ts's is "milldynamics.paramsGroups", a separate panel/preference).
const DEFAULT_GROUP_OPEN: Record<MetricGroupName, boolean> = {
  Drum: false,
  Media: false,
  Grinding: false,
  Slurry: false,
  Solver: false,
};

function readGroupOpenPrefs(): Record<string, boolean> {
  try {
    const raw = localStorage.getItem(GROUPS_STORAGE_KEY);
    if (raw === null) return { ...DEFAULT_GROUP_OPEN };
    const parsed = JSON.parse(raw) as Record<string, unknown>;
    const result: Record<string, boolean> = { ...DEFAULT_GROUP_OPEN };
    for (const key of Object.keys(result)) {
      if (typeof parsed[key] === "boolean") result[key] = parsed[key] as boolean;
    }
    return result;
  } catch {
    return { ...DEFAULT_GROUP_OPEN };
  }
}

function writeGroupOpenPrefs(prefs: Record<string, boolean>): void {
  try {
    localStorage.setItem(GROUPS_STORAGE_KEY, JSON.stringify(prefs));
  } catch {
    // ignore (private browsing / disabled storage)
  }
}

/**
 * `onAutoReportPause` is called synchronously (before the PDF is built) when the armed "Report end
 * time" is reached during a running sim -- it is main.ts's job to actually pause the worker
 * (`send({type: "pause"})`) and reflect that in the toolbar's play/pause icon, mirroring the
 * toolbar's own pause path; this module only owns the aggregation/report-building/DOM side of the
 * feature (see report/aggregator.ts, report/pdf.ts).
 */
export function createMetricsPanel(
  container: HTMLElement,
  onAutoReportPause: () => void,
  getSnapshotPng?: () => string | null,
): MetricsPanel {
  const groupOpenPrefs = readGroupOpenPrefs();

  const header = document.createElement("div");
  header.className = "metrics-panel-header";
  const title = document.createElement("h2");
  title.textContent = "Metrics";
  const exportButton = document.createElement("button");
  exportButton.type = "button";
  exportButton.className = "metrics-export";
  exportButton.textContent = "Export CSV";
  exportButton.disabled = true;
  header.appendChild(title);
  header.appendChild(exportButton);
  container.appendChild(header);

  // --- Report (report/aggregator.ts + report/pdf.ts) ----------------------------------------
  // A small subsection right below the header/CSV button: an optional "Report end time" that
  // arms an auto-pause-and-download, plus a manual "Generate report now" that always works off
  // whatever the aggregator has collected so far. See docs on `onAutoReportPause` above and
  // `reportArmed`/`reportEndTimeS` below for the arm/disarm rules.
  const reportSection = document.createElement("section");
  reportSection.className = "metrics-report";
  const reportHeading = document.createElement("h3");
  reportHeading.textContent = "Report";
  reportSection.appendChild(reportHeading);

  function makeNumberRow(id: string, labelText: string, value: string): HTMLInputElement {
    const row = document.createElement("div");
    row.className = "metrics-report-row";
    const label = document.createElement("label");
    label.textContent = labelText;
    label.htmlFor = id;
    const input = document.createElement("input");
    input.type = "number";
    input.id = id;
    input.className = "metrics-report-input";
    input.min = "0";
    input.step = "0.5";
    input.value = value;
    row.appendChild(label);
    row.appendChild(input);
    reportSection.appendChild(row);
    return input;
  }
  const reportWarmupInput = makeNumberRow("metrics-report-warmup", "Warm-up (s)", String(DEFAULT_REPORT_WARMUP_S));
  const reportEndTimeInput = makeNumberRow("metrics-report-end-time", "Report end (s)", String(DEFAULT_REPORT_END_S));
  reportEndTimeInput.placeholder = "off";

  const autoRow = document.createElement("div");
  autoRow.className = "metrics-report-row";
  const autoLabel = document.createElement("label");
  autoLabel.htmlFor = "metrics-report-auto";
  autoLabel.textContent = "Auto pause & download at end";
  const autoCheckbox = document.createElement("input");
  autoCheckbox.type = "checkbox";
  autoCheckbox.id = "metrics-report-auto";
  autoCheckbox.checked = false;
  autoRow.appendChild(autoLabel);
  autoRow.appendChild(autoCheckbox);
  reportSection.appendChild(autoRow);

  const revNote = document.createElement("div");
  revNote.className = "metrics-note";
  reportSection.appendChild(revNote);

  const generateReportButton = document.createElement("button");
  generateReportButton.type = "button";
  generateReportButton.className = "metrics-export metrics-report-generate";
  generateReportButton.textContent = GENERATE_LABEL;
  reportSection.appendChild(generateReportButton);

  const reportNote = document.createElement("div");
  reportNote.className = "metrics-note";
  reportNote.textContent =
    "PDF summary (stats, time series, snapshot, impact spectrum) over the window [warm-up, end]. " +
    "Samples outside the window are ignored. Generating before the end uses [warm-up, now]. " +
    "End = 0 or empty disables the end.";
  reportSection.appendChild(reportNote);

  container.appendChild(reportSection);

  const reportAggregator = new ReportAggregator();
  let reportWarmupS = DEFAULT_REPORT_WARMUP_S;
  let reportEndTimeS = DEFAULT_REPORT_END_S;
  reportAggregator.setWindow(reportWarmupS, reportEndTimeS);
  let reportArmed = true;
  let windowComplete = false;
  let lastAggregatedMetrics: AppState["metrics"] = null;
  let latestState: AppState | null = null;
  let latestRpm: number | null = null;

  function refreshReportUi(): void {
    let note: string;
    if (reportEndTimeS <= 0) note = "No end time set.";
    else if (latestRpm === null) note = "";
    else {
      const revs = (latestRpm * Math.max(0, reportEndTimeS - reportWarmupS)) / 60;
      note = `≈ ${revs.toFixed(1)} revolutions in the window`;
    }
    if (revNote.textContent !== note) revNote.textContent = note;
    const label = windowComplete ? "Report window complete – download report" : GENERATE_LABEL;
    if (generateReportButton.textContent !== label) generateReportButton.textContent = label;
    generateReportButton.classList.toggle("is-complete", windowComplete);
  }

  function applyWindowInputs(): void {
    const w = Number(reportWarmupInput.value);
    const e = Number(reportEndTimeInput.value);
    reportWarmupS = Number.isFinite(w) && w > 0 ? w : 0;
    reportEndTimeS = Number.isFinite(e) && e > 0 ? e : 0;
    reportAggregator.setWindow(reportWarmupS, reportEndTimeS);
    // Any window edit re-arms the end trigger (raising the end after a completed window resumes it).
    reportArmed = true;
    windowComplete = latestState !== null && reportAggregator.isPastEnd(latestState.simTime);
    refreshReportUi();
  }

  function generateReport(state: AppState): void {
    const doc = buildReportPdf({
      params: state.params,
      aggregator: reportAggregator,
      impactEnergyHistogram: reportAggregator.windowHistogram() ?? state.metrics?.impact_energy_histogram,
      simTimeAtGeneration: state.simTime,
      snapshotPng: getSnapshotPng ? getSnapshotPng() : null,
    });
    downloadReportPdf(doc, state.simTime);
  }

  reportWarmupInput.addEventListener("input", applyWindowInputs);
  reportEndTimeInput.addEventListener("input", applyWindowInputs);

  generateReportButton.addEventListener("click", () => {
    if (latestState) generateReport(latestState);
  });

  const valueEls = new Map<string, HTMLSpanElement>();
  const sparklineEls = new Map<string, HTMLCanvasElement>();
  const rowEls = new Map<string, HTMLElement>();
  const groupDetailsEls = new Map<MetricGroupName, HTMLDetailsElement>();
  const groupChipEls = new Map<MetricGroupName, HTMLSpanElement>();
  let histogramCanvas: HTMLCanvasElement | null = null;
  let histogramWrapEl: HTMLElement | null = null;
  let coarseGrainingWarnEl: HTMLElement | null = null;
  let lastWarnedIds: string[] = [];

  /** Builds one `.metric-row` for `spec`: label + value spans, and an optional sparkline canvas,
   * registering all three into the maps above. Shared by the key-results block and the
   * collapsible groups below so a row's DOM is only ever built in one place -- each metric still
   * renders exactly once, just in one of two locations depending on `spec.key`. */
  function buildRow(spec: MetricSpec, extraClass?: string): HTMLElement {
    const row = document.createElement("div");
    row.className = extraClass ? `metric-row ${extraClass}` : "metric-row";

    const label = document.createElement("span");
    label.className = "metric-label";
    label.append(spec.unit ? `${spec.label} (${spec.unit})` : spec.label);
    const helpIcon = createHelpIcon(spec.help);
    if (helpIcon) label.appendChild(helpIcon);

    const value = document.createElement("span");
    value.className = "metric-value";
    value.textContent = "-";

    row.appendChild(label);
    row.appendChild(value);

    if (spec.sparkline) {
      const canvas = document.createElement("canvas");
      canvas.className = "metric-sparkline";
      row.appendChild(canvas);
      sparklineEls.set(spec.id, canvas);
    }

    valueEls.set(spec.id, value);
    rowEls.set(spec.id, row);
    return row;
  }

  // --- Key results ------------------------------------------------------------------------------
  const keySection = document.createElement("section");
  keySection.className = "metrics-key";
  const keyHeading = document.createElement("h3");
  keyHeading.textContent = "Key results";
  keySection.appendChild(keyHeading);

  for (const id of KEY_METRIC_IDS) {
    const spec = METRIC_SPECS.find((s) => s.id === id);
    if (!spec) continue; // guarded by specs.test.ts: every KEY_METRIC_IDS entry must have a spec
    keySection.appendChild(buildRow(spec, "is-key"));
  }

  // Solver-health roll-up: one line summarising the five WARN_THRESHOLDS diagnostics (below)
  // instead of making the reader open three separately-collapsed groups to find them. A button
  // (not a plain row) so clicking it can open whichever group(s) currently hold a warned row.
  const healthRow = document.createElement("button");
  healthRow.type = "button";
  healthRow.className = "metric-row is-key is-health metric-health-link";
  const healthLabel = document.createElement("span");
  healthLabel.className = "metric-label";
  healthLabel.append("Solver health");
  const healthHelp = createHelpIcon(
    "Roll-up of the solver diagnostics (ball overlap, sub-step displacement, compression error, coupling clamp hits). OK means none is over its warning threshold; click to open the offending group.",
  );
  if (healthHelp) healthLabel.appendChild(healthHelp);
  const healthValue = document.createElement("span");
  healthValue.className = "metric-value";
  healthValue.textContent = "OK";
  healthRow.appendChild(healthLabel);
  healthRow.appendChild(healthValue);
  healthRow.addEventListener("click", () => {
    let firstDetails: HTMLDetailsElement | null = null;
    for (const id of lastWarnedIds) {
      const spec = METRIC_SPECS.find((s) => s.id === id);
      if (!spec) continue;
      const details = groupDetailsEls.get(spec.group);
      if (!details) continue;
      details.open = true;
      groupOpenPrefs[spec.group] = true;
      firstDetails ??= details;
    }
    if (firstDetails) {
      writeGroupOpenPrefs(groupOpenPrefs);
      firstDetails.scrollIntoView({ block: "center" });
    }
  });
  keySection.appendChild(healthRow);

  container.appendChild(keySection);

  // --- Collapsible groups -------------------------------------------------------------------------
  for (const group of METRIC_GROUPS) {
    const details = document.createElement("details");
    details.className = "metrics-group";
    details.open = groupOpenPrefs[group] ?? false;
    details.addEventListener("toggle", () => {
      groupOpenPrefs[group] = details.open;
      writeGroupOpenPrefs(groupOpenPrefs);
    });
    groupDetailsEls.set(group, details);

    const summary = document.createElement("summary");
    summary.textContent = group;
    if (group === "Grinding") {
      // Keeps the fact that coarse-graining is active (and by how much) discoverable while the
      // group is collapsed -- otherwise it's only visible after opening the group and reading the
      // banner below, same reasoning as hud.ts's media-substitution badge.
      const chip = document.createElement("span");
      chip.className = "metrics-summary-chip";
      chip.hidden = true;
      summary.appendChild(chip);
      groupChipEls.set(group, chip);
    }
    details.appendChild(summary);

    for (const spec of METRIC_SPECS.filter((s) => s.group === group && !s.key)) {
      details.appendChild(buildRow(spec));
    }

    if (group === "Grinding") {
      const warnEl = document.createElement("div");
      warnEl.className = "metrics-warn";
      warnEl.hidden = true;
      details.appendChild(warnEl);
      coarseGrainingWarnEl = warnEl;

      const histWrap = document.createElement("div");
      histWrap.className = "impact-histogram-wrap";
      const histLabel = document.createElement("div");
      histLabel.className = "impact-histogram-label";
      histLabel.textContent = "Impact energy distribution";
      const canvas = document.createElement("canvas");
      canvas.className = "impact-histogram";
      histWrap.appendChild(histLabel);
      histWrap.appendChild(canvas);
      details.appendChild(histWrap);
      histogramCanvas = canvas;
      histogramWrapEl = histWrap;
    }

    container.appendChild(details);
  }

  const history = new MetricsHistory(600, 0.1);
  const csvColumns: CsvColumn[] = METRIC_SPECS.filter((s) => s.csv).map((s) => ({ id: s.id, label: s.label, unit: s.unit }));

  let lastSimTime = 0;

  exportButton.addEventListener("click", () => {
    const csv = history.toCsv(csvColumns);
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    try {
      const a = document.createElement("a");
      a.href = url;
      a.download = `milldynamics-metrics-${lastSimTime.toFixed(1)}s.csv`;
      a.click();
    } finally {
      URL.revokeObjectURL(url);
    }
  });

  let lastRenderMs = -Infinity;

  function update(state: AppState, fps: number, nowMs: number): void {
    const mill = millOf(state.params);
    const ctx: MetricContext = {
      metrics: state.metrics,
      simTime: state.simTime,
      fps,
      achievedTimeScale: state.achievedTimeScale,
      subStepsPerSecondAchieved: state.subStepsPerSecondAchieved,
      subStepsPerSecondRequired: state.subStepsPerSecondRequired,
      rpm: mill ? rpmOf(mill) : null,
      percentCritical: mill ? percentCriticalOf(mill) : null,
      criticalSpeedRpm: mill ? criticalSpeedRpm(mill.diameter_m) : null,
      slurryEnabled: slurryEnabledOf(state.params),
    };

    const values: Record<string, number | null> = {};
    const reasons: Record<string, string | null> = {};
    for (const spec of METRIC_SPECS) {
      const raw = spec.value(ctx);
      const reason = spec.unavailable ? spec.unavailable(raw, ctx) : null;
      reasons[spec.id] = reason;
      values[spec.id] = reason !== null ? null : raw;
    }

    // Sampling must happen every frame (MetricsHistory self-throttles by sim time); only the DOM
    // paint below is throttled by wall-clock time. A metric currently `unavailable` (e.g.
    // collision_rate while coarse-grained) is sampled as `null`, so its sparkline/CSV column shows
    // a gap/empty cell rather than a value that can't be compared to a real mill.
    history.push(state.simTime, values);
    lastSimTime = state.simTime;
    exportButton.disabled = history.length === 0;
    latestState = state;

    // Feed the report aggregator once per *new* metrics object (worker.ts throttles how often
    // `state.metrics` actually gets a fresh value, well below render rate -- see FrameMessage's
    // doc comment on `metrics`/`fluidSurface`), not once per render frame, so a metrics snapshot
    // that hasn't changed since last frame isn't double-counted into the running mean/variance.
    if (state.metrics && state.metrics !== lastAggregatedMetrics) {
      reportAggregator.add(state.simTime, values, state.metrics.impact_energy_histogram);
      lastAggregatedMetrics = state.metrics;
    }

    // Window end: the first time sim time passes the end, mark the window complete; with the
    // auto checkbox on, also pause + generate + download once (re-armed by reset()/window edits).
    latestRpm = mill ? rpmOf(mill) : null;
    if (reportArmed && reportEndTimeS > 0 && reportAggregator.isPastEnd(state.simTime)) {
      reportArmed = false;
      windowComplete = true;
      if (autoCheckbox.checked) {
        onAutoReportPause();
        generateReport(state);
      }
    }

    if (nowMs - lastRenderMs < 100) return;
    lastRenderMs = nowMs;
    refreshReportUi();

    for (const spec of METRIC_SPECS) {
      const raw = values[spec.id];
      const value = raw === undefined ? null : raw;
      const reason = reasons[spec.id] ?? null;
      const span = valueEls.get(spec.id);
      if (span) {
        const text = reason ?? (value === null ? "-" : spec.format ? spec.format(value) : String(value));
        if (span.textContent !== text) span.textContent = text;
      }

      const row = rowEls.get(spec.id);
      if (row) {
        row.hidden = spec.hideWhenUnavailable === true && reason !== null;
      }

      const warnThreshold = WARN_THRESHOLDS.get(spec.id);
      if (warnThreshold !== undefined && row) {
        const warn = value !== null && value > warnThreshold;
        row.classList.toggle("is-warn", warn);
      }
    }

    // Generic guard: hide a group's `<details>` if every non-key row in it ended up hidden above
    // (nothing left to show), so it never dangles open above an empty group. Key-block rows are
    // excluded from this check -- they aren't rendered inside the group at all.
    for (const group of METRIC_GROUPS) {
      const details = groupDetailsEls.get(group);
      if (!details) continue;
      const groupSpecs = METRIC_SPECS.filter((s) => s.group === group && !s.key);
      details.hidden = groupSpecs.length > 0 && groupSpecs.every((s) => rowEls.get(s.id)?.hidden === true);
    }

    for (const spec of METRIC_SPECS) {
      if (!spec.sparkline) continue;
      const canvas = sparklineEls.get(spec.id);
      if (!canvas) continue;
      drawSparkline(canvas, history.seriesFor(spec.id));
    }

    const kReason = coarseGrainingReason(ctx);
    if (coarseGrainingWarnEl) {
      coarseGrainingWarnEl.hidden = kReason === null;
      if (kReason !== null && state.metrics) {
        const k = state.metrics.coarse_graining_factor.toFixed(2);
        const nTrue = state.metrics.true_ball_count.toFixed(0);
        coarseGrainingWarnEl.textContent =
          `Coarse-graining is active (k = ${k}). Collision rate and the impact energy distribution ` +
          `are hidden -- impact counts scale as 1/k^2 and per-impact energies as k^2, so both would ` +
          `be misleading. Raise Max balls above ${nTrue} to show them.`;
      }
    }

    if (histogramWrapEl) {
      histogramWrapEl.hidden = kReason !== null;
    }
    if (histogramCanvas && kReason === null) {
      drawImpactHistogram(histogramCanvas, state.metrics?.impact_energy_histogram);
    }

    const grindingChip = groupChipEls.get("Grinding");
    if (grindingChip) {
      grindingChip.hidden = kReason === null;
      if (kReason !== null && state.metrics) {
        grindingChip.textContent = `k = ${state.metrics.coarse_graining_factor.toFixed(2)}`;
      }
    }

    // Solver-health roll-up: same WARN_THRESHOLDS rule as the per-row `is-warn` highlight above,
    // just collected into the key block's one line instead of read off five separately-collapsed
    // rows. `lastWarnedIds` is read by the row's own click handler (above) to open the right
    // group(s).
    const warnedSpecs = METRIC_SPECS.filter((s) => {
      const threshold = WARN_THRESHOLDS.get(s.id);
      if (threshold === undefined) return false;
      const v = values[s.id];
      return v !== null && v !== undefined && v > threshold;
    });
    lastWarnedIds = warnedSpecs.map((s) => s.id);
    const isHealthWarn = warnedSpecs.length > 0;
    const healthText = isHealthWarn ? `⚠ ${warnedSpecs.map((s) => s.label).join(", ")}` : "OK";
    if (healthValue.textContent !== healthText) healthValue.textContent = healthText;
    healthRow.classList.toggle("is-warn", isHealthWarn);
  }

  return {
    update,
    reset(): void {
      history.reset();
      // Fresh aggregation window on every simulation reset (worker "init"/"ready" cycle, see
      // main.ts) -- otherwise a reset would silently mix pre-/post-reset samples into one report,
      // and a previously-fired auto-trigger would stay permanently disarmed.
      reportAggregator.reset();
      reportArmed = true;
      windowComplete = false;
      lastAggregatedMetrics = null;
      refreshReportUi();
    },
  };
}
