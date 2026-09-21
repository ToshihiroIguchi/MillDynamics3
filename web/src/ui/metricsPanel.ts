// Right-side metrics panel (docs/PLAN.md ss4.4): grouped read-outs of every `MetricSpec`
// (metrics/specs.ts), a rolling history feeding a few sparklines plus a CSV export, and a small
// impact-energy histogram chart. All the data plumbing (spec table, history ring buffer, sparkline
// renderer, `Metrics` type) lives elsewhere; this file only builds/updates the DOM.

import { coarseGrainingReason, METRIC_GROUPS, METRIC_SPECS, type MetricContext, type MetricGroupName } from "../metrics/specs";
import { MetricsHistory, type CsvColumn } from "../metrics/history";
import type { ImpactEnergyHistogram } from "../metrics/types";
import { criticalSpeedRpm, percentCriticalOf, rpmOf } from "../params/derived";
import type { ParamsJson } from "../protocol";
import type { AppState } from "../state";
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

  ctx.fillStyle = "#9aa3b0";
  ctx.font = "9px system-ui, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(`${maxCount.toFixed(1)}/s`, padLeft, 10);

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
// threshold-0 case of that same rule, not a different rule.
const COUPLING_CLAMP_HITS_WARN_THRESHOLD = 0;
const SUBSTEP_DISPLACEMENT_WARN_THRESHOLD = 0.5;
const MAX_COMPRESSION_ERROR_WARN_THRESHOLD = 0.05;
const MAX_BALL_OVERLAP_WARN_THRESHOLD = 0.5;

const WARN_THRESHOLDS = new Map<string, number>([
  ["coupling_clamp_hits", COUPLING_CLAMP_HITS_WARN_THRESHOLD],
  ["substep_displacement", SUBSTEP_DISPLACEMENT_WARN_THRESHOLD], // confirmed against metrics/specs.ts
  ["max_compression_error", MAX_COMPRESSION_ERROR_WARN_THRESHOLD],
  ["max_ball_overlap", MAX_BALL_OVERLAP_WARN_THRESHOLD],
]);

export function createMetricsPanel(container: HTMLElement): MetricsPanel {
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

  const valueEls = new Map<string, HTMLSpanElement>();
  const sparklineEls = new Map<string, HTMLCanvasElement>();
  const rowEls = new Map<string, HTMLElement>();
  const groupHeadingEls = new Map<MetricGroupName, HTMLHeadingElement>();
  let histogramCanvas: HTMLCanvasElement | null = null;
  let histogramWrapEl: HTMLElement | null = null;
  let coarseGrainingWarnEl: HTMLElement | null = null;

  for (const group of METRIC_GROUPS) {
    const section = document.createElement("section");
    const h3 = document.createElement("h3");
    h3.textContent = group;
    section.appendChild(h3);
    groupHeadingEls.set(group, h3);

    for (const spec of METRIC_SPECS.filter((s) => s.group === group)) {
      const row = document.createElement("div");
      row.className = "metric-row";

      const label = document.createElement("span");
      label.className = "metric-label";
      label.textContent = spec.unit ? `${spec.label} (${spec.unit})` : spec.label;

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

      section.appendChild(row);
      valueEls.set(spec.id, value);
      rowEls.set(spec.id, row);
    }

    container.appendChild(section);

    if (group === "Grinding") {
      const warnEl = document.createElement("div");
      warnEl.className = "metrics-warn";
      warnEl.hidden = true;
      container.appendChild(warnEl);
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
      container.appendChild(histWrap);
      histogramCanvas = canvas;
      histogramWrapEl = histWrap;
    }
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

    if (nowMs - lastRenderMs < 100) return;
    lastRenderMs = nowMs;

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

    // Generic guard: hide a group's heading if every row in it ended up hidden above (nothing left
    // to head), so it never dangles above an empty group.
    for (const group of METRIC_GROUPS) {
      const heading = groupHeadingEls.get(group);
      if (!heading) continue;
      const groupSpecs = METRIC_SPECS.filter((s) => s.group === group);
      heading.hidden = groupSpecs.length > 0 && groupSpecs.every((s) => rowEls.get(s.id)?.hidden === true);
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
  }

  return {
    update,
    reset(): void {
      history.reset();
    },
  };
}
