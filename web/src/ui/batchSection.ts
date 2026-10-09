// "Repeat runs (seeds)" section of the Grinding proxy (3D) group (docs/METRICS.md). Runs the SAME
// Simulation with the SAME applied params as the main screen, headless in a dedicated worker, once
// per seed (1..N), then shows mean +- standard error of the proxy rows, ratios against the stored
// reference, and downloads. All logic lives in the pure modules under batch/; this file only builds
// and updates the DOM.

import { buildBatchCsv, buildBatchJson, type BatchRecord } from "../batch/exporters";
import {
  estimateSeconds,
  formatDuration,
  refusalReason,
  seedList,
  summarizeParams,
  validateInputs,
  type BatchInputs,
  type LiveSnapshot,
} from "../batch/plan";
import { formatPercent, formatSig, rowInfo } from "../batch/rows";
import { BatchCancelled, BatchRunner } from "../batch/runner";
import { buildSensitivityConditions, buildSensitivityRows, exceedsSensitivity } from "../batch/sensitivity";
import {
  combineSeeds,
  DRIFT_WARN_PER_REV,
  parseSeedResult,
  RATIO_ERROR_WARN,
  steadyState,
  type Combined,
  type SeedResult,
} from "../batch/stats";
import { buildTable, referenceFromBatch } from "../batch/table";
import { formatRatio, settingsMismatch, settingsOf, type ProxyReference } from "../metrics/reference";
import type { ParamsJson } from "../protocol";

export interface BatchSectionDeps {
  /** Current applied state of the main simulation, or null before it is ready. */
  getLive(): LiveSnapshot | null;
  getReference(): ProxyReference | null;
  /** Stores a new reference (persisted by the caller) -- the same path as "Set reference". */
  setReference(ref: ProxyReference): void;
}

export interface BatchSection {
  element: HTMLElement;
  /** Called from the metrics panel's periodic repaint. */
  update(): void;
}

interface FinishedRun {
  params: ParamsJson;
  inputs: BatchInputs;
  results: SeedResult[];
  combined: Combined;
  thresholds: readonly number[];
  k: number;
}

function el<K extends keyof HTMLElementTagNameMap>(tag: K, className?: string, text?: string): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  if (className) e.className = className;
  if (text !== undefined) e.textContent = text;
  return e;
}

function download(filename: string, mime: string, content: string): void {
  const url = URL.createObjectURL(new Blob([content], { type: mime }));
  try {
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    a.click();
  } finally {
    URL.revokeObjectURL(url);
  }
}

function labeledInput(parent: HTMLElement, id: string, label: string, value: string, step: string, min: string): HTMLInputElement {
  const row = el("div", "metrics-report-row");
  const l = el("label", undefined, label);
  l.htmlFor = id;
  const input = el("input", "metrics-report-input");
  input.type = "number";
  input.id = id;
  input.min = min;
  input.step = step;
  input.value = value;
  row.append(l, input);
  parent.appendChild(row);
  return input;
}

const HELP_TEXT =
  "Why seeds: the charge is chaotic, so one run is one sample; the spread between seeds is the " +
  "statistical error of every row. Why whole revolutions: the motion is periodic per revolution, so a " +
  "partial revolution biases the mean. window_* rows are exact window means from the cumulative " +
  "energies (not from the smoothed live values). Ratios between conditions are only trustworthy where " +
  "their relative error is small.";

export function createBatchSection(deps: BatchSectionDeps): BatchSection {
  const runner = new BatchRunner();
  let busy = false;
  let cancelled = false;
  let finished: FinishedRun | null = null;
  let renderedRef: ProxyReference | null | undefined;

  const root = el("section", "batch-section");
  root.appendChild(el("h4", "batch-heading", "Repeat runs (seeds)"));
  const help = el("div", "metrics-note", HELP_TEXT);
  root.appendChild(help);

  const seedsInput = labeledInput(root, "batch-seeds", "Seeds (1..N)", "3", "1", "1");
  const warmupInput = labeledInput(root, "batch-warmup", "Warm-up revolutions", "4", "1", "0");
  const measureInput = labeledInput(root, "batch-measure", "Measured revolutions", "4", "1", "1");

  const summaryEl = el("div", "batch-summary");
  const noteEl = el("div", "metrics-note");
  noteEl.textContent =
    "The run uses exactly the settings of the main screen (applied parameters, including slurry if " +
    "enabled, which makes it slow). Pause the main simulation while it runs to free CPU.";
  const hintEl = el("div", "metrics-note");
  root.append(summaryEl, noteEl, hintEl);

  const controls = el("div", "metrics-ref-controls");
  const runBtn = el("button", "metrics-export", "Run seeds");
  runBtn.type = "button";
  const cancelBtn = el("button", "metrics-export", "Cancel");
  cancelBtn.type = "button";
  cancelBtn.disabled = true;
  const sensBtn = el("button", "metrics-export", "Sensitivity");
  sensBtn.type = "button";
  sensBtn.title =
    "Re-run seed 1 with ball-ball/ball-wall friction and restitution each scaled by 0.8 and 1.2 (9 runs including the base)";
  controls.append(runBtn, cancelBtn, sensBtn);
  root.appendChild(controls);

  const progressEl = el("div", "batch-progress");
  const errorEl = el("div", "metrics-warn");
  errorEl.hidden = true;
  root.append(progressEl, errorEl);

  const resultEl = el("div", "batch-result");
  root.appendChild(resultEl);
  const sensEl = el("div", "batch-result");
  root.appendChild(sensEl);

  function showError(text: string | null): void {
    errorEl.hidden = text === null;
    if (text !== null) errorEl.textContent = text;
  }

  function setBusy(b: boolean): void {
    busy = b;
    runBtn.disabled = b;
    sensBtn.disabled = b;
    cancelBtn.disabled = !b;
    for (const i of [seedsInput, warmupInput, measureInput]) i.disabled = b;
  }

  async function runOne(paramsJson: string, seed: number, inputs: BatchInputs): Promise<SeedResult> {
    const text = await runner.run(paramsJson, seed, inputs.warmup, inputs.measure);
    return parseSeedResult(text);
  }

  /** Common start: validation, refusal, snapshot of the applied params. Returns null (error shown) on refusal. */
  function begin(seedsText: string): { inputs: BatchInputs; live: LiveSnapshot; params: ParamsJson; k: number } | null {
    showError(null);
    const inputs = validateInputs(seedsText, warmupInput.value, measureInput.value);
    if (typeof inputs === "string") {
      showError(inputs);
      return null;
    }
    const live = deps.getLive();
    const refusal = refusalReason(live);
    if (refusal !== null || !live || !live.params) {
      showError(refusal ?? "The simulation is not ready yet.");
      return null;
    }
    return { inputs, live, params: structuredClone(live.params), k: live.coarseGrainingFactor ?? 1 };
  }

  function handleFailure(err: unknown): void {
    if (err instanceof BatchCancelled) {
      progressEl.textContent = "Cancelled.";
      return;
    }
    progressEl.textContent = "Failed.";
    showError(`Run failed: ${err instanceof Error ? err.message : String(err)}`);
  }

  runBtn.addEventListener("click", () => {
    const start = begin(seedsInput.value);
    if (!start || busy) return;
    const { inputs, live, params, k } = start;
    const paramsJson = JSON.stringify(params);
    const seeds = seedList(inputs.seeds);
    cancelled = false;
    setBusy(true);
    resultEl.replaceChildren();
    sensEl.replaceChildren();
    finished = null;
    void (async () => {
      try {
        const results: SeedResult[] = [];
        for (const [i, seed] of seeds.entries()) {
          progressEl.textContent = `seed ${i + 1}/${seeds.length}`;
          results.push(await runOne(paramsJson, seed, inputs));
          if (cancelled) throw new BatchCancelled();
        }
        const combined = combineSeeds(results);
        const resultK = Math.max(...results.map((r) => r.coarse_graining_factor));
        finished = { params, inputs, results, combined, thresholds: live.thresholdsScaled, k: Math.max(k, resultK) };
        progressEl.textContent = `Done: ${seeds.length} seed(s).`;
        renderResult();
      } catch (err) {
        handleFailure(err);
      } finally {
        setBusy(false);
      }
    })();
  });

  cancelBtn.addEventListener("click", () => {
    cancelled = true;
    runner.cancel();
  });

  sensBtn.addEventListener("click", () => {
    const start = begin("1");
    if (!start || busy) return;
    const { inputs, live, params } = start;
    let conditions;
    try {
      conditions = buildSensitivityConditions(params);
    } catch (err) {
      showError(err instanceof Error ? err.message : String(err));
      return;
    }
    cancelled = false;
    setBusy(true);
    sensEl.replaceChildren();
    void (async () => {
      try {
        const results: SeedResult[] = [];
        for (const [i, c] of conditions.entries()) {
          progressEl.textContent = `sensitivity ${i + 1}/${conditions.length}: ${c.label}`;
          results.push(await runOne(JSON.stringify(c.params), 1, inputs));
          if (cancelled) throw new BatchCancelled();
        }
        progressEl.textContent = `Sensitivity done (${conditions.length} runs, seed 1).`;
        renderSensitivity(
          conditions.map((c) => c.label),
          results,
          live.thresholdsScaled,
        );
      } catch (err) {
        handleFailure(err);
      } finally {
        setBusy(false);
      }
    })();
  });

  function cell(tr: HTMLElement, text: string, className?: string, title?: string): HTMLTableCellElement {
    const td = el("td", className, text);
    if (title) td.title = title;
    tr.appendChild(td);
    return td;
  }

  function renderResult(): void {
    resultEl.replaceChildren();
    renderedRef = deps.getReference();
    const run = finished;
    if (!run) return;
    const ref = renderedRef;
    const table = buildTable(run.combined, run.thresholds, ref);

    const info = el("div", "metrics-note");
    info.textContent =
      `${run.combined.seeds} seed(s) [${run.results.map((r) => r.seed).join(", ")}], ` +
      `${run.inputs.warmup} warm-up + ${run.inputs.measure} measured revolution(s), window ${run.results[0]?.window_s.toFixed(2)} s. ` +
      `Settings: ${summarizeParams(run.params)}. Mean ± standard error (n-1).`;
    resultEl.appendChild(info);

    const ss = steadyState(run.results);
    const steady = el("div", ss.warn ? "metrics-note is-warn-note" : "metrics-note");
    steady.textContent =
      `Steady state: wall-power drift ${(ss.mean * 100).toFixed(1)} %/rev (mean over seeds; mean |drift| ` +
      `${(ss.meanAbs * 100).toFixed(1)} %, max ${(ss.maxAbs * 100).toFixed(1)} %).` +
      (ss.warn
        ? ` Above ${DRIFT_WARN_PER_REV * 100} %/rev: the charge is still evolving, lengthen the warm-up.`
        : "");
    resultEl.appendChild(steady);

    if (ref) {
      const diff = settingsMismatch(ref.settings, settingsOf(run.params, run.k));
      const refNote = el("div", diff.length > 0 ? "metrics-note is-warn-note" : "metrics-note");
      refNote.textContent =
        `Ratio = batch mean / reference${ref.fromBatch ? " (reference is itself a batch mean: both errors propagated)" : " (reference treated as exact)"}. ` +
        `Ratios with relative error above ${RATIO_ERROR_WARN * 100} % are marked "!" and should not be read.` +
        (diff.length > 0
          ? ` Warning: numerical settings differ from the reference run (${diff.join(", ")}); ratios are less reliable.`
          : "");
      resultEl.appendChild(refNote);
    }

    const wrap = el("div", "batch-table-wrap");
    const tbl = el("table", "batch-table");
    const head = el("tr");
    for (const h of ["Row", "Mean ± SE", "Rel. SE", ...(ref ? ["× ref"] : [])]) head.appendChild(el("th", undefined, h));
    tbl.appendChild(head);
    for (const row of table) {
      const tr = el("tr");
      cell(tr, row.info.unit ? `${row.info.label} (${row.info.unit})` : row.info.label, "batch-label", row.info.name);
      cell(tr, `${formatSig(row.mean)} ± ${formatSig(row.se, 2)}`);
      cell(tr, formatPercent(row.relSe));
      if (ref) {
        if (row.ratio === null) cell(tr, "-");
        else {
          const err = row.ratioRelErr === null ? "" : ` ± ${formatPercent(row.ratioRelErr, 0)}`;
          const td = cell(
            tr,
            `${formatRatio(row.ratio)}${err}${row.ratioUncertain ? " !" : ""}`,
            row.ratioUncertain ? "batch-uncertain" : undefined,
            row.ratioUncertain ? "Relative error of the ratio exceeds 10 %: not trustworthy" : `Reference: ${formatSig(row.refValue ?? NaN)}`,
          );
          td.classList.add("batch-ratio");
        }
      }
      tbl.appendChild(tr);
    }
    wrap.appendChild(tbl);
    resultEl.appendChild(wrap);

    const actions = el("div", "metrics-ref-controls");
    const useBtn = el("button", "metrics-export", "Use as reference");
    useBtn.type = "button";
    useBtn.addEventListener("click", () => {
      deps.setReference(referenceFromBatch(run.combined, run.params, run.k));
      renderResult();
    });
    const csvBtn = el("button", "metrics-export", "Download CSV");
    csvBtn.type = "button";
    const jsonBtn = el("button", "metrics-export", "Download JSON");
    jsonBtn.type = "button";
    const record = (): BatchRecord => ({
      params: run.params,
      warmupRevolutions: run.inputs.warmup,
      measureRevolutions: run.inputs.measure,
      results: run.results,
      combined: run.combined,
      thresholdsScaled: run.thresholds,
      reference: deps.getReference(),
      createdAt: new Date().toISOString(),
    });
    csvBtn.addEventListener("click", () => download("milldynamics-repeat-runs.csv", "text/csv", buildBatchCsv(record())));
    jsonBtn.addEventListener("click", () =>
      download("milldynamics-repeat-runs.json", "application/json", buildBatchJson(record())),
    );
    actions.append(useBtn, csvBtn, jsonBtn);
    resultEl.appendChild(actions);
  }

  function renderSensitivity(labels: string[], results: SeedResult[], thresholds: readonly number[]): void {
    sensEl.replaceChildren();
    const rows = buildSensitivityRows(results);
    const note = el("div", "metrics-note");
    note.textContent =
      "Sensitivity to the contact coefficients (one seed, so it includes seed noise of the same order as the " +
      "relative SE above). Each cell is the change versus the base run; changes above 10 % are highlighted. " +
      "Ratios between conditions are only trustworthy where these changes are small.";
    sensEl.appendChild(note);
    const wrap = el("div", "batch-table-wrap");
    const tbl = el("table", "batch-table");
    const head = el("tr");
    head.appendChild(el("th", undefined, "Row"));
    head.appendChild(el("th", undefined, "Base"));
    for (const l of labels.slice(1)) head.appendChild(el("th", undefined, l.replace("friction", "fric.").replace("restitution", "rest.")));
    tbl.appendChild(head);
    for (const row of rows) {
      const info = rowInfo(row.name, thresholds);
      const tr = el("tr");
      cell(tr, info.unit ? `${info.label} (${info.unit})` : info.label, "batch-label", row.name);
      cell(tr, formatSig(row.base));
      for (const c of row.changes) {
        cell(
          tr,
          c === null ? "-" : `${c >= 0 ? "+" : ""}${(c * 100).toFixed(0)} %`,
          exceedsSensitivity(c) ? "batch-uncertain" : undefined,
        );
      }
      tbl.appendChild(tr);
    }
    wrap.appendChild(tbl);
    sensEl.appendChild(wrap);
  }

  function update(): void {
    const live = deps.getLive();
    const text = `Applied settings: ${summarizeParams(live?.params ?? null)}`;
    if (summaryEl.textContent !== text) summaryEl.textContent = text;
    const seedsN = Number(seedsInput.value);
    const w = Number(warmupInput.value);
    const m = Number(measureInput.value);
    let hint = "";
    if (live && Number.isInteger(seedsN) && seedsN >= 1 && Number.isFinite(w) && Number.isFinite(m) && m >= 1) {
      const est = estimateSeconds(live.params, { warmup: w, measure: m }, seedsN, live.subStepsPerSecondAchieved);
      hint =
        est === null
          ? "Duration unknown (the main simulation is not running, so no sub-step rate is available). Expect it to be slow, especially with slurry."
          : `Rough estimate at the main screen's current sub-step rate: ${formatDuration(est)} (${seedsN} run(s)); slower if the main simulation keeps running.`;
    }
    if (hintEl.textContent !== hint) hintEl.textContent = hint;
    if (!busy) {
      const refusal = refusalReason(live);
      if (refusal !== null && errorEl.hidden) {
        // Shown passively through the hint until the user presses Run; do not overwrite run errors.
        runBtn.title = refusal;
      } else runBtn.title = "";
    }
    if (finished && renderedRef !== deps.getReference()) renderResult();
  }

  return { element: root, update };
}
