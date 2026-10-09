// "Scale-up study" panel: compares grinding-proxy metrics RELATIVELY between a lab (reference)
// condition and candidate large-scale conditions. All numerics run in scaleup/scaleupWorker.ts;
// pure logic (params building, ratios, CSV, extrapolation) lives in scaleup/logic.ts.

import type { ParamsJson } from "../protocol";
import { drawLineChart, seriesColor, type Series } from "../scaleup/charts";
import { CancelledError, ScaleupClient } from "../scaleup/client";
import {
  ballCount,
  ballCountAtRatio,
  binCentres,
  buildParams,
  buildScaleupConfig,
  combineEnsembles,
  type Condition,
  DEFAULT_PACKING_FRACTION_2D,
  DEFAULT_STUDY_CONFIG,
  ddRatio,
  defaultCondition,
  ensembleWarnings,
  extrapolate,
  formatNumber,
  getMetric,
  histogramRates,
  intervalsCsv,
  isSensitive,
  MAX_BALLS,
  MIN_BALLS,
  metricDefs,
  metricLabelForKey,
  metricsCsv,
  primaryMetricKeys,
  ratioWithError,
  recordedDuration,
  seedList,
  sensitivityChange,
  steadyIndicators,
  type StudyConfig,
  validateCondition,
  validateStudyConfig,
  DISPLACEMENT_WARN,
  DRIFT_WARN_REL_PER_REV,
  RATIO_ERROR_NOTE_REL,
} from "../scaleup/logic";
import { CLASS_NAMES, type EnsembleResult, type SensitivityCase, type SweepResult } from "../scaleup/types";

function el<K extends keyof HTMLElementTagNameMap>(tag: K, cls?: string, text?: string): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function button(label: string, onClick: () => void, cls = "su-btn"): HTMLButtonElement {
  const b = el("button", cls, label);
  b.type = "button";
  b.addEventListener("click", onClick);
  return b;
}

function downloadText(filename: string, text: string, mime: string): void {
  const url = URL.createObjectURL(new Blob([text], { type: mime }));
  try {
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    a.click();
  } finally {
    URL.revokeObjectURL(url);
  }
}

function numInput(value: number, step: string, onChange: (v: number) => void, min?: number): HTMLInputElement {
  const i = el("input", "su-num");
  i.type = "number";
  i.step = step;
  if (min !== undefined) i.min = String(min);
  i.value = String(value);
  i.addEventListener("input", () => onChange(i.value === "" ? Number.NaN : Number(i.value)));
  return i;
}

function labelled(text: string, control: HTMLElement, title?: string): HTMLLabelElement {
  const l = el("label", "su-field");
  l.append(el("span", "su-field-label", text), control);
  if (title) l.title = title;
  return l;
}

/** Columns of the condition table. `get`/`set` map a numeric Condition field. */
interface NumCol {
  head: string;
  title: string;
  key: keyof Omit<Condition, "name">;
  step: string;
}

const COLS: NumCol[] = [
  { head: "Drum D (mm)", title: "Drum inner diameter", key: "drumDiameterMm", step: "1" },
  { head: "Ball d (mm)", title: "Ball diameter", key: "ballDiameterMm", step: "0.1" },
  { head: "% crit.", title: "Drum speed as % of critical speed", key: "percentCritical", step: "1" },
  { head: "Fill J", title: "Ball bed fill fraction J", key: "fillFraction", step: "0.01" },
  { head: "Density", title: "Ball density (kg/m3)", key: "ballDensity", step: "100" },
  { head: "Lifters", title: "Lifter count (0 = smooth wall); lifter height = d, base width = 1.5 d, top width = d", key: "lifterCount", step: "1" },
  { head: "e bb", title: "Restitution ball-ball", key: "restitutionBallBall", step: "0.05" },
  { head: "e bw", title: "Restitution ball-wall", key: "restitutionBallWall", step: "0.05" },
  { head: "mu bb", title: "Friction ball-ball", key: "frictionBallBall", step: "0.05" },
  { head: "mu bw", title: "Friction ball-wall", key: "frictionBallWall", step: "0.05" },
  { head: "Visc. (Pa s)", title: "Lubrication-cushion viscosity between balls; 0 = dry", key: "viscosityPaS", step: "0.1" },
];

interface StudyResult {
  names: string[];
  conditions: Condition[];
  params: ParamsJson[];
  config: StudyConfig;
  ens: EnsembleResult[];
}

export interface ScaleupPanel {
  el: HTMLElement;
  isVisible(): boolean;
  setVisible(visible: boolean): void;
}

/** Creates the (initially hidden) full-screen scale-up panel. `onClosed` fires when the user closes it from inside. */
export function createScaleupPanel(onClosed: () => void): ScaleupPanel {
  const root = el("div", "su-overlay");
  root.hidden = true;
  root.setAttribute("role", "dialog");
  root.setAttribute("aria-label", "Scale-up study");

  let client: ScaleupClient | null = null;
  let baseParams: ParamsJson | null = null;
  let packing = DEFAULT_PACKING_FRACTION_2D;
  let conditions: Condition[] = [defaultCondition("Lab"), { ...defaultCondition("Candidate 1"), drumDiameterMm: 200, ballDiameterMm: 4 }];
  const cfg: StudyConfig = { ...DEFAULT_STUDY_CONFIG };
  let busy = false;
  let result: StudyResult | null = null;
  let sweepResult: { baseName: string; ratios: number[]; data: SweepResult; config: StudyConfig; params: ParamsJson } | null = null;
  let sensResult: { name: string; rel: number; cases: SensitivityCase[]; params: ParamsJson } | null = null;

  function getClient(): ScaleupClient {
    if (!client) {
      client = new ScaleupClient();
      client.ready
        .then((r) => {
          baseParams = JSON.parse(r.defaultParamsJson) as ParamsJson;
          const pf = (baseParams.media as { packing_fraction_2d?: number } | undefined)?.packing_fraction_2d;
          if (typeof pf === "number") packing = pf;
          refreshDerived();
        })
        .catch(() => undefined);
    }
    return client;
  }

  async function getBase(): Promise<ParamsJson> {
    const c = getClient();
    if (!baseParams) await c.ready;
    if (!baseParams) throw new Error("default parameters unavailable");
    return baseParams;
  }

  // ---- header -----------------------------------------------------------------------------
  const header = el("div", "su-header");
  header.append(el("h2", "su-title", "Scale-up study"));
  const closeBtn = button("Close", () => {
    setVisible(false);
    onClosed();
  }, "su-btn su-close");
  header.append(closeBtn);
  const intro = el(
    "p",
    "su-note",
    "Relative comparison only: absolute values carry the 2D-model bias and are not trusted. Compare each candidate to the first (reference) condition through the ratio columns. Grinding is not simulated; the metrics are impact statistics (stress energies, frequencies) and power. The slurry is not resolved (dry, optionally with a lubrication-cushion viscosity between balls).",
  );

  // ---- conditions -------------------------------------------------------------------------
  const condSection = el("section", "su-section");
  condSection.append(el("h3", "su-h3", "Conditions"));
  const condTableWrap = el("div", "su-scroll");
  const condMsg = el("div", "su-msg");
  const condButtons = el("div", "su-row");
  condButtons.append(
    button("Add candidate", () => {
      const ref = conditions[0] ?? defaultCondition("Lab");
      conditions.push({ ...ref, name: `Candidate ${conditions.length}` });
      renderConditions();
    }),
  );
  condSection.append(condTableWrap, condButtons, condMsg);

  let derivedCells: { n: HTMLElement; dd: HTMLElement; row: HTMLElement }[] = [];

  function renderConditions(): void {
    condTableWrap.replaceChildren();
    derivedCells = [];
    const table = el("table", "su-table su-cond");
    const head = el("tr");
    head.append(el("th", undefined, "#"), el("th", undefined, "Name"));
    for (const c of COLS) {
      const th = el("th", undefined, c.head);
      th.title = c.title;
      head.append(th);
    }
    head.append(el("th", undefined, "Balls"), el("th", undefined, "D/d"), el("th"));
    const thead = el("thead");
    thead.append(head);
    const tbody = el("tbody");
    conditions.forEach((cond, idx) => {
      const tr = el("tr");
      tr.append(el("td", undefined, idx === 0 ? "ref" : String(idx)));
      const nameTd = el("td");
      const nameIn = el("input", "su-text");
      nameIn.type = "text";
      nameIn.value = cond.name;
      nameIn.addEventListener("input", () => {
        cond.name = nameIn.value;
        refreshSelects();
      });
      nameTd.append(nameIn);
      tr.append(nameTd);
      for (const col of COLS) {
        const td = el("td");
        td.append(
          numInput(cond[col.key], col.step, (v) => {
            cond[col.key] = v;
            refreshDerived();
          }, 0),
        );
        tr.append(td);
      }
      const nTd = el("td", "su-derived");
      const ddTd = el("td", "su-derived");
      tr.append(nTd, ddTd);
      const rm = el("td");
      if (idx > 0) {
        rm.append(
          button("Remove", () => {
            conditions.splice(idx, 1);
            renderConditions();
          }, "su-btn su-small"),
        );
      }
      tr.append(rm);
      derivedCells.push({ n: nTd, dd: ddTd, row: tr });
      tbody.append(tr);
    });
    table.append(thead, tbody);
    condTableWrap.append(table);
    refreshDerived();
    refreshSelects();
  }

  /** Validation problems of the whole input state (empty = can run). */
  function allProblems(): string[] {
    const out: string[] = [];
    conditions.forEach((c, i) => {
      for (const p of validateCondition(c, packing)) out.push(`${c.name || `condition ${i + 1}`}: ${p}`);
    });
    for (const p of validateStudyConfig(cfg)) out.push(p);
    return out;
  }

  function refreshDerived(): void {
    conditions.forEach((c, i) => {
      const cells = derivedCells[i];
      if (!cells) return;
      const n = ballCount(c, packing);
      cells.n.textContent = Number.isFinite(n) ? String(Math.round(n)) : "-";
      cells.dd.textContent = formatNumber(ddRatio(c), 3);
      const bad = n > MAX_BALLS || n < MIN_BALLS || validateCondition(c, packing).length > 0;
      cells.row.classList.toggle("su-bad", bad);
    });
    const problems = allProblems();
    condMsg.textContent = problems.length ? problems.join(" | ") : "";
    condMsg.classList.toggle("su-error", problems.length > 0);
    updateRunState();
  }

  // ---- global config ----------------------------------------------------------------------
  const cfgSection = el("section", "su-section");
  cfgSection.append(el("h3", "su-h3", "Run settings"));
  const cfgRow = el("div", "su-row su-wrap");
  const cfgFields: [string, keyof StudyConfig, string, string][] = [
    ["Seeds", "seedCount", "1", "Ensemble size: number of random initial placements (mean and standard error are over seeds)"],
    ["Warm-up revolutions", "warmupRevolutions", "1", "Revolutions simulated before recording (charge settling)"],
    ["Measured revolutions", "measureRevolutions", "1", "Revolutions recorded; use a whole number"],
    ["Samples / revolution", "samplesPerRevolution", "1", "Recording intervals per revolution (time-series resolution)"],
    ["dt ratio", "dtRatio", "0.01", "Sub-step as a multiple of sqrt(d/g)"],
    ["DEM iterations", "demIterations", "1", "Contact-solver iterations per sub-step"],
  ];
  for (const [label, key, step, title] of cfgFields) {
    cfgRow.append(
      labelled(label, numInput(cfg[key], step, (v) => {
        cfg[key] = v;
        refreshDerived();
      }, 0), title),
    );
  }
  cfgSection.append(cfgRow);

  // ---- run bar ----------------------------------------------------------------------------
  const runBar = el("div", "su-row su-wrap su-runbar");
  const runBtn = button("Run study", () => void runStudy(), "su-btn su-primary");
  const cancelBtn = button("Cancel", () => cancelRun(), "su-btn");
  cancelBtn.disabled = true;
  const progress = el("span", "su-progress");
  const runMsg = el("div", "su-msg");
  runBar.append(runBtn, cancelBtn, progress);

  function updateRunState(): void {
    const problems = allProblems();
    runBtn.disabled = busy || problems.length > 0;
    runBtn.title = problems.length ? problems.join("\n") : "";
    cancelBtn.disabled = !busy;
    sweepRunBtn.disabled = busy || problems.length > 0;
    sensRunBtn.disabled = busy || problems.length > 0;
  }

  function setBusy(b: boolean, text = ""): void {
    busy = b;
    progress.textContent = text;
    updateRunState();
  }

  function cancelRun(): void {
    client?.cancel();
  }

  function reportError(target: HTMLElement, err: unknown, context: string): void {
    if (err instanceof CancelledError) {
      target.textContent = "Cancelled.";
      target.classList.remove("su-error");
      return;
    }
    const message = err instanceof Error ? err.message : String(err);
    target.textContent = `${context}: ${message}`;
    target.classList.add("su-error");
  }

  async function runStudy(): Promise<void> {
    runMsg.textContent = "";
    runMsg.classList.remove("su-error");
    const snapshot = conditions.map((c) => ({ ...c }));
    const config = { ...cfg };
    const seeds = seedList(config.seedCount);
    setBusy(true, "starting worker...");
    try {
      const base = await getBase();
      const params: ParamsJson[] = [];
      const ens: EnsembleResult[] = [];
      for (let ci = 0; ci < snapshot.length; ci++) {
        const cond = snapshot[ci];
        if (!cond) continue;
        const p = buildParams(base, cond);
        params.push(p);
        const paramsJson = JSON.stringify(p);
        const parts: EnsembleResult[] = [];
        for (let si = 0; si < seeds.length; si++) {
          setBusy(true, `condition ${ci + 1}/${snapshot.length} (${cond.name}), seed ${si + 1}/${seeds.length}`);
          const seed = seeds[si] ?? 1;
          const text = await getClient().run(paramsJson, JSON.stringify(buildScaleupConfig(config, [seed])));
          parts.push(JSON.parse(text) as EnsembleResult);
        }
        ens.push(combineEnsembles(parts));
      }
      result = { names: snapshot.map((c) => c.name), conditions: snapshot, params, config, ens };
      runMsg.textContent = "Done.";
      renderResults();
    } catch (err) {
      reportError(runMsg, err, "Run failed");
    } finally {
      setBusy(false);
    }
  }

  // ---- results ----------------------------------------------------------------------------
  const resultsSection = el("section", "su-section");
  resultsSection.hidden = true;

  const warnBox = el("div", "su-warnings");
  const tableWrap = el("div", "su-scroll");
  const steadyWrap = el("div", "su-scroll");
  const chartControls = el("div", "su-row su-wrap");
  const classSel = el("select", "su-select");
  CLASS_NAMES.forEach((n, i) => {
    const o = el("option", undefined, n);
    o.value = String(i);
    classSel.append(o);
  });
  const logYBox = el("input");
  logYBox.type = "checkbox";
  logYBox.checked = true;
  const thrSel = el("select", "su-select");
  chartControls.append(
    labelled("Impact class", classSel),
    labelled("Log y (histogram)", logYBox),
    labelled("Energy threshold", thrSel),
  );
  const legend = el("div", "su-legend");
  const chartGrid = el("div", "su-charts");
  const mkChart = (title: string): HTMLCanvasElement => {
    const box = el("div", "su-chart");
    box.append(el("div", "su-chart-title", title));
    const cv = el("canvas", "su-canvas");
    box.append(cv);
    chartGrid.append(box);
    return cv;
  };
  const histCanvas = mkChart("Stress-energy histogram (SE / m g d), first seed");
  const powerCanvas = mkChart("Wall power per ball vs revolutions, first seed");
  const energyCanvas = mkChart("Cumulative effective energy per ball vs revolutions, first seed");
  const exportRow = el("div", "su-row su-wrap");
  exportRow.append(
    button("Download CSV (metrics)", () => {
      if (result) downloadText("scaleup-metrics.csv", metricsCsv(csvConds(result)), "text/csv");
    }),
    button("Download CSV (intervals)", () => {
      if (result) downloadText("scaleup-intervals.csv", intervalsCsv(csvConds(result)), "text/csv");
    }),
    button("Download JSON", () => {
      const r = result;
      if (!r) return;
      const payload = {
        generated: new Date().toISOString(),
        note: "Relative scale-up proxy study; params are the exact Params JSON passed to the solver per condition.",
        study_config: r.config,
        conditions: r.conditions.map((c, i) => ({
          condition: c,
          params: r.params[i],
          wasm_config_template: buildScaleupConfig(r.config, seedList(r.config.seedCount)),
          result: r.ens[i],
        })),
        sweep: sweepResult,
        sensitivity: sensResult,
      };
      downloadText("scaleup-results.json", JSON.stringify(payload, null, 1), "application/json");
    }),
  );
  resultsSection.append(
    el("h3", "su-h3", "Results"),
    warnBox,
    tableWrap,
    el("h4", "su-h4", "Steady-state indicators"),
    steadyWrap,
    el("h4", "su-h4", "Charts"),
    chartControls,
    legend,
    chartGrid,
    exportRow,
  );

  function csvConds(r: StudyResult): { name: string; ens: EnsembleResult }[] {
    return r.names.map((name, i) => ({ name, ens: r.ens[i] as EnsembleResult }));
  }

  function fmtMeanSe(m: { mean: number; se: number } | null): string {
    if (!m) return "-";
    return `${formatNumber(m.mean, 4)} ± ${formatNumber(m.se, 2)}`;
  }

  function renderResults(): void {
    if (!result) return;
    const r = result;
    resultsSection.hidden = false;
    const ref = r.ens[0];
    if (!ref) return;

    // warnings
    warnBox.replaceChildren();
    r.ens.forEach((e, i) => {
      for (const w of ensembleWarnings(r.names[i] ?? `condition ${i + 1}`, e)) {
        warnBox.append(el("div", "su-warning", w.text));
      }
    });

    // metric table
    const table = el("table", "su-table");
    const head = el("tr");
    head.append(el("th", undefined, "Metric"));
    r.names.forEach((n, i) => {
      head.append(el("th", undefined, i === 0 ? `${n} (reference)` : n));
      if (i > 0) head.append(el("th", undefined, `${n} / reference`));
    });
    const thead = el("thead");
    thead.append(head);
    const tbody = el("tbody");
    for (const def of metricDefs(ref.first.summary)) {
      const tr = el("tr");
      const label = el("td", "su-metric");
      label.append(el("div", undefined, def.label), el("div", "su-help", def.help));
      tr.append(label);
      const refM = getMetric(ref, def.key);
      r.ens.forEach((e, i) => {
        const m = getMetric(e, def.key);
        tr.append(el("td", "su-val", fmtMeanSe(m)));
        if (i > 0) {
          const td = el("td", "su-val su-ratio");
          if (m && refM) {
            const q = ratioWithError(m, refM);
            td.textContent = Number.isFinite(q.ratio) ? `${formatNumber(q.ratio, 3)} ± ${formatNumber(q.err, 2)}` : "-";
            if (q.noisy) {
              td.classList.add("su-noisy");
              const note = el("div", "su-help", `ratio error ${Number.isFinite(q.relErr) ? (q.relErr * 100).toFixed(0) : "?"}% > ${RATIO_ERROR_NOTE_REL * 100}%: add seeds`);
              td.append(note);
            }
          } else {
            td.textContent = "-";
          }
          tr.append(td);
        }
      });
      tbody.append(tr);
    }
    table.append(thead, tbody);
    tableWrap.replaceChildren(table);

    // steady-state table
    const st = el("table", "su-table");
    const sh = el("tr");
    for (const h of ["Condition", "Balls", "D/d", "Drift (/rev, worst seed)", "Power rel. std", "Max sub-step displacement / d"]) sh.append(el("th", undefined, h));
    const sHead = el("thead");
    sHead.append(sh);
    const sBody = el("tbody");
    r.ens.forEach((e, i) => {
      const ind = steadyIndicators(e);
      const tr = el("tr");
      tr.append(
        el("td", undefined, r.names[i] ?? ""),
        el("td", "su-val", String(e.first.summary.n_balls)),
        el("td", "su-val", formatNumber(e.first.summary.dd_ratio, 3)),
        el("td", ind.driftMax > DRIFT_WARN_REL_PER_REV ? "su-val su-noisy" : "su-val", `${(ind.driftMax * 100).toFixed(2)} %`),
        el("td", "su-val", `${(ind.relStdMax * 100).toFixed(1)} %`),
        el("td", ind.displacementMax > DISPLACEMENT_WARN ? "su-val su-noisy" : "su-val", ind.displacementMax.toFixed(3)),
      );
      sBody.append(tr);
    });
    st.append(sHead, sBody);
    steadyWrap.replaceChildren(st);

    // threshold selector
    const prevThr = thrSel.value;
    thrSel.replaceChildren();
    ref.first.summary.thresholds_j.forEach((thr, i) => {
      const rel = thr / ref.first.summary.scale_energy_j;
      const o = el("option", undefined, `thr${i}: SE > ${formatNumber(rel, 3)} m g d`);
      o.value = String(i);
      thrSel.append(o);
    });
    if (prevThr && Number(prevThr) < ref.first.summary.thresholds_j.length) thrSel.value = prevThr;

    legend.replaceChildren();
    r.names.forEach((n, i) => {
      const item = el("span", "su-legend-item");
      const sw = el("span", "su-swatch");
      sw.style.background = seriesColor(i);
      item.append(sw, document.createTextNode(n));
      legend.append(item);
    });
    drawCharts();
  }

  function drawCharts(): void {
    if (!result || root.hidden) return;
    const r = result;
    const cls = Number(classSel.value) || 0;
    const thr = Number(thrSel.value) || 0;
    const hist: Series[] = r.ens.map((e, i) => ({
      name: r.names[i] ?? "",
      xs: binCentres(e.first.se_bin_edges),
      ys: histogramRates(e.first.histogram_total, cls, e.first.summary.n_balls, recordedDuration(e.first.intervals)),
    }));
    drawLineChart(histCanvas, hist, {
      logX: true,
      logY: logYBox.checked,
      xLabel: "SE / (m g d)",
      yLabel: "impacts / ball / s / bin",
    });
    const power: Series[] = r.ens.map((e, i) => ({
      name: r.names[i] ?? "",
      xs: e.first.intervals.map((iv) => iv.revolutions),
      ys: e.first.intervals.map((iv) => iv.wall_power_w),
    }));
    drawLineChart(powerCanvas, power, { xLabel: "revolutions", yLabel: "wall power (W/ball)" });
    const energy: Series[] = r.ens.map((e, i) => ({
      name: r.names[i] ?? "",
      xs: e.first.intervals.map((iv) => iv.revolutions),
      ys: e.first.intervals.map((iv) => iv.cum_effective_energy_j[thr] ?? Number.NaN),
    }));
    drawLineChart(energyCanvas, energy, { xLabel: "revolutions", yLabel: `cum. effective energy thr${thr} (J/ball)` });
  }
  classSel.addEventListener("change", drawCharts);
  logYBox.addEventListener("change", drawCharts);
  thrSel.addEventListener("change", drawCharts);
  window.addEventListener("resize", drawCharts);
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", drawCharts);

  // ---- matching helper ---------------------------------------------------------------------
  const matchSection = el("section", "su-section");
  matchSection.append(el("h3", "su-h3", "Matching ball diameter"));
  const matchRow = el("div", "su-row su-wrap");
  let matchTargetMm = 600;
  let matchedD: number | null = null;
  const matchOut = el("span", "su-derived");
  const matchAdd = button("Add as candidate", () => {
    const ref = conditions[0];
    if (!ref || matchedD === null) return;
    conditions.push({ ...ref, name: `Match D=${matchTargetMm} mm`, drumDiameterMm: matchTargetMm, ballDiameterMm: matchedD });
    renderConditions();
  });
  matchAdd.disabled = true;
  async function updateMatch(): Promise<void> {
    const ref = conditions[0];
    matchedD = null;
    matchAdd.disabled = true;
    if (!ref || !(matchTargetMm > 0) || !(ref.drumDiameterMm > 0) || !(ref.ballDiameterMm > 0)) {
      matchOut.textContent = "-";
      return;
    }
    try {
      const d = await getClient().matchingBallDiameter(ref.ballDiameterMm / 1000, ref.drumDiameterMm / 1000, matchTargetMm / 1000);
      matchedD = Number((d * 1000).toPrecision(4));
      matchOut.textContent = `ball diameter ${matchedD} mm (D/d = ${formatNumber(matchTargetMm / matchedD, 3)}, ${Math.round(ballCountAtRatio(ref.fillFraction, matchTargetMm / matchedD, packing))} balls)`;
      matchAdd.disabled = false;
    } catch (err) {
      matchOut.textContent = err instanceof Error ? err.message : String(err);
    }
  }
  matchRow.append(
    labelled("Target drum diameter (mm)", numInput(matchTargetMm, "10", (v) => {
      matchTargetMm = v;
      void updateMatch();
    }, 0)),
    button("Compute", () => void updateMatch()),
    matchOut,
    matchAdd,
  );
  matchSection.append(
    el("p", "su-note", "Keeps the drop-type stress energy SE ~ rho d^3 g D of the reference at the same % critical speed: d = d_ref (D_ref / D_target)^(1/3). A starting point for a candidate, not a result."),
    matchRow,
  );

  // ---- D/d sweep ---------------------------------------------------------------------------
  const sweepSection = el("section", "su-section");
  sweepSection.append(el("h3", "su-h3", "D/d sweep"));
  const sweepBaseSel = el("select", "su-select");
  let sweepRatiosText = "20,30,45,60,90";
  const sweepRatiosIn = el("input", "su-text");
  sweepRatiosIn.value = sweepRatiosText;
  sweepRatiosIn.addEventListener("input", () => {
    sweepRatiosText = sweepRatiosIn.value;
  });
  const sweepRunBtn = button("Run sweep", () => void runSweep(), "su-btn su-primary");
  const sweepMsg = el("div", "su-msg");
  const sweepOut = el("div");
  sweepSection.append(
    el("p", "su-note", "Holds the base condition's ball diameter and every dimensionless setting fixed and varies the drum diameter (D = ratio x d) at constant % critical. Fitted power laws allow extrapolating to a larger D/d without simulating its ball count."),
    (() => {
      const row = el("div", "su-row su-wrap");
      row.append(labelled("Base condition", sweepBaseSel), labelled("D/d ratios (comma separated)", sweepRatiosIn), sweepRunBtn);
      return row;
    })(),
    sweepMsg,
    sweepOut,
  );

  function parseRatios(text: string): number[] {
    return text
      .split(/[\s,;]+/)
      .filter((s) => s !== "")
      .map(Number);
  }

  async function runSweep(): Promise<void> {
    sweepMsg.textContent = "";
    sweepMsg.classList.remove("su-error");
    const cond = conditions[Number(sweepBaseSel.value) || 0];
    if (!cond) return;
    const ratios = parseRatios(sweepRatiosText);
    if (ratios.length < 2 || ratios.some((x) => !Number.isFinite(x) || x <= 0)) {
      sweepMsg.textContent = "Enter at least two positive ratios.";
      sweepMsg.classList.add("su-error");
      return;
    }
    const bad = ratios.filter((x) => ballCountAtRatio(cond.fillFraction, x, packing) > MAX_BALLS || ballCountAtRatio(cond.fillFraction, x, packing) < MIN_BALLS);
    if (bad.length) {
      sweepMsg.textContent = `D/d = ${bad.join(", ")} would give fewer than ${MIN_BALLS} or more than ${MAX_BALLS} balls (J = ${cond.fillFraction}).`;
      sweepMsg.classList.add("su-error");
      return;
    }
    setBusy(true, `sweep of ${ratios.length} ratios x ${cfg.seedCount} seeds (no finer progress available)`);
    try {
      const base = await getBase();
      const params = buildParams(base, cond);
      const config = { ...cfg };
      const text = await getClient().sweep(JSON.stringify(params), JSON.stringify(buildScaleupConfig(config, seedList(config.seedCount))), ratios);
      sweepResult = { baseName: cond.name, ratios, data: JSON.parse(text) as SweepResult, config, params };
      sweepMsg.textContent = "Done.";
      renderSweep();
    } catch (err) {
      reportError(sweepMsg, err, "Sweep failed");
    } finally {
      setBusy(false);
    }
  }

  let extrapMetric = "";
  let extrapTarget = 200;

  function renderSweep(): void {
    sweepOut.replaceChildren();
    if (!sweepResult) return;
    const { data, ratios } = sweepResult;
    const fitNames = data.fits.map((f) => f.name);
    const pts = el("table", "su-table");
    const h = el("tr");
    h.append(el("th", undefined, "D/d"), el("th", undefined, "Balls"));
    for (const n of fitNames) h.append(el("th", undefined, metricLabelForKey(n)));
    const th = el("thead");
    th.append(h);
    const tb = el("tbody");
    for (const p of data.points) {
      const tr = el("tr");
      tr.append(el("td", "su-val", formatNumber(p.dd_ratio, 4)), el("td", "su-val", String(p.n_balls)));
      for (const n of fitNames) tr.append(el("td", "su-val", fmtMeanSe(getMetric(p.ensemble, n))));
      tb.append(tr);
    }
    pts.append(th, tb);

    const fits = el("table", "su-table");
    const fh = el("tr");
    for (const t of ["Metric", "Exponent", "Expected", "R^2", "Prefactor"]) fh.append(el("th", undefined, t));
    const fth = el("thead");
    fth.append(fh);
    const ftb = el("tbody");
    for (const f of data.fits) {
      const tr = el("tr");
      tr.append(
        el("td", undefined, metricLabelForKey(f.name)),
        el("td", "su-val", f.exponent.toFixed(3)),
        el("td", "su-val", f.expected_exponent === null ? "-" : String(f.expected_exponent)),
        el("td", f.r2 < 0.95 ? "su-val su-noisy" : "su-val", f.r2.toFixed(4)),
        el("td", "su-val", formatNumber(f.prefactor, 4)),
      );
      ftb.append(tr);
    }
    fits.append(fth, ftb);

    // extrapolation calculator
    const calc = el("div", "su-row su-wrap");
    const sel = el("select", "su-select");
    for (const f of data.fits) {
      const o = el("option", undefined, metricLabelForKey(f.name));
      o.value = f.name;
      sel.append(o);
    }
    if (!data.fits.some((f) => f.name === extrapMetric)) extrapMetric = data.fits[0]?.name ?? "";
    sel.value = extrapMetric;
    const out = el("span", "su-derived");
    const warn = el("div");
    const update = (): void => {
      warn.replaceChildren();
      const fit = data.fits.find((f) => f.name === extrapMetric);
      if (!fit) {
        out.textContent = "-";
        return;
      }
      const ex = extrapolate(fit, extrapTarget, ratios);
      out.textContent = Number.isFinite(ex.value) ? `predicted ${formatNumber(ex.value, 4)}` : "-";
      for (const w of ex.warnings) warn.append(el("div", "su-warning", w));
    };
    sel.addEventListener("change", () => {
      extrapMetric = sel.value;
      update();
    });
    calc.append(
      labelled("Metric", sel),
      labelled("Target D/d", numInput(extrapTarget, "10", (v) => {
        extrapTarget = v;
        update();
      }, 0)),
      out,
    );
    update();
    sweepOut.append(el("h4", "su-h4", "Sweep points"), wrapScroll(pts), el("h4", "su-h4", "Power-law fits (y = prefactor x D/d^exponent)"), wrapScroll(fits), el("h4", "su-h4", "Extrapolation"), calc, warn);
  }

  function wrapScroll(node: HTMLElement): HTMLElement {
    const w = el("div", "su-scroll");
    w.append(node);
    return w;
  }

  // ---- sensitivity ---------------------------------------------------------------------------
  const sensSection = el("section", "su-section");
  sensSection.append(el("h3", "su-h3", "Contact sensitivity"));
  const sensRunBtn = button("Run sensitivity (+-20 %)", () => void runSensitivity(), "su-btn su-primary");
  const sensMsg = el("div", "su-msg");
  const sensOut = el("div");
  sensSection.append(
    el("p", "su-note", "Re-runs the reference condition with restitution and friction (ball-ball, ball-wall) scaled by +-20 % to show how much each metric depends on the assumed contact values. Changes above 10 % are highlighted: ratios between conditions are only trustworthy where this is small."),
    (() => {
      const row = el("div", "su-row su-wrap");
      row.append(sensRunBtn);
      return row;
    })(),
    sensMsg,
    sensOut,
  );

  async function runSensitivity(): Promise<void> {
    sensMsg.textContent = "";
    sensMsg.classList.remove("su-error");
    const cond = conditions[0];
    if (!cond) return;
    const rel = 0.2;
    setBusy(true, `sensitivity of "${cond.name}": 9 runs x ${cfg.seedCount} seeds (no finer progress available)`);
    try {
      const base = await getBase();
      const params = buildParams(base, cond);
      const text = await getClient().sensitivity(JSON.stringify(params), JSON.stringify(buildScaleupConfig(cfg, seedList(cfg.seedCount))), rel);
      sensResult = { name: cond.name, rel, cases: JSON.parse(text) as SensitivityCase[], params };
      sensMsg.textContent = "Done.";
      renderSensitivity();
    } catch (err) {
      reportError(sensMsg, err, "Sensitivity failed");
    } finally {
      setBusy(false);
    }
  }

  function renderSensitivity(): void {
    sensOut.replaceChildren();
    if (!sensResult) return;
    const [base, ...rest] = sensResult.cases;
    if (!base) return;
    const table = el("table", "su-table");
    const h = el("tr");
    h.append(el("th", undefined, "Metric (change vs base)"));
    for (const c of rest) h.append(el("th", undefined, c.label));
    const th = el("thead");
    th.append(h);
    const tb = el("tbody");
    for (const key of primaryMetricKeys(base.scalar_names)) {
      const tr = el("tr");
      tr.append(el("td", undefined, metricLabelForKey(key)));
      for (const c of rest) {
        const pct = sensitivityChange(base, c, key);
        const td = el("td", isSensitive(pct) ? "su-val su-noisy" : "su-val", Number.isFinite(pct) ? `${pct >= 0 ? "+" : ""}${pct.toFixed(1)} %` : "-");
        tr.append(td);
      }
      tb.append(tr);
    }
    table.append(th, tb);
    sensOut.append(wrapScroll(table));
  }

  // ---- selects that list condition names -----------------------------------------------------
  function refreshSelects(): void {
    const prev = sweepBaseSel.value;
    sweepBaseSel.replaceChildren();
    conditions.forEach((c, i) => {
      const o = el("option", undefined, c.name || `condition ${i + 1}`);
      o.value = String(i);
      sweepBaseSel.append(o);
    });
    if (prev !== "" && Number(prev) < conditions.length) sweepBaseSel.value = prev;
  }

  // ---- assemble ---------------------------------------------------------------------------
  const body = el("div", "su-body");
  body.append(intro, condSection, cfgSection, runBar, runMsg, resultsSection, matchSection, sweepSection, sensSection);
  root.append(header, body);
  renderConditions();
  window.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && !root.hidden) {
      setVisible(false);
      onClosed();
    }
  });

  function setVisible(visible: boolean): void {
    root.hidden = !visible;
    if (visible) {
      getClient();
      requestAnimationFrame(drawCharts);
      void updateMatch();
    }
  }

  return { el: root, isVisible: () => !root.hidden, setVisible };
}
