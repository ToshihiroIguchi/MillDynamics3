import { FrameRequestGate } from "./frameGate";
import { CanvasRenderer, type LiftersRenderState } from "./render/canvas";
import type { FrameMessage, MainToWorkerMessage, ParamsJson, WorkerToMainMessage } from "./protocol";
import { paramsChangeRequiresReset } from "./params/schema";
import { createInitialState } from "./state";
import { createHud } from "./ui/hud";
import { createMetricsPanel } from "./ui/metricsPanel";
import { createParamsPanel } from "./ui/paramsPanel";
import { createToolbar } from "./ui/toolbar";
import { MediaRecorderCanvasRecorder } from "./video/recorder";

const state = createInitialState();

const canvasOrNull = document.querySelector<HTMLCanvasElement>("#scene");
if (!canvasOrNull) {
  throw new Error("Missing #scene canvas element");
}
const canvas: HTMLCanvasElement = canvasOrNull;
const sceneWrapOrNull = document.querySelector<HTMLElement>(".scene-wrap");
if (!sceneWrapOrNull) {
  throw new Error("Missing .scene-wrap element");
}
const sceneWrap: HTMLElement = sceneWrapOrNull;
const renderer = new CanvasRenderer(canvas);
const recorder = new MediaRecorderCanvasRecorder();

function resizeToWrap(): void {
  const { clientWidth: w, clientHeight: h } = sceneWrap;
  if (w <= 0 || h <= 0) return;
  renderer.resize(w, h);
}
resizeToWrap();
new ResizeObserver(resizeToWrap).observe(sceneWrap);

const worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
const gate = new FrameRequestGate();

function send(message: MainToWorkerMessage): void {
  worker.postMessage(message);
}

/** Reads `mill.diameter_m` out of the (currently opaque) params blob; see protocol.ts. */
function radiusMFromParams(params: ParamsJson | null): number {
  const mill = params?.mill as { diameter_m?: number } | undefined;
  return (mill?.diameter_m ?? 1.0) / 2;
}

const NO_LIFTERS: LiftersRenderState = { count: 0, heightM: 0, baseWidthM: 0, topWidthM: 0, phaseDeg: 0 };

/** Reads `lifters` out of the (currently opaque) params blob; see protocol.ts. */
function liftersFromParams(params: ParamsJson | null): LiftersRenderState {
  const lifters = params?.lifters as
    | { count?: number; height_m?: number; base_width_m?: number; top_width_m?: number; phase_deg?: number }
    | undefined;
  if (!lifters) return NO_LIFTERS;
  return {
    count: lifters.count ?? 0,
    heightM: lifters.height_m ?? 0,
    baseWidthM: lifters.base_width_m ?? 0,
    topWidthM: lifters.top_width_m ?? 0,
    phaseDeg: lifters.phase_deg ?? 0,
  };
}

function applyFrame(msg: FrameMessage): void {
  state.drumAngle = msg.drumAngle;
  state.simTime = msg.simTime;
  state.achievedTimeScale = msg.achievedTimeScale;
  state.subStepsPerSecondAchieved = msg.subStepsPerSecondAchieved;
  state.subStepsPerSecondRequired = msg.subStepsPerSecondRequired;
  state.ballPositions = msg.ballPositions;
  state.ballOrientations = msg.ballOrientations;
  state.ballRadiusM = msg.ballRadiusM;
  state.fluidPositions = msg.fluidPositions;
  state.fluidDye = msg.fluidDye;
  // `fluidSurface`/`metrics` are throttled well below render rate (worker.ts's
  // `SLOW_UPDATE_INTERVAL_MS`) and arrive as `undefined` on frames that didn't recompute them --
  // keep showing the last received value rather than clearing it.
  if (msg.fluidSurface) state.fluidSurface = msg.fluidSurface;
  if (msg.metrics) state.metrics = msg.metrics;
}

worker.onmessage = (event: MessageEvent<WorkerToMainMessage>) => {
  const msg = event.data;
  switch (msg.type) {
    case "ready":
      state.params = msg.params;
      metricsPanel.reset();
      paramsPanel.setParams(msg.params);
      break;
    case "frame":
      applyFrame(msg);
      if (msg.requestId !== undefined) gate.settle(msg.requestId);
      break;
    case "frameSkipped":
      gate.settle(msg.requestId);
      break;
    case "error":
      console.error("[mill-worker]", msg.message);
      paramsPanel.showError(msg.message);
      break;
  }
};

function togglePlay(): boolean {
  state.running = !state.running;
  send({ type: state.running ? "play" : "pause" });
  if (!state.running) gate.resetClock();
  return state.running;
}

// Apply hot-applies via "setParams" (Simulation::set_params, keeps t/drum angle/ball & fluid
// population) unless the change actually needs a reset (params/schema.ts's
// `paramsChangeRequiresReset`, e.g. drum diameter or media geometry), in which case it falls back
// to "init" like the Reset button always does. `state.params` is updated locally right away on the
// hot path since the worker sends no reply to "setParams" to update it from (unlike "init"'s
// "ready" message) -- the value being applied is already authoritative, constructed by the panel
// from the form the user just submitted.
const paramsPanelEl = document.querySelector<HTMLElement>("#params-panel");
if (!paramsPanelEl) {
  throw new Error("Missing #params-panel element");
}
const paramsPanel = createParamsPanel(paramsPanelEl, (params) => {
  const resetCause = state.params ? paramsChangeRequiresReset(state.params, params) : "initial load";
  if (resetCause) {
    console.info(`[params] resetting simulation: "${resetCause}" changed`);
    gate.resetClock();
    send({ type: "init", params });
  } else {
    console.info("[params] hot-applying (no reset)");
    state.params = params;
    send({ type: "setParams", params });
    // No "ready" reply follows a hot "setParams" (unlike "init"), so the panel's status indicator
    // would otherwise stay stuck on "Edited -- not applied" even though the apply succeeded --
    // restore it to "Applied" here, mirroring the "ready" handler above.
    paramsPanel.setParams(params);
  }
});

const PANEL_STORAGE_KEY = "milldynamics.panel";
const PARAMS_PANEL_STORAGE_KEY = "milldynamics.paramsPanel";

function readPanelVisiblePref(): boolean {
  try {
    const raw = localStorage.getItem(PANEL_STORAGE_KEY);
    if (raw === null) return true;
    return JSON.parse(raw) === true;
  } catch {
    return true;
  }
}

function setPanelCollapsed(collapsed: boolean): void {
  document.body.classList.toggle("panel-collapsed", collapsed);
}

const initialPanelVisible = readPanelVisiblePref();
setPanelCollapsed(!initialPanelVisible);

function togglePanel(): boolean {
  const wasCollapsed = document.body.classList.contains("panel-collapsed");
  const nowVisible = wasCollapsed;
  setPanelCollapsed(!nowVisible);
  try {
    localStorage.setItem(PANEL_STORAGE_KEY, JSON.stringify(nowVisible));
  } catch {
    // ignore (private browsing / disabled storage)
  }
  return nowVisible;
}

function readParamsPanelVisiblePref(): boolean {
  try {
    const raw = localStorage.getItem(PARAMS_PANEL_STORAGE_KEY);
    if (raw === null) return true;
    return JSON.parse(raw) === true;
  } catch {
    return true;
  }
}

function setParamsPanelCollapsed(collapsed: boolean): void {
  document.body.classList.toggle("params-collapsed", collapsed);
}

const initialParamsVisible = readParamsPanelVisiblePref();
setParamsPanelCollapsed(!initialParamsVisible);

function toggleParamsPanel(): boolean {
  const wasCollapsed = document.body.classList.contains("params-collapsed");
  const nowVisible = wasCollapsed;
  setParamsPanelCollapsed(!nowVisible);
  try {
    localStorage.setItem(PARAMS_PANEL_STORAGE_KEY, JSON.stringify(nowVisible));
  } catch {
    // ignore (private browsing / disabled storage)
  }
  return nowVisible;
}

/** Starts/stops video/recorder.ts's `CanvasRecorder` on `#scene`, independent of simulation
 * running/reset state (recording just keeps capturing whatever gets rendered, including a reset
 * transient, mirroring how a real screen recording would behave -- see video/recorder.ts's doc
 * comment). `stop()`'s promise isn't awaited by the toolbar (its own click handler needs the new
 * state synchronously to swap the icon immediately), so a rejection is caught here instead of
 * being left to become an unhandled rejection. */
function toggleRecord(): boolean {
  if (recorder.isRecording()) {
    state.isRecording = false;
    state.recordingStartSimTime = null;
    void recorder.stop(state.simTime).catch((err) => console.error("[recorder] stop failed", err));
    return false;
  }
  recorder.start(canvas, state.simTime);
  state.isRecording = recorder.isRecording();
  state.recordingStartSimTime = state.isRecording ? state.simTime : null;
  return state.isRecording;
}

const toolbar = createToolbar({
  onTogglePlay: togglePlay,
  onStep: () => send({ type: "step" }),
  onReset: () => {
    gate.resetClock();
    send({ type: "init", params: state.params ?? undefined });
  },
  onToggleParams: toggleParamsPanel,
  initialParamsVisible,
  onTogglePanel: togglePanel,
  initialPanelVisible,
  onToggleRecord: toggleRecord,
  recordingSupported: recorder.isSupported(),
});
sceneWrap.appendChild(toolbar.el);

const hud = createHud();
sceneWrap.appendChild(hud.el);

/** Pauses the simulation programmatically (not via the toolbar button's own click), mirroring
 * `togglePlay`'s "pause" branch exactly -- same `state.running` update, same `send({type:
 * "pause"})`, same `gate.resetClock()` -- but also syncing the toolbar's icon by hand since
 * `toolbar.ts`'s own click handler (the only other place the icon updates) never runs here. Used
 * by the metrics panel's report auto-trigger (ui/metricsPanel.ts's `onAutoReportPause`) so a run
 * that reaches its configured "Report end time" pauses and reflects that in the UI before the PDF
 * downloads. A no-op if the sim is already paused. */
function pauseForAutoReport(): void {
  if (!state.running) return;
  state.running = false;
  send({ type: "pause" });
  gate.resetClock();
  toolbar.setRunning(false);
}

const panelEl = document.querySelector<HTMLElement>("#metrics-panel");
if (!panelEl) {
  throw new Error("Missing #metrics-panel element");
}
const metricsPanel = createMetricsPanel(panelEl, pauseForAutoReport);

const SPACE_GUARD_SELECTOR = "input, select, textarea, summary, button, [contenteditable]";

window.addEventListener("keydown", (event) => {
  if (
    event.code === "Space" &&
    !(event.target instanceof HTMLElement && event.target.closest(SPACE_GUARD_SELECTOR))
  ) {
    event.preventDefault();
    togglePlay();
  }
});

send({ type: "init" });

let lastTimeMs: number | null = null;
let smoothedFps = 60;
function frameLoop(nowMs: number): void {
  if (lastTimeMs !== null) {
    const wallDt = (nowMs - lastTimeMs) / 1000;
    if (wallDt > 0) {
      smoothedFps = smoothedFps * 0.9 + (1 / wallDt) * 0.1;
    }
    if (state.running) {
      const issued = gate.tryIssue(nowMs, wallDt);
      if (issued) {
        send({ type: "requestFrame", wallDt: issued.wallDt, requestId: issued.requestId });
      }
    }
  }
  lastTimeMs = nowMs;
  renderer.render({
    radiusM: radiusMFromParams(state.params),
    drumAngle: state.drumAngle,
    lifters: liftersFromParams(state.params),
    ballPositions: state.ballPositions,
    ballOrientations: state.ballOrientations,
    ballRadiusM: state.ballRadiusM,
    fluidPositions: state.fluidPositions,
    fluidDye: state.fluidDye,
    fluidSurface: state.fluidSurface,
  });
  recorder.captureFrame(); // no-op when not recording, see video/recorder.ts
  hud.update(state, smoothedFps);
  metricsPanel.update(state, smoothedFps, nowMs);
  requestAnimationFrame(frameLoop);
}
requestAnimationFrame(frameLoop);
