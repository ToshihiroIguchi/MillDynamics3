import { CanvasRenderer, type LiftersRenderState } from "./render/canvas";
import type { FrameMessage, MainToWorkerMessage, ParamsJson, WorkerToMainMessage } from "./protocol";
import { createInitialState } from "./state";
import { createHud } from "./ui/hud";
import { createMetricsPanel } from "./ui/metricsPanel";
import { createParamsModal } from "./ui/paramsModal";
import { createToolbar } from "./ui/toolbar";

const state = createInitialState();

const canvas = document.querySelector<HTMLCanvasElement>("#scene");
if (!canvas) {
  throw new Error("Missing #scene canvas element");
}
const sceneWrapOrNull = document.querySelector<HTMLElement>(".scene-wrap");
if (!sceneWrapOrNull) {
  throw new Error("Missing .scene-wrap element");
}
const sceneWrap: HTMLElement = sceneWrapOrNull;
const renderer = new CanvasRenderer(canvas);

function resizeToWrap(): void {
  const { clientWidth: w, clientHeight: h } = sceneWrap;
  if (w <= 0 || h <= 0) return;
  renderer.resize(w, h);
}
resizeToWrap();
new ResizeObserver(resizeToWrap).observe(sceneWrap);

const worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });

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
  state.ballPositions = msg.ballPositions;
  state.ballOrientations = msg.ballOrientations;
  state.ballRadiusM = msg.ballRadiusM;
  state.fluidPositions = msg.fluidPositions;
  state.fluidDye = msg.fluidDye;
  state.fluidSurface = msg.fluidSurface;
  state.metrics = msg.metrics;
}

worker.onmessage = (event: MessageEvent<WorkerToMainMessage>) => {
  const msg = event.data;
  switch (msg.type) {
    case "ready":
      state.params = msg.params;
      metricsPanel.reset();
      break;
    case "frame":
      applyFrame(msg);
      break;
    case "error":
      console.error("[mill-worker]", msg.message);
      break;
  }
};

function togglePlay(): boolean {
  state.running = !state.running;
  send({ type: state.running ? "play" : "pause" });
  return state.running;
}

// v1 simplification (see params/schema.ts): both "Apply" and "Reset" fully reconstruct the
// simulation via an "init" message; there is no partial "hot" apply yet (M5).
const paramsModal = createParamsModal((params) => {
  send({ type: "init", params });
});

const PANEL_STORAGE_KEY = "milldynamics.panel";

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

const toolbarEl = createToolbar({
  onTogglePlay: togglePlay,
  onStep: () => send({ type: "step" }),
  onReset: () => send({ type: "init", params: state.params ?? undefined }),
  onOpenParams: () => paramsModal.open(state.params ?? {}),
  onTogglePanel: togglePanel,
  initialPanelVisible,
});
sceneWrap.appendChild(toolbarEl);

const hud = createHud();
sceneWrap.appendChild(hud.el);

const panelEl = document.querySelector<HTMLElement>("#metrics-panel");
if (!panelEl) {
  throw new Error("Missing #metrics-panel element");
}
const metricsPanel = createMetricsPanel(panelEl);

window.addEventListener("keydown", (event) => {
  if (event.code === "Space" && !(event.target instanceof HTMLElement && event.target.closest("dialog"))) {
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
      send({ type: "requestFrame", wallDt });
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
  hud.update(state, smoothedFps);
  metricsPanel.update(state, smoothedFps, nowMs);
  requestAnimationFrame(frameLoop);
}
requestAnimationFrame(frameLoop);
