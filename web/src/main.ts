import { CanvasRenderer } from "./render/canvas";
import type { FrameMessage, MainToWorkerMessage, ParamsJson, WorkerToMainMessage } from "./protocol";
import { createInitialState } from "./state";

const state = createInitialState();

const canvas = document.querySelector<HTMLCanvasElement>("#scene");
if (!canvas) {
  throw new Error("Missing #scene canvas element");
}
const renderer = new CanvasRenderer(canvas);

function resizeToViewport(): void {
  renderer.resize(window.innerWidth, window.innerHeight);
}
resizeToViewport();
window.addEventListener("resize", resizeToViewport);

const worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });

function send(message: MainToWorkerMessage): void {
  worker.postMessage(message);
}

/** Reads `mill.diameter_m` out of the (currently opaque) params blob; see protocol.ts. */
function radiusMFromParams(params: ParamsJson | null): number {
  const mill = params?.mill as { diameter_m?: number } | undefined;
  return (mill?.diameter_m ?? 1.0) / 2;
}

function applyFrame(msg: FrameMessage): void {
  state.drumAngle = msg.drumAngle;
  state.simTime = msg.simTime;
  state.achievedTimeScale = msg.achievedTimeScale;
  state.ballPositions = msg.ballPositions;
  state.ballOrientations = msg.ballOrientations;
  state.ballRadiusM = msg.ballRadiusM;
}

worker.onmessage = (event: MessageEvent<WorkerToMainMessage>) => {
  const msg = event.data;
  switch (msg.type) {
    case "ready":
      state.params = msg.params;
      break;
    case "frame":
      applyFrame(msg);
      break;
    case "error":
      console.error("[mill-worker]", msg.message);
      break;
  }
};

send({ type: "init" });

let lastTimeMs: number | null = null;
function frameLoop(nowMs: number): void {
  if (lastTimeMs !== null) {
    const wallDt = (nowMs - lastTimeMs) / 1000;
    if (state.running) {
      send({ type: "requestFrame", wallDt });
    }
  }
  lastTimeMs = nowMs;
  renderer.render({
    radiusM: radiusMFromParams(state.params),
    drumAngle: state.drumAngle,
    ballPositions: state.ballPositions,
    ballOrientations: state.ballOrientations,
    ballRadiusM: state.ballRadiusM,
  });
  requestAnimationFrame(frameLoop);
}
requestAnimationFrame(frameLoop);
