// Owns the wasm `Simulation` instance and drives it from `requestFrame` messages sent by the
// main thread's rAF loop (see src/main.ts and docs/PLAN.md ss4.1).
//
// This file intentionally avoids the "webworker" TS lib (which conflicts with "DOM", used by
// src/main.ts, in a single tsconfig program): `self` is narrowed to exactly the shape used below
// instead of relying on ambient `DedicatedWorkerGlobalScope` types.

import init, { Simulation, setPanicHook } from "./wasm/mill_wasm.js";
import type { FrameMessage, MainToWorkerMessage, ParamsJson, WorkerToMainMessage } from "./protocol";

interface WorkerScope {
  onmessage: ((event: MessageEvent<MainToWorkerMessage>) => void) | null;
  postMessage: (message: WorkerToMainMessage, transfer: Transferable[]) => void;
}
const scope = self as unknown as WorkerScope;

let sim: Simulation | null = null;
let running = true;

/** Largest wall-clock dt accepted from a single requestFrame (e.g. after a backgrounded tab
 * resumes), so a long stall doesn't get integrated as one huge, potentially unstable sub-step. */
const MAX_FRAME_DT = 1 / 15;

function post(message: WorkerToMainMessage, transfer: Transferable[] = []): void {
  scope.postMessage(message, transfer);
}

function buildFrame(sim: Simulation, achievedTimeScale: number): FrameMessage {
  return {
    type: "frame",
    drumAngle: sim.drumAngle(),
    simTime: sim.simTime(),
    achievedTimeScale,
    ballPositions: sim.ballPositions(),
    ballOrientations: sim.ballOrientations(),
    ballRadiusM: sim.ballRadiusM(),
    fluidPositions: sim.fluidPositions(),
    fluidDye: sim.fluidDye(),
    fluidSurface: sim.fluidSurface(),
  };
}

function postFrame(sim: Simulation, achievedTimeScale: number): void {
  const frame = buildFrame(sim, achievedTimeScale);
  post(frame, [
    frame.ballPositions.buffer,
    frame.ballOrientations.buffer,
    frame.fluidPositions.buffer,
    frame.fluidDye.buffer,
    frame.fluidSurface.buffer,
  ]);
}

async function boot(initialParams?: ParamsJson): Promise<void> {
  try {
    await init();
    setPanicHook();
    const json = initialParams ? JSON.stringify(initialParams) : undefined;
    sim = new Simulation(json);
    const params = JSON.parse(sim.paramsJson()) as ParamsJson;
    post({ type: "ready", params });
  } catch (err) {
    post({ type: "error", message: String(err) });
  }
}

scope.onmessage = (event) => {
  const msg = event.data;
  switch (msg.type) {
    case "init":
      void boot(msg.params);
      break;
    case "setParams":
      if (sim) {
        try {
          sim.setParams(JSON.stringify(msg.params));
        } catch (err) {
          post({ type: "error", message: String(err) });
        }
      }
      break;
    case "play":
      running = true;
      break;
    case "pause":
      running = false;
      break;
    case "step":
      if (sim) {
        // One "step" advances by one nominal 60 Hz frame's worth of sim time, split into
        // simulation.substeps fixed sub-steps by mill-core (docs/PLAN.md ss3.2/4.1).
        sim.step(1 / 60);
        postFrame(sim, 1);
      }
      break;
    case "requestFrame": {
      if (!sim || !running) return;
      // Clamp so a long stall (e.g. a backgrounded tab) isn't integrated as one huge sub-step.
      // Frame-budget throttling for sustained slow-frame handling lands with PBF (M3/M4).
      const dt = Math.min(msg.wallDt, MAX_FRAME_DT);
      sim.step(dt);
      const achievedTimeScale = msg.wallDt > 0 ? dt / msg.wallDt : 1;
      postFrame(sim, achievedTimeScale);
      break;
    }
  }
};
