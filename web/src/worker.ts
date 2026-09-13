// Owns the wasm `Simulation` instance and drives it from `requestFrame` messages sent by the
// main thread's rAF loop (see src/main.ts and docs/PLAN.md ss4.1).
//
// This file intentionally avoids the "webworker" TS lib (which conflicts with "DOM", used by
// src/main.ts, in a single tsconfig program): `self` is narrowed to exactly the shape used below
// instead of relying on ambient `DedicatedWorkerGlobalScope` types.

import init, { Simulation, setPanicHook } from "./wasm/mill_wasm.js";
import type { MainToWorkerMessage, ParamsJson, WorkerToMainMessage } from "./protocol";

interface WorkerScope {
  onmessage: ((event: MessageEvent<MainToWorkerMessage>) => void) | null;
  postMessage: (message: WorkerToMainMessage) => void;
}
const scope = self as unknown as WorkerScope;

let sim: Simulation | null = null;
let running = true;

function post(message: WorkerToMainMessage): void {
  scope.postMessage(message);
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
        // Through M0, one "step" advances by one nominal 60 Hz frame's worth of sim time; this
        // becomes a fixed-substep (1/240 s) advance once DEM/PBF land (M1/M3).
        sim.step(1 / 60);
        post({
          type: "frame",
          drumAngle: sim.drumAngle(),
          simTime: sim.simTime(),
          achievedTimeScale: 1,
        });
      }
      break;
    case "requestFrame": {
      if (!sim || !running) return;
      // M0: advance simulation time directly by wall-clock time (no fixed-step accumulator or
      // frame-budget throttling yet -- both land with DEM/PBF in M1/M3, see docs/PLAN.md ss4.1).
      sim.step(msg.wallDt);
      post({
        type: "frame",
        drumAngle: sim.drumAngle(),
        simTime: sim.simTime(),
        achievedTimeScale: 1,
      });
      break;
    }
  }
};
