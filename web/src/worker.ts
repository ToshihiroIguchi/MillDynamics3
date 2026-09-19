// Owns the wasm `Simulation` instance and drives it from `requestFrame` messages sent by the
// main thread's rAF loop (see src/main.ts and docs/PLAN.md ss4.1).
//
// This file intentionally avoids the "webworker" TS lib (which conflicts with "DOM", used by
// src/main.ts, in a single tsconfig program): `self` is narrowed to exactly the shape used below
// instead of relying on ambient `DedicatedWorkerGlobalScope` types.
//
// Fixed-sub-step accumulator (docs/PLAN.md ss4.1): `Simulation::fixedSubDt()` is a *constant*
// sub-step size (`1 / (60 * simulation.substeps)`), independent of both the actual wall-clock
// frame rate and `simulation.time_scale`. This replaced an earlier, buggy design where the
// sub-step size was `wallDt * timeScale` -- i.e. a slow frame or a high `time_scale` produced a
// larger sub-step, which is exactly the kind of oversized step that risks DEM tunnelling (a ball
// moving further than its own radius within one sub-step and skipping a contact). Now `time_scale`
// only controls how much simulated time we try to *catch up* per unit of wall-clock time; the
// sub-step size itself never changes. `pendingSimTime` below is the accumulator that reconciles
// the two: each wall-clock frame adds `wallDt * timeScale` sim-seconds to it, and we drain it in
// fixed `subDt`-sized chunks (bounded by `frame_budget_ms` of wall-clock work per call) via
// `sim.stepFixed()`.

import init, { Simulation, setPanicHook } from "./wasm/mill_wasm.js";
import type { Metrics } from "./metrics/types";
import type { FrameMessage, MainToWorkerMessage, ParamsJson, WorkerToMainMessage } from "./protocol";

interface WorkerScope {
  onmessage: ((event: MessageEvent<MainToWorkerMessage>) => void) | null;
  postMessage: (message: WorkerToMainMessage, transfer: Transferable[]) => void;
}
const scope = self as unknown as WorkerScope;

let sim: Simulation | null = null;
let running = true;
/** `simulation.time_scale` from the current params (params.rs, default 1.0): how much simulated
 * time to advance per unit of wall-clock time (i.e. how many fixed sub-steps to run to catch up),
 * NOT the size of each sub-step -- that's `sim.fixedSubDt()`, which is constant. Set in `boot()`. */
let timeScale = 1;
/** `simulation.frame_budget_ms` from the current params (params.rs, default 12): the wall-clock
 * time budget (ms) `drainPendingSimTime` may spend running `stepFixed()` calls within a single
 * `requestFrame`/`step` message, so a backlog of pending sim time is drained gradually across
 * several rendered frames rather than blocking one frame for as long as it takes to fully catch
 * up. Set in `boot()`. */
let frameBudgetMs = 12;
/** Accumulated simulated time (seconds) not yet advanced via `stepFixed()`. Grows by
 * `wallDt * timeScale` per rendered frame and drains in `sim.fixedSubDt()`-sized chunks; see the
 * module comment above and `drainPendingSimTime` below. */
let pendingSimTime = 0;

/** Largest wall-clock dt accepted from a single requestFrame (e.g. after a backgrounded tab
 * resumes), so a long stall doesn't get folded into `pendingSimTime` as one huge catch-up demand
 * that would otherwise take many frames (or, before the frame-budget cap below, one very long
 * frame) to drain. */
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
    metrics: JSON.parse(sim.metricsJson()) as Metrics,
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

/**
 * Drains `pendingSimTime` in fixed `subDt`-sized chunks via `sim.stepFixed()`, resetting the
 * per-frame diagnostics first (`resetFrameStats()`) since `stepFixed()` itself doesn't reset them.
 * Bounded by `frameBudgetMs` of wall-clock time spent in *this* call, checked once per sub-step so
 * a slow machine can't blow past the budget by more than one sub-step's worth of work.
 *
 * Also guards against `pendingSimTime` growing unboundedly when the sim is chronically slower than
 * real time (e.g. a very high `time_scale`, or a machine too slow to drain even one frame's worth
 * of backlog within `frame_budget_ms`): anything left over past a small cap is dropped rather than
 * queued up, so a stretch of slowness can't later dump a huge, destabilizing "catch-up" burst once
 * conditions improve. The cap is a handful of sub-steps (or 1 sim-second, whichever is larger, so
 * it doesn't fire spuriously when `substeps` is large and `subDt` tiny).
 *
 * Returns the number of `stepFixed()` calls actually made, which callers use to compute a
 * genuinely *measured* `achievedTimeScale` (it can now read below `timeScale` when the frame
 * budget or catch-up couldn't keep up -- unlike the old code, where it was algebraically pinned to
 * `timeScale`).
 */
function drainPendingSimTime(sim: Simulation, subDt: number): number {
  sim.resetFrameStats();
  const budgetStart = performance.now();
  let stepsRun = 0;
  while (pendingSimTime >= subDt) {
    sim.stepFixed();
    pendingSimTime -= subDt;
    stepsRun += 1;
    if (performance.now() - budgetStart >= frameBudgetMs) break;
  }
  const runawayCap = Math.max(subDt * 8, 1);
  if (pendingSimTime > runawayCap) {
    pendingSimTime = 0;
  }
  return stepsRun;
}

async function boot(initialParams?: ParamsJson): Promise<void> {
  try {
    await init();
    setPanicHook();
    const json = initialParams ? JSON.stringify(initialParams) : undefined;
    sim = new Simulation(json);
    const params = JSON.parse(sim.paramsJson()) as ParamsJson;
    const simulation = params.simulation as { time_scale?: number; frame_budget_ms?: number } | undefined;
    timeScale = simulation?.time_scale ?? 1;
    frameBudgetMs = simulation?.frame_budget_ms ?? 12;
    pendingSimTime = 0;
    post({ type: "ready", params });
    // Send an immediate snapshot so the drum/balls/slurry reflect the new params right away, even
    // while paused -- without this the canvas keeps showing the previous sim's frame until Play or
    // Step is pressed.
    postFrame(sim, 1);
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
    case "step": {
      if (sim) {
        // Manual single-frame-advance: same catch-up mechanics as requestFrame below, for exactly
        // one nominal 60 Hz frame's worth of sim time, so the Step button and the running loop
        // behave identically.
        const nominalWallDt = 1 / 60;
        pendingSimTime += nominalWallDt * timeScale;
        const subDt = sim.fixedSubDt();
        const stepsRun = drainPendingSimTime(sim, subDt);
        const achievedTimeScale = (stepsRun * subDt) / nominalWallDt;
        postFrame(sim, achievedTimeScale);
      }
      break;
    }
    case "requestFrame": {
      if (!sim || !running) return;
      // Clamp so a long stall (e.g. a backgrounded tab) isn't folded into pendingSimTime as one
      // huge catch-up demand.
      const wallDt = Math.min(msg.wallDt, MAX_FRAME_DT);
      pendingSimTime += wallDt * timeScale;
      const subDt = sim.fixedSubDt();
      const stepsRun = drainPendingSimTime(sim, subDt);
      const achievedTimeScale = msg.wallDt > 0 ? (stepsRun * subDt) / msg.wallDt : timeScale;
      postFrame(sim, achievedTimeScale);
      break;
    }
  }
};
