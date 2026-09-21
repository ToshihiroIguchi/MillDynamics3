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

import init, { defaultParamsJson, Simulation, setPanicHook } from "./wasm/mill_wasm.js";
import type { Metrics } from "./metrics/types";
import { QUALITY_PRESETS } from "./params/presets";
import type { FrameMessage, MainToWorkerMessage, ParamsJson, WorkerToMainMessage } from "./protocol";
import { withPath } from "./params/schema";

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

/**
 * Minimum wall-clock interval (ms) between `fluidSurface()`/`metricsJson()` recomputations --
 * both are expensive relative to a render frame (marching squares over a 128x128 grid; a full
 * toe/shoulder/pool/mixing-index pass over every ball and fluid particle) and neither needs to be
 * fresher than this to look smooth or read correctly: the free-surface *outline* barely moves
 * frame-to-frame at typical sim speeds, and the metrics panel already throttles its own DOM paint
 * to 100 ms (ui/metricsPanel.ts). Before this, both ran unconditionally every rendered frame --
 * one contributor to the achieved-time-scale shortfall the HUD reports (ui/hud.ts).
 */
const SLOW_UPDATE_INTERVAL_MS = 1000 / 15;
let lastSlowUpdateMs = Number.NEGATIVE_INFINITY;

function buildFrame(
  sim: Simulation,
  achievedTimeScale: number,
  subStepsPerSecondAchieved: number,
  subStepsPerSecondRequired: number,
): FrameMessage {
  const frame: FrameMessage = {
    type: "frame",
    drumAngle: sim.drumAngle(),
    simTime: sim.simTime(),
    achievedTimeScale,
    subStepsPerSecondAchieved,
    subStepsPerSecondRequired,
    ballPositions: sim.ballPositions(),
    ballOrientations: sim.ballOrientations(),
    ballRadiusM: sim.ballRadiusM(),
    fluidPositions: sim.fluidPositions(),
    fluidDye: sim.fluidDye(),
  };
  const now = performance.now();
  if (now - lastSlowUpdateMs >= SLOW_UPDATE_INTERVAL_MS) {
    lastSlowUpdateMs = now;
    frame.fluidSurface = sim.fluidSurface();
    frame.metrics = JSON.parse(sim.metricsJson()) as Metrics;
  }
  return frame;
}

function postFrame(
  sim: Simulation,
  achievedTimeScale: number,
  subStepsPerSecondAchieved: number,
  subStepsPerSecondRequired: number,
): void {
  const frame = buildFrame(sim, achievedTimeScale, subStepsPerSecondAchieved, subStepsPerSecondRequired);
  const transfer: Transferable[] = [
    frame.ballPositions.buffer,
    frame.ballOrientations.buffer,
    frame.fluidPositions.buffer,
    frame.fluidDye.buffer,
  ];
  if (frame.fluidSurface) transfer.push(frame.fluidSurface.buffer);
  post(frame, transfer);
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

/** Re-reads the two params fields the worker itself caches in plain JS locals (not read through
 * wasm on every `stepFixed()` the way the Rust-side solver reads its own params) -- needed after
 * both `boot()` and a hot `setParams`, since the latter would otherwise leave a stale `timeScale`/
 * `frameBudgetMs` driving `drainPendingSimTime`'s catch-up math even though `sim`'s own params.rs
 * state already reflects the change. */
function syncCachedParams(params: ParamsJson): void {
  const simulation = params.simulation as { time_scale?: number; frame_budget_ms?: number } | undefined;
  timeScale = simulation?.time_scale ?? 1;
  frameBudgetMs = simulation?.frame_budget_ms ?? 12;
}

/** id of the preset (params/presets.ts) a fresh app load with no explicit params starts at -- see
 * `boot()`'s use of it below. Not mill-core's own `SimulationParams::default()`
 * (`max_balls = 600`, `resolution = 40`): that Default is deliberately left at this project's
 * original, best-fidelity values since native tests/benches (docs/PERF.md) compare against it
 * directly, but it achieves well under 1x real time (see docs/PERF.md, presets.ts's own doc
 * comment) -- a poor first impression for a fresh page load with nothing else to compare against.
 * "accuracy" here would reproduce that Default's own values; "realtime" is deliberately not that.
 */
const INITIAL_PRESET_ID = "realtime";

async function boot(initialParams?: ParamsJson): Promise<void> {
  try {
    await init();
    setPanicHook();
    let params = initialParams;
    if (!params) {
      const preset = QUALITY_PRESETS.find((p) => p.id === INITIAL_PRESET_ID);
      const rustDefaults = JSON.parse(defaultParamsJson()) as ParamsJson;
      params = preset
        ? withPath(withPath(rustDefaults, "simulation.max_balls", preset.maxBalls), "simulation.resolution", preset.resolution)
        : rustDefaults;
    }
    const json = JSON.stringify(params);
    sim = new Simulation(json);
    const effectiveParams = JSON.parse(sim.paramsJson()) as ParamsJson;
    syncCachedParams(effectiveParams);
    pendingSimTime = 0;
    // Force the very next `buildFrame` to recompute fluidSurface/metrics rather than skipping
    // them per `SLOW_UPDATE_INTERVAL_MS` -- otherwise a reset that lands within that window of
    // the previous simulation's last update would send `undefined` for both, leaving the old
    // (now-replaced) simulation's stale surface/metrics on screen.
    lastSlowUpdateMs = Number.NEGATIVE_INFINITY;
    post({ type: "ready", params: effectiveParams });
    // Send an immediate snapshot so the drum/balls/slurry reflect the new params right away, even
    // while paused -- without this the canvas keeps showing the previous sim's frame until Play or
    // Step is pressed. No sub-steps actually ran, so report 0 achieved against the real
    // requirement rather than claiming a fictitious 1x.
    postFrame(sim, 1, 0, 1 / sim.fixedSubDt());
  } catch (err) {
    // Leaves `sim` as `null`, but this is not a permanent freeze: "init" messages are not gated
    // on `sim` being non-null, so the toolbar's Reset button (main.ts, wired to `{ type: "init",
    // params: state.params ?? undefined }`) retries `boot()` and can recover. Say so, since the
    // banner is otherwise the only thing the user sees and nothing else tells them Reset helps.
    post({ type: "error", message: `${String(err)} Press Reset to retry.` });
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
          syncCachedParams(msg.params);
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
        postFrame(sim, achievedTimeScale, stepsRun / nominalWallDt, 1 / subDt);
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
      const achievedSubStepsPerSecond = msg.wallDt > 0 ? stepsRun / msg.wallDt : 0;
      postFrame(sim, achievedTimeScale, achievedSubStepsPerSecond, 1 / subDt);
      break;
    }
  }
};
