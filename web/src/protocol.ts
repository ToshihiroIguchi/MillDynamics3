// Typed messages exchanged between the main thread and the simulation worker.
// See docs/PLAN.md ss4.1 for the runtime architecture this implements.

import type { Metrics } from "./metrics/types";

/**
 * Parameters as they cross the worker boundary: JSON-serializable, matching mill-core's `Params`
 * (see crates/mill-core/src/params.rs). The web app treats this as an opaque blob through M0;
 * typed field access lands with the parameters modal and params/schema.ts in M1+.
 */
export type ParamsJson = Record<string, unknown>;

export interface InitMessage {
  type: "init";
  /** Params to start with, or omitted to use mill-core's defaults. */
  params?: ParamsJson;
}

export interface SetParamsMessage {
  type: "setParams";
  /**
   * Params to hot-apply in place (`Simulation::set_params`, no reset of drum angle/sim
   * time/ball/fluid population). The sender (main.ts, via `params/schema.ts`'s
   * `paramsChangeRequiresReset`) is responsible for routing a change that actually needs a reset
   * through `InitMessage` instead -- by the time a message reaches here it is assumed safe to
   * apply live.
   */
  params: ParamsJson;
}

export interface PlayMessage {
  type: "play";
}

export interface PauseMessage {
  type: "pause";
}

export interface StepMessage {
  type: "step";
}

export interface RequestFrameMessage {
  type: "requestFrame";
  /** Wall-clock time elapsed since the previous requestFrame was *sent*, in seconds (not since the
   * previous rAF tick -- those can now differ because a tick may be skipped when at the in-flight
   * cap; see frameGate.ts). */
  wallDt: number;
  /** Identifies this request so the worker's reply (`frame` or `frameSkipped`) can be matched back
   * to it by `frameGate.ts`'s `FrameRequestGate`. */
  requestId: number;
}

export type MainToWorkerMessage =
  | InitMessage
  | SetParamsMessage
  | PlayMessage
  | PauseMessage
  | StepMessage
  | RequestFrameMessage;

export interface ReadyMessage {
  type: "ready";
  params: ParamsJson;
}

export interface FrameMessage {
  type: "frame";
  drumAngle: number;
  simTime: number;
  /** Wall-clock-time-to-sim-time ratio actually achieved this frame (see docs/PLAN.md ss4.1). */
  achievedTimeScale: number;
  /**
   * Fixed sub-steps per second this frame actually ran (`stepsRun / wallDt`), and how many a full
   * `time_scale = 1` real-time rate needs (`1 / fixedSubDt()`, a constant for the current params).
   * `achievedTimeScale` is exactly their ratio; these two are surfaced separately so the HUD can
   * show *why* it's short (e.g. "62 / 480 sub-steps/s"), not just the ratio.
   */
  subStepsPerSecondAchieved: number;
  subStepsPerSecondRequired: number;
  /** Ball center positions, flattened as [x0, y0, x1, y1, ...] (m), transferred (not copied). */
  ballPositions: Float32Array;
  /** Ball orientations (radians), one per ball, same order as ballPositions. */
  ballOrientations: Float32Array;
  ballRadiusM: number;
  /** Fluid (slurry) particle positions, flattened as [x0, y0, x1, y1, ...] (m). */
  fluidPositions: Float32Array;
  /** Fluid dye tracer values ([0, 1]), one per particle, same order as fluidPositions. */
  fluidDye: Float32Array;
  /**
   * Free-surface contour(s), flattened as [n_polys, len_0, x, y, ..., len_1, ...] (docs/PLAN.md
   * ss3.5). `undefined` when this frame did not recompute it -- marching squares over a 128x128
   * grid is expensive enough that the worker throttles it well below render rate (see worker.ts's
   * `SLOW_UPDATE_INTERVAL_MS`); the consumer should keep showing the last value it received.
   */
  fluidSurface?: Float32Array;
  /**
   * Derived metrics (toe/shoulder, slurry pool extent, mixing index, grinding/solver diagnostics,
   * debug checks, docs/PLAN.md ss3.5), parsed once by the worker (see worker.ts's `buildFrame`)
   * from `Simulation::metrics_json()`'s JSON text into the typed shape in metrics/types.ts.
   * `undefined` when this frame did not recompute it, same throttling and same "keep the last
   * value" contract as `fluidSurface`.
   */
  metrics?: Metrics;
  /** Set only when this frame answers a gated `requestFrame` message (not set for the `step`
   * message's frame, nor for the initial snapshot frame posted at the end of `boot()`); see
   * frameGate.ts's `FrameRequestGate.settle`. */
  requestId?: number;
}

export interface ErrorMessage {
  type: "error";
  message: string;
}

/** Reply to a `requestFrame` that could not be answered with a `frame` (sim not ready/running, or
 * an internal error while stepping/building the frame) -- still frees the request's in-flight slot
 * in `frameGate.ts`'s `FrameRequestGate` so the gate can never deadlock waiting for a reply that
 * would otherwise never come. */
export interface FrameSkippedMessage {
  type: "frameSkipped";
  requestId: number;
}

export type WorkerToMainMessage = ReadyMessage | FrameMessage | ErrorMessage | FrameSkippedMessage;
