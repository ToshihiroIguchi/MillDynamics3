// Typed messages exchanged between the main thread and the simulation worker.
// See docs/PLAN.md ss4.1 for the runtime architecture this implements.

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
  params: ParamsJson;
  /** Whether this requires a full simulation reset (true) or can be applied live (false). */
  reset: boolean;
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
  /** Wall-clock time elapsed since the previous requestFrame, in seconds. */
  wallDt: number;
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
  /** Ball center positions, flattened as [x0, y0, x1, y1, ...] (m), transferred (not copied). */
  ballPositions: Float32Array;
  /** Ball orientations (radians), one per ball, same order as ballPositions. */
  ballOrientations: Float32Array;
  ballRadiusM: number;
  /** Fluid (slurry) particle positions, flattened as [x0, y0, x1, y1, ...] (m). */
  fluidPositions: Float32Array;
  /** Fluid dye tracer values ([0, 1]), one per particle, same order as fluidPositions. */
  fluidDye: Float32Array;
  /** Free-surface contour(s), flattened as [n_polys, len_0, x, y, ..., len_1, ...] (docs/PLAN.md ss3.5). */
  fluidSurface: Float32Array;
}

export interface ErrorMessage {
  type: "error";
  message: string;
}

export type WorkerToMainMessage = ReadyMessage | FrameMessage | ErrorMessage;
