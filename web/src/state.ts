import type { Metrics } from "./metrics/types";
import type { ParamsJson } from "./protocol";

/** App-level state, updated from worker messages and read by the render loop. */
export interface AppState {
  params: ParamsJson | null;
  running: boolean;
  drumAngle: number;
  simTime: number;
  achievedTimeScale: number;
  /** See protocol.ts's `FrameMessage.subStepsPerSecondAchieved`/`subStepsPerSecondRequired`. */
  subStepsPerSecondAchieved: number;
  subStepsPerSecondRequired: number;
  ballPositions: Float32Array;
  ballOrientations: Float32Array;
  ballRadiusM: number;
  fluidPositions: Float32Array;
  fluidDye: Float32Array;
  fluidSurface: Float32Array;
  /** Derived metrics (docs/PLAN.md ss3.5); see protocol.ts's `FrameMessage.metrics`. `null`
   * before the first frame arrives. */
  metrics: Metrics | null;
  /** Whether video/recorder.ts's `CanvasRecorder` is currently recording -- drives the toolbar's
   * Record button icon (ui/toolbar.ts) and the HUD's recording indicator (ui/hud.ts). Recording is
   * independent of `running`/simulation reset: it just keeps capturing whatever gets rendered. */
  isRecording: boolean;
  /** `simTime` at the moment recording started, or `null` while not recording. Used (instead of
   * wall-clock time) to show elapsed recording duration in the HUD, consistent with this app's
   * sim-time-based accounting elsewhere. */
  recordingStartSimTime: number | null;
}

export function createInitialState(): AppState {
  return {
    params: null,
    running: true,
    drumAngle: 0,
    simTime: 0,
    achievedTimeScale: 1,
    subStepsPerSecondAchieved: 0,
    subStepsPerSecondRequired: 0,
    ballPositions: new Float32Array(0),
    ballOrientations: new Float32Array(0),
    ballRadiusM: 0,
    fluidPositions: new Float32Array(0),
    fluidDye: new Float32Array(0),
    fluidSurface: new Float32Array(0),
    metrics: null,
    isRecording: false,
    recordingStartSimTime: null,
  };
}
