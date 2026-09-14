import type { ParamsJson } from "./protocol";

/** App-level state, updated from worker messages and read by the render loop. */
export interface AppState {
  params: ParamsJson | null;
  running: boolean;
  drumAngle: number;
  simTime: number;
  achievedTimeScale: number;
  ballPositions: Float32Array;
  ballOrientations: Float32Array;
  ballRadiusM: number;
  fluidPositions: Float32Array;
  fluidDye: Float32Array;
  fluidSurface: Float32Array;
}

export function createInitialState(): AppState {
  return {
    params: null,
    running: true,
    drumAngle: 0,
    simTime: 0,
    achievedTimeScale: 1,
    ballPositions: new Float32Array(0),
    ballOrientations: new Float32Array(0),
    ballRadiusM: 0,
    fluidPositions: new Float32Array(0),
    fluidDye: new Float32Array(0),
    fluidSurface: new Float32Array(0),
  };
}
