// On-screen HUD: sim time, rpm/%Nc, fps, achieved time scale. Everything else (toe/shoulder,
// pool extent, mixing index, grinding/solver diagnostics, ...) lives in the metrics panel
// (ui/metricsPanel.ts, docs/PLAN.md ss4.4) now that it exists.

import type { ParamsJson } from "../protocol";
import { criticalSpeedRpm, percentCriticalOf, rpmOf } from "../params/derived";
import type { AppState } from "../state";

export interface Hud {
  el: HTMLElement;
  update(state: AppState, fps: number): void;
}

function millOf(params: ParamsJson | null): { diameter_m: number; speed_mode: string; speed_value: number } | null {
  const mill = params?.mill as { diameter_m?: number; speed_mode?: string; speed_value?: number } | undefined;
  if (!mill || mill.diameter_m === undefined || mill.speed_mode === undefined || mill.speed_value === undefined) {
    return null;
  }
  return { diameter_m: mill.diameter_m, speed_mode: mill.speed_mode, speed_value: mill.speed_value };
}

/**
 * Threshold, as a fraction of `1.0x`, below which the achieved rate reads as "SLOW MOTION" rather
 * than "REAL TIME". Not exactly `1.0`: `achievedTimeScale` is itself noisy frame-to-frame (a
 * single frame's `stepsRun` is an integer count of a possibly-large `fixedSubDt()`, see
 * worker.ts), so requiring *exactly* `>= 1.0` would flicker between states on an otherwise-healthy
 * run. `0.98` absorbs that noise without hiding a real, sustained shortfall.
 */
const REAL_TIME_THRESHOLD = 0.98;

export function createHud(): Hud {
  const el = document.createElement("div");
  el.className = "hud";

  return {
    el,
    update(state, fps) {
      const mill = millOf(state.params);
      const ballCount = state.ballPositions.length / 2;
      const isRealTime = state.achievedTimeScale >= REAL_TIME_THRESHOLD;
      // Deliberately not a subtle number: a simulator whose physics runs slower than the clock
      // it's being watched against must say so plainly, not bury it in a "speed: 0.12x" figure
      // easy to read as a setting rather than a shortfall (docs/PERF.md).
      const speedLine = isRealTime
        ? `REAL TIME (${state.achievedTimeScale.toFixed(2)}x)`
        : `SLOW MOTION ${state.achievedTimeScale.toFixed(2)}x -- ${state.subStepsPerSecondAchieved.toFixed(0)}/${state.subStepsPerSecondRequired.toFixed(0)} sub-steps/s`;
      const lines = [`t = ${state.simTime.toFixed(1)} s`, `balls: ${ballCount}`, `fps: ${fps.toFixed(0)}`, speedLine];
      if (mill) {
        lines.splice(1, 0, `${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc, Nc=${criticalSpeedRpm(mill.diameter_m).toFixed(1)})`);
      }
      el.textContent = lines.join("  |  ");
      el.classList.toggle("hud-slow-motion", !isRealTime);
    },
  };
}
