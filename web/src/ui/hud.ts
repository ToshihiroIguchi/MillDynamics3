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

export function createHud(): Hud {
  const el = document.createElement("div");
  el.className = "hud";

  return {
    el,
    update(state, fps) {
      const mill = millOf(state.params);
      const ballCount = state.ballPositions.length / 2;
      const lines = [`t = ${state.simTime.toFixed(1)} s`, `balls: ${ballCount}`, `fps: ${fps.toFixed(0)}`, `speed: ${state.achievedTimeScale.toFixed(2)}x`];
      if (mill) {
        lines.splice(1, 0, `${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc, Nc=${criticalSpeedRpm(mill.diameter_m).toFixed(1)})`);
      }
      el.textContent = lines.join("  |  ");
    },
  };
}
