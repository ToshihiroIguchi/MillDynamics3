// On-screen HUD: sim time, rpm/%Nc, fps, achieved time scale. Everything else (toe/shoulder,
// pool extent, mixing index, grinding/solver diagnostics, ...) lives in the metrics panel
// (ui/metricsPanel.ts, docs/PLAN.md ss4.4) now that it exists.

import { effectiveMedia, criticalSpeedRpm, percentCriticalOf, rpmOf } from "../params/derived";
import type { ParamsJson } from "../protocol";
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

/** Subset of `media`/`simulation.max_balls` this HUD needs to mirror `Params::effective_media`
 * (see params/derived.ts's `effectiveMedia`) -- same defensive-cast pattern as `millOf` above. */
function mediaSubstitutionOf(
  params: ParamsJson | null,
): { trueDiameterM: number; diameterM: number; ballCount: number; scaleFactor: number } | null {
  const mill = params?.mill as { diameter_m?: number } | undefined;
  const media = params?.media as
    | { ball_diameter_m?: number; fill_fraction?: number; packing_fraction_2d?: number; density_kg_m3?: number }
    | undefined;
  const simulation = params?.simulation as { max_balls?: number } | undefined;
  if (
    mill?.diameter_m === undefined ||
    media?.ball_diameter_m === undefined ||
    media.fill_fraction === undefined ||
    media.packing_fraction_2d === undefined ||
    media.density_kg_m3 === undefined ||
    simulation?.max_balls === undefined
  ) {
    return null;
  }
  const eff = effectiveMedia(
    mill.diameter_m,
    {
      ball_diameter_m: media.ball_diameter_m,
      fill_fraction: media.fill_fraction,
      packing_fraction_2d: media.packing_fraction_2d,
      density_kg_m3: media.density_kg_m3,
    },
    simulation.max_balls,
  );
  return { trueDiameterM: eff.trueDiameterM, diameterM: eff.diameterM, ballCount: eff.ballCount, scaleFactor: eff.scaleFactor };
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

  const statsLine = document.createElement("div");
  statsLine.className = "hud-line";
  el.appendChild(statsLine);

  // Speed indicator: a graphical dot + proportional bar alongside the existing text, so a
  // sustained shortfall reads at a glance rather than requiring the number to be parsed (user
  // feedback: a text-only "0.37x" is easy to skim past).
  const speedRow = document.createElement("div");
  speedRow.className = "hud-speed";
  const speedDot = document.createElement("span");
  speedDot.className = "hud-speed-dot";
  const speedBar = document.createElement("span");
  speedBar.className = "hud-speed-bar";
  const speedBarFill = document.createElement("span");
  speedBarFill.className = "hud-speed-bar-fill";
  speedBar.appendChild(speedBarFill);
  const speedText = document.createElement("span");
  speedText.className = "hud-speed-text";
  speedRow.append(speedDot, speedBar, speedText);
  el.appendChild(speedRow);

  // Media-substitution badge: only shown once coarse-graining is actually active (the params
  // panel's Derived section shows this too, but that panel can be collapsed -- this keeps the
  // fact that the *simulated* ball size differs from the entered one visible on the main screen
  // itself, per user feedback).
  const mediaBadge = document.createElement("div");
  mediaBadge.className = "hud-media";
  mediaBadge.hidden = true;
  el.appendChild(mediaBadge);

  return {
    el,
    update(state, fps) {
      const mill = millOf(state.params);
      const ballCount = state.ballPositions.length / 2;
      const isRealTime = state.achievedTimeScale >= REAL_TIME_THRESHOLD;

      const lines = [`t = ${state.simTime.toFixed(1)} s`, `balls: ${ballCount}`, `fps: ${fps.toFixed(0)}`];
      if (mill) {
        lines.splice(1, 0, `${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc, Nc=${criticalSpeedRpm(mill.diameter_m).toFixed(1)})`);
      }
      statsLine.textContent = lines.join("  |  ");

      // Deliberately not a subtle number: a simulator whose physics runs slower than the clock
      // it's being watched against must say so plainly, not bury it in a "speed: 0.12x" figure
      // easy to read as a setting rather than a shortfall (docs/PERF.md).
      speedText.textContent = isRealTime
        ? `REAL TIME (${state.achievedTimeScale.toFixed(2)}x)`
        : `SLOW MOTION ${state.achievedTimeScale.toFixed(2)}x -- ${state.subStepsPerSecondAchieved.toFixed(0)}/${state.subStepsPerSecondRequired.toFixed(0)} sub-steps/s`;
      const fillPct = Math.max(0, Math.min(1, state.achievedTimeScale)) * 100;
      speedBarFill.style.width = `${fillPct}%`;
      speedRow.classList.toggle("hud-speed-slow", !isRealTime);
      el.classList.toggle("hud-slow-motion", !isRealTime);

      const media = mediaSubstitutionOf(state.params);
      if (media && media.scaleFactor > 1) {
        mediaBadge.hidden = false;
        mediaBadge.textContent =
          `media: ${(media.trueDiameterM * 1000).toFixed(1)} mm -> ${(media.diameterM * 1000).toFixed(1)} mm ` +
          `(substituted, x${media.scaleFactor.toFixed(2)}, N=${media.ballCount})`;
      } else {
        mediaBadge.hidden = true;
      }
    },
  };
}
