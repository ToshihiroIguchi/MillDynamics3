// On-screen HUD (v1: sim time, rpm/%Nc, ball count, fps, achieved time scale). Extended with
// slurry/mixing metrics in M3/M5 (see docs/PLAN.md ss4.4).

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

/** Shape of the JSON-encoded `mill_core::metrics::Metrics` (crates/mill-core/src/metrics.rs), as
 * relevant to the HUD -- only the fields this file actually reads. */
interface MetricsJson {
  toe_angle_rad: number | null;
  shoulder_angle_rad: number | null;
  pool_angle_min_rad: number | null;
  pool_angle_max_rad: number | null;
  mixing_index: number | null;
  drum_omega_rad_s: number;
}

/**
 * Converts an internal `atan2`-from-`+x` angle (radians, `[0, 2*pi)`, CCW convention) to the
 * mill-literature "angle from vertical" convention (0 deg at 12 o'clock, increasing clockwise for
 * a counter-clockwise-rotating drum, flipped for a clockwise-rotating one). Mirrors
 * `mill_core::metrics::to_vertical_degrees` (crates/mill-core/src/metrics.rs) exactly -- see its
 * doc comment for the derivation.
 */
function toVerticalDegrees(atan2Rad: number, omega: number): number {
  const sign = omega >= 0 ? -1 : 1;
  const verticalRad = sign * (atan2Rad - Math.PI / 2);
  const deg = (verticalRad * 180) / Math.PI;
  return ((deg % 360) + 360) % 360;
}

/** Plain radians-to-degrees, wrapped into `[0, 360)` -- no "from vertical" conversion (used for
 * the pool extent, which the HUD reports in the raw internal angle convention). */
function toDegrees(atan2Rad: number): number {
  const deg = (atan2Rad * 180) / Math.PI;
  return ((deg % 360) + 360) % 360;
}

/** Formats one HUD text segment from the worker's opaque `metricsJson` string (protocol.ts's
 * `FrameMessage.metricsJson`), or `null` if there is nothing to show yet (e.g. before the first
 * frame). Missing/unavailable individual metrics (e.g. a centrifuged charge has no toe/shoulder)
 * render as "-". */
function formatMetricsSegment(metricsJson: string): string | null {
  if (!metricsJson) return null;
  let m: MetricsJson;
  try {
    m = JSON.parse(metricsJson) as MetricsJson;
  } catch {
    return null;
  }
  const toe = m.toe_angle_rad === null ? "-" : `${toVerticalDegrees(m.toe_angle_rad, m.drum_omega_rad_s).toFixed(1)}°`;
  const shoulder =
    m.shoulder_angle_rad === null ? "-" : `${toVerticalDegrees(m.shoulder_angle_rad, m.drum_omega_rad_s).toFixed(1)}°`;
  const pool =
    m.pool_angle_min_rad === null || m.pool_angle_max_rad === null
      ? "-"
      : `${toDegrees(m.pool_angle_min_rad).toFixed(0)}-${toDegrees(m.pool_angle_max_rad).toFixed(0)}°`;
  const mixing = m.mixing_index === null ? "-" : m.mixing_index.toFixed(2);
  return `toe/shoulder: ${toe}/${shoulder}  pool: ${pool}  mix: ${mixing}`;
}

export function createHud(): Hud {
  const el = document.createElement("div");
  el.className = "hud";
  document.body.appendChild(el);

  return {
    el,
    update(state, fps) {
      const mill = millOf(state.params);
      const ballCount = state.ballPositions.length / 2;
      const lines = [`t = ${state.simTime.toFixed(1)} s`, `balls: ${ballCount}`, `fps: ${fps.toFixed(0)}`, `speed: ${state.achievedTimeScale.toFixed(2)}x`];
      if (mill) {
        lines.splice(1, 0, `${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc, Nc=${criticalSpeedRpm(mill.diameter_m).toFixed(1)})`);
      }
      const metricsSegment = formatMetricsSegment(state.metricsJson);
      if (metricsSegment) {
        lines.push(metricsSegment);
      }
      el.textContent = lines.join("  |  ");
    },
  };
}
