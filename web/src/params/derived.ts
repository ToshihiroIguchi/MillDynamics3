// Small, display-only mirrors of the mill-core formulas in crates/mill-core/src/params.rs
// (MillParams::critical_speed_rpm/rpm/percent_critical). Duplicated here deliberately: these are
// tiny, stable, well-documented formulas (docs/PLAN.md ss3.1), so mirroring them in TypeScript for
// instant UI feedback is lower-risk than round-tripping through the worker/wasm for every keypress
// in the parameters modal. mill-core remains the source of truth for anything that affects the
// actual simulation (this file is read-only display).

export interface MillLike {
  diameter_m: number;
  speed_mode: string; // "rpm" | "percent_critical"
  speed_value: number;
}

/** Critical speed (rpm): Nc = 42.3 / sqrt(D), D in meters. */
export function criticalSpeedRpm(diameterM: number): number {
  return 42.3 / Math.sqrt(diameterM);
}

export function rpmOf(mill: MillLike): number {
  return mill.speed_mode === "rpm" ? mill.speed_value : (criticalSpeedRpm(mill.diameter_m) * mill.speed_value) / 100;
}

export function percentCriticalOf(mill: MillLike): number {
  return mill.speed_mode === "percent_critical" ? mill.speed_value : (100 * mill.speed_value) / criticalSpeedRpm(mill.diameter_m);
}
