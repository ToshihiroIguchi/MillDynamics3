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

/** Subset of `MediaParams` (crates/mill-core/src/params.rs) that the formulas below need. */
export interface MediaLike {
  ball_diameter_m: number;
  fill_fraction: number;
  packing_fraction_2d: number;
  density_kg_m3: number;
}

/** Mirrors `Params::true_ball_count` (params.rs:392-400). `diameterM` is the mill's drum diameter. */
export function trueBallCount(diameterM: number, media: MediaLike): number {
  const rTrue = media.ball_diameter_m * 0.5;
  if (rTrue <= 0) return 0;
  const radiusM = diameterM * 0.5;
  const drumArea = Math.PI * radiusM * radiusM;
  return (media.fill_fraction * media.packing_fraction_2d * drumArea) / (Math.PI * rTrue * rTrue);
}

/** Mirrors `Params::effective_media` (params.rs:416-441). `diameterM` is the mill's drum diameter. */
export function effectiveMedia(
  diameterM: number,
  media: MediaLike,
  maxBalls: number,
): { trueDiameterM: number; diameterM: number; densityKgM3: number; ballCount: number; scaleFactor: number } {
  const trueDiameterM = media.ball_diameter_m;
  const densityKgM3 = media.density_kg_m3;
  const nTrue = trueBallCount(diameterM, media);
  if (nTrue > maxBalls && maxBalls > 0) {
    const scaleFactor = Math.sqrt(nTrue / maxBalls);
    const ballCount = Math.max(1, Math.round(nTrue / (scaleFactor * scaleFactor)));
    return { trueDiameterM, diameterM: trueDiameterM * scaleFactor, densityKgM3, ballCount, scaleFactor };
  }
  return { trueDiameterM, diameterM: trueDiameterM, densityKgM3, ballCount: Math.max(0, Math.round(nTrue)), scaleFactor: 1 };
}

/**
 * Mirrors the target particle count computed by `FluidParticles::seed_lattice` (pbf.rs:388-396).
 * Note `slurryFillFraction` is `slurry.fill_fraction`, a different field from `media.fill_fraction`.
 */
export function fluidParticleCountEstimate(radiusM: number, resolution: number, slurryFillFraction: number): number {
  const dx = radiusM / Math.max(1, resolution);
  const targetArea = slurryFillFraction * Math.PI * radiusM * radiusM;
  return Math.round(Math.max(0, targetArea / (dx * dx)));
}

/**
 * Conservative static estimate of dem.rs's `max_substep_displacement_over_diameter` metric
 * (dem.rs:352-356), using the drum's nominal rim speed (rpm -> rad/s -> tangential speed at the
 * wall) as an upper bound on ball speed, since this is a params-only estimate with no live
 * per-ball velocities available. This is NOT what dem.rs actually measures (dem.rs uses each
 * ball's own velocity, dominated by rotation and cascading, not just rim speed); the live
 * `Metrics.max_substep_displacement_over_diameter` (surfaced in the metrics panel) is the
 * authoritative value once the sim is running.
 */
export function substepDisplacementOverDiameter(
  mill: MillLike,
  timeScale: number,
  substeps: number,
  effectiveDiameterM: number,
): number {
  const omega = (rpmOf(mill) * 2 * Math.PI) / 60;
  const rimSpeedMS = omega * (mill.diameter_m / 2);
  const subDt = timeScale / 60 / Math.max(1, substeps);
  return effectiveDiameterM > 0 ? (rimSpeedMS * subDt) / effectiveDiameterM : 0;
}
