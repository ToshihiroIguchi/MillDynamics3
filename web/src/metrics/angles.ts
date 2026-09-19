// Angle-convention helpers shared by the HUD and the metrics panel. See
// crates/mill-core/src/metrics.rs's module doc comment for the internal convention every
// `*_angle_rad` field in metrics/types.ts's `Metrics` uses: `atan2(y, x)` in `[0, 2*pi)`,
// measured counter-clockwise from `+x` in the world frame.

/**
 * Converts an internal `atan2`-from-`+x` angle (radians, `[0, 2*pi)`, CCW convention) to the
 * mill-literature "angle from vertical" convention (0 deg at 12 o'clock, increasing clockwise for
 * a counter-clockwise-rotating drum, flipped for a clockwise-rotating one). Mirrors
 * `mill_core::metrics::to_vertical_degrees` (crates/mill-core/src/metrics.rs) exactly -- see its
 * doc comment for the derivation.
 */
export function toVerticalDegrees(atan2Rad: number, omega: number): number {
  const sign = omega >= 0 ? -1 : 1;
  const verticalRad = sign * (atan2Rad - Math.PI / 2);
  const deg = (verticalRad * 180) / Math.PI;
  return ((deg % 360) + 360) % 360;
}

/** Plain radians-to-degrees, wrapped into `[0, 360)` -- no "from vertical" conversion (used for
 * the pool extent, which is reported in the raw internal angle convention). */
export function toDegrees(atan2Rad: number): number {
  const deg = (atan2Rad * 180) / Math.PI;
  return ((deg % 360) + 360) % 360;
}
