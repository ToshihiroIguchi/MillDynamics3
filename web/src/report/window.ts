// Pure conversion of the report-window inputs (revolutions) to the aggregator's seconds.

export interface ReportWindowSeconds {
  /** False when the drum speed is unknown or not positive: the window is then unset. */
  valid: boolean;
  warmupS: number;
  /** Absolute sim time at which the window ends; 0 = no end. */
  endS: number;
}

/**
 * Converts warm-up and measurement lengths in drum revolutions to absolute sim seconds at `rpm`.
 * The window is [warmup, warmup + measure]. Non-finite / non-positive `measureRev` means no end
 * (`endS` 0); negative or non-finite `warmupRev` counts as 0. With `rpm` null/NaN/<= 0 the window
 * is unset (warm-up 0, no end, `valid` false).
 */
export function windowSecondsFromRevs(warmupRev: number, measureRev: number, rpm: number | null): ReportWindowSeconds {
  if (rpm === null || !Number.isFinite(rpm) || rpm <= 0) return { valid: false, warmupS: 0, endS: 0 };
  const warm = Number.isFinite(warmupRev) && warmupRev > 0 ? warmupRev : 0;
  const meas = Number.isFinite(measureRev) && measureRev > 0 ? measureRev : 0;
  const secPerRev = 60 / rpm;
  return { valid: true, warmupS: warm * secPerRev, endS: meas > 0 ? (warm + meas) * secPerRev : 0 };
}
