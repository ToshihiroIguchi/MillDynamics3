// Smoothing for the speed badge (ui/hud.ts). The raw per-frame achieved rate is noisy (a frame's step
// count is an integer number of fixed sub-steps, see worker.ts), so the badge shows an exponential
// moving average, refreshed at a low rate, with hysteresis on the REAL TIME / SLOW MOTION switch.

/** Smoothed values the badge shows; the three rates come from the same average so they stay consistent. */
export interface SpeedReading {
  timeScale: number;
  subStepsAchieved: number;
  subStepsRequired: number;
  slow: boolean;
}

export interface SpeedSmootherOptions {
  /** Time constant of the moving average, in ms. */
  tauMs?: number;
  /** Minimum time between published readings, in ms. */
  refreshMs?: number;
  /** Go to SLOW MOTION when the smoothed time scale falls below this fraction of 1.0x. */
  slowBelow?: number;
  /** Return to REAL TIME only once the smoothed time scale reaches this (> slowBelow: hysteresis). */
  realAtOrAbove?: number;
}

export class SpeedSmoother {
  private readonly tauMs: number;
  private readonly refreshMs: number;
  private readonly slowBelow: number;
  private readonly realAtOrAbove: number;
  private scale = 0;
  private achieved = 0;
  private required = 0;
  private lastSampleMs: number | null = null;
  private lastPublishMs = Number.NEGATIVE_INFINITY;
  private slow = false;

  constructor(options: SpeedSmootherOptions = {}) {
    this.tauMs = options.tauMs ?? 1000;
    this.refreshMs = options.refreshMs ?? 500;
    this.slowBelow = options.slowBelow ?? 0.98;
    this.realAtOrAbove = options.realAtOrAbove ?? 0.99;
  }

  /**
   * Feeds one frame's raw values. Returns a reading when one is due (the first sample, then at most
   * every `refreshMs`), otherwise `null` and the badge keeps showing the previous reading.
   */
  push(nowMs: number, timeScale: number, subStepsAchieved: number, subStepsRequired: number): SpeedReading | null {
    if (this.lastSampleMs === null) {
      this.scale = timeScale;
      this.achieved = subStepsAchieved;
      this.required = subStepsRequired;
    } else {
      const alpha = 1 - Math.exp(-Math.max(0, nowMs - this.lastSampleMs) / this.tauMs);
      this.scale += alpha * (timeScale - this.scale);
      this.achieved += alpha * (subStepsAchieved - this.achieved);
      this.required += alpha * (subStepsRequired - this.required);
    }
    this.lastSampleMs = nowMs;

    if (nowMs - this.lastPublishMs < this.refreshMs) return null;
    this.lastPublishMs = nowMs;
    if (this.slow) {
      if (this.scale >= this.realAtOrAbove) this.slow = false;
    } else if (this.scale < this.slowBelow) {
      this.slow = true;
    }
    return {
      timeScale: this.scale,
      subStepsAchieved: this.achieved,
      subStepsRequired: this.required,
      slow: this.slow,
    };
  }
}
