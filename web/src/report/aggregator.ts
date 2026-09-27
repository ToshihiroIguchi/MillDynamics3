// Online (Welford's algorithm) mean/variance/min/max aggregator feeding the PDF report
// (report/pdf.ts, ui/metricsPanel.ts). Unlike metrics/history.ts's `MetricsHistory` -- a bounded
// ring buffer sized for a 60s sparkline/CSV window -- this tracks running statistics per metric id
// in O(1) memory regardless of how long the aggregation window (a whole simulation run) lasts, so
// it never needs to store individual samples.

/** Aggregated statistics for one metric id over the current window. */
export interface MetricStats {
  /** Number of finite samples this metric received (may differ from `ReportAggregator.samples`
   * when a metric was `null`/unavailable on some frames, e.g. collision_rate while coarse-grained). */
  count: number;
  mean: number;
  /** Sample standard deviation (Bessel-corrected, n-1 denominator); `0` for a single sample. */
  std: number;
  min: number;
  max: number;
}

interface Accumulator {
  count: number;
  mean: number;
  m2: number;
  min: number;
  max: number;
}

/**
 * Feed one `add(simTime, values)` call per throttled metrics update; read back per-id stats with
 * `statsFor`, or `null` if that id never received a finite value (including the whole-aggregator
 * "no data yet" case -- `add` was never called, or every value fed for that id was null/NaN). Never
 * throws on zero samples: every accessor has a documented "no data" answer instead.
 */
export class ReportAggregator {
  private accumulators = new Map<string, Accumulator>();
  private startSimTime: number | null = null;
  private endSimTime: number | null = null;
  private sampleCount = 0;

  /** Records one set of metric values (id -> value) at `simTime`. `null`/`undefined`/non-finite
   * values are skipped per-id (that id's stats simply don't advance this call) rather than
   * poisoning its running mean/variance with a placeholder. */
  add(simTime: number, values: Record<string, number | null | undefined>): void {
    if (this.startSimTime === null) this.startSimTime = simTime;
    this.endSimTime = simTime;
    this.sampleCount += 1;

    for (const [id, value] of Object.entries(values)) {
      if (value === null || value === undefined || !Number.isFinite(value)) continue;
      let acc = this.accumulators.get(id);
      if (!acc) {
        acc = { count: 0, mean: 0, m2: 0, min: value, max: value };
        this.accumulators.set(id, acc);
      }
      acc.count += 1;
      const delta = value - acc.mean;
      acc.mean += delta / acc.count;
      const delta2 = value - acc.mean;
      acc.m2 += delta * delta2;
      acc.min = Math.min(acc.min, value);
      acc.max = Math.max(acc.max, value);
    }
  }

  /** Clears all accumulated statistics and the aggregation window -- called on simulation reset
   * (ui/metricsPanel.ts's `reset()`) to start a fresh window rather than mixing pre-/post-reset
   * samples. */
  reset(): void {
    this.accumulators.clear();
    this.startSimTime = null;
    this.endSimTime = null;
    this.sampleCount = 0;
  }

  /** Stats for one metric id, or `null` if it never received a finite value. */
  statsFor(id: string): MetricStats | null {
    const acc = this.accumulators.get(id);
    if (!acc || acc.count === 0) return null;
    const variance = acc.count > 1 ? acc.m2 / (acc.count - 1) : 0;
    return { count: acc.count, mean: acc.mean, std: Math.sqrt(variance), min: acc.min, max: acc.max };
  }

  /** Every metric id that has received at least one finite value so far. */
  ids(): string[] {
    return [...this.accumulators.keys()];
  }

  /** Sim time of the first `add()` call since construction/reset, or `null` before any call. */
  get startTime(): number | null {
    return this.startSimTime;
  }

  /** Sim time of the most recent `add()` call, or `null` before any call. */
  get endTime(): number | null {
    return this.endSimTime;
  }

  /** Total number of `add()` calls since construction/reset (not per-metric; see `MetricStats.count`
   * for how many of those calls actually carried a finite value for a given id). */
  get samples(): number {
    return this.sampleCount;
  }
}
