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

/** One decimated time-series sample. */
export interface SeriesPoint {
  t: number;
  v: number;
}

/** Upper bound on stored points per metric id (see `ReportAggregator.seriesFor`). */
export const MAX_SERIES_POINTS = 1000;

/** Sample window `[warmupS, endS]` in sim seconds. `endS <= 0` (or non-finite) means "no end". */
export interface ReportWindow {
  warmupS: number;
  endS: number;
}

interface Series {
  points: SeriesPoint[];
  /** Keep every `stride`-th accepted sample; doubles each time the store fills up. */
  stride: number;
  seen: number;
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
  private window: ReportWindow = { warmupS: 0, endS: 0 };
  private series = new Map<string, Series>();
  private histEdges: number[] | null = null;
  private histSums: number[] | null = null;
  private histLastCounts: number[] | null = null;
  private histDuration = 0;
  private lastAcceptedTime: number | null = null;

  /** Sets the sample window. Samples outside `[warmupS, endS]` are ignored by `add`. Does not
   * clear already-collected data. */
  setWindow(warmupS: number, endS: number): void {
    this.window = {
      warmupS: Number.isFinite(warmupS) && warmupS > 0 ? warmupS : 0,
      endS: Number.isFinite(endS) && endS > 0 ? endS : 0,
    };
  }

  getWindow(): ReportWindow {
    return { ...this.window };
  }

  /** True once `simTime` is past a configured window end. */
  isPastEnd(simTime: number): boolean {
    return this.window.endS > 0 && simTime > this.window.endS;
  }

  /** Records one set of metric values (id -> value) at `simTime`. `null`/`undefined`/non-finite
   * values are skipped per-id (that id's stats simply don't advance this call) rather than
   * poisoning its running mean/variance with a placeholder. */
  add(
    simTime: number,
    values: Record<string, number | null | undefined>,
    histogram?: { bin_edges_j: number[]; counts_per_s: number[] } | null,
  ): void {
    if (simTime < this.window.warmupS) return;
    if (this.window.endS > 0 && simTime > this.window.endS) return;
    this.addHistogram(simTime, histogram);
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
      this.pushSeries(id, simTime, value);
    }
  }

  private pushSeries(id: string, t: number, v: number): void {
    let s = this.series.get(id);
    if (!s) {
      s = { points: [], stride: 1, seen: 0 };
      this.series.set(id, s);
    }
    const keep = s.seen % s.stride === 0;
    s.seen += 1;
    if (!keep) return;
    s.points.push({ t, v });
    if (s.points.length >= MAX_SERIES_POINTS) {
      s.points = s.points.filter((_, i) => i % 2 === 0);
      s.stride *= 2;
    }
  }

  private addHistogram(simTime: number, h: { bin_edges_j: number[]; counts_per_s: number[] } | null | undefined): void {
    const dt = this.lastAcceptedTime === null ? 0 : simTime - this.lastAcceptedTime;
    this.lastAcceptedTime = simTime;
    if (!h || h.counts_per_s.length === 0) return;
    if (!this.histSums || this.histSums.length !== h.counts_per_s.length) {
      this.histSums = new Array<number>(h.counts_per_s.length).fill(0);
      this.histDuration = 0;
    }
    this.histEdges = h.bin_edges_j.slice();
    this.histLastCounts = h.counts_per_s.slice();
    if (dt > 0) {
      for (let i = 0; i < h.counts_per_s.length; i++) this.histSums[i]! += h.counts_per_s[i]! * dt;
      this.histDuration += dt;
    }
  }

  /** Decimated (<= `MAX_SERIES_POINTS`) time series for one metric id; empty if none. */
  seriesFor(id: string): readonly SeriesPoint[] {
    return this.series.get(id)?.points ?? [];
  }

  /** Window-mean impact histogram (sum of counts_per_s x dt / total dt per bin), or `null` if no
   * histogram was ever fed. With a zero-length window (single sample) falls back to that sample. */
  windowHistogram(): { bin_edges_j: number[]; counts_per_s: number[] } | null {
    if (!this.histSums || !this.histEdges) return null;
    const counts =
      this.histDuration > 0
        ? this.histSums.map((x) => x / this.histDuration)
        : (this.histLastCounts ?? this.histSums).slice();
    return { bin_edges_j: this.histEdges.slice(), counts_per_s: counts };
  }

  /** Clears all accumulated statistics and the aggregation window -- called on simulation reset
   * (ui/metricsPanel.ts's `reset()`) to start a fresh window rather than mixing pre-/post-reset
   * samples. */
  reset(): void {
    this.accumulators.clear();
    this.startSimTime = null;
    this.endSimTime = null;
    this.sampleCount = 0;
    this.series.clear();
    this.histEdges = null;
    this.histSums = null;
    this.histLastCounts = null;
    this.histDuration = 0;
    this.lastAcceptedTime = null;
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
