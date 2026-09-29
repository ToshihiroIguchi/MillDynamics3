import { describe, expect, it } from "vitest";
import { MAX_SERIES_POINTS, ReportAggregator } from "../../src/report/aggregator";

describe("ReportAggregator", () => {
  it("reports no data (null stats, null window, 0 samples) before any add() call", () => {
    const agg = new ReportAggregator();
    expect(agg.statsFor("a")).toBeNull();
    expect(agg.startTime).toBeNull();
    expect(agg.endTime).toBeNull();
    expect(agg.samples).toBe(0);
    expect(agg.ids()).toEqual([]);
  });

  it("matches hand-computed mean/std/min/max for a small synthetic series", () => {
    // Series: 1, 2, 3, 4, 5. mean = 3, sample variance (n-1) = 2.5, std = sqrt(2.5).
    const agg = new ReportAggregator();
    const series = [1, 2, 3, 4, 5];
    series.forEach((v, i) => agg.add(i * 0.1, { a: v }));

    const stats = agg.statsFor("a");
    expect(stats).not.toBeNull();
    expect(stats!.count).toBe(5);
    expect(stats!.mean).toBeCloseTo(3, 10);
    expect(stats!.std).toBeCloseTo(Math.sqrt(2.5), 10);
    expect(stats!.min).toBe(1);
    expect(stats!.max).toBe(5);
  });

  it("reports std = 0 for a single sample (no Bessel-correction division by zero)", () => {
    const agg = new ReportAggregator();
    agg.add(0, { a: 42 });
    const stats = agg.statsFor("a");
    expect(stats).not.toBeNull();
    expect(stats!.count).toBe(1);
    expect(stats!.mean).toBe(42);
    expect(stats!.std).toBe(0);
    expect(stats!.min).toBe(42);
    expect(stats!.max).toBe(42);
  });

  it("tracks the aggregation window (start/end sim time) and total sample count across multiple ids", () => {
    const agg = new ReportAggregator();
    agg.add(1.0, { a: 1, b: 10 });
    agg.add(2.5, { a: 2, b: 20 });
    agg.add(4.0, { a: 3, b: 30 });

    expect(agg.startTime).toBe(1.0);
    expect(agg.endTime).toBe(4.0);
    expect(agg.samples).toBe(3);
    expect(agg.ids().sort()).toEqual(["a", "b"]);
    expect(agg.statsFor("b")!.mean).toBeCloseTo(20, 10);
  });

  it("skips null/undefined/non-finite values per-id without affecting the window or other ids", () => {
    const agg = new ReportAggregator();
    agg.add(0, { a: 1, b: null });
    agg.add(1, { a: null, b: 5 });
    agg.add(2, { a: 3, b: undefined });
    agg.add(3, { a: NaN, b: 7 });

    // 4 add() calls -> 4 samples in the window, even though each id only got some of them.
    expect(agg.samples).toBe(4);
    expect(agg.startTime).toBe(0);
    expect(agg.endTime).toBe(3);

    const a = agg.statsFor("a");
    expect(a!.count).toBe(2); // 1, 3
    expect(a!.mean).toBeCloseTo(2, 10);

    const b = agg.statsFor("b");
    expect(b!.count).toBe(2); // 5, 7
    expect(b!.mean).toBeCloseTo(6, 10);
  });

  it("returns null for an id that only ever received non-finite values", () => {
    const agg = new ReportAggregator();
    agg.add(0, { a: NaN });
    agg.add(1, { a: null });
    expect(agg.statsFor("a")).toBeNull();
    expect(agg.samples).toBe(2);
  });

  it("clears accumulators and the window on reset(), starting a fresh window on the next add()", () => {
    const agg = new ReportAggregator();
    agg.add(0, { a: 1 });
    agg.add(1, { a: 2 });
    agg.reset();

    expect(agg.statsFor("a")).toBeNull();
    expect(agg.startTime).toBeNull();
    expect(agg.endTime).toBeNull();
    expect(agg.samples).toBe(0);
    expect(agg.ids()).toEqual([]);

    agg.add(5, { a: 100 });
    expect(agg.startTime).toBe(5);
    expect(agg.endTime).toBe(5);
    expect(agg.samples).toBe(1);
    expect(agg.statsFor("a")!.mean).toBe(100);
  });
});

describe("ReportAggregator windowing", () => {
  it("ignores samples before warm-up and after the end", () => {
    const agg = new ReportAggregator();
    agg.setWindow(2, 7);
    for (let t = 0; t <= 10; t += 1) agg.add(t, { a: t });
    expect(agg.startTime).toBe(2);
    expect(agg.endTime).toBe(7);
    expect(agg.samples).toBe(6); // t = 2..7
    expect(agg.statsFor("a")!.mean).toBeCloseTo(4.5, 10);
    expect(agg.isPastEnd(7)).toBe(false);
    expect(agg.isPastEnd(7.1)).toBe(true);
  });

  it("end = 0 means no end, and reset() keeps the window but clears data", () => {
    const agg = new ReportAggregator();
    agg.setWindow(1, 0);
    agg.add(0.5, { a: 1 });
    agg.add(100, { a: 2 });
    expect(agg.samples).toBe(1);
    agg.reset();
    expect(agg.samples).toBe(0);
    expect(agg.seriesFor("a")).toEqual([]);
    expect(agg.getWindow()).toEqual({ warmupS: 1, endS: 0 });
  });
});

describe("ReportAggregator series decimation", () => {
  it("never stores more than MAX_SERIES_POINTS and keeps the full time span", () => {
    const agg = new ReportAggregator();
    const n = 20000;
    for (let i = 0; i < n; i++) agg.add(i * 0.01, { a: i });
    const pts = agg.seriesFor("a");
    expect(pts.length).toBeLessThanOrEqual(MAX_SERIES_POINTS);
    expect(pts.length).toBeGreaterThan(MAX_SERIES_POINTS / 2 - 1);
    expect(pts[0]!.t).toBe(0);
    expect(pts[pts.length - 1]!.t).toBeGreaterThan(n * 0.01 * 0.95);
    // Exact stats are unaffected by decimation.
    expect(agg.statsFor("a")!.count).toBe(n);
  });
});

describe("ReportAggregator window-mean histogram", () => {
  it("weights each frame's counts_per_s by its dt", () => {
    const agg = new ReportAggregator();
    const edges = [1, 10, 100];
    agg.add(0, { a: 1 }, { bin_edges_j: edges, counts_per_s: [100, 100] }); // dt 0: no weight
    agg.add(1, { a: 1 }, { bin_edges_j: edges, counts_per_s: [10, 0] }); // dt 1
    agg.add(4, { a: 1 }, { bin_edges_j: edges, counts_per_s: [20, 40] }); // dt 3
    const h = agg.windowHistogram()!;
    expect(h.bin_edges_j).toEqual(edges);
    expect(h.counts_per_s[0]).toBeCloseTo((10 * 1 + 20 * 3) / 4, 10);
    expect(h.counts_per_s[1]).toBeCloseTo((0 * 1 + 40 * 3) / 4, 10);
  });

  it("is null with no histogram, and falls back to the sample for one frame", () => {
    const agg = new ReportAggregator();
    expect(agg.windowHistogram()).toBeNull();
    agg.add(0, {}, { bin_edges_j: [1, 2], counts_per_s: [5] });
    expect(agg.windowHistogram()!.counts_per_s).toEqual([5]);
  });
});
