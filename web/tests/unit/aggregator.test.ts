import { describe, expect, it } from "vitest";
import { ReportAggregator } from "../../src/report/aggregator";

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
