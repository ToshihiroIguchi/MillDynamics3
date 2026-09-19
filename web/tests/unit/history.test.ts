import { describe, expect, it } from "vitest";
import { MetricsHistory } from "../../src/metrics/history";

describe("MetricsHistory", () => {
  it("throttles pushes by simulation time", () => {
    const history = new MetricsHistory(100, 0.1);
    history.push(0, { a: 1 });
    history.push(0.05, { a: 2 }); // too soon after the first sample, dropped
    expect(history.length).toBe(1);
    history.push(0.1, { a: 3 }); // interval elapsed, lands
    expect(history.length).toBe(2);
    expect(history.seriesFor("a")).toEqual([1, 3]);
  });

  it("evicts the oldest sample once capacity is exceeded", () => {
    const history = new MetricsHistory(3, 0);
    history.push(0, { x: 10 });
    history.push(1, { x: 11 });
    history.push(2, { x: 12 });
    history.push(3, { x: 13 });
    expect(history.length).toBe(3);
    expect(history.seriesFor("x")).toEqual([11, 12, 13]);
  });

  it("clears all samples and the throttle high-water mark on reset", () => {
    const history = new MetricsHistory(100, 0.1);
    history.push(0, { a: 1 });
    history.push(1, { a: 2 });
    history.reset();
    expect(history.length).toBe(0);
    history.push(0, { a: 3 });
    expect(history.length).toBe(1);
  });

  it("clears history when sim time moves backward (unannounced reset)", () => {
    const history = new MetricsHistory(100, 0.1);
    history.push(0, { a: 1 });
    history.push(1, { a: 2 });
    history.push(0.5, { a: 3 }); // simTime < lastSampledSimTime -> auto-clear
    expect(history.length).toBe(1);
    expect(history.seriesFor("a")).toEqual([3]);
  });

  it("renders CSV with labelled headers, formatted sim time, and empty cells for null values", () => {
    const history = new MetricsHistory(100, 0);
    history.push(1.23456, { a: 1, b: null });
    const csv = history.toCsv([
      { id: "a", label: "A", unit: "m" },
      { id: "b", label: "B" },
    ]);
    const lines = csv.split("\n");
    expect(lines[0]).toBe("sim_time_s,A (m),B");
    expect(lines).toHaveLength(2);
    const row = lines[1];
    expect(row).toBeDefined();
    const cells = (row ?? "").split(",");
    expect(cells[0]).toBe("1.235");
    expect(cells[1]).toBe("1");
    expect(cells[2]).toBe("");
  });
});
