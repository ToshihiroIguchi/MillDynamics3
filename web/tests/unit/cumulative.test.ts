import { describe, expect, it } from "vitest";
import { energyCumulativeFraction } from "../../src/report/pdf";

describe("energyCumulativeFraction", () => {
  it("is monotone, ends at 1, and weights bins by rate x geometric-mean energy", () => {
    const edges = [1, 10, 100];
    const cum = energyCumulativeFraction(edges, [10, 1])!;
    const w0 = 10 * Math.sqrt(10);
    const w1 = 1 * Math.sqrt(1000);
    expect(cum[0]).toBeCloseTo(w0 / (w0 + w1), 10);
    expect(cum[1]).toBeCloseTo(1, 10);
  });

  it("returns null for empty or all-zero histograms", () => {
    expect(energyCumulativeFraction([], [])).toBeNull();
    expect(energyCumulativeFraction([1, 10], [0])).toBeNull();
  });
});
