import { describe, expect, it } from "vitest";
import { windowSecondsFromRevs } from "../../src/report/window";

describe("windowSecondsFromRevs", () => {
  it("converts revolutions to seconds at the given rpm", () => {
    const w = windowSecondsFromRevs(4, 8, 100);
    expect(w.valid).toBe(true);
    expect(w.warmupS).toBeCloseTo(2.4, 12);
    expect(w.endS).toBeCloseTo(7.2, 12);
  });
  it("treats zero / non-finite measure as no end", () => {
    expect(windowSecondsFromRevs(4, 0, 60).endS).toBe(0);
    expect(windowSecondsFromRevs(4, Number.NaN, 60).endS).toBe(0);
    expect(windowSecondsFromRevs(4, 0, 60).warmupS).toBeCloseTo(4, 12);
  });
  it("clamps negative warm-up to zero", () => {
    const w = windowSecondsFromRevs(-3, 2, 30);
    expect(w.warmupS).toBe(0);
    expect(w.endS).toBeCloseTo(4, 12);
  });
  it("is unset without a positive rpm", () => {
    for (const rpm of [null, 0, -5, Number.NaN]) {
      expect(windowSecondsFromRevs(4, 8, rpm)).toEqual({ valid: false, warmupS: 0, endS: 0 });
    }
  });
});
