import { describe, expect, it } from "vitest";
import { toDegrees, toVerticalDegrees } from "../../src/metrics/angles";

describe("toVerticalDegrees", () => {
  it("maps atan2Rad = PI/2 (world up) to 0 degrees regardless of rotation sign", () => {
    expect(toVerticalDegrees(Math.PI / 2, 1)).toBeCloseTo(0, 5);
    expect(toVerticalDegrees(Math.PI / 2, -1)).toBeCloseTo(0, 5);
  });

  it("maps atan2Rad = 0 to 90 degrees for a positive (CCW) omega", () => {
    // sign = -1, verticalRad = -1 * (0 - PI/2) = PI/2 -> 90 deg.
    expect(toVerticalDegrees(0, 5)).toBeCloseTo(90, 5);
  });

  it("mirrors to 270 degrees for a negative (CW) omega at the same atan2Rad", () => {
    // sign = 1, verticalRad = 1 * (0 - PI/2) = -PI/2 -> -90 deg, wrapped to 270.
    expect(toVerticalDegrees(0, -5)).toBeCloseTo(270, 5);
  });

  it("wraps a result near 2*PI back into [0, 360)", () => {
    const result = toVerticalDegrees(2 * Math.PI - 1e-9, 1);
    expect(result).toBeGreaterThanOrEqual(0);
    expect(result).toBeLessThan(360);
  });
});

describe("toDegrees", () => {
  it("converts 0 rad to 0 deg", () => {
    expect(toDegrees(0)).toBe(0);
  });

  it("converts PI rad to 180 deg", () => {
    expect(toDegrees(Math.PI)).toBeCloseTo(180, 5);
  });

  it("wraps a negative angle into [0, 360)", () => {
    expect(toDegrees(-Math.PI / 2)).toBeCloseTo(270, 5);
  });
});
