import { describe, expect, it } from "vitest";
import { isFrameDue } from "../../src/video/recorder";

describe("isFrameDue", () => {
  it("always captures the first frame", () => {
    expect(isFrameDue(0, null)).toBe(true);
  });

  it("waits one frame interval of simulated time between frames", () => {
    expect(isFrameDue(0.02, 0, 30)).toBe(false);
    expect(isFrameDue(1 / 30, 0, 30)).toBe(true);
    expect(isFrameDue(0.5, 0.4, 30)).toBe(true);
  });

  it("captures nothing while simulated time is not advancing (paused)", () => {
    expect(isFrameDue(1.0, 1.0, 30)).toBe(false);
  });
});
