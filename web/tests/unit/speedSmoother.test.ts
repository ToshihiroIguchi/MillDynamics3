import { describe, expect, it } from "vitest";
import { SpeedSmoother } from "../../src/ui/speedSmoother";

describe("SpeedSmoother", () => {
  it("publishes the first sample immediately, then at most once per refresh interval", () => {
    const s = new SpeedSmoother({ refreshMs: 500 });
    expect(s.push(0, 1, 960, 960)).not.toBeNull();
    expect(s.push(16, 1, 960, 960)).toBeNull();
    expect(s.push(499, 1, 960, 960)).toBeNull();
    expect(s.push(500, 1, 960, 960)).not.toBeNull();
  });

  it("smooths a noisy rate: figures stay near the mean instead of following each frame", () => {
    const s = new SpeedSmoother({ tauMs: 1000, refreshMs: 0 });
    let last = null;
    for (let t = 0; t <= 20000; t += 16) {
      const noisy = t % 32 === 0 ? 0.12 : 0.19;
      last = s.push(t, noisy, noisy * 960, 960) ?? last;
    }
    expect(last!.timeScale).toBeGreaterThan(0.14);
    expect(last!.timeScale).toBeLessThan(0.17);
    expect(last!.subStepsAchieved / last!.subStepsRequired).toBeCloseTo(last!.timeScale, 6);
    expect(last!.slow).toBe(true);
  });

  it("converges to a sustained new rate", () => {
    const s = new SpeedSmoother({ tauMs: 1000, refreshMs: 0 });
    s.push(0, 1, 960, 960);
    let last = null;
    for (let t = 16; t <= 10000; t += 16) last = s.push(t, 0.5, 480, 960) ?? last;
    expect(last!.timeScale).toBeCloseTo(0.5, 2);
  });

  it("does not flip between states on jitter inside the hysteresis band", () => {
    const s = new SpeedSmoother({ tauMs: 1, refreshMs: 0, slowBelow: 0.98, realAtOrAbove: 0.99 });
    expect(s.push(0, 1.0, 960, 960)!.slow).toBe(false);
    expect(s.push(10, 0.985, 945, 960)!.slow).toBe(false);
    expect(s.push(20, 0.97, 931, 960)!.slow).toBe(true);
    expect(s.push(30, 0.985, 945, 960)!.slow).toBe(true);
    expect(s.push(40, 0.995, 955, 960)!.slow).toBe(false);
  });
});
