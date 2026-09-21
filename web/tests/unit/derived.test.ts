import { describe, expect, it } from "vitest";
import {
  criticalSpeedRpm,
  effectiveMedia,
  fluidParticleCountEstimate,
  interstitialFilling,
  substepDisplacementOverDiameter,
  trueBallCount,
  type MediaLike,
  type MillLike,
} from "../../src/params/derived";

// Mirrors mill-core's defaults (crates/mill-core/src/params.rs): MillParams::diameter_m = 1.0,
// MediaParams::ball_diameter_m = 0.010, fill_fraction = 0.30, packing_fraction_2d = 0.82,
// SimulationParams::max_balls = 600 (the `maxBalls` args passed explicitly below, 2000/5000,
// pre-date that default change and are unaffected by it).
const diameterM = 1.0;
const media: MediaLike = { ball_diameter_m: 0.010, fill_fraction: 0.3, packing_fraction_2d: 0.82, density_kg_m3: 6000 };

describe("trueBallCount", () => {
  it("matches the exact default true ball count (2460)", () => {
    expect(trueBallCount(diameterM, media)).toBeCloseTo(2460, 0);
  });
});

describe("effectiveMedia", () => {
  it("coarse-grains when the true count exceeds maxBalls", () => {
    const eff = effectiveMedia(diameterM, media, 2000);
    expect(eff.scaleFactor).toBeCloseTo(Math.sqrt(1.23), 5);
    expect(eff.diameterM).toBeCloseTo(0.01 * Math.sqrt(1.23), 6);
    expect(eff.ballCount).toBe(2000);
    expect(eff.trueDiameterM).toBe(0.01);
  });

  it("leaves media unscaled when the true count is already within maxBalls", () => {
    const eff = effectiveMedia(diameterM, media, 5000);
    expect(eff.scaleFactor).toBe(1);
    expect(eff.diameterM).toBe(eff.trueDiameterM);
    expect(eff.diameterM).toBe(0.01);
    expect(eff.ballCount).toBeCloseTo(2460, 0);
  });

  it("coarse-grains at the project's actual shipped defaults (1.0 m drum, 10 mm balls, max_balls = 600)", () => {
    const eff = effectiveMedia(1.0, media, 600);
    expect(eff.scaleFactor).toBeGreaterThan(1);
    expect(eff.scaleFactor).toBeCloseTo(Math.sqrt(2460 / 600), 2);
  });
});

describe("interstitialFilling", () => {
  it("matches the hand-computed U formula", () => {
    // U = slurry_fill / (media_fill * (1 - packing_fraction_2d)) = 0.35 / (0.30 * 0.18) ~= 6.48
    expect(interstitialFilling(0.35, media)).toBeCloseTo(0.35 / (0.3 * 0.18), 6);
  });

  it("scales linearly with slurry fill fraction", () => {
    expect(interstitialFilling(0.7, media)).toBeCloseTo(2 * interstitialFilling(0.35, media), 10);
  });

  it("returns 0 when the void volume fraction is zero (media.fill_fraction = 0)", () => {
    const zeroFillMedia: MediaLike = { ...media, fill_fraction: 0 };
    expect(interstitialFilling(0.35, zeroFillMedia)).toBe(0);
  });

  it("returns 0 when the void volume fraction is zero (packing_fraction_2d = 1)", () => {
    const fullyPackedMedia: MediaLike = { ...media, packing_fraction_2d: 1 };
    expect(interstitialFilling(0.35, fullyPackedMedia)).toBe(0);
  });
});

describe("fluidParticleCountEstimate", () => {
  it("matches the hand-computed target-area / dx^2 formula", () => {
    const radiusM = 0.5;
    const resolution = 40;
    const slurryFillFraction = 0.15;
    const dx = radiusM / resolution;
    const targetArea = slurryFillFraction * Math.PI * radiusM * radiusM;
    const expected = Math.round(targetArea / (dx * dx));
    expect(fluidParticleCountEstimate(radiusM, resolution, slurryFillFraction)).toBe(expected);
  });
});

describe("substepDisplacementOverDiameter", () => {
  const mill: MillLike = { diameter_m: 1.0, speed_mode: "rpm", speed_value: 30 };

  it("returns a finite, order-of-magnitude-sane positive value for normal inputs", () => {
    const result = substepDisplacementOverDiameter(mill, 4, 0.01109);
    expect(result).toBeGreaterThan(0);
    expect(result).toBeLessThan(10);
  });

  it("returns exactly 0 for a degenerate zero effective diameter", () => {
    expect(substepDisplacementOverDiameter(mill, 4, 0)).toBe(0);
  });

  it("halves when substeps doubles, since subDt is now fixed and independent of time_scale", () => {
    const base = substepDisplacementOverDiameter(mill, 4, 0.01109);
    const doubled = substepDisplacementOverDiameter(mill, 8, 0.01109);
    expect(doubled).toBeCloseTo(base / 2, 10);
  });
});

describe("criticalSpeedRpm", () => {
  it("matches the mill-core formula for a 1 m drum", () => {
    expect(criticalSpeedRpm(1.0)).toBeCloseTo(42.3, 1);
  });
});
