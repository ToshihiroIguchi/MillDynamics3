import { describe, expect, it } from "vitest";
import { pickSupportedMimeType } from "../../src/video/recorder";

describe("pickSupportedMimeType", () => {
  it("picks the first supported candidate", () => {
    const candidates = ["video/webm;codecs=vp9", "video/webm;codecs=vp8", "video/webm"];
    const isTypeSupported = (type: string) => type === "video/webm;codecs=vp8" || type === "video/webm";
    expect(pickSupportedMimeType(candidates, isTypeSupported)).toBe("video/webm;codecs=vp8");
  });

  it("returns the only supported candidate when it's first in the list", () => {
    const candidates = ["video/webm;codecs=vp9", "video/webm;codecs=vp8"];
    const isTypeSupported = (type: string) => type === "video/webm;codecs=vp9";
    expect(pickSupportedMimeType(candidates, isTypeSupported)).toBe("video/webm;codecs=vp9");
  });

  it("returns null when none of the candidates are supported", () => {
    const candidates = ["video/webm;codecs=vp9", "video/webm;codecs=vp8", "video/webm"];
    expect(pickSupportedMimeType(candidates, () => false)).toBeNull();
  });

  it("returns null for an empty candidate list", () => {
    expect(pickSupportedMimeType([], () => true)).toBeNull();
  });
});
