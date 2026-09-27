import { describe, expect, it } from "vitest";
import { FrameRequestGate } from "../../src/frameGate";

describe("FrameRequestGate", () => {
  it("allows up to maxInFlight requests, then withholds until a settle frees a slot", () => {
    const gate = new FrameRequestGate(2);
    const first = gate.tryIssue(0, 1 / 60);
    const second = gate.tryIssue(16, 1 / 60);
    expect(first).not.toBeNull();
    expect(second).not.toBeNull();
    expect(gate.tryIssue(32, 1 / 60)).toBeNull();

    gate.settle(first!.requestId);
    const third = gate.tryIssue(48, 1 / 60);
    expect(third).not.toBeNull();
    expect(gate.tryIssue(64, 1 / 60)).toBeNull();
  });

  it("respects a custom maxInFlight and defaults to 2", () => {
    const gate = new FrameRequestGate();
    expect(gate.tryIssue(0, 1 / 60)).not.toBeNull();
    expect(gate.tryIssue(1, 1 / 60)).not.toBeNull();
    expect(gate.tryIssue(2, 1 / 60)).toBeNull();

    const single = new FrameRequestGate(1);
    expect(single.tryIssue(0, 1 / 60)).not.toBeNull();
    expect(single.tryIssue(1, 1 / 60)).toBeNull();
  });

  it("settle() on an unknown or already-settled id is a no-op", () => {
    const gate = new FrameRequestGate(1);
    const req = gate.tryIssue(0, 1 / 60);
    expect(req).not.toBeNull();

    gate.settle(9999); // unknown id: no-op, must not free the one real in-flight slot
    expect(gate.tryIssue(1, 1 / 60)).toBeNull();

    gate.settle(req!.requestId);
    gate.settle(req!.requestId); // already-settled: second call is a no-op
    expect(gate.tryIssue(2, 1 / 60)).not.toBeNull();
  });

  it("falls back to rafDtS for wallDt on the very first tryIssue", () => {
    const gate = new FrameRequestGate();
    const req = gate.tryIssue(1000, 0.0167);
    expect(req).not.toBeNull();
    expect(req!.wallDt).toBe(0.0167);
  });

  it("falls back to rafDtS right after resetClock()", () => {
    const gate = new FrameRequestGate();
    gate.tryIssue(0, 1 / 60);
    gate.resetClock();
    const req = gate.tryIssue(500, 0.02);
    expect(req).not.toBeNull();
    expect(req!.wallDt).toBe(0.02);
  });

  it("computes wallDt as send-to-send elapsed time on a subsequent tryIssue", () => {
    const gate = new FrameRequestGate(10);
    gate.tryIssue(1000, 1 / 60);
    const req = gate.tryIssue(1250, 1 / 60);
    expect(req).not.toBeNull();
    expect(req!.wallDt).toBeCloseTo((1250 - 1000) / 1000, 10);
  });

  it("resetClock() does not clear outstanding in-flight ids", () => {
    const gate = new FrameRequestGate(2);
    const first = gate.tryIssue(0, 1 / 60);
    expect(first).not.toBeNull();
    gate.resetClock();
    // Only 1 of 2 slots is free; the cap still applies to the id issued before resetClock.
    const second = gate.tryIssue(16, 1 / 60);
    expect(second).not.toBeNull();
    expect(gate.tryIssue(32, 1 / 60)).toBeNull();

    // Settling the pre-reset id still frees its slot normally.
    gate.settle(first!.requestId);
    expect(gate.tryIssue(48, 1 / 60)).not.toBeNull();
  });
});
