// Credit-based flow control for `requestFrame` messages sent from main.ts's rAF loop to
// worker.ts, so the worker's FIFO mailbox never accumulates an unbounded backlog of unanswered
// `requestFrame` messages (which would also starve `pause`/`setParams`/`init`, stuck behind that
// backlog). See docs/PLAN.md ss4.1.
//
// This class is pure and has no DOM/worker dependencies, so it can be unit-tested directly (see
// tests/unit/frameGate.test.ts).

const DEFAULT_MAX_IN_FLIGHT = 2;

export class FrameRequestGate {
  private readonly maxInFlight: number;
  private readonly inFlight = new Set<number>();
  private nextRequestId = 0;
  private lastIssueMs: number | null = null;

  constructor(maxInFlight: number = DEFAULT_MAX_IN_FLIGHT) {
    this.maxInFlight = maxInFlight;
  }

  /** Called once per rAF tick when the sim should be running. `nowMs` is the current high-res
   * timestamp; `rafDtS` is this tick's wall-clock delta in seconds (from the caller's own rAF
   * timing), used as a fallback wallDt when there's no prior send to measure from. Returns the
   * {requestId, wallDt} to send as a requestFrame message, or null if already at the in-flight cap
   * (caller should send nothing this tick). */
  tryIssue(nowMs: number, rafDtS: number): { requestId: number; wallDt: number } | null {
    if (this.inFlight.size >= this.maxInFlight) return null;
    const wallDt = this.lastIssueMs !== null ? (nowMs - this.lastIssueMs) / 1000 : rafDtS;
    this.lastIssueMs = nowMs;
    const requestId = this.nextRequestId;
    this.nextRequestId += 1;
    this.inFlight.add(requestId);
    return { requestId, wallDt };
  }

  /** Call when a `frame` (with a requestId) or `frameSkipped` message arrives, to free up an
   * in-flight slot. Unknown/already-settled ids are ignored (no-op), which matters after
   * resetClock() is called mid-flight -- stale replies for ids issued before a pause/reset must
   * not corrupt state. */
  settle(requestId: number): void {
    this.inFlight.delete(requestId);
  }

  /** Call on pause and on any reset (Reset button, or a reset-path Apply that sends "init").
   * Forgets the last-issue wall-clock timestamp so the NEXT tryIssue's wallDt doesn't include the
   * paused/reset interval -- it should fall back to that call's rafDtS instead, as if issuing for
   * the first time. Does NOT touch outstanding in-flight ids; those still get settled normally
   * when their replies arrive (the FIFO worker will still answer them), so credits return to the
   * pool naturally rather than needing to be force-cleared. */
  resetClock(): void {
    this.lastIssueMs = null;
  }
}
