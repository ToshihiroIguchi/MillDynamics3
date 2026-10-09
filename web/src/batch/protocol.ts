// Messages exchanged between the main thread and the repeat-runs worker (batchWorker.ts).
// The worker runs `stressRunJson` (crates/mill-wasm), which is a blocking, headless run of the SAME
// Simulation with the SAME params as the main screen (only the seed is replaced). One seed per call.

export interface BatchRunMessage {
  type: "run";
  /** Correlates the reply with this request. */
  id: number;
  /** Full applied params of the main simulation, serialised. */
  paramsJson: string;
  seed: number;
  /** Whole revolutions simulated before measuring (>= 0). */
  warmup: number;
  /** Whole revolutions averaged over (integer >= 1). */
  measure: number;
}

export interface BatchResultMessage {
  type: "result";
  id: number;
  /** `SeedResult` serialised by Rust (see stats.ts's `SeedResult`). */
  resultJson: string;
}

export interface BatchErrorMessage {
  type: "error";
  id: number;
  /** A thrown wasm error (usually a string) or init failure, as text. Never swallowed. */
  error: string;
}

export type BatchMainToWorker = BatchRunMessage;
export type BatchWorkerToMain = BatchResultMessage | BatchErrorMessage;
