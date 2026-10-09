// Typed messages between the main thread and the scale-up worker (scaleupWorker.ts). Every request
// carries an id that its reply echoes; a failure of any wasm call (the wasm functions throw
// strings) comes back as an "error" reply, never swallowed.

export interface RunRequest {
  type: "run";
  id: number;
  paramsJson: string;
  configJson: string;
}

export interface SweepRequest {
  type: "sweep";
  id: number;
  paramsJson: string;
  configJson: string;
  /** Drum-to-ball diameter ratios. */
  ratios: number[];
}

export interface SensitivityRequest {
  type: "sensitivity";
  id: number;
  paramsJson: string;
  configJson: string;
  rel: number;
}

export interface MatchingRequest {
  type: "matching";
  id: number;
  dLabM: number;
  drumLabM: number;
  drumTargetM: number;
}

export interface IntervalsCsvRequest {
  type: "intervalsCsv";
  id: number;
  paramsJson: string;
  configJson: string;
}

export type ScaleupRequest = RunRequest | SweepRequest | SensitivityRequest | MatchingRequest | IntervalsCsvRequest;

/** Distributes `Omit` over a union so each member keeps its own fields. */
export type WithoutId<T> = T extends unknown ? Omit<T, "id"> : never;

export interface ReadyReply {
  type: "ready";
  /** `defaultParamsJson()` text. */
  defaultParamsJson: string;
  binEdges: number[];
}

export interface TextReply {
  type: "text";
  id: number;
  /** JSON text (run/sweep/sensitivity) or CSV text (intervalsCsv). */
  text: string;
}

export interface NumberReply {
  type: "number";
  id: number;
  value: number;
}

export interface ErrorReply {
  type: "error";
  /** Request id, or -1 for a failure outside any request (e.g. wasm init). */
  id: number;
  message: string;
}

export type ScaleupReply = ReadyReply | TextReply | NumberReply | ErrorReply;
