// Dedicated worker for repeat runs (see batch/protocol.ts). Blocks inside `stressRunJson` for the
// whole run, so it must never share a thread with the main simulation worker. Started lazily by
// batch/runner.ts; cancelling terminates this worker and creates a fresh one on the next run.

import init, { setPanicHook, stressRunJson } from "../wasm/mill_wasm.js";
import type { BatchMainToWorker, BatchWorkerToMain } from "./protocol";

interface BatchScope {
  onmessage: ((event: MessageEvent<BatchMainToWorker>) => void) | null;
  postMessage: (message: BatchWorkerToMain) => void;
}
const scope = self as unknown as BatchScope;

let ready: Promise<void> | null = null;
function ensureInit(): Promise<void> {
  ready ??= init().then(() => {
    setPanicHook();
  });
  return ready;
}

scope.onmessage = (event) => {
  const msg = event.data;
  if (msg.type !== "run") return;
  void (async () => {
    try {
      await ensureInit();
      const resultJson = stressRunJson(msg.paramsJson, msg.seed, msg.warmup, msg.measure);
      scope.postMessage({ type: "result", id: msg.id, resultJson });
    } catch (err) {
      // wasm-bindgen throws plain strings for Result::Err; Error objects for panics/init failures.
      scope.postMessage({ type: "error", id: msg.id, error: err instanceof Error ? err.message : String(err) });
    }
  })();
};
