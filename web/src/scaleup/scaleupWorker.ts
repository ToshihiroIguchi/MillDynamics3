// Dedicated worker for the blocking scale-up wasm calls (crates/mill-core/src/scaleup.rs), kept
// off the main thread and separate from the live-simulation worker (src/worker.ts). Avoids the
// "webworker" TS lib for the same reason worker.ts does: `self` is narrowed to the shape used.

import init, {
  defaultParamsJson,
  scaleupBinEdges,
  scaleupIntervalsCsv,
  scaleupMatchingBallDiameter,
  scaleupRunJson,
  scaleupSensitivityJson,
  scaleupSweepJson,
  setPanicHook,
} from "../wasm/mill_wasm.js";
import type { ScaleupReply, ScaleupRequest } from "./protocol";

interface WorkerScope {
  onmessage: ((event: MessageEvent<ScaleupRequest>) => void) | null;
  postMessage: (message: ScaleupReply) => void;
}
const scope = self as unknown as WorkerScope;

function errorText(err: unknown): string {
  if (typeof err === "string") return err;
  if (err instanceof Error) return err.message;
  return String(err);
}

function handle(req: ScaleupRequest): ScaleupReply {
  switch (req.type) {
    case "run":
      return { type: "text", id: req.id, text: scaleupRunJson(req.paramsJson, req.configJson) };
    case "sweep":
      return {
        type: "text",
        id: req.id,
        text: scaleupSweepJson(req.paramsJson, req.configJson, JSON.stringify(req.ratios)),
      };
    case "sensitivity":
      return {
        type: "text",
        id: req.id,
        text: scaleupSensitivityJson(req.paramsJson, req.configJson, req.rel),
      };
    case "matching":
      return {
        type: "number",
        id: req.id,
        value: scaleupMatchingBallDiameter(req.dLabM, req.drumLabM, req.drumTargetM),
      };
    case "intervalsCsv":
      return { type: "text", id: req.id, text: scaleupIntervalsCsv(req.paramsJson, req.configJson) };
  }
}

const ready = init().then(() => {
  setPanicHook();
  scope.postMessage({
    type: "ready",
    defaultParamsJson: defaultParamsJson(),
    binEdges: Array.from(scaleupBinEdges()),
  });
});

ready.catch((err) => scope.postMessage({ type: "error", id: -1, message: `wasm init failed: ${errorText(err)}` }));

scope.onmessage = (event) => {
  const req = event.data;
  // Requests sent before init finished simply wait for it.
  void ready.then(() => {
    try {
      scope.postMessage(handle(req));
    } catch (err) {
      scope.postMessage({ type: "error", id: req.id, message: errorText(err) });
    }
  }, () => undefined);
};
