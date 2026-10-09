// Main-thread handle to the scale-up worker: promise-based request/reply with ids, and `cancel()`
// which terminates the worker (the wasm call is blocking, so this is the only way to stop it) and
// starts a fresh one.

import type { ScaleupReply, ScaleupRequest, WithoutId } from "./protocol";

export class CancelledError extends Error {
  constructor() {
    super("cancelled");
    this.name = "CancelledError";
  }
}

export interface ScaleupReady {
  defaultParamsJson: string;
  binEdges: number[];
}

interface Pending {
  resolve: (reply: ScaleupReply) => void;
  reject: (err: Error) => void;
}

export class ScaleupClient {
  private worker!: Worker;
  private nextId = 1;
  private pending = new Map<number, Pending>();
  /** Resolves when the current worker's wasm has initialised. */
  ready!: Promise<ScaleupReady>;
  private readyResolve!: (r: ScaleupReady) => void;
  private readyReject!: (e: Error) => void;

  constructor() {
    this.spawn();
  }

  private spawn(): void {
    this.ready = new Promise<ScaleupReady>((resolve, reject) => {
      this.readyResolve = resolve;
      this.readyReject = reject;
    });
    // Avoid an unhandled rejection when nobody awaits `ready` before a cancel.
    this.ready.catch(() => undefined);
    this.worker = new Worker(new URL("./scaleupWorker.ts", import.meta.url), { type: "module" });
    this.worker.onmessage = (event: MessageEvent<ScaleupReply>) => this.onReply(event.data);
    this.worker.onerror = (event) => this.failAll(new Error(event.message || "scale-up worker error"));
  }

  private failAll(err: Error): void {
    this.readyReject(err);
    for (const p of this.pending.values()) p.reject(err);
    this.pending.clear();
  }

  private onReply(reply: ScaleupReply): void {
    if (reply.type === "ready") {
      this.readyResolve({ defaultParamsJson: reply.defaultParamsJson, binEdges: reply.binEdges });
      return;
    }
    if (reply.type === "error" && reply.id < 0) {
      this.failAll(new Error(reply.message));
      return;
    }
    const p = this.pending.get(reply.id);
    if (!p) return;
    this.pending.delete(reply.id);
    if (reply.type === "error") p.reject(new Error(reply.message));
    else p.resolve(reply);
  }

  private request(req: WithoutId<ScaleupRequest>): Promise<ScaleupReply> {
    const id = this.nextId++;
    return new Promise<ScaleupReply>((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.worker.postMessage({ ...req, id } as ScaleupRequest);
    });
  }

  private async text(req: WithoutId<ScaleupRequest>): Promise<string> {
    const r = await this.request(req);
    if (r.type !== "text") throw new Error(`unexpected reply type ${r.type}`);
    return r.text;
  }

  run(paramsJson: string, configJson: string): Promise<string> {
    return this.text({ type: "run", paramsJson, configJson });
  }

  sweep(paramsJson: string, configJson: string, ratios: number[]): Promise<string> {
    return this.text({ type: "sweep", paramsJson, configJson, ratios });
  }

  sensitivity(paramsJson: string, configJson: string, rel: number): Promise<string> {
    return this.text({ type: "sensitivity", paramsJson, configJson, rel });
  }

  intervalsCsv(paramsJson: string, configJson: string): Promise<string> {
    return this.text({ type: "intervalsCsv", paramsJson, configJson });
  }

  async matchingBallDiameter(dLabM: number, drumLabM: number, drumTargetM: number): Promise<number> {
    const r = await this.request({ type: "matching", dLabM, drumLabM, drumTargetM });
    if (r.type !== "number") throw new Error(`unexpected reply type ${r.type}`);
    return r.value;
  }

  /** Terminates the worker, rejects every pending request with `CancelledError`, starts a new worker. */
  cancel(): void {
    this.worker.terminate();
    const err = new CancelledError();
    this.readyReject(err);
    for (const p of this.pending.values()) p.reject(err);
    this.pending.clear();
    this.spawn();
  }
}
