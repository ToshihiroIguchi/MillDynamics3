// Main-thread side of the repeat-runs worker: lazy start, one request at a time, cancel by
// terminating the worker (a blocking wasm call cannot be interrupted any other way).

import type { BatchWorkerToMain } from "./protocol";

export class BatchCancelled extends Error {
  constructor() {
    super("cancelled");
    this.name = "BatchCancelled";
  }
}

interface Pending {
  resolve: (resultJson: string) => void;
  reject: (err: Error) => void;
}

export class BatchRunner {
  private worker: Worker | null = null;
  private nextId = 1;
  private pending = new Map<number, Pending>();

  private ensureWorker(): Worker {
    if (this.worker) return this.worker;
    const w = new Worker(new URL("./batchWorker.ts", import.meta.url), { type: "module" });
    w.onmessage = (event: MessageEvent<BatchWorkerToMain>) => {
      const msg = event.data;
      const p = this.pending.get(msg.id);
      if (!p) return;
      this.pending.delete(msg.id);
      if (msg.type === "result") p.resolve(msg.resultJson);
      else p.reject(new Error(msg.error));
    };
    w.onerror = (event) => {
      this.failAll(new Error(event.message || "batch worker crashed"));
      this.terminate();
    };
    this.worker = w;
    return w;
  }

  /** Runs one seed. Rejects with the wasm error text, or `BatchCancelled` after `cancel()`. */
  run(paramsJson: string, seed: number, warmup: number, measure: number): Promise<string> {
    const worker = this.ensureWorker();
    const id = this.nextId++;
    return new Promise<string>((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      worker.postMessage({ type: "run", id, paramsJson, seed, warmup, measure });
    });
  }

  cancel(): void {
    this.failAll(new BatchCancelled());
    this.terminate();
  }

  private failAll(err: Error): void {
    const all = [...this.pending.values()];
    this.pending.clear();
    for (const p of all) p.reject(err);
  }

  private terminate(): void {
    this.worker?.terminate();
    this.worker = null;
  }
}
