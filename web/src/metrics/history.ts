// Ring-buffer-style sample history for the metrics panel's sparklines and CSV export
// (ui/metricsPanel.ts). Samples are taken at a fixed *simulation*-time cadence (not wall-clock
// time), so a paused/slow-running sim doesn't over-sample and a fast-forwarded one doesn't
// under-sample relative to what actually happened in the simulated mill.

export interface CsvColumn {
  id: string;
  label: string;
  unit?: string;
}

interface Sample {
  simTime: number;
  values: Record<string, number | null>;
}

export class MetricsHistory {
  private samples: Sample[] = [];
  private lastSampledSimTime = -Infinity;

  constructor(
    private readonly capacity: number,
    private readonly sampleIntervalS: number,
  ) {}

  /** Records one sample if at least `sampleIntervalS` of *simulation* time has elapsed since the
   * last recorded sample (or this is the first sample, or `simTime` moved backward, e.g. a
   * reset -- treated the same as a fresh start). `values` should already be keyed by the same
   * metric ids `seriesFor`/`toCsv` are called with. */
  push(simTime: number, values: Record<string, number | null>): void {
    if (simTime < this.lastSampledSimTime) {
      // Sim time went backward (a reset without an explicit `reset()` call) -- start over rather
      // than mixing pre-/post-reset samples into one history.
      this.samples = [];
      this.lastSampledSimTime = -Infinity;
    }
    if (simTime - this.lastSampledSimTime < this.sampleIntervalS) return;
    this.lastSampledSimTime = simTime;
    this.samples.push({ simTime, values });
    if (this.samples.length > this.capacity) {
      this.samples.shift();
    }
  }

  reset(): void {
    this.samples = [];
    this.lastSampledSimTime = -Infinity;
  }

  /** The recorded values for one metric id, oldest first, `NaN` for a missing/null sample (so
   * callers can treat the return value as a plain numeric series and let `NaN` render as a gap). */
  seriesFor(id: string): number[] {
    return this.samples.map((s) => {
      const v = s.values[id];
      return v == null ? NaN : v;
    });
  }

  get length(): number {
    return this.samples.length;
  }

  /** Renders the full history as CSV: `sim_time_s` first, then one column per entry in
   * `columns`, labelled `"label (unit)"` (or just `label` when `unit` is omitted). A missing/null
   * value renders as an empty cell. */
  toCsv(columns: CsvColumn[]): string {
    const header = ["sim_time_s", ...columns.map((c) => (c.unit ? `${c.label} (${c.unit})` : c.label))];
    const lines = [header.join(",")];
    for (const sample of this.samples) {
      const cells = [sample.simTime.toFixed(3)];
      for (const col of columns) {
        const v = sample.values[col.id];
        cells.push(v == null || Number.isNaN(v) ? "" : String(v));
      }
      lines.push(cells.join(","));
    }
    return lines.join("\n");
  }
}
