// Minimal, dependency-free sparkline renderer for the metrics panel (ui/metricsPanel.ts). No
// axes/labels -- just a compact trend line plus a dot marking the latest value, matching this
// app's existing fixed dark palette (styles.css).

const CSS_WIDTH = 120;
const CSS_HEIGHT = 28;
const PAD = 2;

/** (Re)sizes `canvas`'s backing store for the sparkline's fixed CSS footprint, accounting for
 * device pixel ratio, only when it doesn't already match (avoids thrashing the backing store on
 * every redraw). */
function ensureSized(canvas: HTMLCanvasElement): CanvasRenderingContext2D | null {
  const dpr = window.devicePixelRatio || 1;
  const wantW = Math.round(CSS_WIDTH * dpr);
  const wantH = Math.round(CSS_HEIGHT * dpr);
  if (canvas.width !== wantW || canvas.height !== wantH) {
    canvas.width = wantW;
    canvas.height = wantH;
    canvas.style.width = `${CSS_WIDTH}px`;
    canvas.style.height = `${CSS_HEIGHT}px`;
  }
  const ctx = canvas.getContext("2d");
  if (!ctx) return null;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return ctx;
}

export interface SparklineOptions {
  /** Fixed y-axis min/max; omit either to auto-range from the finite values in `values`. */
  min?: number;
  max?: number;
}

/** Draws a single-stroke trend line for `values` (oldest first; `NaN`/non-finite entries render
 * as a gap in the line) into `canvas`, with a dot at the last finite value. No-op (blank) if
 * fewer than two finite values are present. */
export function drawSparkline(canvas: HTMLCanvasElement, values: number[], options: SparklineOptions = {}): void {
  const ctx = ensureSized(canvas);
  if (!ctx) return;
  ctx.clearRect(0, 0, CSS_WIDTH, CSS_HEIGHT);

  const finite = values.filter((v) => Number.isFinite(v));
  if (finite.length < 2) return;

  let min = options.min ?? Math.min(...finite);
  let max = options.max ?? Math.max(...finite);
  if (max - min < 1e-9) {
    // A perfectly flat series still gets a visible centred line rather than collapsing to a
    // divide-by-zero-flavoured NaN plot.
    min -= 0.5;
    max += 0.5;
  }

  const n = values.length;
  const plotW = CSS_WIDTH - 2 * PAD;
  const plotH = CSS_HEIGHT - 2 * PAD;
  const xAt = (i: number) => PAD + (n <= 1 ? 0 : (i / (n - 1)) * plotW);
  const yAt = (v: number) => CSS_HEIGHT - PAD - ((v - min) / (max - min)) * plotH;

  ctx.strokeStyle = "#3a82a8";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  let penDown = false;
  let lastFiniteIdx = -1;
  values.forEach((v, i) => {
    if (!Number.isFinite(v)) {
      penDown = false;
      return;
    }
    const x = xAt(i);
    const y = yAt(v);
    if (!penDown) {
      ctx.moveTo(x, y);
      penDown = true;
    } else {
      ctx.lineTo(x, y);
    }
    lastFiniteIdx = i;
  });
  ctx.stroke();

  if (lastFiniteIdx >= 0) {
    const lastValue = values[lastFiniteIdx];
    if (lastValue !== undefined) {
      ctx.beginPath();
      ctx.arc(xAt(lastFiniteIdx), yAt(lastValue), 2, 0, Math.PI * 2);
      ctx.fillStyle = "#e8a33d";
      ctx.fill();
    }
  }
}
