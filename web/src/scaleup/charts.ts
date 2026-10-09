// Minimal dependency-free canvas line chart for the scale-up panel (see also ui/sparkline.ts).
// Text/grid colours are read from the panel's CSS variables at draw time, so light/dark themes
// work; series colours are a fixed palette that is legible on both backgrounds.

export interface Series {
  name: string;
  xs: number[];
  ys: number[];
}

export interface ChartOptions {
  logX?: boolean;
  logY?: boolean;
  xLabel: string;
  yLabel: string;
  height?: number;
}

/** Okabe-Ito based palette (colour-blind safe, readable on light and dark backgrounds). */
export const SERIES_COLORS = ["#0072b2", "#e69f00", "#009e73", "#cc79a7", "#d55e00", "#56b4e9", "#f0e442", "#999999"];

export function seriesColor(i: number): string {
  return SERIES_COLORS[i % SERIES_COLORS.length] ?? "#999999";
}

const PAD_L = 58;
const PAD_R = 12;
const PAD_T = 10;
const PAD_B = 34;

function cssVar(el: HTMLElement, name: string, fallback: string): string {
  const v = getComputedStyle(el).getPropertyValue(name).trim();
  return v || fallback;
}

/** "Nice" tick values for a linear axis. */
export function niceTicks(min: number, max: number, count = 5): number[] {
  if (!(max > min)) return [min];
  const raw = (max - min) / count;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = (norm < 1.5 ? 1 : norm < 3.5 ? 2 : norm < 7.5 ? 5 : 10) * mag;
  const out: number[] = [];
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-9; v += step) out.push(v);
  return out;
}

function fmtTick(v: number): string {
  if (v === 0) return "0";
  const a = Math.abs(v);
  return a >= 1e4 || a < 1e-2 ? v.toExponential(0) : String(Number(v.toPrecision(3)));
}

/** Draws the chart into `canvas`, sized to its CSS width. Series with no finite positive points on a log axis are skipped. */
export function drawLineChart(canvas: HTMLCanvasElement, series: Series[], opts: ChartOptions): void {
  const dpr = window.devicePixelRatio || 1;
  const cssW = Math.max(240, canvas.clientWidth || canvas.parentElement?.clientWidth || 480);
  const cssH = opts.height ?? 220;
  if (canvas.width !== Math.round(cssW * dpr) || canvas.height !== Math.round(cssH * dpr)) {
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
  }
  canvas.style.height = `${cssH}px`;
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, cssH);

  const text = cssVar(canvas, "--su-text", "#b7bfcc");
  const grid = cssVar(canvas, "--su-grid", "#2a313b");

  const usable = (x: number, y: number) =>
    Number.isFinite(x) && Number.isFinite(y) && (!opts.logX || x > 0) && (!opts.logY || y > 0);
  const pts = series.map((s) => s.xs.map((x, i) => [x, s.ys[i] ?? Number.NaN] as const).filter(([x, y]) => usable(x, y)));
  const all = pts.flat();
  ctx.font = "11px system-ui, sans-serif";
  ctx.fillStyle = text;
  if (all.length === 0) {
    ctx.fillText("no data", PAD_L, cssH / 2);
    return;
  }
  const tx = (v: number) => (opts.logX ? Math.log10(v) : v);
  const ty = (v: number) => (opts.logY ? Math.log10(v) : v);
  let x0 = Math.min(...all.map((p) => tx(p[0])));
  let x1 = Math.max(...all.map((p) => tx(p[0])));
  let y0 = Math.min(...all.map((p) => ty(p[1])));
  let y1 = Math.max(...all.map((p) => ty(p[1])));
  if (!opts.logY) y0 = Math.min(y0, 0);
  if (opts.logX) {
    x0 = Math.floor(x0);
    x1 = Math.ceil(x1);
  }
  if (opts.logY) {
    y0 = Math.floor(y0);
    y1 = Math.ceil(y1);
  }
  if (x1 - x0 < 1e-12) x1 = x0 + 1;
  if (y1 - y0 < 1e-12) y1 = y0 + 1;
  const plotW = cssW - PAD_L - PAD_R;
  const plotH = cssH - PAD_T - PAD_B;
  const px = (v: number) => PAD_L + ((tx(v) - x0) / (x1 - x0)) * plotW;
  const py = (v: number) => PAD_T + plotH - ((ty(v) - y0) / (y1 - y0)) * plotH;

  const ticksFor = (lo: number, hi: number, log: boolean): number[] => {
    if (!log) return niceTicks(lo, hi);
    const out: number[] = [];
    for (let e = Math.ceil(lo); e <= Math.floor(hi); e++) out.push(Math.pow(10, e));
    return out;
  };
  ctx.strokeStyle = grid;
  ctx.lineWidth = 1;
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (const v of ticksFor(y0, y1, !!opts.logY)) {
    const y = py(v);
    ctx.beginPath();
    ctx.moveTo(PAD_L, y);
    ctx.lineTo(cssW - PAD_R, y);
    ctx.stroke();
    ctx.fillText(fmtTick(v), PAD_L - 5, y);
  }
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  for (const v of ticksFor(x0, x1, !!opts.logX)) {
    const x = px(v);
    ctx.beginPath();
    ctx.moveTo(x, PAD_T);
    ctx.lineTo(x, PAD_T + plotH);
    ctx.stroke();
    ctx.fillText(fmtTick(v), x, PAD_T + plotH + 4);
  }
  ctx.fillText(opts.xLabel, PAD_L + plotW / 2, cssH - 14);
  ctx.save();
  ctx.translate(11, PAD_T + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = "middle";
  ctx.fillText(opts.yLabel, 0, 0);
  ctx.restore();

  pts.forEach((line, i) => {
    ctx.strokeStyle = seriesColor(i);
    ctx.lineWidth = 1.8;
    ctx.beginPath();
    line.forEach(([x, y], k) => (k === 0 ? ctx.moveTo(px(x), py(y)) : ctx.lineTo(px(x), py(y))));
    ctx.stroke();
  });
}
