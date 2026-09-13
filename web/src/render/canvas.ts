// Minimal Canvas 2D renderer. Through M0 this draws only the drum wall and a rotation marker;
// balls/lifters/surface/fluid/dye layers are added in M1/M2/M3 (see docs/PLAN.md ss4.2).

export interface DrumRenderState {
  radiusM: number;
  drumAngle: number;
}

export class CanvasRenderer {
  private readonly ctx: CanvasRenderingContext2D;

  constructor(private readonly canvas: HTMLCanvasElement) {
    const ctx = canvas.getContext("2d");
    if (!ctx) {
      throw new Error("2D canvas context is not available");
    }
    this.ctx = ctx;
  }

  /** Resizes the backing store for the given CSS pixel size, accounting for device pixel ratio. */
  resize(widthCss: number, heightCss: number): void {
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(widthCss * dpr);
    this.canvas.height = Math.round(heightCss * dpr);
    this.canvas.style.width = `${widthCss}px`;
    this.canvas.style.height = `${heightCss}px`;
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  render(state: DrumRenderState): void {
    const { ctx, canvas } = this;
    const dpr = window.devicePixelRatio || 1;
    const widthCss = canvas.width / dpr;
    const heightCss = canvas.height / dpr;

    ctx.fillStyle = "#12161c";
    ctx.fillRect(0, 0, widthCss, heightCss);

    const cx = widthCss / 2;
    const cy = heightCss / 2;
    const pxPerM = (Math.min(widthCss, heightCss) * 0.9) / (2 * state.radiusM);
    const radiusPx = state.radiusM * pxPerM;

    // Drum wall.
    ctx.beginPath();
    ctx.arc(cx, cy, radiusPx, 0, Math.PI * 2);
    ctx.strokeStyle = "#5b6472";
    ctx.lineWidth = 3;
    ctx.stroke();

    // Rotation marker: a radial line fixed to the wall at drumAngle, so its motion visually
    // confirms rotation direction and speed.
    // Math convention is +y up; canvas is +y down, hence the sign flip on the y component.
    const markerX = cx + Math.cos(state.drumAngle) * radiusPx;
    const markerY = cy - Math.sin(state.drumAngle) * radiusPx;
    ctx.beginPath();
    ctx.moveTo(cx, cy);
    ctx.lineTo(markerX, markerY);
    ctx.strokeStyle = "#e8a33d";
    ctx.lineWidth = 2;
    ctx.stroke();
  }
}
