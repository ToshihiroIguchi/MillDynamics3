// Canvas 2D renderer: drum wall, lifters, rotation marker, and the ball (grinding media)
// population. Surface/fluid/dye layers are added in M3 (see docs/PLAN.md ss4.2).

export interface LiftersRenderState {
  count: number;
  heightM: number;
  baseWidthM: number;
  topWidthM: number;
  phaseDeg: number;
}

export interface DrumRenderState {
  radiusM: number;
  drumAngle: number;
  /** Lifters rotate rigidly with the drum; `count <= 0` (the default) means a smooth wall. */
  lifters: LiftersRenderState;
  /** Ball centers, flattened as [x0, y0, x1, y1, ...] (m), in the same world frame as radiusM. */
  ballPositions: Float32Array;
  /** Ball orientations (radians), one per ball, same order as ballPositions. */
  ballOrientations: Float32Array;
  ballRadiusM: number;
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

    this.renderLifters(state, cx, cy, pxPerM);
    this.renderBalls(state, cx, cy, pxPerM);

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

  /** Draws each lifter as a trapezoid rotating with the drum. Vertex layout mirrors
   * `lifter_cross_section_sdf` in crates/mill-core/src/geometry.rs exactly (base at the wall,
   * tapering to a narrower or wider top at `radiusM - heightM`), converted from the lifter's own
   * (radial, tangential) frame into world coordinates. */
  private renderLifters(state: DrumRenderState, cx: number, cy: number, pxPerM: number): void {
    const { lifters, radiusM } = state;
    if (!lifters || lifters.count <= 0) return;
    const { ctx } = this;

    const rTop = radiusM - lifters.heightM;
    const baseHalf = lifters.baseWidthM / 2;
    const topHalf = lifters.topWidthM / 2;
    const phaseRad = (lifters.phaseDeg * Math.PI) / 180;

    ctx.fillStyle = "#6b7280";
    ctx.strokeStyle = "#3a4048";
    ctx.lineWidth = 1;
    for (let i = 0; i < lifters.count; i++) {
      const theta = phaseRad + (i * 2 * Math.PI) / lifters.count + state.drumAngle;
      const urX = Math.cos(theta);
      const urY = Math.sin(theta);
      const utX = -urY;
      const utY = urX;

      const worldVerts: Array<[number, number]> = [
        [urX * radiusM - utX * baseHalf, urY * radiusM - utY * baseHalf],
        [urX * radiusM + utX * baseHalf, urY * radiusM + utY * baseHalf],
        [urX * rTop + utX * topHalf, urY * rTop + utY * topHalf],
        [urX * rTop - utX * topHalf, urY * rTop - utY * topHalf],
      ];

      ctx.beginPath();
      worldVerts.forEach(([wx, wy], idx) => {
        // Same +y-up (world) to +y-down (canvas) flip used throughout this renderer.
        const px = cx + wx * pxPerM;
        const py = cy - wy * pxPerM;
        if (idx === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      });
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    }
  }

  private renderBalls(state: DrumRenderState, cx: number, cy: number, pxPerM: number): void {
    const { ctx } = this;
    const { ballPositions, ballOrientations, ballRadiusM } = state;
    const n = Math.min(ballPositions.length / 2, ballOrientations.length);
    if (n <= 0) return;

    const ballRadiusPx = Math.max(ballRadiusM * pxPerM, 1);

    ctx.fillStyle = "#9aa3b0";
    ctx.strokeStyle = "#3a4048";
    ctx.lineWidth = 1;
    for (let i = 0; i < n; i++) {
      const worldX = ballPositions[i * 2] ?? 0;
      const worldY = ballPositions[i * 2 + 1] ?? 0;
      // Same +y-up (world) to +y-down (canvas) flip as the rotation marker above.
      const px = cx + worldX * pxPerM;
      const py = cy - worldY * pxPerM;

      ctx.beginPath();
      ctx.arc(px, py, ballRadiusPx, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();

      // Spin indicator: a short radial line from the ball's center, at its own orientation, so
      // rolling due to friction is visually confirmable (per-ball theta, world +y-up flipped).
      const theta = ballOrientations[i] ?? 0;
      const spinX = px + Math.cos(theta) * ballRadiusPx;
      const spinY = py - Math.sin(theta) * ballRadiusPx;
      ctx.beginPath();
      ctx.moveTo(px, py);
      ctx.lineTo(spinX, spinY);
      ctx.strokeStyle = "#e0522f";
      ctx.stroke();
      ctx.strokeStyle = "#3a4048";
    }
  }
}
