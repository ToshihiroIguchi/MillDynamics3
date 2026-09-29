// Client-side "record the simulation as a video" feature, mirroring report/ module structure.
//
// Frames are encoded with WebCodecs (via mediabunny, which also muxes the WebM container) rather
// than `MediaRecorder`: MediaRecorder stamps frames with the *wall clock*, so a simulation running
// slower (or faster, via `simulation.time_scale`) than real time would produce a video whose
// duration does not match the simulated time. Here every frame is stamped with its *simulation*
// time (`simTime - simTimeAtStart`), so the video always plays back in real-time physical scale:
// 10 s of simulated time == a 10 s video, however long the computation took. Paused periods
// (simTime not advancing) add no frames and therefore no video time.
//
// Frames are captured at most `VIDEO_FPS` per simulated second, from the same canvas the user is
// watching (one call to `captureFrame` per rendered frame).
//
// Container/codec: MP4 (H.264) preferred for compatibility; WebM (VP9/VP8/AV1) if H.264 cannot be
// encoded.
// If WebCodecs is unavailable the feature disables itself (see `isSupported`) rather than throwing.

import {
  BufferTarget,
  CanvasSource,
  Mp4OutputFormat,
  Output,
  WebMOutputFormat,
  getFirstEncodableVideoCodec,
  type VideoCodec,
} from "mediabunny";

/** Video frame rate, in frames per simulated second. */
export const VIDEO_FPS = 30;

interface Container {
  codecs: VideoCodec[];
  extension: "mp4" | "webm";
  mime: string;
  format: () => Mp4OutputFormat | WebMOutputFormat;
}

/** Ordered by preference: MP4/H.264 plays everywhere (Windows Media Player, PowerPoint, iOS,
 * browsers); WebM is the fallback for browsers that cannot encode H.264 (e.g. some Firefox builds). */
const CONTAINERS: Container[] = [
  { codecs: ["avc"], extension: "mp4", mime: "video/mp4", format: () => new Mp4OutputFormat() },
  { codecs: ["vp9", "vp8", "av1"], extension: "webm", mime: "video/webm", format: () => new WebMOutputFormat() },
];

/** Pure: whether a frame is due, given the simulated time elapsed since recording started and the
 * timestamp of the last captured frame (`null` before the first). The first frame is always due;
 * later ones once at least one frame interval of simulated time has passed. Unit-tested. */
export function isFrameDue(simElapsed: number, lastFrameTime: number | null, fps: number = VIDEO_FPS): boolean {
  if (lastFrameTime === null) return true;
  return simElapsed - lastFrameTime >= 1 / fps - 1e-9;
}

/** Downloads `blob` as `filename`, mirroring ui/metricsPanel.ts's CSV export and
 * report/pdf.ts's `downloadReportPdf` idiom exactly: an object URL, a synthetic `<a download>`
 * click, then revoke. */
function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  try {
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    a.click();
  } finally {
    URL.revokeObjectURL(url);
  }
}

export interface CanvasRecorder {
  /** Feature-detects WebCodecs. Callers (ui/toolbar.ts) disable the Record button entirely rather
   * than letting `start` fail when this is `false`. */
  isSupported(): boolean;
  isRecording(): boolean;
  /** Begins recording `canvas`. No-op if already recording or if `isSupported()` is `false`. */
  start(canvas: HTMLCanvasElement, simTimeAtStart: number): void;
  /** Pushes the canvas's current contents as a frame stamped with `simTime`, if a frame is due.
   * No-op (cheap boolean check) if not currently recording, so callers can call this
   * unconditionally once per rendered frame. */
  captureFrame(simTime: number): void;
  /** Finalizes the recording, builds the Blob, and triggers its download. Resolves once the
   * download has been triggered. No-op (resolves immediately) if not currently recording. */
  stop(simTimeAtStop: number): Promise<void>;
}

/** Default `CanvasRecorder` implementation, backed by WebCodecs + mediabunny. See the module doc
 * comment above for the design rationale. */
export class WebCodecsCanvasRecorder implements CanvasRecorder {
  private recording = false;
  private simTimeAtStart = 0;
  private lastFrameTime: number | null = null;
  private output: Output<Mp4OutputFormat | WebMOutputFormat, BufferTarget> | null = null;
  private container: Container = CONTAINERS[0];
  private source: CanvasSource | null = null;
  /** Set while the encoder is still accepting a frame (backpressure): frames are dropped, not
   * queued, so a slow encoder never stalls the render loop. */
  private busy = false;
  /** Resolves when `start`'s async setup (codec probe, `output.start()`) has finished. */
  private ready: Promise<void> = Promise.resolve();
  private pending: Promise<void> = Promise.resolve();

  isSupported(): boolean {
    return typeof VideoEncoder !== "undefined" && typeof VideoFrame !== "undefined";
  }

  isRecording(): boolean {
    return this.recording;
  }

  start(canvas: HTMLCanvasElement, simTimeAtStart: number): void {
    if (this.recording || !this.isSupported()) return;
    this.recording = true;
    this.simTimeAtStart = simTimeAtStart;
    this.lastFrameTime = null;
    this.busy = true; // frames are dropped until setup completes

    this.ready = (async () => {
      let chosen: { container: Container; codec: VideoCodec } | null = null;
      for (const container of CONTAINERS) {
        const codec = await getFirstEncodableVideoCodec(container.codecs, {
          width: canvas.width,
          height: canvas.height,
        });
        if (codec) {
          chosen = { container, codec };
          break;
        }
      }
      if (!chosen) throw new Error("no encodable video codec");
      const { container, codec } = chosen;
      this.container = container;
      const output = new Output({ format: container.format(), target: new BufferTarget() });
      const source = new CanvasSource(canvas, { codec, bitrate: 6_000_000 });
      output.addVideoTrack(source, { frameRate: VIDEO_FPS });
      await output.start();
      this.output = output;
      this.source = source;
      this.busy = false;
    })();
    this.ready.catch((err) => {
      console.error("[recorder] setup failed", err);
      this.recording = false;
    });
  }

  captureFrame(simTime: number): void {
    if (!this.recording || this.busy || !this.source) return;
    const elapsed = Math.max(0, simTime - this.simTimeAtStart);
    if (!isFrameDue(elapsed, this.lastFrameTime)) return;
    this.lastFrameTime = elapsed;
    this.busy = true;
    this.pending = this.source
      .add(elapsed, 1 / VIDEO_FPS)
      .catch((err) => console.error("[recorder] frame encode failed", err))
      .finally(() => {
        this.busy = false;
      });
  }

  async stop(simTimeAtStop: number): Promise<void> {
    if (!this.recording) return;
    this.recording = false;
    try {
      await this.ready;
      await this.pending;
      const { output, source } = this;
      if (!output || !source) return;
      source.close();
      await output.finalize();
      const buffer = output.target.buffer;
      if (buffer) {
        const { mime, extension } = this.container;
        downloadBlob(new Blob([buffer], { type: mime }), `milldynamics-recording-${simTimeAtStop.toFixed(1)}s.${extension}`);
      }
    } finally {
      this.output = null;
      this.source = null;
      this.busy = false;
    }
  }
}
