// Client-side "record the simulation as a video" feature, mirroring report/ module structure.
//
// Uses `HTMLCanvasElement.captureStream()` + `MediaRecorder` (both native browser APIs, zero new
// dependencies) instead of an offline re-simulation/ffmpeg.wasm approach: this records exactly the
// frames the user is actually watching, which is correct both for normal playback and for
// `simulation.time_scale` running faster/slower than real time -- the video shows what was
// rendered, not a separately-computed deterministic replay.
//
// Frame-accurate capture: the canvas's stream is created in "manual" mode (`captureStream(0)`, see
// MDN), so it only advances when `track.requestFrame()` is called explicitly (from main.ts's
// `frameLoop`, right after `renderer.render(...)`) rather than being sampled on MediaRecorder's own
// wall-clock timer. This keeps the recorded video's frame timing tied to what was actually drawn,
// independent of `requestAnimationFrame` jitter.
//
// Container/codec: WebM, VP9 preferred with VP8 fallback -- natively supported by every browser
// this project targets (Chrome/Firefox/Edge) with no codec license concerns (unlike MP4/H.264).
// The candidate list + `MediaRecorder.isTypeSupported` check is a small ordered list, tried in
// order; if none are supported the feature disables itself (see `isSupported`) rather than
// throwing.

/** Ordered by preference: VP9 first (better compression), VP8 fallback, then an unspecified codec
 * within the WebM container as a last resort. */
export const MIME_TYPE_CANDIDATES: string[] = ["video/webm;codecs=vp9", "video/webm;codecs=vp8", "video/webm"];

/** Pure (no DOM/MediaRecorder access) so it can be unit-tested with a fake `isTypeSupported` --
 * see tests/unit/recorder.test.ts. Returns the first candidate `isTypeSupported` accepts, or
 * `null` if none are (including an empty candidate list). */
export function pickSupportedMimeType(candidates: string[], isTypeSupported: (type: string) => boolean): string | null {
  for (const candidate of candidates) {
    if (isTypeSupported(candidate)) return candidate;
  }
  return null;
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
  /** Feature-detects `captureStream` + `MediaRecorder` + at least one supported mime type from
   * `MIME_TYPE_CANDIDATES`. Callers (ui/toolbar.ts) disable the Record button entirely rather than
   * letting `start` throw when this is `false`. */
  isSupported(): boolean;
  isRecording(): boolean;
  /** Begins recording `canvas`. No-op if already recording or if `isSupported()` is `false`. */
  start(canvas: HTMLCanvasElement, simTimeAtStart: number): void;
  /** Pushes one frame onto the recording. No-op (cheap boolean check) if not currently recording,
   * so callers can call this unconditionally once per rendered frame. */
  captureFrame(): void;
  /** Finalizes the recording, builds the Blob, and triggers its download. Resolves once the
   * download has been triggered. No-op (resolves immediately) if not currently recording. */
  stop(simTimeAtStop: number): Promise<void>;
}

/** Default `CanvasRecorder` implementation, backed by `HTMLCanvasElement.captureStream()` +
 * `MediaRecorder`. See the module doc comment above for the design rationale. */
export class MediaRecorderCanvasRecorder implements CanvasRecorder {
  private recording = false;
  private mediaRecorder: MediaRecorder | null = null;
  // `captureStream()`'s video track is a `CanvasCaptureMediaStreamTrack` at runtime (the subtype
  // that carries `requestFrame()` for manual-mode capture); `MediaStream.getVideoTracks()` is
  // typed to return the broader `MediaStreamTrack[]`, so this narrows it explicitly.
  private track: CanvasCaptureMediaStreamTrack | null = null;
  private stream: MediaStream | null = null;
  private chunks: Blob[] = [];
  private mimeType: string | null = null;

  isSupported(): boolean {
    if (typeof MediaRecorder === "undefined") return false;
    if (typeof HTMLCanvasElement === "undefined" || typeof HTMLCanvasElement.prototype.captureStream !== "function") {
      return false;
    }
    return pickSupportedMimeType(MIME_TYPE_CANDIDATES, (type) => MediaRecorder.isTypeSupported(type)) !== null;
  }

  isRecording(): boolean {
    return this.recording;
  }

  start(canvas: HTMLCanvasElement, _simTimeAtStart: number): void {
    if (this.recording || !this.isSupported()) return;

    const mimeType = pickSupportedMimeType(MIME_TYPE_CANDIDATES, (type) => MediaRecorder.isTypeSupported(type));
    if (!mimeType) return; // isSupported() already checked this, but guard defensively regardless.

    // 0 fps = "manual" mode: the stream only advances on an explicit `track.requestFrame()` call
    // (see `captureFrame` below), not on MediaRecorder's own wall-clock sampling.
    const stream = canvas.captureStream(0);
    const track = stream.getVideoTracks()[0] as CanvasCaptureMediaStreamTrack | undefined;
    if (!track) return;

    const mediaRecorder = new MediaRecorder(stream, { mimeType });
    this.chunks = [];
    mediaRecorder.ondataavailable = (event: BlobEvent) => {
      if (event.data.size > 0) this.chunks.push(event.data);
    };
    mediaRecorder.start();

    this.mediaRecorder = mediaRecorder;
    this.track = track;
    this.stream = stream;
    this.mimeType = mimeType;
    this.recording = true;
  }

  captureFrame(): void {
    if (!this.recording || !this.track) return;
    this.track.requestFrame();
  }

  stop(simTimeAtStop: number): Promise<void> {
    if (!this.recording || !this.mediaRecorder) return Promise.resolve();

    const mediaRecorder = this.mediaRecorder;
    const mimeType = this.mimeType ?? "video/webm";
    const stream = this.stream;

    return new Promise((resolve) => {
      mediaRecorder.onstop = () => {
        const blob = new Blob(this.chunks, { type: mimeType });
        downloadBlob(blob, `milldynamics-recording-${simTimeAtStop.toFixed(1)}s.webm`);

        stream?.getTracks().forEach((t) => t.stop());
        this.mediaRecorder = null;
        this.track = null;
        this.stream = null;
        this.chunks = [];
        this.recording = false;
        resolve();
      };
      mediaRecorder.stop();
    });
  }
}
