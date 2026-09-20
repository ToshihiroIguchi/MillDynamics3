// Quality presets: pre-measured (`simulation.max_balls`, `simulation.resolution`) pairs trading
// ball/fluid particle count for wall-clock cost, per docs/PERF.md's methodology. Both fields are
// `resetRequired` (schema.ts), so choosing a preset only edits the form -- the user still clicks
// Apply to commit it (same flow as any other field edit).
//
// Measured in-browser (this machine, WASM release build, default drum/media/slurry otherwise,
// `lifters.count = 0`, steady cascading after ~10 s of sim time; see docs/PERF.md for the full
// methodology and raw readings, including the exact "fluid spacing / ball diameter" readout this
// same panel shows live -- `dx / effective_ball_diameter`, equivalently `h / (2 * diameter)`).
// "Achieved speed" is real, not a target -- Realtime is the only tier that reaches >= 1.0x. All
// three tiers stay under the panel's own "a fluid particle is wider than a ball" warning (that
// ratio is < 1 for each), but all three are still coarse relative to what a visually smooth *thin*
// wetting film specifically needs -- that is a different, stricter bar the warning threshold does
// not check, and none of these three tiers clears it. That trade-off is real and disclosed here,
// not hidden: removing it needs either a fundamentally faster fluid solver (the M6 performance
// pass this project has not yet done) or much more aggressive ball coarse-graining than any of
// these three tiers use.
export interface QualityPreset {
  id: string;
  label: string;
  maxBalls: number;
  resolution: number;
  /** One-line, honest summary of this tier's actual measured trade-off (not a marketing claim). */
  note: string;
}

export const QUALITY_PRESETS: QualityPreset[] = [
  {
    id: "realtime",
    label: "Realtime (default)",
    maxBalls: 150,
    resolution: 15,
    note: "~1.0-1.1x achieved speed. Fewest, largest coarse-grained balls and coarsest slurry of the three tiers -- a thin wetting film is not resolved smoothly.",
  },
  {
    id: "balanced",
    label: "Balanced",
    maxBalls: 300,
    resolution: 25,
    note: "~0.7-0.75x achieved speed (noticeably slower than real time). More balls and finer slurry than Realtime; wetting resolution improves but a thin film is still coarse.",
  },
  {
    id: "accuracy",
    label: "Accuracy",
    maxBalls: 600,
    resolution: 40,
    note: "~0.3-0.4x achieved speed (clearly slow motion). This project's original default ball/fluid counts -- the finest of the three tiers, and still under-resolves a thin wetting film.",
  },
];
