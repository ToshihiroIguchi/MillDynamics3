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
// not hidden: removing it needs either a fundamentally faster fluid solver -- the 2026-09-27 M6
// solver pass (docs/PERF.md) sped up the existing PBF/coupling implementation substantially
// (constant-factor wins: a counting-sort spatial grid, WASM SIMD, removed per-substep
// allocations) without changing its algorithmic complexity or resolution/quality trade-off, so
// this specific fidelity gap remains -- or much more aggressive ball coarse-graining than any of
// these three tiers use.
//
// As of the 2 mm ball in a 63 mm drum (see `crates/mill-core/src/params.rs`'s
// `shipped_defaults_need_no_coarse_graining` test) `N_true` (~244) sits under every tier's
// `maxBalls`, so no tier coarse-grains -- only `resolution` (fluid particle count) differs. The
// notes below give native (not browser) ms/frame from `examples/perf_probe` with the water default
// (2026-10-08) and the power error from `examples/coarse_graining_probe --reference`; the older
// "achieved speed" browser figures in docs/PERF.md pre-date that default.
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
    label: "Fast (default)",
    maxBalls: 300,
    resolution: 25,
    note: "Native 67 ms/frame (water default). Mill power is about +5 % off a converged reference at this fluid resolution (+19 % at the former 15, +10 % at 20); a thin wetting film is not resolved.",
  },
  {
    id: "balanced",
    label: "Balanced",
    maxBalls: 300,
    resolution: 30,
    note: "Native 72 ms/frame (1143 fluid particles). Mill power within about +-3 % of a converged reference at the 2 mm / 63 mm default with water (docs/VERIFICATION.md); a thin film is still coarse. Coarse-graining the balls (k > 1) is not within +-3 %: about +10-20 % power at k = 2.",
  },
  {
    id: "accuracy",
    label: "Accuracy",
    maxBalls: 600,
    resolution: 50,
    note: "Native 202 ms/frame (3174 fluid particles). The finest tier: mill power within about +-3 % of the converged reference (-1 % at this resolution). Same ball population as every tier at the current default.",
  },
];
