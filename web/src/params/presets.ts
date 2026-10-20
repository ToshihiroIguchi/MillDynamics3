// Quality presets: pre-measured (`simulation.max_balls`, `simulation.resolution`) pairs trading
// ball/fluid particle count for wall-clock cost, per docs/PERF.md's methodology. Both fields are
// `resetRequired` (schema.ts), so choosing a preset only edits the form -- the user still clicks
// Apply to commit it (same flow as any other field edit).
//
// Measured in-browser (this machine, WASM release build, default drum/media/slurry otherwise,
// `lifters.count = 0`, steady cascading after ~10 s of sim time; see docs/PERF.md for the full
// methodology and raw readings, including the exact "fluid spacing / ball diameter" readout this
// same panel shows live -- `dx / effective_ball_diameter`, equivalently `h / (2 * diameter)`).
// "Achieved speed" is real, not a target -- see each tier's own ms/frame note below. All
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
// Coarse-graining is Off by default and no preset changes the mode except `impact_stats`, which
// sets it to "off" explicitly; `maxBalls` only matters if the user later picks Auto. At the
// 2 mm ball in a 63 mm drum `N_true` (~244) is small anyway, so only `resolution` (fluid
// particle count) differs between tiers. The
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
  /** Optional: `simulation.substeps`. Presets that omit it leave the field untouched. */
  substeps?: number;
  /** Optional: `simulation.dem_iterations`. Presets that omit it leave the field untouched. */
  demIterations?: number;
  /** Optional: `simulation.coarse_graining_mode`. Presets that omit it leave the field untouched. */
  coarseGrainingMode?: string;
}

/** One form field a preset sets: the schema path and the value as the form input would hold it. */
export interface PresetFieldValue {
  path: string;
  value: string;
}

/** Every form path any preset may write; edits to these re-check whether the selected preset still matches. */
export const PRESET_FIELD_PATHS: ReadonlySet<string> = new Set([
  "simulation.max_balls",
  "simulation.resolution",
  "simulation.substeps",
  "simulation.dem_iterations",
  "simulation.coarse_graining_mode",
]);

/**
 * The (path, value) pairs a preset writes into the form. The two required fields come first;
 * optional fields appear only when the preset defines them, so existing presets yield exactly
 * the original two entries.
 */
export function presetFieldValues(preset: QualityPreset): PresetFieldValue[] {
  const out: PresetFieldValue[] = [
    { path: "simulation.max_balls", value: String(preset.maxBalls) },
    { path: "simulation.resolution", value: String(preset.resolution) },
  ];
  if (preset.substeps !== undefined) out.push({ path: "simulation.substeps", value: String(preset.substeps) });
  if (preset.demIterations !== undefined) out.push({ path: "simulation.dem_iterations", value: String(preset.demIterations) });
  if (preset.coarseGrainingMode !== undefined) {
    out.push({ path: "simulation.coarse_graining_mode", value: preset.coarseGrainingMode });
  }
  return out;
}

/**
 * True when every field the preset sets currently holds the preset's value. `read` returns the
 * form's current string value for a path (or undefined when the field is missing). Numeric
 * fields compare numerically, the select field compares as a string.
 */
export function presetMatchesFields(preset: QualityPreset, read: (path: string) => string | undefined): boolean {
  return presetFieldValues(preset).every(({ path, value }) => {
    const current = read(path);
    if (current === undefined) return false;
    return path === "simulation.coarse_graining_mode" ? current === value : Number(current) === Number(value);
  });
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
    note: "Native 72 ms/frame (1143 fluid particles). Mill power within about +-3 % of a converged reference at the 2 mm / 63 mm default with water (checked against a finer-resolution reference run); a thin film is still coarse. Coarse-graining is Off; if chosen explicitly it is not within +-3 % (about +10-20 % power at k = 2).",
  },
  {
    id: "accuracy",
    label: "Accuracy",
    maxBalls: 600,
    resolution: 50,
    note: "Native 202 ms/frame (3174 fluid particles). The finest tier: mill power within about +-3 % of the converged reference (-1 % at this resolution). Same true-size ball population as every tier at the current default (coarse-graining Off).",
  },
  {
    id: "impact_stats",
    label: "Impact statistics (slow)",
    maxBalls: 600,
    resolution: 25,
    substeps: 32,
    demIterations: 12,
    coarseGrainingMode: "off",
    note: "About 4x the sub-steps of the default (32 vs 8, plus 12 ball-solver iterations), so the simulation runs several times slower than real time. The Grinding proxy impact statistics are only reasonably converged this way: impact counts still drift roughly +5 % per doubling of sub-steps and mean stress energies move 3-5 %, so compare runs only with identical settings. Coarse-graining is off (true-size balls), so the ball count can be large for small balls (the panel warns above 3000 and blocks Apply above 30000 true balls).",
  },
];
