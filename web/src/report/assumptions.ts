// Pure text builders stating the model basis and numerical settings of a run, shared by the PDF
// report, the live-metrics CSV and the batch exports so every output carries the same caveats.
// Read params with the same defensive-cast idiom as report/pdf.ts's `readRunConfig`: every field
// is optional and falls back to its `Params::default()` value.

import type { ParamsJson } from "../protocol";

/** Implied 3D solid fraction of the one-ball-thick slab basis: phi_3D = (2/3) * packing_fraction_2d. */
export function impliedSolidFraction3d(params: ParamsJson | null): number {
  const media = params?.media as { packing_fraction_2d?: number } | undefined;
  return (2 / 3) * (media?.packing_fraction_2d ?? 0.82);
}

export function modelBasisLines(params: ParamsJson | null): string[] {
  const phi2d = (params?.media as { packing_fraction_2d?: number } | undefined)?.packing_fraction_2d ?? 0.82;
  const phi3d = impliedSolidFraction3d(params);
  return [
    "2D cross-section model: values marked /m are per metre of drum length and must not be scaled by a real drum length.",
    "Stress energy SE uses the mass of the 3D sphere each disc stands for.",
    `Per-ball powers assume a one-ball-thick slab (depth = ball diameter), i.e. an implied 3D solid fraction phi_3D = (2/3) x packing_fraction_2d = ${phi3d.toFixed(3)} (packing_fraction_2d = ${phi2d.toFixed(2)}).`,
    "Ratios between conditions are independent of the depth basis.",
    "Impact statistics require coarse-graining k = 1.",
    "Contact coefficients are uncalibrated: use results for relative comparison between conditions.",
  ];
}

export function numericsLines(params: ParamsJson | null): string[] {
  const sim = params?.simulation as
    | {
        substeps?: number;
        dem_iterations?: number;
        pbf_iterations?: number;
        resolution?: number;
        coarse_graining_mode?: string;
        coarse_graining_k?: number;
        seed?: number;
      }
    | undefined;
  const media = params?.media as
    | {
        friction_ball_ball?: number;
        friction_ball_wall?: number;
        restitution_ball_ball?: number;
        restitution_ball_wall?: number;
      }
    | undefined;
  const mode = sim?.coarse_graining_mode ?? "off";
  return [
    `Sub-steps per frame: ${sim?.substeps ?? 8}`,
    `Ball solver iterations: ${sim?.dem_iterations ?? 2}`,
    `Fluid solver iterations: ${sim?.pbf_iterations ?? 3}`,
    `Fluid resolution: ${sim?.resolution ?? 40}`,
    `Coarse-graining mode: ${mode}${mode === "manual" ? `, k = ${sim?.coarse_graining_k ?? 1}` : ""}`,
    `Friction ball-ball / ball-wall: ${media?.friction_ball_ball ?? 0.25} / ${media?.friction_ball_wall ?? 0.35}`,
    `Restitution ball-ball / ball-wall: ${media?.restitution_ball_ball ?? 0.7} / ${media?.restitution_ball_wall ?? 0.5}`,
    `Random seed: ${sim?.seed ?? 1}`,
  ];
}
