// Media material presets: pre-measured `media.density_kg_m3` values for a few common grinding
// media materials, driving the "Media material" dropdown convenience in the params panel (mirrors
// presets.ts's QUALITY_PRESETS pattern for simulation.max_balls/resolution).
//
// These are representative/typical densities for each material family, not a specific certified
// grade's datasheet value -- real density varies with alloy/purity/manufacturing process, as noted
// per-entry below. Pick the closest match and fall back to the raw "Media density" number input
// (the "(custom)" option) for a known exact value.
export interface MediaMaterial {
  id: string;
  label: string;
  densityKgM3: number;
  /** One-line note on the value's provenance/caveat. */
  note: string;
}

export const MEDIA_MATERIALS: MediaMaterial[] = [
  {
    id: "ysz",
    label: "Yttria-stabilized zirconia (YSZ)",
    densityKgM3: 6000,
    note: "Matches this project's default media density (models ZrO2 ceramic); see docs/PARAMETERS.md.",
  },
  {
    id: "stainless_steel",
    label: "Stainless steel",
    densityKgM3: 7700,
    note: "Representative figure -- real stainless steel grinding media vary ~7700-7900 kg/m3 by grade/alloy.",
  },
  {
    id: "alumina",
    label: "Alumina (Al2O3)",
    densityKgM3: 3600,
    note: "Representative figure -- real alumina media vary ~3400-3900 kg/m3 depending on alumina purity grade.",
  },
  {
    id: "mullite",
    label: "Mullite",
    densityKgM3: 3000,
    note: "Representative figure -- real mullite media vary ~2800-3200 kg/m3.",
  },
];
