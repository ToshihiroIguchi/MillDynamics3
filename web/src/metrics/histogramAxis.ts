// Axis helpers for the impact stress-energy histogram. Its bin edges are log-spaced in the
// dimensionless SE / (m g d) but reported in joules (scaled by the 3D ball's m g d), so they do
// NOT land on powers of ten; ticks are therefore placed every `BINS_PER_DECADE` edges instead.

export const BINS_PER_DECADE = 4;

/** Indices of the bin edges to label: the first edge and every `BINS_PER_DECADE`-th one after it
 * (one decade apart). Non-positive/non-finite edges are skipped. */
export function tickEdgeIndices(edges: readonly number[]): number[] {
  const out: number[] = [];
  for (let i = 0; i < edges.length; i += BINS_PER_DECADE) {
    const e = edges[i];
    if (e !== undefined && Number.isFinite(e) && e > 0) out.push(i);
  }
  return out;
}

/** Compact joule label such as "1.0e-6". */
export function formatEdgeJ(edge: number): string {
  return edge.toExponential(1);
}
