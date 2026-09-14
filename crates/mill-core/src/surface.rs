//! Free-surface extraction for rendering.
//!
//! Splats fluid particles onto a scalar occupancy field on a `G x G` grid spanning the drum's
//! bounding box, masks that field with the drum's wall SDF so no bulk-interior value survives
//! outside the wall (or inside a lifter bar), then runs marching squares at half the field's peak
//! (bulk-interior) value to produce closed polylines approximating the slurry free surface,
//! lightly smoothed (one pass of Chaikin corner-cutting) for a less blocky outline, and finally
//! projects any resulting contour point that still landed outside the wall back onto it along the
//! wall normal. See docs/PLAN.md ss3.5.

use std::collections::HashMap;

use glam::Vec2;

use crate::geometry::Drum;
use crate::pbf::FluidParticles;

/// Grid points per axis (docs/PLAN.md ss3.5: "G = 128"); the grid spans `[-R, R]` on each axis
/// with `(GRID_SIZE - 1)` marching-squares cells per axis.
const GRID_SIZE: usize = 128;
/// Splat kernel radius as a multiple of the fluid's own kernel radius `h` (docs/PLAN.md ss3.5).
const SPLAT_RADIUS_FACTOR: f32 = 1.5;
/// Contour threshold as a fraction of the field's peak (bulk-interior) value.
const THRESHOLD_FRACTION: f32 = 0.5;

/// A closed (or, rarely, open if the fluid touches the grid boundary) polyline approximating one
/// connected piece of the free surface.
pub type Polygon = Vec<Vec2>;

/// Smooth, compactly-supported splat kernel (unnormalized -- only the threshold-relative shape
/// matters here, unlike [`crate::pbf`]'s physically-normalized kernels).
fn splat_kernel(r2: f32, radius: f32) -> f32 {
    let radius2 = radius * radius;
    if r2 >= radius2 {
        return 0.0;
    }
    let t = 1.0 - r2 / radius2;
    t * t * t
}

/// Builds the scalar occupancy field on a `GRID_SIZE x GRID_SIZE` grid of points spanning
/// `[-drum.radius_m, drum.radius_m]` on each axis, by splatting every fluid particle with
/// [`splat_kernel`] of radius `SPLAT_RADIUS_FACTOR * fluid.h`, then masking out (zeroing) every
/// grid point that falls outside the wall (or inside a lifter bar) per `drum`'s SDF, evaluated in
/// the drum-local frame via `drum_angle` -- this is what keeps the marching-squares contour from
/// bulging through the wall in the first place; see [`extract_surface`] for the second
/// (contour-point) pass that mops up the remaining sub-cell overshoot.
fn build_field(fluid: &FluidParticles, drum: &Drum, drum_angle: f32) -> Vec<f32> {
    let drum_radius_m = drum.radius_m;
    let mut field = vec![0.0f32; GRID_SIZE * GRID_SIZE];
    if fluid.is_empty() || drum_radius_m <= 0.0 {
        return field;
    }
    let cell_size = 2.0 * drum_radius_m / (GRID_SIZE as f32 - 1.0);
    let splat_radius = SPLAT_RADIUS_FACTOR * fluid.h;

    for &p in &fluid.x {
        // Grid-index bounding box of this particle's splat footprint.
        let gi_min = grid_index(p.x - splat_radius, drum_radius_m, cell_size);
        let gi_max = grid_index(p.x + splat_radius, drum_radius_m, cell_size);
        let gj_min = grid_index(p.y - splat_radius, drum_radius_m, cell_size);
        let gj_max = grid_index(p.y + splat_radius, drum_radius_m, cell_size);
        for gj in gj_min..=gj_max {
            for gi in gi_min..=gi_max {
                let gp = grid_point(gi, gj, drum_radius_m, cell_size);
                let r2 = (gp - p).length_squared();
                field[gj * GRID_SIZE + gi] += splat_kernel(r2, splat_radius);
            }
        }
    }

    // Mask with the drum SDF: any grid point outside the wall (or inside a lifter bar) cannot be
    // fluid, however much splat weight landed on it. Skip points already at zero (the bulk of the
    // grid) to keep this cheap.
    let to_local = Vec2::from_angle(-drum_angle);
    for gj in 0..GRID_SIZE {
        for gi in 0..GRID_SIZE {
            let idx = gj * GRID_SIZE + gi;
            if field[idx] == 0.0 {
                continue;
            }
            let p_world = grid_point(gi, gj, drum_radius_m, cell_size);
            let p_local = rotate(p_world, to_local);
            if drum.sdf(p_local) < 0.0 {
                field[idx] = 0.0;
            }
        }
    }

    field
}

/// Rotates a 2D vector by a unit vector representing `(cos, sin)` of the rotation angle (same
/// convention as the private helper of the same name in [`crate::geometry`]).
fn rotate(v: Vec2, unit: Vec2) -> Vec2 {
    Vec2::new(v.x * unit.x - v.y * unit.y, v.x * unit.y + v.y * unit.x)
}

/// Clamps a world coordinate to the nearest in-bounds grid index along one axis.
fn grid_index(world: f32, drum_radius_m: f32, cell_size: f32) -> usize {
    let raw = ((world + drum_radius_m) / cell_size).round();
    raw.clamp(0.0, (GRID_SIZE - 1) as f32) as usize
}

fn grid_point(gi: usize, gj: usize, drum_radius_m: f32, cell_size: f32) -> Vec2 {
    Vec2::new(
        -drum_radius_m + gi as f32 * cell_size,
        -drum_radius_m + gj as f32 * cell_size,
    )
}

/// Identifies one edge of the grid (between two adjacent grid points), independent of which
/// marching-squares cell references it -- the key that lets [`extract_surface`] chain segments
/// from neighbouring cells into continuous polylines using exact (integer) matching.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum EdgeId {
    /// Between grid points `(i, j)` and `(i + 1, j)`.
    Horizontal(i32, i32),
    /// Between grid points `(i, j)` and `(i, j + 1)`.
    Vertical(i32, i32),
}

/// Linearly interpolates the threshold-crossing point along `edge`, given the field and grid
/// geometry. Always interpolates from the edge's lower-index endpoint to its higher-index
/// endpoint, so both marching-squares cells sharing this edge compute the identical point.
fn interp_edge(
    edge: EdgeId,
    field: &[f32],
    threshold: f32,
    drum_radius_m: f32,
    cell_size: f32,
) -> Vec2 {
    let ((gi0, gj0), (gi1, gj1)) = match edge {
        EdgeId::Horizontal(i, j) => ((i as usize, j as usize), (i as usize + 1, j as usize)),
        EdgeId::Vertical(i, j) => ((i as usize, j as usize), (i as usize, j as usize + 1)),
    };
    let v0 = field[gj0 * GRID_SIZE + gi0];
    let v1 = field[gj1 * GRID_SIZE + gi1];
    let t = if (v1 - v0).abs() > 1e-9 {
        ((threshold - v0) / (v1 - v0)).clamp(0.0, 1.0)
    } else {
        0.5
    };
    let p0 = grid_point(gi0, gj0, drum_radius_m, cell_size);
    let p1 = grid_point(gi1, gj1, drum_radius_m, cell_size);
    p0.lerp(p1, t)
}

/// Marching-squares segments for one cell with corner bit-mask `case` (bit0=bottom-left,
/// bit1=bottom-right, bit2=top-right, bit3=top-left; set if the corner's field value exceeds the
/// threshold), as pairs of the cell's four edges (0=bottom, 1=right, 2=top, 3=left). The two
/// ambiguous cases (5, 10) are resolved with a fixed (non-adaptive) choice; a minor, documented
/// simplification.
fn cell_edge_pairs(case: u8) -> &'static [(u8, u8)] {
    match case {
        0 | 15 => &[],
        1 | 14 => &[(3, 0)],
        2 | 13 => &[(0, 1)],
        3 | 12 => &[(3, 1)],
        4 | 11 => &[(1, 2)],
        6 | 9 => &[(0, 2)],
        7 | 8 => &[(3, 2)],
        5 => &[(3, 0), (1, 2)],
        10 => &[(0, 1), (3, 2)],
        _ => unreachable!("case is a 4-bit value, 0..=15"),
    }
}

fn cell_edge_id(cell_i: usize, cell_j: usize, edge_side: u8) -> EdgeId {
    let (i, j) = (cell_i as i32, cell_j as i32);
    match edge_side {
        0 => EdgeId::Horizontal(i, j),     // bottom
        1 => EdgeId::Vertical(i + 1, j),   // right
        2 => EdgeId::Horizontal(i, j + 1), // top
        3 => EdgeId::Vertical(i, j),       // left
        _ => unreachable!("edge_side is 0..=3"),
    }
}

/// Extracts the free-surface contour(s) of `fluid` as closed (or, rarely, open) polylines, in the
/// same world frame as `fluid.x`, against `drum` at its current `drum_angle`. See the module doc
/// comment for the algorithm, including the wall-SDF field mask and the final contour-point
/// projection that together keep every returned point on or inside the wall (and outside any
/// lifter bar).
pub fn extract_surface(fluid: &FluidParticles, drum: &Drum, drum_angle: f32) -> Vec<Polygon> {
    let drum_radius_m = drum.radius_m;
    if fluid.is_empty() || drum_radius_m <= 0.0 {
        return Vec::new();
    }
    let field = build_field(fluid, drum, drum_angle);
    let peak = field.iter().cloned().fold(0.0f32, f32::max);
    if peak <= 0.0 {
        return Vec::new();
    }
    let threshold = THRESHOLD_FRACTION * peak;
    let cell_size = 2.0 * drum_radius_m / (GRID_SIZE as f32 - 1.0);

    // Pass 1: collect every segment (as a pair of EdgeIds) from every marching-squares cell.
    let mut segments: Vec<(EdgeId, EdgeId)> = Vec::new();
    for cj in 0..GRID_SIZE - 1 {
        for ci in 0..GRID_SIZE - 1 {
            let corner =
                |di: usize, dj: usize| field[(cj + dj) * GRID_SIZE + (ci + di)] > threshold;
            let case = corner(0, 0) as u8
                | (corner(1, 0) as u8) << 1
                | (corner(1, 1) as u8) << 2
                | (corner(0, 1) as u8) << 3;
            for &(a, b) in cell_edge_pairs(case) {
                segments.push((cell_edge_id(ci, cj, a), cell_edge_id(ci, cj, b)));
            }
        }
    }
    if segments.is_empty() {
        return Vec::new();
    }

    // Pass 2: chain segments sharing an EdgeId into polylines.
    let mut by_edge: HashMap<EdgeId, Vec<usize>> = HashMap::new();
    for (idx, &(a, b)) in segments.iter().enumerate() {
        by_edge.entry(a).or_default().push(idx);
        by_edge.entry(b).or_default().push(idx);
    }
    let mut visited = vec![false; segments.len()];
    let mut polygons: Vec<Vec<EdgeId>> = Vec::new();

    for start in 0..segments.len() {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let (a0, b0) = segments[start];
        let mut chain = vec![a0, b0];
        // Extend forward from b0, then (if the chain didn't already close) backward from a0.
        extend_chain(&mut chain, &segments, &by_edge, &mut visited, true);
        extend_chain(&mut chain, &segments, &by_edge, &mut visited, false);
        // A chain that closed on itself has its first EdgeId re-appended at the end (the closing
        // segment's "other" endpoint); drop the duplicate so the closed-polygon index wrap in
        // chaikin_smooth_closed doesn't produce a spurious zero-length seam edge. Contours that
        // reach the grid boundary instead (rare: fluid always stays within the drum radius, well
        // inside the grid) stay open here and are still treated as closed by the smoothing step --
        // a known, minor simplification.
        if chain.len() > 1 && chain[0] == chain[chain.len() - 1] {
            chain.pop();
        }
        polygons.push(chain);
    }

    polygons
        .into_iter()
        .map(|chain| {
            let pts: Vec<Vec2> = chain
                .iter()
                .map(|&e| interp_edge(e, &field, threshold, drum_radius_m, cell_size))
                .collect();
            let smoothed = chaikin_smooth_closed(&pts);
            smoothed
                .into_iter()
                .map(|p| project_inside_wall(p, drum, drum_angle))
                .collect()
        })
        .collect()
}

/// Projects `p` back onto the wall (moving it along the outward normal) if it lies outside the
/// wall or inside a lifter bar, per `drum.sdf_world`. This mops up the sub-cell overshoot that
/// linear interpolation across a threshold-straddling cell (and Chaikin's lerp of such points) can
/// still leave outside the field mask applied in [`build_field`] -- up to one grid cell, since the
/// mask only guarantees grid *points* are correctly classified, not the cell interior between
/// them.
fn project_inside_wall(p: Vec2, drum: &Drum, drum_angle: f32) -> Vec2 {
    let (d, normal) = drum.sdf_world(p, drum_angle);
    if d < 0.0 {
        p - d * normal
    } else {
        p
    }
}

/// Follows unvisited segments sharing an endpoint with the current end (`forward`) or start
/// (`!forward`) of `chain`, appending (or prepending) their other endpoint, until no more
/// unvisited segment continues the chain.
fn extend_chain(
    chain: &mut Vec<EdgeId>,
    segments: &[(EdgeId, EdgeId)],
    by_edge: &HashMap<EdgeId, Vec<usize>>,
    visited: &mut [bool],
    forward: bool,
) {
    loop {
        let tip = if forward {
            *chain.last().unwrap()
        } else {
            chain[0]
        };
        let Some(candidates) = by_edge.get(&tip) else {
            return;
        };
        let Some(&next_idx) = candidates.iter().find(|&&idx| !visited[idx]) else {
            return;
        };
        visited[next_idx] = true;
        let (a, b) = segments[next_idx];
        let other = if a == tip { b } else { a };
        if forward {
            chain.push(other);
        } else {
            chain.insert(0, other);
        }
    }
}

/// One pass of Chaikin corner-cutting on a closed polygon (docs/PLAN.md ss3.5: "Chaikin smoothing
/// x1"). No-op for fewer than 3 points.
fn chaikin_smooth_closed(pts: &[Vec2]) -> Vec<Vec2> {
    if pts.len() < 3 {
        return pts.to_vec();
    }
    let n = pts.len();
    let mut out = Vec::with_capacity(n * 2);
    for i in 0..n {
        let p0 = pts[i];
        let p1 = pts[(i + 1) % n];
        out.push(p0.lerp(p1, 0.25));
        out.push(p0.lerp(p1, 0.75));
    }
    out
}

/// Flattens polygons for cross-boundary (e.g. wasm) export as `[n_polys, len_0, x, y, ..., len_1,
/// x, y, ...]` (docs/PLAN.md ss3.5).
pub fn flatten_polygons(polys: &[Polygon]) -> Vec<f32> {
    let mut out = vec![polys.len() as f32];
    for poly in polys {
        out.push(poly.len() as f32);
        for p in poly {
            out.push(p.x);
            out.push(p.y);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::SlurryParams;

    #[test]
    fn splat_kernel_is_zero_outside_radius_and_positive_inside() {
        assert!(splat_kernel(0.0, 1.0) > 0.0);
        assert!(splat_kernel(0.5, 1.0) > 0.0);
        assert_eq!(splat_kernel(1.0, 1.0), 0.0);
        assert_eq!(splat_kernel(2.0, 1.0), 0.0);
    }

    #[test]
    fn no_fluid_gives_no_surface() {
        let fluid = FluidParticles::seed_lattice(
            &SlurryParams {
                fill_fraction: 0.0,
                ..SlurryParams::default()
            },
            0.5,
            20,
            &[],
            0.0,
        );
        let drum = crate::geometry::Drum::new(0.5, 0.0, crate::params::LiftersParams::default());
        assert!(extract_surface(&fluid, &drum, 0.0).is_empty());
    }

    #[test]
    fn a_settled_puddle_has_one_closed_contour_of_plausible_area() {
        let slurry = SlurryParams {
            fill_fraction: 0.2,
            ..SlurryParams::default()
        };
        let radius_m = 0.5;
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 24, &[], 0.0);

        // Let it settle briefly (still drum) so the shape is representative, not the raw lattice.
        let drum = crate::geometry::Drum::new(
            radius_m,
            0.0,
            crate::params::LiftersParams {
                count: 0,
                ..crate::params::LiftersParams::default()
            },
        );
        for _ in 0..120 {
            fluid.step(&drum, 0.0, &slurry, 3, 1.0 / 240.0);
        }

        let polys = extract_surface(&fluid, &drum, 0.0);
        assert_eq!(
            polys.len(),
            1,
            "expected exactly one connected puddle contour, got {}",
            polys.len()
        );
        let poly = &polys[0];
        assert!(
            poly.len() >= 8,
            "contour has suspiciously few points: {}",
            poly.len()
        );

        // Shoelace area of the contour should be in the right ballpark for the requested fill
        // fraction (loose tolerance: this is a coarse marching-squares contour, not the particle
        // area itself).
        let area = shoelace_area(poly);
        let expected_area = slurry.fill_fraction * std::f32::consts::PI * radius_m * radius_m;
        let rel_err = (area - expected_area).abs() / expected_area;
        assert!(
            rel_err < 0.35,
            "contour area {area} far from expected {expected_area} (rel_err={rel_err})"
        );

        // Every contour point should lie within the drum (wall-SDF mask + contour-point
        // projection should keep this tight, not just "close").
        for &p in poly {
            assert!(
                p.length() <= radius_m + 1e-4,
                "contour point outside the drum: {p:?}"
            );
        }
    }

    #[test]
    fn spinning_drum_keeps_flung_fluid_contour_inside_the_wall() {
        let slurry = SlurryParams {
            fill_fraction: 0.2,
            ..SlurryParams::default()
        };
        let radius_m = 0.5;
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 24, &[], 0.0);

        // Fast-spinning smooth drum: fling the puddle against the wall for ~1s of sim time.
        let omega = 6.0;
        let drum = crate::geometry::Drum::new(
            radius_m,
            omega,
            crate::params::LiftersParams {
                count: 0,
                ..crate::params::LiftersParams::default()
            },
        );
        let dt = 1.0 / 240.0;
        let n_steps = (1.0 / dt) as usize;
        let mut drum_angle = 0.0f32;
        for _ in 0..n_steps {
            fluid.step(&drum, drum_angle, &slurry, 3, dt);
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
        }

        let polys = extract_surface(&fluid, &drum, drum_angle);
        assert!(!polys.is_empty(), "expected at least one contour");
        for poly in &polys {
            for &p in poly {
                assert!(
                    p.length() <= radius_m + 1e-4,
                    "contour point outside the drum after spin-up: {p:?}"
                );
            }
        }
    }

    #[test]
    fn lifters_never_appear_inside_the_extracted_contour() {
        let slurry = SlurryParams {
            fill_fraction: 0.2,
            ..SlurryParams::default()
        };
        let radius_m = 0.5;
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 24, &[], 0.0);

        let drum = crate::geometry::Drum::new(
            radius_m,
            0.0,
            crate::params::LiftersParams {
                count: 4,
                ..crate::params::LiftersParams::default()
            },
        );
        for _ in 0..120 {
            fluid.step(&drum, 0.0, &slurry, 3, 1.0 / 240.0);
        }

        let polys = extract_surface(&fluid, &drum, 0.0);
        for poly in &polys {
            for &p in poly {
                let (d, _) = drum.sdf_world(p, 0.0);
                assert!(
                    d >= -1e-4,
                    "contour point inside a lifter bar or outside the wall: {p:?} (d={d})"
                );
            }
        }
    }

    fn shoelace_area(poly: &[Vec2]) -> f32 {
        let n = poly.len();
        let mut sum = 0.0;
        for i in 0..n {
            let a = poly[i];
            let b = poly[(i + 1) % n];
            sum += a.x * b.y - b.x * a.y;
        }
        (sum * 0.5).abs()
    }

    #[test]
    fn flatten_polygons_round_trips_counts_and_points() {
        let polys = vec![
            vec![
                Vec2::new(0.0, 0.0),
                Vec2::new(1.0, 0.0),
                Vec2::new(1.0, 1.0),
            ],
            vec![Vec2::new(-1.0, -1.0), Vec2::new(-2.0, -2.0)],
        ];
        let flat = flatten_polygons(&polys);
        assert_eq!(flat[0], 2.0); // n_polys
        assert_eq!(flat[1], 3.0); // len of poly 0
        assert_eq!(&flat[2..8], &[0.0, 0.0, 1.0, 0.0, 1.0, 1.0]);
        assert_eq!(flat[8], 2.0); // len of poly 1
        assert_eq!(&flat[9..13], &[-1.0, -1.0, -2.0, -2.0]);
    }
}
