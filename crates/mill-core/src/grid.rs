//! Uniform-grid spatial hash for neighbour search.
//!
//! Shared broadphase structure used by both the ball solver ([`crate::dem`]) and the PBF fluid
//! solver ([`crate::pbf`], M3): a hash map from integer cell coordinates to the point indices
//! inside that cell, cell size chosen per-solver (ball diameter for DEM, kernel radius `h` for
//! PBF), rebuilt each sub-step. See docs/PLAN.md ss3.2/3.3 for how each solver uses it.

use std::collections::HashMap;

use glam::Vec2;

/// A uniform grid over a set of 2D points, used to enumerate candidate neighbour pairs (or a
/// point's neighbours) without an all-pairs O(n^2) scan.
pub struct UniformGrid {
    cell_size: f32,
    cells: HashMap<(i32, i32), Vec<u32>>,
}

impl UniformGrid {
    /// Builds a grid bucketing every point in `points` by its cell of side `cell_size`.
    ///
    /// `cell_size` must be at least as large as the interaction range being queried (e.g. the sum
    /// of two particle radii for contact detection, or the kernel radius for SPH-style
    /// neighbourhoods): [`UniformGrid::for_each_candidate_pair`] only considers same-and-adjacent
    /// cells, so any true neighbour pair farther apart than `cell_size` would be missed.
    pub fn build(points: &[Vec2], cell_size: f32) -> Self {
        assert!(
            cell_size.is_finite() && cell_size > 0.0,
            "cell_size must be positive"
        );
        let mut cells: HashMap<(i32, i32), Vec<u32>> = HashMap::new();
        for (i, &p) in points.iter().enumerate() {
            cells
                .entry(Self::cell_key(p, cell_size))
                .or_default()
                .push(i as u32);
        }
        Self { cell_size, cells }
    }

    fn cell_key(p: Vec2, cell_size: f32) -> (i32, i32) {
        (
            (p.x / cell_size).floor() as i32,
            (p.y / cell_size).floor() as i32,
        )
    }

    /// Calls `f(i, j)` (with `i < j`) exactly once for every pair of point indices whose cells are
    /// the same or adjacent (a 3x3 neighbourhood) — i.e. every pair that could plausibly be within
    /// `cell_size` of each other. Callers must still check the actual distance; this only narrows
    /// an O(n^2) scan down to spatially-nearby candidates.
    pub fn for_each_candidate_pair<F: FnMut(u32, u32)>(&self, mut f: F) {
        // Only a "forward" half of the 8-neighbourhood plus the cell itself is visited from each
        // cell, so every unordered pair of distinct adjacent cells is covered exactly once
        // (the mirror-image offsets are what the *other* cell would use to reach this one).
        const FORWARD_OFFSETS: [(i32, i32); 4] = [(1, 0), (1, 1), (0, 1), (-1, 1)];

        for (&(cx, cy), indices) in &self.cells {
            // Pairs within this cell.
            for a in 0..indices.len() {
                for b in (a + 1)..indices.len() {
                    let (i, j) = (indices[a], indices[b]);
                    f(i.min(j), i.max(j));
                }
            }
            // Pairs against each forward-neighbouring cell.
            for (dx, dy) in FORWARD_OFFSETS {
                let Some(other) = self.cells.get(&(cx + dx, cy + dy)) else {
                    continue;
                };
                for &i in indices {
                    for &j in other {
                        f(i.min(j), i.max(j));
                    }
                }
            }
        }
    }

    /// Calls `f(j)` for every point index in the same or an adjacent cell to `p` (a 3x3
    /// neighbourhood around `p`'s own cell). Used to query neighbours of a point that is not
    /// itself necessarily one of the grid's indexed points (e.g. a ball centre when broad-phasing
    /// against fluid particles).
    pub fn for_each_near<F: FnMut(u32)>(&self, p: Vec2, mut f: F) {
        let (cx, cy) = Self::cell_key(p, self.cell_size);
        for dy in -1..=1 {
            for dx in -1..=1 {
                let Some(indices) = self.cells.get(&(cx + dx, cy + dy)) else {
                    continue;
                };
                for &i in indices {
                    f(i);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_all_pairs_within_cell_size_no_duplicates_no_self_pairs() {
        // A small cluster where every point is mutually within cell_size of every other, plus one
        // far-away outlier that must not be paired with anything.
        let points = vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(0.05, 0.0),
            Vec2::new(0.0, 0.05),
            Vec2::new(100.0, 100.0), // outlier, far from the cluster
        ];
        let grid = UniformGrid::build(&points, 0.2);

        let mut pairs = Vec::new();
        grid.for_each_candidate_pair(|i, j| pairs.push((i, j)));

        // No self-pairs, no duplicates, and always reported as (min, max).
        for &(i, j) in &pairs {
            assert!(i < j);
        }
        let mut deduped = pairs.clone();
        deduped.sort_unstable();
        deduped.dedup();
        assert_eq!(
            pairs.len(),
            deduped.len(),
            "pair reported more than once: {pairs:?}"
        );

        // The three clustered points (0,1,2) must all appear pairwise; the outlier (3) must not.
        for expected in [(0, 1), (0, 2), (1, 2)] {
            assert!(
                pairs.contains(&expected),
                "missing expected pair {expected:?} in {pairs:?}"
            );
        }
        assert!(!pairs.iter().any(|&(i, j)| i == 3 || j == 3));
    }

    #[test]
    fn candidate_pairs_are_a_superset_of_true_neighbours() {
        // Randomized-ish grid of points; every pair genuinely within cell_size must be reported
        // as a candidate (the grid is allowed false positives, never false negatives).
        let cell_size = 1.0;
        let points: Vec<Vec2> = (0..40)
            .map(|i| {
                let t = i as f32 * 0.37;
                Vec2::new(t.sin() * 3.0, (t * 1.3).cos() * 3.0)
            })
            .collect();
        let grid = UniformGrid::build(&points, cell_size);

        let mut candidates = std::collections::HashSet::new();
        grid.for_each_candidate_pair(|i, j| {
            candidates.insert((i, j));
        });

        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                if points[i].distance(points[j]) <= cell_size {
                    assert!(
                        candidates.contains(&(i as u32, j as u32)),
                        "true neighbour pair ({i},{j}) missing from candidates"
                    );
                }
            }
        }
    }

    #[test]
    fn for_each_near_finds_points_in_adjacent_cells() {
        let points = vec![
            Vec2::new(0.9, 0.0),
            Vec2::new(-0.9, 0.0),
            Vec2::new(50.0, 50.0),
        ];
        let grid = UniformGrid::build(&points, 1.0);

        let mut found = Vec::new();
        grid.for_each_near(Vec2::ZERO, |i| found.push(i));
        found.sort_unstable();

        assert_eq!(found, vec![0, 1]);
    }

    #[test]
    fn empty_grid_has_no_pairs_and_no_neighbours() {
        let grid = UniformGrid::build(&[], 1.0);
        let mut count = 0;
        grid.for_each_candidate_pair(|_, _| count += 1);
        assert_eq!(count, 0);

        let mut near_count = 0;
        grid.for_each_near(Vec2::ZERO, |_| near_count += 1);
        assert_eq!(near_count, 0);
    }
}
