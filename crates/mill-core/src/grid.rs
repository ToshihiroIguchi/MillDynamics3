//! Uniform-grid spatial hash for neighbour search.
//!
//! Shared broadphase structure used by both the ball solver ([`crate::dem`]) and the PBF fluid
//! solver ([`crate::pbf`], M3): cell coordinates map to the point indices inside that cell, cell
//! size chosen per-solver (ball diameter for DEM, kernel radius `h` for PBF), rebuilt each
//! sub-step. See docs/PLAN.md ss3.2/3.3 for how each solver uses it.
//!
//! **`BTreeMap`, not `HashMap`.** `std::collections::HashMap`'s default hasher is randomly seeded
//! per process, so its iteration order (used by [`UniformGrid::for_each_candidate_pair`]) is not
//! reproducible run-to-run even for byte-identical input -- silently breaking this crate's stated
//! "every run is reproducible from `Params` + `seed`" convention (CLAUDE.md) at the level of exact
//! floating-point summation order (contact/kernel contributions get accumulated in a different
//! sequence each run, which f32 addition is not associative under). `BTreeMap` costs `O(log n)`
//! instead of `O(1)` per cell lookup, negligible next to the per-cell candidate-pair work this
//! grid exists to bound in the first place, in exchange for a fully deterministic iteration order.

use std::collections::BTreeMap;

use glam::Vec2;

/// A uniform grid over a set of 2D points, used to enumerate candidate neighbour pairs (or a
/// point's neighbours) without an all-pairs O(n^2) scan.
pub struct UniformGrid {
    cell_size: f32,
    cells: BTreeMap<(i32, i32), Vec<u32>>,
}

impl UniformGrid {
    /// Builds a grid bucketing every point in `points` by its cell of side `cell_size`. A
    /// non-finite point (NaN or +/-infinity in either coordinate) is skipped entirely -- it is
    /// never inserted into any cell, so it can never be reported as a candidate neighbour of
    /// anything (see [`UniformGrid::cell_key`]'s doc comment for why simply computing *some* cell
    /// key for it, instead of excluding it, is not safe). Its index is otherwise left unused,
    /// which every caller already treats as "no candidates" for that point.
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
        let mut cells: BTreeMap<(i32, i32), Vec<u32>> = BTreeMap::new();
        for (i, &p) in points.iter().enumerate() {
            if !p.is_finite() {
                continue;
            }
            cells
                .entry(Self::cell_key(p, cell_size))
                .or_default()
                .push(i as u32);
        }
        Self { cell_size, cells }
    }

    /// Integer cell coordinates of `p`. A non-finite `p` is never passed here by
    /// [`UniformGrid::build`] (it skips such points) or [`UniformGrid::for_each_near`] (it
    /// returns early for one) -- but if it ever were, `as i32`'s saturating float-to-int cast
    /// would silently fold `NaN` to `(0, 0)` (a real, populated cell, not a harmless out-of-range
    /// one) and +/-infinity to `i32::{MAX, MIN}` (which, added to a same-sign offset one cell
    /// over, would *overflow* rather than simply miss). Neither failure mode is a safe "no
    /// neighbours" default, which is why both callers guard against non-finite input themselves
    /// rather than relying on this function to degrade gracefully.
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
            // Pairs against each forward-neighbouring cell. `build` never inserts a non-finite
            // point, but a merely huge *finite* coordinate (a badly diverged but not yet
            // infinite simulation state) can still floor/cast to a cell key at the `i32`
            // boundary; a plain `cx + dx` there would overflow (debug builds panic, release
            // builds silently wrap), so the offset uses `saturating_add`. That alone is not
            // sufficient right at the boundary, though: `i32::MAX.saturating_add(1) ==
            // i32::MAX`, so a "forward" offset can degenerate to the *same* cell as `(cx, cy)`
            // -- the explicit `!=` check below skips that case, since this cell's own points are
            // already paired against each other just above and must not be paired against
            // themselves a second time as bogus zero-distance self-pairs.
            for (dx, dy) in FORWARD_OFFSETS {
                let key = (cx.saturating_add(dx), cy.saturating_add(dy));
                if key == (cx, cy) {
                    continue;
                }
                let Some(other) = self.cells.get(&key) else {
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
    /// against fluid particles). A non-finite `p` is treated as having no neighbours at all
    /// (`f` is never called) -- see [`UniformGrid::cell_key`]'s doc comment for why a non-finite
    /// coordinate cannot be allowed to produce *some* cell key here (`NaN` would alias a real,
    /// populated cell at the origin, not miss cleanly).
    pub fn for_each_near<F: FnMut(u32)>(&self, p: Vec2, mut f: F) {
        if !p.is_finite() {
            return;
        }
        let (cx, cy) = Self::cell_key(p, self.cell_size);
        // Collect the up-to-9 candidate cell keys before querying, deduplicating with a linear
        // scan (cheap at this size): right at the `i32` boundary (see `cell_key`'s doc comment),
        // `saturating_add` can fold two or three distinct offsets onto the same key, which would
        // otherwise report that cell's points more than once to `f`.
        let mut keys: [(i32, i32); 9] = [(0, 0); 9];
        let mut n_keys = 0;
        for dy in -1..=1 {
            for dx in -1..=1 {
                let key = (cx.saturating_add(dx), cy.saturating_add(dy));
                if !keys[..n_keys].contains(&key) {
                    keys[n_keys] = key;
                    n_keys += 1;
                }
            }
        }
        for &key in &keys[..n_keys] {
            let Some(indices) = self.cells.get(&key) else {
                continue;
            };
            for &i in indices {
                f(i);
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

    #[test]
    fn for_each_near_does_not_overflow_on_a_non_finite_query_point() {
        // Regression: `cell_key`'s `as i32` cast saturates a non-finite coordinate to
        // `i32::{MAX, MIN}`; a plain `cx + dx` on that value then overflows (debug builds panic).
        // Every non-finite value that can appear in this crate's f32 state (see e.g.
        // `crate::pbf::FluidParticles::step_coupled`'s step-0 sanitize comment) must be safe here.
        let points = vec![Vec2::new(0.9, 0.0), Vec2::new(-0.9, 0.0)];
        let grid = UniformGrid::build(&points, 1.0);
        for p in [
            Vec2::splat(f32::INFINITY),
            Vec2::splat(f32::NEG_INFINITY),
            Vec2::new(f32::INFINITY, 0.0),
            Vec2::splat(f32::NAN),
        ] {
            let mut count = 0;
            grid.for_each_near(p, |_| count += 1);
            assert_eq!(
                count, 0,
                "expected no neighbours for non-finite point {p:?}"
            );
        }
    }

    #[test]
    fn build_and_candidate_pairs_do_not_overflow_when_a_point_is_non_finite() {
        // Same regression as `for_each_near_does_not_overflow_on_a_non_finite_query_point`, but
        // for a non-finite point stored *in* the grid (e.g. an already-diverged ball position)
        // rather than only used as a query -- `for_each_candidate_pair` walks every stored cell
        // key, including one saturated to `i32::MAX`, so it needs the same guard.
        let points = vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(0.05, 0.0),
            Vec2::splat(f32::INFINITY),
        ];
        let grid = UniformGrid::build(&points, 0.2);
        let mut pairs = Vec::new();
        grid.for_each_candidate_pair(|i, j| pairs.push((i, j)));
        // The two finite points are still found as neighbours; the non-finite one pairs with
        // nothing (its saturated cell has no real neighbours).
        assert!(pairs.contains(&(0, 1)));
        assert!(!pairs.iter().any(|&(i, j)| i == 2 || j == 2));
    }

    #[test]
    fn candidate_pair_iteration_order_is_deterministic_across_rebuilds() {
        // Regression: this crate is documented (CLAUDE.md) as "every run is reproducible from
        // Params + seed", which requires not just the same *set* of candidate pairs but the same
        // *order* every time (contact/kernel contributions are accumulated with f32 addition,
        // which is not associative, so a different summation order can shift the result at the
        // bit level). A plain `HashMap`'s default hasher is randomly seeded per process, so its
        // iteration order is not reproducible run-to-run for byte-identical input -- this is why
        // `UniformGrid` uses a `BTreeMap` (see the module doc comment). Rebuilding the identical
        // grid many times within *this* process and comparing the exact candidate-pair sequence
        // each time is a reasonable proxy for that cross-process guarantee, since the ordering
        // only depends on the map's key comparison, not on any process-specific hasher state.
        let points: Vec<Vec2> = (0..60)
            .map(|i| {
                let t = i as f32 * 0.53;
                Vec2::new(t.sin() * 2.0, (t * 0.7).cos() * 2.0)
            })
            .collect();

        let first_run: Vec<(u32, u32)> = {
            let grid = UniformGrid::build(&points, 0.3);
            let mut pairs = Vec::new();
            grid.for_each_candidate_pair(|i, j| pairs.push((i, j)));
            pairs
        };
        for _ in 0..5 {
            let grid = UniformGrid::build(&points, 0.3);
            let mut pairs = Vec::new();
            grid.for_each_candidate_pair(|i, j| pairs.push((i, j)));
            assert_eq!(
                pairs, first_run,
                "candidate-pair order changed across an identical rebuild"
            );
        }
    }
}
