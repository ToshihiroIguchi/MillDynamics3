//! Uniform-grid spatial hash for neighbour search.
//!
//! Shared broadphase structure used by both the ball solver ([`crate::dem`]) and the PBF fluid
//! solver ([`crate::pbf`], M3): cell coordinates map to the point indices inside that cell, cell
//! size chosen per-solver (ball diameter for DEM, kernel radius `h` for PBF), rebuilt each
//! sub-step. See docs/PLAN.md ss3.2/3.3 for how each solver uses it.
//!
//! **Counting-sort CSR, not a tree.** An earlier version stored `(i32, i32) -> Vec<u32>` in a
//! `BTreeMap` specifically for its always-sorted iteration order (`std::collections::HashMap`'s
//! default hasher is randomly seeded per process, so its iteration order is not reproducible
//! run-to-run even for byte-identical input -- silently breaking this crate's stated "every run is
//! reproducible from `Params` + seed" convention, CLAUDE.md, at the level of exact floating-point
//! summation order). That correctness property does not actually require a tree: this grid is
//! rebuilt from scratch every sub-step (never incrementally updated), so the same guarantee --
//! candidate pairs and per-point neighbour lists visited in a fixed, `(cx, cy)`-then-point-index
//! order -- is obtained instead with a single counting sort into flat `Vec<u32>` buffers
//! ([`GridLayout::Dense`]), with `O(1)` cell lookups and one allocation for the whole point set
//! instead of one small `Vec` per occupied cell. `docs/PERF.md`'s 2026-09-27 "M6 solver pass"
//! section measured this against the `BTreeMap` version this replaced.
//!
//! [`GridLayout::Sparse`] is the fallback for a pathologically large bounding box (e.g. one point
//! at a diverged-but-still-finite coordinate far from otherwise-tight data): sorting `(key, index)`
//! pairs and grouping by key gives the identical order via a stable sort instead of a dense array
//! sized to the (here, huge) bounding box. Both layouts are covered by the same determinism test
//! below (`candidate_pair_iteration_order_is_deterministic_across_rebuilds`).

use glam::Vec2;

/// Above this many addressable cells (bounding-box width * height) relative to the point count, a
/// dense linear-index array would spend far more memory and empty-cell iteration than the points
/// it actually holds warrant -- [`UniformGrid::build`] falls back to [`GridLayout::Sparse`]
/// instead. `2` (not a larger, more "generous" multiple) was chosen empirically, not
/// theoretically: the DEM ball population's own occupied-cell fraction, specifically, is well
/// below the naive "about one point per cell" assumption a larger budget implicitly relies on --
/// `cell_size` there is sized to one ball *diameter*, but a bounding box is a square while the
/// drum's fill area is a circle inside it (a `pi/4 ~= 0.785` factor on its own), and
/// `media.fill_fraction` (routinely 0.25-0.35) leaves most of the interior area itself unoccupied
/// too. `docs/PERF.md`'s 2026-09-27 "M6 solver pass" section measured `cargo bench`'s `dem_step`
/// directly against this constant: the original, more generous `8` (chosen without this empirical
/// check) cost `dem_step` a measurable few percent of the dense array's own wasted-empty-cell
/// iteration relative to `2`, for no offsetting benefit at the population sizes this crate
/// actually runs (a large point count is exactly where `2` protects against the same problem this
/// budget exists for in the first place).
const DENSE_CELLS_PER_POINT_BUDGET: u64 = 2;
/// Floor added to [`DENSE_CELLS_PER_POINT_BUDGET`]'s per-point allowance, so a small point set
/// (where `2 * n` alone would be a tiny, easily-exceeded budget) still gets the dense path for any
/// reasonably compact bounding box.
const DENSE_CELLS_BUDGET_FLOOR: u64 = 512;

/// How [`UniformGrid`] addresses its `indices` buffer's per-cell ranges. See the module doc
/// comment for why there are two layouts and when each is chosen.
enum GridLayout {
    /// Cell `(cx, cy)` maps to linear index `k = (cx - min_cx) * height + (cy - min_cy)`, and
    /// `cell_start[k]..cell_start[k + 1]` slices [`UniformGrid::indices`] for that cell.
    /// Iterating `k` ascending visits cells in exactly the `(cx, cy)` lexicographic order a
    /// `BTreeMap<(i32, i32), _>` would (first ascending `cx`, then ascending `cy` within it) --
    /// `k`'s definition is precisely that row-major order. `cell_start` has `width * height + 1`
    /// entries (the usual CSR trailing sentinel).
    Dense {
        min_cx: i32,
        min_cy: i32,
        width: i64,
        height: i64,
        cell_start: Vec<u32>,
    },
    /// Fallback ([`DENSE_CELLS_PER_POINT_BUDGET`]) for a pathologically large bounding box: only
    /// the cells that actually contain a point are stored, as `keys` (sorted ascending, unique,
    /// same `(cx, cy)` tuple order as `Dense`'s row-major `k`) with `cell_start` parallel to it
    /// (`keys.len() + 1` entries, same CSR convention). A cell's range is located by binary search
    /// instead of direct indexing.
    Sparse {
        keys: Vec<(i32, i32)>,
        cell_start: Vec<u32>,
    },
}

/// A uniform grid over a set of 2D points, used to enumerate candidate neighbour pairs (or a
/// point's neighbours) without an all-pairs O(n^2) scan.
pub struct UniformGrid {
    cell_size: f32,
    layout: GridLayout,
    /// Point indices, grouped by cell (see [`GridLayout`]); within one cell's range, indices are
    /// in ascending original-point-index order (both build paths use a stable sort/counting sort
    /// that preserves it) -- the same per-cell order the former `BTreeMap<_, Vec<u32>>`'s push-in-
    /// point-order `Vec` produced.
    indices: Vec<u32>,
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

        // Pass 1: cell key for every finite point, kept in original point-index order -- both
        // `build_dense`'s counting sort and `build_sparse`'s stable sort rely on that order to
        // reproduce the former `BTreeMap` version's per-cell insertion order.
        let mut keyed: Vec<((i32, i32), u32)> = Vec::with_capacity(points.len());
        let mut min_cx = i32::MAX;
        let mut min_cy = i32::MAX;
        let mut max_cx = i32::MIN;
        let mut max_cy = i32::MIN;
        for (i, &p) in points.iter().enumerate() {
            if !p.is_finite() {
                continue;
            }
            let key = Self::cell_key(p, cell_size);
            min_cx = min_cx.min(key.0);
            min_cy = min_cy.min(key.1);
            max_cx = max_cx.max(key.0);
            max_cy = max_cy.max(key.1);
            keyed.push((key, i as u32));
        }

        if keyed.is_empty() {
            return Self {
                cell_size,
                layout: GridLayout::Sparse {
                    keys: Vec::new(),
                    cell_start: vec![0],
                },
                indices: Vec::new(),
            };
        }

        // i64 throughout: `max_cx - min_cx` etc. as i32 arithmetic could overflow at the
        // `cell_key` saturating-cast boundary (see that method's doc comment) even though no
        // individual coordinate does.
        let width = max_cx as i64 - min_cx as i64 + 1;
        let height = max_cy as i64 - min_cy as i64 + 1;
        let total_cells = (width as u64).saturating_mul(height as u64);
        let budget = DENSE_CELLS_PER_POINT_BUDGET
            .saturating_mul(keyed.len() as u64)
            .saturating_add(DENSE_CELLS_BUDGET_FLOOR);

        if total_cells <= budget {
            Self::build_dense(keyed, cell_size, min_cx, min_cy, width, height)
        } else {
            Self::build_sparse(keyed, cell_size)
        }
    }

    /// Dense-layout path of [`UniformGrid::build`]: a two-pass counting sort (histogram, then
    /// prefix-sum, then a stable fill) into `total_cells + 1`-length `cell_start` / `keyed.len()`-
    /// length `indices` buffers -- `O(n + total_cells)`, no per-cell allocation.
    fn build_dense(
        keyed: Vec<((i32, i32), u32)>,
        cell_size: f32,
        min_cx: i32,
        min_cy: i32,
        width: i64,
        height: i64,
    ) -> Self {
        let total_cells = (width * height) as usize;
        let mut cell_start = vec![0u32; total_cells + 1];
        for &(key, _) in &keyed {
            let k = Self::dense_index(key, min_cx, min_cy, height);
            cell_start[k + 1] += 1;
        }
        for i in 1..cell_start.len() {
            cell_start[i] += cell_start[i - 1];
        }
        // Scratch write cursor, seeded from the just-computed start offsets: walking `keyed` in
        // its original (point-index-ascending) order and appending at each cell's cursor is a
        // stable counting sort, so a cell's slice of `indices` ends up point-index-ascending too
        // -- exactly the former `BTreeMap<_, Vec<u32>>`'s per-cell push order.
        let mut cursor = cell_start.clone();
        let mut indices = vec![0u32; keyed.len()];
        for &(key, idx) in &keyed {
            let k = Self::dense_index(key, min_cx, min_cy, height);
            indices[cursor[k] as usize] = idx;
            cursor[k] += 1;
        }
        Self {
            cell_size,
            layout: GridLayout::Dense {
                min_cx,
                min_cy,
                width,
                height,
                cell_start,
            },
            indices,
        }
    }

    /// Sparse-layout path of [`UniformGrid::build`] (see [`GridLayout::Sparse`]'s doc comment for
    /// when this is chosen): sort `(key, index)` pairs by key -- `Vec::sort_by_key` is a stable
    /// sort, so points sharing a cell keep their relative (point-index-ascending) order, the same
    /// per-cell order [`UniformGrid::build_dense`]'s counting sort produces -- then group into
    /// `keys`/`cell_start`.
    fn build_sparse(mut keyed: Vec<((i32, i32), u32)>, cell_size: f32) -> Self {
        keyed.sort_by_key(|&(key, _)| key);
        let mut keys = Vec::new();
        let mut cell_start = vec![0u32];
        let mut indices = Vec::with_capacity(keyed.len());
        let mut i = 0;
        while i < keyed.len() {
            let key = keyed[i].0;
            let mut j = i;
            while j < keyed.len() && keyed[j].0 == key {
                indices.push(keyed[j].1);
                j += 1;
            }
            keys.push(key);
            cell_start.push(j as u32);
            i = j;
        }
        Self {
            cell_size,
            layout: GridLayout::Sparse { keys, cell_start },
            indices,
        }
    }

    /// [`GridLayout::Dense`]'s linear cell index for `key`, given that layout's own bounding-box
    /// origin/height. Callers already know `key` falls within the dense bounding box (it was
    /// either derived from one of the points used to size it, or already range-checked by
    /// [`UniformGrid::cell_range`]).
    fn dense_index(key: (i32, i32), min_cx: i32, min_cy: i32, height: i64) -> usize {
        let dcx = key.0 as i64 - min_cx as i64;
        let dcy = key.1 as i64 - min_cy as i64;
        (dcx * height + dcy) as usize
    }

    /// The `(start, end)` range into [`UniformGrid::indices`] for cell `key`, or `None` if that
    /// cell is out of range (`Dense`, a key outside the bounding box) or simply empty (either
    /// layout) -- both cases mean "no point in that cell", matching what a `BTreeMap::get` on an
    /// absent key would have meant in the layout this replaced.
    fn cell_range(&self, key: (i32, i32)) -> Option<(usize, usize)> {
        match &self.layout {
            GridLayout::Dense {
                min_cx,
                min_cy,
                width,
                height,
                cell_start,
            } => {
                let dcx = key.0 as i64 - *min_cx as i64;
                let dcy = key.1 as i64 - *min_cy as i64;
                if dcx < 0 || dcy < 0 || dcx >= *width || dcy >= *height {
                    return None;
                }
                let k = (dcx * *height + dcy) as usize;
                let start = cell_start[k] as usize;
                let end = cell_start[k + 1] as usize;
                if start == end {
                    None
                } else {
                    Some((start, end))
                }
            }
            GridLayout::Sparse { keys, cell_start } => {
                let idx = keys.binary_search(&key).ok()?;
                Some((cell_start[idx] as usize, cell_start[idx + 1] as usize))
            }
        }
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

        // A cell's own pairs, then its forward-neighbour pairs -- shared between both layouts
        // below, given that cell's `(cx, cy)` and its already-looked-up `self_indices` slice.
        let visit_cell = |cx: i32, cy: i32, self_indices: &[u32], f: &mut F| {
            for a in 0..self_indices.len() {
                for b in (a + 1)..self_indices.len() {
                    let (i, j) = (self_indices[a], self_indices[b]);
                    f(i.min(j), i.max(j));
                }
            }
            for (dx, dy) in FORWARD_OFFSETS {
                // `build` never inserts a non-finite point, but a merely huge *finite* coordinate
                // (a badly diverged but not yet infinite simulation state) can still floor/cast to
                // a cell key at the `i32` boundary; a plain `cx + dx` there would overflow (debug
                // builds panic, release builds silently wrap), so the offset uses
                // `saturating_add`. That alone is not sufficient right at the boundary, though:
                // `i32::MAX.saturating_add(1) == i32::MAX`, so a "forward" offset can degenerate
                // to the *same* cell as `(cx, cy)` -- the explicit `==` check below skips that
                // case, since this cell's own points are already paired against each other just
                // above and must not be paired against themselves a second time as bogus
                // zero-distance self-pairs.
                let key = (cx.saturating_add(dx), cy.saturating_add(dy));
                if key == (cx, cy) {
                    continue;
                }
                let Some((ostart, oend)) = self.cell_range(key) else {
                    continue;
                };
                let other_indices = &self.indices[ostart..oend];
                for &i in self_indices {
                    for &j in other_indices {
                        f(i.min(j), i.max(j));
                    }
                }
            }
        };

        match &self.layout {
            GridLayout::Dense {
                min_cx,
                min_cy,
                width,
                height,
                cell_start,
            } => {
                let total_cells = (*width * *height) as usize;
                for k in 0..total_cells {
                    let start = cell_start[k] as usize;
                    let end = cell_start[k + 1] as usize;
                    if start == end {
                        continue;
                    }
                    let cx = *min_cx + (k as i64 / *height) as i32;
                    let cy = *min_cy + (k as i64 % *height) as i32;
                    visit_cell(cx, cy, &self.indices[start..end], &mut f);
                }
            }
            GridLayout::Sparse { keys, cell_start } => {
                for idx in 0..keys.len() {
                    let (cx, cy) = keys[idx];
                    let start = cell_start[idx] as usize;
                    let end = cell_start[idx + 1] as usize;
                    visit_cell(cx, cy, &self.indices[start..end], &mut f);
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
            let Some((start, end)) = self.cell_range(key) else {
                continue;
            };
            for &i in &self.indices[start..end] {
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
        // bit level). This grid's ordering depends only on its own layout/build logic (a
        // deterministic counting sort or stable sort, see the module doc comment), not on any
        // process-specific hasher state, so rebuilding the identical grid many times within this
        // process and comparing the exact candidate-pair sequence each time is a reasonable proxy
        // for that cross-process guarantee.
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

    #[test]
    fn sparse_fallback_matches_dense_layout_pair_order_and_set() {
        // Forces `GridLayout::Sparse` (a huge bounding box relative to the point count -- one
        // point at a diverged-but-finite coordinate far from an otherwise tight cluster) and
        // checks it against the same points run through the dense path (small `cell_size` here
        // keeps the *dense* bounding box, in cell units, under the budget for this same point
        // set, isolating the comparison to layout choice rather than also changing which cells
        // are occupied). The two layouts must agree exactly: same pairs, same order.
        let mut points: Vec<Vec2> = (0..20)
            .map(|i| {
                let t = i as f32 * 0.41;
                Vec2::new(t.sin() * 0.5, (t * 0.9).cos() * 0.5)
            })
            .collect();
        points.push(Vec2::new(1.0e6, 1.0e6)); // forces a huge bounding box => Sparse

        let sparse_grid = UniformGrid::build(&points, 0.2);
        assert!(
            matches!(sparse_grid.layout, GridLayout::Sparse { .. }),
            "test setup should have forced the sparse fallback"
        );
        let mut sparse_pairs = Vec::new();
        sparse_grid.for_each_candidate_pair(|i, j| sparse_pairs.push((i, j)));

        // Same points, minus the outlier, through the ordinary dense path.
        let dense_points = &points[..points.len() - 1];
        let dense_grid = UniformGrid::build(dense_points, 0.2);
        assert!(
            matches!(dense_grid.layout, GridLayout::Dense { .. }),
            "test setup should have kept the outlier-free set on the dense path"
        );
        let mut dense_pairs = Vec::new();
        dense_grid.for_each_candidate_pair(|i, j| dense_pairs.push((i, j)));

        // The outlier point (index `points.len() - 1`) cannot be within `cell_size` of anything,
        // so it contributes no pairs -- the sparse grid's pairs (over all 21 points) must exactly
        // match the dense grid's pairs (over the 20 non-outlier points).
        assert_eq!(sparse_pairs, dense_pairs);
    }
}
