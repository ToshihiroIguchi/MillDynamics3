//! Flat (CSR) neighbour lists: fluid-fluid (mutual) and fluid-boundary.

use glam::Vec2;

use crate::grid::UniformGrid;

/// CSR neighbour lists: `of(i)` is `nbrs[offsets[i]..offsets[i + 1]]`.
pub struct Csr {
    pub offsets: Vec<u32>,
    pub nbrs: Vec<u32>,
}

impl Csr {
    #[inline]
    pub fn range(&self, i: usize) -> std::ops::Range<usize> {
        self.offsets[i] as usize..self.offsets[i + 1] as usize
    }

    fn from_edges(n: usize, edges: &[(u32, u32)]) -> Self {
        let mut offsets = vec![0u32; n + 1];
        for &(src, _) in edges {
            offsets[src as usize + 1] += 1;
        }
        for k in 1..offsets.len() {
            offsets[k] += offsets[k - 1];
        }
        let mut cursor = offsets.clone();
        let mut nbrs = vec![0u32; edges.len()];
        for &(src, dst) in edges {
            let pos = cursor[src as usize] as usize;
            nbrs[pos] = dst;
            cursor[src as usize] += 1;
        }
        Self { offsets, nbrs }
    }
}

/// Mutual fluid-fluid neighbours within `support` (stable counting sort, so every per-particle
/// sum below visits neighbours in a fixed, deterministic order).
pub fn build_ff(x: &[Vec2], support: f32) -> Csr {
    let grid = UniformGrid::build(x, support);
    let mut edges: Vec<(u32, u32)> = Vec::new();
    grid.for_each_candidate_pair(|i, j| {
        if (x[i as usize] - x[j as usize]).length_squared() < support * support {
            edges.push((i, j));
            edges.push((j, i));
        }
    });
    Csr::from_edges(x.len(), &edges)
}
