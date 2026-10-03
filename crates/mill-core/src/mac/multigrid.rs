//! Multigrid-preconditioned CG for the weighted Poisson operator
//! `(A p)_c = sum_faces w_f (p_c - p_nb)` on a cut-cell grid (pure Neumann: singular, the mean is
//! projected out).
//!
//! Coarse levels aggregate 2x2 cells with piecewise-constant transfer; summing the fine face
//! weights across each coarse face makes every coarse operator the exact Galerkin product, so the
//! V-cycle is a symmetric positive semi-definite preconditioner for CG.

use super::embedded::CircleDomain;

const SMOOTH_OMEGA: f64 = 0.8;
const PRE_POST_SWEEPS: usize = 2;
const COARSEST_SWEEPS: usize = 60;
const MIN_DIAG: f64 = 1e-12;
/// Over-correction of the aggregated coarse-grid correction (standard for piecewise-constant
/// transfer; symmetric, so the V-cycle stays a valid CG preconditioner).
const COARSE_SCALE: f64 = 1.8;

/// One grid level of the hierarchy.
#[derive(Clone, Debug)]
pub struct Level {
    pub n: usize,
    pub wx: Vec<f64>,
    pub wy: Vec<f64>,
    /// Extra diagonal per cell (mass term, Dirichlet links); zero everywhere = singular Neumann.
    pub extra: Vec<f64>,
    /// `sum w + extra` per cell; cells at or below `MIN_DIAG` are decoupled and ignored.
    pub diag: Vec<f64>,
}

impl Level {
    fn from_weights(n: usize, wx: Vec<f64>, wy: Vec<f64>, extra: Vec<f64>) -> Self {
        let mut diag = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                diag[i + n * j] = wx[i + (n + 1) * j]
                    + wx[i + 1 + (n + 1) * j]
                    + wy[i + n * j]
                    + wy[i + n * (j + 1)]
                    + extra[i + n * j];
            }
        }
        Self {
            n,
            wx,
            wy,
            extra,
            diag,
        }
    }

    fn active(&self, c: usize) -> bool {
        self.diag[c] > MIN_DIAG
    }

    /// `out = A x`.
    pub fn apply(&self, x: &[f64], out: &mut [f64]) {
        let n = self.n;
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if !self.active(c) {
                    out[c] = 0.0;
                    continue;
                }
                let xc = x[c];
                let (wl, wr) = (self.wx[i + (n + 1) * j], self.wx[i + 1 + (n + 1) * j]);
                let (wd, wu) = (self.wy[i + n * j], self.wy[i + n * (j + 1)]);
                let mut s = self.diag[c] * xc;
                if wl > 0.0 {
                    s -= wl * x[c - 1];
                }
                if wr > 0.0 {
                    s -= wr * x[c + 1];
                }
                if wd > 0.0 {
                    s -= wd * x[c - n];
                }
                if wu > 0.0 {
                    s -= wu * x[c + n];
                }
                out[c] = s;
            }
        }
    }

    fn jacobi(&self, b: &[f64], x: &mut [f64], tmp: &mut [f64], sweeps: usize) {
        for _ in 0..sweeps {
            self.apply(x, tmp);
            for c in 0..x.len() {
                if self.active(c) {
                    x[c] += SMOOTH_OMEGA * (b[c] - tmp[c]) / self.diag[c];
                }
            }
        }
    }

    fn coarsen(&self) -> Level {
        let (n, nc) = (self.n, self.n / 2);
        let mut wx = vec![0.0; (nc + 1) * nc];
        let mut wy = vec![0.0; nc * (nc + 1)];
        for jc in 0..nc {
            for ic in 0..=nc {
                let i = 2 * ic;
                wx[ic + (nc + 1) * jc] =
                    self.wx[i + (n + 1) * 2 * jc] + self.wx[i + (n + 1) * (2 * jc + 1)];
            }
        }
        for jc in 0..=nc {
            for ic in 0..nc {
                let j = 2 * jc;
                wy[ic + nc * jc] = self.wy[2 * ic + n * j] + self.wy[2 * ic + 1 + n * j];
            }
        }
        let mut extra = vec![0.0; nc * nc];
        for jc in 0..nc {
            for ic in 0..nc {
                extra[ic + nc * jc] = self.extra[2 * ic + n * 2 * jc]
                    + self.extra[2 * ic + 1 + n * 2 * jc]
                    + self.extra[2 * ic + n * (2 * jc + 1)]
                    + self.extra[2 * ic + 1 + n * (2 * jc + 1)];
            }
        }
        Level::from_weights(nc, wx, wy, extra)
    }
}

/// Convergence report of one solve.
#[derive(Clone, Copy, Debug, Default)]
pub struct SolveStats {
    pub iterations: u32,
    pub residual: f64,
}

pub struct PoissonSolver {
    levels: Vec<Level>,
    /// No extra diagonal anywhere: pure Neumann, the mean is projected out.
    singular: bool,
}

impl PoissonSolver {
    pub fn new(domain: &CircleDomain) -> Self {
        Self::from_weights(
            domain.n,
            domain.wx.clone(),
            domain.wy.clone(),
            vec![0.0; domain.n * domain.n],
        )
    }

    /// General operator `(A x)_c = extra_c x_c + sum_f w_f (x_c - x_nb)` on an `n x n` grid
    /// (`wx`: `(n+1) x n` vertical faces, `wy`: `n x (n+1)` horizontal faces).
    pub fn from_weights(n: usize, wx: Vec<f64>, wy: Vec<f64>, extra: Vec<f64>) -> Self {
        let singular = extra.iter().all(|&e| e == 0.0);
        let mut levels = vec![Level::from_weights(n, wx, wy, extra)];
        while levels.last().unwrap().n % 2 == 0 && levels.last().unwrap().n > 4 {
            let c = levels.last().unwrap().coarsen();
            levels.push(c);
        }
        Self { levels, singular }
    }

    pub fn fine(&self) -> &Level {
        &self.levels[0]
    }

    pub fn level_count(&self) -> usize {
        self.levels.len()
    }

    fn project_mean(&self, level: usize, v: &mut [f64]) {
        if !self.singular {
            return;
        }
        let lv = &self.levels[level];
        let (mut s, mut cnt) = (0.0, 0usize);
        for (c, x) in v.iter().enumerate() {
            if lv.active(c) {
                s += x;
                cnt += 1;
            }
        }
        if cnt == 0 {
            return;
        }
        let m = s / cnt as f64;
        for (c, x) in v.iter_mut().enumerate() {
            if lv.active(c) {
                *x -= m;
            }
        }
    }

    /// One V-cycle `x = B r` (starting from `x = 0`).
    fn vcycle(&self, level: usize, r: &[f64], x: &mut [f64]) {
        let lv = &self.levels[level];
        let len = r.len();
        x.fill(0.0);
        let mut tmp = vec![0.0; len];
        if level + 1 == self.levels.len() {
            lv.jacobi(r, x, &mut tmp, COARSEST_SWEEPS);
            return;
        }
        lv.jacobi(r, x, &mut tmp, PRE_POST_SWEEPS);
        lv.apply(x, &mut tmp);
        let n = lv.n;
        let nc = n / 2;
        let mut rc = vec![0.0; nc * nc];
        for jc in 0..nc {
            for ic in 0..nc {
                let mut s = 0.0;
                for (di, dj) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                    let c = 2 * ic + di + n * (2 * jc + dj);
                    s += r[c] - tmp[c];
                }
                rc[ic + nc * jc] = s;
            }
        }
        let mut ec = vec![0.0; nc * nc];
        self.vcycle(level + 1, &rc, &mut ec);
        for jc in 0..nc {
            for ic in 0..nc {
                let e = COARSE_SCALE * ec[ic + nc * jc];
                for (di, dj) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                    x[2 * ic + di + n * (2 * jc + dj)] += e;
                }
            }
        }
        lv.jacobi(r, x, &mut tmp, PRE_POST_SWEEPS);
    }

    /// Solves `A p = b` to a relative residual `tol` (zero-mean `p`); `b` need not be compatible
    /// (its mean over active cells is removed).
    pub fn solve(&self, b: &[f64], p: &mut [f64], tol: f64, max_iter: u32) -> SolveStats {
        let lv = &self.levels[0];
        let len = b.len();
        let mut r = b.to_vec();
        self.project_mean(0, &mut r);
        p.fill(0.0);
        let bnorm = norm(&r).max(1e-300);
        let mut z = vec![0.0; len];
        self.vcycle(0, &r, &mut z);
        self.project_mean(0, &mut z);
        let mut d = z.clone();
        let mut rz = dot(&r, &z);
        let mut ad = vec![0.0; len];
        let mut stats = SolveStats::default();
        for it in 0..max_iter {
            stats.iterations = it;
            stats.residual = norm(&r) / bnorm;
            if stats.residual <= tol {
                break;
            }
            lv.apply(&d, &mut ad);
            let alpha = rz / dot(&d, &ad);
            for c in 0..len {
                p[c] += alpha * d[c];
                r[c] -= alpha * ad[c];
            }
            self.vcycle(0, &r, &mut z);
            self.project_mean(0, &mut z);
            let rz_new = dot(&r, &z);
            let beta = rz_new / rz;
            rz = rz_new;
            for c in 0..len {
                d[c] = z[c] + beta * d[c];
            }
            stats.iterations = it + 1;
            stats.residual = norm(&r) / bnorm;
        }
        self.project_mean(0, p);
        stats
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
