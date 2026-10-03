//! Variational (symmetric) discretisation of the viscous / yield stress on the staggered cut-cell
//! grid.
//!
//! The strain rate `B u + b0` (deviatoric, 2D: `xx`, `xy`, `yy` at every quadrature point) is one
//! sparse linear map of the face velocities plus wall data; the stress force is exactly its
//! transpose, `F = -B^T W g tau`, so that the implicit viscous operator `B^T W g B` and the explicit
//! stress terms of the augmented-Lagrangian iteration are adjoint to each other (a stress cannot
//! wind up in the null space of the divergence). Free surfaces need no special condition: points
//! that touch air carry no stress (the natural condition of the energy).
//!
//! Quadrature points: cell centres (`xx`, `yy` from the cell, `xy` averaged from the four corners)
//! and cell corners (`xy` from the corner, `xx`, `yy` averaged from the four cells), each with half
//! the open area of the cell. Walls enter through ghost velocities extrapolated through the
//! no-slip wall point.

use super::staggered::StaggeredMesh;
use super::viscous::{FluidGrid, Link};

/// Linear form `sum coef * x[col] + constant` in the unknown vector (`u` nodes then `v` nodes).
#[derive(Clone, Default)]
struct Lin {
    terms: Vec<(usize, f64)>,
    c: f64,
}

impl Lin {
    fn add(&mut self, other: &Lin, scale: f64) {
        self.terms
            .extend(other.terms.iter().map(|&(k, a)| (k, a * scale)));
        self.c += other.c * scale;
    }
}

pub struct ViscousOperator {
    n: usize,
    /// Rows `3 p + comp` (`xx`, `xy`, `yy`) of the points `p < n^2` (cells) and `p >= n^2`
    /// (corners); empty for invalid points.
    rows: Vec<Vec<(usize, f64)>>,
    b0: Vec<f64>,
    /// `weight / dx^2` of each point, zero for invalid points.
    pub weight: Vec<f64>,
}

/// Pairing of the stress components: `lambda : D = xx xx + 2 xy xy + yy yy`.
const PAIR: [f64; 3] = [1.0, 2.0, 1.0];

impl ViscousOperator {
    /// `cell_ok[c]`: the cell carries liquid; `fraction[c]`: its open area fraction.
    /// `active_u/v`: unknown face nodes; `wall(x, y)`: wall velocity.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        sm: &StaggeredMesh,
        active_u: &[bool],
        active_v: &[bool],
        cell_ok: &[bool],
        fraction: &[f64],
        wall: &dyn Fn(f64, f64) -> (f64, f64),
    ) -> Self {
        let n = sm.n();
        let nn = n * n;
        let inv = 1.0 / sm.dx();
        // Velocity at a node as a linear form (unknown, wall-aware ghost, or none).
        let node = |lat: usize, idx: usize| -> Option<Lin> {
            let (g, active): (&FluidGrid, &[bool]) = if lat == 0 {
                (&sm.gu, active_u)
            } else {
                (&sm.gv, active_v)
            };
            let col = |c: usize| c + lat * nn;
            if active[idx] {
                return Some(Lin {
                    terms: vec![(col(idx), 1.0)],
                    c: 0.0,
                });
            }
            if g.fluid[idx] {
                return None;
            }
            let (i, j) = (idx % n, idx / n);
            // Active neighbour along an axis whose link towards `idx` is a wall.
            for (k, (di, dj)) in [(1isize, 0isize), (-1, 0), (0, 1), (0, -1)]
                .iter()
                .enumerate()
            {
                let (ni, nj) = (i as isize + di, j as isize + dj);
                if ni < 0 || nj < 0 || ni as usize >= n || nj as usize >= n {
                    continue;
                }
                let c = ni as usize + n * nj as usize;
                if !(g.fluid[c] && active[c]) {
                    continue;
                }
                // Direction from c to idx: opposite of (di, dj).
                let dir = match k {
                    0 => 0, // idx is to the left of c (c = idx + 1)
                    1 => 1,
                    2 => 2,
                    _ => 3,
                };
                if let Link::Boundary { theta, bx, by } = g.links[c][dir] {
                    let w = wall(bx, by);
                    let wall_value = if lat == 0 { w.0 } else { w.1 };
                    let opposite = match g.links[c][dir ^ 1] {
                        Link::Fluid(o) if active[o] => Some(o),
                        _ => None,
                    };
                    return Some(match opposite {
                        Some(o) if theta < 0.25 => {
                            let f = (1.0 - theta) / (1.0 + theta);
                            Lin {
                                terms: vec![(col(o), -f)],
                                c: wall_value * (1.0 + f),
                            }
                        }
                        _ => {
                            let f = (1.0 - theta) / theta;
                            Lin {
                                terms: vec![(col(c), -f)],
                                c: wall_value * (1.0 + f),
                            }
                        }
                    });
                }
            }
            None
        };
        let diff = |a: &Lin, b: &Lin| {
            let mut out = Lin::default();
            out.add(a, 1.0);
            out.add(b, -1.0);
            out
        };
        let mut rows = vec![Vec::new(); 6 * nn];
        let mut b0 = vec![0.0; 6 * nn];
        let mut weight = vec![0.0; 2 * nn];
        // Per-cell deviatoric `xx` form and per-corner `xy` form.
        let mut xx_form: Vec<Option<Lin>> = vec![None; nn];
        let mut xy_form: Vec<Option<Lin>> = vec![None; nn];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                if cell_ok[c] && fraction[c] > 0.0 {
                    let forms = (node(0, c + 1), node(0, c), node(1, c + n), node(1, c));
                    if let (Some(ur), Some(ul), Some(vt), Some(vb)) = forms {
                        let mut f = Lin::default();
                        f.add(&diff(&ur, &ul), 0.5 * inv);
                        f.add(&diff(&vt, &vb), -0.5 * inv);
                        xx_form[c] = Some(f);
                    }
                }
                // Corner `c` (lower-left of cell c).
                let forms = (node(0, c), node(0, c - n), node(1, c), node(1, c - 1));
                if let (Some(u0), Some(u1), Some(v0), Some(v1)) = forms {
                    let mut f = Lin::default();
                    f.add(&diff(&u0, &u1), 0.5 * inv);
                    f.add(&diff(&v0, &v1), 0.5 * inv);
                    xy_form[c] = Some(f);
                }
            }
        }
        let mut put = |p: usize, comps: [Lin; 3], w: f64| {
            for (k, lin) in comps.into_iter().enumerate() {
                rows[3 * p + k] = lin.terms;
                b0[3 * p + k] = lin.c;
            }
            weight[p] = w;
        };
        for j in 2..n - 2 {
            for i in 2..n - 2 {
                let c = i + n * j;
                // Cell point: own `xx`, shear averaged over the corners that exist.
                let corners = [c, c + 1, c + n, c + n + 1];
                if let Some(xx) = &xx_form[c] {
                    let have: Vec<usize> = corners
                        .iter()
                        .copied()
                        .filter(|&k| xy_form[k].is_some())
                        .collect();
                    if !have.is_empty() {
                        let mut xy = Lin::default();
                        for &k in &have {
                            xy.add(
                                xy_form[k].as_ref().expect("checked"),
                                1.0 / have.len() as f64,
                            );
                        }
                        let mut yy = Lin::default();
                        yy.add(xx, -1.0);
                        put(c, [xx.clone(), xy, yy], 0.5 * fraction[c]);
                    }
                }
                // Corner point: own shear, normal parts averaged over the cells that exist.
                let cells = [c - n - 1, c - n, c - 1, c];
                if let Some(xy) = &xy_form[c] {
                    let have: Vec<usize> = cells
                        .iter()
                        .copied()
                        .filter(|&k| xx_form[k].is_some())
                        .collect();
                    let mut xx = Lin::default();
                    for &k in &have {
                        xx.add(
                            xx_form[k].as_ref().expect("checked"),
                            1.0 / have.len() as f64,
                        );
                    }
                    let mut yy = Lin::default();
                    yy.add(&xx, -1.0);
                    let w = 0.25 * cells.iter().map(|&k| fraction[k]).sum::<f64>();
                    if w > 0.0 {
                        put(nn + c, [xx, xy.clone(), yy], 0.5 * w);
                    }
                }
            }
        }
        Self {
            n,
            rows,
            b0,
            weight,
        }
    }

    pub fn points(&self) -> usize {
        2 * self.n * self.n
    }

    pub fn valid(&self, p: usize) -> bool {
        self.weight[p] > 0.0
    }

    /// Strain rate `[xx, xy, yy]` at every point including the wall data (zero where invalid).
    pub fn strain(&self, x: &[f64]) -> Vec<[f64; 3]> {
        (0..self.points())
            .map(|p| {
                let mut out = [0.0; 3];
                if self.valid(p) {
                    for (k, o) in out.iter_mut().enumerate() {
                        let row = 3 * p + k;
                        *o = self.b0[row]
                            + self.rows[row].iter().map(|&(c, a)| a * x[c]).sum::<f64>();
                    }
                }
                out
            })
            .collect()
    }

    /// `B^T W g z` for a stress-like field `z` (zero rows for invalid points).
    #[allow(clippy::needless_range_loop)]
    pub fn transpose(&self, z: &[[f64; 3]], out: &mut [f64]) {
        out.fill(0.0);
        for p in 0..self.points() {
            if !self.valid(p) {
                continue;
            }
            for k in 0..3 {
                let s = self.weight[p] * PAIR[k] * z[p][k];
                if s == 0.0 {
                    continue;
                }
                for &(c, a) in &self.rows[3 * p + k] {
                    out[c] += a * s;
                }
            }
        }
    }

    /// `B_lin x` (without wall data).
    fn strain_linear(&self, x: &[f64]) -> Vec<[f64; 3]> {
        (0..self.points())
            .map(|p| {
                let mut out = [0.0; 3];
                if self.valid(p) {
                    for (k, o) in out.iter_mut().enumerate() {
                        *o = self.rows[3 * p + k]
                            .iter()
                            .map(|&(c, a)| a * x[c])
                            .sum::<f64>();
                    }
                }
                out
            })
            .collect()
    }

    /// `B^T W g B_lin x`.
    pub fn apply(&self, x: &[f64], out: &mut [f64]) {
        self.transpose(&self.strain_linear(x), out);
    }

    /// `B^T W g b0r` where `b0r` is the wall-data part of the strain scaled by `scale`, added to
    /// the stress `s` first: returns `B^T W g (s + scale * b0)`.
    pub fn transpose_with_wall(&self, s: &[[f64; 3]], scale: f64, out: &mut [f64]) {
        let z: Vec<[f64; 3]> = (0..self.points())
            .map(|p| {
                let mut t = s[p];
                if self.valid(p) {
                    for (k, tk) in t.iter_mut().enumerate() {
                        *tk += scale * self.b0[3 * p + k];
                    }
                }
                t
            })
            .collect();
        self.transpose(&z, out);
    }
}

/// Preconditioned conjugate gradients for `A x = b` on the entries with `mask`; `a` applies the
/// operator, `m` the preconditioner. Returns the iteration count and the final relative residual.
pub fn pcg(
    a: &dyn Fn(&[f64], &mut [f64]),
    m: &dyn Fn(&[f64], &mut [f64]),
    mask: &[bool],
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_iter: usize,
) -> (usize, f64) {
    let len = b.len();
    let dot =
        |p: &[f64], q: &[f64]| -> f64 { (0..len).filter(|&i| mask[i]).map(|i| p[i] * q[i]).sum() };
    let mut ax = vec![0.0; len];
    a(x, &mut ax);
    let mut r: Vec<f64> = (0..len)
        .map(|i| if mask[i] { b[i] - ax[i] } else { 0.0 })
        .collect();
    let bnorm = dot(b, b).sqrt().max(1e-300);
    let mut z = vec![0.0; len];
    m(&r, &mut z);
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut res = dot(&r, &r).sqrt() / bnorm;
    let mut it = 0;
    while it < max_iter && res > tol {
        a(&p, &mut ax);
        let alpha = rz / dot(&p, &ax).max(1e-300);
        for i in 0..len {
            if mask[i] {
                x[i] += alpha * p[i];
                r[i] -= alpha * ax[i];
            }
        }
        m(&r, &mut z);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz.max(1e-300);
        rz = rz_new;
        for i in 0..len {
            p[i] = if mask[i] { z[i] + beta * p[i] } else { 0.0 };
        }
        res = dot(&r, &r).sqrt() / bnorm;
        it += 1;
    }
    (it, res)
}
