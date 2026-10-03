//! Cell-centred level set on the MAC grid for the free surface: liquid where `psi < 0`.
//!
//! * advection `psi_t + u . grad psi = 0` with WENO5 upwind derivatives and TVD-RK3,
//! * PDE reinitialisation to a signed distance with the Russo-Smereka sub-cell fix (the zero
//!   contour is held fixed in the cells it crosses),
//! * liquid volume from a smoothed Heaviside and a global volume correction that shifts `psi`.
//!
//! The field is defined on the whole grid (also inside solids) so the contact line is free to
//! slide; the box border is far from the interface by contract, indices are clamped there.

use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub struct LevelSet {
    pub n: usize,
    pub dx: f64,
    pub half: f64,
    pub psi: Vec<f64>,
    /// Fraction of each cell that is part of the domain (1 everywhere unless set); the liquid
    /// volume and its correction only count the domain part.
    pub weight: Vec<f64>,
}

/// WENO5 derivative from five consecutive one-sided differences (Jiang-Peng / Osher-Fedkiw).
fn weno5(v1: f64, v2: f64, v3: f64, v4: f64, v5: f64) -> f64 {
    let p1 = v1 / 3.0 - 7.0 * v2 / 6.0 + 11.0 * v3 / 6.0;
    let p2 = -v2 / 6.0 + 5.0 * v3 / 6.0 + v4 / 3.0;
    let p3 = v3 / 3.0 + 5.0 * v4 / 6.0 - v5 / 6.0;
    let s1 = 13.0 / 12.0 * (v1 - 2.0 * v2 + v3).powi(2) + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3).powi(2);
    let s2 = 13.0 / 12.0 * (v2 - 2.0 * v3 + v4).powi(2) + 0.25 * (v2 - v4).powi(2);
    let s3 = 13.0 / 12.0 * (v3 - 2.0 * v4 + v5).powi(2) + 0.25 * (3.0 * v3 - 4.0 * v4 + v5).powi(2);
    let vmax = [v1, v2, v3, v4, v5]
        .iter()
        .fold(0.0f64, |m, v| m.max(v * v));
    let eps = 1e-6 * vmax + 1e-99;
    let a1 = 0.1 / (s1 + eps).powi(2);
    let a2 = 0.6 / (s2 + eps).powi(2);
    let a3 = 0.3 / (s3 + eps).powi(2);
    (a1 * p1 + a2 * p2 + a3 * p3) / (a1 + a2 + a3)
}

impl LevelSet {
    /// Samples `f(x, y)` at the cell centres (a signed distance or any function with the right sign).
    pub fn new(n: usize, half: f64, f: impl Fn(f64, f64) -> f64) -> Self {
        let dx = 2.0 * half / n as f64;
        let mut psi = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                psi[i + n * j] = f(-half + (i as f64 + 0.5) * dx, -half + (j as f64 + 0.5) * dx);
            }
        }
        Self {
            n,
            dx,
            half,
            psi,
            weight: vec![1.0; n * n],
        }
    }

    pub fn centre(&self, i: usize, j: usize) -> (f64, f64) {
        (
            -self.half + (i as f64 + 0.5) * self.dx,
            -self.half + (j as f64 + 0.5) * self.dx,
        )
    }

    fn at(&self, q: &[f64], i: isize, j: isize) -> f64 {
        let n = self.n as isize;
        q[(i.clamp(0, n - 1) + n * j.clamp(0, n - 1)) as usize]
    }

    /// Upwind WENO5 `(psi_x, psi_y)` at cell `(i, j)` for the velocity `(u, v)`.
    fn upwind_gradient(&self, q: &[f64], i: isize, j: isize, u: f64, v: f64) -> (f64, f64) {
        let inv = 1.0 / self.dx;
        let dxs = |k: isize| (self.at(q, k, j) - self.at(q, k - 1, j)) * inv;
        let dys = |k: isize| (self.at(q, i, k) - self.at(q, i, k - 1)) * inv;
        let gx = if u > 0.0 {
            weno5(dxs(i - 2), dxs(i - 1), dxs(i), dxs(i + 1), dxs(i + 2))
        } else {
            weno5(dxs(i + 3), dxs(i + 2), dxs(i + 1), dxs(i), dxs(i - 1))
        };
        let gy = if v > 0.0 {
            weno5(dys(j - 2), dys(j - 1), dys(j), dys(j + 1), dys(j + 2))
        } else {
            weno5(dys(j + 3), dys(j + 2), dys(j + 1), dys(j), dys(j - 1))
        };
        (gx, gy)
    }

    fn advection_rhs(&self, q: &[f64], uc: &[f64], vc: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut out = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if uc[c] == 0.0 && vc[c] == 0.0 {
                    continue;
                }
                let (gx, gy) = self.upwind_gradient(q, i as isize, j as isize, uc[c], vc[c]);
                out[c] = -(uc[c] * gx + vc[c] * gy);
            }
        }
        out
    }

    /// One TVD-RK3 advection step with cell-centre velocities `(uc, vc)`.
    pub fn advect(&mut self, uc: &[f64], vc: &[f64], dt: f64) {
        let q0 = self.psi.clone();
        let l0 = self.advection_rhs(&q0, uc, vc);
        let q1: Vec<f64> = q0.iter().zip(&l0).map(|(q, l)| q + dt * l).collect();
        let l1 = self.advection_rhs(&q1, uc, vc);
        let q2: Vec<f64> = (0..q0.len())
            .map(|c| 0.75 * q0[c] + 0.25 * (q1[c] + dt * l1[c]))
            .collect();
        let l2 = self.advection_rhs(&q2, uc, vc);
        for c in 0..q0.len() {
            self.psi[c] = q0[c] / 3.0 + 2.0 / 3.0 * (q2[c] + dt * l2[c]);
        }
    }

    /// Averages face velocities (`u[i + n j]` at the left face of cell `(i, j)`, `v[i + n j]` at
    /// the bottom face) to cell centres.
    pub fn centre_velocity(&self, u: &[f64], v: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = self.n;
        let mut uc = vec![0.0; n * n];
        let mut vc = vec![0.0; n * n];
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                let c = i + n * j;
                uc[c] = 0.5 * (u[c] + u[c + 1]);
                vc[c] = 0.5 * (v[c] + v[c + n]);
            }
        }
        (uc, vc)
    }

    /// PDE reinitialisation towards `|grad psi| = 1` with the interface held fixed (Russo-Smereka).
    pub fn reinitialize(&mut self, iterations: usize) {
        let n = self.n;
        let dx = self.dx;
        let dtau = 0.3 * dx;
        let psi0 = self.psi.clone();
        let sgn: Vec<f64> = psi0.iter().map(|&p| p / (p * p + dx * dx).sqrt()).collect();
        // Interface cells: sign change with a 4-neighbour. D = distance estimate to the contour.
        let mut dist = vec![f64::NAN; n * n];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                let crosses = [c - 1, c + 1, c - n, c + n]
                    .iter()
                    .any(|&k| psi0[k] * psi0[c] < 0.0);
                if crosses || psi0[c] == 0.0 {
                    let (xc, yc) = self.centre(i, j);
                    dist[c] = self.contour_distance(&psi0, xc, yc, psi0[c]);
                }
            }
        }
        for _ in 0..iterations {
            let q = self.psi.clone();
            for j in 1..n - 1 {
                for i in 1..n - 1 {
                    let c = i + n * j;
                    if !dist[c].is_nan() {
                        self.psi[c] = q[c] - dtau / dx * (sgn[c].signum() * q[c].abs() - dist[c]);
                        continue;
                    }
                    let s = sgn[c];
                    let (dxm, dxp) = ((q[c] - q[c - 1]) / dx, (q[c + 1] - q[c]) / dx);
                    let (dym, dyp) = ((q[c] - q[c - n]) / dx, (q[c + n] - q[c]) / dx);
                    let g2 = if s > 0.0 {
                        dxm.max(0.0).powi(2).max((-dxp).max(0.0).powi(2))
                            + dym.max(0.0).powi(2).max((-dyp).max(0.0).powi(2))
                    } else {
                        (-dxm).max(0.0).powi(2).max(dxp.max(0.0).powi(2))
                            + (-dym).max(0.0).powi(2).max(dyp.max(0.0).powi(2))
                    };
                    self.psi[c] = q[c] - dtau * s * (g2.sqrt() - 1.0);
                }
            }
        }
    }

    /// Tensor cubic Lagrange interpolant of `q` at `(x, y)` with its gradient.
    fn interpolate(&self, q: &[f64], x: f64, y: f64) -> (f64, f64, f64) {
        let fx = (x + self.half) / self.dx - 0.5;
        let fy = (y + self.half) / self.dx - 0.5;
        let (i0, j0) = (fx.floor(), fy.floor());
        let (sx, sy) = (fx - i0, fy - j0);
        let weights = |s: f64| {
            (
                [
                    -s * (s - 1.0) * (s - 2.0) / 6.0,
                    (s + 1.0) * (s - 1.0) * (s - 2.0) / 2.0,
                    -(s + 1.0) * s * (s - 2.0) / 2.0,
                    (s + 1.0) * s * (s - 1.0) / 6.0,
                ],
                [
                    -(3.0 * s * s - 6.0 * s + 2.0) / 6.0,
                    (3.0 * s * s - 4.0 * s - 1.0) / 2.0,
                    -(3.0 * s * s - 2.0 * s - 2.0) / 2.0,
                    (3.0 * s * s - 1.0) / 6.0,
                ],
            )
        };
        let (wx, dwx) = weights(sx);
        let (wy, dwy) = weights(sy);
        let (mut v, mut gx, mut gy) = (0.0, 0.0, 0.0);
        for (b, (wyb, dwyb)) in wy.iter().zip(&dwy).enumerate() {
            for (a, (wxa, dwxa)) in wx.iter().zip(&dwx).enumerate() {
                let val = self.at(
                    q,
                    i0 as isize + a as isize - 1,
                    j0 as isize + b as isize - 1,
                );
                v += wxa * wyb * val;
                gx += dwxa * wyb * val;
                gy += wxa * dwyb * val;
            }
        }
        (v, gx / self.dx, gy / self.dx)
    }

    /// Signed distance from `(xc, yc)` to the zero contour of the cubic interpolant of `q`
    /// (Newton projection onto the contour); sign of `sign_of`.
    fn contour_distance(&self, q: &[f64], xc: f64, yc: f64, sign_of: f64) -> f64 {
        let (mut x, mut y) = (xc, yc);
        for _ in 0..8 {
            let (v, gx, gy) = self.interpolate(q, x, y);
            let g2 = gx * gx + gy * gy;
            if g2 < 1e-24 {
                break;
            }
            let (ddx, ddy) = (v * gx / g2, v * gy / g2);
            x -= ddx;
            y -= ddy;
            if ddx.abs() + ddy.abs() < 1e-10 * self.dx {
                break;
            }
        }
        let d = ((x - xc).powi(2) + (y - yc).powi(2)).sqrt();
        if sign_of < 0.0 {
            -d
        } else {
            d
        }
    }

    fn smooth_heaviside(psi: f64, eps: f64) -> f64 {
        if psi <= -eps {
            0.0
        } else if psi >= eps {
            1.0
        } else {
            0.5 * (1.0 + psi / eps + (PI * psi / eps).sin() / PI)
        }
    }

    /// Liquid area `sum (1 - H(psi)) dx^2` (smoothed Heaviside of half-width `1.5 dx`).
    pub fn volume(&self) -> f64 {
        let eps = 1.5 * self.dx;
        self.psi
            .iter()
            .zip(&self.weight)
            .map(|(&p, &w)| w * (1.0 - Self::smooth_heaviside(p, eps)))
            .sum::<f64>()
            * self.dx
            * self.dx
    }

    /// Interface length `sum delta(psi) |grad psi| dx^2` (for the volume correction).
    fn perimeter(&self) -> f64 {
        let n = self.n;
        let eps = 1.5 * self.dx;
        let mut s = 0.0;
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                let p = self.psi[c];
                if p.abs() < eps && self.weight[c] > 0.0 {
                    let gx = (self.psi[c + 1] - self.psi[c - 1]) / (2.0 * self.dx);
                    let gy = (self.psi[c + n] - self.psi[c - n]) / (2.0 * self.dx);
                    let delta = 0.5 / eps * (1.0 + (PI * p / eps).cos());
                    s += self.weight[c] * delta * (gx * gx + gy * gy).sqrt();
                }
            }
        }
        s * self.dx * self.dx
    }

    /// Shifts `psi` uniformly so that the liquid area equals `target`.
    pub fn correct_volume(&mut self, target: f64) {
        let per = self.perimeter();
        if per > 0.0 {
            let shift = (target - self.volume()) / per;
            for p in &mut self.psi {
                *p -= shift;
            }
        }
    }
}
