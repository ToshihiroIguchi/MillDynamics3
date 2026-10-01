//! Implicit Newtonian viscosity: pairwise-central (Monaghan/Cleary) form, solved by conjugate
//! gradient, with the moving wall as boundary pairs.
//!
//! `f_ij = -C_ij e_ij (e_ij . (v_i - v_j))`, `C_ij = 8 mu m_i m_j / (rho_i rho_j) |W'(r)| r /
//! (r^2 + eta^2)` (`8 = 2 (d + 2)` for `d = 2`; reproduces `mu laplacian(v)` on a quadratic
//! profile). Central pair forces conserve linear and angular momentum exactly and vanish for a
//! rigid rotation. Wall pairs use `rho_b = rest_density`, `m_j -> psi_b`, scaled by
//! `slurry.wall_no_slip`.

use glam::Vec2;

use super::{account_wall_pair, kinetic_energy_f64, Fluid, FluidStats, Work};
use crate::params::SlurryParams;

const CG_TOLERANCE: f32 = 1e-5;
const CG_MAX_ITERATIONS: u32 = 200;
/// `eta^2 = ETA_FACTOR * H^2` regularises `1 / r^2` for coincident particles.
const ETA_FACTOR: f32 = 0.01;

struct System<'a> {
    work: &'a Work,
    c_ff: Vec<f32>,
    e_ff: Vec<Vec2>,
    c_fb: Vec<f32>,
    e_fb: Vec<Vec2>,
    dt: f32,
}

impl System<'_> {
    /// `out = (I + dt L) v`.
    fn apply(&self, v: &[Vec2], out: &mut [Vec2]) {
        for i in 0..v.len() {
            let mut s = Vec2::ZERO;
            for k in self.work.ff.range(i) {
                let j = self.work.ff.nbrs[k] as usize;
                let e = self.e_ff[k];
                s += self.c_ff[k] * e * e.dot(v[i] - v[j]);
            }
            for k in self.work.fb.range(i) {
                let e = self.e_fb[k];
                s += self.c_fb[k] * e * e.dot(v[i]);
            }
            out[i] = v[i] + self.dt * s;
        }
    }
}

fn dot(a: &[Vec2], b: &[Vec2]) -> f32 {
    a.iter().zip(b).map(|(&x, &y)| x.dot(y)).sum()
}

impl Fluid {
    pub(super) fn viscosity_solve(
        &mut self,
        work: &Work,
        slurry: &SlurryParams,
        dt: f32,
        stats: &mut FluidStats,
    ) {
        let mu = slurry.viscosity_pa_s.max(0.0);
        if mu <= 0.0 || self.is_empty() {
            return;
        }
        let beta = slurry.wall_no_slip.clamp(0.0, 1.0);
        let n = self.len();
        let m = self.particle_mass;
        let rho0 = self.rest_density;
        let eta2 = ETA_FACTOR * self.h * self.h;

        let mut c_ff = vec![0.0f32; work.ff.nbrs.len()];
        let mut e_ff = vec![Vec2::ZERO; work.ff.nbrs.len()];
        let mut c_fb = vec![0.0f32; work.fb.nbrs.len()];
        let mut e_fb = vec![Vec2::ZERO; work.fb.nbrs.len()];
        for i in 0..n {
            let rho_i = work.rho[i].max(1e-6);
            for k in work.ff.range(i) {
                let j = work.ff.nbrs[k] as usize;
                let delta = self.x[i] - self.x[j];
                let r = delta.length();
                if r <= 1e-9 {
                    continue;
                }
                let rho_j = work.rho[j].max(1e-6);
                c_ff[k] = 8.0 * mu * m / (rho_i * rho_j) * self.kernel.dw_dr(r).abs() * r
                    / (r * r + eta2);
                e_ff[k] = delta / r;
            }
            if beta > 0.0 {
                for k in work.fb.range(i) {
                    let b = work.fb.nbrs[k] as usize;
                    let delta = self.x[i] - work.bw[b];
                    let r = delta.length();
                    if r <= 1e-9 {
                        continue;
                    }
                    c_fb[k] = beta * 8.0 * mu * self.boundary.psi[b] / (rho_i * rho0)
                        * self.kernel.dw_dr(r).abs()
                        * r
                        / (r * r + eta2);
                    e_fb[k] = delta / r;
                }
            }
        }
        let sys = System {
            work,
            c_ff,
            e_ff,
            c_fb,
            e_fb,
            dt,
        };

        // Right-hand side: incoming velocity plus the wall's drag toward its own velocity.
        let mut rhs = self.v.clone();
        for (i, r) in rhs.iter_mut().enumerate() {
            for k in work.fb.range(i) {
                let b = work.fb.nbrs[k] as usize;
                let e = sys.e_fb[k];
                *r += dt * sys.c_fb[k] * e * e.dot(work.vb[b]);
            }
        }

        let ke_before = kinetic_energy_f64(&self.v, m);
        let b_norm = dot(&rhs, &rhs).sqrt();
        if b_norm > 0.0 && b_norm.is_finite() {
            let target = CG_TOLERANCE * b_norm;
            let mut ax = vec![Vec2::ZERO; n];
            sys.apply(&self.v, &mut ax);
            let mut r: Vec<Vec2> = rhs.iter().zip(&ax).map(|(&b, &a)| b - a).collect();
            let mut p = r.clone();
            let mut ap = vec![Vec2::ZERO; n];
            let mut rs_old = dot(&r, &r);
            let mut iters = 0u32;
            if rs_old.sqrt() > target {
                loop {
                    if iters >= CG_MAX_ITERATIONS {
                        stats.viscosity_cap_hits += 1;
                        break;
                    }
                    sys.apply(&p, &mut ap);
                    let p_ap = dot(&p, &ap);
                    if p_ap <= 0.0 || !p_ap.is_finite() {
                        break;
                    }
                    let alpha = rs_old / p_ap;
                    for i in 0..n {
                        self.v[i] += alpha * p[i];
                        r[i] -= alpha * ap[i];
                    }
                    iters += 1;
                    let rs_new = dot(&r, &r);
                    if rs_new.sqrt() <= target {
                        break;
                    }
                    let beta_cg = rs_new / rs_old;
                    for i in 0..n {
                        p[i] = r[i] + beta_cg * p[i];
                    }
                    rs_old = rs_new;
                }
            }
            stats.viscosity_iterations += iters;
        }
        stats.ke_delta_viscosity_j += kinetic_energy_f64(&self.v, m) - ke_before;

        // Wall reaction from the converged velocities: `J_ib = -m dt c_ib e (e . (v_i - v_b))`.
        if beta > 0.0 {
            for i in 0..n {
                for k in work.fb.range(i) {
                    let b = work.fb.nbrs[k] as usize;
                    let e = sys.e_fb[k];
                    let j = -(m * dt * sys.c_fb[k]) * e * e.dot(self.v[i] - work.vb[b]);
                    account_wall_pair(stats, j, work.bw[b], work.vb[b], false);
                }
            }
        }
    }
}
