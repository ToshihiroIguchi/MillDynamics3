//! Implicit Newtonian viscosity: pairwise-central (Monaghan/Cleary) form, solved by
//! preconditioned conjugate gradient, with the solids (drum wall, balls) as smooth coupling
//! tensors and the balls' own velocities as unknowns.
//!
//! Fluid pairs: `f_ij = -C_ij e (e . (v_i - v_j))`, `C_ij = 8 mu m_i m_j / (rho_i rho_j) |W'(r)| r /
//! (r^2 + eta^2)` (`8 = 2 (d + 2)` for `d = 2`; reproduces `mu laplacian(v)` on a quadratic
//! profile, corrected by the lattice factor in [`Fluid::viscosity_norm`]). Central pair forces
//! conserve linear and angular momentum exactly and vanish for a rigid rotation. Solids: the same
//! pair form integrated over the solid region, i.e. the force per unit mass `-(8 mu / rho_i) T
//! (v_i - v_S)` with `T` from [`super::solid`], scaled by `slurry.wall_no_slip` /
//! `slurry.ball_no_slip`.
//!
//! **Balls are unknowns of the same system.** A ball's viscous relaxation rate (`~ 4 pi mu / m`,
//! ~8e3 /s at the defaults) times the step is ~8, so an explicit (staggered) coupling overshoots;
//! instead `(M + dt K) u' = M u` is solved for `u = (v_i, V_B, omega_B)` with `M = diag(m, m_B,
//! I_B)` and `K` the (symmetric positive semi-definite) Hessian of the dissipation `1/2 sum m k
//! Delta^T T Delta`, `Delta = v_i - V_B - omega_B perp(x_i - x_B)`. Linear and angular momentum of
//! fluid + balls is conserved exactly by construction. With `bodies_prescribed` the balls keep
//! their velocity (they act like the wall), which the Taylor-Couette verification uses.

use glam::Vec2;

use super::solid::{SolidSample, ETA_FACTOR, WALL};
use super::{account_solid, kinetic_energy_f64, Fluid, FluidStats, Work};
use crate::params::SlurryParams;

const CG_TOLERANCE: f32 = 1e-5;
const CG_MAX_ITERATIONS: u32 = 200;

/// Unknowns: fluid velocities, ball velocities and ball spins.
#[derive(Clone)]
struct State {
    v: Vec<Vec2>,
    bv: Vec<Vec2>,
    bw: Vec<f32>,
}

impl State {
    fn zeros(n: usize, nb: usize) -> Self {
        Self {
            v: vec![Vec2::ZERO; n],
            bv: vec![Vec2::ZERO; nb],
            bw: vec![0.0; nb],
        }
    }

    fn dot(&self, o: &State) -> f32 {
        self.v.iter().zip(&o.v).map(|(a, b)| a.dot(*b)).sum::<f32>()
            + self
                .bv
                .iter()
                .zip(&o.bv)
                .map(|(a, b)| a.dot(*b))
                .sum::<f32>()
            + self.bw.iter().zip(&o.bw).map(|(a, b)| a * b).sum::<f32>()
    }
}

#[inline]
fn perp(l: Vec2) -> Vec2 {
    Vec2::new(-l.y, l.x)
}

struct System<'a> {
    work: &'a Work,
    m: f32,
    m_ball: f32,
    i_ball: f32,
    c_ff: Vec<f32>,
    e_ff: Vec<Vec2>,
    /// Per solid sample: `8 mu beta / rho_i`.
    k_solid: Vec<f32>,
    /// Balls are unknowns (false: their velocity is fixed like the wall's).
    dof: bool,
    dt: f32,
}

impl System<'_> {
    /// Velocity of the solid sample's body at the particle *as a function of the unknowns*: the
    /// ball's `V + omega perp(l)` when balls are unknowns, zero otherwise (fixed velocities live
    /// in the right-hand side so the operator stays linear).
    #[inline]
    fn solid_velocity(&self, s: &SolidSample, u: &State) -> Vec2 {
        if self.dof && s.owner != WALL {
            let o = s.owner as usize;
            u.bv[o] + u.bw[o] * perp(s.lever)
        } else {
            Vec2::ZERO
        }
    }

    /// `out = (M + dt K) u`.
    fn apply(&self, u: &State, out: &mut State) {
        let dt = self.dt;
        let (m, m_ball, i_ball) = (self.m, self.m_ball, self.i_ball);
        for (o, (v, w)) in out.bv.iter_mut().zip(out.bw.iter_mut()).enumerate() {
            *v = m_ball * u.bv[o];
            *w = i_ball * u.bw[o];
        }
        for i in 0..u.v.len() {
            let mut s = Vec2::ZERO;
            for k in self.work.ff.range(i) {
                let j = self.work.ff.nbrs[k] as usize;
                let e = self.e_ff[k];
                s += self.c_ff[k] * e * e.dot(u.v[i] - u.v[j]);
            }
            for idx in self.work.solids_of(i) {
                let sol = &self.work.solids[idx];
                let force = self.k_solid[idx] * sol.apply_t(u.v[i] - self.solid_velocity(sol, u));
                s += force;
                if self.dof && sol.owner != WALL {
                    let o = sol.owner as usize;
                    out.bv[o] -= dt * m * force;
                    out.bw[o] -= dt * m * perp(sol.lever).dot(force);
                }
            }
            out.v[i] = m * (u.v[i] + dt * s);
        }
    }
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
        let beta_wall = slurry.wall_no_slip.clamp(0.0, 1.0);
        let beta_ball = slurry.ball_no_slip.clamp(0.0, 1.0);
        let n = self.len();
        let nb = self.bodies.x.len();
        let m = self.particle_mass;
        let eta2 = ETA_FACTOR * self.h * self.h;
        let dof = nb > 0 && !self.bodies_prescribed && self.ball_mass > 0.0;
        let nb_unknown = if dof { nb } else { 0 };

        let mut c_ff = vec![0.0f32; work.ff.nbrs.len()];
        let mut e_ff = vec![Vec2::ZERO; work.ff.nbrs.len()];
        let mut k_solid = vec![0.0f32; work.solids.len()];
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
                c_ff[k] = self.viscosity_norm * 8.0 * mu * m / (rho_i * rho_j)
                    * self.kernel.dw_dr(r).abs()
                    * r
                    / (r * r + eta2);
                e_ff[k] = delta / r;
            }
            for idx in work.solids_of(i) {
                let beta = if work.solids[idx].owner == WALL {
                    beta_wall
                } else {
                    beta_ball
                };
                k_solid[idx] = 8.0 * mu * beta / rho_i;
            }
        }
        let sys = System {
            work,
            m,
            m_ball: self.ball_mass,
            i_ball: self.ball_inertia,
            c_ff,
            e_ff,
            k_solid,
            dof,
            dt,
        };

        // Initial state and right-hand side `M u0 + dt m k T v_S` for solids whose velocity is
        // not an unknown (the wall; balls when prescribed).
        let mut u = State::zeros(n, nb_unknown);
        u.v.copy_from_slice(&self.v);
        if dof {
            u.bv.copy_from_slice(&self.bodies.v);
            u.bw.copy_from_slice(&self.bodies.omega);
        }
        let mut rhs = State::zeros(n, nb_unknown);
        for i in 0..n {
            let mut r = self.v[i];
            for idx in work.solids_of(i) {
                let s = &work.solids[idx];
                if !(dof && s.owner != WALL) {
                    r += dt * sys.k_solid[idx] * s.apply_t(s.vel);
                }
            }
            rhs.v[i] = m * r;
        }
        for o in 0..nb_unknown {
            rhs.bv[o] = self.ball_mass * self.bodies.v[o];
            rhs.bw[o] = self.ball_inertia * self.bodies.omega[o];
        }

        // Jacobi preconditioner: diagonal of `M + dt K` per component.
        let mut dg = State::zeros(n, nb_unknown);
        for i in 0..n {
            let mut d = Vec2::splat(m);
            for k in work.ff.range(i) {
                let e = sys.e_ff[k];
                d += dt * m * sys.c_ff[k] * Vec2::new(e.x * e.x, e.y * e.y);
            }
            for idx in work.solids_of(i) {
                let s = &work.solids[idx];
                let kk = dt * m * sys.k_solid[idx];
                let txx = s.b + (s.a - s.b) * s.n.x * s.n.x;
                let tyy = s.b + (s.a - s.b) * s.n.y * s.n.y;
                d += kk * Vec2::new(txx, tyy);
                if dof && s.owner != WALL {
                    let o = s.owner as usize;
                    dg.bv[o] += kk * Vec2::new(txx, tyy);
                    let p = perp(s.lever);
                    let tp = s.apply_t(p);
                    dg.bw[o] += kk * p.dot(tp);
                }
            }
            dg.v[i] = d;
        }
        for o in 0..nb_unknown {
            dg.bv[o] += Vec2::splat(self.ball_mass);
            dg.bw[o] += self.ball_inertia;
        }
        let precond = |r: &State, z: &mut State| {
            for i in 0..r.v.len() {
                z.v[i] = Vec2::new(r.v[i].x / dg.v[i].x, r.v[i].y / dg.v[i].y);
            }
            for o in 0..r.bv.len() {
                z.bv[o] = Vec2::new(r.bv[o].x / dg.bv[o].x, r.bv[o].y / dg.bv[o].y);
                z.bw[o] = r.bw[o] / dg.bw[o];
            }
        };

        let ke_before = kinetic_energy_f64(&self.v, m);
        let b_norm = rhs.dot(&rhs).sqrt();
        if b_norm > 0.0 && b_norm.is_finite() {
            let target = CG_TOLERANCE * b_norm;
            let mut au = State::zeros(n, nb_unknown);
            sys.apply(&u, &mut au);
            let mut r = rhs.clone();
            let sub = |a: &mut State, b: &State| {
                for (x, y) in a.v.iter_mut().zip(&b.v) {
                    *x -= *y;
                }
                for (x, y) in a.bv.iter_mut().zip(&b.bv) {
                    *x -= *y;
                }
                for (x, y) in a.bw.iter_mut().zip(&b.bw) {
                    *x -= *y;
                }
            };
            sub(&mut r, &au);
            let mut z = State::zeros(n, nb_unknown);
            precond(&r, &mut z);
            let mut p = z.clone();
            let mut ap = State::zeros(n, nb_unknown);
            let mut rz_old = r.dot(&z);
            let mut iters = 0u32;
            if r.dot(&r).sqrt() > target {
                loop {
                    if iters >= CG_MAX_ITERATIONS {
                        stats.viscosity_cap_hits += 1;
                        break;
                    }
                    sys.apply(&p, &mut ap);
                    let p_ap = p.dot(&ap);
                    if p_ap <= 0.0 || !p_ap.is_finite() {
                        break;
                    }
                    let alpha = rz_old / p_ap;
                    for i in 0..n {
                        u.v[i] += alpha * p.v[i];
                        r.v[i] -= alpha * ap.v[i];
                    }
                    for o in 0..nb_unknown {
                        u.bv[o] += alpha * p.bv[o];
                        r.bv[o] -= alpha * ap.bv[o];
                        u.bw[o] += alpha * p.bw[o];
                        r.bw[o] -= alpha * ap.bw[o];
                    }
                    iters += 1;
                    if r.dot(&r).sqrt() <= target {
                        break;
                    }
                    precond(&r, &mut z);
                    let rz_new = r.dot(&z);
                    let beta_cg = rz_new / rz_old;
                    for i in 0..n {
                        p.v[i] = z.v[i] + beta_cg * p.v[i];
                    }
                    for o in 0..nb_unknown {
                        p.bv[o] = z.bv[o] + beta_cg * p.bv[o];
                        p.bw[o] = z.bw[o] + beta_cg * p.bw[o];
                    }
                    rz_old = rz_new;
                }
            }
            stats.viscosity_iterations += iters;
            self.v.copy_from_slice(&u.v);
            if dof {
                self.bodies.v.copy_from_slice(&u.bv);
                self.bodies.omega.copy_from_slice(&u.bw);
            }
        }
        stats.ke_delta_viscosity_j += kinetic_energy_f64(&self.v, m) - ke_before;

        // Solid reaction from the converged velocities: `J = -m dt k T (v_i - v_S)`.
        for i in 0..n {
            for idx in work.solids_of(i) {
                let s = &work.solids[idx];
                let vs = if dof && s.owner != WALL {
                    let o = s.owner as usize;
                    self.bodies.v[o] + self.bodies.omega[o] * perp(s.lever)
                } else {
                    s.vel
                };
                let j = -(m * dt * sys.k_solid[idx]) * s.apply_t(self.v[i] - vs);
                account_solid(stats, &mut self.ball_acc, self.x[i], s, j, false);
            }
        }
        // The balls' velocities already include this impulse.
        self.mark_ball_impulses_applied();
    }
}
