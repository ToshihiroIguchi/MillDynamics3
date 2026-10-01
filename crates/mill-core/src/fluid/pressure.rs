//! DFSPH divergence-free and constant-density pressure solves (Jacobi, accumulated `kappa`).
//!
//! Pressure acceleration `a_i = -sum_j m (kappa_i + kappa_j) grad W_ij - kappa_i sum_b psi_b grad
//! W_ib`, `kappa = p / rho^2`. Driving particle `i`'s own error `e_i` to zero with only `kappa_i`
//! gives `kappa_i = e_i / (dt^2 D_i)` (`e_i` a density change over `dt`), which is iterated
//! Jacobi-style. `kappa` is accumulated and clamped at `>= 0` (no suction), so an over-corrected
//! particle can relax back toward zero pressure.

use glam::Vec2;

use super::{account_solid, kinetic_energy_f64, Fluid, FluidStats, Work};

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Mode {
    Divergence,
    /// Constant-density solve, warm-started from the stored `kappa`.
    Density,
    /// A second density solve in the same step (after viscosity): continues from the
    /// accumulated `kappa` without re-applying it.
    DensityRefine,
}

/// Density solve stops once the mean relative error of the pressure-active particles is below
/// this (0.01 %).
const DENSITY_TOLERANCE: f32 = 1e-4;
/// ... and the worst particle below this (0.5 %): a mean criterion alone lets a few wall
/// particles stay badly over-compressed.
const DENSITY_MAX_ERROR: f32 = 5e-3;
const DENSITY_MIN_ITERATIONS: u32 = 2;
const DENSITY_MAX_ITERATIONS: u32 = 200;
/// Divergence solve: mean relative density change per internal step.
const DIVERGENCE_TOLERANCE: f32 = 5e-4;
const DIVERGENCE_MAX_ERROR: f32 = 5e-3;
const DIVERGENCE_MIN_ITERATIONS: u32 = 1;
const DIVERGENCE_MAX_ITERATIONS: u32 = 100;
/// Jacobi under-relaxation of the per-iteration pressure increment.
const RELAXATION: f32 = 0.5;
/// Particles below this fraction of rest density (free surface) are skipped by the divergence
/// solve.
const DIVERGENCE_MIN_DENSITY_RATIO: f32 = 0.9;

impl Fluid {
    pub(super) fn pressure_solve(
        &mut self,
        work: &Work,
        mode: Mode,
        dt: f32,
        stats: &mut FluidStats,
    ) {
        let n = self.len();
        let m = self.particle_mass;
        let rho0 = self.rest_density;
        let (min_it, max_it, tol, tol_max) = match mode {
            Mode::Density | Mode::DensityRefine => (
                DENSITY_MIN_ITERATIONS,
                DENSITY_MAX_ITERATIONS,
                DENSITY_TOLERANCE,
                DENSITY_MAX_ERROR,
            ),
            Mode::Divergence => (
                DIVERGENCE_MIN_ITERATIONS,
                DIVERGENCE_MAX_ITERATIONS,
                DIVERGENCE_TOLERANCE,
                DIVERGENCE_MAX_ERROR,
            ),
        };
        let ke_before = kinetic_energy_f64(&self.v, m);
        let mut kappa = match mode {
            Mode::Density | Mode::DensityRefine => std::mem::take(&mut self.kappa),
            Mode::Divergence => vec![0.0f32; n],
        };
        let mut kinc = vec![0.0f32; n];
        if mode == Mode::Density {
            kinc.copy_from_slice(&kappa);
            self.apply_kappa(work, &kinc, dt, stats);
        }

        let inv_dt2 = 1.0 / (dt * dt);
        let mut iters = 0u32;
        loop {
            let mut err_sum = 0.0f64;
            let mut err_max = 0.0f32;
            for i in 0..n {
                let mut rho_dot = 0.0f32;
                for k in work.ff.range(i) {
                    let j = work.ff.nbrs[k] as usize;
                    rho_dot += m * (self.v[i] - self.v[j]).dot(work.grad_ff[k]);
                }
                for idx in work.solids_of(i) {
                    let sol = &work.solids[idx];
                    rho_dot += rho0 * (self.v[i] - sol.vel).dot(sol.grad);
                }
                let e = match mode {
                    Mode::Density | Mode::DensityRefine => work.rho[i] + dt * rho_dot - rho0,
                    Mode::Divergence => {
                        if work.rho[i] >= DIVERGENCE_MIN_DENSITY_RATIO * rho0 {
                            dt * rho_dot
                        } else {
                            0.0
                        }
                    }
                };
                // A particle with no accumulated pressure can only be pushed apart; one that
                // already carries pressure may be relaxed back (down to zero).
                let e = if kappa[i] > 0.0 { e } else { e.max(0.0) };
                let mut inc = RELAXATION * e * work.inv_d[i] * inv_dt2;
                if kappa[i] + inc < 0.0 {
                    inc = -kappa[i];
                }
                kinc[i] = inc;
                let rel = e.abs() / rho0;
                err_sum += rel as f64;
                err_max = err_max.max(rel);
            }
            let err_mean = (err_sum / n as f64) as f32;
            if mode != Mode::Divergence {
                stats.mean_density_error = err_mean;
                stats.max_density_error = err_max;
            }
            if iters >= min_it && err_mean <= tol && err_max <= tol_max {
                break;
            }
            if iters >= max_it {
                match mode {
                    Mode::Density | Mode::DensityRefine => stats.density_cap_hits += 1,
                    Mode::Divergence => stats.divergence_cap_hits += 1,
                }
                break;
            }
            self.apply_kappa(work, &kinc, dt, stats);
            for (k, inc) in kappa.iter_mut().zip(&kinc) {
                *k += *inc;
            }
            iters += 1;
        }

        match mode {
            Mode::Density | Mode::DensityRefine => {
                stats.density_iterations += iters;
                self.kappa = kappa;
            }
            Mode::Divergence => stats.divergence_iterations += iters,
        }
        stats.ke_delta_pressure_j += kinetic_energy_f64(&self.v, m) - ke_before;
        self.apply_ball_impulses();
    }

    /// Applies the pressure acceleration for the per-particle pressure variable `kinc` (Jacobi: all
    /// accelerations from the same `kinc`, then added) and accounts the wall reaction.
    fn apply_kappa(&mut self, work: &Work, kinc: &[f32], dt: f32, stats: &mut FluidStats) {
        let m = self.particle_mass;
        let rho0 = self.rest_density;
        let mut dv = vec![Vec2::ZERO; self.len()];
        for (i, dvi) in dv.iter_mut().enumerate() {
            let ki = kinc[i];
            let xi = self.x[i];
            let mut a = Vec2::ZERO;
            for k in work.ff.range(i) {
                let j = work.ff.nbrs[k] as usize;
                a -= m * (ki + kinc[j]) * work.grad_ff[k];
            }
            let mut dvb = Vec2::ZERO;
            if ki != 0.0 {
                for idx in work.solids_of(i) {
                    let sol = &work.solids[idx];
                    let dvib = -(dt * ki * rho0) * sol.grad;
                    dvb += dvib;
                    account_solid(stats, &mut self.ball_acc, xi, sol, m * dvib, true);
                }
            }
            *dvi = dt * a + dvb;
        }
        for (v, d) in self.v.iter_mut().zip(&dv) {
            *v += *d;
        }
    }
}
