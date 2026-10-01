//! DFSPH (divergence-free SPH) slurry solver with Akinci boundary particles for the drum wall.
//!
//! Replaces the PBF solver ([`crate::pbf`]) in the planned rebuild (docs/VERIFICATION.md): PBF has no
//! solid wall density (particles pile up at the wall), reconstructs velocity from position
//! projection (creates and destroys energy) and applies no-slip only to wall-touching particles.
//! This module is fluid + drum wall only; balls are added in a later phase.
//!
//! # Algorithm (one internal step; Bender & Koschier 2015/2017 ordering)
//! 1. Neighbour search (fluid-fluid and fluid-boundary), density `rho_i = sum_j m W_ij + sum_b
//!    psi_b W_ib`, and the pressure-solver denominator `D_i = |sum_j m grad W_ij + sum_b psi_b
//!    grad W_ib|^2 + sum_j m^2 |grad W_ij|^2`.
//! 2. Divergence-free solve on `v` (removes positive volume change rate for particles with a full
//!    neighbourhood).
//! 3. Gravity and Akinci cohesion (surface tension), explicit.
//! 4. Implicit viscosity, `(I + dt L) v = v*` by conjugate gradient. `L` is built from the
//!    pairwise-central Monaghan/Cleary form (`f_ij = -C_ij e_ij e_ij . (v_i - v_j)`, `C_ij = 8 mu m^2
//!    / (rho_i rho_j) |W'| r / (r^2 + eta^2)`), so linear and angular momentum are conserved
//!    exactly and a rigid rotation feels no viscous force. The wall enters the same system as
//!    boundary pairs moving with `omega x r`, scaled by `slurry.wall_no_slip`.
//! 5. Constant-density solve (`rho*_i = rho_i + dt * d rho_i / dt` driven to `rest_density`,
//!    pressure `>= 0`, so a free surface exerts no suction), `kappa` accumulated and warm-started
//!    from the previous step. Pressure is mirrored onto the boundary particles (exact reaction).
//! 6. Advect; a hard SDF backstop catches any particle that crosses the wall (counted).
//!
//! Each call to [`Fluid::step`] is split into `n` equal internal steps chosen from a CFL bound
//! evaluated from the state at the start of the call (deterministic, bounded by
//! [`MAX_INTERNAL_STEPS`]).
//!
//! # Wall accounting
//! Every pressure and viscous interaction with the wall is a central pair force, so the impulse
//! the wall delivers to the fluid is summed per boundary particle; `wall_work_j = sum J . v_b`
//! equals `omega * wall_angular_impulse` and fluid angular momentum changes by exactly that
//! angular impulse plus gravity's torque.

mod boundary;
mod kernel;
mod neighbors;
mod pressure;
mod viscosity;

use glam::Vec2;

use crate::geometry::Drum;
use crate::params::{DyePattern, SlurryParams};
use boundary::Boundary;
use kernel::{cohesion_kernel, Kernel};
use neighbors::{build_fb, build_ff, Csr};

const GRAVITY: f32 = -9.81;
/// Number of polar-angle bins in [`FluidStats::wall_normal_impulse_by_angle`].
pub const WALL_IMPULSE_BINS: usize = 32;
/// Upper bound on internal steps per [`Fluid::step`] call.
pub const MAX_INTERNAL_STEPS: u32 = 32;
/// CFL factor: a particle may move at most this fraction of `dx` per internal step.
const CFL_FACTOR: f32 = 0.4;
/// Reference surface tension / cohesion constants, shared with the PBF solver so
/// `slurry.surface_tension_n_m` keeps its calibrated meaning.
const REFERENCE_SURFACE_TENSION_N_M: f32 = 0.072;
const COHESION_ACCEL_FACTOR: f32 = 2.0;
/// Floor on the pressure denominator `D_i`, as a fraction of its interior-lattice value.
const D_FLOOR_FRACTION: f32 = 0.05;

/// Per-call solver diagnostics. Energies in J per metre of depth, impulses in N*s per metre.
#[derive(Debug, Clone, Copy, Default)]
pub struct FluidStats {
    pub internal_steps: u32,
    /// Times the CFL bound asked for more than [`MAX_INTERNAL_STEPS`].
    pub cfl_cap_hits: u32,
    pub density_iterations: u32,
    pub divergence_iterations: u32,
    pub viscosity_iterations: u32,
    pub density_cap_hits: u32,
    pub divergence_cap_hits: u32,
    pub viscosity_cap_hits: u32,
    /// Mean / max relative density error (`max(0, rho* - rho0) / rho0`) at the last density solve.
    pub mean_density_error: f32,
    pub max_density_error: f32,
    /// Particles the hard SDF backstop had to move out of the wall.
    pub wall_backstop_hits: u32,
    /// Wall -> fluid impulse and its angular impulse about the drum centre.
    pub wall_impulse: [f64; 2],
    pub wall_angular_impulse: f64,
    /// The viscous (wall-friction) part of [`Self::wall_impulse`].
    pub wall_viscous_impulse: [f64; 2],
    /// Pressure impulse along the inward wall normal, binned by the boundary particle's polar
    /// angle (`atan2(y, x)` over `[-pi, pi)`).
    pub wall_normal_impulse_by_angle: [f64; WALL_IMPULSE_BINS],
    /// Work the wall does on the fluid (`sum J . v_b`): pressure and viscous parts.
    pub wall_work_j: f64,
    pub wall_pressure_work_j: f64,
    pub wall_viscous_work_j: f64,
    /// Kinetic-energy change of each stage; `pe_delta_j` is the potential-energy change over the
    /// advect stage.
    pub ke_delta_gravity_j: f64,
    pub ke_delta_cohesion_j: f64,
    pub ke_delta_viscosity_j: f64,
    pub ke_delta_pressure_j: f64,
    pub ke_delta_backstop_j: f64,
    pub pe_delta_j: f64,
}

/// Per-internal-step working data (neighbours, kernel gradients, density, denominators).
pub(crate) struct Work {
    pub ff: Csr,
    pub fb: Csr,
    /// `grad_i W_ij` per directed fluid-fluid edge / `grad_i W_ib` per fluid-boundary edge.
    pub grad_ff: Vec<Vec2>,
    pub grad_fb: Vec<Vec2>,
    /// Boundary particles' world positions and velocities.
    pub bw: Vec<Vec2>,
    pub vb: Vec<Vec2>,
    pub rho: Vec<f32>,
    pub inv_d: Vec<f32>,
}

pub struct Fluid {
    pub x: Vec<Vec2>,
    pub v: Vec<Vec2>,
    pub dye: Vec<f32>,
    pub particle_mass: f32,
    /// Lattice spacing `dx = R / resolution`.
    pub dx: f32,
    /// Kernel support radius `H = 2 dx`.
    pub h: f32,
    pub rest_density: f32,
    kernel: Kernel,
    boundary: Boundary,
    /// Accumulated density-solve pressure variable `kappa = p / rho^2` (warm start).
    kappa: Vec<f32>,
    row_h: f32,
    d_floor: f32,
    /// Drum angle at the end of the last step (for [`Fluid::densities`]).
    angle: f32,
}

pub(crate) fn kinetic_energy_f64(v: &[Vec2], mass: f32) -> f64 {
    0.5 * mass as f64 * v.iter().map(|v| v.length_squared() as f64).sum::<f64>()
}

impl Fluid {
    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Seeds a hexagonal lattice (spacing `dx = R / resolution`) filling `slurry.fill_fraction` of
    /// the drum's area from the bottom up, keeping at least `row_h / 2` clear of the wall.
    /// `drum_angle` is the drum's current rotation.
    pub fn new(slurry: &SlurryParams, drum: &Drum, drum_angle: f32, resolution: u32) -> Self {
        let r = drum.radius_m;
        let dx = r / resolution.max(1) as f32;
        let support = 2.0 * dx;
        let row_h = dx * (3.0f32.sqrt() * 0.5);
        let kernel = Kernel::new(support);
        let rho0 = slurry.density_kg_m3;

        // Interior lattice sums set the particle mass so a perfect lattice reads exactly `rho0`.
        let (mut rho_unit, mut grad_sq_unit) = (0.0f32, 0.0f32);
        for row in -4i32..=4 {
            let off = if row.rem_euclid(2) == 0 { 0.0 } else { 0.5 };
            for col in -4i32..=4 {
                let p = Vec2::new((col as f32 + off) * dx, row as f32 * row_h);
                let d = p.length();
                rho_unit += kernel.w(d);
                grad_sq_unit += kernel.grad(p, d).length_squared();
            }
        }
        let mass = rho0 / rho_unit;
        let d_floor = D_FLOOR_FRACTION * mass * mass * grad_sq_unit;

        let target_area = slurry.fill_fraction * std::f32::consts::PI * r * r;
        let target_count = (target_area / (dx * row_h)).round().max(0.0) as usize;
        let mut x: Vec<Vec2> = Vec::with_capacity(target_count);
        let mut dye: Vec<f32> = Vec::with_capacity(target_count);
        let clearance = 0.5 * row_h;
        if target_count > 0 {
            let n_rows = (2.0 * r / row_h).ceil() as i32 + 1;
            let n_cols = (r / dx).ceil() as i32 + 1;
            'rows: for row in -n_rows / 2..=n_rows / 2 {
                let y = row as f32 * row_h;
                let off = if row.rem_euclid(2) == 0 { 0.0 } else { 0.5 };
                for col in -n_cols..=n_cols {
                    let p = Vec2::new((col as f32 + off) * dx, y);
                    if drum.sdf_world(p, drum_angle).0 < clearance {
                        continue;
                    }
                    dye.push(dye_value(slurry.dye_pattern, p));
                    x.push(p);
                    if x.len() >= target_count {
                        break 'rows;
                    }
                }
            }
        }

        let n = x.len();
        let boundary = Boundary::build(drum, dx, row_h, support, mass);
        Self {
            x,
            v: vec![Vec2::ZERO; n],
            dye,
            particle_mass: mass,
            dx,
            h: support,
            rest_density: rho0,
            kernel,
            boundary,
            kappa: vec![0.0; n],
            row_h,
            d_floor,
            angle: drum_angle,
        }
    }

    /// Fluid mechanical energy (J per metre): kinetic + `m g y`.
    pub fn mechanical_energy_j(&self) -> f64 {
        let potential: f64 = self.x.iter().map(|p| p.y as f64).sum::<f64>()
            * self.particle_mass as f64
            * GRAVITY.abs() as f64;
        kinetic_energy_f64(&self.v, self.particle_mass) + potential
    }

    /// Per-particle SPH density **including the boundary term**, at the drum angle the last
    /// [`Fluid::step`] ended on.
    pub fn densities(&self) -> Vec<f32> {
        let n = self.len();
        if n == 0 {
            return Vec::new();
        }
        let (bw, _) = self.boundary_world(self.angle, 0.0);
        let to_local = Vec2::from_angle(-self.angle);
        let ff = build_ff(&self.x, self.h);
        let fb = build_fb(&self.x, to_local, &self.boundary, &bw, self.h);
        (0..n)
            .map(|i| {
                let mut rho = self.particle_mass * self.kernel.w(0.0);
                for k in ff.range(i) {
                    let j = ff.nbrs[k] as usize;
                    rho += self.particle_mass * self.kernel.w((self.x[i] - self.x[j]).length());
                }
                for k in fb.range(i) {
                    let b = fb.nbrs[k] as usize;
                    rho += self.boundary.psi[b] * self.kernel.w((self.x[i] - bw[b]).length());
                }
                rho
            })
            .collect()
    }

    fn boundary_world(&self, angle: f32, omega: f32) -> (Vec<Vec2>, Vec<Vec2>) {
        let rot = Vec2::from_angle(angle);
        let bw: Vec<Vec2> = self.boundary.local.iter().map(|&p| rot.rotate(p)).collect();
        let vb = bw
            .iter()
            .map(|p| omega * Vec2::new(-p.y, p.x))
            .collect::<Vec<_>>();
        (bw, vb)
    }

    /// Advances the fluid by `dt` against `drum` (at `drum_angle` at the start of the call).
    pub fn step(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        slurry: &SlurryParams,
        dt: f32,
    ) -> FluidStats {
        let mut stats = FluidStats::default();
        let n = self.len();
        if n == 0 || dt <= 0.0 {
            self.angle = drum_angle + drum.omega * dt.max(0.0);
            return stats;
        }
        if !self.boundary.matches(drum) {
            self.boundary = Boundary::build(drum, self.dx, self.row_h, self.h, self.particle_mass);
        }
        for (x, v) in self.x.iter_mut().zip(self.v.iter_mut()) {
            if !x.is_finite() || !v.is_finite() {
                *x = Vec2::ZERO;
                *v = Vec2::ZERO;
            }
        }

        let v_max = self
            .v
            .iter()
            .map(|v| v.length())
            .fold(0.0f32, f32::max)
            .max(drum.omega.abs() * drum.radius_m)
            + GRAVITY.abs() * dt;
        let wanted = (dt * v_max / (CFL_FACTOR * self.dx)).ceil().max(1.0) as u32;
        let n_int = if wanted > MAX_INTERNAL_STEPS {
            stats.cfl_cap_hits = 1;
            MAX_INTERNAL_STEPS
        } else {
            wanted
        };
        let dt_int = dt / n_int as f32;
        stats.internal_steps = n_int;

        let mut angle = drum_angle;
        for _ in 0..n_int {
            self.internal_step(drum, angle, slurry, dt_int, &mut stats);
            angle += drum.omega * dt_int;
        }
        self.angle = angle;
        stats
    }

    fn internal_step(
        &mut self,
        drum: &Drum,
        angle: f32,
        slurry: &SlurryParams,
        dt: f32,
        stats: &mut FluidStats,
    ) {
        let work = self.prepare(angle, drum.omega);
        let m = self.particle_mass;

        // Gravity and cohesion, then the density solve: pressure (not wall friction) must carry
        // the weight, so it acts on the raw gravity kick before viscosity does.
        let ke0 = kinetic_energy_f64(&self.v, m);
        for v in self.v.iter_mut() {
            v.y += GRAVITY * dt;
        }
        let ke1 = kinetic_energy_f64(&self.v, m);
        stats.ke_delta_gravity_j += ke1 - ke0;
        self.apply_cohesion(&work, slurry, dt);
        stats.ke_delta_cohesion_j += kinetic_energy_f64(&self.v, m) - ke1;

        self.pressure_solve(&work, pressure::Mode::Density, dt, stats);
        self.viscosity_solve(&work, slurry, dt, stats);
        // The stiff implicit viscosity damps the decompressing velocity the first solve just
        // produced; a second density solve (after it, so nothing damps it) restores the density.
        // Measured: this ordering is the only one of three tried (viscosity first; pressure,
        // viscosity, divergence; incremental pressure-correction) whose rotating-drum power is
        // stable under halving/quartering `dt` and has a quiet hydrostatic state.
        self.pressure_solve(&work, pressure::Mode::DensityRefine, dt, stats);

        // Advect, then the hard wall backstop against the drum at the end of the step.
        let angle_next = angle + drum.omega * dt;
        let pe0: f64 = self.x.iter().map(|p| p.y as f64).sum();
        for (x, v) in self.x.iter_mut().zip(&self.v) {
            *x += *v * dt;
        }
        let ke_pre = kinetic_energy_f64(&self.v, m);
        for (x, v) in self.x.iter_mut().zip(self.v.iter_mut()) {
            let (d, normal) = drum.sdf_world(*x, angle_next);
            if d < 0.0 {
                *x += -d * normal;
                let rel = *v - drum.wall_velocity(*x);
                let vn = rel.dot(normal);
                if vn < 0.0 {
                    *v -= vn * normal;
                }
                stats.wall_backstop_hits += 1;
            }
        }
        stats.ke_delta_backstop_j += kinetic_energy_f64(&self.v, m) - ke_pre;
        let pe1: f64 = self.x.iter().map(|p| p.y as f64).sum();
        stats.pe_delta_j += (pe1 - pe0) * m as f64 * GRAVITY.abs() as f64;
    }

    /// Neighbour search, kernel gradients, density and pressure denominators at the current
    /// positions and drum angle.
    fn prepare(&self, angle: f32, omega: f32) -> Work {
        let n = self.len();
        let m = self.particle_mass;
        let (bw, vb) = self.boundary_world(angle, omega);
        let to_local = Vec2::from_angle(-angle);
        let ff = build_ff(&self.x, self.h);
        let fb = build_fb(&self.x, to_local, &self.boundary, &bw, self.h);
        let mut grad_ff = vec![Vec2::ZERO; ff.nbrs.len()];
        let mut grad_fb = vec![Vec2::ZERO; fb.nbrs.len()];
        let mut rho = vec![0.0f32; n];
        let mut inv_d = vec![0.0f32; n];
        for i in 0..n {
            let mut rho_i = m * self.kernel.w(0.0);
            let mut g_sum = Vec2::ZERO;
            let mut sq = 0.0f32;
            for k in ff.range(i) {
                let j = ff.nbrs[k] as usize;
                let delta = self.x[i] - self.x[j];
                let r = delta.length();
                rho_i += m * self.kernel.w(r);
                let g = self.kernel.grad(delta, r);
                grad_ff[k] = g;
                g_sum += m * g;
                sq += m * m * g.length_squared();
            }
            for k in fb.range(i) {
                let b = fb.nbrs[k] as usize;
                let delta = self.x[i] - bw[b];
                let r = delta.length();
                let psi = self.boundary.psi[b];
                rho_i += psi * self.kernel.w(r);
                let g = self.kernel.grad(delta, r);
                grad_fb[k] = g;
                g_sum += psi * g;
            }
            rho[i] = rho_i;
            inv_d[i] = 1.0 / (g_sum.length_squared() + sq).max(self.d_floor);
        }
        Work {
            ff,
            fb,
            grad_ff,
            grad_fb,
            bw,
            vb,
            rho,
            inv_d,
        }
    }

    /// Akinci fluid-fluid cohesion as an explicit velocity kick (pairwise central, so momentum is
    /// conserved).
    fn apply_cohesion(&mut self, work: &Work, slurry: &SlurryParams, dt: f32) {
        let accel = (slurry.surface_tension_n_m.max(0.0) / REFERENCE_SURFACE_TENSION_N_M)
            * COHESION_ACCEL_FACTOR
            * GRAVITY.abs();
        if accel <= 0.0 {
            return;
        }
        let m = self.particle_mass;
        let rho0 = self.rest_density;
        let mut dv = vec![Vec2::ZERO; self.len()];
        for (i, dvi) in dv.iter_mut().enumerate() {
            let mut pull = Vec2::ZERO;
            for k in work.ff.range(i) {
                let j = work.ff.nbrs[k] as usize;
                let delta = self.x[i] - self.x[j];
                let r = delta.length();
                if r <= 1e-9 {
                    continue;
                }
                pull -= (m * cohesion_kernel(r, self.h) / rho0) * (delta / r);
            }
            *dvi = accel * pull * dt;
        }
        for (v, d) in self.v.iter_mut().zip(&dv) {
            *v += *d;
        }
    }
}

fn dye_value(pattern: DyePattern, p: Vec2) -> f32 {
    match pattern {
        DyePattern::LeftRight => {
            if p.x < 0.0 {
                0.0
            } else {
                1.0
            }
        }
        DyePattern::TopBottom => {
            if p.y < 0.0 {
                0.0
            } else {
                1.0
            }
        }
        DyePattern::None => 0.0,
    }
}

/// Accumulates one wall-pair impulse `j` (wall -> fluid) acting at boundary particle `b`.
pub(crate) fn account_wall_pair(
    stats: &mut FluidStats,
    j: Vec2,
    bw: Vec2,
    vb: Vec2,
    pressure: bool,
) {
    stats.wall_impulse[0] += j.x as f64;
    stats.wall_impulse[1] += j.y as f64;
    stats.wall_angular_impulse += (bw.x * j.y - bw.y * j.x) as f64;
    let work = j.dot(vb) as f64;
    stats.wall_work_j += work;
    if pressure {
        stats.wall_pressure_work_j += work;
        let len = bw.length();
        if len > 1e-9 {
            let inward = -bw / len;
            let angle = bw.y.atan2(bw.x);
            let bin = (((angle + std::f32::consts::PI) / std::f32::consts::TAU)
                * WALL_IMPULSE_BINS as f32) as usize;
            stats.wall_normal_impulse_by_angle[bin.min(WALL_IMPULSE_BINS - 1)] +=
                j.dot(inward) as f64;
        }
    } else {
        stats.wall_viscous_impulse[0] += j.x as f64;
        stats.wall_viscous_impulse[1] += j.y as f64;
        stats.wall_viscous_work_j += work;
    }
}

#[cfg(test)]
mod tests;
