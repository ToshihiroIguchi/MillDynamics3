//! DFSPH (divergence-free SPH) slurry solver with smooth solid boundaries.
//!
//! Replaces the PBF solver ([`crate::pbf`]) in the planned rebuild (docs/VERIFICATION.md): PBF has no
//! solid wall density (particles pile up at the wall), reconstructs velocity from position
//! projection (creates and destroys energy) and applies no-slip only to wall-touching particles.
//!
//! # Solid boundaries
//! The drum wall, lifters and balls are not particle sets but smooth solids ([`solid`]): the kernel
//! integrals of the solid region, `Gamma` (volume fraction) and `T` (viscous coupling tensor), are
//! tabulated once per body radius. A discrete boundary-particle layer was tried first and failed
//! the Taylor-Couette torque gate (docs/VERIFICATION.md): its rows, gap to the first fluid row and
//! staggering enter the viscous coupling at O(1) and do not converge cleanly.
//!
//! # Algorithm (one internal step)
//! 1. Neighbour search (fluid-fluid), solid samples per particle, density `rho_i = m sum_j W_ij +
//!    rho0 sum_S Gamma_S`, and the pressure denominator `D_i = |sum_j m grad W_ij + rho0 sum_S
//!    grad Gamma_S|^2 + sum_j m^2 |grad W_ij|^2`.
//! 2. Gravity and Akinci cohesion (surface tension), explicit.
//! 3. Constant-density solve (`rho*_i = rho_i + dt * d rho_i / dt` driven to `rest_density`,
//!    pressure `>= 0`, so a free surface exerts no suction), `kappa` accumulated and warm-started
//!    from the previous step; pressure acts on the solids with exact reaction.
//! 4. Implicit viscosity, `(I + dt L) v = v*` by conjugate gradient, `L` from the pairwise-central
//!    Monaghan/Cleary form (`f_ij = -C_ij e e . (v_i - v_j)`, `C_ij = 8 mu m^2 / (rho_i rho_j)
//!    |W'| r / (r^2 + eta^2)`): momentum and angular momentum are conserved exactly and a rigid
//!    rotation feels no viscous force. Solids enter as `-(8 mu / rho_i) T (v_i - v_S)`.
//! 5. A second density solve (viscosity damps the decompressing velocity of the first; the
//!    ordering was chosen by measurement, see [`Fluid::internal_step`]).
//! 6. Advect; hard backstops catch any particle that crosses the wall or enters a ball (counted).
//!
//! Each call to [`Fluid::step`] is split into `n` equal internal steps chosen from a CFL bound
//! evaluated from the state at the start of the call (deterministic, bounded by
//! [`MAX_INTERNAL_STEPS`]).
//!
//! # Accounting
//! Every pressure and viscous interaction with a solid is a central force, so the impulse a solid
//! delivers to the fluid is summed per (particle, solid) sample; `wall_work_j = sum J . v_S`
//! equals `omega * wall_angular_impulse`, and a ball receives `-J` with torque `-cross(x_i - x_B, J)`.

mod kernel;
mod neighbors;
mod pressure;
mod solid;
mod viscosity;

use glam::Vec2;

use crate::coupling::CouplingImpulses;
use crate::dem::Balls;
use crate::geometry::Drum;
use crate::grid::UniformGrid;
use crate::params::{DyePattern, SlurryParams};
use kernel::{cohesion_kernel, Kernel};
use neighbors::{build_ff, Csr};
use solid::{BallModel, SolidSample, WallModel, WALL};

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
/// Largest summed solid volume fraction a fluid particle may sit at (see the pore handling in
/// the advect stage of [`Fluid::internal_step`]).
const GAMMA_MAX: f32 = 0.85;
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
    /// Particles the hard SDF backstop had to move out of the wall / out of a ball.
    pub wall_backstop_hits: u32,
    pub ball_backstop_hits: u32,
    /// Work the balls do on the fluid (`sum J . v_b` over ball pairs, mirror of `wall_work_j`).
    pub ball_work_j: f64,

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

/// Per-internal-step working data (neighbours, kernel gradients, density, denominators, solids).
pub(crate) struct Work {
    pub ff: Csr,
    /// `grad_i W_ij` per directed fluid-fluid edge.
    pub grad_ff: Vec<Vec2>,
    /// Solid samples of particle `i`: `solids[solid_off[i]..solid_off[i + 1]]`.
    pub solid_off: Vec<u32>,
    pub solids: Vec<SolidSample>,
    pub rho: Vec<f32>,
    pub inv_d: Vec<f32>,
}

impl Work {
    #[inline]
    pub fn solids_of(&self, i: usize) -> std::ops::Range<usize> {
        self.solid_off[i] as usize..self.solid_off[i + 1] as usize
    }
}

/// Private copy of the balls' kinematics the fluid sub-cycles with.
struct Bodies {
    x: Vec<Vec2>,
    v: Vec<Vec2>,
    theta: Vec<f32>,
    omega: Vec<f32>,
}

/// Per-ball impulse the fluid delivered to the ball this call (linear, and angular about the
/// ball centre).
pub(crate) struct BallAcc {
    pub impulse: Vec<Vec2>,
    pub angular: Vec<f32>,
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
    /// Multiplier on gravity (1 by default; a verification hook, e.g. to run a Taylor-Couette
    /// case with and without gravity).
    pub gravity_scale: f32,
    /// Upper bound on the viscous diffusion number `dt nu / dx^2` of one internal step (the
    /// pressure-viscosity splitting error is first order in it); `f32::INFINITY` (default) leaves
    /// the step count to the CFL bound alone.
    pub max_viscous_number: f32,
    /// Cap on internal steps per call (default [`MAX_INTERNAL_STEPS`]).
    pub max_internal_steps: u32,
    /// Treat the balls' velocities as prescribed (they act like the wall); otherwise they are
    /// unknowns of the implicit viscous solve and receive the pressure impulses.
    pub bodies_prescribed: bool,
    /// Lattice correction of the fluid-fluid viscous coefficient (1 / the lattice's effective
    /// viscosity ratio on a quadratic profile).
    pub(crate) viscosity_norm: f32,
    ball_mass: f32,
    ball_inertia: f32,
    /// Ball impulses already folded into the ball copy's velocities.
    applied: BallAcc,
    kernel: Kernel,
    wall: WallModel,
    drum: Drum,
    /// Accumulated density-solve pressure variable `kappa = p / rho^2` (warm start).
    kappa: Vec<f32>,
    d_floor: f32,
    /// Drum angle at the end of the last step (for [`Fluid::densities`]).
    angle: f32,
    bodies: Bodies,
    ball_model: BallModel,
    pub(crate) ball_acc: BallAcc,
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
        Self::new_with_balls(slurry, drum, drum_angle, resolution, &[], 0.0)
    }

    /// As [`Fluid::new`], skipping lattice sites that overlap a ball (`existing_balls` centres,
    /// common `ball_radius`).
    pub fn new_with_balls(
        slurry: &SlurryParams,
        drum: &Drum,
        drum_angle: f32,
        resolution: u32,
        existing_balls: &[Vec2],
        ball_radius: f32,
    ) -> Self {
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
                    if existing_balls
                        .iter()
                        .any(|&b| (p - b).length() < ball_radius + clearance)
                    {
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

        // The particle mass is fixed by the requested fluid volume, not by the lattice cell: the
        // fluid mass must be `rho0 * fill * area` exactly. A mass even ~1 % too large for the
        // space available makes the density solve push particles into the solid's kernel layer
        // (the smooth boundary has no discrete rows to stop them), shrinking the effective gap and
        // inflating every viscous torque; the count of lattice sites is only commensurate with the
        // domain by accident, so the mass absorbs the mismatch (a few percent at most).
        let mass = if !x.is_empty() {
            rho0 * target_area / x.len() as f32
        } else {
            mass
        };
        let n = x.len();
        let wall = WallModel::build(drum, &kernel);
        Self {
            x,
            v: vec![Vec2::ZERO; n],
            dye,
            particle_mass: mass,
            dx,
            h: support,
            rest_density: rho0,
            gravity_scale: 1.0,
            max_viscous_number: f32::INFINITY,
            max_internal_steps: MAX_INTERNAL_STEPS,
            bodies_prescribed: false,
            viscosity_norm: lattice_viscosity_norm(&kernel, dx, row_h, mass, rho0),
            ball_mass: 0.0,
            ball_inertia: 0.0,
            applied: BallAcc {
                impulse: Vec::new(),
                angular: Vec::new(),
            },
            kernel,
            wall,
            drum: *drum,
            kappa: vec![0.0; n],
            d_floor,
            angle: drum_angle,
            bodies: Bodies {
                x: Vec::new(),
                v: Vec::new(),
                theta: Vec::new(),
                omega: Vec::new(),
            },
            ball_model: BallModel::empty(),
            ball_acc: BallAcc {
                impulse: Vec::new(),
                angular: Vec::new(),
            },
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
        let (soff, solids) = self.sample_solids(&self.drum, self.angle);
        let ff = build_ff(&self.x, self.h);
        (0..n)
            .map(|i| {
                let mut rho = self.particle_mass * self.kernel.w(0.0);
                for k in ff.range(i) {
                    let j = ff.nbrs[k] as usize;
                    rho += self.particle_mass * self.kernel.w((self.x[i] - self.x[j]).length());
                }
                for s in &solids[soff[i] as usize..soff[i + 1] as usize] {
                    rho += self.rest_density * s.gamma;
                }
                rho
            })
            .collect()
    }

    /// Solid samples (wall first, then nearby balls in index order) for every fluid particle.
    fn sample_solids(&self, drum: &Drum, angle: f32) -> (Vec<u32>, Vec<SolidSample>) {
        let n = self.x.len();
        let mut off = Vec::with_capacity(n + 1);
        off.push(0u32);
        let mut samples: Vec<SolidSample> = Vec::new();
        let balls_present = !self.bodies.x.is_empty() && self.ball_model.radius > 0.0;
        let grid = if balls_present {
            Some(UniformGrid::build(
                &self.bodies.x,
                self.ball_model.reach().max(1e-6),
            ))
        } else {
            None
        };
        let mut nearby: Vec<u32> = Vec::new();
        for i in 0..n {
            let x = self.x[i];
            if let Some(s) = self.wall.sample(x, drum, angle) {
                if s.gamma > 0.0 {
                    samples.push(s);
                }
            }
            if let Some(grid) = &grid {
                nearby.clear();
                grid.for_each_near(x, |b| nearby.push(b));
                nearby.sort_unstable();
                for &b in &nearby {
                    let b = b as usize;
                    if let Some(s) = self.ball_model.sample(
                        b as u32,
                        x,
                        self.bodies.x[b],
                        self.bodies.v[b],
                        self.bodies.omega[b],
                    ) {
                        if s.gamma > 0.0 {
                            samples.push(s);
                        }
                    }
                }
            }
            off.push(samples.len() as u32);
        }
        (off, samples)
    }

    /// Advances the fluid by `dt` against `drum` (at `drum_angle` at the start of the call), with
    /// no balls.
    pub fn step(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        slurry: &SlurryParams,
        dt: f32,
    ) -> FluidStats {
        self.step_with_balls(drum, drum_angle, slurry, dt, &Balls::empty())
            .1
    }

    /// As [`Fluid::step`], with `balls` as moving rigid boundaries. The balls' kinematics are
    /// copied and sub-cycled with the fluid; the impulse (and angular impulse about the ball
    /// centre) the fluid delivers to each ball over the call is returned for
    /// [`crate::dem::DemState::step_with_external_forces`] to apply.
    pub fn step_with_balls(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        slurry: &SlurryParams,
        dt: f32,
        balls: &Balls,
    ) -> (CouplingImpulses, FluidStats) {
        let mut stats = FluidStats::default();
        let nb = balls.len();
        let n = self.len();
        self.bodies = Bodies {
            x: balls.x.clone(),
            v: balls.v.clone(),
            theta: balls.theta.clone(),
            omega: balls.omega.clone(),
        };
        if nb > 0 && (self.ball_model.radius - balls.radius).abs() > 1e-9 {
            self.ball_model = BallModel::build(balls.radius, &self.kernel);
        }
        self.ball_acc = BallAcc {
            impulse: vec![Vec2::ZERO; nb],
            angular: vec![0.0; nb],
        };
        self.applied = BallAcc {
            impulse: vec![Vec2::ZERO; nb],
            angular: vec![0.0; nb],
        };
        self.ball_mass = balls.mass;
        self.ball_inertia = balls.inertia;
        let finish = |acc: &BallAcc| CouplingImpulses {
            impulses: acc.impulse.clone(),
            angular_impulses: acc.angular.clone(),
            clamp_hits: 0,
            fluid_momentum_change: Vec2::ZERO,
        };
        if n == 0 || dt <= 0.0 {
            self.angle = drum_angle + drum.omega * dt.max(0.0);
            return (finish(&self.ball_acc), stats);
        }
        if !self.wall.matches(drum) {
            self.wall = WallModel::build(drum, &self.kernel);
        }
        self.drum = *drum;
        for (x, v) in self.x.iter_mut().zip(self.v.iter_mut()) {
            if !x.is_finite() || !v.is_finite() {
                *x = Vec2::ZERO;
                *v = Vec2::ZERO;
            }
        }

        let ball_speed = self
            .bodies
            .v
            .iter()
            .zip(&self.bodies.omega)
            .map(|(v, w)| v.length() + w.abs() * balls.radius)
            .filter(|s| s.is_finite())
            .fold(0.0f32, f32::max);
        let v_max = self
            .v
            .iter()
            .map(|v| v.length())
            .fold(0.0f32, f32::max)
            .max(drum.omega.abs() * drum.radius_m)
            .max(ball_speed)
            + GRAVITY.abs() * dt;
        let wanted_cfl = (dt * v_max / (CFL_FACTOR * self.dx)).ceil().max(1.0);
        let nu = slurry.viscosity_pa_s.max(0.0) / self.rest_density;
        let wanted_visc = if self.max_viscous_number.is_finite() && self.max_viscous_number > 0.0 {
            (dt * nu / (self.dx * self.dx * self.max_viscous_number))
                .ceil()
                .max(1.0)
        } else {
            1.0
        };
        let wanted = wanted_cfl.max(wanted_visc) as u32;
        let n_int = if wanted > self.max_internal_steps {
            stats.cfl_cap_hits = 1;
            self.max_internal_steps
        } else {
            wanted
        };
        let dt_int = dt / n_int as f32;
        stats.internal_steps = n_int;

        let mut angle = drum_angle;
        for _ in 0..n_int {
            self.internal_step(drum, angle, slurry, dt_int, balls.radius, &mut stats);
            angle += drum.omega * dt_int;
        }
        self.angle = angle;
        (finish(&self.ball_acc), stats)
    }

    fn internal_step(
        &mut self,
        drum: &Drum,
        angle: f32,
        slurry: &SlurryParams,
        dt: f32,
        ball_radius: f32,
        stats: &mut FluidStats,
    ) {
        let work = self.prepare(drum, angle);
        let m = self.particle_mass;

        // Gravity and cohesion, then the density solve: pressure (not wall friction) must carry
        // the weight, so it acts on the raw gravity kick before viscosity does.
        let ke0 = kinetic_energy_f64(&self.v, m);
        for v in self.v.iter_mut() {
            v.y += GRAVITY * self.gravity_scale * dt;
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
        // The balls advance first (their copy moved with the fluid-updated velocities), then any
        // fluid particle left inside a ball is moved back to its surface.
        for b in 0..self.bodies.x.len() {
            let v = self.bodies.v[b];
            let w = self.bodies.omega[b];
            self.bodies.x[b] += v * dt;
            self.bodies.theta[b] += w * dt;
        }
        if !self.bodies.x.is_empty() && ball_radius > 0.0 {
            let centres: Vec<Vec2> = self.bodies.x.clone();
            let grid = UniformGrid::build(&centres, (2.0 * ball_radius).max(1e-6));
            for (x, v) in self.x.iter_mut().zip(self.v.iter_mut()) {
                grid.for_each_near(*x, |b| {
                    let b = b as usize;
                    let delta = *x - centres[b];
                    let d = delta.length();
                    if d < ball_radius {
                        let normal = if d > 1e-9 { delta / d } else { Vec2::X };
                        *x = centres[b] + normal * ball_radius;
                        let surface = self.bodies.v[b]
                            + self.bodies.omega[b] * Vec2::new(-normal.y, normal.x) * ball_radius;
                        let vn = (*v - surface).dot(normal);
                        if vn < 0.0 {
                            *v -= vn * normal;
                        }
                        stats.ball_backstop_hits += 1;
                    }
                });
            }
        }
        // Pores narrower than the kernel (a ball against the wall or another ball) cannot hold a
        // fluid particle: the density solve would squeeze it out at an unphysical speed. A
        // particle where the solids' summed volume fraction exceeds [`GAMMA_MAX`] is moved down
        // the gradient to the boundary of that region and loses its inward velocity.
        if !self.bodies.x.is_empty() {
            let (off, solids) = self.sample_solids(drum, angle_next);
            for i in 0..self.x.len() {
                let range = off[i] as usize..off[i + 1] as usize;
                if range.len() < 2 {
                    continue; // a single solid is handled by the exact backstops above
                }
                let (mut gamma, mut grad, mut vel) = (0.0f32, Vec2::ZERO, Vec2::ZERO);
                for sol in &solids[range.clone()] {
                    gamma += sol.gamma;
                    grad += sol.grad;
                    vel += sol.vel;
                }
                let gl = grad.length();
                if gamma > GAMMA_MAX && gl > 1e-9 {
                    let n = grad / gl;
                    self.x[i] -= n * ((gamma - GAMMA_MAX) / gl).min(0.5 * self.dx);
                    let vn = (self.v[i] - vel / range.len() as f32).dot(n);
                    if vn > 0.0 {
                        self.v[i] -= vn * n;
                    }
                    stats.ball_backstop_hits += 1;
                }
            }
        }
        stats.ke_delta_backstop_j += kinetic_energy_f64(&self.v, m) - ke_pre;
        let pe1: f64 = self.x.iter().map(|p| p.y as f64).sum();
        stats.pe_delta_j += (pe1 - pe0) * m as f64 * GRAVITY.abs() as f64;
    }

    /// Neighbour search, solid samples, kernel gradients, density and pressure denominators at the
    /// current positions and drum angle.
    fn prepare(&self, drum: &Drum, angle: f32) -> Work {
        let n = self.len();
        let m = self.particle_mass;
        let rho0 = self.rest_density;
        let (solid_off, solids) = self.sample_solids(drum, angle);
        let ff = build_ff(&self.x, self.h);
        let mut grad_ff = vec![Vec2::ZERO; ff.nbrs.len()];
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
            for s in &solids[solid_off[i] as usize..solid_off[i + 1] as usize] {
                rho_i += rho0 * s.gamma;
                g_sum += rho0 * s.grad;
            }
            rho[i] = rho_i;
            inv_d[i] = 1.0 / (g_sum.length_squared() + sq).max(self.d_floor);
        }
        Work {
            ff,
            grad_ff,
            solid_off,
            solids,
            rho,
            inv_d,
        }
    }

    /// Folds ball impulses accumulated since the last sync (the pressure stages') into the ball
    /// copy's velocities; a no-op for prescribed balls.
    pub(super) fn apply_ball_impulses(&mut self) {
        if self.bodies_prescribed || self.ball_mass <= 0.0 {
            return;
        }
        for b in 0..self.bodies.x.len() {
            let dp = self.ball_acc.impulse[b] - self.applied.impulse[b];
            let dl = self.ball_acc.angular[b] - self.applied.angular[b];
            self.bodies.v[b] += dp / self.ball_mass;
            if self.ball_inertia > 0.0 {
                self.bodies.omega[b] += dl / self.ball_inertia;
            }
        }
        self.mark_ball_impulses_applied();
    }

    pub(super) fn mark_ball_impulses_applied(&mut self) {
        self.applied.impulse.copy_from_slice(&self.ball_acc.impulse);
        self.applied.angular.copy_from_slice(&self.ball_acc.angular);
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

/// `1 / (effective viscosity / nu)` of the pairwise-central operator on the seeded hex lattice
/// for a quadratic velocity profile (0.965 for the shipped kernel: the `eta^2` regularisation and
/// the finite sum), applied to the fluid-fluid coefficient so the bulk reproduces `mu`.
fn lattice_viscosity_norm(kernel: &Kernel, dx: f32, row_h: f32, mass: f32, rho0: f32) -> f32 {
    let eta2 = solid::ETA_FACTOR * kernel.support * kernel.support;
    let mut acc_x = 0.0f32;
    for row in -5i32..=5 {
        let off = if row.rem_euclid(2) == 0 { 0.0 } else { 0.5 };
        for col in -5i32..=5 {
            let pj = Vec2::new((col as f32 + off) * dx, row as f32 * row_h);
            let r = pj.length();
            if r < 1e-9 || r >= kernel.support {
                continue;
            }
            let e = -pj / r;
            // u = (y^2 / 2, 0): u_i - u_j at the origin is -pj.y^2 / 2
            let c = 8.0 * mass / (rho0 * rho0) * kernel.dw_dr(r).abs() * r / (r * r + eta2) * rho0;
            acc_x += -c * e.x * e.x * (-0.5 * pj.y * pj.y) * 1.0;
        }
    }
    // acc_x is `a_x / nu` for `mu = rho0 nu`; exact is 1.
    if acc_x > 0.1 {
        1.0 / acc_x
    } else {
        1.0
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

/// Accumulates one (particle, solid) impulse `j` (solid -> fluid) on the particle at `x`: wall
/// samples into the wall statistics (impulse, angular impulse about the drum centre, work, pressure
/// profile by the particle's polar angle); ball samples into the ball's impulse and angular impulse
/// about its centre (Newton's third law) and `ball_work_j`.
pub(crate) fn account_solid(
    stats: &mut FluidStats,
    balls: &mut BallAcc,
    x: Vec2,
    s: &SolidSample,
    j: Vec2,
    pressure: bool,
) {
    if s.owner != WALL {
        let o = s.owner as usize;
        balls.impulse[o] -= j;
        balls.angular[o] -= s.lever.x * j.y - s.lever.y * j.x;
        stats.ball_work_j += j.dot(s.vel) as f64;
        return;
    }
    stats.wall_impulse[0] += j.x as f64;
    stats.wall_impulse[1] += j.y as f64;
    stats.wall_angular_impulse += (x.x * j.y - x.y * j.x) as f64;
    let work_j = j.dot(s.vel) as f64;
    stats.wall_work_j += work_j;
    if pressure {
        stats.wall_pressure_work_j += work_j;
        let len = x.length();
        if len > 1e-9 {
            let inward = -x / len;
            let angle = x.y.atan2(x.x);
            let bin = (((angle + std::f32::consts::PI) / std::f32::consts::TAU)
                * WALL_IMPULSE_BINS as f32) as usize;
            stats.wall_normal_impulse_by_angle[bin.min(WALL_IMPULSE_BINS - 1)] +=
                j.dot(inward) as f64;
        }
    } else {
        stats.wall_viscous_impulse[0] += j.x as f64;
        stats.wall_viscous_impulse[1] += j.y as f64;
        stats.wall_viscous_work_j += work_j;
    }
}

#[cfg(test)]
mod tests;
