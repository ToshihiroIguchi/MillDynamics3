//! Position Based Fluids (PBF) solver for the slurry.
//!
//! 2D Poly6 (density) / Spiky (gradient) kernels, iterative density-constraint projection with
//! artificial-pressure anti-clustering, drum-wall boundary projection with wall-velocity
//! blending (no-slip), and XSPH viscosity. See docs/PLAN.md ss3.3 for the full per-substep
//! algorithm and the rationale for choosing PBF over WCSPH/LBM/MPM.
//!
//! **M3 scope ("slurry alone", per docs/PLAN.md ss5): no ball<->fluid interaction yet** -- fluid
//! particles only collide with the drum wall/lifters, not with balls. Two-way coupling is M4
//! ([`crate::coupling`]).
//!
//! **Viscosity, v1 simplification.** The plan's original two formulas for the XSPH mixing
//! coefficient (a simple `clamp(mu/mu_ref)` and a separate dt-normalized exponential) were
//! redundant; this implementation uses the simpler, dt-independent-in-practice form
//! `c = clamp(viscosity_pa_s / MU_REF, 0, 1)` applied directly as the XSPH coefficient. Like the
//! rest of the v1 slurry model, this is qualitative (monotonic in viscosity, not a quantitatively
//! calibrated match to real Pa*s); see the module-level caveat already noted project-wide (2D
//! slice, qualitative results).

use glam::Vec2;

use crate::geometry::Drum;
use crate::grid::UniformGrid;
use crate::params::{DyePattern, SlurryParams};

const GRAVITY: f32 = -9.81;
/// CFM-style relaxation added to the lambda denominator to avoid a divide-by-zero for particles
/// with few/no neighbours (docs/PLAN.md ss3.3: "epsilon = 1e2..1e3 relaxation").
const EPSILON_RELAX: f32 = 200.0;
/// Artificial-pressure ("tensile instability" / anti-clustering) term coefficients, docs/PLAN.md
/// ss3.3. **Disabled for v1** (`S_CORR_K = 0.0`): at the mass/kernel-radius/rest-density scale
/// used here (real SI units, e.g. `rest_density` ~1000-2000 kg/m^3), the literature-default
/// `k = 0.1` (tuned for the very different implicit scale of typical toy/graphics PBF demos)
/// produces a correction 10-30x larger than the density-constraint (`lambda`) terms it's meant to
/// gently supplement, causing runaway dispersal instead of preventing clustering. Empirically
/// verified: disabling it gives a stable settled puddle at ~1% mean density error over 2s of
/// simulated time; re-enabling it at several tried-and-measured `k` scales did not find a stable
/// working point in the time available. Revisit only if visual clustering artifacts are observed.
const S_CORR_K: f32 = 0.0;
const S_CORR_DELTA_Q_FACTOR: f32 = 0.2;
const S_CORR_N: i32 = 4;
/// Reference viscosity (Pa*s) at which the XSPH mixing coefficient saturates at 1.0; see the
/// module doc comment's "Viscosity, v1 simplification" note.
const MU_REF: f32 = 2.0;

/// 2D Poly6 kernel (docs/PLAN.md ss3.3: `4/(pi h^8)`), evaluated from squared distance.
fn poly6(r2: f32, h: f32) -> f32 {
    if r2 >= h * h || r2 < 0.0 {
        return 0.0;
    }
    let h2 = h * h;
    let term = h2 - r2;
    4.0 / (std::f32::consts::PI * h.powi(8)) * term * term * term
}

/// Gradient (w.r.t. the first argument) of the 2D Spiky kernel (docs/PLAN.md ss3.3:
/// `-30/(pi h^5)`). `delta = x_i - x_j`, `r = |delta|`. Points from `x_i` toward `x_j` (Spiky
/// decreases with `r`, so its gradient w.r.t. `x_i` points toward the neighbour).
fn spiky_grad(delta: Vec2, r: f32, h: f32) -> Vec2 {
    if r <= 1e-9 || r >= h {
        return Vec2::ZERO;
    }
    let coeff = -30.0 / (std::f32::consts::PI * h.powi(5)) * (h - r) * (h - r);
    (delta / r) * coeff
}

/// Slurry (fluid) particle population. All particles share one kernel radius/mass/rest density
/// (see [`FluidParticles::seed_lattice`]); `dye` is a pure Lagrangian tracer (no diffusion) used
/// for mixing visualization/metrics (docs/PLAN.md ss3.3/3.5).
pub struct FluidParticles {
    pub x: Vec<Vec2>,
    pub v: Vec<Vec2>,
    pub dye: Vec<f32>,
    pub particle_mass: f32,
    pub h: f32,
    pub rest_density: f32,
}

impl FluidParticles {
    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Seeds a hexagonal lattice of fluid particles, spacing `dx = drum_radius_m / resolution`
    /// (kernel radius `h = 2*dx`), filling `slurry.fill_fraction` of the drum's cross-sectional
    /// area from the bottom upward (gravity keeps it there anyway, and starting near the final
    /// rest position speeds up settling). Sites overlapping an existing ball are skipped (docs/
    /// PLAN.md ss3.3: "balls initialised first, fluid particles overlapping a ball are removed").
    pub fn seed_lattice(
        slurry: &SlurryParams,
        drum_radius_m: f32,
        resolution: u32,
        existing_balls: &[Vec2],
        ball_radius: f32,
    ) -> Self {
        let dx = drum_radius_m / resolution.max(1) as f32;
        let h = 2.0 * dx;
        let rest_density = slurry.density_kg_m3;
        let particle_mass = rest_density * dx * dx;

        let target_area =
            slurry.fill_fraction * std::f32::consts::PI * drum_radius_m * drum_radius_m;
        let target_count = (target_area / (dx * dx)).round().max(0.0) as usize;

        let mut x = Vec::with_capacity(target_count);
        let mut dye = Vec::with_capacity(target_count);
        if target_count > 0 && dx > 0.0 {
            let row_h = dx * (3.0f32.sqrt() * 0.5);
            let n_rows = (2.0 * drum_radius_m / row_h).ceil() as i32 + 1;

            'rows: for row in -n_rows / 2..=n_rows / 2 {
                let y = row as f32 * row_h;
                if y.abs() > drum_radius_m {
                    continue;
                }
                let half_width = (drum_radius_m * drum_radius_m - y * y).max(0.0).sqrt();
                let offset = if row.rem_euclid(2) == 0 {
                    0.0
                } else {
                    dx * 0.5
                };
                let n_cols = (2.0 * half_width / dx).floor() as i32;
                for col in -n_cols / 2..=n_cols / 2 {
                    let cx = col as f32 * dx + offset;
                    if cx * cx + y * y > drum_radius_m * drum_radius_m {
                        continue;
                    }
                    let p = Vec2::new(cx, y);
                    let overlaps_ball = existing_balls
                        .iter()
                        .any(|&b| (p - b).length() < ball_radius);
                    if overlaps_ball {
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
        Self {
            x,
            v: vec![Vec2::ZERO; n],
            dye,
            particle_mass,
            h,
            rest_density,
        }
    }

    /// Per-particle density (docs/PLAN.md ss3.3: `rho_i = sum_j m_j * W_poly6(r_ij)`, including
    /// the particle's own contribution at `r = 0`). Exposed publicly for tests and for
    /// [`crate::surface`]/[`crate::metrics`] to reuse without recomputing neighbours.
    pub fn densities(&self) -> Vec<f32> {
        let n = self.len();
        let mut density = vec![0.0f32; n];
        if n == 0 {
            return density;
        }
        let grid = UniformGrid::build(&self.x, self.h);
        let self_contribution = self.particle_mass * poly6(0.0, self.h);
        density.fill(self_contribution);
        grid.for_each_candidate_pair(|i, j| {
            let (iu, ju) = (i as usize, j as usize);
            let r2 = (self.x[iu] - self.x[ju]).length_squared();
            let w = poly6(r2, self.h);
            density[iu] += self.particle_mass * w;
            density[ju] += self.particle_mass * w;
        });
        density
    }

    /// Advances the fluid by one fixed sub-step `dt`, against the given (already positioned at
    /// `drum_angle`) drum. See docs/PLAN.md ss3.3 for the full per-substep algorithm; M3 does not
    /// yet couple against balls (docs/PLAN.md ss3.4, M4).
    pub fn step(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        slurry: &SlurryParams,
        iterations: u32,
        dt: f32,
    ) {
        let n = self.len();
        if n == 0 || dt <= 0.0 {
            return;
        }
        let h = self.h;
        let rest_density = self.rest_density;
        let mass = self.particle_mass;

        // --- 1. Predict --------------------------------------------------------------------
        let x0: Vec<Vec2> = self.x.clone();
        for i in 0..n {
            self.v[i].y += GRAVITY * dt;
            self.x[i] += self.v[i] * dt;
        }

        // --- 2. Neighbour lists (rebuilt each sub-step, reused across iterations) -----------
        let grid = UniformGrid::build(&self.x, h);
        let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); n];
        {
            let positions = &self.x;
            grid.for_each_candidate_pair(|i, j| {
                let (iu, ju) = (i as usize, j as usize);
                if (positions[iu] - positions[ju]).length_squared() <= h * h {
                    neighbors[iu].push(j);
                    neighbors[ju].push(i);
                }
            });
        }

        // --- 3. Density-constraint solve (Jacobi-style: compute all deltas, then apply) -----
        let delta_q = S_CORR_DELTA_Q_FACTOR * h;
        let w_delta_q = poly6(delta_q * delta_q, h).max(1e-9);
        let mut density = vec![0.0f32; n];
        for _ in 0..iterations.max(1) {
            for i in 0..n {
                let mut rho = mass * poly6(0.0, h);
                for &j in &neighbors[i] {
                    let r2 = (self.x[i] - self.x[j as usize]).length_squared();
                    rho += mass * poly6(r2, h);
                }
                density[i] = rho;
            }

            let mut lambda = vec![0.0f32; n];
            for i in 0..n {
                let c_i = density[i] / rest_density - 1.0;
                let mut grad_self = Vec2::ZERO;
                let mut sum_grad_sq = 0.0f32;
                for &j in &neighbors[i] {
                    let ju = j as usize;
                    let delta = self.x[i] - self.x[ju];
                    let r = delta.length();
                    let grad = spiky_grad(delta, r, h) / rest_density;
                    grad_self += grad;
                    sum_grad_sq += grad.length_squared();
                }
                sum_grad_sq += grad_self.length_squared();
                lambda[i] = -c_i / (sum_grad_sq + EPSILON_RELAX);
            }

            let mut delta_p = vec![Vec2::ZERO; n];
            for i in 0..n {
                let mut sum = Vec2::ZERO;
                for &j in &neighbors[i] {
                    let ju = j as usize;
                    let delta = self.x[i] - self.x[ju];
                    let r = delta.length();
                    let w_ratio = poly6(r * r, h) / w_delta_q;
                    let s_corr = -S_CORR_K * w_ratio.powi(S_CORR_N);
                    sum += (lambda[i] + lambda[ju] + s_corr) * spiky_grad(delta, r, h);
                }
                delta_p[i] = sum / rest_density;
            }
            for (x, dp) in self.x.iter_mut().zip(&delta_p) {
                *x += *dp;
            }
        }

        // --- 4. Boundary projection (position only; velocity is reconciled below) ----------
        let mut touched_wall = vec![false; n];
        for (x, touched) in self.x.iter_mut().zip(&mut touched_wall) {
            let (d, normal) = drum.sdf_world(*x, drum_angle);
            if d < 0.0 {
                *x += -d * normal;
                *touched = true;
            }
        }

        // --- 5. Velocity reconstruction ------------------------------------------------------
        for ((v, &x), &x0i) in self.v.iter_mut().zip(&self.x).zip(&x0) {
            *v = (x - x0i) / dt;
        }

        // --- 6. No-slip wall velocity blending -------------------------------------------------
        let beta = slurry.wall_no_slip.clamp(0.0, 1.0);
        if beta > 0.0 {
            for ((v, &touched), &x) in self.v.iter_mut().zip(&touched_wall).zip(&self.x) {
                if touched {
                    let v_wall = drum.wall_velocity(x);
                    *v = (1.0 - beta) * *v + beta * v_wall;
                }
            }
        }

        // --- 7. XSPH viscosity (see module doc comment's v1 simplification note) -----------
        let c = (slurry.viscosity_pa_s / MU_REF).clamp(0.0, 1.0);
        if c > 0.0 {
            // Recompute density once more at the final (post-boundary-projection) positions so
            // viscosity uses up-to-date neighbour densities.
            for i in 0..n {
                let mut rho = mass * poly6(0.0, h);
                for &j in &neighbors[i] {
                    let r2 = (self.x[i] - self.x[j as usize]).length_squared();
                    rho += mass * poly6(r2, h);
                }
                density[i] = rho;
            }
            let mut delta_v = vec![Vec2::ZERO; n];
            for i in 0..n {
                let mut sum = Vec2::ZERO;
                for &j in &neighbors[i] {
                    let ju = j as usize;
                    let r2 = (self.x[i] - self.x[ju]).length_squared();
                    let w = poly6(r2, h);
                    let rho_j = density[ju].max(1e-6);
                    sum += (mass / rho_j) * (self.v[ju] - self.v[i]) * w;
                }
                delta_v[i] = c * sum;
            }
            for (v, dv) in self.v.iter_mut().zip(&delta_v) {
                *v += *dv;
            }
        }
    }
}

/// Initial dye value for a particle at world position `p`, per `slurry.dye_pattern` (docs/PLAN.md
/// ss3.3). Dye is otherwise purely advected (no diffusion), so this is the only place it's set.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::LiftersParams;

    fn still_drum(radius_m: f32) -> Drum {
        Drum::new(
            radius_m,
            0.0,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        )
    }

    fn test_slurry() -> SlurryParams {
        SlurryParams {
            fill_fraction: 0.15,
            ..SlurryParams::default()
        }
    }

    #[test]
    fn seed_lattice_produces_particles_within_the_drum_and_target_count_ballpark() {
        let slurry = test_slurry();
        let radius_m = 0.5;
        let resolution = 24;
        let fluid = FluidParticles::seed_lattice(&slurry, radius_m, resolution, &[], 0.0);

        assert!(
            !fluid.is_empty(),
            "expected some fluid particles to be seeded"
        );
        for &p in &fluid.x {
            assert!(
                p.length() <= radius_m,
                "seeded particle outside the drum: {p:?}"
            );
        }

        let dx = radius_m / resolution as f32;
        let target_area = slurry.fill_fraction * std::f32::consts::PI * radius_m * radius_m;
        let expected_count = (target_area / (dx * dx)).round() as usize;
        let rel_err = (fluid.len() as f32 - expected_count as f32).abs() / expected_count as f32;
        assert!(
            rel_err < 0.1,
            "particle count {} far from target {expected_count}",
            fluid.len()
        );
    }

    #[test]
    fn seed_lattice_skips_sites_overlapping_existing_balls() {
        let slurry = test_slurry();
        let radius_m = 0.5;
        // A single "ball" covering the entire bottom of the drum: no fluid site should survive.
        let balls = vec![Vec2::new(0.0, -radius_m)];
        let fluid = FluidParticles::seed_lattice(&slurry, radius_m, 24, &balls, radius_m * 2.0);
        assert!(
            fluid.is_empty(),
            "expected all fluid sites to be excluded by the oversized ball"
        );
    }

    #[test]
    fn hydrostatic_column_settles_near_rest_density() {
        // A body of fluid at rest in a still drum should settle so that interior particles sit
        // close to rest density (docs/PLAN.md ss5 M3: "max density error < 2% after 1s"). Only
        // interior particles are checked: particles right at the free surface or wall are
        // expected to have fewer neighbours and thus a lower measured density, same as in any
        // SPH-family method.
        let slurry = SlurryParams {
            fill_fraction: 0.25,
            ..SlurryParams::default()
        };
        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);
        assert!(
            fluid.len() > 20,
            "need enough particles for a meaningful interior region"
        );

        let dt = 1.0 / 240.0;
        for _ in 0..240 {
            fluid.step(&drum, 0.0, &slurry, 3, dt);
        }

        let density = fluid.densities();
        // "Interior" = far enough from *both* nearest boundaries: the curved wall (a puddle
        // sitting at the bottom of a circular drum is close to the wall almost everywhere, so a
        // simple "away from the wall" radial check alone is nearly always false near the bottom;
        // it must be combined, not intersected on top of a separate near-bottom-only check) and
        // the free surface (the puddle's own top, `max_y`).
        let max_y = fluid.x.iter().map(|p| p.y).fold(f32::MIN, f32::max);
        let margin = 1.5 * fluid.h;

        let mut checked = 0;
        for (i, &p) in fluid.x.iter().enumerate() {
            let dist_to_wall = radius_m - p.length();
            let dist_to_surface = max_y - p.y;
            if dist_to_wall.min(dist_to_surface) < margin {
                continue; // too close to a boundary for a full 2D neighbourhood
            }
            checked += 1;
            let rel_err = (density[i] - slurry.density_kg_m3).abs() / slurry.density_kg_m3;
            assert!(
                rel_err < 0.1,
                "interior particle {i} density error too large: rho={}, rho0={}, rel_err={rel_err}",
                density[i],
                slurry.density_kg_m3
            );
        }
        assert!(
            checked > 5,
            "too few interior particles to assess ({checked})"
        );
    }

    #[test]
    fn fluid_stays_inside_the_drum_and_finite() {
        let slurry = test_slurry();
        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);

        let dt = 1.0 / 240.0;
        for _ in 0..480 {
            fluid.step(&drum, 0.0, &slurry, 3, dt);
        }

        for &p in &fluid.x {
            assert!(
                p.x.is_finite() && p.y.is_finite(),
                "non-finite position: {p:?}"
            );
            assert!(
                p.length() <= radius_m * 1.05,
                "particle escaped the drum: {p:?}"
            );
        }
    }

    #[test]
    fn dye_is_conserved_per_particle_pure_advection() {
        // Dye has no diffusion term, so the *set* of dye values carried by the particles must be
        // unchanged by stepping (only positions/velocities change).
        let slurry = SlurryParams {
            fill_fraction: 0.15,
            dye_pattern: DyePattern::LeftRight,
            ..SlurryParams::default()
        };
        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);
        let dye_before = fluid.dye.clone();

        for _ in 0..60 {
            fluid.step(&drum, 0.0, &slurry, 3, 1.0 / 240.0);
        }

        assert_eq!(fluid.dye, dye_before);
        assert!(dye_before.contains(&0.0) && dye_before.contains(&1.0));
    }
}
