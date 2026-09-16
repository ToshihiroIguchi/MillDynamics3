//! Position Based Fluids (PBF) solver for the slurry.
//!
//! 2D Poly6 (density) / Spiky (gradient) kernels, iterative density-constraint projection with
//! artificial-pressure anti-clustering, drum-wall boundary projection with wall-velocity
//! blending (no-slip), and XSPH viscosity. See docs/PLAN.md ss3.3 for the full per-substep
//! algorithm and the rationale for choosing PBF over WCSPH/LBM/MPM.
//!
//! Two-way ball<->fluid coupling ([`FluidParticles::step_coupled`], docs/PLAN.md ss3.4) projects
//! fluid particles out of overlapping balls and exchanges viscous momentum near a ball's surface,
//! accumulating the reaction as a [`crate::coupling::CouplingImpulses`] for [`crate::dem`] to
//! apply. [`FluidParticles::step`] is the ball-free special case (equivalent to an empty ball list).
//!
//! **Viscosity, v1 simplification.** The XSPH mixing coefficient is a smooth saturating curve in
//! `sqrt(viscosity_pa_s)` (see [`xsph_coefficient`]), applied directly as the XSPH coefficient.
//! Like the rest of the v1 slurry model, this is qualitative (monotonic in viscosity, not a
//! quantitatively calibrated match to real Pa*s -- see [`xsph_coefficient`]'s "range caveat");
//! see the module-level caveat already noted project-wide (2D slice, qualitative results).
//!
//! An earlier version used a hard-clamped linear mapping (`clamp(viscosity_pa_s / 2.0, 0, 1)`)
//! that saturated at 2 Pa*s: every viscosity from 2 Pa*s upward produced bit-identical output,
//! so a UI control allowing up to (say) 200 Pa*s would have moved a slider that did nothing
//! past 2. [`xsph_coefficient`]'s curve never hard-saturates, so raising the UI's range stays
//! meaningful across the whole advertised span.

use glam::Vec2;

use crate::coupling::CouplingImpulses;
use crate::dem::Balls;
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
/// Viscosity (Pa*s) at which the XSPH mixing coefficient reaches half of its (asymptotic, never
/// attained) saturation value of 1.0. See [`xsph_coefficient`].
const MU_HALF_SATURATION_PA_S: f32 = 15.0;

/// XSPH mixing coefficient for a given dynamic viscosity (Pa*s): `c = sqrt(mu) / (sqrt(mu) +
/// sqrt(MU_HALF_SATURATION_PA_S))`, strictly monotonic on `[0, inf)` into `[0, 1)`. `c(0) = 0`,
/// `c(MU_HALF_SATURATION_PA_S) = 0.5`, and `c < 1` for every finite input -- so the XSPH update
/// (step 7) stays a strict convex combination for *any* non-negative viscosity, without a hard
/// clamp (a hard clamp used to be the stability guard; the curve's own asymptote now is).
///
/// Replaces an earlier `clamp(mu / 2.0, 0, 1)`, which saturated at 2 Pa*s -- see the module doc
/// comment's "Viscosity, v1 simplification" note for why that made a wide UI range meaningless.
///
/// Why `sqrt(mu)` and not `mu`: a rational (Michaelis-Menten-style) curve is required at all --
/// an exponential `1 - exp(-mu/S)` has a single scale, so keeping e.g. 100 vs. 200 Pa*s visually
/// distinguishable forces an `S` large enough to collapse the *default* slurry viscosity (0.5
/// Pa*s) down to `c ~ 0.005`, an unhelpfully weak effect. Taking the rational curve in
/// `sqrt(mu)` instead of plain `mu` widens the low end -- where real slurries mostly live -- at
/// negligible cost to the top-end separation: with `MU_HALF_SATURATION_PA_S = 15.0`,
/// `c(0.5) = 0.154` (vs. `0.032` for the same curve in plain `mu`), `c(100) = 0.870`,
/// `c(200) = 0.930` (a 6 percentage-point gap, clearly distinguishable, vs. the old scheme's
/// zero gap across the entire 2-200 Pa*s range).
///
/// **Range caveat.** XSPH's effective kinematic viscosity is `~ c * h^2 / dt`, linear in `c`
/// and therefore bounded by `h^2 / dt` regardless of the mapping: doubling the parameter from
/// 100 to 200 Pa*s changes the effective viscosity by roughly 9%, not 2x. The mapping is
/// *ordinal* across the whole advertised range (a higher number is always measurably thicker),
/// not a calibrated physical match -- consistent with this module's "qualitative, not
/// quantitatively calibrated" caveat. A genuinely wider dynamic range needs multi-pass XSPH or
/// the planned M7 implicit viscosity (docs/PLAN.md ss3.3), not a different curve; not pursued
/// here.
fn xsph_coefficient(viscosity_pa_s: f32) -> f32 {
    let s = viscosity_pa_s.max(0.0).sqrt();
    s / (s + MU_HALF_SATURATION_PA_S.sqrt())
}

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

/// 2D cross product (the z-component of the 3D cross product of `(a, 0)` and `(b, 0)`), used for
/// angular impulse = lever_arm x linear impulse.
fn cross2(a: Vec2, b: Vec2) -> f32 {
    a.x * b.y - a.y * b.x
}

/// Clamp magnitude for the coupling impulse applied to a ball over one sub-step (docs/PLAN.md
/// ss3.4: equivalent to a force of `F_CLAMP_G_MULTIPLE * m_b * g` sustained for one sub-step,
/// `F_clamp * dt`), for stability against transient large overlaps.
///
/// An earlier version of this project gave balls a 3D-sphere mass while fluid particles used a
/// unit-depth 2D mass, making a fluid particle hundreds of times heavier than a coarse-grained
/// ball; `F_CLAMP_G_MULTIPLE` was then tuned down to `3x` purely to survive that mismatch (see
/// git history). Balls are now unit-depth discs too ([`crate::dem::ball_mass`]), so the coupling
/// forces (overlap push, viscous drag, buoyancy, docs/PLAN.md ss3.4) are dimensionally consistent
/// and no longer need an artificially tight ceiling: `20x` is a pure numerical-stability backstop
/// against a transient large overlap, not a value that shapes normal behaviour. Persistent
/// clamping (see [`crate::coupling::CouplingImpulses::clamp_hits`]) after the charge has settled
/// indicates a real problem (e.g. an under-resolved fluid), not an expected steady state.
const F_CLAMP_G_MULTIPLE: f32 = 20.0;

fn impulse_clamp(ball_mass: f32, dt: f32) -> f32 {
    F_CLAMP_G_MULTIPLE * ball_mass * 9.81 * dt
}

/// Safety factor on the physically-attainable slurry speed, used by
/// [`FluidParticles::step_coupled`]'s fluid speed clamp (step 7.5): `v_max =
/// FLUID_SPEED_SAFETY_FACTOR * (|omega| * radius_m + sqrt(4 * g * radius_m))` -- the drum's
/// wall speed plus the free-fall speed across the full drum diameter. `5.0` is far above
/// anything the slurry reaches in a physically meaningful state (e.g. ~6 m/s attainable vs. a
/// ~30 m/s clamp at this project's defaults), so this is a pure stability/visibility backstop,
/// never a model parameter users tune.
///
/// Deliberately *not* a CFL-style `h / dt` bound: that is coupled to `simulation.resolution`
/// and at high resolution would sit below the physically-attainable speed, clipping real motion
/// instead of blow-ups.
const FLUID_SPEED_SAFETY_FACTOR: f32 = 5.0;

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
    /// `drum_angle`) drum, with no ball interaction. Equivalent to
    /// [`FluidParticles::step_coupled`] with an empty ball list. See docs/PLAN.md ss3.3 for the
    /// full per-substep algorithm.
    pub fn step(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        slurry: &SlurryParams,
        iterations: u32,
        dt: f32,
    ) {
        let _ = self.step_coupled(drum, drum_angle, slurry, iterations, dt, &Balls::empty());
    }

    /// As [`FluidParticles::step`], but also two-way couples against `balls` (docs/PLAN.md ss3.4,
    /// [`crate::coupling`]): fluid particles overlapping a ball are projected to its surface and
    /// fluid near a ball's surface exchanges viscous momentum with it, both accumulated into the
    /// returned [`CouplingImpulses`] (clamped) for
    /// [`crate::dem::DemState::step_with_external_forces`] to apply. `balls` is treated as a fixed
    /// boundary for this call -- it advances separately, afterwards, using the returned impulses.
    pub fn step_coupled(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        slurry: &SlurryParams,
        iterations: u32,
        dt: f32,
        balls: &Balls,
    ) -> CouplingImpulses {
        let n = self.len();
        let mut coupling = CouplingImpulses::zeros(balls.len());
        if n == 0 || dt <= 0.0 {
            return coupling;
        }

        // --- 0. Sanitize state carried in from the previous sub-step ------------------------
        // Given finite input this solver provably produces finite output (position only ever
        // advances via `x += v * dt` in step 1; every later correction is individually bounded
        // -- step 3's by the CFM-relaxed lambda denominator, step 3.5's by `0.5 * balls.radius`,
        // step 4's by the drum radius). So this can only fire on a value that entered from
        // outside this struct -- e.g. a ball position, read directly in steps 3.5/6.5, which
        // `dem.rs` does not itself guard. Resetting to the drum centre at rest is the cheapest
        // way to stop one bad particle from silently manufacturing a *finite-looking* bogus
        // ball impulse in step 3.5 (see the `is_nan` guard there) and poisoning the rendered
        // surface field (`crate::surface`).
        for (x, v) in self.x.iter_mut().zip(self.v.iter_mut()) {
            if !x.is_finite() || !v.is_finite() {
                *x = Vec2::ZERO;
                *v = Vec2::ZERO;
            }
        }

        let h = self.h;
        let rest_density = self.rest_density;
        let mass = self.particle_mass;
        let ball_grid = if balls.is_empty() {
            None
        } else {
            Some(UniformGrid::build(&balls.x, (2.0 * balls.radius).max(1e-6)))
        };

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

        // --- 3.5 Ball overlap projection (docs/PLAN.md ss3.4 step 1): push fluid particles out
        // of any ball they've penetrated, accumulating the Newton's-third-law reaction on that
        // ball. `balls` is a fixed boundary for this call (see this method's doc comment).
        if let Some(ball_grid) = &ball_grid {
            let contact_radius = balls.radius + 0.25 * h; // balls.radius + 0.5*dx (dx = h/2)
            for i in 0..n {
                let mut nearby: Vec<u32> = Vec::new();
                ball_grid.for_each_near(self.x[i], |b| nearby.push(b));
                for b in nearby {
                    let b = b as usize;
                    let delta = self.x[i] - balls.x[b];
                    let dist = delta.length();
                    if dist.is_nan() || dist >= contact_radius {
                        // `>=` alone is false for NaN, and `f32::min` below ignores NaN
                        // (returns the other operand) -- together they would turn a
                        // non-finite position into a finite, plausible, permanently-repeating
                        // reaction impulse on the ball. (Step 0 should make this unreachable
                        // in practice; this is the second layer.)
                        continue;
                    }
                    let n_hat = if dist > 1e-9 { delta / dist } else { Vec2::X };
                    // Cap the correction at half the ball's radius: when a fluid particle starts
                    // (or, at coarse fluid resolution relative to ball size, persistently ends up)
                    // deep inside `contact_radius`, projecting it out in a single sub-step would
                    // otherwise produce an outsized reaction impulse, standard practice for
                    // contact solvers (bounding the per-step correction, not just its downstream
                    // force/impulse).
                    let push_mag = (contact_radius - dist).min(0.5 * balls.radius);
                    let push = push_mag * n_hat;
                    self.x[i] += push;
                    // `push` is a position correction; the fluid particle's resulting momentum
                    // change (impulse) is `mass * (push / dt)` (velocity is reconstructed as
                    // displacement over dt, see step 5 below). Newton's third law gives the ball
                    // the exact opposite impulse -- see CouplingImpulses's doc comment for why
                    // this is expressed as an impulse, not a force (no extra `/ dt`).
                    let reaction_impulse = -(mass * push) / dt;
                    coupling.impulses[b] += reaction_impulse;
                    coupling.fluid_momentum_change -= reaction_impulse;
                    let lever = n_hat * balls.radius; // approximate contact point on the ball's surface
                    coupling.angular_impulses[b] += cross2(lever, reaction_impulse);
                }
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

        // --- 6.5 Ball viscous drag (docs/PLAN.md ss3.4 step 2): fluid within `h` of a ball's
        // surface blends toward the ball's surface velocity over this sub-step's `dt`, with the
        // blend factor derived from the disc's own viscous relaxation time -- an
        // order-of-magnitude Stokes-regime closure for drag on a disc immersed in a viscous
        // fluid, `tau = rho_ball * r^2 / (4 * mu)` (`mu` the *physical* slurry viscosity, not
        // the qualitative XSPH coefficient `c` used by step 7's fluid-fluid mixing below):
        // `beta = ball_no_slip * (1 - exp(-dt / tau))`, so a highly viscous fluid (small `tau`)
        // reaches full no-slip within one sub-step while an inviscid fluid (`mu -> 0`,
        // `tau -> inf`) applies essentially no drag. The momentum removed from the fluid is
        // added to that ball as an impulse (see CouplingImpulses's doc comment).
        let mu = slurry.viscosity_pa_s.max(0.0);
        if let Some(ball_grid) = &ball_grid {
            let beta_b = slurry.ball_no_slip.clamp(0.0, 1.0);
            let rho_ball = if balls.radius > 0.0 {
                balls.mass / (std::f32::consts::PI * balls.radius * balls.radius)
            } else {
                0.0
            };
            let factor = if mu > 1e-6 && rho_ball > 0.0 {
                let tau = (rho_ball * balls.radius * balls.radius / (4.0 * mu)).max(1e-9);
                beta_b * (1.0 - (-dt / tau).exp())
            } else {
                0.0
            };
            if factor > 0.0 {
                for i in 0..n {
                    let mut nearby: Vec<u32> = Vec::new();
                    ball_grid.for_each_near(self.x[i], |b| nearby.push(b));
                    for b in nearby {
                        let b = b as usize;
                        let r_vec = self.x[i] - balls.x[b];
                        let dist = r_vec.length();
                        if dist.is_nan() || dist >= balls.radius + h {
                            // See the identical guard in step 3.5 above for why `is_nan()` is
                            // checked explicitly rather than relying on `>=` alone.
                            continue;
                        }
                        let v_surf = balls.v[b] + balls.omega[b] * Vec2::new(-r_vec.y, r_vec.x);
                        let dv = factor * (v_surf - self.v[i]);
                        self.v[i] += dv;
                        // `dv` is already a velocity change, so the fluid's momentum change
                        // (impulse) is `mass * dv` directly -- no `/ dt` (unlike the
                        // position-derived overlap-projection impulse above).
                        let reaction_impulse = -(mass * dv);
                        coupling.impulses[b] += reaction_impulse;
                        coupling.fluid_momentum_change -= reaction_impulse;
                        coupling.angular_impulses[b] += cross2(r_vec, reaction_impulse);
                    }
                }
            }
        }

        // --- 6.6 Buoyancy (docs/PLAN.md ss3.4 step 3): balls are typically sub-resolution
        // relative to the fluid spacing and do not contribute to the density-constraint sum
        // (step 3), so buoyancy is not an emergent effect here -- it is modelled directly. Each
        // ball samples the local fluid density with the same Poly6 kernel PBF already uses
        // (reusing the fluid's own neighbour grid `grid`, built in step 2; ball positions are
        // fixed for the whole of this fluid sub-step, so querying it with a ball's position is
        // the same kind of "neighbour set from slightly-stale positions, kernel evaluated at
        // current positions" approximation step 7 already relies on) and receives an Archimedes
        // buoyant impulse `-rho_eff * (pi r^2) * g_vec * dt`, `rho_eff = min(rho_local,
        // rest_density)` -- clamping to rest density avoids over-buoyancy from a locally
        // compacted pocket, and tapers smoothly to zero as a ball approaches the free surface
        // (lower sampled density there) rather than an on/off cutoff. The reaction is applied
        // immediately as a velocity change split across the contributing fluid particles,
        // weighted by each one's share of `rho_local` -- consistent with step 6.5's direct
        // velocity-based exchange (this must run after step 5's `self.v` reconstruction, not
        // before, or the reconstruction would discard it).
        if !balls.is_empty() && balls.radius > 0.0 {
            let area = std::f32::consts::PI * balls.radius * balls.radius;
            for b in 0..balls.len() {
                if !balls.x[b].is_finite() {
                    // See the identical guard rationale in steps 3.5/6.5 above: a non-finite ball
                    // position must never manufacture a plausible-looking buoyant impulse (or, via
                    // `UniformGrid::cell_key`'s NaN-to-zero saturation, spuriously borrow real
                    // fluid neighbours near the grid origin).
                    continue;
                }
                let mut nearby: Vec<u32> = Vec::new();
                grid.for_each_near(balls.x[b], |j| nearby.push(j));
                if nearby.is_empty() {
                    continue;
                }
                let mut rho_local = 0.0f32;
                let mut weights: Vec<(usize, f32)> = Vec::with_capacity(nearby.len());
                for j in nearby {
                    let ju = j as usize;
                    let r2 = (balls.x[b] - self.x[ju]).length_squared();
                    let w = mass * poly6(r2, h);
                    if w > 0.0 {
                        rho_local += w;
                        weights.push((ju, w));
                    }
                }
                if rho_local <= 0.0 {
                    continue;
                }
                let rho_eff = rho_local.min(rest_density);
                // Acts through the ball's centroid, so it contributes no angular impulse.
                let impulse_on_ball = Vec2::new(0.0, -rho_eff * area * GRAVITY * dt);
                coupling.impulses[b] += impulse_on_ball;
                for (ju, w) in weights {
                    let frac = w / rho_local;
                    let dv = -(impulse_on_ball * frac) / mass;
                    self.v[ju] += dv;
                    coupling.fluid_momentum_change += mass * dv;
                }
            }
        }

        // --- 7. XSPH viscosity (see module doc comment's v1 simplification note) -----------
        let c = xsph_coefficient(slurry.viscosity_pa_s);
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
                let mut w_sum = 0.0f32;
                for &j in &neighbors[i] {
                    let ju = j as usize;
                    let r2 = (self.x[i] - self.x[ju]).length_squared();
                    let w = poly6(r2, h);
                    let rho_j = density[ju].max(1e-6);
                    sum += (mass / rho_j) * (self.v[ju] - self.v[i]) * w;
                    w_sum += (mass / rho_j) * w;
                }
                // Normalizing by `w_sum` (clamped to >= 1 so the well-relaxed case, where
                // `w_sum ~ 0.68`, is unaffected) keeps this a strict convex combination toward
                // the neighbour average even for the pathological density ratios where the
                // unnormalized weight sum could otherwise exceed 1 -- tightening the "usually a
                // contraction" stability argument to "always a contraction" for any `c < 1`.
                delta_v[i] = c * sum / w_sum.max(1.0);
            }
            for (v, dv) in self.v.iter_mut().zip(&delta_v) {
                *v += *dv;
            }
        }

        // --- 7.5 Fluid speed clamp (stability backstop) --------------------------------------
        // Step 5 reconstructs `v = (x - x0) / dt` from *position* corrections, so a particle
        // persistently squeezed -- e.g. trapped between a centrifuged ball layer and the wall
        // at high rotation speed -- gains `correction / dt` of speed every sub-step, with
        // nothing to stop it until step 4's wall projection caps it at roughly
        // `2 * radius_m / dt` (~240 m/s at this project's defaults) -- finite, but ~40x the
        // physically attainable speed, and a fluid field running that fast shreds the rendered
        // free surface into incoherent noise (`crate::surface`). Note: momentum removed here is
        // deliberately *not* returned to any ball (unlike step 6.5's exchange) -- a small
        // conservation violation confined to states that are already unphysical, the right
        // trade for an operation that must never itself inject energy.
        let v_max = FLUID_SPEED_SAFETY_FACTOR
            * (drum.omega.abs() * drum.radius_m + (4.0 * GRAVITY.abs() * drum.radius_m).sqrt());
        for v in self.v.iter_mut() {
            let speed = v.length();
            if !speed.is_finite() {
                *v = Vec2::ZERO;
            } else if speed > v_max {
                *v *= v_max / speed;
            }
        }

        // --- 8. Clamp accumulated coupling impulses for stability (docs/PLAN.md ss3.4) --------
        if !balls.is_empty() {
            let clamp = impulse_clamp(balls.mass, dt);
            let angular_clamp = clamp * balls.radius;
            for (impulse, angular) in coupling
                .impulses
                .iter_mut()
                .zip(&mut coupling.angular_impulses)
            {
                let mag = impulse.length();
                if mag > clamp && mag > 1e-9 {
                    *impulse *= clamp / mag;
                    coupling.clamp_hits += 1;
                }
                *angular = angular.clamp(-angular_clamp, angular_clamp);
            }
        }
        coupling
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

    #[test]
    fn viscosity_coefficient_is_strictly_monotonic_over_the_whole_ui_range() {
        let mut prev = -1.0f32;
        for &mu in &[
            0.0, 0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 30.0, 50.0, 100.0, 150.0, 200.0, 1.0e6,
        ] {
            let c = xsph_coefficient(mu);
            assert!(
                c > prev,
                "not strictly increasing at mu={mu}: {c} <= {prev}"
            );
            assert!((0.0..1.0).contains(&c), "c({mu}) = {c} outside [0, 1)");
            prev = c;
        }
    }

    #[test]
    fn viscosity_coefficient_separates_the_top_of_the_ui_range() {
        // Regression: the old `clamp(mu / 2.0, 0, 1)` mapping made every viscosity from 2 Pa*s
        // upward produce an identical coefficient of exactly 1.0, so the UI's viscosity slider
        // did nothing above 2 Pa*s.
        let (c100, c200) = (xsph_coefficient(100.0), xsph_coefficient(200.0));
        assert!(
            c200 - c100 > 0.05,
            "100 vs 200 Pa*s barely differ: {c100} -> {c200}"
        );
        assert!(
            xsph_coefficient(2.0) < 0.5,
            "2 Pa*s must no longer be near-saturated"
        );
    }

    #[test]
    fn viscosity_coefficient_keeps_the_default_slurry_meaningfully_viscous() {
        // Widening the advertised range must not make the *default* (0.5 Pa*s) look inviscid.
        let c = xsph_coefficient(SlurryParams::default().viscosity_pa_s);
        assert!(
            c > 0.1,
            "default viscosity gives a negligible XSPH coefficient: {c}"
        );
    }

    #[test]
    fn viscosity_coefficient_matches_its_documented_half_saturation_point() {
        assert!((xsph_coefficient(MU_HALF_SATURATION_PA_S) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn higher_viscosity_damps_a_high_frequency_perturbation_more_across_the_whole_range() {
        // XSPH smooths toward the local neighbourhood average, so a spatially checkerboarded
        // (the highest frequency this particle spacing can resolve) velocity perturbation
        // should be damped more by higher viscosity -- including from 100 to 200 Pa*s, which
        // the old hard-clamped mapping made bit-identical. A *linear* shear profile would not
        // distinguish any viscosity here: a symmetric averaging kernel reproduces a linear
        // field exactly, so only a genuinely high-frequency pattern exercises the coefficient.
        //
        // Settling dynamics themselves depend weakly on viscosity, so settling independently at
        // each candidate `mu` (as an earlier version of this test did) lets that per-mu
        // divergence -- not the injected signal -- dominate the measured residual, non-
        // monotonically. Instead, settle *once* at a fixed baseline viscosity to get one common
        // pre-perturbation state, then re-run the perturbation-and-decay phase from an identical
        // clone of that state at each candidate `mu` -- isolating mu's effect purely to the
        // decay phase.
        let radius_m = 0.5;
        let resolution = 20u32;
        let dx = radius_m / resolution as f32;
        let baseline_slurry = SlurryParams {
            fill_fraction: 0.2,
            viscosity_pa_s: 0.5,
            ..SlurryParams::default()
        };
        let drum = still_drum(radius_m);
        let mut baseline =
            FluidParticles::seed_lattice(&baseline_slurry, radius_m, resolution, &[], 0.0);
        for _ in 0..240 {
            baseline.step(&drum, 0.0, &baseline_slurry, 3, 1.0 / 240.0);
        }
        // Checkerboard the x-velocity by spatial column parity, so neighbours in space carry
        // opposite-sign perturbations -- the highest frequency this particle spacing resolves.
        const PERTURBATION_MPS: f32 = 2.0;
        for i in 0..baseline.len() {
            let col = (baseline.x[i].x / dx).round() as i64;
            baseline.v[i].x += if col % 2 == 0 {
                PERTURBATION_MPS
            } else {
                -PERTURBATION_MPS
            };
        }

        let residual_after_stepping = |mu: f32| -> f32 {
            let slurry = SlurryParams {
                viscosity_pa_s: mu,
                ..baseline_slurry
            };
            let mut fluid = FluidParticles {
                x: baseline.x.clone(),
                v: baseline.v.clone(),
                dye: baseline.dye.clone(),
                particle_mass: baseline.particle_mass,
                h: baseline.h,
                rest_density: baseline.rest_density,
            };
            for _ in 0..20 {
                fluid.step(&drum, 0.0, &slurry, 3, 1.0 / 240.0);
            }
            fluid.v.iter().map(|v| v.x * v.x).sum::<f32>() / fluid.len() as f32
        };

        let viscosities = [0.5f32, 2.0, 10.0, 50.0, 100.0, 200.0];
        let residuals: Vec<f32> = viscosities
            .iter()
            .map(|&mu| residual_after_stepping(mu))
            .collect();
        for w in residuals.windows(2) {
            assert!(
                w[1] < w[0],
                "viscosity stopped mattering across {viscosities:?}: residuals {residuals:?}"
            );
        }
    }

    #[test]
    fn an_injected_extreme_velocity_is_bounded_within_one_sub_step() {
        let slurry = test_slurry();
        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);
        assert!(!fluid.is_empty());
        fluid.v[0] = Vec2::new(1.0e6, 0.0);

        let dt = 1.0 / 240.0;
        fluid.step(&drum, 0.0, &slurry, 3, dt);

        let v_max = FLUID_SPEED_SAFETY_FACTOR
            * (drum.omega.abs() * drum.radius_m + (4.0 * GRAVITY.abs() * drum.radius_m).sqrt());
        for &v in &fluid.v {
            assert!(v.is_finite(), "non-finite velocity after clamping: {v:?}");
            assert!(
                v.length() <= v_max * 1.001,
                "speed {} exceeds the computed clamp {v_max}",
                v.length()
            );
        }
    }

    #[test]
    fn a_non_finite_fluid_particle_is_reset_and_never_reaches_the_ball_coupling() {
        use crate::params::EffectiveMedia;

        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let slurry = SlurryParams {
            fill_fraction: 0.0, // inject the one fluid particle manually, below
            ..SlurryParams::default()
        };
        let dt = 1.0 / 240.0;

        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);
        assert!(fluid.is_empty());
        fluid.x.push(Vec2::new(f32::NAN, 0.0));
        fluid.v.push(Vec2::ZERO);
        fluid.dye.push(0.0);

        // One ball, placed well away from the drum centre -- where the sanitize step resets a
        // non-finite particle to -- so a legitimate, unrelated overlap can't happen and any
        // impulse the ball receives must have come from the (supposedly sanitized) particle.
        let effective = EffectiveMedia {
            true_diameter_m: 0.05,
            diameter_m: 0.05,
            density_kg_m3: 7800.0,
            ball_count: 1,
            scale_factor: 1.0,
        };
        let mut dem = crate::dem::DemState::new(&effective, radius_m, 1);
        dem.balls.x[0] = Vec2::new(0.35, 0.0);
        dem.balls.v[0] = Vec2::ZERO;
        dem.balls.omega[0] = 0.0;

        let impulses = fluid.step_coupled(&drum, 0.0, &slurry, 3, dt, &dem.balls);

        assert!(
            fluid.x[0].is_finite() && fluid.v[0].is_finite(),
            "non-finite fluid particle survived the sanitize step: x={:?} v={:?}",
            fluid.x[0],
            fluid.v[0]
        );
        assert!(fluid.x[0].length() <= radius_m * 1.05);
        assert!(
            impulses.impulses[0].length() < 1e-6,
            "a non-finite fluid particle leaked a reaction impulse into the ball solver: {:?}",
            impulses.impulses[0]
        );
    }

    #[test]
    fn a_non_finite_ball_position_does_not_contaminate_the_fluid() {
        use crate::params::EffectiveMedia;

        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let slurry = test_slurry();
        let dt = 1.0 / 240.0;

        let effective = EffectiveMedia {
            true_diameter_m: 0.05,
            diameter_m: 0.05,
            density_kg_m3: 7800.0,
            ball_count: 4,
            scale_factor: 1.0,
        };
        let mut dem = crate::dem::DemState::new(&effective, radius_m, 1);
        dem.balls.x[0] = Vec2::splat(f32::INFINITY);
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);

        let _ = fluid.step_coupled(&drum, 0.0, &slurry, 3, dt, &dem.balls);

        for (&x, &v) in fluid.x.iter().zip(&fluid.v) {
            assert!(
                x.is_finite() && v.is_finite(),
                "fluid contaminated by a non-finite ball position: x={x:?} v={v:?}"
            );
        }
    }
}
