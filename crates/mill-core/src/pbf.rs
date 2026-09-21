//! Position Based Fluids (PBF) solver for the slurry.
//!
//! 2D Poly6 (density) / Spiky (gradient) kernels, iterative one-sided density-constraint
//! projection, drum-wall boundary projection with wall-velocity blending (no-slip), and implicit
//! Newtonian viscosity. See docs/PLAN.md ss3.3 for the full per-substep algorithm and the
//! rationale for choosing PBF over WCSPH/LBM/MPM.
//!
//! Two-way ball<->fluid coupling ([`FluidParticles::step_coupled`], docs/PLAN.md ss3.4) pushes
//! fluid particles out of overlapping balls, exchanges viscous momentum near a ball's surface, and
//! (step 3.6) applies Akinci-style fluid-fluid cohesion and fluid-boundary adhesion -- the only
//! attraction anywhere in this solver, giving the slurry a genuine (if empirically calibrated)
//! contact angle instead of the density constraint's incidental unclamped-attraction side effect
//! an earlier version relied on. Every ball-side reaction accumulates into a
//! [`crate::coupling::CouplingImpulses`] for [`crate::dem`] to apply. [`FluidParticles::step`] is
//! the ball-free special case (equivalent to an empty ball list).
//!
//! **Viscosity.** Step 7 discretises the Newtonian viscous term `(mu/rho) * laplacian(v)` with the
//! Morris (1997) SPH viscous Laplacian ([`morris_weights`]) and solves the resulting backward-Euler
//! system `(I + dt*L) v_new = v` implicitly by conjugate gradient ([`solve_implicit_viscosity`]),
//! `L` being the (symmetric positive semi-definite) graph Laplacian `(Lv)_i = sum_j c_ij (v_i -
//! v_j)`. Unlike an explicit scheme, this is unconditionally stable for any `viscosity_pa_s` at
//! this project's fixed sub-step, and `viscosity_pa_s` enters the physics directly (Pa*s) rather
//! than through a separately-calibrated qualitative coefficient -- an earlier version used
//! explicit XSPH mixing with a hand-tuned saturating coefficient curve instead (see git history),
//! whose effective kinematic viscosity was bounded by `~h^2/dt` regardless of the input and so
//! could not represent a genuinely low-viscosity (near-water) slurry. This is still a 2D-slice
//! model (see the project-wide qualitative-results caveat), and `nu = mu/rho` uses the fluid's own
//! *measured* SPH density at each particle, not `slurry.density_kg_m3` directly -- the two agree
//! closely in a well-relaxed interior but diverge at a low-density boundary (free surface, sparse
//! neighbourhood), where the resulting locally elevated `nu` is itself physically reasonable (a
//! thin/rarefied film should be *more*, not less, resistant to shearing relative to its own mass).

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

/// Step 3's density-constraint solve treats `simulation.pbf_iterations` as a floor, not a fixed
/// count: once that many passes have run, it keeps going -- up to this many further passes -- only
/// while the worst-particle compression error (`c_i`, see step 3's loop) still exceeds
/// [`COMPRESSION_CONVERGENCE_TOLERANCE`]. A violent lifter/cataracting impact can inject enough
/// velocity in one sub-step that the fixed default (3) under-converges, producing the compression-
/// error spikes documented in docs/PHYSICS.md; escalating only when actually needed keeps the
/// common (already-converged) case at its prior cost.
const ADAPTIVE_MAX_EXTRA_ITERATIONS: u32 = 6;
/// Compression-error convergence target for step 3's escalation, matching the bound this project's
/// own quiescent-state regression test already asserts (`metrics.rs`,
/// `settled_puddle_compression_error_is_small`).
const COMPRESSION_CONVERGENCE_TOLERANCE: f32 = 0.02;

/// Step 3.5's overlap-push standoff, as a multiple of `h` beyond the ball's own radius
/// (`contact_radius = balls.radius + OVERLAP_PUSH_MARGIN * h`). A small numerical margin only
/// (previously `0.25`, which at this project's defaults forced roughly half a particle spacing of
/// permanently fluid-free gap around every ball -- a geometric 180 degree contact angle baked into
/// the discretisation, independent of any wetting parameter). See step 3.6 for the cohesion/
/// adhesion terms that give this project's slurry an actual (if empirically calibrated) contact
/// angle instead.
const OVERLAP_PUSH_MARGIN: f32 = 0.05;

/// Reference surface tension (real water, N/m) at which `COHESION_ACCEL_FACTOR * |GRAVITY|` of
/// peak fluid-fluid cohesive pull applies; `slurry.surface_tension_n_m` scales linearly from
/// there. Like the wider SPH surface-tension literature -- Akinci, Akinci & Teschner ("Versatile
/// Surface Tension and Adhesion for SPH Fluids", 2013) does not give a closed-form mapping from
/// its own cohesion coefficient to a macroscopic N/m value either, since the mapping is inherently
/// resolution-dependent (a coarser or finer fluid lattice needs a different coefficient for the
/// same emergent surface tension) -- this is a calibrated reference point at this project's
/// current default resolution, not a first-principles unit conversion.
const REFERENCE_SURFACE_TENSION_N_M: f32 = 0.072;
/// Peak fluid-fluid cohesion acceleration, as a multiple of `|GRAVITY|`, at
/// [`REFERENCE_SURFACE_TENSION_N_M`]. Same role as [`ADHESION_ACCEL_FACTOR`] for the fluid-
/// boundary term: comparable to (not dominant over) gravity, so real cohesion is visible without
/// becoming the largest term in the coupling budget.
const COHESION_ACCEL_FACTOR: f32 = 2.0;
/// Peak fluid-boundary adhesion acceleration, as a multiple of `|GRAVITY|`, at
/// `slurry.wettability = 1`. Chosen to be comparable to (not dominant over) gravity/buoyancy, so a
/// fully-wetting slurry visibly clings to media without the term becoming the largest force in the
/// coupling budget; bounded further by the existing per-ball impulse clamp (step 8,
/// `F_CLAMP_G_MULTIPLE`) alongside every other coupling mechanism.
const ADHESION_ACCEL_FACTOR: f32 = 2.0;

/// Relative-residual stopping tolerance for [`solve_implicit_viscosity`]'s conjugate-gradient
/// solve: iterate until `|r| <= VISCOSITY_CG_TOLERANCE * |b|`, `b` the right-hand side (the
/// incoming velocity field).
const VISCOSITY_CG_TOLERANCE: f32 = 1e-3;
/// Iteration budget for [`solve_implicit_viscosity`]. The operator is well-conditioned (a
/// graph Laplacian scaled by `dt`, itself small at this project's fixed sub-step), so this is a
/// generous ceiling in practice, not a value tuned tight against real convergence needs.
const VISCOSITY_CG_MAX_ITERS: u32 = 50;
/// Regularisation `eta^2 = VISCOSITY_ETA_FACTOR * h^2` added to the denominator of
/// [`morris_weights`], preventing a blow-up as `|x_ij| -> 0` (docs/PLAN.md ss3.3's own
/// `EPSILON_RELAX` plays the analogous role for the density-constraint solve). `0.01` is the
/// standard choice in the SPH viscosity literature (Morris 1997, Monaghan 1992's artificial
/// viscosity).
const VISCOSITY_ETA_FACTOR: f32 = 0.01;

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

/// Akinci-style cohesion kernel (Akinci, Akinci & Teschner 2013, eq. 2), adapted to 2D. Used for
/// fluid-fluid cohesion directly, and (restricted to its outer half, see [`adhesion_kernel`]) as
/// the shape fluid-boundary adhesion reuses. Vanishes at `r = 0` and `r = h`, peaks near
/// `r ~= 0.6h`; unlike a density kernel ([`poly6`]) it goes slightly *negative* for `r < h/2` (the
/// paper's own construction), so two particles very close together get a small locally-repulsive
/// contribution instead of being pulled all the way to coincidence -- the force itself is always
/// attractive-by-construction (applied as `-strength * kernel(r) * direction_toward_neighbour`),
/// so a negative kernel value here means a locally repulsive contribution to that same force.
///
/// Normalised for **2D** (`2*pi * integral_0^h kernel(r) r dr == 1`, matching this module's
/// [`poly6`]/[`spiky_grad`] convention of integrating to exactly 1 over the disc of radius `h`) by
/// direct integration -- the published 3D coefficient (`32/(pi h^9)`) does not carry over to 2D by
/// simply dropping a power of `h`, the same way every other kernel in this module is independently
/// normalised for 2D rather than reusing its 3D literature constant.
fn cohesion_kernel(r: f32, h: f32) -> f32 {
    if r <= 0.0 || r >= h {
        return 0.0;
    }
    let coeff = 35840.0 / (209.0 * std::f32::consts::PI * h.powi(8));
    let term = (h - r).powi(3) * r.powi(3);
    if r <= 0.5 * h {
        coeff * (2.0 * term - h.powi(6) / 64.0)
    } else {
        coeff * term
    }
}

/// Fluid-boundary adhesion kernel: the same cubic shape as [`cohesion_kernel`]'s outer half
/// (`r` in `(h/2, h)`), independently re-normalised over that restricted domain
/// (`2*pi * integral_{h/2}^h kernel(r) r dr == 1`) rather than reusing the original paper's
/// fractional-power adhesion formula (`0.007 h^-3.25 * (-4r^2/h + 6r - 2h)^0.25`) -- that
/// formula's non-integer exponents have no clean closed-form 2D re-derivation, and this kernel
/// needs no property the simpler cubic shape (zero at both ends of `(h/2, h)`, single smooth hump,
/// exact to normalise) lacks for this solver's purpose: pulling fluid onto a ball's surface from
/// moderate range, tapering smoothly to zero right at contact (already the overlap push's job,
/// step 3.5) and at the outer cutoff.
fn adhesion_kernel(r: f32, h: f32) -> f32 {
    if r <= 0.5 * h || r >= h {
        return 0.0;
    }
    let coeff = 35840.0 / (163.0 * std::f32::consts::PI * h.powi(8));
    coeff * (h - r).powi(3) * r.powi(3)
}

/// Boundary-particle positions for Akinci-style fluid-boundary adhesion: each ball's
/// circumference sampled at roughly `dx` (the fluid lattice spacing) arc-length spacing, at least
/// 6 points regardless of how small the ball is. Recomputed fresh each sub-step from the balls'
/// *current* centres (cheap: a fixed small point count per ball, no allocation carried across
/// sub-steps) rather than cached, since balls move every sub-step. Points for ball `b` occupy a
/// contiguous run of the returned `Vec` (ball order, then angle order), so a caller that already
/// knows the per-ball point count can slice directly instead of filtering by ball index.
fn ball_boundary_particles(balls: &Balls, dx: f32) -> Vec<Vec2> {
    if balls.is_empty() || balls.radius <= 0.0 || dx <= 0.0 {
        return Vec::new();
    }
    let count_per_ball = boundary_points_per_ball(balls.radius, dx);
    let mut position = Vec::with_capacity(balls.len() * count_per_ball as usize);
    for b in 0..balls.len() {
        for k in 0..count_per_ball {
            let angle = std::f32::consts::TAU * (k as f32) / (count_per_ball as f32);
            let offset = Vec2::new(angle.cos(), angle.sin()) * balls.radius;
            position.push(balls.x[b] + offset);
        }
    }
    position
}

/// Number of boundary points [`ball_boundary_particles`] places on one ball's circumference:
/// roughly `dx` arc-length spacing, floored at 6 so even a very small (relative to `dx`) ball
/// gets a reasonable angular sampling.
fn boundary_points_per_ball(ball_radius: f32, dx: f32) -> u32 {
    let circumference = std::f32::consts::TAU * ball_radius;
    ((circumference / dx).round() as u32).max(6)
}

/// 2D cross product (the z-component of the 3D cross product of `(a, 0)` and `(b, 0)`), used for
/// angular impulse = lever_arm x linear impulse.
fn cross2(a: Vec2, b: Vec2) -> f32 {
    a.x * b.y - a.y * b.x
}

/// Builds each particle's neighbour-index list (mutual: `j` in `neighbors[i]` implies `i` in
/// `neighbors[j]`) from `grid` (a [`UniformGrid`] already built over `x` with cell size `h`),
/// keeping only candidate pairs actually within `h`. Factored out of [`FluidParticles::step_coupled`]
/// (its step 2) so the viscosity solver's tests below can drive it in isolation without a full
/// coupled sub-step.
fn build_neighbor_lists(x: &[Vec2], grid: &UniformGrid, h: f32) -> Vec<Vec<u32>> {
    let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); x.len()];
    grid.for_each_candidate_pair(|i, j| {
        let (iu, ju) = (i as usize, j as usize);
        if (x[iu] - x[ju]).length_squared() <= h * h {
            neighbors[iu].push(j);
            neighbors[ju].push(i);
        }
    });
    neighbors
}

/// Real inner product `<a, b> = sum_i a_i . b_i` treating a `&[Vec2]` field as one flat `2n`-real
/// vector, used throughout [`solve_implicit_viscosity`]'s conjugate-gradient iteration.
fn dot(a: &[Vec2], b: &[Vec2]) -> f32 {
    a.iter().zip(b).map(|(&ai, &bi)| ai.dot(bi)).sum()
}

/// Morris (1997) SPH viscous-Laplacian weights, laid out parallel to `neighbors` (`out[i][k]` is
/// the weight of `neighbors[i][k]`):
///
/// ```text
/// c_ij = m_j * (mu_i + mu_j) / (rho_i * rho_j)
///        * |x_ij . grad_W_ij| / (|x_ij|^2 + eta^2), eta^2 = VISCOSITY_ETA_FACTOR * h^2
/// ```
///
/// so that `(Lv)_i = sum_j c_ij * (v_i - v_j)` (see [`laplacian_apply`]) discretises `-(mu/rho) *
/// laplacian(v)` at particle `i`. This crate's fluid is single-phase (`mu_i = mu_j = mu` for every
/// particle), so the `(mu_i + mu_j)` term is simply `2*mu` here, but the two-sided form is kept
/// for clarity against the cited derivation. `c_ij >= 0` for every pair (each factor is
/// non-negative), so the resulting `L` is a graph Laplacian: symmetric (`c_ij` depends only on the
/// unordered pair) and positive semi-definite -- the property [`solve_implicit_viscosity`]'s
/// conjugate-gradient solve relies on.
fn morris_weights(
    x: &[Vec2],
    neighbors: &[Vec<u32>],
    density: &[f32],
    mass: f32,
    mu: f32,
    h: f32,
) -> Vec<Vec<f32>> {
    let eta2 = VISCOSITY_ETA_FACTOR * h * h;
    neighbors
        .iter()
        .enumerate()
        .map(|(i, js)| {
            let rho_i = density[i].max(1e-6);
            js.iter()
                .map(|&j| {
                    let ju = j as usize;
                    let rho_j = density[ju].max(1e-6);
                    let delta = x[i] - x[ju];
                    let r = delta.length();
                    let grad = spiky_grad(delta, r, h);
                    let numerator = (mass * 2.0 * mu / (rho_i * rho_j)) * delta.dot(grad).abs();
                    numerator / (delta.length_squared() + eta2)
                })
                .collect()
        })
        .collect()
}

/// Matrix-free application of the graph Laplacian `(Lv)_i = sum_j c_ij * (v_i - v_j)`, `weights`
/// laid out parallel to `neighbors` (see [`morris_weights`]).
fn laplacian_apply(neighbors: &[Vec<u32>], weights: &[Vec<f32>], v: &[Vec2]) -> Vec<Vec2> {
    neighbors
        .iter()
        .zip(weights)
        .enumerate()
        .map(|(i, (js, cs))| {
            let mut sum = Vec2::ZERO;
            for (&j, &c) in js.iter().zip(cs) {
                sum += c * (v[i] - v[j as usize]);
            }
            sum
        })
        .collect()
}

/// Applies the implicit-viscosity system operator `A = I + dt*L` to `v`.
fn apply_viscosity_system(
    neighbors: &[Vec<u32>],
    weights: &[Vec<f32>],
    dt: f32,
    v: &[Vec2],
) -> Vec<Vec2> {
    let lv = laplacian_apply(neighbors, weights, v);
    v.iter().zip(&lv).map(|(&vi, &lvi)| vi + dt * lvi).collect()
}

/// Solves `(I + dt*L) v_new = v` in place by conjugate gradient, warm-started from `v` itself
/// (the incoming velocity is usually already close to the diffused result at this project's
/// small fixed sub-step), and returns the iteration count actually used (0 if the right-hand
/// side was already within tolerance, or degenerate). `A = I + dt*L` is symmetric positive
/// definite whenever `L` is symmetric positive semi-definite (true of [`morris_weights`]'s
/// output, docs/PLAN.md ss3.3) and `dt > 0`, since `L`'s eigenvalues are `>= 0` and `I` shifts
/// them to `>= 1`.
///
/// The 2D vector system is solved directly over `Vec2`, with the real inner product `<a, b> =
/// sum_i a_i . b_i` ([`dot`]), rather than as two independent scalar systems. This is ordinary CG
/// on the full `2n`-dimensional real vector space `A` and `v` actually live in (`L`, and hence
/// `A`, acts identically and independently on the x- and y-components, so this is mathematically
/// equivalent to running the scalar algorithm twice -- just without the bookkeeping of splitting
/// and re-merging the two components).
fn solve_implicit_viscosity(
    neighbors: &[Vec<u32>],
    weights: &[Vec<f32>],
    v: &mut [Vec2],
    dt: f32,
) -> u32 {
    // The right-hand side `b` *is* the incoming velocity, and `x` is warm-started to the same
    // value, so neither needs its own copy: `v` plays both roles and is updated in place.
    let b_norm = dot(v, v).sqrt();
    if b_norm <= 0.0 || !b_norm.is_finite() {
        // All-zero (nothing to diffuse) or non-finite (step 0 should have made this unreachable;
        // running CG on it would only spread the NaN across the whole field).
        return 0;
    }
    let target = VISCOSITY_CG_TOLERANCE * b_norm;

    let ax = apply_viscosity_system(neighbors, weights, dt, v);
    let mut r: Vec<Vec2> = v.iter().zip(&ax).map(|(&bi, &axi)| bi - axi).collect();
    let mut p = r.clone();
    let mut rs_old = dot(&r, &r);
    if rs_old.sqrt() <= target {
        return 0;
    }

    for iter in 1..=VISCOSITY_CG_MAX_ITERS {
        let ap = apply_viscosity_system(neighbors, weights, dt, &p);
        let p_ap = dot(&p, &ap);
        if p_ap <= 0.0 || !p_ap.is_finite() {
            // Unreachable for a genuinely SPD operator; bail out rather than divide by ~0 if f32
            // rounding on a degenerate (e.g. fully disconnected) neighbourhood ever produces it.
            return iter;
        }
        let alpha = rs_old / p_ap;
        for i in 0..v.len() {
            v[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rs_new = dot(&r, &r);
        if rs_new.sqrt() <= target {
            return iter;
        }
        let beta = rs_new / rs_old;
        for i in 0..v.len() {
            p[i] = r[i] + beta * p[i];
        }
        rs_old = rs_new;
    }
    VISCOSITY_CG_MAX_ITERS
}

/// Mean shear rate `gamma_dot = sqrt(2 D:D)` over all particles, `D = sym(grad v)` the symmetric
/// part of the SPH velocity-gradient tensor `grad_v_i = sum_j (m_j/rho_j) * (v_j - v_i) (x)
/// grad_W_ij` (standard SPH gradient estimator, the same neighbour-averaging form the density
/// constraint (step 3) and viscosity ([`morris_weights`]) already use). Used for
/// [`crate::metrics`]'s mean-shear-rate readout and (future work) a Bingham/Herschel-Bulkley
/// effective viscosity. `0.0` if there are no particles or every particle's gradient is
/// degenerate (no neighbours).
fn mean_shear_rate(
    x: &[Vec2],
    v: &[Vec2],
    neighbors: &[Vec<u32>],
    density: &[f32],
    mass: f32,
    h: f32,
) -> f32 {
    let n = x.len();
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0f32;
    let mut count = 0u32;
    for i in 0..n {
        let (mut gxx, mut gxy, mut gyx, mut gyy) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
        for &j in &neighbors[i] {
            let ju = j as usize;
            let delta = x[i] - x[ju];
            let r = delta.length();
            let grad = spiky_grad(delta, r, h);
            let rho_j = density[ju].max(1e-6);
            let coeff = mass / rho_j;
            let dv = v[ju] - v[i];
            gxx += coeff * dv.x * grad.x;
            gxy += coeff * dv.x * grad.y;
            gyx += coeff * dv.y * grad.x;
            gyy += coeff * dv.y * grad.y;
        }
        let dxy = 0.5 * (gxy + gyx);
        let d_contraction = gxx * gxx + gyy * gyy + 2.0 * dxy * dxy;
        let gamma_dot = (2.0 * d_contraction).sqrt();
        if gamma_dot.is_finite() {
            sum += gamma_dot;
            count += 1;
        }
    }
    if count == 0 {
        0.0
    } else {
        sum / count as f32
    }
}

/// Per-sub-step fluid solver diagnostics, returned alongside (or instead of, for
/// [`FluidParticles::step`]) [`crate::coupling::CouplingImpulses`]. Surfaced by
/// [`crate::Simulation`] for [`crate::metrics`] (docs/PLAN.md ss3.5).
#[derive(Debug, Clone, Copy, Default)]
pub struct FluidStepStats {
    /// Conjugate-gradient iterations [`solve_implicit_viscosity`] used this sub-step. A healthy
    /// run stays comfortably under [`VISCOSITY_CG_MAX_ITERS`]; persistently hitting the ceiling
    /// would mean the viscosity solve is not actually converging.
    pub viscosity_iterations: u32,
    /// Mean shear rate (1/s) over the fluid population this sub-step, from [`mean_shear_rate`].
    pub mean_shear_rate_per_s: f32,
    /// Copy of this sub-step's [`crate::coupling::CouplingImpulses::clamp_hits`], surfaced here
    /// too (alongside the other solver diagnostics `step_coupled` already returns) so
    /// [`crate::Simulation`] does not need a third return value from `step_coupled` just for one
    /// `u32`.
    pub coupling_clamp_hits: u32,
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
/// and no longer need an artificially tight ceiling.
///
/// `60x` (raised from an intermediate `20x`) accounts for a real difference in how the two
/// coupling mechanisms this clamp bounds scale. Step 3.5's overlap-push impulse is now itself
/// momentum-bounded (see that step's doc comment: it exchanges only the momentum needed to
/// arrest a particle's actual approach velocity, not `mass * push / dt` for the *position*
/// correction's full, dt-scaled magnitude as an earlier version did) and so rarely approaches
/// this clamp at all any more. Step 6.5's drag correction (see that step's own doc comment) is a
/// *velocity*-matching relaxation, `impulse = mass * a * (v_bar - v_b)` with `a <= 1` bounded but
/// otherwise `dt`-independent once the relaxation factor `beta` saturates (`dt >= tau`, i.e.
/// `viscosity_pa_s` above roughly 10 Pa*s at this project's ball sizes) -- physically correct
/// (achieving no-slip within one sub-step is what "highly viscous" means, and can legitimately
/// need a large single-sub-step velocity change when a small ball is embedded in a much larger
/// entrained fluid mass moving fast relative to it), so this clamp still needs to stay generous
/// enough for that mechanism specifically, not just for step 3.5's (now much smaller) needs.
/// `60x` keeps clamping rare (single-digit percent of ball-substeps, see
/// `coupling::tests::cascading_charge_keeps_coupling_clamp_hits_rare_once_settled`) at
/// `viscosity_pa_s` up to 200 while still catching genuine numerical blow-ups; persistent
/// clamping after the charge has settled still indicates a real problem, not an expected
/// steady state.
const F_CLAMP_G_MULTIPLE: f32 = 60.0;

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
        // Area each particle represents on the hexagonal lattice this function seeds: columns are
        // `dx` apart, rows `row_h = dx*sqrt(3)/2` apart, so the per-site cell is `dx * row_h`, and
        // the lattice's number density is `1/(dx*row_h) = 1.1547/dx^2`. Deriving both the mass and
        // the site count from that same cell area keeps three things consistent at once: total mass
        // is `rest_density * target_area`, the seeded footprint *is* `target_area` (i.e. the
        // `fill_fraction` the user asked for), and the Poly6 sum over the seeded lattice already
        // reads `rest_density` to within 0.1% -- so the run no longer starts with a compression
        // transient. An earlier version used the square-lattice cell `dx*dx` for both while seeding
        // hexagonally, which left every run ~15% over-dense at `t = 0`.
        let row_h = dx * (3.0f32.sqrt() * 0.5);
        let cell_area = dx * row_h;
        let particle_mass = rest_density * cell_area;

        let target_area =
            slurry.fill_fraction * std::f32::consts::PI * drum_radius_m * drum_radius_m;
        let target_count = (target_area / cell_area).round().max(0.0) as usize;

        let mut x = Vec::with_capacity(target_count);
        let mut dye = Vec::with_capacity(target_count);
        if target_count > 0 && dx > 0.0 {
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
    ) -> FluidStepStats {
        self.step_coupled(drum, drum_angle, slurry, iterations, dt, &Balls::empty())
            .1
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
    ) -> (CouplingImpulses, FluidStepStats) {
        let n = self.len();
        let mut coupling = CouplingImpulses::zeros(balls.len());
        if n == 0 || dt <= 0.0 {
            return (coupling, FluidStepStats::default());
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
        // Step 3.5's own query radius (`contact_radius`); step 3.6's cohesion/adhesion terms
        // build and query their own grids separately (fluid-fluid neighbours already come from
        // step 2's `grid`/`neighbors`, and fluid-boundary adhesion queries that same `grid` from
        // each boundary point, see step 3.6 below).
        let ball_grid = if balls.is_empty() {
            None
        } else {
            let contact_radius = balls.radius + OVERLAP_PUSH_MARGIN * h;
            Some(UniformGrid::build(
                &balls.x,
                (2.0 * balls.radius).max(contact_radius).max(1e-6),
            ))
        };

        // --- 1. Predict --------------------------------------------------------------------
        let x0: Vec<Vec2> = self.x.clone();
        for i in 0..n {
            self.v[i].y += GRAVITY * dt;
            self.x[i] += self.v[i] * dt;
        }

        // --- 2. Neighbour lists (rebuilt each sub-step, reused across iterations) -----------
        let grid = UniformGrid::build(&self.x, h);
        let neighbors = build_neighbor_lists(&self.x, &grid, h);

        // --- 3. Density-constraint solve (Jacobi-style: compute all deltas, then apply) -----
        // `iterations` is a floor here, not a fixed count -- see [`ADAPTIVE_MAX_EXTRA_ITERATIONS`].
        let mut density = vec![0.0f32; n];
        let base_iterations = iterations.max(1);
        let mut iter_count = 0u32;
        loop {
            for i in 0..n {
                let mut rho = mass * poly6(0.0, h);
                for &j in &neighbors[i] {
                    let r2 = (self.x[i] - self.x[j as usize]).length_squared();
                    rho += mass * poly6(r2, h);
                }
                density[i] = rho;
            }

            let mut lambda = vec![0.0f32; n];
            let mut max_c = 0.0f32;
            for i in 0..n {
                // Pressure-only constraint: `C_i = max(0, rho_i/rest_density - 1)`, so only
                // *over*-dense particles are projected apart and a density-deficient one is left
                // alone. Macklin & Muller's original PBF leaves `C_i` unclamped and relies on an
                // artificial-pressure term (`s_corr`) to suppress the resulting tensile
                // instability; an earlier version of this crate had the unclamped form with that
                // remedy disabled (`k = 0`, never found stable at this project's SI parameter
                // scale), which is an inconsistent combination. Unclamped, a deficient particle
                // gets `lambda_i > 0`, and since `spiky_grad`'s coefficient is negative its
                // `delta_p` points *toward* its neighbours -- an attraction strongest exactly
                // where the kernel is truncated (free surface, thin films, isolated splash),
                // which is what balled slurry up into mid-air blobs before this constraint was
                // clamped. This is not this solver's fluid-fluid attraction mechanism any more
                // than an unclamped constraint's incidental side effect ever should have been --
                // step 3.6's cohesion term (`slurry.surface_tension_n_m`) is the real one, and is
                // explicit, bounded, and tunable independently of the incompressibility solve.
                let c_i = (density[i] / rest_density - 1.0).max(0.0);
                max_c = max_c.max(c_i);
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
                    sum += (lambda[i] + lambda[ju]) * spiky_grad(delta, r, h);
                }
                delta_p[i] = sum / rest_density;
            }
            for (x, dp) in self.x.iter_mut().zip(&delta_p) {
                *x += *dp;
            }

            iter_count += 1;
            let converged = max_c <= COMPRESSION_CONVERGENCE_TOLERANCE;
            if (iter_count >= base_iterations && converged)
                || iter_count >= base_iterations + ADAPTIVE_MAX_EXTRA_ITERATIONS
            {
                break;
            }
        }

        // --- 3.5 Ball overlap projection (docs/PLAN.md ss3.4 step 1): push fluid particles out
        // of any ball they've penetrated, accumulating the Newton's-third-law reaction on that
        // ball. `balls` is a fixed boundary for this call (see this method's doc comment).
        //
        // The position correction (`push`) must satisfy the geometric constraint in full, but the
        // *momentum* it exchanges with the ball is bounded separately, by how much of the
        // particle's actual velocity toward the ball's surface this contact needs to arrest --
        // not by how far the position moved. Without this split, the implicit velocity
        // `push_mag / dt` (step 5 reconstructs velocity from position) diverges as `dt` shrinks or
        // as fluid particles pile up against a ball, independent of any real contact speed -- this
        // was the coupling half of the fluidised-charge energy-injection bug (docs/PHYSICS.md):
        // a purely resting/settled overlap (particle not actually approaching the ball) must
        // exchange no momentum at all, only a position correction. Any part of the push beyond the
        // physically-justified impulse is tracked in `push_velocity_excess` and subtracted back
        // out of the fluid's own velocity after step 5's position-to-velocity reconstruction.
        let mut push_velocity_excess = vec![Vec2::ZERO; n];
        if let Some(ball_grid) = &ball_grid {
            let contact_radius = balls.radius + OVERLAP_PUSH_MARGIN * h;
            // `i` indexes both `self.x`/`self.v` (mutated in place, alongside `push_velocity_excess`)
            // and is looked up via the grid, so a plain iterator/enumerate over one collection
            // doesn't fit cleanly here.
            #[allow(clippy::needless_range_loop)]
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
                    // otherwise produce an outsized *position* correction. This bounds the
                    // geometric correction only; see above for the separate momentum bound.
                    let push_mag = (contact_radius - dist).min(0.5 * balls.radius);
                    let push = push_mag * n_hat;
                    self.x[i] += push;

                    let lever = n_hat * balls.radius; // approximate contact point on the ball's surface
                    let v_surface = balls.v[b] + balls.omega[b] * Vec2::new(-lever.y, lever.x);
                    // How much of the particle's velocity toward the ball this contact needs to
                    // arrest (>= 0; a particle already moving away, or at rest relative to the
                    // ball, needs none).
                    let arrest_speed = (-(self.v[i] - v_surface).dot(n_hat)).max(0.0);
                    let full_impulse_mag = mass * push_mag / dt;
                    let impulse_mag = full_impulse_mag.min(mass * arrest_speed);
                    push_velocity_excess[i] += ((full_impulse_mag - impulse_mag) / mass) * n_hat;

                    // `impulse_mag` is the fluid particle's actual momentum change (impulse) from
                    // this contact; Newton's third law gives the ball the exact opposite -- see
                    // `CouplingImpulses`'s doc comment for why this is expressed as an impulse,
                    // not a force (no extra `/ dt`).
                    let reaction_impulse = -impulse_mag * n_hat;
                    coupling.impulses[b] += reaction_impulse;
                    coupling.fluid_momentum_change += impulse_mag * n_hat;
                    coupling.angular_impulses[b] += cross2(lever, reaction_impulse);
                }
            }
        }

        // --- 3.6 Fluid-fluid cohesion and fluid-boundary adhesion (Akinci, Akinci & Teschner,
        // "Versatile Surface Tension and Adhesion for SPH Fluids", 2013, adapted to 2D and to this
        // solver's impulse-based ball coupling) -- the only attraction anywhere in this solver.
        // Together these give the slurry a genuine contact angle: cohesion (3.6a) holds a fluid
        // film together against gravity, adhesion (3.6b) pulls it onto media, and which one wins
        // where is an emergent competition between two symmetric pairwise forces, not a single
        // ball<->fluid-only knob with no fluid-fluid counterpart (docs/PHYSICS.md ss6.1a's earlier
        // "deliberately not a true Young's-equation contact angle" limitation).

        // --- 3.6a Fluid-fluid cohesion (`slurry.surface_tension_n_m`). A fully interior particle
        // has a (near-)symmetric neighbourhood, so the pairwise pulls cancel and it feels no net
        // force; only a free-surface/thin-film particle, whose neighbourhood is asymmetric, feels
        // a net inward pull -- ordinary surface tension. Reuses `neighbors` (step 2) directly, no
        // new broad-phase pass. Purely a velocity effect (`accel * dt`, like gravity/buoyancy), so
        // it must be deferred past step 5's position-driven velocity reconstruction exactly like
        // step 3.5's `push_velocity_excess` -- applying it to `self.v` here would be silently
        // discarded.
        let mut cohesion_velocity_delta = vec![Vec2::ZERO; n];
        let cohesion_accel = (slurry.surface_tension_n_m.max(0.0) / REFERENCE_SURFACE_TENSION_N_M)
            * COHESION_ACCEL_FACTOR
            * GRAVITY.abs();
        if cohesion_accel > 0.0 {
            for i in 0..n {
                let mut pull = Vec2::ZERO;
                for &j in &neighbors[i] {
                    let ju = j as usize;
                    let delta = self.x[i] - self.x[ju];
                    let r = delta.length();
                    if r <= 1e-9 {
                        continue;
                    }
                    // `mass * kernel(r) / rest_density` is dimensionless (kernel is normalised to
                    // 1/area, so `mass * kernel` is a density and dividing by `rest_density`
                    // cancels units) -- roughly 1 for a particle deep in bulk fluid, tapering to 0
                    // near a free surface, which is what makes `cohesion_accel` a genuine peak
                    // acceleration scale rather than needing its own separate normalisation.
                    let w = mass * cohesion_kernel(r, h) / rest_density;
                    pull -= w * (delta / r);
                }
                cohesion_velocity_delta[i] = cohesion_accel * pull * dt;
            }
        }

        // --- 3.6b Fluid-boundary adhesion (`slurry.wettability`). Each ball's circumference is
        // sampled into boundary particles ([`ball_boundary_particles`]); a fluid particle within
        // the adhesion kernel's support of one is pulled toward it. `h_adhesion` caps the kernel's
        // support radius at `2 * balls.radius` regardless of the fluid's own `h`: at this
        // project's coupling resolution `h` can exceed the ball's own radius several times over,
        // and using it directly (as the shell mechanism this replaces did, before its own
        // separate `ADHESION_RANGE_BALL_RADII` cap) let a distant fluid particle read as a
        // floating clump rather than a clinging film (the mid-air ball+slurry clumps an earlier
        // review reported) -- see docs/PHYSICS.md ss6.1a's history.
        //
        // **Bounded per ball, not per contributing particle.** Mirrors step 6.5's viscous-drag
        // closure exactly, for the identical reason (an amplification-by-neighbour-count bug
        // caught by `coupling::tests::ball_dropped_into_a_pool_decelerates_more_than_a_dry_drop`,
        // fixed there first): a taper-weighted (here, kernel-weighted) mean pull direction and
        // coverage fraction give *one* target closing velocity per ball, mass-ratio-relaxed
        // toward the contributing fluid mass (`attach`, vanishing as that mass -> 0 so a single
        // nearby SPH particle cannot move an immovable ball), then distributed back across the
        // contributing particles proportional to their own weight so Newton's third law holds
        // exactly regardless of how many particles are actually in range. **Balls are processed
        // sequentially (Gauss-Seidel)**, same reason as step 6.5: a fluid particle can be near
        // more than one ball's boundary at once at this project's coupling resolution.
        let mut adhesion_velocity_delta = vec![Vec2::ZERO; n];
        let wettability = slurry.wettability.clamp(0.0, 1.0);
        if wettability > 0.0 && !balls.is_empty() && balls.radius > 0.0 {
            let dx = 0.5 * h; // seed_lattice defines h = 2*dx; see ball_boundary_particles's doc
            let h_adhesion = h.min(2.0 * balls.radius);
            let boundary_pos = ball_boundary_particles(balls, dx);
            let count_per_ball = boundary_points_per_ball(balls.radius, dx) as usize;
            if !boundary_pos.is_empty() && h_adhesion > 0.0 {
                let adhesion_accel = wettability * ADHESION_ACCEL_FACTOR * GRAVITY.abs();
                for b in 0..balls.len() {
                    if !balls.x[b].is_finite() {
                        continue; // see the identical guard rationale in steps 3.5/6.5/6.6
                    }
                    let start = b * count_per_ball;
                    let end = (start + count_per_ball).min(boundary_pos.len());
                    // Merge every boundary point's contribution to the same fluid particle before
                    // it enters `items`, so a particle near more than one of this ball's own
                    // boundary points (routine: they are ~dx apart, well inside `h_adhesion` of
                    // each other at this project's typical resolutions) contributes once with its
                    // combined weight/direction, not once per boundary point.
                    let mut items: Vec<(usize, Vec2, f32)> = Vec::new();
                    for &bp in &boundary_pos[start..end] {
                        let mut nearby: Vec<u32> = Vec::new();
                        grid.for_each_near(bp, |j| nearby.push(j));
                        for j in nearby {
                            let ju = j as usize;
                            let delta = bp - self.x[ju]; // points from the fluid toward the ball
                            let r = delta.length();
                            if !r.is_finite() || r <= 0.5 * h_adhesion || r >= h_adhesion {
                                continue;
                            }
                            let w = mass * adhesion_kernel(r, h_adhesion) / rest_density;
                            if w <= 0.0 {
                                continue;
                            }
                            let n_hat = delta / r;
                            match items.iter_mut().find(|(idx, _, _)| *idx == ju) {
                                Some(entry) => {
                                    entry.1 += w * n_hat;
                                    entry.2 += w;
                                }
                                None => items.push((ju, w * n_hat, w)),
                            }
                        }
                    }
                    if items.is_empty() {
                        continue;
                    }
                    let mut weight_sum = 0.0f32;
                    // Weighted mean direction the *fluid* gets pulled in (toward the ball's
                    // nearby boundary points -- each `dir` is already `w * n_hat`,
                    // `n_hat = (boundary - fluid)/r`).
                    let mut fluid_pull_dir_sum = Vec2::ZERO;
                    for &(_, dir, w) in &items {
                        fluid_pull_dir_sum += dir;
                        weight_sum += w;
                    }
                    if weight_sum <= 1e-9 {
                        continue;
                    }
                    let fluid_pull_len = fluid_pull_dir_sum.length();
                    if fluid_pull_len <= 1e-9 {
                        // Coverage is (near-)symmetric around the ball (e.g. fully submerged) --
                        // no net direction to pull, same as a fully-entrained ball feeling no net
                        // drag (step 6.5).
                        continue;
                    }
                    // The ball's own reaction points the opposite way -- toward its fluid, not
                    // toward its own surface.
                    let ball_pull_dir = -fluid_pull_dir_sum / fluid_pull_len;
                    let coverage = weight_sum.min(1.0);
                    let m_contrib = mass * weight_sum;
                    let attach = m_contrib / (balls.mass + m_contrib);
                    let dv_b = ball_pull_dir * (adhesion_accel * coverage * dt * attach);
                    coupling.impulses[b] += balls.mass * dv_b;
                    // Approximate contact point (a single mean lever arm across every
                    // contributing boundary point, not one lever per pair) -- the same level of
                    // simplification the drag closure (step 6.5) and this mechanism's predecessor
                    // already used for their own aggregate reactions.
                    let lever = ball_pull_dir * balls.radius;
                    coupling.angular_impulses[b] += balls.mass * cross2(lever, dv_b);

                    for (ju, _dir, w) in items {
                        let frac = w / weight_sum;
                        // `-dv_b` points toward the ball (opposite `ball_pull_dir`), i.e. exactly
                        // `fluid_pull_dir_sum`'s direction -- consistent with each pair's own
                        // `n_hat`.
                        let dv_i = -(balls.mass * frac / mass) * dv_b;
                        adhesion_velocity_delta[ju] += dv_i;
                        coupling.fluid_momentum_change += mass * dv_i;
                    }
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
        // Remove the part of step 3.5's overlap push that was purely geometric (see that step's
        // doc comment) so it doesn't silently show up as fluid velocity via the reconstruction
        // above -- velocity is linear in position, so this cleanly isolates just that push's
        // excess regardless of what else (the density constraint, the wall projection) also moved
        // this particle this sub-step.
        for (v, &excess) in self.v.iter_mut().zip(&push_velocity_excess) {
            *v -= excess;
        }
        // Add step 3.6's cohesion/adhesion kicks back in now that reconstruction is done (see
        // that step's deferral comments above).
        for (v, &delta) in self.v.iter_mut().zip(&cohesion_velocity_delta) {
            *v += delta;
        }
        for (v, &delta) in self.v.iter_mut().zip(&adhesion_velocity_delta) {
            *v += delta;
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

        // --- 6.5 Ball viscous drag (docs/PLAN.md ss3.4 step 2, docs/PHYSICS.md ss6.2): each ball
        // relaxes toward its locally-entrained fluid mass via a Stokes-regime relaxation time
        // `tau = rho_ball * r^2 / (4 * mu)`, `beta = ball_no_slip * (1 - exp(-dt / tau))`, same
        // closure as before. The exchange is a **centre-of-mass relaxation between the ball and
        // its entrained fluid mass**, not a per-particle blend: an earlier version applied `beta`
        // independently to *every* nearby fluid particle and let the ball absorb the sum of all
        // those reactions, which amplifies the ball's actual response by the entrained/ball mass
        // ratio (`~9x` at this project's defaults) -- harmless at low viscosity (the amplified
        // response was still small) but an explicit over-relaxation once `beta` itself approaches
        // 1 (mu >= ~10 Pa*s), which the stability clamp (step 8) then masked as "clamping is
        // normal" rather than surfacing as the real bug it was. Weighting every nearby fluid
        // particle `i` by a smooth taper `phi_i = poly6(|r_i|^2, h_c) / poly6(0, h_c)`
        // (`r_i = x_i - x_b`, `h_c = balls.radius + h`, going to zero at the cutoff instead of a
        // top hat) and mass `w_i = particle_mass * phi_i`:
        //
        //   m_ent  = sum_i w_i                              entrained fluid mass
        //   v_bar  = (sum_i w_i * v_i) / m_ent               entrained-mass velocity
        //   a_lin  = beta * m_ent / (balls.mass + m_ent)     <= 1 for every mu, dt (the fix)
        //   dv_b   = a_lin * (v_bar - v_b)
        //
        // and the mirror-image moment-weighted closure for the ball's spin (`i_ent`, `omega_bar`,
        // `a_rot`, `dw_b`, using `balls.inertia` in place of `balls.mass`). `dv_b`/`dw_b` are then
        // the ball's impulse (`impulses[b] += balls.mass * dv_b`, no `/dt`, same convention as
        // before); the fluid's reaction is distributed back across the same weighted neighbours
        // so linear and angular momentum are conserved *exactly* (not just bounded), including
        // the cross term from distributing an off-centre reaction (`m_b * cross2(r_bar, dv_b)`,
        // see docs/PHYSICS.md ss6.2 for the full derivation of why this term is required for exact
        // conservation once the linear reaction is spread non-uniformly over an off-centre `r_bar`
        // -- see [`CouplingImpulses`]'s doc comment for why impulses, not forces, throughout).
        //
        // **Balls are processed sequentially (Gauss-Seidel), not Jacobi.** `self.v` is updated
        // immediately after each ball's contribution, not accumulated separately and applied once
        // at the end: at this project's coupling resolution, `h_c` easily spans several
        // neighbouring balls' worth of packed charge, so a fluid particle commonly sits within
        // more than one ball's entrainment radius at once. Computing every ball's `v_bar` against
        // the *same* stale fluid state and only summing the results onto the fluid afterwards
        // (Jacobi-style) lets several balls' otherwise-individually-bounded (`a_lin <= 1`)
        // corrections stack on the same fluid particle, well beyond what any single ball's own
        // bound was meant to allow -- reproducibly driving the fluid to its speed clamp (step
        // 7.5) even in a perfectly still drum with no rotation to drive it. Sequential
        // application means each ball after the first already sees the fluid's up-to-date state,
        // so its own relaxation is computed (and bounded) against reality rather than against a
        // stale snapshot every other overlapping ball is *also* independently correcting.
        let mu = slurry.viscosity_pa_s.max(0.0);
        if !balls.is_empty() && balls.radius > 0.0 {
            let beta_b = slurry.ball_no_slip.clamp(0.0, 1.0);
            let rho_ball = balls.mass / (std::f32::consts::PI * balls.radius * balls.radius);
            let beta = if mu > 1e-6 && rho_ball > 0.0 {
                let tau = (rho_ball * balls.radius * balls.radius / (4.0 * mu)).max(1e-9);
                beta_b * (1.0 - (-dt / tau).exp())
            } else {
                0.0
            };
            if beta > 0.0 {
                let h_c = balls.radius + h;
                let w_hc0 = poly6(0.0, h_c).max(1e-12);
                // Coupling-only grid, sized to this step's actual query radius (`h_c`), queried
                // ball -> fluid (the reverse direction from the fluid -> ball `ball_grid` above,
                // which is sized for step 3.5's much smaller `contact_radius`). Built once from
                // the pre-loop fluid positions; positions don't change in this step (only
                // velocities do), so the neighbour sets themselves stay valid throughout even
                // though `self.v` is updated in place as the loop progresses.
                let coupling_fluid_grid = UniformGrid::build(&self.x, h_c.max(1e-6));
                let i_b = balls.inertia;
                for b in 0..balls.len() {
                    if !balls.x[b].is_finite() {
                        // See the identical rationale in step 6.6 below: a non-finite ball
                        // position must never manufacture a plausible-looking reaction.
                        continue;
                    }
                    let mut nearby: Vec<u32> = Vec::new();
                    coupling_fluid_grid.for_each_near(balls.x[b], |j| nearby.push(j));
                    if nearby.is_empty() {
                        continue;
                    }
                    let mut items: Vec<(usize, Vec2, f32)> = Vec::with_capacity(nearby.len());
                    let mut m_ent = 0.0f32;
                    let mut r_sum = Vec2::ZERO;
                    let mut v_sum = Vec2::ZERO;
                    let mut i_ent = 0.0f32;
                    let mut l_sum = 0.0f32;
                    for j in nearby {
                        let ju = j as usize;
                        let r_i = self.x[ju] - balls.x[b];
                        let r2 = r_i.length_squared();
                        if !r2.is_finite() || r2 >= h_c * h_c {
                            continue;
                        }
                        let phi = poly6(r2, h_c) / w_hc0;
                        if phi <= 0.0 {
                            continue;
                        }
                        let w = mass * phi;
                        m_ent += w;
                        r_sum += w * r_i;
                        v_sum += w * self.v[ju];
                        i_ent += w * r2;
                        l_sum += w * cross2(r_i, self.v[ju]);
                        items.push((ju, r_i, phi));
                    }
                    if m_ent <= 1e-12 || items.is_empty() {
                        continue;
                    }

                    let r_bar = r_sum / m_ent;
                    let v_bar = v_sum / m_ent;
                    let a_lin = beta * m_ent / (balls.mass + m_ent);
                    let dv_b = a_lin * (v_bar - balls.v[b]);

                    // Angular relaxation, and the tangential field that distributes its reaction
                    // back to the fluid with zero net linear momentum (mean-subtracted around
                    // `p_bar = perp(r_bar)`, so `sum_i w_i * (perp(r_i) - p_bar) == 0` exactly).
                    // Only actually applied when the distribution isn't degenerate (`s` nonzero)
                    // -- otherwise skipped for *both* sides this sub-step, since giving the ball
                    // an angular impulse without a matching fluid reaction would break
                    // conservation.
                    let mut dw_b = 0.0f32;
                    let mut c_rot = 0.0f32;
                    let mut p_bar = Vec2::ZERO;
                    if i_ent > 1e-12 {
                        let omega_bar = l_sum / i_ent;
                        let a_rot = beta * i_ent / (i_b + i_ent);
                        let candidate_dw = a_rot * (omega_bar - balls.omega[b]);
                        p_bar = Vec2::new(-r_bar.y, r_bar.x);
                        // `s` (the exact-conservation normalizer, see below) is
                        // `sum_i mass*phi_i^2*(|r_i|^2 - cross2(r_i, p_bar))`: a positive-definite
                        // term (`|r_i|^2`) minus a cross term that can be comparable in size to it
                        // for perfectly ordinary neighbour geometries -- not just contrived
                        // degenerate ones -- whenever the entrained neighbour count is small
                        // (routine at this project's fine coupling resolution relative to a
                        // ball's size, docs/PHYSICS.md ss6.x). `s` can therefore land anywhere
                        // from comparable to that positive scale down to many orders of magnitude
                        // smaller, essentially at random depending on exact particle placement.
                        // Comparing `|s|` against a *relative* fraction of its own
                        // positive-definite part (rather than a fixed absolute epsilon, which
                        // would only catch an exact-zero-like degeneracy) is what actually
                        // detects "this division is about to amplify without bound" -- a fixed
                        // `1e-9` absolute threshold let this fire routinely and was the source of
                        // a genuine coupling instability (a ball+slurry charge failing to settle
                        // even in a stationary drum, `c_rot` observed in the hundreds to tens of
                        // thousands here).
                        let mut s = 0.0f32;
                        let mut s_scale = 0.0f32;
                        for &(_, r_i, phi) in &items {
                            let perp_i = Vec2::new(-r_i.y, r_i.x);
                            let phi2 = phi * phi;
                            s += mass * phi2 * cross2(r_i, perp_i - p_bar);
                            s_scale += mass * phi2 * r_i.length_squared();
                        }
                        if s.abs() > 0.1 * s_scale {
                            dw_b = candidate_dw;
                            c_rot = -i_b * dw_b / s;
                        }
                    }

                    for &(ju, r_i, phi) in &items {
                        let mut dv_i = -(balls.mass / m_ent) * dv_b * phi;
                        if c_rot != 0.0 {
                            let perp_i = Vec2::new(-r_i.y, r_i.x);
                            dv_i += c_rot * phi * (perp_i - p_bar);
                        }
                        // Applied immediately (Gauss-Seidel), not accumulated for later -- see
                        // this step's doc comment for why.
                        self.v[ju] += dv_i;
                        coupling.fluid_momentum_change += mass * dv_i;
                    }

                    // `dv_b`/`dw_b` are velocity changes, so the ball's impulse is `mass * dv_b`
                    // directly (no `/ dt`, unlike step 3.5's position-derived impulse). The
                    // `balls.mass * cross2(r_bar, dv_b)` term is the reaction torque from
                    // distributing that linear reaction over an off-centre `r_bar` rather than
                    // exactly at the ball's centre -- see this step's doc comment above.
                    coupling.impulses[b] += balls.mass * dv_b;
                    coupling.angular_impulses[b] += i_b * dw_b + balls.mass * cross2(r_bar, dv_b);
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

        // --- 7. Implicit Newtonian viscosity (see module doc comment) -----------------------
        let mu = slurry.viscosity_pa_s.max(0.0);
        let mut viscosity_iterations = 0u32;
        // Recompute density once more at the final (post-boundary-projection) positions, so both
        // the viscosity solve and the shear-rate readout below use up-to-date neighbour densities.
        for i in 0..n {
            let mut rho = mass * poly6(0.0, h);
            for &j in &neighbors[i] {
                let r2 = (self.x[i] - self.x[j as usize]).length_squared();
                rho += mass * poly6(r2, h);
            }
            density[i] = rho;
        }
        if mu > 0.0 {
            let weights = morris_weights(&self.x, &neighbors, &density, mass, mu, h);
            viscosity_iterations = solve_implicit_viscosity(&neighbors, &weights, &mut self.v, dt);
        }
        let mean_shear_rate_per_s =
            mean_shear_rate(&self.x, &self.v, &neighbors, &density, mass, h);

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
        let coupling_clamp_hits = coupling.clamp_hits;
        (
            coupling,
            FluidStepStats {
                viscosity_iterations,
                mean_shear_rate_per_s,
                coupling_clamp_hits,
            },
        )
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
        let row_h = dx * (3.0f32.sqrt() * 0.5);
        let target_area = slurry.fill_fraction * std::f32::consts::PI * radius_m * radius_m;
        let expected_count = (target_area / (dx * row_h)).round() as usize;
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

    /// Builds a still-drum scenario with one ball at rest at the origin, surrounded by a dense
    /// uniform patch of fluid particles all moving at `fluid_velocity`, filling a large enough
    /// area that the entrained mass `m_ent` (within `balls.radius + h` of the ball) is much
    /// greater than the ball's own mass -- isolating step 6.5's ball<->fluid drag exchange in the
    /// regime where its centre-of-mass relaxation (docs/PHYSICS.md ss6.2) should reduce to the
    /// ball simply relaxing toward `fluid_velocity` at rate `beta`, independent of the fluid's own
    /// mass (no wall, no gravity yet applied, no buoyancy since the ball isn't tracked by
    /// `step_coupled`'s buoyancy step until `dt` has actually elapsed once).
    fn ball_in_uniform_fluid_patch(
        ball_radius_m: f32,
        ball_density_kg_m3: f32,
        slurry: &SlurryParams,
        fluid_velocity: Vec2,
    ) -> (Balls, FluidParticles) {
        let mass = crate::dem::ball_mass(2.0 * ball_radius_m, ball_density_kg_m3);
        let inertia = crate::dem::ball_inertia(mass, ball_radius_m);
        let balls = Balls {
            x: vec![Vec2::ZERO],
            v: vec![Vec2::ZERO],
            theta: vec![0.0],
            omega: vec![0.0],
            radius: ball_radius_m,
            mass,
            inertia,
        };

        // Fluid lattice spacing large enough (relative to the ball) that `h_c = ball_radius + h`
        // spans an area whose total fluid mass `m_ent` dominates the ball's own mass -- `m_ent`
        // scales with the *area* under the smooth taper (`~ rest_density * h_c^2`), not with
        // particle count, so this requires a large `h` relative to `ball_radius`, not a fine
        // lattice (`h = 2*dx`, so `dx` itself must be large relative to `ball_radius`).
        let dx = ball_radius_m * 4.0;
        let h = 2.0 * dx;
        let rest_density = slurry.density_kg_m3;
        let particle_mass = rest_density * dx * dx;
        // covers h_c = radius + h with margin
        let half_span = ball_radius_m + h;
        // Exclude the ball's own overlap zone (step 3.5's `contact_radius`) -- a fluid particle
        // placed on top of the ball would trigger a large overlap-push impulse that has nothing
        // to do with the drag mechanism (step 6.5) this test isolates.
        let exclude_radius = 1.01 * (ball_radius_m + 0.25 * h);
        let n = (2.0 * half_span / dx).ceil() as i32;
        let mut x = Vec::new();
        let mut v = Vec::new();
        let mut dye = Vec::new();
        for iy in -n..=n {
            for ix in -n..=n {
                let p = Vec2::new(ix as f32 * dx, iy as f32 * dx);
                let dist = p.length();
                if dist <= half_span && dist >= exclude_radius {
                    x.push(p);
                    v.push(fluid_velocity);
                    dye.push(0.0);
                }
            }
        }
        let fluid = FluidParticles {
            x,
            v,
            dye,
            particle_mass,
            h,
            rest_density,
        };
        (balls, fluid)
    }

    #[test]
    fn ball_drag_matches_two_dimensional_stokes_scaling() {
        // docs/PHYSICS.md ss6.2: in the entrained-mass-dominated, small-`beta` regime (`m_ent >>
        // balls.mass`, `dt << tau`), the centre-of-mass relaxation's ball impulse reduces
        // analytically to `impulse ~= 4*pi*mu*(v_fluid - v_ball)*dt` -- ordinary 2D Stokes drag on
        // a disc, with no dependence on the ball/fluid mass ratio (unlike the removed
        // per-particle blend, which was off by the entrained/ball mass ratio, ~9x at this
        // project's own defaults). This directly checks that calibration, and that increasing
        // viscosity roughly linearly increases the drag response.
        let drum = still_drum(2.0); // large, unused directly (no wall contact at the origin)
        let ball_radius_m = 0.005;
        let ball_density = 6000.0;
        let fluid_velocity = Vec2::new(1.0, 0.0);
        let dt = 1.0 / 240.0;

        let stokes_drag_impulse = |mu: f32| -> f32 {
            let slurry = SlurryParams {
                viscosity_pa_s: mu,
                fill_fraction: 0.0, // unused here (fluid is placed manually)
                wettability: 0.0,   // isolates viscous drag (step 6.5) from adhesion (step 3.6);
                // this test's own `exclude_radius` only clears step 3.5's small overlap zone, not
                // step 3.6's much wider adhesion shell
                ..SlurryParams::default()
            };
            let (mut balls, mut fluid) =
                ball_in_uniform_fluid_patch(ball_radius_m, ball_density, &slurry, fluid_velocity);
            let (impulses, _stats) = fluid.step_coupled(&drum, 0.0, &slurry, 1, dt, &balls);
            balls.v[0] += impulses.impulses[0] / balls.mass;
            impulses.impulses[0].length()
        };

        let mu_low = 0.5f32;
        let mu_high = 5.0f32;
        let impulse_low = stokes_drag_impulse(mu_low);
        let impulse_high = stokes_drag_impulse(mu_high);

        let predicted = |mu: f32| 4.0 * std::f32::consts::PI * mu * fluid_velocity.x * dt;
        for (mu, impulse) in [(mu_low, impulse_low), (mu_high, impulse_high)] {
            let expected = predicted(mu);
            let rel_err = (impulse - expected).abs() / expected;
            assert!(
                rel_err < 0.5,
                "mu={mu}: impulse={impulse} far from 2D-Stokes prediction {expected} \
                 (rel_err={rel_err})"
            );
        }
        assert!(
            impulse_high > impulse_low * 3.0,
            "drag should scale up substantially with viscosity: \
             impulse(mu={mu_low})={impulse_low}, impulse(mu={mu_high})={impulse_high}"
        );
    }

    #[test]
    fn ball_drag_never_overshoots_the_local_fluid_velocity() {
        // Even at a viscosity high enough to saturate `beta` toward 1 within a single sub-step
        // (docs/PHYSICS.md ss6.2), the ball's resulting velocity must never overshoot past the
        // entrained fluid's velocity, and the stability clamp (step 8) should not need to fire --
        // the centre-of-mass relaxation's `a <= 1` bound (the actual A1 fix) should make it
        // unnecessary here, unlike the removed per-particle blend.
        let drum = still_drum(2.0);
        let ball_radius_m = 0.005;
        let ball_density = 6000.0;
        let fluid_velocity = Vec2::new(1.0, 0.0);
        let dt = 1.0 / 240.0;
        let slurry = SlurryParams {
            viscosity_pa_s: 200.0,
            fill_fraction: 0.0,
            wettability: 0.0, // isolates drag (step 6.5) from adhesion (step 3.6)
            ..SlurryParams::default()
        };
        let (balls, mut fluid) =
            ball_in_uniform_fluid_patch(ball_radius_m, ball_density, &slurry, fluid_velocity);

        let (impulses, _stats) = fluid.step_coupled(&drum, 0.0, &slurry, 1, dt, &balls);
        assert_eq!(
            impulses.clamp_hits, 0,
            "a bounded (a<=1) drag relaxation should not need the stability clamp"
        );
        let new_v = balls.v[0] + impulses.impulses[0] / balls.mass;
        assert!(
            new_v.x >= 0.0 && new_v.x <= fluid_velocity.x * 1.001,
            "ball overshot the fluid velocity it was relaxing toward: new_v={new_v:?}, \
             fluid_velocity={fluid_velocity:?}"
        );
    }

    #[test]
    fn coupling_finds_ball_neighbours_beyond_two_radii() {
        // Regression for the `ball_grid`/coupling-fluid-grid cell-size fix (docs/PHYSICS.md
        // ss6.1/6.4): a fluid particle at `r + 0.5*h` from a ball -- farther than `2*radius` at
        // this test's proportions -- must still be found and produce a non-zero drag reaction.
        let drum = still_drum(2.0);
        let ball_radius_m = 0.01;
        let slurry = SlurryParams {
            viscosity_pa_s: 50.0,
            fill_fraction: 0.0,
            ..SlurryParams::default()
        };
        let dx = ball_radius_m * 1.5;
        let h = 2.0 * dx;
        assert!(
            ball_radius_m + 0.5 * h > 2.0 * ball_radius_m,
            "test setup must place the fluid particle beyond 2*radius"
        );
        let mass = crate::dem::ball_mass(2.0 * ball_radius_m, 6000.0);
        let inertia = crate::dem::ball_inertia(mass, ball_radius_m);
        let balls = Balls {
            x: vec![Vec2::ZERO],
            v: vec![Vec2::ZERO],
            theta: vec![0.0],
            omega: vec![0.0],
            radius: ball_radius_m,
            mass,
            inertia,
        };
        let mut fluid = FluidParticles {
            x: vec![Vec2::new(ball_radius_m + 0.5 * h, 0.0)],
            v: vec![Vec2::new(1.0, 0.0)],
            dye: vec![0.0],
            particle_mass: slurry.density_kg_m3 * dx * dx,
            h,
            rest_density: slurry.density_kg_m3,
        };

        let (impulses, _stats) = fluid.step_coupled(&drum, 0.0, &slurry, 1, 1.0 / 240.0, &balls);
        assert!(
            impulses.impulses[0].length() > 1e-9,
            "expected a non-zero drag reaction from a fluid particle beyond 2*radius"
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
    fn viscosity_solve_is_a_no_op_at_zero_viscosity() {
        // mu = 0 must skip the CG solve entirely (no diffusion, no iteration cost) rather than
        // running an expensive no-op solve.
        let slurry = SlurryParams {
            fill_fraction: 0.2,
            viscosity_pa_s: 0.0,
            ..SlurryParams::default()
        };
        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);
        let stats = fluid.step(&drum, 0.0, &slurry, 3, 1.0 / 240.0);
        assert_eq!(stats.viscosity_iterations, 0);
    }

    #[test]
    fn viscosity_keeps_mattering_between_100_and_200_pa_s() {
        // Regression: the old qualitative XSPH coefficient this solver replaced saturated by
        // ~100 Pa*s, so 100 vs 200 Pa*s produced bit-identical decay -- the UI's viscosity
        // slider did nothing in the top half of its advertised range. The physically implicit
        // solver used now should still show a measurable difference this far up the range, even
        // though the unconditionally-stable implicit scheme itself saturates *somewhat* at very
        // large `dt * (viscosity term)` (that is expected and fine -- it's what keeps the solve
        // stable without a vanishingly small sub-step -- the point here is only that it must not
        // be bit-identical the way the old ordinal coefficient was).
        let (rate_100, rate_200) = (
            checkerboard_shear_decay_rate(100.0),
            checkerboard_shear_decay_rate(200.0),
        );
        // Margin lowered from 1.1x to 1.05x when the hexagonal-lattice seeding mass fix (see
        // `FluidParticles::seed_lattice`'s doc comment) shifted the settled lattice's local
        // densities slightly, which the implicit solve's `nu = mu/rho_i` is sensitive to; measured
        // ~1.095x after that fix (was ~1.13x before), still comfortably distinguishing this from
        // the pre-implicit-solver saturation this test guards against (which was bit-identical,
        // i.e. a ratio of exactly 1.0).
        assert!(
            rate_200 > rate_100 * 1.05,
            "100 vs 200 Pa*s should still measurably differ: rate_100={rate_100} rate_200={rate_200}"
        );
    }

    #[test]
    fn viscosity_solve_stays_finite_and_converges_at_high_viscosity_and_rotation() {
        // Far beyond anything the UI permits (omega = 6 rad/s at this drum radius), with the
        // top of the advertised viscosity range, to exercise the CG solve under a genuinely
        // stiff, fast-moving field. Checked on *every* sub-step, not just at the end -- a
        // transient excursion (non-finite velocity, or the solve failing to converge, i.e.
        // hitting VISCOSITY_CG_MAX_ITERS) that later recovers is exactly the kind of event a
        // final-state-only check would miss.
        let slurry = SlurryParams {
            fill_fraction: 0.2,
            viscosity_pa_s: 200.0,
            ..SlurryParams::default()
        };
        let radius_m = 0.5;
        let omega = 6.0;
        let drum = Drum::new(
            radius_m,
            omega,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        );
        let dt = 1.0 / 240.0;
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);
        let mut drum_angle = 0.0f32;
        for s in 0..(2 * 240) {
            let stats = fluid.step(&drum, drum_angle, &slurry, 3, dt);
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
            assert!(
                stats.viscosity_iterations < VISCOSITY_CG_MAX_ITERS,
                "sub-step {s}: viscosity CG hit the iteration ceiling without converging"
            );
            assert!(
                stats.mean_shear_rate_per_s.is_finite() && stats.mean_shear_rate_per_s >= 0.0,
                "sub-step {s}: non-finite or negative mean shear rate: {}",
                stats.mean_shear_rate_per_s
            );
            for (i, &v) in fluid.v.iter().enumerate() {
                assert!(
                    v.is_finite(),
                    "sub-step {s}, particle {i}: non-finite velocity {v:?}"
                );
            }
        }
    }

    /// Shared harness for the checkerboard-perturbation viscosity tests: settles a puddle once
    /// at a fixed baseline viscosity, imposes a spatially-checkerboarded (the highest frequency
    /// this particle spacing can resolve) x-velocity perturbation, then measures the perturbation
    /// energy's exponential decay rate after stepping for `n_steps` at `mu`. Settling once and
    /// re-running the decay phase from an identical clone of that state at each candidate `mu`
    /// (rather than settling independently per `mu`) isolates `mu`'s effect from per-mu settling
    /// divergence -- the same rationale
    /// `higher_viscosity_damps_a_high_frequency_perturbation_more_across_the_whole_range` uses.
    fn checkerboard_shear_decay_rate(mu: f32) -> f32 {
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
        const PERTURBATION_MPS: f32 = 2.0;
        for i in 0..baseline.len() {
            let col = (baseline.x[i].x / dx).round() as i64;
            baseline.v[i].x += if col % 2 == 0 {
                PERTURBATION_MPS
            } else {
                -PERTURBATION_MPS
            };
        }
        let initial_energy: f32 =
            baseline.v.iter().map(|v| v.x * v.x).sum::<f32>() / baseline.len() as f32;

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
        let n_steps = 20;
        let dt = 1.0 / 240.0;
        for _ in 0..n_steps {
            fluid.step(&drum, 0.0, &slurry, 3, dt);
        }
        let final_energy: f32 = fluid.v.iter().map(|v| v.x * v.x).sum::<f32>() / fluid.len() as f32;
        // Perturbation energy decays as exp(-2 * decay_rate * t) (energy ~ amplitude^2); solve
        // for decay_rate from the before/after ratio, clamped so a fully-decayed or noisy
        // final_energy can't produce a nonsensical negative/infinite rate.
        let ratio = (final_energy / initial_energy.max(1e-12)).clamp(1e-6, 1.0);
        -0.5 * ratio.ln() / (n_steps as f32 * dt)
    }

    #[test]
    fn higher_viscosity_damps_a_high_frequency_perturbation_more_across_the_whole_range() {
        // The implicit solver diffuses velocity toward the local neighbourhood average, so a
        // spatially checkerboarded (the highest frequency this particle spacing can resolve)
        // velocity perturbation should be damped more by higher viscosity -- including from 100
        // to 200 Pa*s, which the old hard-clamped XSPH mapping this solver replaced made
        // bit-identical (see `viscosity_keeps_mattering_between_100_and_200_pa_s` for that
        // specific regression). A *linear* shear profile would not distinguish any viscosity
        // here: a symmetric averaging kernel reproduces a linear field exactly, so only a
        // genuinely high-frequency pattern exercises the solver's viscosity-dependence.
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

        let (impulses, _stats) = fluid.step_coupled(&drum, 0.0, &slurry, 3, dt, &dem.balls);

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

    #[test]
    fn cohesion_pulls_two_isolated_fluid_particles_together_only_when_surface_tension_is_positive()
    {
        // Regression for fluid-fluid cohesion (step 3.6a, `slurry.surface_tension_n_m`),
        // replacing the disabled `S_CORR_K` artificial-pressure term as this solver's only
        // fluid-fluid attraction (docs/PHYSICS.md ss5.2's former "no fluid-fluid cohesion/
        // surface-tension term" limitation). Two particles placed near the cohesion kernel's
        // peak (`r ~= 0.6h`) should stay put at `surface_tension_n_m = 0` and be pulled
        // together, symmetrically, when positive.
        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let dt = 1.0 / 240.0;

        let run = |surface_tension_n_m: f32| -> (Vec2, Vec2, f32) {
            let slurry = SlurryParams {
                fill_fraction: 0.0,
                wall_no_slip: 0.0,
                viscosity_pa_s: 0.0,
                surface_tension_n_m,
                ..SlurryParams::default()
            };
            let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 20, &[], 0.0);
            assert!(fluid.is_empty());
            let sep = 0.6 * fluid.h; // near the cohesion kernel's peak
            fluid.x.push(Vec2::new(-0.5 * sep, 0.0));
            fluid.x.push(Vec2::new(0.5 * sep, 0.0));
            fluid.v.push(Vec2::ZERO);
            fluid.v.push(Vec2::ZERO);
            fluid.dye.push(0.0);
            fluid.dye.push(0.0);

            fluid.step(&drum, 0.0, &slurry, 1, dt);
            (fluid.v[0], fluid.v[1], fluid.particle_mass)
        };

        let (v0_off, v1_off, _mass) = run(0.0);
        assert_eq!(
            v0_off.x, 0.0,
            "surface_tension_n_m = 0 must not move particles horizontally"
        );
        assert_eq!(
            v1_off.x, 0.0,
            "surface_tension_n_m = 0 must not move particles horizontally"
        );

        let (v0_on, v1_on, mass) = run(0.072);
        assert!(
            v0_on.x > 1e-6,
            "the left particle should be pulled toward the right one (positive x): v0={v0_on:?}"
        );
        assert!(
            v1_on.x < -1e-6,
            "the right particle should be pulled toward the left one (negative x): v1={v1_on:?}"
        );
        // Newton's third law: cohesion is a symmetric pairwise force between two otherwise-
        // identical particles, so the total x-momentum change must cancel exactly.
        let sum_x = mass * (v0_on.x + v1_on.x);
        assert!(
            sum_x.abs() < 1e-6,
            "momentum not conserved in isolated cohesion contact: v0={v0_on:?} v1={v1_on:?}"
        );
    }
}
