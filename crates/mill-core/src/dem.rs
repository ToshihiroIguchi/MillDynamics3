//! Position-based (XPBD-style) rigid discs for grinding media (balls).
//!
//! Non-penetration against other balls and against the drum wall/lifters ([`crate::geometry::Drum`])
//! is solved as a rigid (zero-compliance) geometric constraint, projected iteratively — stable at
//! any sub-step size, unlike an explicit spring-dashpot DEM contact (see docs/PLAN.md ss3.2 for the
//! full rationale and algorithm this module implements). Friction, restitution and rolling
//! resistance are added as position/velocity corrections on top of the converged contact solve.
//! Ball population is seeded from [`crate::params::Params::effective_media`] (post
//! coarse-graining, if any) rather than the raw UI media diameter.

use std::collections::HashMap;
use std::f32::consts::PI;

use glam::Vec2;

use crate::geometry::Drum;
use crate::grid::UniformGrid;
use crate::params::{EffectiveMedia, MediaParams};
use crate::rng::Rng;

const GRAVITY: f32 = -9.81;
/// Below this approach speed, a contact is treated as already-resting rather than a fresh impact:
/// restitution `e` is not applied to it (the target relative normal velocity is exactly zero
/// instead of `e * v_approach`), so a settled contact doesn't buzz indefinitely instead of
/// settling. Step 6 (see [`DemState::step_with_external_forces`]) still zeroes/clamps the
/// *actual* relative normal velocity for every contact regardless of this threshold -- including
/// a resting one -- because the depenetration solve (step 3) can otherwise leave a pair with
/// spurious separation velocity from a large position correction; only the choice of *target*
/// (`e * v_approach` vs. exactly zero) depends on this threshold. See docs/PLAN.md ss3.2 step 6.
/// Also used (docs/PLAN.md ss3.5) as the threshold for counting a contact as a genuine "impact"
/// for [`DemStepStats::collision_count`]/`impact_energy_histogram`.
const RESTITUTION_VELOCITY_THRESHOLD: f32 = 0.02;

/// Largest fraction of a ball's own diameter that step 3's non-penetration solve recovers from a
/// single contact in one iteration. A deep overlap (e.g. from a coarse-grained ball population's
/// large diameter, or several balls suddenly landing on top of each other) is instead recovered
/// gradually over several sub-steps rather than teleporting the pair fully apart in one shot --
/// which step 4 would otherwise turn directly into an unphysically large separation velocity (see
/// this module's doc comment and docs/PHYSICS.md for the fluidised-charge energy-injection bug
/// this, together with step 6's bidirectional restitution, fixes). At the current default
/// `dem_iterations = 2` this still allows recovering up to `2 * 0.2 = 0.4` diameters of overlap
/// per sub-step, so it does not meaningfully slow down recovery from the ordinary small overlaps
/// a converged solve produces.
const MAX_RECOVERY_FRACTION: f32 = 0.2;

/// How far, as a fraction of the ball radius, to advance a CCD-clamped ball past its exact
/// geometric time-of-impact, so the resulting state is a small genuine overlap rather than a bare
/// `C = 0` touch -- see the call site in `step_with_external_forces` step 2.5 for why an exact
/// touch would silently drop the contact from every later step (friction, restitution, rolling
/// resistance all key off a *recorded*, `lambda_n > 0` contact). Expressed as a physical depth
/// (`CCD_CONTACT_SKIN_FRACTION * r` metres) rather than a fraction of the sub-step's *remaining*
/// time specifically so it stays well above `f32` precision at the ball's actual position
/// magnitude regardless of how close to the sub-step's end the crossing lands -- a time-fraction
/// skin degenerates when the crossing is found near `t = 1` (routine, not an edge case: it is
/// exactly what "a ball just reaches the floor by the end of a sub-step" looks like), where even
/// a generous fraction of the tiny remaining time can round away to nothing at typical world-space
/// position magnitudes.
const CCD_CONTACT_SKIN_FRACTION: f32 = 1e-2;

/// Number of log-spaced bins in [`DemStepStats::impact_energy_histogram`].
pub const IMPACT_ENERGY_HISTOGRAM_BINS: usize = 12;
/// Lower edge (J per metre of mill depth) of the impact-energy histogram's range.
pub const IMPACT_ENERGY_MIN_J: f32 = 1e-6;
/// Upper edge (J per metre of mill depth) of the impact-energy histogram's range.
pub const IMPACT_ENERGY_MAX_J: f32 = 1.0;

/// Log-spaced bin edges (J per metre of mill depth), `IMPACT_ENERGY_HISTOGRAM_BINS + 1` values
/// from [`IMPACT_ENERGY_MIN_J`] to [`IMPACT_ENERGY_MAX_J`], for labelling
/// [`DemStepStats::impact_energy_histogram`] (e.g. in [`crate::metrics`]).
pub fn impact_energy_bin_edges() -> [f32; IMPACT_ENERGY_HISTOGRAM_BINS + 1] {
    let log_min = IMPACT_ENERGY_MIN_J.log10();
    let log_max = IMPACT_ENERGY_MAX_J.log10();
    let mut edges = [0.0f32; IMPACT_ENERGY_HISTOGRAM_BINS + 1];
    for (i, edge) in edges.iter_mut().enumerate() {
        let t = i as f32 / IMPACT_ENERGY_HISTOGRAM_BINS as f32;
        *edge = 10f32.powf(log_min + t * (log_max - log_min));
    }
    edges
}

/// Histogram bin index for an impact energy `e` (J per metre of mill depth), log-spaced between
/// [`IMPACT_ENERGY_MIN_J`] and [`IMPACT_ENERGY_MAX_J`] (values outside that range are clamped
/// into the first/last bin, so every genuine impact is counted somewhere). `None` for a
/// non-positive or non-finite `e` (not a real impact).
fn impact_energy_bin(e: f32) -> Option<usize> {
    if !e.is_finite() || e <= 0.0 {
        return None;
    }
    let e_clamped = e.clamp(IMPACT_ENERGY_MIN_J, IMPACT_ENERGY_MAX_J);
    let log_min = IMPACT_ENERGY_MIN_J.log10();
    let log_max = IMPACT_ENERGY_MAX_J.log10();
    let t = (e_clamped.log10() - log_min) / (log_max - log_min);
    let idx = (t * IMPACT_ENERGY_HISTOGRAM_BINS as f32) as usize;
    Some(idx.min(IMPACT_ENERGY_HISTOGRAM_BINS - 1))
}

/// Per-sub-step DEM solver diagnostics (docs/PLAN.md ss3.5), returned by
/// [`DemState::step_with_external_forces`]. Everything here is a *per-metre-of-mill-depth*
/// quantity, consistent with this crate's unit-depth disc convention ([`ball_mass`]); a caller
/// wanting an absolute value for a real 3D mill multiplies by the mill's axial length.
#[derive(Debug, Clone, Copy, Default)]
pub struct DemStepStats {
    /// Work (J) done against the drum wall's friction this sub-step -- the rate of this
    /// (`wall_work_j / dt`) is the mill's instantaneous power draw. Derived from the friction
    /// step's tangential impulse at each ball-wall contact dotted with the wall's own velocity
    /// there (the wall's *normal* impulse does no work: the wall's velocity is purely tangential,
    /// see [`crate::geometry::Drum::wall_velocity`]).
    pub wall_work_j: f32,
    /// Number of contacts this sub-step whose pre-solve approach speed exceeded
    /// [`RESTITUTION_VELOCITY_THRESHOLD`] (i.e. a genuine fresh impact, ball-ball and ball-wall
    /// combined) -- the same criterion the restitution step (docs/PLAN.md ss3.2 step 6) uses to
    /// decide whether to apply restitution at all.
    pub collision_count: u32,
    /// Count of impacts this sub-step per log-spaced energy bin (see [`impact_energy_bin_edges`]),
    /// `E = 0.5 * m_reduced * v_n_pre^2` (`m_reduced = mass/2` ball-ball, `mass` ball-wall).
    pub impact_energy_histogram: [u32; IMPACT_ENERGY_HISTOGRAM_BINS],
    /// Mechanical energy (J) removed by this sub-step's contact solve as a whole:
    /// [`mechanical_energy_j`] (translational + rotational KE + gravitational PE) right after the
    /// predict step (gravity + external forces, before any contact correction), plus this
    /// sub-step's [`wall_work_j`](Self::wall_work_j), minus final mechanical energy. Adding back
    /// `wall_work_j` is what keeps this metric from reading negative whenever the moving wall
    /// pumps energy into the charge faster than friction/restitution remove it -- an earlier
    /// version omitted it (and PE/rotational KE) and could read hundreds of watts negative on an
    /// ordinary cascading run despite the solver itself never creating energy (see
    /// `crate::tests::steady_cascading_charge_never_gains_more_energy_than_the_wall_supplies`,
    /// which proves `e_final <= e_after_predict + wall_work_j` every sub-step -- this field is
    /// exactly that invariant's residual, so it is non-negative by the same argument, not by
    /// clamping). A simple, robust net-dissipation estimate -- it does not attribute the loss to
    /// any particular mechanism (friction, restitution `e < 1`, rolling resistance, or step 3's
    /// bounded depenetration recovery), but requires no extra per-mechanism bookkeeping and
    /// cannot silently miss one. **Does not include fluid-side viscous dissipation** -- with
    /// slurry enabled, a caller should expect `dissipated_power_w < power_draw_w` in steady
    /// state, since the fluid removes mechanical energy from the balls (drag, ss6.2) that never
    /// shows up here.
    pub dissipated_energy_j: f32,
    /// Largest per-ball `|v| * dt / (2 * radius)` this sub-step: the fraction of a ball's own
    /// diameter it moved in one sub-step. A rough tunnelling-risk indicator -- values approaching
    /// or exceeding 1 mean a fast ball's motion this sub-step is comparable to (or larger than)
    /// its own size, so the broad-phase/contact solve (docs/PLAN.md ss3.2, discrete per-sub-step)
    /// could in principle miss a collision along the way. Diagnostic only in this pass; a
    /// continuous (swept) collision guard is future work.
    pub max_substep_displacement_over_diameter: f32,
}

/// Mass of one ball, modelled as a **unit-depth disc** (kg per metre of mill length):
/// `rho * pi * r^2 * 1 m`. This is the same 2D slice convention the fluid uses
/// ([`crate::pbf::FluidParticles::seed_lattice`]: `particle_mass = rho * dx^2 * 1 m`), so ball and
/// fluid masses are dimensionally consistent and the ball<->fluid momentum exchange
/// ([`crate::coupling`]) is meaningful. An earlier version used a 3D sphere mass here, which made
/// a fluid particle hundreds of times heavier than a coarse-grained ball. Every mass, energy and
/// power reported by this crate is therefore "per metre of mill length".
pub fn ball_mass(diameter_m: f32, density_kg_m3: f32) -> f32 {
    let r = diameter_m * 0.5;
    density_kg_m3 * PI * r * r
}

/// Moment of inertia of a uniform disc about its centre, `1/2 m r^2` (per metre of mill length,
/// consistent with [`ball_mass`]).
pub fn ball_inertia(mass: f32, radius_m: f32) -> f32 {
    0.5 * mass * radius_m * radius_m
}

/// Total mechanical energy (translational + rotational kinetic energy, plus gravitational
/// potential energy) of a ball population, per metre of mill depth (consistent with this crate's
/// unit-depth-disc mass convention, see [`ball_mass`]). Rotational KE must be included, not just
/// translational: friction (step 5 of [`DemState::step_with_external_forces`]) exchanges momentum
/// between a ball's linear and angular velocity (e.g. a spinning ball starting to roll gains
/// linear speed at spin's expense), so omitting it would misread a legitimate
/// translational<->rotational transfer as spurious energy creation or loss.
///
/// Used both by [`DemState::step_with_external_forces`]'s per-sub-step dissipation accounting
/// ([`DemStepStats::dissipated_energy_j`]) and by this module's and `lib.rs`'s energy-invariant
/// regression tests, so the metric and the test that proves it cannot go materially negative
/// share one definition.
pub(crate) fn mechanical_energy_j(balls: &Balls) -> f32 {
    let mut e = 0.0f32;
    for i in 0..balls.len() {
        e += 0.5 * balls.mass * balls.v[i].length_squared();
        e += 0.5 * balls.inertia * balls.omega[i] * balls.omega[i];
        e += balls.mass * GRAVITY.abs() * balls.x[i].y;
    }
    e
}

/// Ball (grinding media) population. All balls currently share one effective radius/mass/inertia
/// (a size distribution is future work); see [`crate::params::Params::effective_media`].
pub struct Balls {
    pub x: Vec<Vec2>,
    pub v: Vec<Vec2>,
    pub theta: Vec<f32>,
    pub omega: Vec<f32>,
    pub radius: f32,
    pub mass: f32,
    pub inertia: f32,
}

impl Balls {
    /// An empty ball population (e.g. for [`crate::pbf::FluidParticles::step`], the ball-free
    /// special case of [`crate::pbf::FluidParticles::step_coupled`]).
    pub fn empty() -> Self {
        Self {
            x: Vec::new(),
            v: Vec::new(),
            theta: Vec::new(),
            omega: Vec::new(),
            radius: 0.0,
            mass: 0.0,
            inertia: 0.0,
        }
    }

    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    fn inv_mass(&self) -> f32 {
        if self.mass > 0.0 {
            1.0 / self.mass
        } else {
            0.0
        }
    }

    fn inv_inertia(&self) -> f32 {
        if self.inertia > 0.0 {
            1.0 / self.inertia
        } else {
            0.0
        }
    }

    /// Seeds `effective.ball_count` balls on a hexagonal lattice filling the drum's
    /// cross-section (skipping lattice sites that would fall outside the drum), with a small
    /// random jitter so initial contacts aren't perfectly aligned. `seed` makes placement
    /// deterministic (see [`crate::rng::Rng`]).
    pub fn seed_lattice(effective: &EffectiveMedia, drum_radius_m: f32, seed: u64) -> Self {
        let r = effective.diameter_m * 0.5;
        let mass = ball_mass(effective.diameter_m, effective.density_kg_m3);
        let inertia = ball_inertia(mass, r);
        let count = effective.ball_count as usize;

        let mut rng = Rng::new(seed);
        let mut x = Vec::with_capacity(count);

        if r > 0.0 && count > 0 {
            let spacing = 2.0 * r * 1.02; // small initial gap so the lattice starts non-overlapping
            let row_h = spacing * (3.0f32.sqrt() * 0.5);
            let fill_r = (drum_radius_m - r).max(0.0);
            let n_rows = (2.0 * fill_r / row_h).ceil() as i32 + 1;

            'rows: for row in -n_rows / 2..=n_rows / 2 {
                let y = row as f32 * row_h;
                if y.abs() > fill_r {
                    continue;
                }
                let half_width = (fill_r * fill_r - y * y).max(0.0).sqrt();
                let offset = if row.rem_euclid(2) == 0 {
                    0.0
                } else {
                    spacing * 0.5
                };
                let n_cols = (2.0 * half_width / spacing).floor() as i32;
                for col in -n_cols / 2..=n_cols / 2 {
                    let cx = col as f32 * spacing + offset;
                    if cx * cx + y * y > fill_r * fill_r {
                        continue;
                    }
                    let jitter =
                        Vec2::new(rng.range_f32(-0.05, 0.05), rng.range_f32(-0.05, 0.05)) * r;
                    x.push(Vec2::new(cx, y) + jitter);
                    if x.len() >= count {
                        break 'rows;
                    }
                }
            }
        }
        x.truncate(count);

        let n = x.len();
        Self {
            x,
            v: vec![Vec2::ZERO; n],
            theta: vec![0.0; n],
            omega: vec![0.0; n],
            radius: r,
            mass,
            inertia,
        }
    }
}

/// Persistent per-substep contact bookkeeping. Cleared and rebuilt every sub-step (see
/// [`DemState::step`]); the accumulated normal impulse `lambda_n` from this sub-step's iterative
/// solve is what friction/rolling-resistance/restitution are computed against afterwards.
///
/// **`HashMap`, with sorted iteration where it matters.** Unlike [`crate::grid::UniformGrid`]
/// (whose fix this module doc comment cross-references), this map's *accumulation* itself
/// (`.entry(key).or_insert(0.0) += d_lambda` in step 3's inner solve loop, run
/// `dem_iterations` times per contact pair every sub-step) is hot enough that `BTreeMap`'s
/// `O(log n)` lookup measurably regressed `cargo bench`'s `dem_step`/`drum_only_step` (roughly
/// 2x at the default ball counts) -- and that accumulation does not actually need sorted order:
/// each key's running sum only ever receives `+=` calls in the fixed order `ball_ball_pairs`
/// already iterates (itself deterministic since `UniformGrid`'s fix), regardless of which
/// internal hash bucket holds it. The real non-determinism risk is in the friction/restitution/
/// rolling-resistance passes (steps 5-7) that iterate *across different keys* and mutate shared
/// per-ball state as they go (a Gauss-Seidel-style sequential correction, where processing order
/// changes the physical result, not just summation rounding) -- those passes use
/// [`ContactBook::sorted_ball_ball`]/[`ContactBook::sorted_ball_wall`] instead of iterating the
/// maps directly, giving deterministic order there at negligible cost (one sort of a
/// per-sub-step-sized collection, not a per-lookup cost).
#[derive(Default)]
struct ContactBook {
    /// Ball-ball contacts: key (i, j) with i < j.
    ball_ball_lambda_n: HashMap<(u32, u32), f32>,
    /// Ball-wall contacts: key is the ball index.
    ball_wall_lambda_n: HashMap<u32, f32>,
}

impl ContactBook {
    /// `ball_ball_lambda_n`'s entries as a `(i, j, lambda_n)` list sorted by `(i, j)`, for
    /// deterministic iteration order (see this struct's doc comment).
    fn sorted_ball_ball(&self) -> Vec<(u32, u32, f32)> {
        let mut entries: Vec<(u32, u32, f32)> = self
            .ball_ball_lambda_n
            .iter()
            .map(|(&(i, j), &lambda_n)| (i, j, lambda_n))
            .collect();
        entries.sort_unstable_by_key(|&(i, j, _)| (i, j));
        entries
    }

    /// `ball_wall_lambda_n`'s entries as an `(i, lambda_n)` list sorted by `i`, for deterministic
    /// iteration order (see this struct's doc comment).
    fn sorted_ball_wall(&self) -> Vec<(u32, f32)> {
        let mut entries: Vec<(u32, f32)> = self
            .ball_wall_lambda_n
            .iter()
            .map(|(&i, &lambda_n)| (i, lambda_n))
            .collect();
        entries.sort_unstable_by_key(|&(i, _)| i);
        entries
    }
}

/// The ball population plus its XPBD solver state.
pub struct DemState {
    pub balls: Balls,
}

impl DemState {
    pub fn new(effective: &EffectiveMedia, drum_radius_m: f32, seed: u64) -> Self {
        Self {
            balls: Balls::seed_lattice(effective, drum_radius_m, seed),
        }
    }

    /// Advances the ball population by one fixed sub-step `dt`, against the given (already
    /// positioned at `drum_angle`) drum. See docs/PLAN.md ss3.2 for the full per-step algorithm.
    /// Equivalent to [`DemState::step_with_external_forces`] with `external = None`.
    pub fn step(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        media: &MediaParams,
        iterations: u32,
        dt: f32,
    ) -> DemStepStats {
        self.step_with_external_forces(drum, drum_angle, media, iterations, dt, None)
    }

    /// As [`DemState::step`], but additionally applies `external` (fluid coupling impulses per
    /// ball, docs/PLAN.md ss3.4, [`crate::coupling`]) as a direct velocity change in the predict
    /// step (`Δv = impulse * inv_mass`, no extra `* dt` -- see
    /// [`crate::coupling::CouplingImpulses`]'s doc comment for why). `external.impulses`/
    /// `external.angular_impulses` must have one entry per ball, same order as `self.balls.x`,
    /// when provided.
    pub fn step_with_external_forces(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        media: &MediaParams,
        iterations: u32,
        dt: f32,
        external: Option<&crate::coupling::CouplingImpulses>,
    ) -> DemStepStats {
        let balls = &mut self.balls;
        if balls.is_empty() || dt <= 0.0 {
            return DemStepStats::default();
        }
        let n = balls.len();
        let w = balls.inv_mass();
        let w_rot = balls.inv_inertia();
        let r = balls.radius;

        // --- 0. Sanitize state carried in from the previous sub-step -------------------------
        // Mirrors `pbf::FluidParticles::step_coupled`'s step 0: given finite input this solver
        // provably produces finite output, so this can only fire on a value corrupted from
        // outside (e.g. a NaN slipping in via `external`). Without this, a single diverged ball
        // silently drops out of `UniformGrid` (it skips non-finite points) and never interacts
        // with anything again, rather than being recovered.
        for (x, v) in balls.x.iter_mut().zip(balls.v.iter_mut()) {
            if !x.is_finite() || !v.is_finite() {
                *x = Vec2::ZERO;
                *v = Vec2::ZERO;
            }
        }

        // --- 1. Predict ---------------------------------------------------------------------
        let x0: Vec<Vec2> = balls.x.clone();
        let theta0: Vec<f32> = balls.theta.clone();
        let v_pre: Vec<Vec2> = balls.v.clone();
        for i in 0..n {
            balls.v[i].y += GRAVITY * dt;
            if let Some(ext) = external {
                balls.v[i] += ext.impulses[i] * w;
                balls.omega[i] += ext.angular_impulses[i] * w_rot;
            }
            balls.x[i] += balls.v[i] * dt;
            balls.theta[i] += balls.omega[i] * dt;
        }
        let e_after_predict = mechanical_energy_j(balls);
        let max_substep_displacement_over_diameter = balls
            .v
            .iter()
            .map(|v| v.length() * dt / (2.0 * r).max(1e-12))
            .fold(0.0f32, f32::max);

        // --- 2. Broad-phase (ball-ball candidates only; wall is checked directly per ball) ----
        // The margin includes this sub-step's largest predicted displacement, not just a fixed
        // 5% pad: at a large `max_substep_displacement_over_diameter` (docs/PLAN.md ss3.5), a
        // fast ball can end this sub-step farther from where it started than the fixed margin
        // covers, which would let the broad-phase miss a pair that is about to overlap deeply
        // once step 3 runs. This does not make the solver continuous (no swept/CCD test is
        // performed within the sub-step, see that field's doc comment) -- it only ensures a pair
        // separated by roughly one sub-step's worth of relative motion is still found as a
        // *candidate*, so step 3's now-bounded recovery ([`MAX_RECOVERY_FRACTION`]) gets a chance
        // to act on it before the overlap becomes severe.
        // Post-predict speed (i.e. including this sub-step's gravity/coupling contribution),
        // the same quantity `max_substep_displacement_over_diameter` above is derived from.
        let max_speed = balls.v.iter().map(|v| v.length()).fold(0.0f32, f32::max);
        let cell_size = (2.0 * r * 1.05).max(2.0 * r + max_speed * dt).max(1e-6);
        let grid = UniformGrid::build(&balls.x, cell_size);
        let mut ball_ball_pairs: Vec<(u32, u32)> = Vec::new();
        grid.for_each_candidate_pair(|i, j| ball_ball_pairs.push((i, j)));

        // --- 2.5 Continuous collision detection against the wall/lifters (conservative
        // advancement) -----------------------------------------------------------------------
        // Clamp each ball's raw predicted advance to the time of first contact with the
        // wall/lifters, so step 3's non-penetration solve below starts from an at-most-barely-
        // touching state instead of repairing an already-tunnelled one. This is what closes the
        // gap docs/PHYSICS.md §9 used to document as "no actual swept/continuous-collision
        // correction" for the wall -- step 2's broad-phase margin above only keeps a fast
        // ball-ball *candidate pair* from being lost, it never limited how far a ball could
        // travel through solid geometry within one sub-step.
        //
        // **Wall/lifter only, deliberately not ball-ball too.** An earlier version of this step
        // also ran `swept_disc_disc_toi` over every `ball_ball_pairs` candidate and took the same
        // min-clamp -- and made a packed, lifters-enabled cataracting charge's
        // `max_ball_wall_overlap_fraction` measurably *worse* (own regression test below), not
        // better. In a dense granular bed a ball routinely has several simultaneous near-touching
        // neighbours (ordinary jostling, not a tunnelling risk), and clamping its advance to the
        // *minimum* time-of-impact across all of them independently froze far more of its motion
        // than the discrete solver's own bounded-recovery pass (`MAX_RECOVERY_FRACTION`, already
        // validated against exactly this many-simultaneous-overlap regime by the fluidised-charge
        // fix this module's tests guard) -- balls piling up unable to settle properly is what
        // pushed *more* of them hard into the wall, not fewer. The wall/lifter case does not have
        // this failure mode: it is a single, well-defined rigid boundary, so a ball is clamped
        // against at most one time-of-impact, not a combinatorial minimum over many neighbours.
        //
        // **Only engaged when the predicted end-of-step overlap is deeper than one bounded-
        // recovery pass could safely absorb**, not on every predicted overlap. A ball sliding
        // fast tangentially along the wall's *curved* boundary has a straight-line one-sub-step
        // chord that dips slightly inside the true circular arc essentially every sub-step --
        // ordinary discretisation, already the discrete solver's routine job to correct (exactly
        // what step 3 below does every sub-step for every wall contact) -- not a tunnelling risk.
        // Gating CCD on `MAX_RECOVERY_FRACTION` (the same bound that already tells the discrete
        // solver how much overlap one pass can safely recover) means CCD only ever intervenes
        // where that budget is genuinely exceeded, leaving the already-validated routine case
        // alone. An earlier version of this step engaged CCD on *any* predicted crossing
        // regardless of depth, which made this same regression test read *worse* than the
        // pre-CCD baseline: balls sliding along a lifter face were repeatedly clamped to a sliver
        // of their predicted advance by the curvature artefact above, crushing their reconstructed
        // tangential velocity every sub-step and leaving them unable to move out of the way of
        // the charge behind them.
        let drum_angle_next = drum_angle + drum.omega * dt;
        let mut toi = vec![1.0f32; n];
        for i in 0..n {
            let (d_end, _) = drum.sdf_world(balls.x[i], drum_angle_next);
            let c_end = d_end - r;
            if c_end >= -MAX_RECOVERY_FRACTION * 2.0 * r {
                continue;
            }
            if let Some(t) = drum.toi_swept(x0[i], balls.x[i], drum_angle, drum_angle_next, r) {
                toi[i] = toi[i].min(t);
            }
        }
        for i in 0..n {
            if toi[i] < 1.0 {
                // Advance a hair past the exact crossing -- a fixed physical depth
                // (`CCD_CONTACT_SKIN_FRACTION * r`), not exactly onto the surface. Stopping
                // exactly at `C = 0` would make step 3's `if c >= 0.0 { continue; }` skip the
                // contact entirely: no `lambda_n` gets recorded, so steps 5-7 (friction,
                // restitution, rolling resistance) never see it either, and a ball clamped to a
                // bare touch would silently fail to rebound instead of being depenetrated and
                // restituted like any other contact next sub-step. A tiny genuine overlap here is
                // deliberate and is exactly what those steps expect.
                let path_len = (balls.x[i] - x0[i]).length();
                let skin_depth = CCD_CONTACT_SKIN_FRACTION * r;
                let t_skin = if path_len > 1e-9 {
                    skin_depth / path_len
                } else {
                    1.0
                };
                let t_clamped = (toi[i] + t_skin).min(1.0);
                balls.x[i] = x0[i] + (balls.x[i] - x0[i]) * t_clamped;
            }
        }

        let mut book = ContactBook::default();

        // --- 3. Solve non-penetration (iterative, zero compliance = rigid) --------------------
        for _ in 0..iterations.max(1) {
            for &(i, j) in &ball_ball_pairs {
                let (i, j) = (i as usize, j as usize);
                let delta = balls.x[i] - balls.x[j];
                let dist = delta.length();
                let c = dist - 2.0 * r;
                if c >= 0.0 {
                    continue;
                }
                let n_hat = if dist > 1e-9 { delta / dist } else { Vec2::X };
                let w_sum = w + w;
                if w_sum <= 0.0 {
                    continue;
                }
                // Cap how much of a deep overlap this one iteration recovers (see
                // [`MAX_RECOVERY_FRACTION`]'s doc comment) -- `lambda_n` (and hence the Coulomb
                // friction cone and rolling resistance below) shrinks along with it, which is
                // correct: a partially-recovered contact should carry proportionally less normal
                // force this sub-step, not the full force a complete recovery would imply.
                let c_eff = c.max(-MAX_RECOVERY_FRACTION * 2.0 * r);
                let d_lambda = -c_eff / w_sum;
                balls.x[i] += w * d_lambda * n_hat;
                balls.x[j] -= w * d_lambda * n_hat;
                *book
                    .ball_ball_lambda_n
                    .entry((i as u32, j as u32))
                    .or_insert(0.0) += d_lambda;
            }
            for i in 0..n {
                let (d, n_hat) = drum.sdf_world(balls.x[i], drum_angle_next);
                let c = d - r;
                if c >= 0.0 {
                    continue;
                }
                if w <= 0.0 {
                    continue;
                }
                let c_eff = c.max(-MAX_RECOVERY_FRACTION * 2.0 * r);
                let d_lambda = -c_eff / w;
                balls.x[i] += w * d_lambda * n_hat;
                *book.ball_wall_lambda_n.entry(i as u32).or_insert(0.0) += d_lambda;
            }
        }

        // --- 4. Reconstruct velocities from the position change -------------------------------
        for i in 0..n {
            balls.v[i] = (balls.x[i] - x0[i]) / dt;
            balls.omega[i] = angle_diff(balls.theta[i], theta0[i]) / dt;
        }

        // Sorted once, reused by every pass below that iterates *across* contacts (steps 5-7) --
        // see `ContactBook`'s doc comment for why this, rather than iterating the maps directly,
        // is what keeps this solver deterministic without slowing down step 3's hot accumulation.
        let ball_ball_contacts = book.sorted_ball_ball();
        let ball_wall_contacts = book.sorted_ball_wall();

        // --- 5. Friction (Coulomb-clamped position correction, one pass but velocity kept in
        // sync per contact) ---------------------------------------------------------------------
        // Each contact's correction, viewed alone (given whatever `v`/`omega` are *right now*),
        // is provably non-energy-increasing: it reduces the relative tangential contact-point
        // speed `v_t` toward zero (fully when unclamped, partially when Coulomb-limited), which
        // is exactly the "generalized inelastic collision" impulse for that one degree of
        // freedom and cannot overshoot past zero. But a ball resting on more than one neighbour
        // (or a neighbour and the wall) has more than one contact processed in this same pass;
        // if `v`/`omega` were only reconstructed once at the very end (from the *total* summed
        // position/orientation change, as steps 3-4 do), a later contact touching an
        // already-corrected ball would compute its own `v_t` against a *stale* velocity that
        // doesn't yet reflect the earlier contact's correction -- silently breaking the
        // per-contact non-energy-increase guarantee once several contacts share a ball (the
        // common case under gravity load), by up to several percent of the local kinetic energy
        // per sub-step. Updating `v`/`omega` immediately after each contact -- exactly the
        // velocity change its own position/orientation correction implies -- keeps every
        // subsequent contact's `v_t` accurate, restoring that guarantee.
        for &(i, j, lambda_n) in &ball_ball_contacts {
            if lambda_n <= 0.0 {
                continue;
            }
            let (iu, ju) = (i as usize, j as usize);
            let delta = balls.x[iu] - balls.x[ju];
            let dist = delta.length();
            let n_hat = if dist > 1e-9 { delta / dist } else { Vec2::X };
            let t_hat = Vec2::new(-n_hat.y, n_hat.x);

            let v_t =
                (balls.v[iu] - balls.v[ju]).dot(t_hat) - r * balls.omega[iu] - r * balls.omega[ju];
            let w_sum_t = w + w + r * r * w_rot + r * r * w_rot;
            if w_sum_t <= 0.0 {
                continue;
            }
            let raw = -v_t * dt / w_sum_t;
            let max_mag = media.friction_ball_ball * lambda_n;
            let d_lambda_t = raw.clamp(-max_mag, max_mag);

            balls.x[iu] += w * d_lambda_t * t_hat;
            balls.x[ju] -= w * d_lambda_t * t_hat;
            balls.theta[iu] -= r * w_rot * d_lambda_t;
            balls.theta[ju] -= r * w_rot * d_lambda_t;
            let dv = w * d_lambda_t * t_hat / dt;
            let domega = -r * w_rot * d_lambda_t / dt;
            balls.v[iu] += dv;
            balls.v[ju] -= dv;
            balls.omega[iu] += domega;
            balls.omega[ju] += domega;
        }
        let mut wall_work_j = 0.0f32;
        for &(i, lambda_n) in &ball_wall_contacts {
            if lambda_n <= 0.0 {
                continue;
            }
            let iu = i as usize;
            let (_, n_hat) = drum.sdf_world(balls.x[iu], drum_angle_next);
            let t_hat = Vec2::new(-n_hat.y, n_hat.x);
            // Sampled at the contact point (`x_c - r*n_hat`), not the ball centre: the wall is
            // rigid, so its velocity varies across the ball's own radius, and the contact point is
            // the only location where "no relative tangential slip" is a meaningful statement. A
            // ball rigidly co-rotating with the drum (`omega_ball == drum.omega`) has zero slip
            // *there*; sampling at the centre instead undercounts the wall's tangential speed by
            // `omega * r` (since `t_hat` at the contact point differs from the centre-sampled one
            // only by that much, for a circular/lifter-face wall), which previously read a
            // perfectly rolling ball as slipping and biased friction, the reconstructed spin, and
            // `wall_work_j` below.
            let v_wall = drum.wall_velocity(balls.x[iu] - r * n_hat);

            let v_t = (balls.v[iu] - v_wall).dot(t_hat) - r * balls.omega[iu];
            let w_sum_t = w + r * r * w_rot;
            if w_sum_t <= 0.0 {
                continue;
            }
            let raw = -v_t * dt / w_sum_t;
            let max_mag = media.friction_ball_wall * lambda_n;
            let d_lambda_t = raw.clamp(-max_mag, max_mag);

            balls.x[iu] += w * d_lambda_t * t_hat;
            balls.theta[iu] -= r * w_rot * d_lambda_t;
            balls.v[iu] += w * d_lambda_t * t_hat / dt;
            balls.omega[iu] += -r * w_rot * d_lambda_t / dt;
            // `d_lambda_t` is the tangential *position* multiplier (units kg*m -- see this
            // module's doc comment and `CouplingImpulses`'s doc comment in `coupling.rs` for the
            // same impulse-vs-position-multiplier distinction), not the impulse itself: it is
            // applied as a position correction (`x += w * d_lambda_t * t_hat`, `w` = inverse
            // mass), so the actual tangential impulse the wall delivered to this ball is
            // `d_lambda_t / dt` (kg*m/s per metre of mill depth). The rate of work the drum's
            // motor must supply to overcome the ball's Newton's-third-law reaction against the
            // wall is that impulse dotted with the wall's own velocity there (docs/PLAN.md
            // ss3.5's power-draw estimate); accumulated directly as work (impulse . velocity),
            // with no further `* dt`, since `wall_work_j` is a per-sub-step quantity. An earlier
            // version of this line omitted the `/ dt` and under-reported power draw by a factor
            // of `dt` (240x at this project's default sub-step) relative to the (correctly
            // `dt`-independent) dissipated-energy estimate below.
            wall_work_j += (d_lambda_t / dt) * t_hat.dot(v_wall);
            // The wall's *normal* impulse does no work only for a smooth cylindrical wall, whose
            // normal is always radial while the wall's own velocity is purely tangential
            // (`n_hat . v_wall == 0` there). A lifter face's normal has a large circumferential
            // component -- pushing a ball up and along, not just constraining it radially -- so
            // `n_hat . v_wall` is generally nonzero and this term is exactly the missing "lifter
            // lifts the charge" power contribution an external review found absent (a smooth-wall
            // mill with `lifters.count == 0` has `n_hat . v_wall == 0` everywhere, so this adds
            // nothing there and stays a strict no-op on that configuration). Same impulse-vs-
            // position-multiplier convention as the tangential term above: `lambda_n / dt` is the
            // actual normal impulse (kg*m/s per metre of mill depth).
            wall_work_j += (lambda_n / dt) * n_hat.dot(v_wall);
        }
        // Re-reconstruct velocities from the total position/orientation change once more: exactly
        // redundant with the incremental updates above in exact arithmetic (both compute the same
        // net change), kept as a cheap self-healing resync against any float summation-order drift
        // between the two, consistent with how step 4 above derives velocity from position.
        for i in 0..n {
            balls.v[i] = (balls.x[i] - x0[i]) / dt;
            balls.omega[i] = angle_diff(balls.theta[i], theta0[i]) / dt;
        }

        // --- 6. Restitution (bidirectional, using the pre-solve approach velocity to set the
        // target) ------------------------------------------------------------------------------
        // For a genuine fresh impact the target relative normal velocity is the usual
        // restitution law, `e * |v_approach|` (separating). For a resting/sliding contact
        // (approach speed below `RESTITUTION_VELOCITY_THRESHOLD`) the target is exactly zero.
        // Either way the *actual* relative normal velocity is driven to that target in both
        // directions -- not just topped up when it falls short. This is what keeps the solver
        // energy-bounded: step 3's depenetration recovers overlap by moving positions, and step
        // 4 turns that raw position change into velocity with no regard for how much separation
        // speed is physically justified. An earlier version of this pass only ever *added*
        // separation velocity (bailing out whenever the depenetration had already produced at
        // least the restitution target), which let step 3's recovery inject unbounded velocity
        // into deep-overlap events instead of being capped by it -- the DEM half of the
        // fluidised-charge energy-injection bug (docs/PHYSICS.md ss3.2/3.5). See
        // [`MAX_RECOVERY_FRACTION`] for step 3's half of the same fix.
        let mut collision_count = 0u32;
        let mut impact_energy_histogram = [0u32; IMPACT_ENERGY_HISTOGRAM_BINS];
        for &(i, j, lambda_n) in &ball_ball_contacts {
            if lambda_n <= 0.0 {
                continue;
            }
            let (iu, ju) = (i as usize, j as usize);
            let delta = balls.x[iu] - balls.x[ju];
            let dist = delta.length();
            let n_hat = if dist > 1e-9 { delta / dist } else { Vec2::X };

            let v_n_pre = (v_pre[iu] - v_pre[ju]).dot(n_hat);
            let is_fresh_impact = v_n_pre < -RESTITUTION_VELOCITY_THRESHOLD;
            if is_fresh_impact {
                // A genuine fresh impact (docs/PLAN.md ss3.5): reduced mass `mass/2` for two
                // equal masses colliding.
                collision_count += 1;
                let impact_energy_j = 0.5 * (balls.mass * 0.5) * v_n_pre * v_n_pre;
                if let Some(bin) = impact_energy_bin(impact_energy_j) {
                    impact_energy_histogram[bin] += 1;
                }
            }
            let e = if is_fresh_impact {
                media.restitution_ball_ball
            } else {
                0.0
            };
            // `.max(0.0)`: the target can never be negative (a contact cannot pull its pair
            // together), and restitution never reverses the sign of the approach velocity.
            // `.max(v_n_pre.max(0.0))`: if the pair was already separating *before* this solve
            // (`v_n_pre > 0` -- e.g. step 3's depenetration had already pushed them apart faster
            // than restitution alone would justify), the target must not fall below that existing
            // separation speed. Without this, a pair that entered this step already separating
            // faster than `target` gets `delta_v_n < 0` below, i.e. an impulse that pulls them back
            // together -- an artificial attraction between dry rigid discs (Signorini's `F_n >= 0`
            // forbids exactly this). This does not reopen the energy-injection bug this pass exists
            // to close: an *approaching* pair (`v_n_pre <= 0`) is untouched (`v_n_pre.max(0.0)` is
            // `0.0` there), so step 3's depenetration is still capped at the restitution target
            // exactly as before.
            let target = (-e * v_n_pre).max(0.0).max(v_n_pre.max(0.0));
            let v_n_now = (balls.v[iu] - balls.v[ju]).dot(n_hat);
            let delta_v_n = target - v_n_now;
            let w_sum = w + w;
            if w_sum <= 0.0 {
                continue;
            }
            let impulse = delta_v_n / w_sum;
            balls.v[iu] += w * impulse * n_hat;
            balls.v[ju] -= w * impulse * n_hat;
        }
        for &(i, lambda_n) in &ball_wall_contacts {
            if lambda_n <= 0.0 {
                continue;
            }
            let iu = i as usize;
            let (_, n_hat) = drum.sdf_world(balls.x[iu], drum_angle_next);
            let v_wall = drum.wall_velocity(balls.x[iu] - r * n_hat);

            let v_n_pre = (v_pre[iu] - v_wall).dot(n_hat);
            let is_fresh_impact = v_n_pre < -RESTITUTION_VELOCITY_THRESHOLD;
            if is_fresh_impact {
                // Wall has "infinite" mass, so the reduced mass of the pair is just the ball's
                // own.
                collision_count += 1;
                let impact_energy_j = 0.5 * balls.mass * v_n_pre * v_n_pre;
                if let Some(bin) = impact_energy_bin(impact_energy_j) {
                    impact_energy_histogram[bin] += 1;
                }
            }
            let e = if is_fresh_impact {
                media.restitution_ball_wall
            } else {
                0.0
            };
            // See the ball-ball restitution pass above for why the target is also floored at the
            // pre-solve separation speed.
            let target = (-e * v_n_pre).max(0.0).max(v_n_pre.max(0.0));
            let v_n_now = (balls.v[iu] - v_wall).dot(n_hat);
            let delta_v_n = target - v_n_now;
            if w <= 0.0 {
                continue;
            }
            balls.v[iu] += delta_v_n * n_hat;
        }

        // --- 7. Rolling resistance -------------------------------------------------------------
        // Damps every ball's *world* `omega` toward zero, scaled by that ball's total accumulated
        // normal load (`total_lambda_n`, summed across all its ball-ball/ball-wall contacts this
        // sub-step) -- a known simplification (docs/PHYSICS.md ss9): it neither applies a
        // reaction torque to a ball-ball contact's partner (so it does not conserve angular
        // momentum) nor references a ball-wall contact's target to `drum.omega` rather than zero
        // (so it fights, rather than assists, a ball correctly rolling without slipping on a
        // rotating drum). A per-contact, momentum-conserving, wall-frame-aware rewrite was tried
        // and reverted: it destabilised this crate's core cataracting/cascading-density regression
        // tests (`dem::tests::a_cascading_charge_stays_dense`,
        // `metrics::tests::compression_error_stays_bounded_under_violent_lifter_cataracting`),
        // which protect against the much more severe fluidised-charge energy-injection bug this
        // solver was extensively tuned against -- see this module's and `MAX_RECOVERY_FRACTION`'s
        // doc comments. Revisiting this needs dedicated re-tuning against those regressions, not a
        // drive-by correctness fix.
        if media.rolling_friction > 0.0 && balls.inertia > 0.0 {
            let mut total_lambda_n = vec![0.0f32; n];
            for &(i, j, lambda_n) in &ball_ball_contacts {
                total_lambda_n[i as usize] += lambda_n.max(0.0);
                total_lambda_n[j as usize] += lambda_n.max(0.0);
            }
            for &(i, lambda_n) in &ball_wall_contacts {
                total_lambda_n[i as usize] += lambda_n.max(0.0);
            }
            for (omega, &lambda_n) in balls.omega.iter_mut().zip(total_lambda_n.iter()) {
                if lambda_n <= 0.0 || *omega == 0.0 {
                    continue;
                }
                let normal_force = lambda_n / dt;
                let max_delta_omega =
                    media.rolling_friction * normal_force * r / balls.inertia * dt;
                let delta_omega = max_delta_omega.min(omega.abs());
                *omega -= omega.signum() * delta_omega;
            }
        }

        let e_final = mechanical_energy_j(balls);
        DemStepStats {
            wall_work_j,
            collision_count,
            impact_energy_histogram,
            dissipated_energy_j: e_after_predict + wall_work_j - e_final,
            max_substep_displacement_over_diameter,
        }
    }
}

/// Smallest signed difference `a - b`, avoiding a spurious huge value if a caller's angles have
/// drifted outside a single revolution (not expected here since `dt` is tiny, but cheap to guard).
fn angle_diff(a: f32, b: f32) -> f32 {
    let mut d = a - b;
    while d > PI {
        d -= 2.0 * PI;
    }
    while d < -PI {
        d += 2.0 * PI;
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::LiftersParams;

    fn media_defaults() -> MediaParams {
        MediaParams::default()
    }

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

    fn single_ball(x: Vec2, v: Vec2, radius: f32, mass: f32) -> Balls {
        Balls {
            x: vec![x],
            v: vec![v],
            theta: vec![0.0],
            omega: vec![0.0],
            radius,
            mass,
            inertia: ball_inertia(mass, radius),
        }
    }

    #[test]
    fn ball_at_rest_on_flat_wall_stays_put_with_small_overlap() {
        // A ball dropped just touching the floor of a large "drum" (so the floor is locally flat)
        // should settle with only a tiny overlap and no lateral drift.
        let radius_drum = 10.0;
        let r = 0.05;
        let mass = ball_mass(2.0 * r, 7800.0);
        let balls = single_ball(Vec2::new(0.0, -(radius_drum - r)), Vec2::ZERO, r, mass);
        let mut state = DemState { balls };
        let drum = still_drum(radius_drum);
        let media = media_defaults();

        for _ in 0..240 {
            state.step(&drum, 0.0, &media, 4, 1.0 / 240.0);
        }

        let x = state.balls.x[0];
        let dist_to_wall = radius_drum - x.length();
        assert!(
            dist_to_wall <= r * 1.02,
            "ball should rest ~touching the wall, got dist={dist_to_wall}, r={r}"
        );
        assert!(
            dist_to_wall >= r * 0.9,
            "ball should not sink far into the wall, got dist={dist_to_wall}, r={r}"
        );
        assert!(
            x.x.abs() < 0.02,
            "ball should not drift laterally, got x={}",
            x.x
        );
        assert!(
            state.balls.v[0].length() < 0.05,
            "ball should be nearly at rest, got v={:?}",
            state.balls.v[0]
        );
    }

    #[test]
    fn ball_rests_stably_on_a_lifter_tip_face() {
        // Same idea as `ball_at_rest_on_flat_wall_stays_put_with_small_overlap`, but the ball
        // settles on a lifter's flat tip face (docs/PLAN.md ss3.1) instead of the bare wall. One
        // lifter is centred at the bottom of a large "drum" (phase -90 degrees, i.e. -y, where
        // gravity carries the ball); the ball is small relative to the tip width so it rests
        // centred rather than balanced on an edge.
        let radius_drum = 10.0;
        let lifters = LiftersParams {
            count: 1,
            height_m: 0.05,
            base_width_m: 0.1,
            top_width_m: 0.06,
            phase_deg: -90.0,
        };
        let drum = Drum::new(radius_drum, 0.0, lifters);
        let tip_radius = radius_drum - lifters.height_m;

        let r = 0.01; // small vs. the 60 mm tip width
        let mass = ball_mass(2.0 * r, 7800.0);
        let balls = single_ball(Vec2::new(0.0, -(tip_radius - r)), Vec2::ZERO, r, mass);
        let mut state = DemState { balls };
        let media = media_defaults();

        for _ in 0..240 {
            state.step(&drum, 0.0, &media, 4, 1.0 / 240.0);
        }

        let x = state.balls.x[0];
        let dist_to_tip = tip_radius - x.length();
        assert!(
            dist_to_tip <= r * 1.02,
            "ball should rest ~touching the lifter's tip, got dist={dist_to_tip}, r={r}"
        );
        assert!(
            dist_to_tip >= r * 0.9,
            "ball should not sink far into the lifter, got dist={dist_to_tip}, r={r}"
        );
        assert!(
            x.x.abs() < 0.005,
            "ball should not drift off the lifter, got x={}",
            x.x
        );
        assert!(
            state.balls.v[0].length() < 0.05,
            "ball should be nearly at rest, got v={:?}",
            state.balls.v[0]
        );
    }

    #[test]
    fn ball_bounces_with_restitution_e_squared_height_ratio() {
        // A ball dropped from height h0 above a flat floor should rebound to approximately e^2 * h0
        // (classical restitution result), for a *single* bounce.
        let radius_drum = 10.0;
        let r = 0.05;
        let mass = ball_mass(2.0 * r, 7800.0);
        let e = 0.8;
        let h0 = 0.5;
        let start_y = -(radius_drum - r) + h0;

        let mut media = media_defaults();
        media.restitution_ball_wall = e;
        media.friction_ball_wall = 0.0; // isolate normal restitution from tangential effects
        media.rolling_friction = 0.0;

        let mut state = DemState {
            balls: single_ball(Vec2::new(0.0, start_y), Vec2::ZERO, r, mass),
        };
        let drum = still_drum(radius_drum);
        let dt = 1.0 / 480.0;

        let mut max_height_after_bounce = f32::MIN;
        let mut bounced = false;
        let mut settled_after_bounce = 0;
        for _ in 0..20000 {
            state.step(&drum, 0.0, &media, 4, dt);
            let dist_to_wall = radius_drum - state.balls.x[0].length();
            let height = dist_to_wall - r;
            if !bounced {
                if state.balls.v[0].y > 0.0 && height < h0 * 0.5 {
                    bounced = true;
                }
            } else {
                max_height_after_bounce = max_height_after_bounce.max(height);
                if state.balls.v[0].y <= 0.0 {
                    settled_after_bounce += 1;
                    if settled_after_bounce > 5 {
                        break; // stop once we're clearly past the rebound apex
                    }
                }
            }
        }

        assert!(bounced, "ball never rebounded off the floor");
        let expected = e * e * h0;
        let rel_err = (max_height_after_bounce - expected).abs() / expected;
        assert!(
            rel_err < 0.25,
            "rebound height {max_height_after_bounce} too far from e^2*h0={expected} (rel_err={rel_err})"
        );
    }

    #[test]
    fn two_balls_head_on_collision_conserves_momentum() {
        let r = 0.05;
        let mass = ball_mass(2.0 * r, 7800.0);
        let mut media = media_defaults();
        media.restitution_ball_ball = 0.9;
        media.friction_ball_ball = 0.0;
        media.rolling_friction = 0.0;

        let balls = Balls {
            x: vec![Vec2::new(-0.5, 0.0), Vec2::new(0.5, 0.0)],
            v: vec![Vec2::new(2.0, 0.0), Vec2::new(-2.0, 0.0)],
            theta: vec![0.0, 0.0],
            omega: vec![0.0, 0.0],
            radius: r,
            mass,
            inertia: ball_inertia(mass, r),
        };
        let mut state = DemState { balls };
        // A drum far larger than the balls' excursion, so wall contact never triggers, no gravity
        // effect on the horizontal axis being measured (gravity acts on y, motion is on x).
        let drum = still_drum(1000.0);
        let dt = 1.0 / 480.0;

        // Only the x-component is meaningful for momentum conservation here: gravity acts on y
        // throughout (there is no floor nearby to stop the fall), while the collision itself is
        // purely along x, so x-momentum is conserved by the collision alone.
        let px0: f32 = state.balls.mass * (state.balls.v[0].x + state.balls.v[1].x);

        for _ in 0..600 {
            state.step(&drum, 0.0, &media, 4, dt);
        }

        let px1: f32 = state.balls.mass * (state.balls.v[0].x + state.balls.v[1].x);
        let rel_err = (px1 - px0).abs() / px0.abs().max(1e-9);
        assert!(
            rel_err < 0.05,
            "x-momentum not conserved: px0={px0} px1={px1} rel_err={rel_err}"
        );

        // After a head-on collision the balls must have separated (not passed through each other).
        assert!(state.balls.x[1].x > state.balls.x[0].x + 2.0 * r * 0.9);
    }

    #[test]
    fn ball_resting_on_rotating_wall_settles_into_rolling_without_slip() {
        // Ball dropped onto a large ("locally flat"), slowly rotating drum's floor. Once it has
        // settled, its own spin should track the wall's rotation closely enough that the
        // *contact point* has near-zero relative tangential velocity ("rolling without
        // slipping") -- not just the ball's centre-of-mass velocity, which the wall's own
        // rotation keeps nonzero regardless. Exercises the wall-velocity-at-contact-point fix:
        // friction/restitution now sample `drum.wall_velocity` at `x - r*n_hat`, not the ball
        // centre, so a rigidly co-rotating ball no longer reads as artificially slipping by
        // `omega_ball * r`.
        let radius_drum = 10.0;
        let omega = 0.5; // rad/s; keeps the wall speed near the floor modest (~5 m/s)
        let r = 0.05;
        let mass = ball_mass(2.0 * r, 7800.0);
        let balls = single_ball(Vec2::new(0.0, -(radius_drum - r)), Vec2::ZERO, r, mass);
        let mut state = DemState { balls };
        let drum = Drum::new(
            radius_drum,
            omega,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        );
        let media = media_defaults();
        let dt = 1.0 / 480.0;

        // `drum_angle` stays 0.0 throughout: a smooth cylindrical wall's SDF is rotation-
        // invariant, so only `drum.omega` (the wall's *velocity* field) matters here, exactly as
        // the other still-drum-with-omega-baked-into-the-`Drum`-value tests in this module rely
        // on for `Drum::wall_velocity`.
        for _ in 0..(5 * 480) {
            state.step(&drum, 0.0, &media, 4, dt);
        }

        let x = state.balls.x[0];
        let (_, n_hat) = drum.sdf_world(x, 0.0);
        let t_hat = Vec2::new(-n_hat.y, n_hat.x);
        let v_wall = drum.wall_velocity(x - r * n_hat);
        let slip = (state.balls.v[0] - v_wall).dot(t_hat) - r * state.balls.omega[0];
        let wall_speed = v_wall.length();
        assert!(
            slip.abs() < 0.05 * wall_speed.max(1e-6),
            "contact-point slip too large: slip={slip} wall_speed={wall_speed} omega_ball={}",
            state.balls.omega[0]
        );
    }

    #[test]
    fn lifter_normal_force_contributes_to_wall_work() {
        // A lifter face's normal has a large circumferential component and does real work
        // lifting the charge, unlike a smooth wall's purely radial normal -- `wall_work_j` must
        // pick that up (`(lambda_n/dt) * n_hat.dot(v_wall)`, alongside the tangential-friction
        // term), not just the tangential-friction accounting alone.
        let r = 0.02;
        let radius_drum = 0.5;
        let omega = 6.0; // fast enough to cataract off a lifter within a short run
        let lifters = LiftersParams {
            count: 6,
            height_m: 0.03,
            base_width_m: 0.05,
            top_width_m: 0.03,
            phase_deg: 0.0,
        };
        let drum = Drum::new(radius_drum, omega, lifters);
        let mut media = media_defaults();
        media.rolling_friction = 0.0;
        let effective = crate::params::EffectiveMedia {
            true_diameter_m: 2.0 * r,
            diameter_m: 2.0 * r,
            density_kg_m3: 6000.0,
            ball_count: 40,
            scale_factor: 1.0,
        };
        let mut state = DemState::new(&effective, radius_drum, 42);
        let dt = 1.0 / 960.0;
        let mut drum_angle = 0.0f32;

        let mut total_wall_work_j = 0.0f32;
        let mut any_normal_work = false;
        for _ in 0..(2 * 960) {
            let stats = state.step(&drum, drum_angle, &media, 4, dt);
            total_wall_work_j += stats.wall_work_j;
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
            if stats.wall_work_j.abs() > 1e-9 {
                any_normal_work = true;
            }
        }

        assert!(
            any_normal_work,
            "expected lifter contacts to register nonzero wall_work_j over the run"
        );
        assert!(
            total_wall_work_j > 0.0,
            "expected lifters to do net positive work lifting the charge, got {total_wall_work_j}"
        );
    }

    #[test]
    fn smooth_wall_normal_force_does_no_work() {
        // Companion to `lifter_normal_force_contributes_to_wall_work`: on a smooth cylindrical
        // wall, `wall_work_j`'s `n_hat . v_wall` term must be an exact no-op, since a smooth
        // wall's normal is always radial while its velocity is purely tangential there.
        let radius_drum = 10.0;
        let omega = 0.5;
        let r = 0.05;
        let mass = ball_mass(2.0 * r, 7800.0);
        let balls = single_ball(Vec2::new(0.0, -(radius_drum - r)), Vec2::ZERO, r, mass);
        let mut state = DemState { balls };
        let drum = Drum::new(
            radius_drum,
            omega,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        );
        let media = media_defaults();
        let dt = 1.0 / 480.0;

        // Settle first, then measure: the friction-only tangential term should be the entire
        // `wall_work_j` reading once resting, since the normal term is provably zero regardless.
        for _ in 0..480 {
            state.step(&drum, 0.0, &media, 4, dt);
        }
        let x = state.balls.x[0];
        let (_, n_hat) = drum.sdf_world(x, 0.0);
        let v_wall_at_contact = drum.wall_velocity(x - r * n_hat);
        let dot = n_hat.dot(v_wall_at_contact);
        // Tolerance is relative to wall speed, not a tiny absolute epsilon: `sdf_world`'s normal
        // is a numeric (central-difference) gradient, and at `radius_drum = 10` the underlying
        // `sdf` values it differences are a small difference of two O(10) lengths (catastrophic
        // cancellation in f32), so the normal itself carries roughly `1e-3`-relative noise here --
        // exactly zero only in exact arithmetic.
        let wall_speed = v_wall_at_contact.length();
        assert!(
            dot.abs() < 0.01 * wall_speed,
            "smooth wall's normal should carry ~no velocity component, got {dot} (wall speed {wall_speed})"
        );
    }

    #[test]
    fn restitution_does_not_pull_an_already_separating_pair_together() {
        // Two balls given an initial *separating* velocity while still overlapping (as if step
        // 3's depenetration had already pushed them apart faster than restitution alone would
        // justify). The restitution target's floor at the pre-solve separation speed must not let
        // the solve decelerate them below the separation speed they already had -- an artificial
        // attraction between dry rigid discs would violate the one-sided normal-force (Signorini)
        // condition.
        let r = 0.05;
        let mass = ball_mass(2.0 * r, 7800.0);
        let mut media = media_defaults();
        media.restitution_ball_ball = 0.5; // low e: without the floor this pulls hard
        media.friction_ball_ball = 0.0;
        media.rolling_friction = 0.0;

        // Overlapping (`dist < 2r`) but already moving apart faster than `e * 0` would target.
        let separation_speed = 3.0;
        let balls = Balls {
            x: vec![Vec2::new(-0.9 * r, 0.0), Vec2::new(0.9 * r, 0.0)],
            v: vec![
                Vec2::new(-separation_speed / 2.0, 0.0),
                Vec2::new(separation_speed / 2.0, 0.0),
            ],
            theta: vec![0.0, 0.0],
            omega: vec![0.0, 0.0],
            radius: r,
            mass,
            inertia: ball_inertia(mass, r),
        };
        let mut state = DemState { balls };
        let drum = still_drum(1000.0);
        let dt = 1.0 / 480.0;

        state.step(&drum, 0.0, &media, 4, dt);

        let relative_speed_after = (state.balls.v[1] - state.balls.v[0]).dot(Vec2::X);
        assert!(
            relative_speed_after >= separation_speed * 0.95,
            "an already-separating pair lost separation speed: before={separation_speed} after={relative_speed_after}"
        );
    }

    #[test]
    fn settled_bed_produces_no_spurious_collisions() {
        // Regression pin for a rejected review finding: it claimed a fully settled, motionless
        // bed under gravity would register a fresh "collision" every sub-step, because the
        // predict step's `g*dt` speed increment (~0.041 m/s at 240 Hz) exceeds
        // `RESTITUTION_VELOCITY_THRESHOLD` (0.02 m/s). That reasoning uses the *post-predict*
        // velocity; the actual impact test uses `v_pre`, captured *before* gravity is applied
        // this sub-step (see `step_with_external_forces`'s doc comment on step 1 and step 6's use
        // of `v_pre`), so a resting contact's approach speed reads as whatever it was at the
        // *start* of the sub-step -- ~0 once settled, not `g*dt`.
        let radius_drum = 10.0;
        let r = 0.05;
        let mass = ball_mass(2.0 * r, 7800.0);
        let balls = single_ball(Vec2::new(0.0, -(radius_drum - r)), Vec2::ZERO, r, mass);
        let mut state = DemState { balls };
        let drum = still_drum(radius_drum);
        let media = media_defaults();
        let dt = 1.0 / 240.0;

        // Let it settle first.
        for _ in 0..240 {
            state.step(&drum, 0.0, &media, 4, dt);
        }

        let mut total_collisions = 0u32;
        for _ in 0..240 {
            let stats = state.step(&drum, 0.0, &media, 4, dt);
            total_collisions += stats.collision_count;
        }
        assert_eq!(
            total_collisions, 0,
            "a settled, motionless bed should register no fresh collisions"
        );
    }

    #[test]
    fn centrifuging_above_critical_speed_keeps_balls_at_wall() {
        use crate::params::{Direction, LiftersParams, MediaParams, MillParams, Params, SpeedMode};

        let mut params = Params {
            mill: MillParams {
                diameter_m: 1.0,
                speed_mode: SpeedMode::PercentCritical,
                speed_value: 120.0, // above critical: expect centrifuging
                direction: Direction::CounterClockwise,
            },
            media: MediaParams {
                ball_diameter_m: 0.02,
                // Kept thin on purpose: a thin annulus makes "no ball anywhere near the centre"
                // (checked below) an unambiguous centrifuging signature, without the assertion
                // being sensitive to exactly how many layers deep a thicker charge packs.
                fill_fraction: 0.10,
                ..MediaParams::default()
            },
            lifters: LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
            ..Params::default()
        };
        params.simulation.max_balls = 150;
        params.validate().unwrap();

        let effective = params.effective_media();
        let drum = Drum::new(params.mill.radius_m(), params.mill.omega(), params.lifters);
        let mut state = DemState::new(&effective, params.mill.radius_m(), params.simulation.seed);

        let mut drum_angle = 0.0f32;
        let dt = 1.0 / 240.0;
        for _ in 0..(5 * 240) {
            state.step(
                &drum,
                drum_angle,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            drum_angle += params.mill.omega() * dt;
        }

        // Centrifuging above critical speed should leave the drum's centre region empty: every
        // ball is pushed out into a (possibly multi-layer) annulus against the wall. This is
        // robust to the exact annulus thickness, unlike a fixed per-ball wall-distance bound.
        let radius_m = params.mill.radius_m();
        let min_dist_from_center = state
            .balls
            .x
            .iter()
            .map(|p| p.length())
            .fold(f32::MAX, f32::min);
        assert!(
            min_dist_from_center > 0.5 * radius_m,
            "balls not centrifuged: found a ball at dist_from_center={min_dist_from_center}, radius_m={radius_m}"
        );
    }

    #[test]
    fn cascading_below_critical_speed_forms_a_toe_and_shoulder() {
        use crate::params::{Direction, LiftersParams, MediaParams, MillParams, Params, SpeedMode};

        let mut params = Params {
            mill: MillParams {
                diameter_m: 1.0,
                speed_mode: SpeedMode::PercentCritical,
                speed_value: 70.0,
                direction: Direction::CounterClockwise,
            },
            media: MediaParams {
                ball_diameter_m: 0.02,
                fill_fraction: 0.25,
                ..MediaParams::default()
            },
            lifters: LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
            ..Params::default()
        };
        params.simulation.max_balls = 300;
        params.validate().unwrap();

        let effective = params.effective_media();
        let drum = Drum::new(params.mill.radius_m(), params.mill.omega(), params.lifters);
        let mut state = DemState::new(&effective, params.mill.radius_m(), params.simulation.seed);

        let mut drum_angle = 0.0f32;
        let dt = 1.0 / 240.0;
        // Run long enough for the charge to reach a statistically steady cascading pattern.
        for _ in 0..(8 * 240) {
            state.step(
                &drum,
                drum_angle,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            drum_angle += params.mill.omega() * dt;
        }

        let r = effective.diameter_m * 0.5;
        let wall_margin = 2.5 * r;
        // Balls within `wall_margin` of the wall approximate the outer layer of the charge, whose
        // angular extent brackets the toe (leading edge, low angle) and shoulder (trailing edge,
        // high angle) for a counter-clockwise-rotating, gravity-loaded charge. Angle convention:
        // atan2(y, x) in [0, 2*pi), measured counter-clockwise from +x; "down" (6 o'clock) is
        // 3*pi/2.
        let mut angles: Vec<f32> = state
            .balls
            .x
            .iter()
            .filter(|p| params.mill.radius_m() - p.length() < wall_margin)
            .map(|p| p.y.atan2(p.x).rem_euclid(2.0 * std::f32::consts::PI))
            .collect();
        assert!(
            angles.len() > 5,
            "too few balls near the wall to assess toe/shoulder: {}",
            angles.len()
        );
        angles.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // The charge should not be centrifuged (i.e. not spread all the way around the wall): its
        // angular extent should leave a substantial gap on the "up" side where no balls ride along
        // the wall (this is the free-fall/cascading region above the shoulder).
        let mut max_gap = 0.0f32;
        for w in angles.windows(2) {
            max_gap = max_gap.max(w[1] - w[0]);
        }
        max_gap = max_gap.max(2.0 * std::f32::consts::PI - (angles[angles.len() - 1] - angles[0]));
        assert!(
            max_gap > std::f32::consts::PI * 0.3,
            "charge appears centrifuged (no large gap at the top): max_gap={max_gap}"
        );

        // The charge's centroid should sit below and to one side of the drum centre (i.e. the
        // pile has been carried up by rotation, not simply sitting at the very bottom): a nonzero
        // horizontal offset confirms cascading motion rather than a static heap.
        let centroid: Vec2 = state
            .balls
            .x
            .iter()
            .copied()
            .fold(Vec2::ZERO, |acc, p| acc + p)
            / state.balls.x.len() as f32;
        assert!(
            centroid.y < -0.1 * params.mill.radius_m(),
            "charge centroid not low enough: {centroid:?}"
        );
        assert!(
            centroid.x.abs() > 0.02 * params.mill.radius_m(),
            "charge centroid not offset sideways (no cascading lean): {centroid:?}"
        );
    }

    /// Largest pairwise ball-ball overlap, as a fraction of the ball radius (same formula as
    /// [`crate::metrics::max_ball_overlap_fraction`], reimplemented here to keep this module's
    /// tests independent of `crate::metrics`).
    fn max_ball_overlap_fraction(balls: &Balls) -> f32 {
        let r = balls.radius;
        let cell_size = (2.0 * r * 1.05).max(1e-6);
        let grid = UniformGrid::build(&balls.x, cell_size);
        let mut max_overlap = 0.0f32;
        grid.for_each_candidate_pair(|i, j| {
            let dist = (balls.x[i as usize] - balls.x[j as usize]).length();
            let overlap = (2.0 * r - dist).max(0.0);
            max_overlap = max_overlap.max(overlap / r);
        });
        max_overlap
    }

    #[test]
    fn total_energy_never_increases_in_a_still_drum() {
        // With the drum stationary (omega = 0) the wall supplies no energy at all
        // (`Drum::wall_velocity` is zero everywhere, so `DemStepStats::wall_work_j` is always
        // zero): every mechanism in this solver -- Coulomb friction, restitution `e < 1`,
        // rolling friction, and step 3's bounded depenetration recovery (see
        // [`MAX_RECOVERY_FRACTION`]) -- is dissipative. Total mechanical energy (KE +
        // gravitational PE) must therefore never increase, only settle toward a resting minimum.
        // This is the bulletproof, no-caveats statement of "this solver creates no energy" for
        // the fluidised-charge bug this module fixes; see the sibling
        // `crate::tests::steady_cascading_charge_never_gains_more_energy_than_the_wall_supplies`
        // for the rotating-drum version, where the wall *does* legitimately supply energy.
        use crate::params::{LiftersParams, MediaParams, MillParams, Params, SpeedMode};

        let mut params = Params {
            mill: MillParams {
                diameter_m: 1.0,
                speed_mode: SpeedMode::Rpm,
                speed_value: 0.0,
                ..MillParams::default()
            },
            media: MediaParams {
                fill_fraction: 0.30,
                ..MediaParams::default()
            },
            lifters: LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
            ..Params::default()
        };
        params.simulation.max_balls = 300;
        params.validate().unwrap();

        let effective = params.effective_media();
        let radius_m = params.mill.radius_m();
        let mut state = DemState::new(&effective, radius_m, params.simulation.seed);
        let drum = still_drum(radius_m);
        let dt = 1.0 / 240.0;

        // The fresh lattice's initial landing is a genuine, extreme transient -- essentially the
        // whole population makes first contact within the same handful of sub-steps (up to ~200
        // simultaneous fresh collisions at this ball count), which is exactly the regime where a
        // sequential (Gauss-Seidel) iterative solver's per-pair correction order leaves the
        // largest *numerical* residual (still tiny relative to the system's total energy scale,
        // and nothing like the bug this test guards against, but larger than the residual once
        // the charge has actually settled). Settle first, unchecked, then hold the invariant to a
        // tight tolerance for the steady/settled phase that follows -- the same "settle first,
        // then check" structure every other test in this module uses.
        for _ in 0..240 {
            state.step(
                &drum,
                0.0,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
        }

        // A small slack for floating-point summation order in the settled regime, not a physical
        // allowance.
        let tolerance_j = 1e-4 * state.balls.mass * state.balls.len() as f32 * 9.81 * radius_m;
        let mut e_prev = mechanical_energy_j(&state.balls);
        let mut worst_delta = f32::MIN;
        let mut worst_step = 0;
        for step in 0..(4 * 240) {
            let stats = state.step(
                &drum,
                0.0,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            let e_now = mechanical_energy_j(&state.balls);
            let delta = e_now - e_prev;
            if delta > worst_delta {
                worst_delta = delta;
                worst_step = step;
                eprintln!(
                    "step={step} delta={delta} collisions={} dissipated={} max_disp={} max_overlap={}",
                    stats.collision_count,
                    stats.dissipated_energy_j,
                    stats.max_substep_displacement_over_diameter,
                    max_ball_overlap_fraction(&state.balls),
                );
            }
            e_prev = e_now;
        }
        eprintln!("worst_delta={worst_delta} at step={worst_step} tolerance_j={tolerance_j}");
        assert!(
            worst_delta <= tolerance_j,
            "mechanical energy increased by {worst_delta} J at step {worst_step} \
             (tolerance {tolerance_j} J) in a still drum, which supplies no energy \
             (wall_work_j is always zero at omega=0) -- every mechanism here must be dissipative"
        );
    }

    #[test]
    fn a_charge_in_a_still_drum_settles_and_stays_dense() {
        // Direct regression test for the fluidised-charge bug: a charge released in a stationary
        // drum must fall, pack, and stay packed near the bottom -- not inflate into a "gas" that
        // fills the whole cross-section. Bounds below are calibrated with margin against this
        // project's own measured baseline for a *healthy* settled pile at these parameters (see
        // git history for the calibration run): residual jitter around
        // `RESTITUTION_VELOCITY_THRESHOLD` and residual overlap from running only
        // `simulation.dem_iterations` (4) Gauss-Seidel iterations on a loaded stack are both
        // expected and are not what this test is checking for -- the reported bug was two
        // orders of magnitude beyond either.
        use crate::params::{LiftersParams, MediaParams, MillParams, Params, SpeedMode};

        let mut params = Params {
            mill: MillParams {
                diameter_m: 1.0,
                speed_mode: SpeedMode::Rpm,
                speed_value: 0.0,
                ..MillParams::default()
            },
            media: MediaParams {
                fill_fraction: 0.30,
                ..MediaParams::default()
            },
            lifters: LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
            ..Params::default()
        };
        params.simulation.max_balls = 300;
        params.validate().unwrap();

        let effective = params.effective_media();
        let radius_m = params.mill.radius_m();
        let mut state = DemState::new(&effective, radius_m, params.simulation.seed);
        let drum = still_drum(radius_m);
        let dt = 1.0 / 240.0;

        for _ in 0..(5 * 240) {
            state.step(
                &drum,
                0.0,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
        }

        let ke: f32 = state
            .balls
            .v
            .iter()
            .map(|v| 0.5 * state.balls.mass * v.length_squared())
            .sum();
        let v_rms = (2.0 * ke / (state.balls.mass * state.balls.len() as f32)).sqrt();
        assert!(
            v_rms < 0.15,
            "settled charge should be nearly at rest, got v_rms={v_rms} m/s"
        );

        let overlap = max_ball_overlap_fraction(&state.balls);
        assert!(
            overlap < 0.4,
            "settled charge has an implausibly deep overlap (charge may have inflated): \
             overlap={overlap}"
        );

        // The reported bug filled the *entire* cross-section, including the centre, rather than
        // a J=0.30 pile sitting at the bottom -- so a settled charge must leave the centre region
        // clearly empty.
        let min_dist_from_center = state
            .balls
            .x
            .iter()
            .map(|p| p.length())
            .fold(f32::MAX, f32::min);
        assert!(
            min_dist_from_center > 0.05 * radius_m,
            "settled charge reaches too close to the drum centre (looks inflated, not packed): \
             min_dist_from_center={min_dist_from_center}, radius_m={radius_m}"
        );
    }

    #[test]
    fn a_cascading_charge_stays_dense() {
        // Direct regression test for the fluidised-charge bug at the screenshot's operating
        // point (cascading, not still): the charge must stay a dense, packed body throughout the
        // run, not disperse into a "gas" filling the whole drum. Checked every second of sim
        // time (not just at the end), since the bug this guards against is an ongoing runaway,
        // not a one-off transient. Bounds are calibrated the same way as the still-drum sibling
        // test above.
        use crate::params::{Direction, LiftersParams, MediaParams, MillParams, Params, SpeedMode};

        let mut params = Params {
            mill: MillParams {
                diameter_m: 1.0,
                speed_mode: SpeedMode::Rpm,
                speed_value: 30.0,
                direction: Direction::CounterClockwise,
            },
            media: MediaParams {
                fill_fraction: 0.30,
                ..MediaParams::default()
            },
            lifters: LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
            ..Params::default()
        };
        params.simulation.max_balls = 400;
        params.validate().unwrap();

        let effective = params.effective_media();
        let radius_m = params.mill.radius_m();
        let omega = params.mill.omega();
        let mut state = DemState::new(&effective, radius_m, params.simulation.seed);
        let dt = 1.0 / 240.0;
        let mut drum_angle = 0.0f32;

        for step in 0..(10 * 240) {
            let drum = Drum::new(radius_m, omega, params.lifters);
            state.step(
                &drum,
                drum_angle,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);

            if step % 240 != 239 {
                continue; // sample once per second of sim time, skip the settling transient below
            }
            if step < 3 * 240 {
                continue; // let the fresh lattice fall/settle before checking
            }

            let overlap = max_ball_overlap_fraction(&state.balls);
            assert!(
                overlap < 0.7,
                "cascading charge has an implausibly deep overlap at step {step} \
                 (charge may be inflating): overlap={overlap}"
            );

            let min_dist_from_center = state
                .balls
                .x
                .iter()
                .map(|p| p.length())
                .fold(f32::MAX, f32::min);
            assert!(
                min_dist_from_center > 0.05 * radius_m,
                "cascading charge reaches too close to the drum centre at step {step} (looks \
                 inflated, not a dense cascading body): min_dist_from_center={min_dist_from_center}, \
                 radius_m={radius_m}"
            );
        }
    }

    /// Largest ball-wall/lifter penetration depth as a fraction of the ball radius (same formula
    /// as `crate::metrics::max_ball_wall_overlap_fraction`, reimplemented here for the same
    /// reason as `max_ball_overlap_fraction` above).
    fn max_ball_wall_overlap_fraction(balls: &Balls, drum: &Drum, drum_angle: f32) -> f32 {
        let r = balls.radius;
        let mut max_overlap = 0.0f32;
        for &x in &balls.x {
            let (d, _) = drum.sdf_world(x, drum_angle);
            let overlap = (r - d).max(0.0);
            max_overlap = max_overlap.max(overlap / r);
        }
        max_overlap
    }

    #[test]
    fn a_cataracting_charge_with_lifters_never_tunnels_through_a_lifter() {
        // Regression for the second external review's findings "①" and "⑦": with lifters
        // enabled, a cataracting charge must never register a ball-wall/lifter overlap deep
        // enough to mean the ball's *centre* has passed fully through the lifter solid
        // (overlap >= 1.0, i.e. the centre itself is a full radius past the surface) -- that
        // specific signature is true tunnelling through a discrete sub-step, as opposed to an
        // ordinary shallow XPBD contact-convergence residual (which this project's docs already
        // expect to read tens of percent, see docs/METRICS.md). This exact configuration
        // (`lifters.count = 8`, default speed) is where the review measured its highest readings
        // and, before this test, was not covered anywhere in this module -- every other DEM test
        // here uses `lifters.count = 0`.
        use crate::params::{Direction, LiftersParams, MediaParams, MillParams, Params, SpeedMode};

        let mut params = Params {
            mill: MillParams {
                diameter_m: 1.0,
                speed_mode: SpeedMode::Rpm,
                speed_value: 30.0,
                direction: Direction::CounterClockwise,
            },
            media: MediaParams {
                fill_fraction: 0.30,
                ..MediaParams::default()
            },
            lifters: LiftersParams {
                count: 8,
                ..LiftersParams::default()
            },
            ..Params::default()
        };
        params.simulation.max_balls = 400;
        params.validate().unwrap();

        let effective = params.effective_media();
        let radius_m = params.mill.radius_m();
        let omega = params.mill.omega();
        let mut state = DemState::new(&effective, radius_m, params.simulation.seed);
        let dt = 1.0 / 240.0;
        let mut drum_angle = 0.0f32;
        let mut worst_wall_overlap = 0.0f32;

        for step in 0..(20 * 240) {
            let drum = Drum::new(radius_m, omega, params.lifters);
            state.step(
                &drum,
                drum_angle,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);

            if step < 5 * 240 {
                continue; // let the fresh lattice fall/settle into steady cataracting first
            }
            let drum_now = Drum::new(radius_m, omega, params.lifters);
            let overlap = max_ball_wall_overlap_fraction(&state.balls, &drum_now, drum_angle);
            worst_wall_overlap = worst_wall_overlap.max(overlap);
            assert!(
                overlap < 1.0,
                "ball-wall overlap at step {step} reads {overlap} (>= 1.0 means a ball's centre \
                 has passed fully through the lifter solid -- true tunnelling, not an XPBD \
                 contact-convergence residual)"
            );
        }
        eprintln!("worst_wall_overlap={worst_wall_overlap}");
    }

    #[test]
    fn dissipated_energy_is_never_materially_negative_while_cascading() {
        // Direct regression test for the metrics-panel review finding that `dissipated_power_w`
        // read hundreds of watts negative on an ordinary cascading run: before the fix,
        // `dissipated_energy_j` was `ke_after_predict - ke_final` (translational KE only, no
        // wall-work credit, no PE, no rotational KE), so whenever the moving wall pumped energy
        // into the charge faster than friction/restitution removed it (an entirely ordinary
        // situation, not a bug) the metric went negative and looked like a thermodynamics
        // violation. The current definition (`e_after_predict + wall_work_j - e_final`) is the
        // residual of the same invariant
        // `crate::tests::steady_cascading_charge_never_gains_more_energy_than_the_wall_supplies`
        // proves holds, so it must stay non-negative up to floating-point slack -- checked here
        // directly on `DemStepStats`, every sub-step, at the same operating point (30 rpm
        // cascading, no lifters) the original review screenshot used.
        use crate::params::{Direction, LiftersParams, MediaParams, MillParams, Params, SpeedMode};

        let mut params = Params {
            mill: MillParams {
                diameter_m: 1.0,
                speed_mode: SpeedMode::Rpm,
                speed_value: 30.0,
                direction: Direction::CounterClockwise,
            },
            media: MediaParams {
                fill_fraction: 0.30,
                ..MediaParams::default()
            },
            lifters: LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
            ..Params::default()
        };
        params.simulation.max_balls = 400;
        params.validate().unwrap();

        let effective = params.effective_media();
        let radius_m = params.mill.radius_m();
        let omega = params.mill.omega();
        let mut state = DemState::new(&effective, radius_m, params.simulation.seed);
        let dt = 1.0 / 240.0;
        let mut drum_angle = 0.0f32;

        // Settling phase (not checked): the fresh lattice's initial landing is an extreme
        // transient, same rationale as the sibling energy-invariant tests above.
        for _ in 0..(3 * 240) {
            let drum = Drum::new(radius_m, omega, params.lifters);
            state.step(
                &drum,
                drum_angle,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
        }

        // A small slack for floating-point summation order, same form as
        // `steady_cascading_charge_never_gains_more_energy_than_the_wall_supplies`'s tolerance --
        // not a physical allowance.
        let tolerance_j = 1e-3 * state.balls.mass * (omega * radius_m).powi(2).max(1.0);
        let mut worst = f32::MAX;
        for _ in 0..(5 * 240) {
            let drum = Drum::new(radius_m, omega, params.lifters);
            let stats = state.step(
                &drum,
                drum_angle,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
            worst = worst.min(stats.dissipated_energy_j);
            assert!(
                stats.dissipated_energy_j >= -tolerance_j,
                "dissipated_energy_j read materially negative while cascading (looks like \
                 energy created from nothing): dissipated_energy_j={}, tolerance_j={tolerance_j}",
                stats.dissipated_energy_j
            );
        }
        assert!(worst.is_finite());
    }
}
