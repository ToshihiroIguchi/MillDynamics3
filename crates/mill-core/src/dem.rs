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
/// Below this approach speed, a contact is treated as already-resting rather than a fresh impact,
/// so restitution is not (re-)applied every sub-step (which would otherwise make resting contacts
/// buzz indefinitely instead of settling). See docs/PLAN.md ss3.2 step 6.
const RESTITUTION_VELOCITY_THRESHOLD: f32 = 0.02;

fn ball_mass(diameter_m: f32, density_kg_m3: f32) -> f32 {
    let r = diameter_m * 0.5;
    density_kg_m3 * (4.0 / 3.0) * PI * r * r * r
}

fn ball_inertia(mass: f32, radius_m: f32) -> f32 {
    0.4 * mass * radius_m * radius_m // 2/5 m r^2, solid sphere
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
#[derive(Default)]
struct ContactBook {
    /// Ball-ball contacts: key (i, j) with i < j.
    ball_ball_lambda_n: HashMap<(u32, u32), f32>,
    /// Ball-wall contacts: key is the ball index.
    ball_wall_lambda_n: HashMap<u32, f32>,
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
    ) {
        self.step_with_external_forces(drum, drum_angle, media, iterations, dt, None);
    }

    /// As [`DemState::step`], but additionally applies `external` (fluid coupling forces/torques
    /// per ball, docs/PLAN.md ss3.4, [`crate::coupling`]) as an extra acceleration in the predict
    /// step. `external.forces`/`external.torques` must have one entry per ball, same order as
    /// `self.balls.x`, when provided.
    pub fn step_with_external_forces(
        &mut self,
        drum: &Drum,
        drum_angle: f32,
        media: &MediaParams,
        iterations: u32,
        dt: f32,
        external: Option<&crate::coupling::CouplingForces>,
    ) {
        let balls = &mut self.balls;
        if balls.is_empty() || dt <= 0.0 {
            return;
        }
        let n = balls.len();
        let w = balls.inv_mass();
        let w_rot = balls.inv_inertia();
        let r = balls.radius;

        // --- 1. Predict ---------------------------------------------------------------------
        let x0: Vec<Vec2> = balls.x.clone();
        let theta0: Vec<f32> = balls.theta.clone();
        let v_pre: Vec<Vec2> = balls.v.clone();
        for i in 0..n {
            balls.v[i].y += GRAVITY * dt;
            if let Some(ext) = external {
                balls.v[i] += ext.forces[i] * w * dt;
                balls.omega[i] += ext.torques[i] * w_rot * dt;
            }
            balls.x[i] += balls.v[i] * dt;
            balls.theta[i] += balls.omega[i] * dt;
        }

        // --- 2. Broad-phase (ball-ball candidates only; wall is checked directly per ball) ----
        let cell_size = (2.0 * r * 1.05).max(1e-6);
        let grid = UniformGrid::build(&balls.x, cell_size);
        let mut ball_ball_pairs: Vec<(u32, u32)> = Vec::new();
        grid.for_each_candidate_pair(|i, j| ball_ball_pairs.push((i, j)));

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
                let d_lambda = -c / w_sum;
                balls.x[i] += w * d_lambda * n_hat;
                balls.x[j] -= w * d_lambda * n_hat;
                *book
                    .ball_ball_lambda_n
                    .entry((i as u32, j as u32))
                    .or_insert(0.0) += d_lambda;
            }
            for i in 0..n {
                let (d, n_hat) = drum.sdf_world(balls.x[i], drum_angle);
                let c = d - r;
                if c >= 0.0 {
                    continue;
                }
                if w <= 0.0 {
                    continue;
                }
                let d_lambda = -c / w;
                balls.x[i] += w * d_lambda * n_hat;
                *book.ball_wall_lambda_n.entry(i as u32).or_insert(0.0) += d_lambda;
            }
        }

        // --- 4. Reconstruct velocities from the position change -------------------------------
        for i in 0..n {
            balls.v[i] = (balls.x[i] - x0[i]) / dt;
            balls.omega[i] = angle_diff(balls.theta[i], theta0[i]) / dt;
        }

        // --- 5. Friction (single pass, Coulomb-clamped position correction) -------------------
        for (&(i, j), &lambda_n) in &book.ball_ball_lambda_n {
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
        }
        for (&i, &lambda_n) in &book.ball_wall_lambda_n {
            if lambda_n <= 0.0 {
                continue;
            }
            let iu = i as usize;
            let (_, n_hat) = drum.sdf_world(balls.x[iu], drum_angle);
            let t_hat = Vec2::new(-n_hat.y, n_hat.x);
            let v_wall = drum.wall_velocity(balls.x[iu]);

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
        }
        // Re-reconstruct velocities once more now that friction has adjusted positions/orientation.
        for i in 0..n {
            balls.v[i] = (balls.x[i] - x0[i]) / dt;
            balls.omega[i] = angle_diff(balls.theta[i], theta0[i]) / dt;
        }

        // --- 6. Restitution (single pass, using the pre-solve approach velocity) -------------
        for (&(i, j), &lambda_n) in &book.ball_ball_lambda_n {
            if lambda_n <= 0.0 {
                continue;
            }
            let (iu, ju) = (i as usize, j as usize);
            let delta = balls.x[iu] - balls.x[ju];
            let dist = delta.length();
            let n_hat = if dist > 1e-9 { delta / dist } else { Vec2::X };

            let v_n_pre = (v_pre[iu] - v_pre[ju]).dot(n_hat);
            if v_n_pre >= -RESTITUTION_VELOCITY_THRESHOLD {
                continue;
            }
            let v_n_now = (balls.v[iu] - balls.v[ju]).dot(n_hat);
            let target = -media.restitution_ball_ball * v_n_pre;
            let delta_v_n = target - v_n_now;
            let w_sum = w + w;
            if w_sum <= 0.0 || delta_v_n <= 0.0 {
                continue;
            }
            let impulse = delta_v_n / w_sum;
            balls.v[iu] += w * impulse * n_hat;
            balls.v[ju] -= w * impulse * n_hat;
        }
        for (&i, &lambda_n) in &book.ball_wall_lambda_n {
            if lambda_n <= 0.0 {
                continue;
            }
            let iu = i as usize;
            let (_, n_hat) = drum.sdf_world(balls.x[iu], drum_angle);
            let v_wall = drum.wall_velocity(balls.x[iu]);

            let v_n_pre = (v_pre[iu] - v_wall).dot(n_hat);
            if v_n_pre >= -RESTITUTION_VELOCITY_THRESHOLD {
                continue;
            }
            let v_n_now = (balls.v[iu] - v_wall).dot(n_hat);
            let target = -media.restitution_ball_wall * v_n_pre;
            let delta_v_n = target - v_n_now;
            if w <= 0.0 || delta_v_n <= 0.0 {
                continue;
            }
            balls.v[iu] += delta_v_n * n_hat;
        }

        // --- 7. Rolling resistance -------------------------------------------------------------
        if media.rolling_friction > 0.0 && balls.inertia > 0.0 {
            let mut total_lambda_n = vec![0.0f32; n];
            for (&(i, j), &lambda_n) in &book.ball_ball_lambda_n {
                total_lambda_n[i as usize] += lambda_n.max(0.0);
                total_lambda_n[j as usize] += lambda_n.max(0.0);
            }
            for (&i, &lambda_n) in &book.ball_wall_lambda_n {
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
}
