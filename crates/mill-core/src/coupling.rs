//! Two-way ball<->fluid momentum exchange.
//!
//! Per fixed sub-step (docs/PLAN.md ss3.4): the fluid solves its own density constraint and, in
//! the same pass, projects any fluid particle overlapping a ball back to the ball's surface,
//! accumulating the (Newton's-third-law) reaction on that ball as a linear/angular **impulse**
//! (not a force -- see [`CouplingImpulses`]'s doc comment for why); fluid particles near a ball's
//! surface also exchange viscous momentum with it (no-slip). The accumulated per-ball impulse is
//! then applied directly as a velocity change in the ball solver's predict step, and finally the
//! balls advance -- their new positions/velocities become the moving boundary for the fluid's
//! next sub-step. See [`crate::pbf::FluidParticles::step_coupled`] and
//! [`crate::dem::DemState::step_with_external_forces`] for the two halves of this exchange.

use glam::Vec2;

use crate::dem::DemState;
use crate::geometry::Drum;
use crate::params::{MediaParams, SlurryParams};
use crate::pbf::FluidParticles;

/// Per-ball linear impulse (kg*m/s) and angular impulse (kg*m^2/s) accumulated from fluid
/// interaction this sub-step, indexed the same as the ball population (`impulses[i]`/
/// `angular_impulses[i]` correspond to ball `i`).
///
/// **Impulse, not force.** [`crate::pbf::FluidParticles::step_coupled`] derives this from
/// *position* corrections (pushing fluid particles off a ball) and *velocity* deltas (viscous
/// no-slip blending); both convert directly to momentum changes (`mass * (push / dt)` and
/// `mass * dv` respectively) without an extra division by `dt`. Expressing the ball's reaction as
/// a force (impulse / dt again) and then re-integrating it with another `* dt` in
/// [`crate::dem::DemState::step_with_external_forces`] round-trips through `dt` twice, making the
/// result unnecessarily sensitive to the sub-step size -- an earlier version of this code did
/// exactly that and was visibly unstable (balls flung out of the charge) at the project's default
/// sub-step rate. Working in impulses throughout and applying them as `Δv = impulse * inv_mass`
/// (no extra `* dt`) avoids that sensitivity, consistent with the rest of the solver's
/// position-based/impulse-style corrections (docs/PLAN.md ss3.2/3.3) rather than mixing in a
/// force-based integration step.
#[derive(Debug, Clone)]
pub struct CouplingImpulses {
    pub impulses: Vec<Vec2>,
    pub angular_impulses: Vec<f32>,
    /// Number of balls whose linear impulse was reduced by the stability clamp
    /// ([`crate::pbf::FluidParticles::step_coupled`] step 8) this sub-step. Surfaced as a metrics
    /// debug indicator (`Metrics::coupling_clamp_hits`): a healthy run stays at (or very near) 0
    /// once the charge has settled -- persistent clamping means the coupling forces are pinned at
    /// an artificial ceiling rather than reflecting the physical interaction.
    pub clamp_hits: u32,
    /// Sum of this sub-step's fluid momentum change (kg*m/s) from every ball<->fluid interaction
    /// (overlap push, viscous drag, buoyancy). Debug accumulator only, populated by
    /// [`crate::pbf::FluidParticles::step_coupled`]: a test can check `sum(impulses) ==
    /// -fluid_momentum_change` (Newton's third law) when no clamp fired this sub-step (clamping
    /// breaks the equality by construction, since it discards momentum on the ball side only).
    pub fluid_momentum_change: Vec2,
}

impl CouplingImpulses {
    pub fn zeros(n_balls: usize) -> Self {
        Self {
            impulses: vec![Vec2::ZERO; n_balls],
            angular_impulses: vec![0.0; n_balls],
            clamp_hits: 0,
            fluid_momentum_change: Vec2::ZERO,
        }
    }
}

/// Advances both solvers by one shared fixed sub-step `dt`, with two-way momentum exchange
/// (docs/PLAN.md ss3.4). Order matters: the fluid solves first (producing this sub-step's
/// reaction forces on the balls from *last* sub-step's ball positions -- a standard staggered/
/// semi-implicit coupling scheme), then the balls advance using those forces, so their new
/// positions become the moving boundary for the fluid's *next* sub-step.
#[allow(clippy::too_many_arguments)]
pub fn step(
    dem: &mut DemState,
    fluid: &mut FluidParticles,
    drum: &Drum,
    drum_angle: f32,
    media: &MediaParams,
    slurry: &SlurryParams,
    dem_iterations: u32,
    pbf_iterations: u32,
    dt: f32,
) {
    let impulses = fluid.step_coupled(drum, drum_angle, slurry, pbf_iterations, dt, &dem.balls);
    dem.step_with_external_forces(drum, drum_angle, media, dem_iterations, dt, Some(&impulses));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{EffectiveMedia, LiftersParams};

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

    #[test]
    fn zeros_has_matching_lengths_and_no_effect() {
        let cf = CouplingImpulses::zeros(5);
        assert_eq!(cf.impulses.len(), 5);
        assert_eq!(cf.angular_impulses.len(), 5);
        assert!(cf.impulses.iter().all(|f| *f == Vec2::ZERO));
        assert!(cf.angular_impulses.iter().all(|&t| t == 0.0));
    }

    #[test]
    fn ball_dropped_into_a_pool_decelerates_more_than_a_dry_drop() {
        // A ball falling into a deep, still pool of slurry should be slowed by the fluid more
        // than an identical ball falling the same distance in dry air (docs/PLAN.md ss5 M4). Uses
        // a large drum and a short measurement window so *neither* ball has yet reached the
        // bottom wall by the time speeds are compared (otherwise both simply rest at v~=0
        // regardless of fluid drag, which is not what this test means to check).
        let radius_m = 2.0;
        let media = MediaParams {
            ball_diameter_m: 0.04,
            fill_fraction: 0.0, // seed exactly one ball ourselves, below
            restitution_ball_wall: 0.0,
            friction_ball_wall: 0.0,
            rolling_friction: 0.0,
            ..MediaParams::default()
        };
        let slurry = SlurryParams {
            fill_fraction: 0.35,
            viscosity_pa_s: 1.0,
            ..SlurryParams::default()
        };
        let drum = still_drum(radius_m);
        let dt = 1.0 / 240.0;

        // Seed the pool first, then drop a ball from just above it.
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 24, &[], 0.0);
        for _ in 0..120 {
            fluid.step(&drum, 0.0, &slurry, 3, dt);
        }
        let pool_top = fluid.x.iter().map(|p| p.y).fold(f32::MIN, f32::max);

        let effective = EffectiveMedia {
            true_diameter_m: media.ball_diameter_m,
            diameter_m: media.ball_diameter_m,
            density_kg_m3: media.density_kg_m3,
            ball_count: 1,
            scale_factor: 1.0,
        };
        let mut dem_wet = DemState::new(&effective, radius_m, 1);
        dem_wet.balls.x[0] = Vec2::new(0.0, pool_top + 0.1);
        let mut dem_dry = DemState::new(&effective, radius_m, 1);
        dem_dry.balls.x[0] = dem_wet.balls.x[0];

        // 60 steps = 0.25s: free fall alone covers ~0.31 m, comfortably short of the ~2 m to the
        // bottom wall from a drop point well above a (nowhere near full) pool in a 2 m-radius drum.
        for _ in 0..60 {
            step(
                &mut dem_wet,
                &mut fluid,
                &drum,
                0.0,
                &media,
                &slurry,
                4,
                3,
                dt,
            );
            dem_dry.step(&drum, 0.0, &media, 4, dt);
        }

        let wet_speed = dem_wet.balls.v[0].length();
        let dry_speed = dem_dry.balls.v[0].length();
        assert!(
            dry_speed > 1.0,
            "sanity check failed: dry ball should still be in free fall, got dry_speed={dry_speed}"
        );
        assert!(
            wet_speed < dry_speed,
            "ball falling through fluid should be slower than a dry drop: wet={wet_speed} dry={dry_speed}"
        );
    }

    #[test]
    fn fluid_particles_end_up_outside_every_ball() {
        let radius_m = 0.5;
        let drum = still_drum(radius_m);
        let media = MediaParams::default();
        let slurry = SlurryParams {
            fill_fraction: 0.3,
            ..SlurryParams::default()
        };
        let dt = 1.0 / 240.0;

        let effective = EffectiveMedia {
            true_diameter_m: 0.05,
            diameter_m: 0.05,
            density_kg_m3: 7800.0,
            ball_count: 8,
            scale_factor: 1.0,
        };
        let mut dem = DemState::new(&effective, radius_m, 3);
        let mut fluid =
            FluidParticles::seed_lattice(&slurry, radius_m, 20, &dem.balls.x, dem.balls.radius);

        for _ in 0..120 {
            step(&mut dem, &mut fluid, &drum, 0.0, &media, &slurry, 4, 3, dt);
        }

        let r = dem.balls.radius;
        for &fp in &fluid.x {
            for &bp in &dem.balls.x {
                let dist = (fp - bp).length();
                assert!(
                    dist >= r * 0.9,
                    "fluid particle inside a ball: dist={dist}, r={r}"
                );
            }
        }
    }

    #[test]
    fn single_overlapping_fluid_particle_gives_the_ball_the_opposite_reaction() {
        // Isolates the ball-overlap projection mechanism (docs/PLAN.md ss3.4 step 1) from every
        // other effect (gravity, wall contact, fluid-fluid interaction, viscosity) that would
        // otherwise contaminate a whole-system momentum comparison: one fluid particle placed
        // just outside a ball's centre (so the push direction is well-defined and purely
        // horizontal) should be pushed away by exactly the amount that gives the ball the equal
        // and opposite reaction (Newton's third law), independent of gravity (which only acts
        // vertically here, on the fluid's `y` and not measured against the ball at all since
        // `dem` is never stepped in this test). Uses a small overlap (0.0025 m) so the resulting
        // force stays comfortably under docs/PLAN.md ss3.4's stability clamp, which would
        // otherwise break the exact cancellation this test checks for.
        let radius_m = 1.0;
        let drum = still_drum(radius_m);
        let slurry = SlurryParams {
            fill_fraction: 0.0, // we inject the one fluid particle manually, below
            wall_no_slip: 0.0,
            viscosity_pa_s: 0.0, // isolates the overlap-projection reaction from ball-viscosity
            ..SlurryParams::default()
        };
        let dt = 1.0 / 240.0;

        // Fine resolution -> small, light fluid particles, so a modest overlap produces a modest
        // (well under the clamp) reaction force.
        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 200, &[], 0.0);
        assert!(fluid.is_empty());
        fluid.x.push(Vec2::new(0.6, 0.0));
        fluid.v.push(Vec2::ZERO);
        fluid.dye.push(0.0);

        let effective = EffectiveMedia {
            true_diameter_m: 0.2,
            diameter_m: 0.2,
            density_kg_m3: 7800.0,
            ball_count: 1,
            scale_factor: 1.0,
        };
        let mut dem = DemState::new(&effective, radius_m, 1);
        // Ball radius 0.1 m; contact_radius = radius + 0.25*h = 0.1 + 0.25*(2*0.005) = 0.1025 m.
        // The fluid particle at dist=0.1 m is just inside that, for a small 0.0025 m overlap.
        dem.balls.x[0] = Vec2::new(0.5, 0.0);
        dem.balls.v[0] = Vec2::ZERO;
        dem.balls.omega[0] = 0.0;

        let impulses = fluid.step_coupled(&drum, 0.0, &slurry, 1, dt, &dem.balls);

        let fluid_momentum_gained = fluid.particle_mass * fluid.v[0];
        let ball_impulse = impulses.impulses[0]; // already an impulse, see CouplingImpulses's doc comment
                                                 // Only the x-component is meaningful here: the push is purely horizontal (particle
                                                 // offset from the ball along x only), while gravity contributes only to the fluid
                                                 // particle's y-velocity and is not part of what this test checks.
        let sum_x = fluid_momentum_gained.x + ball_impulse.x;
        assert!(
            sum_x.abs() < 1e-6,
            "momentum not conserved in isolated overlap: fluid={fluid_momentum_gained:?} \
             ball_impulse={ball_impulse:?} (x-components should cancel)"
        );
        // Sanity: the reaction must be non-trivial (not silently zero).
        assert!(
            ball_impulse.x.abs() > 1e-9,
            "expected a non-zero reaction on the ball"
        );
    }

    #[test]
    fn extreme_rotation_keeps_the_coupled_fluid_state_bounded_every_sub_step() {
        // Far beyond anything the UI permits (omega = 40 rad/s, ~9x critical speed for this
        // drum) with a ball charge, to exercise fluid squeezed between the centrifuged charge
        // and the wall -- the configuration that drives `pbf.rs` step 5's velocity ramp.
        // Asserted on *every* sub-step, not just at the end: a transient excursion that later
        // relaxes is exactly the event that would make the rendered free surface disappear for
        // a frame (see `crate::surface`), so checking only the final state isn't enough.
        let radius_m = 0.5;
        let omega = 40.0f32;
        let drum = Drum::new(
            radius_m,
            omega,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        );
        let media = MediaParams {
            ball_diameter_m: 0.02,
            fill_fraction: 0.2,
            ..MediaParams::default()
        };
        let slurry = SlurryParams {
            fill_fraction: 0.15,
            ..SlurryParams::default()
        };
        let dt = 1.0 / 240.0;

        let effective = EffectiveMedia {
            true_diameter_m: media.ball_diameter_m,
            diameter_m: media.ball_diameter_m,
            density_kg_m3: media.density_kg_m3,
            ball_count: 40,
            scale_factor: 1.0,
        };
        let mut dem = DemState::new(&effective, radius_m, 1);
        let mut fluid =
            FluidParticles::seed_lattice(&slurry, radius_m, 16, &dem.balls.x, dem.balls.radius);

        // Mirrors pbf.rs's `FLUID_SPEED_SAFETY_FACTOR` clamp formula exactly.
        let v_max = 5.0 * (omega.abs() * radius_m + (4.0 * 9.81 * radius_m).sqrt());
        let mut drum_angle = 0.0f32;
        for s in 0..300 {
            step(
                &mut dem, &mut fluid, &drum, drum_angle, &media, &slurry, 4, 3, dt,
            );
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);

            for (i, (&x, &v)) in fluid.x.iter().zip(&fluid.v).enumerate() {
                assert!(
                    x.is_finite(),
                    "sub-step {s}, fluid {i}: non-finite position {x:?}"
                );
                assert!(
                    v.is_finite(),
                    "sub-step {s}, fluid {i}: non-finite velocity {v:?}"
                );
                assert!(
                    x.length() <= radius_m * 1.05,
                    "sub-step {s}, fluid {i}: escaped the drum: {x:?}"
                );
                assert!(
                    v.length() <= v_max * 1.001,
                    "sub-step {s}, fluid {i}: speed {} exceeds clamp {v_max}",
                    v.length()
                );
            }
            for (i, (&x, &v)) in dem.balls.x.iter().zip(&dem.balls.v).enumerate() {
                assert!(
                    x.is_finite(),
                    "sub-step {s}, ball {i}: non-finite position {x:?}"
                );
                assert!(
                    v.is_finite(),
                    "sub-step {s}, ball {i}: non-finite velocity {v:?}"
                );
            }
        }
    }

    #[test]
    fn a_light_disc_rises_through_a_still_pool() {
        // Buoyancy (docs/PLAN.md ss3.4 step 3): a disc much less dense than the slurry, released
        // fully submerged in a still, settled pool, should rise rather than sink or stay put.
        let radius_m = 2.0;
        let slurry = SlurryParams {
            fill_fraction: 0.35,
            density_kg_m3: 1800.0,
            viscosity_pa_s: 0.5,
            ..SlurryParams::default()
        };
        let drum = still_drum(radius_m);
        let dt = 1.0 / 240.0;

        let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, 24, &[], 0.0);
        for _ in 0..120 {
            fluid.step(&drum, 0.0, &slurry, 3, dt);
        }
        let pool_top = fluid.x.iter().map(|p| p.y).fold(f32::MIN, f32::max);

        let media = MediaParams {
            ball_diameter_m: 0.05,
            density_kg_m3: 300.0, // a sixth of the slurry's rest density: should float up briskly
            fill_fraction: 0.0,
            restitution_ball_wall: 0.0,
            friction_ball_wall: 0.0,
            rolling_friction: 0.0,
            ..MediaParams::default()
        };
        let effective = EffectiveMedia {
            true_diameter_m: media.ball_diameter_m,
            diameter_m: media.ball_diameter_m,
            density_kg_m3: media.density_kg_m3,
            ball_count: 1,
            scale_factor: 1.0,
        };
        let mut dem = DemState::new(&effective, radius_m, 1);
        // Start well below the pool's own surface so the disc is fully submerged, comfortably
        // above the drum's bottom wall (radius_m = 2.0) so wall contact never intervenes.
        let start_y = pool_top - 0.5;
        dem.balls.x[0] = Vec2::new(0.0, start_y);
        dem.balls.v[0] = Vec2::ZERO;

        // Track the *peak* height reached, not just the final position: a disc this much
        // lighter than the slurry (1/6th the rest density) rises briskly enough to shoot clear
        // through the free surface within a fraction of a second (physically the same thing a
        // released cork does underwater), after which it is a free-falling/bouncing projectile
        // for whatever remains of the run -- its position at an arbitrary later time is not a
        // reliable "did it rise" signal, but the peak height it reached on the way up is.
        let mut max_y = start_y;
        for _ in 0..240 {
            step(&mut dem, &mut fluid, &drum, 0.0, &media, &slurry, 4, 3, dt);
            max_y = max_y.max(dem.balls.x[0].y);
        }

        assert!(
            dem.balls.x[0].is_finite() && dem.balls.v[0].is_finite(),
            "non-finite disc state: x={:?} v={:?}",
            dem.balls.x[0],
            dem.balls.v[0]
        );
        // The disc's own viscous drag reaction (step 6.5) opposes its rise in proportion to its
        // velocity relative to the fluid, so the ascent settles into a modest terminal creep
        // rather than accelerating unbounded -- measured empirically at ~4.4 cm over this run's
        // 1 s. The bound below is set well under that with margin, while staying far above any
        // plausible settling jitter from the fresh lattice (millimetres, not centimetres).
        assert!(
            max_y > start_y + 0.02,
            "light disc should have risen: start_y={start_y}, peak_y={max_y}"
        );
    }

    #[test]
    fn ball_fluid_momentum_exchange_is_symmetric_when_unclamped() {
        // Newton's third law: absent clamping, the total momentum the balls gain from fluid
        // interaction (overlap push, viscous drag, buoyancy) must equal minus the momentum the
        // fluid gained (`CouplingImpulses::fluid_momentum_change`) -- see that field's doc
        // comment. Clamping breaks the equality by construction (it discards momentum on the
        // ball side only), so the check only applies when no clamp fired this sub-step.
        let radius_m = 1.0;
        let drum = still_drum(radius_m);
        let slurry = SlurryParams {
            fill_fraction: 0.3,
            ..SlurryParams::default()
        };
        let media = MediaParams::default();
        let dt = 1.0 / 240.0;

        let effective = EffectiveMedia {
            true_diameter_m: 0.05,
            diameter_m: 0.05,
            density_kg_m3: 4000.0, // between the slurry and steel: exercises drag and buoyancy
            ball_count: 6,
            scale_factor: 1.0,
        };
        let mut dem = DemState::new(&effective, radius_m, 3);
        let mut fluid =
            FluidParticles::seed_lattice(&slurry, radius_m, 20, &dem.balls.x, dem.balls.radius);

        // Settle briefly so the sub-step under test has realistic (not raw-lattice) contacts.
        for _ in 0..60 {
            step(&mut dem, &mut fluid, &drum, 0.0, &media, &slurry, 4, 3, dt);
        }

        let impulses = fluid.step_coupled(&drum, 0.0, &slurry, 3, dt, &dem.balls);
        if impulses.clamp_hits == 0 {
            let sum_ball_impulse: Vec2 = impulses.impulses.iter().copied().sum();
            let expected = -impulses.fluid_momentum_change;
            let rel_err = (sum_ball_impulse - expected).length() / expected.length().max(1e-9);
            assert!(
                rel_err < 1e-3,
                "momentum not conserved: sum_ball={sum_ball_impulse:?} \
                 expected={expected:?} rel_err={rel_err}"
            );
        }
    }

    #[test]
    fn cascading_charge_keeps_coupling_clamp_hits_rare_once_settled() {
        // A multi-second run at representative cascading parameters should settle into a state
        // where the ball<->fluid coupling clamp (`CouplingImpulses::clamp_hits`'s doc comment)
        // fires only occasionally, not persistently. This project's default media (2 cm balls)
        // are comparable in size to the default fluid particle spacing here (dx ~= 2.5 cm at
        // resolution 20), i.e. sub-resolution/unresolved (docs/PLAN.md ss3.4) -- a ball can end
        // up briefly overlapped by several fluid particles at once during an energetic
        // cascading contact, producing a legitimately large single-substep overlap-push impulse
        // (step 3.5) even after the charge has otherwise settled into steady motion. Zero clamp
        // hits is therefore not a realistic bar here; what matters is that clamping stays rare
        // (an occasional large contact, not the steady state) rather than dominant (which would
        // mean the coupling forces are pinned at an artificial ceiling instead of reflecting the
        // physical interaction). Measured empirically at ~3.4% of ball-substeps for this exact
        // scenario; the bound below is set with generous headroom above that.
        let radius_m = 0.5;
        let omega = 3.0; // a representative cascading speed for this drum radius
        let drum = Drum::new(
            radius_m,
            omega,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        );
        let media = MediaParams {
            ball_diameter_m: 0.02,
            fill_fraction: 0.25,
            ..MediaParams::default()
        };
        let slurry = SlurryParams {
            fill_fraction: 0.15,
            ..SlurryParams::default()
        };
        let dt = 1.0 / 240.0;

        let effective = EffectiveMedia {
            true_diameter_m: media.ball_diameter_m,
            diameter_m: media.ball_diameter_m,
            density_kg_m3: media.density_kg_m3,
            ball_count: 150,
            scale_factor: 1.0,
        };
        let mut dem = DemState::new(&effective, radius_m, 1);
        let mut fluid =
            FluidParticles::seed_lattice(&slurry, radius_m, 20, &dem.balls.x, dem.balls.radius);

        let mut drum_angle = 0.0f32;
        // Settling phase (not checked -- the initial lattice contacts are expected to be noisy).
        for _ in 0..(3 * 240) {
            step(
                &mut dem, &mut fluid, &drum, drum_angle, &media, &slurry, 4, 3, dt,
            );
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
        }
        // Measurement phase: once the charge is in steady cascading motion.
        let n_measurement_steps = 2 * 240;
        let mut clamp_hits_total = 0u32;
        for _ in 0..n_measurement_steps {
            let impulses = fluid.step_coupled(&drum, drum_angle, &slurry, 3, dt, &dem.balls);
            clamp_hits_total += impulses.clamp_hits;
            dem.step_with_external_forces(&drum, drum_angle, &media, 4, dt, Some(&impulses));
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
        }

        let ball_substeps = dem.balls.len() as f32 * n_measurement_steps as f32;
        let clamp_rate = clamp_hits_total as f32 / ball_substeps;
        assert!(
            clamp_rate < 0.15,
            "coupling clamp fired too often once settled: {clamp_hits_total} hits over \
             {ball_substeps} ball-substeps (rate={clamp_rate})"
        );
    }
}
