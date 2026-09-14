//! mill-core: 2D ball-mill DEM + PBF slurry simulation core.
//!
//! Pure Rust, no WASM dependencies -- kept independently testable/benchmarkable via `cargo test`
//! and `cargo bench`. The `wasm-bindgen` wrapper lives in the sibling `mill-wasm` crate.
//!
//! Module responsibilities (see docs/PLAN.md for the full design):
//! - [`params`]: parameter groups, defaults, validation, coarse-graining derivation.
//! - [`geometry`]: drum wall SDF (circle + optional lifters), normals, wall velocity.
//! - [`grid`]: uniform-grid spatial hash shared by the ball solver and PBF neighbour search (M3).
//! - [`dem`]: position-based (XPBD-style) rigid discs for grinding media.
//! - [`pbf`]: Position Based Fluids solver for the slurry (M3).
//! - [`coupling`]: two-way ball<->fluid momentum exchange.
//! - [`surface`]: free-surface extraction (marching squares) for rendering (M3).
//! - [`metrics`]: toe/shoulder angles, slurry level, mixing index, power draw (M3+).
//! - [`rng`]: small deterministic PRNG for reproducible particle-lattice jitter.

pub mod coupling;
pub mod dem;
pub mod geometry;
pub mod grid;
pub mod metrics;
pub mod params;
pub mod pbf;
pub mod rng;
pub mod surface;

use dem::DemState;
use geometry::Drum;
pub use params::{EffectiveMedia, Params};
use pbf::FluidParticles;

/// Fixed simulation sub-step shared by the ball solver and (once it lands, M3) the PBF solver:
/// 1/240 s. `simulation.substeps` sub-steps are taken per nominal 60 Hz rendered frame at 1x time
/// scale (docs/PLAN.md ss3.2/3.3).
pub const FIXED_DT: f32 = 1.0 / 240.0;

/// The simulation instance. This is the single type the `mill-wasm` wrapper (and, indirectly, the
/// web worker's fixed-step accumulator, docs/PLAN.md ss4.1) drives via [`Simulation::step`].
pub struct Simulation {
    params: Params,
    /// Drum rotation angle (radians), kept wrapped into `[0, 2*pi)`.
    drum_angle: f32,
    /// Total simulated time (s).
    sim_time: f64,
    dem: DemState,
    fluid: FluidParticles,
}

impl Simulation {
    /// Creates a new simulation, validating `params` first, seeding the ball population from its
    /// (possibly coarse-grained) effective media (see [`Params::effective_media`]), and -- if
    /// `params.slurry.enabled` -- seeding the fluid (sites overlapping a ball are skipped, docs/
    /// PLAN.md ss3.3). Balls and fluid interact two-way every sub-step, see [`coupling`].
    pub fn new(params: Params) -> Result<Self, String> {
        params.validate()?;
        let radius_m = params.mill.radius_m();
        let effective = params.effective_media();
        let dem = DemState::new(&effective, radius_m, params.simulation.seed);
        let fluid = if params.slurry.enabled {
            FluidParticles::seed_lattice(
                &params.slurry,
                radius_m,
                params.simulation.resolution,
                &dem.balls.x,
                dem.balls.radius,
            )
        } else {
            FluidParticles::seed_lattice(
                &params::SlurryParams {
                    fill_fraction: 0.0,
                    ..params.slurry
                },
                radius_m,
                params.simulation.resolution,
                &[],
                0.0,
            )
        };
        Ok(Self {
            params,
            drum_angle: 0.0,
            sim_time: 0.0,
            dem,
            fluid,
        })
    }

    pub fn params(&self) -> &Params {
        &self.params
    }

    /// Replaces the parameters in place, without resetting `drum_angle`/`sim_time`/the ball
    /// population. Used for "hot" parameter updates that should not restart the simulation
    /// (docs/PLAN.md ss4.1); the worker is responsible for deciding whether a given change instead
    /// needs a full reset (i.e. constructing a new [`Simulation`]) instead of calling this. In
    /// particular, changes to `mill.diameter_m`/`media`/`simulation.seed`/`simulation.max_balls`
    /// take effect in derived quantities (e.g. `omega`) immediately but do **not** reseed or
    /// resize the already-created ball population -- callers must reset for those.
    pub fn set_params(&mut self, params: Params) -> Result<(), String> {
        params.validate()?;
        self.params = params;
        Ok(())
    }

    pub fn drum_angle(&self) -> f32 {
        self.drum_angle
    }

    pub fn sim_time(&self) -> f64 {
        self.sim_time
    }

    /// The current ball population (positions, velocities, orientation, uniform radius/mass).
    pub fn balls(&self) -> &dem::Balls {
        &self.dem.balls
    }

    /// The current fluid (slurry) particle population. Empty when `params.slurry.enabled` was
    /// `false` at construction time.
    pub fn fluid(&self) -> &FluidParticles {
        &self.fluid
    }

    /// Extracts the fluid's free-surface contour(s) at the current state (docs/PLAN.md ss3.5).
    /// Recomputed on demand each call -- callers should not call this more often than needed for
    /// rendering.
    pub fn fluid_surface(&self) -> Vec<surface::Polygon> {
        surface::extract_surface(&self.fluid, self.params.mill.radius_m())
    }

    /// Advances the simulation by `dt` seconds of simulation time, split into
    /// `simulation.substeps` fixed sub-steps (docs/PLAN.md ss3.2/4.1).
    pub fn step(&mut self, dt: f32) {
        let omega = self.params.mill.omega();
        let radius_m = self.params.mill.radius_m();
        let lifters = self.params.lifters;
        let n_substeps = self.params.simulation.substeps.max(1);
        let dem_iterations = self.params.simulation.dem_iterations;
        let sub_dt = dt / n_substeps as f32;

        let pbf_iterations = self.params.simulation.pbf_iterations;

        for _ in 0..n_substeps {
            let drum = Drum::new(radius_m, omega, lifters);
            coupling::step(
                &mut self.dem,
                &mut self.fluid,
                &drum,
                self.drum_angle,
                &self.params.media,
                &self.params.slurry,
                dem_iterations,
                pbf_iterations,
                sub_dt,
            );
            self.drum_angle = (self.drum_angle + omega * sub_dt).rem_euclid(std::f32::consts::TAU);
            self.sim_time += sub_dt as f64;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Direction, SpeedMode};

    #[test]
    fn new_rejects_invalid_params() {
        let mut params = Params::default();
        params.mill.diameter_m = -1.0;
        assert!(Simulation::new(params).is_err());
    }

    #[test]
    fn step_advances_drum_angle_at_correct_rate() {
        let mut params = Params::default();
        params.mill.speed_mode = SpeedMode::Rpm;
        params.mill.speed_value = 60.0; // 1 rev/s => omega = 2*pi rad/s
        params.mill.direction = Direction::CounterClockwise;
        params.media.fill_fraction = 0.0; // isolate drum kinematics from ball dynamics
        let mut sim = Simulation::new(params).unwrap();
        sim.step(0.25); // quarter turn
        assert!((sim.drum_angle() - std::f32::consts::FRAC_PI_2).abs() < 1e-3);
    }

    #[test]
    fn step_wraps_angle_into_0_tau() {
        let mut params = Params::default();
        params.mill.speed_mode = SpeedMode::Rpm;
        params.mill.speed_value = 60.0;
        params.media.fill_fraction = 0.0;
        let mut sim = Simulation::new(params).unwrap();
        sim.step(1.5); // 1.5 revolutions
        assert!(sim.drum_angle() >= 0.0 && sim.drum_angle() < std::f32::consts::TAU);
    }

    #[test]
    fn clockwise_direction_decreases_angle_then_wraps() {
        let mut params = Params::default();
        params.mill.speed_mode = SpeedMode::Rpm;
        params.mill.speed_value = 60.0;
        params.mill.direction = Direction::Clockwise;
        params.media.fill_fraction = 0.0;
        let mut sim = Simulation::new(params).unwrap();
        sim.step(0.25);
        // A quarter turn clockwise from 0 wraps to three-quarters of a turn in [0, tau).
        let expected = std::f32::consts::TAU - std::f32::consts::FRAC_PI_2;
        assert!((sim.drum_angle() - expected).abs() < 1e-3);
    }

    #[test]
    fn sim_time_accumulates_regardless_of_omega() {
        let mut params = Params::default();
        params.mill.speed_value = 0.0;
        params.media.fill_fraction = 0.0;
        let mut sim = Simulation::new(params).unwrap();
        sim.step(0.4);
        sim.step(0.6);
        // `dt` is f32, so the accumulated f64 sim_time inherits f32-level precision per step.
        assert!((sim.sim_time() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn set_params_validates_before_applying() {
        let mut sim = Simulation::new(Params::default()).unwrap();
        let mut bad = Params::default();
        bad.media.fill_fraction = 5.0;
        assert!(sim.set_params(bad).is_err());
        // The original (valid) params must remain in effect after a rejected update.
        assert!(
            (sim.params().media.fill_fraction - Params::default().media.fill_fraction).abs() < 1e-9
        );
    }

    #[test]
    fn new_seeds_the_expected_number_of_balls() {
        let params = Params::default();
        let expected = params.effective_media().ball_count as usize;
        let sim = Simulation::new(params).unwrap();
        assert_eq!(sim.balls().len(), expected);
    }

    #[test]
    fn step_moves_balls_under_gravity_and_keeps_them_inside_the_drum() {
        let mut params = Params::default();
        params.mill.speed_value = 0.0; // isolate settling under gravity from wall-driven motion
        params.simulation.max_balls = 100;
        let mut sim = Simulation::new(params).unwrap();
        let start: Vec<_> = sim.balls().x.clone();

        for _ in 0..120 {
            sim.step(1.0 / 60.0);
        }

        assert_ne!(
            sim.balls().x,
            start,
            "balls should have moved under gravity"
        );
        let r = sim.balls().radius;
        let radius_m = params.mill.radius_m();
        for &p in &sim.balls().x {
            assert!(
                p.length() + r <= radius_m * 1.05,
                "ball escaped the drum: {p:?}"
            );
        }
    }

    #[test]
    fn default_simulation_seeds_fluid_and_extracts_a_surface_after_stepping() {
        let mut params = Params::default();
        // Keep it small so this integration test runs quickly.
        params.simulation.max_balls = 100;
        let mut sim = Simulation::new(params).unwrap();
        assert!(
            !sim.fluid().x.is_empty(),
            "default params have slurry.enabled = true"
        );

        for _ in 0..30 {
            sim.step(1.0 / 60.0);
        }

        for &p in &sim.fluid().x {
            assert!(
                p.x.is_finite() && p.y.is_finite(),
                "non-finite fluid position: {p:?}"
            );
        }
        let surface = sim.fluid_surface();
        assert!(
            !surface.is_empty(),
            "expected at least one free-surface contour"
        );
    }
}
