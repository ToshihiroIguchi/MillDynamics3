//! mill-core: 2D ball-mill DEM + PBF slurry simulation core.
//!
//! Pure Rust, no WASM dependencies -- kept independently testable/benchmarkable via `cargo test`
//! and `cargo bench`. The `wasm-bindgen` wrapper lives in the sibling `mill-wasm` crate.
//!
//! Module responsibilities (see docs/PLAN.md for the full design):
//! - [`params`]: parameter groups, defaults, validation, coarse-graining derivation.
//! - [`geometry`]: drum wall SDF (circle + optional lifters), normals, wall velocity.
//! - [`grid`]: uniform-grid spatial hash shared by DEM and PBF neighbour search (M1/M3).
//! - [`dem`]: soft-sphere DEM for grinding media (M1).
//! - [`pbf`]: Position Based Fluids solver for the slurry (M3).
//! - [`coupling`]: two-way ball<->fluid momentum exchange (M4).
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

pub use params::{EffectiveMedia, Params};

/// Fixed simulation sub-step used by the PBF/DEM solvers once they land (M1+): 1/240 s. Both
/// `dem` and `pbf` are designed around the same nominal sub-step rate (docs/PLAN.md ss3.2/3.3).
pub const FIXED_DT: f32 = 1.0 / 240.0;

/// The simulation instance. This is the single type the `mill-wasm` wrapper (and, indirectly, the
/// web worker's fixed-step accumulator, docs/PLAN.md ss4.1) drives via [`Simulation::step`].
///
/// Through milestone M0 this only tracks drum kinematics (no media or slurry yet); DEM/PBF
/// sub-stepping is added inside `step` in later milestones without changing this type's public
/// API.
pub struct Simulation {
    params: Params,
    /// Drum rotation angle (radians), kept wrapped into `[0, 2*pi)`.
    drum_angle: f32,
    /// Total simulated time (s).
    sim_time: f64,
}

impl Simulation {
    /// Creates a new simulation, validating `params` first.
    pub fn new(params: Params) -> Result<Self, String> {
        params.validate()?;
        Ok(Self {
            params,
            drum_angle: 0.0,
            sim_time: 0.0,
        })
    }

    pub fn params(&self) -> &Params {
        &self.params
    }

    /// Replaces the parameters in place, without resetting `drum_angle`/`sim_time`. Used for
    /// "hot" parameter updates that should not restart the simulation (docs/PLAN.md ss4.1); the
    /// worker is responsible for deciding whether a given change needs a full reset (i.e.
    /// constructing a new [`Simulation`]) instead of calling this.
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

    /// Advances the simulation by `dt` seconds of simulation time.
    ///
    /// Through M0 this only advances drum kinematics; DEM/PBF sub-stepping is added in later
    /// milestones without changing this signature.
    pub fn step(&mut self, dt: f32) {
        let omega = self.params.mill.omega();
        self.drum_angle = (self.drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
        self.sim_time += dt as f64;
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
        let mut sim = Simulation::new(params).unwrap();
        sim.step(0.25); // quarter turn
        assert!((sim.drum_angle() - std::f32::consts::FRAC_PI_2).abs() < 1e-3);
    }

    #[test]
    fn step_wraps_angle_into_0_tau() {
        let mut params = Params::default();
        params.mill.speed_mode = SpeedMode::Rpm;
        params.mill.speed_value = 60.0;
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
}
