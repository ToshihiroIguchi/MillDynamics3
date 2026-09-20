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

/// Time constant (s) for the exponential moving average [`Simulation`] applies to the per-sub-step
/// grinding diagnostics (power draw, torque, collision rate, dissipated power, impact-energy
/// histogram, docs/PLAN.md ss3.5). A single fixed 1/240 s sub-step is far too short/noisy to read
/// directly; ~1 s smooths that out while still tracking a real change (e.g. a speed change)
/// within about a second.
const GRINDING_STATS_EMA_TAU_S: f32 = 1.0;

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
    /// Fluid solver diagnostics (viscosity CG iterations, mean shear rate) from the most
    /// recently completed sub-step. `Default` (all zero) before the first [`Simulation::step`]
    /// call.
    fluid_stats: pbf::FluidStepStats,
    /// EMA-smoothed grinding diagnostics (docs/PLAN.md ss3.5), updated every sub-step in
    /// [`Simulation::step`] with time constant [`GRINDING_STATS_EMA_TAU_S`]. All zero before the
    /// first sub-step.
    power_draw_w: f32,
    torque_nm: f32,
    collision_rate_per_s: f32,
    dissipated_power_w: f32,
    impact_energy_counts_per_s: [f32; dem::IMPACT_ENERGY_HISTOGRAM_BINS],
    /// Largest [`dem::DemStepStats::max_substep_displacement_over_diameter`] seen since the last
    /// call to [`Simulation::step`] (reset at the start of each call, not EMA-smoothed -- a
    /// tunnelling-risk spike is exactly the kind of transient an average would hide).
    max_substep_displacement_over_diameter: f32,
    /// Sum of [`pbf::FluidStepStats::coupling_clamp_hits`] over every sub-step of the most recent
    /// [`Simulation::step`] call (reset at the start of each call, not EMA-smoothed -- see
    /// [`coupling::CouplingImpulses::clamp_hits`]'s doc comment for why this is meant to read as
    /// "did this happen at all just now", not a smoothed rate).
    coupling_clamp_hits: u32,
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
            fluid_stats: pbf::FluidStepStats::default(),
            power_draw_w: 0.0,
            torque_nm: 0.0,
            collision_rate_per_s: 0.0,
            dissipated_power_w: 0.0,
            impact_energy_counts_per_s: [0.0; dem::IMPACT_ENERGY_HISTOGRAM_BINS],
            max_substep_displacement_over_diameter: 0.0,
            coupling_clamp_hits: 0,
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

    /// Conjugate-gradient iterations the implicit viscosity solve used on the most recently
    /// completed sub-step ([`pbf::FluidStepStats::viscosity_iterations`]). `0` before the first
    /// [`Simulation::step`] call.
    pub fn viscosity_iterations(&self) -> u32 {
        self.fluid_stats.viscosity_iterations
    }

    /// Mean shear rate (1/s) over the fluid population on the most recently completed sub-step
    /// ([`pbf::FluidStepStats::mean_shear_rate_per_s`]). `0.0` before the first
    /// [`Simulation::step`] call.
    pub fn mean_shear_rate_per_s(&self) -> f32 {
        self.fluid_stats.mean_shear_rate_per_s
    }

    /// EMA-smoothed mill power draw (W per metre of mill depth): the rate of work the drum's
    /// motor must supply against the charge's frictional resistance ([`dem::DemStepStats::wall_work_j`]).
    pub fn power_draw_w(&self) -> f32 {
        self.power_draw_w
    }

    /// EMA-smoothed drum torque (N*m per metre of mill depth), `power_draw_w / omega`. `0.0`
    /// while the drum is not rotating (torque is undefined, not infinite, at `omega = 0`).
    pub fn torque_nm(&self) -> f32 {
        self.torque_nm
    }

    /// EMA-smoothed ball-ball + ball-wall collision rate (impacts per second per metre of mill
    /// depth), see [`dem::DemStepStats::collision_count`].
    pub fn collision_rate_per_s(&self) -> f32 {
        self.collision_rate_per_s
    }

    /// EMA-smoothed net mechanical-energy dissipation rate (W per metre of mill depth) from the
    /// DEM contact solve, see [`dem::DemStepStats::dissipated_energy_j`] -- net of the wall's own
    /// work input, so this is non-negative in practice (the raw per-sub-step value is clamped to
    /// `>= 0` before feeding the EMA, guarding only against float round-off in a near-zero
    /// sub-step; the invariant the metric relies on already bounds it there). Excludes fluid-side
    /// viscous dissipation.
    pub fn dissipated_power_w(&self) -> f32 {
        self.dissipated_power_w
    }

    /// EMA-smoothed impact rate (impacts per second per metre of mill depth) in each log-spaced
    /// energy bin (see [`dem::impact_energy_bin_edges`]).
    pub fn impact_energy_counts_per_s(&self) -> &[f32; dem::IMPACT_ENERGY_HISTOGRAM_BINS] {
        &self.impact_energy_counts_per_s
    }

    /// Largest [`dem::DemStepStats::max_substep_displacement_over_diameter`] over every sub-step
    /// of the most recent [`Simulation::step`] call (not EMA-smoothed, see that field's struct
    /// doc comment for why).
    pub fn max_substep_displacement_over_diameter(&self) -> f32 {
        self.max_substep_displacement_over_diameter
    }

    /// Sum of [`pbf::FluidStepStats::coupling_clamp_hits`] over every sub-step of the most recent
    /// [`Simulation::step`] call (not EMA-smoothed, see
    /// [`coupling::CouplingImpulses::clamp_hits`]'s doc comment for why).
    pub fn coupling_clamp_hits(&self) -> u32 {
        self.coupling_clamp_hits
    }

    /// Extracts the fluid's free-surface contour(s) at the current state (docs/PLAN.md ss3.5).
    /// Recomputed on demand each call -- callers should not call this more often than needed for
    /// rendering.
    pub fn fluid_surface(&self) -> Vec<surface::Polygon> {
        let drum = Drum::new(
            self.params.mill.radius_m(),
            self.params.mill.omega(),
            self.params.lifters,
        );
        surface::extract_surface(&self.fluid, &drum, self.drum_angle)
    }

    /// Computes derived metrics (toe/shoulder, pool extent, mixing index, grinding/solver
    /// diagnostics, debug checks, docs/PLAN.md ss3.5) at the current state. The geometric/
    /// physical fields (toe/shoulder, pool, mixing, ...) are recomputed on demand from the
    /// current ball/fluid state each call; the grinding/solver diagnostics (power draw,
    /// collision histogram, coarse-graining, ...) are this instance's own accumulated state at
    /// call time (see [`metrics::Metrics`]'s doc comment). Callers should not call this more
    /// often than needed (e.g. once per rendered HUD update).
    pub fn metrics(&self) -> metrics::Metrics {
        let drum = Drum::new(
            self.params.mill.radius_m(),
            self.params.mill.omega(),
            self.params.lifters,
        );
        let mut m = metrics::compute(&self.dem.balls, &self.fluid, &drum, self.drum_angle);

        let effective = self.params.effective_media();
        m.power_draw_w = self.power_draw_w;
        m.torque_nm = self.torque_nm;
        m.collision_rate_per_s = self.collision_rate_per_s;
        m.dissipated_power_w = self.dissipated_power_w;
        m.impact_energy_histogram = metrics::ImpactEnergyHistogram {
            bin_edges_j: dem::impact_energy_bin_edges().to_vec(),
            counts_per_s: self.impact_energy_counts_per_s.to_vec(),
        };
        m.coupling_clamp_hits = self.coupling_clamp_hits;
        m.effective_ball_diameter_m = effective.diameter_m;
        m.simulated_ball_count = effective.ball_count;
        m.true_ball_count = self.params.true_ball_count().round().max(0.0) as u32;
        m.coarse_graining_factor = effective.scale_factor;
        m.max_substep_displacement_over_diameter = self.max_substep_displacement_over_diameter;
        m.mean_shear_rate_per_s = self.fluid_stats.mean_shear_rate_per_s;
        m.viscosity_solver_iterations = self.fluid_stats.viscosity_iterations;
        m
    }

    /// The fixed sub-step size (s) this simulation's `simulation.substeps` implies at the
    /// project's nominal 60 Hz target frame rate: `1 / (60 * substeps)`. [`FIXED_DT`] is this
    /// value only at `substeps = 4`; the current default is `substeps = 8`
    /// (`params::SimulationParams::default`), so `FIXED_DT` (used only by `benches/step.rs` and
    /// docs, never by the solver itself) is `1/240 s` while this method's actual default-params
    /// value is `1/480 s`. This method gives the correct value for whatever `substeps` this
    /// instance's params currently specify. A caller driving a fixed-sub-step accumulator loop
    /// (docs/PLAN.md ss4.1: "run `step_fixed` until sim time catches up with wall time, capped by
    /// a frame budget") should use this rather than hard-coding [`FIXED_DT`], since it stays
    /// correct if `substeps` is changed from the default.
    pub fn fixed_sub_dt(&self) -> f32 {
        1.0 / (60.0 * self.params.simulation.substeps.max(1) as f32)
    }

    /// Advances the simulation by exactly one fixed sub-step ([`Simulation::fixed_sub_dt`]).
    /// Unlike [`Simulation::step`], this does **not** reset
    /// [`Simulation::max_substep_displacement_over_diameter`]/[`Simulation::coupling_clamp_hits`]
    /// -- a caller running several `step_fixed` calls in a row to catch up with wall time within
    /// one rendered frame (docs/PLAN.md ss4.1) should call [`Simulation::reset_frame_stats`] once
    /// at the start of that frame instead, so those "peak/any-hit this frame" diagnostics
    /// (see their own doc comments for why) cover the whole frame rather than just its last
    /// sub-step.
    pub fn step_fixed(&mut self) {
        let sub_dt = self.fixed_sub_dt();
        self.advance_one_sub_step(sub_dt);
    }

    /// Resets the per-frame diagnostics [`Simulation::max_substep_displacement_over_diameter`]
    /// and [`Simulation::coupling_clamp_hits`] to zero. [`Simulation::step`] does this itself
    /// once per call; a caller instead driving [`Simulation::step_fixed`] in an accumulator loop
    /// (docs/PLAN.md ss4.1) should call this once per rendered frame, before that frame's
    /// `step_fixed` calls, so the two calling conventions report the same "since the last
    /// [`Simulation::metrics`]-worthy checkpoint" semantics.
    pub fn reset_frame_stats(&mut self) {
        self.max_substep_displacement_over_diameter = 0.0;
        self.coupling_clamp_hits = 0;
    }

    /// Advances the simulation by `dt` seconds of simulation time, split into
    /// `simulation.substeps` fixed sub-steps (docs/PLAN.md ss3.2/4.1).
    pub fn step(&mut self, dt: f32) {
        let n_substeps = self.params.simulation.substeps.max(1);
        let sub_dt = dt / n_substeps as f32;

        self.reset_frame_stats();
        for _ in 0..n_substeps {
            self.advance_one_sub_step(sub_dt);
        }
    }

    /// The shared per-sub-step advance both [`Simulation::step`] and [`Simulation::step_fixed`]
    /// are built from: one call to [`coupling::step`] at sub-step size `sub_dt`, folding its
    /// diagnostics into this instance's accumulated/EMA-smoothed state and advancing
    /// `drum_angle`/`sim_time`.
    fn advance_one_sub_step(&mut self, sub_dt: f32) {
        let omega = self.params.mill.omega();
        let radius_m = self.params.mill.radius_m();
        let lifters = self.params.lifters;
        let dem_iterations = self.params.simulation.dem_iterations;
        let pbf_iterations = self.params.simulation.pbf_iterations;

        let drum = Drum::new(radius_m, omega, lifters);
        let (fluid_stats, dem_stats) = coupling::step(
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
        self.coupling_clamp_hits += fluid_stats.coupling_clamp_hits;
        self.fluid_stats = fluid_stats;
        self.update_grinding_stats(&dem_stats, sub_dt, omega);
        self.drum_angle = (self.drum_angle + omega * sub_dt).rem_euclid(std::f32::consts::TAU);
        self.sim_time += sub_dt as f64;
    }

    /// Folds one sub-step's [`dem::DemStepStats`] into the EMA-smoothed grinding diagnostics
    /// (time constant [`GRINDING_STATS_EMA_TAU_S`]), except
    /// `max_substep_displacement_over_diameter`, which tracks this call's peak directly (see its
    /// struct doc comment).
    fn update_grinding_stats(&mut self, stats: &dem::DemStepStats, sub_dt: f32, omega: f32) {
        let alpha = 1.0 - (-sub_dt / GRINDING_STATS_EMA_TAU_S).exp();

        let power_instant = stats.wall_work_j / sub_dt;
        self.power_draw_w += alpha * (power_instant - self.power_draw_w);

        let torque_instant = if omega.abs() > 1e-6 {
            power_instant / omega
        } else {
            0.0
        };
        self.torque_nm += alpha * (torque_instant - self.torque_nm);

        let collision_rate_instant = stats.collision_count as f32 / sub_dt;
        self.collision_rate_per_s += alpha * (collision_rate_instant - self.collision_rate_per_s);

        // Clamped to >= 0: `dissipated_energy_j` is provably non-negative up to floating-point
        // summation-order slack (see its doc comment), so this guards only against that slack
        // surfacing as a visible negative wattage, not against a real solver defect -- a
        // materially negative raw value would still be caught by
        // `dem::tests`/`tests::dissipated_energy_is_never_materially_negative_while_cascading`.
        let dissipated_power_instant = stats.dissipated_energy_j.max(0.0) / sub_dt;
        self.dissipated_power_w += alpha * (dissipated_power_instant - self.dissipated_power_w);

        for (ema, &count) in self
            .impact_energy_counts_per_s
            .iter_mut()
            .zip(&stats.impact_energy_histogram)
        {
            let instant = count as f32 / sub_dt;
            *ema += alpha * (instant - *ema);
        }

        self.max_substep_displacement_over_diameter = self
            .max_substep_displacement_over_diameter
            .max(stats.max_substep_displacement_over_diameter);
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
    fn fixed_sub_dt_matches_60hz_over_substeps() {
        let mut params = Params::default();
        params.simulation.substeps = 8;
        let sim = Simulation::new(params).unwrap();
        assert!((sim.fixed_sub_dt() - 1.0 / 480.0).abs() < 1e-9);
    }

    #[test]
    fn step_fixed_called_substeps_times_matches_one_step_call_at_1x() {
        // `step(1/60)` (a nominal 60 Hz frame at 1x time scale) must advance the simulation
        // identically to calling `step_fixed()` `substeps` times -- the two calling conventions
        // (docs/PLAN.md ss4.1) are meant to be equivalent at that specific `dt`, not just
        // similar, so a worker migrating from one to the other doesn't change simulation
        // behaviour.
        let mut params = Params::default();
        params.media.fill_fraction = 0.0; // isolate kinematics from ball/fluid dynamics
        params.slurry.enabled = false;
        let n_substeps = params.simulation.substeps;

        let mut via_step = Simulation::new(params).unwrap();
        via_step.step(1.0 / 60.0);

        let mut via_fixed = Simulation::new(params).unwrap();
        via_fixed.reset_frame_stats();
        for _ in 0..n_substeps {
            via_fixed.step_fixed();
        }

        assert!((via_step.sim_time() - via_fixed.sim_time()).abs() < 1e-9);
        assert!((via_step.drum_angle() - via_fixed.drum_angle()).abs() < 1e-6);
    }

    #[test]
    fn reset_frame_stats_zeroes_the_per_frame_diagnostics() {
        let mut params = Params::default();
        params.simulation.max_balls = 100;
        let mut sim = Simulation::new(params).unwrap();
        for _ in 0..30 {
            sim.step_fixed();
        }
        assert!(sim.max_substep_displacement_over_diameter() > 0.0);
        sim.reset_frame_stats();
        assert_eq!(sim.max_substep_displacement_over_diameter(), 0.0);
        assert_eq!(sim.coupling_clamp_hits(), 0);
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

    #[test]
    fn steady_cascading_charge_never_gains_more_energy_than_the_wall_supplies() {
        // Replaces a looser "power draw and dissipation are within a couple of orders of
        // magnitude of each other" check (see git history) with the actual invariant this
        // crate's energy accounting must satisfy: total mechanical energy (KE + gravitational
        // PE) can only *increase* by as much as the wall's motor supplies via friction
        // (`dem::DemStepStats::wall_work_j`, the same quantity the "Power draw" metric is an EMA
        // of) -- every other mechanism in the contact solve (Coulomb friction, restitution
        // `e < 1`, rolling friction, and step 3's bounded depenetration recovery, see
        // `dem::MAX_RECOVERY_FRACTION`) is dissipative and can only remove mechanical energy,
        // never add it. Checked every sub-step, not just at the end, and against `DemState`
        // directly (bypassing `Simulation`'s EMA-smoothed `power_draw_w`/`dissipated_power_w`,
        // whose smoothing and net-KE-change definitions respectively make them unsuitable for an
        // exact per-sub-step accounting check like this one).
        //
        // This is the direct regression test for this project's fluidised-charge
        // energy-injection bug (docs/PHYSICS.md): before the fix, this invariant was violated by
        // orders of magnitude on essentially every sub-step (step 3's uncapped depenetration
        // recovery, turned into velocity by step 4, with nothing removing the excess).
        use crate::dem::DemState;
        use crate::geometry::Drum;
        use crate::params::{Direction, LiftersParams, MediaParams, MillParams, SpeedMode};

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
        params.simulation.max_balls = 150;
        params.validate().unwrap();

        let effective = params.effective_media();
        let radius_m = params.mill.radius_m();
        let omega = params.mill.omega();
        let mut dem = DemState::new(&effective, radius_m, params.simulation.seed);
        let dt = 1.0 / 240.0;
        let mut drum_angle = 0.0f32;

        // Shares its definition with `DemStepStats::dissipated_energy_j`'s own accounting (see
        // `dem::mechanical_energy_j`'s doc comment) -- translational + rotational KE + PE, not
        // just translational KE, so a legitimate friction-driven spin<->translation transfer
        // isn't misread as spurious energy creation.
        let mechanical_energy_j = |dem: &DemState| -> f32 { dem::mechanical_energy_j(&dem.balls) };

        // Settling phase (not checked): the fresh lattice's initial fall/pile-up is not
        // representative of the invariant this test cares about.
        for _ in 0..(3 * 240) {
            let drum = Drum::new(radius_m, omega, params.lifters);
            dem.step(
                &drum,
                drum_angle,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
        }

        // Measurement phase: the invariant, checked every sub-step.
        let mut e_prev = mechanical_energy_j(&dem);
        let mut cumulative_wall_work_j = 0.0f32;
        // A tiny slack for floating-point summation order, not a physical allowance: measured
        // empirically, the invariant holds with wide margin (worst observed
        // `delta_e - wall_work_j` is comfortably negative, never positive) at this project's
        // defaults.
        let tolerance_j = 1e-3 * dem.balls.mass * (omega * radius_m).powi(2).max(1.0);
        for _ in 0..(5 * 240) {
            let drum = Drum::new(radius_m, omega, params.lifters);
            let stats = dem.step(
                &drum,
                drum_angle,
                &params.media,
                params.simulation.dem_iterations,
                dt,
            );
            drum_angle = (drum_angle + omega * dt).rem_euclid(std::f32::consts::TAU);
            cumulative_wall_work_j += stats.wall_work_j;

            let e_now = mechanical_energy_j(&dem);
            let delta_e = e_now - e_prev;
            assert!(
                delta_e <= stats.wall_work_j + tolerance_j,
                "mechanical energy increased more than the wall's own work this sub-step \
                 (energy injected from nowhere): delta_e={delta_e} wall_work_j={}",
                stats.wall_work_j
            );
            e_prev = e_now;
        }

        assert!(
            cumulative_wall_work_j > 0.0,
            "expected the wall to do positive net work while cascading, got {cumulative_wall_work_j}"
        );
    }

    #[test]
    fn centrifuging_draws_less_power_than_cascading() {
        // A well-known operational signature of centrifuging (the failure mode docs/PLAN.md
        // ss3.2's cascading-vs-centrifuging distinction warns about): once the charge locks to
        // the wall and co-rotates, the mill draws meaningfully less power than when it is
        // actively cascading -- there is far less ongoing relative sliding at the wall to
        // overcome, and none of the grinding action a real mill is run for. This qualitative
        // ordering (rather than an absolute "power should approach zero" claim) is what real
        // mill operators actually use to detect an over-speed centrifuging mill, and is more
        // robust here than pinning down exactly how close settled centrifuging power gets to
        // zero for this particular ball count/size/friction combination.
        use crate::params::{Direction, LiftersParams, MediaParams, MillParams, SpeedMode};

        fn settled_power_draw(percent_critical: f32, fill_fraction: f32) -> f32 {
            let mut params = Params {
                mill: MillParams {
                    diameter_m: 1.0,
                    speed_mode: SpeedMode::PercentCritical,
                    speed_value: percent_critical,
                    direction: Direction::CounterClockwise,
                },
                media: MediaParams {
                    ball_diameter_m: 0.02,
                    fill_fraction,
                    ..MediaParams::default()
                },
                lifters: LiftersParams {
                    count: 0,
                    ..LiftersParams::default()
                },
                slurry: params::SlurryParams {
                    enabled: false,
                    ..params::SlurryParams::default()
                },
                ..Params::default()
            };
            params.simulation.max_balls = 150;
            params.validate().unwrap();

            let mut sim = Simulation::new(params).unwrap();
            let dt = 1.0 / 60.0;
            for _ in 0..(8 * 60) {
                sim.step(dt);
            }
            sim.power_draw_w()
        }

        let cascading_power = settled_power_draw(70.0, 0.25);
        let centrifuging_power = settled_power_draw(130.0, 0.10);
        assert!(
            cascading_power > 0.0,
            "expected positive power draw while cascading, got {cascading_power}"
        );
        assert!(
            centrifuging_power < cascading_power,
            "centrifuging should draw less power than cascading: \
             cascading={cascading_power} centrifuging={centrifuging_power}"
        );
    }

    #[test]
    fn dissipated_power_approaches_power_draw_in_a_settled_dry_charge() {
        // Regression for the same review finding as
        // `dem::tests::dissipated_energy_is_never_materially_negative_while_cascading`, checked
        // at the `Simulation`-level EMA'd metric instead of the raw per-sub-step value: once a
        // dry (no slurry) cascading charge reaches a statistically steady state, its total
        // mechanical energy is on average constant, so the wall's power input (`power_draw_w`,
        // energy the motor supplies) and the contact solve's net dissipation
        // (`dissipated_power_w`) must track each other -- this is exactly the check that would
        // have caught the original defect (which instead separated the two by a
        // `wall_work_j`-sized negative offset, see `dem::DemStepStats::dissipated_energy_j`'s
        // doc comment).
        use crate::params::{Direction, LiftersParams, MediaParams, MillParams, SpeedMode};

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
            slurry: params::SlurryParams {
                enabled: false,
                ..params::SlurryParams::default()
            },
            ..Params::default()
        };
        params.simulation.max_balls = 150;
        params.validate().unwrap();

        let mut sim = Simulation::new(params).unwrap();
        let dt = 1.0 / 60.0;
        // Settle well past the EMA time constant (`GRINDING_STATS_EMA_TAU_S = 1.0` s) so both
        // readings reflect steady state, not the initial-landing transient.
        for _ in 0..(10 * 60) {
            sim.step(dt);
        }

        let power_draw = sim.power_draw_w();
        let dissipated = sim.dissipated_power_w();
        assert!(
            power_draw > 0.0,
            "expected positive power draw for a settled cascading charge, got {power_draw}"
        );
        assert!(
            dissipated > 0.0,
            "expected positive dissipated power for a settled cascading charge, got {dissipated}"
        );
        let ratio = dissipated / power_draw;
        assert!(
            (0.5..=1.5).contains(&ratio),
            "dissipated_power_w should track power_draw_w within a loose band once a dry \
             charge is statistically steady (energy in ~= energy out): \
             power_draw_w={power_draw}, dissipated_power_w={dissipated}, ratio={ratio}"
        );
    }

    #[test]
    fn metrics_surfaces_grinding_and_coarse_graining_fields_after_stepping() {
        // `Simulation::metrics()` must fill in every field `metrics::compute` itself cannot
        // derive (see that struct's doc comment): ball-count/coarse-graining fields from
        // `Params`, and the EMA-smoothed grinding/solver diagnostics from this instance's own
        // accumulated state. Default params coarse-grain heavily (see
        // `params::tests::effective_media_coarse_grains_when_true_count_is_large`), so this also
        // exercises that path.
        let mut params = Params::default();
        params.simulation.max_balls = 150;
        let mut sim = Simulation::new(params).unwrap();

        for _ in 0..(3 * 60) {
            sim.step(1.0 / 60.0);
        }
        let m = sim.metrics();

        assert!(
            m.true_ball_count > m.simulated_ball_count,
            "expected coarse-graining to have kicked in: true={} sim={}",
            m.true_ball_count,
            m.simulated_ball_count
        );
        assert!(m.coarse_graining_factor > 1.0);
        assert!(m.effective_ball_diameter_m > 0.0);
        assert_eq!(m.simulated_ball_count, sim.balls().len() as u32);
        assert_eq!(m.fluid_particle_count, sim.fluid().len() as u32);

        assert!(m.power_draw_w.is_finite());
        assert!(m.torque_nm.is_finite());
        assert!(m.collision_rate_per_s.is_finite() && m.collision_rate_per_s >= 0.0);
        assert!(m.dissipated_power_w.is_finite());
        assert!(m.mean_shear_rate_per_s.is_finite() && m.mean_shear_rate_per_s >= 0.0);

        assert_eq!(
            m.impact_energy_histogram.bin_edges_j.len(),
            dem::IMPACT_ENERGY_HISTOGRAM_BINS + 1
        );
        assert_eq!(
            m.impact_energy_histogram.counts_per_s.len(),
            dem::IMPACT_ENERGY_HISTOGRAM_BINS
        );
        assert!(
            m.impact_energy_histogram
                .bin_edges_j
                .windows(2)
                .all(|w| w[1] > w[0]),
            "bin edges should be strictly increasing"
        );

        // Round-trips through JSON cleanly (this is exactly the wasm/frontend boundary).
        let json = serde_json::to_string(&m).expect("Metrics should serialize");
        assert!(json.contains("power_draw_w"));
        assert!(json.contains("impact_energy_histogram"));
    }
}
