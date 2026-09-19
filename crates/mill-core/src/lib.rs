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

    /// EMA-smoothed net kinetic-energy dissipation rate (W per metre of mill depth) from the DEM
    /// contact solve, see [`dem::DemStepStats::dissipated_energy_j`].
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

        self.max_substep_displacement_over_diameter = 0.0;
        self.coupling_clamp_hits = 0;
        for _ in 0..n_substeps {
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

        let dissipated_power_instant = stats.dissipated_energy_j / sub_dt;
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
    fn steady_cascading_charge_has_a_plausible_power_energy_balance() {
        // Power draw (work done against wall friction, docs/PLAN.md ss3.5) and net dissipation
        // (kinetic energy removed by the contact solve) both being positive and finite is a
        // basic sanity check on this instrumentation: energy flows in only through wall friction
        // and gravity (net gravitational work is zero over a full toe-shoulder cycle once the
        // charge's average height stops changing) and leaves only through dissipation, so both
        // must be positive once the charge is genuinely cascading, and dissipation cannot run
        // away to some absurd multiple of the input.
        //
        // This deliberately does *not* assert they are numerically close: the XPBD contact
        // solver (docs/PLAN.md ss3.2) resolves many simultaneous, tightly-packed contacts in the
        // toe with a modest fixed iteration count (`dem_iterations = 4`, not run to full
        // convergence every sub-step) and treats any contact below
        // `RESTITUTION_VELOCITY_THRESHOLD` as inelastic by construction (no compensating
        // restitution impulse) -- both are well-known sources of extra *numerical* damping in
        // position-based/iterative rigid-contact solvers, on top of the physical friction and
        // restitution losses `wall_work_j` alone accounts for. Measured empirically at this
        // project's defaults: dissipation runs roughly an order of magnitude above wall power
        // draw, which is why the bound below is wide rather than a tight ratio.
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
        // Settle into steady cascading motion for several EMA time constants
        // (GRINDING_STATS_EMA_TAU_S = 1 s) so the smoothed readings have converged.
        for _ in 0..(6 * 60) {
            sim.step(dt);
        }

        let power = sim.power_draw_w();
        let dissipated = sim.dissipated_power_w();
        assert!(
            power > 0.0,
            "expected positive power draw once cascading, got {power}"
        );
        assert!(
            dissipated > 0.0,
            "expected positive dissipation once cascading, got {dissipated}"
        );
        let ratio = power / dissipated;
        assert!(
            (0.01..1.0).contains(&ratio),
            "power draw and dissipation should be positive and within a couple of orders of \
             magnitude of each other, not disconnected or with dissipation *below* input power: \
             power={power} dissipated={dissipated} ratio={ratio}"
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
