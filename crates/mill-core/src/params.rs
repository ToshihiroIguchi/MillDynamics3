//! Simulation parameters.
//!
//! All quantities here are SI (meters, kilograms, seconds, Pascal-seconds, radians) regardless of
//! what unit the UI displays them in (see `web/src/params/schema.ts` for UI-facing units and
//! conversions). `Display`-only toggles (what to draw) are a rendering concern and intentionally
//! do not live here.

use serde::{Deserialize, Serialize};

/// How the mill's rotation speed is specified.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpeedMode {
    /// `speed_value` is an absolute rotation speed in rpm.
    Rpm,
    /// `speed_value` is a percentage of critical speed (e.g. `70.0` for 70% Nc).
    PercentCritical,
}

/// Drum rotation direction, viewed from the standard (right-handed, +x right / +y up) 2D frame.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Direction {
    CounterClockwise,
    Clockwise,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MillParams {
    /// Drum inner diameter (m).
    pub diameter_m: f32,
    pub speed_mode: SpeedMode,
    /// Interpretation depends on `speed_mode`: rpm, or % of critical speed.
    pub speed_value: f32,
    pub direction: Direction,
}

impl Default for MillParams {
    fn default() -> Self {
        Self {
            diameter_m: 1.0,
            speed_mode: SpeedMode::Rpm,
            // ~70% Nc for the default 1.0 m drum (Nc = 42.3 rpm), expressed directly in rpm.
            speed_value: 30.0,
            direction: Direction::CounterClockwise,
        }
    }
}

impl MillParams {
    /// Critical speed in rpm: the rotation speed at which media theoretically centrifuges
    /// against the wall (Nc = 42.3 / sqrt(D), D in meters).
    pub fn critical_speed_rpm(&self) -> f32 {
        42.3 / self.diameter_m.sqrt()
    }

    /// Rotation speed in rpm, regardless of `speed_mode`.
    pub fn rpm(&self) -> f32 {
        match self.speed_mode {
            SpeedMode::Rpm => self.speed_value,
            SpeedMode::PercentCritical => self.critical_speed_rpm() * self.speed_value / 100.0,
        }
    }

    /// Rotation speed as a percentage of critical speed, regardless of `speed_mode`.
    pub fn percent_critical(&self) -> f32 {
        match self.speed_mode {
            SpeedMode::PercentCritical => self.speed_value,
            SpeedMode::Rpm => 100.0 * self.speed_value / self.critical_speed_rpm(),
        }
    }

    /// Angular velocity magnitude (rad/s).
    pub fn omega_magnitude(&self) -> f32 {
        self.rpm() * 2.0 * std::f32::consts::PI / 60.0
    }

    /// Signed angular velocity (rad/s): positive is counter-clockwise.
    pub fn omega(&self) -> f32 {
        match self.direction {
            Direction::CounterClockwise => self.omega_magnitude(),
            Direction::Clockwise => -self.omega_magnitude(),
        }
    }

    pub fn radius_m(&self) -> f32 {
        self.diameter_m * 0.5
    }

    pub fn validate(&self) -> Result<(), String> {
        if !(self.diameter_m > 0.05 && self.diameter_m < 20.0) {
            return Err("mill.diameter_m must be in (0.05, 20.0)".into());
        }
        if !(self.speed_value.is_finite() && self.speed_value >= 0.0) {
            return Err("mill.speed_value must be finite and >= 0".into());
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MediaParams {
    /// Ball diameter (m). A future size-distribution table can extend this.
    pub ball_diameter_m: f32,
    /// Fraction of the drum's cross-sectional area occupied by the charge *including voids*
    /// (the conventional "J" ball-filling fraction). The solid disc area actually seeded is
    /// `fill_fraction * packing_fraction_2d` of the drum area, see
    /// [`Params::effective_media`].
    pub fill_fraction: f32,
    /// 2D areal packing fraction of the settled disc charge (solid area / charge footprint area),
    /// used to convert `fill_fraction` into a disc count. Random close packing of equal discs is
    /// ~0.82; the hexagonal upper bound is `pi / (2 sqrt 3) ~= 0.907`.
    pub packing_fraction_2d: f32,
    /// Media density (kg/m^3). Balls are modelled as unit-depth discs (see
    /// [`crate::dem::ball_mass`]), so this is also the areal mass per metre of mill length.
    pub density_kg_m3: f32,
    pub restitution_ball_ball: f32,
    pub restitution_ball_wall: f32,
    pub friction_ball_ball: f32,
    pub friction_ball_wall: f32,
    /// Rolling-resistance coefficient (dimensionless torque coefficient).
    pub rolling_friction: f32,
}

impl Default for MediaParams {
    // Defaults model yttria-stabilized zirconia (ZrO2) grinding media, a common ceramic ball-mill
    // charge: density ~6.0 g/cm^3, harder and more elastic (higher restitution, lower friction)
    // than the steel media this used to default to. 10 mm is a typical lab/bench tumbling-mill
    // ball size (as opposed to sub-mm beads used in attritor/bead mills).
    fn default() -> Self {
        Self {
            ball_diameter_m: 0.010,
            fill_fraction: 0.30,
            packing_fraction_2d: 0.82,
            density_kg_m3: 6000.0,
            restitution_ball_ball: 0.7,
            restitution_ball_wall: 0.5,
            friction_ball_ball: 0.25,
            friction_ball_wall: 0.35,
            rolling_friction: 0.01,
        }
    }
}

impl MediaParams {
    pub fn validate(&self) -> Result<(), String> {
        if !(self.ball_diameter_m > 0.0005 && self.ball_diameter_m < 1.0) {
            return Err("media.ball_diameter_m must be in (0.0005, 1.0)".into());
        }
        if !(0.0..=0.9).contains(&self.fill_fraction) {
            return Err("media.fill_fraction must be in [0, 0.9]".into());
        }
        if !(0.5..=0.907).contains(&self.packing_fraction_2d) {
            return Err("media.packing_fraction_2d must be in [0.5, 0.907]".into());
        }
        if !(self.density_kg_m3.is_finite() && self.density_kg_m3 > 0.0) {
            return Err("media.density_kg_m3 must be > 0".into());
        }
        for (name, e) in [
            ("restitution_ball_ball", self.restitution_ball_ball),
            ("restitution_ball_wall", self.restitution_ball_wall),
        ] {
            if !(0.0..=1.0).contains(&e) {
                return Err(format!("media.{name} must be in [0, 1]"));
            }
        }
        for (name, mu) in [
            ("friction_ball_ball", self.friction_ball_ball),
            ("friction_ball_wall", self.friction_ball_wall),
            ("rolling_friction", self.rolling_friction),
        ] {
            if !(mu >= 0.0 && mu.is_finite()) {
                return Err(format!("media.{name} must be >= 0"));
            }
        }
        Ok(())
    }
}

/// Slurry rheology model. Only `Newtonian` is implemented through M6; `Bingham` is reserved for
/// the M7 fidelity extension (see docs/PLAN.md).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Rheology {
    Newtonian,
    Bingham,
}

/// Initial dye tracer pattern, used to visualize and measure mixing.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DyePattern {
    LeftRight,
    TopBottom,
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SlurryParams {
    pub enabled: bool,
    /// Fraction of the drum's cross-sectional area filled with slurry.
    pub fill_fraction: f32,
    pub density_kg_m3: f32,
    pub viscosity_pa_s: f32,
    pub rheology: Rheology,
    /// Bingham yield stress (Pa), used only when `rheology == Bingham`.
    pub yield_stress_pa: f32,
    /// No-slip blend factor at the drum wall/lifters, in [0, 1] (1 = full no-slip).
    pub wall_no_slip: f32,
    /// No-slip blend factor at ball surfaces, in [0, 1] (1 = full no-slip).
    pub ball_no_slip: f32,
    /// Ball<->slurry wettability, in [0, 1]. `0` = non-wetting (fluid is never pulled toward a
    /// ball's surface, only ever pushed off it -- a 180 degree contact angle). `1` = strongly
    /// wetting. Scales [`crate::pbf::FluidParticles::step_coupled`]'s Akinci-style boundary-
    /// adhesion term. Together with `surface_tension_n_m` this gives a genuine (if empirically
    /// calibrated, not first-principles-derived -- see that term's own doc comment) contact angle:
    /// wetting pulls fluid onto media, cohesion holds the resulting film together against gravity,
    /// and the two compete exactly as they do physically. Before `surface_tension_n_m` existed,
    /// wettability alone could only fix whether fluid clings to media at all, not how (see
    /// docs/PHYSICS.md ss6.1a's history).
    pub wettability: f32,
    /// Fluid-fluid surface tension (N/m), driving an Akinci-style pairwise cohesion force. `0`
    /// disables it. Real water is ~0.072 N/m (this field's suggested default); like the wider SPH
    /// surface-tension literature this coefficient is a calibrated proportionality to the physical
    /// value, not a first-principles unit conversion -- see [`crate::pbf`]'s
    /// `REFERENCE_SURFACE_TENSION_N_M`/`COHESION_ACCEL_FACTOR` doc comments for why no such
    /// conversion exists even in the original method (Akinci, Akinci & Teschner, "Versatile
    /// Surface Tension and Adhesion for SPH Fluids", 2013). Replaces the disabled `S_CORR_K`
    /// artificial-pressure term as this solver's only fluid-fluid attraction.
    ///
    /// `#[serde(default = ...)]`, not a bare `#[serde(default)]`: an older client's JSON
    /// (predating this field, e.g. the web UI before it exposes a control for this) must not
    /// silently deserialize to `0.0` (cohesion off) while `wettability` still defaults to a
    /// positive value -- that specific combination is exactly the "adhesion with nothing holding
    /// the resulting film together" imbalance this field exists to fix (see `wettability`'s doc
    /// comment). Falling back to this struct's own intended default (`0.072`) instead keeps an
    /// older client's behaviour physically coherent until it is updated to send the field itself.
    #[serde(default = "default_surface_tension_n_m")]
    pub surface_tension_n_m: f32,
    pub dye_pattern: DyePattern,
}

/// `#[serde(default = ...)]` target for [`SlurryParams::surface_tension_n_m`] -- see that field's
/// doc comment for why this must match `SlurryParams::default()`'s own value rather than `f32`'s
/// bare `0.0`.
fn default_surface_tension_n_m() -> f32 {
    0.072
}

impl Default for SlurryParams {
    fn default() -> Self {
        Self {
            enabled: true,
            // Above `MediaParams::fill_fraction`'s default (0.30) so the liquid's free surface
            // clears the top of the settled ball heap at rest, instead of leaving media exposed
            // above the slurry pool (see docs/PARAMETERS.md's Slurry table).
            fill_fraction: 0.35,
            density_kg_m3: 1800.0,
            viscosity_pa_s: 50.0,
            rheology: Rheology::Newtonian,
            yield_stress_pa: 0.0,
            wall_no_slip: 1.0,
            ball_no_slip: 1.0,
            // Real mineral slurries wet ceramic/steel grinding media reasonably well; 0.6 gives a
            // visible clinging film without dominating the momentum budget (see pbf.rs's
            // ADHESION_ACCEL_FACTOR, scaled by gravity, for the bound).
            wettability: 0.6,
            // Real water's surface tension; see this field's doc comment for why the mapping into
            // this solver's internal cohesion strength is calibrated, not derived.
            surface_tension_n_m: default_surface_tension_n_m(),
            dye_pattern: DyePattern::LeftRight,
        }
    }
}

impl SlurryParams {
    pub fn validate(&self) -> Result<(), String> {
        if !(0.0..=0.9).contains(&self.fill_fraction) {
            return Err("slurry.fill_fraction must be in [0, 0.9]".into());
        }
        if !(self.density_kg_m3.is_finite() && self.density_kg_m3 > 0.0) {
            return Err("slurry.density_kg_m3 must be > 0".into());
        }
        if !(self.viscosity_pa_s >= 0.0 && self.viscosity_pa_s.is_finite()) {
            return Err("slurry.viscosity_pa_s must be >= 0".into());
        }
        if !(self.yield_stress_pa >= 0.0 && self.yield_stress_pa.is_finite()) {
            return Err("slurry.yield_stress_pa must be >= 0".into());
        }
        if !(self.surface_tension_n_m >= 0.0 && self.surface_tension_n_m.is_finite()) {
            return Err("slurry.surface_tension_n_m must be >= 0".into());
        }
        for (name, b) in [
            ("wall_no_slip", self.wall_no_slip),
            ("ball_no_slip", self.ball_no_slip),
            ("wettability", self.wettability),
        ] {
            if !(0.0..=1.0).contains(&b) {
                return Err(format!("slurry.{name} must be in [0, 1]"));
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LiftersParams {
    /// Number of evenly-spaced lifter bars. `0` (the default) means a perfectly smooth wall.
    pub count: u32,
    pub height_m: f32,
    pub base_width_m: f32,
    pub top_width_m: f32,
    /// Angular offset of the first lifter from +x, in degrees.
    pub phase_deg: f32,
}

impl Default for LiftersParams {
    fn default() -> Self {
        Self {
            count: 0,
            height_m: 0.020,
            base_width_m: 0.030,
            top_width_m: 0.020,
            phase_deg: 0.0,
        }
    }
}

impl LiftersParams {
    pub fn validate(&self) -> Result<(), String> {
        if self.count > 0 {
            if !(self.height_m > 0.0 && self.height_m.is_finite()) {
                return Err("lifters.height_m must be > 0".into());
            }
            if !(self.base_width_m > 0.0 && self.base_width_m.is_finite()) {
                return Err("lifters.base_width_m must be > 0".into());
            }
            if !(self.top_width_m >= 0.0 && self.top_width_m.is_finite()) {
                return Err("lifters.top_width_m must be >= 0".into());
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SimulationParams {
    /// Number of fluid particles spanning the drum radius; controls PBF spatial resolution.
    pub resolution: u32,
    /// Fixed sub-steps per rendered frame at 1x time scale (nominal PBF step rate is
    /// `substeps * 60` Hz before frame-budget throttling).
    pub substeps: u32,
    /// PBF density-constraint solver iterations per sub-step. A floor, not a fixed count: the
    /// solver runs additional passes beyond this (bounded separately, see
    /// `pbf::ADAPTIVE_MAX_EXTRA_ITERATIONS`) when the density constraint is still unconverged
    /// after this many, so a violent collision doesn't leave a large compression-error residual.
    pub pbf_iterations: u32,
    /// XPBD non-penetration solver iterations per sub-step for the ball population (see
    /// docs/PLAN.md ss3.2). Analogous to `pbf_iterations`; unlike an explicit DEM time step, this
    /// trades contact-resolution accuracy for cost rather than stability.
    pub dem_iterations: u32,
    /// Wall-clock-to-simulation-time multiplier requested by the user (achieved rate may be
    /// lower and is reported back to the UI, see docs/PLAN.md ss4.1).
    pub time_scale: f32,
    /// Maximum wall-clock milliseconds the worker may spend simulating per animation frame.
    pub frame_budget_ms: f32,
    /// Target upper bound on the number of DEM ball particles actually simulated. If the true
    /// media population (derived from `mill`/`media`) would exceed this, the solver transparently
    /// substitutes a coarser (larger, lighter) effective ball population that preserves total
    /// charge mass and cross-sectional fill area — see [`Params::effective_media`] and
    /// docs/PLAN.md ss3.2 ("Coarse-graining"). The default keeps real-time playback achievable per
    /// the M6 performance targets even at small media diameters.
    pub max_balls: u32,
    pub seed: u64,
}

impl Default for SimulationParams {
    fn default() -> Self {
        Self {
            resolution: 40,
            substeps: 8,
            pbf_iterations: 3,
            dem_iterations: 2,
            time_scale: 1.0,
            frame_budget_ms: 12.0,
            max_balls: 600,
            seed: 1,
        }
    }
}

impl SimulationParams {
    pub fn validate(&self) -> Result<(), String> {
        if !(4..=200).contains(&self.resolution) {
            return Err("simulation.resolution must be in [4, 200]".into());
        }
        if !(1..=16).contains(&self.substeps) {
            return Err("simulation.substeps must be in [1, 16]".into());
        }
        if !(1..=20).contains(&self.pbf_iterations) {
            return Err("simulation.pbf_iterations must be in [1, 20]".into());
        }
        if !(1..=20).contains(&self.dem_iterations) {
            return Err("simulation.dem_iterations must be in [1, 20]".into());
        }
        if !(self.time_scale > 0.0 && self.time_scale.is_finite()) {
            return Err("simulation.time_scale must be > 0".into());
        }
        if !(self.frame_budget_ms > 0.0 && self.frame_budget_ms.is_finite()) {
            return Err("simulation.frame_budget_ms must be > 0".into());
        }
        if !(10..=50_000).contains(&self.max_balls) {
            return Err("simulation.max_balls must be in [10, 50000]".into());
        }
        Ok(())
    }
}

/// The complete set of parameters needed to (re)construct a [`crate::Simulation`].
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct Params {
    pub mill: MillParams,
    pub media: MediaParams,
    pub slurry: SlurryParams,
    pub lifters: LiftersParams,
    pub simulation: SimulationParams,
}

impl Params {
    /// Validates every parameter group, returning the first error found.
    pub fn validate(&self) -> Result<(), String> {
        self.mill.validate()?;
        self.media.validate()?;
        self.slurry.validate()?;
        self.lifters.validate()?;
        self.simulation.validate()?;
        Ok(())
    }

    /// True (uncoarsened) disc count implied by the mill/media parameters:
    /// `N_true = fill_fraction * packing_fraction_2d * drum_area / (pi * r_true^2)`, i.e. the
    /// charge footprint `J * drum_area` times the areal packing fraction, divided by one disc's
    /// area. Fractional (not rounded); see [`Params::effective_media`].
    pub fn true_ball_count(&self) -> f32 {
        let r_true = self.media.ball_diameter_m * 0.5;
        if r_true <= 0.0 {
            return 0.0;
        }
        let drum_area = std::f32::consts::PI * self.mill.radius_m().powi(2);
        (self.media.fill_fraction * self.media.packing_fraction_2d * drum_area)
            / (std::f32::consts::PI * r_true * r_true)
    }

    /// Derives the media population actually simulated, applying coarse-graining (particle
    /// scaling) when the true media population would exceed `simulation.max_balls`.
    ///
    /// The true disc count is [`Params::true_ball_count`]. If `N_true > max_balls`, we substitute
    /// `N_sim` larger discs with scale factor `k = sqrt(N_true / max_balls)`:
    /// `diameter_m = true_diameter_m * k`, so `N_sim = N_true / k^2 ~= max_balls`. Because balls
    /// are unit-depth discs ([`crate::dem::ball_mass`], mass `~ r^2`), preserving the total solid
    /// area automatically preserves the total charge mass, so the density is left unchanged
    /// (`density_kg_m3 == true density`).
    ///
    /// This is a standard coarse-grained DEM approximation (see docs/PLAN.md ss3.2): it preserves
    /// bulk charge mass and footprint area but not the true interstitial void structure or
    /// single-collision statistics at the true particle size. When `N_true <= max_balls`,
    /// `scale_factor` is `1.0` and the true media parameters are returned unchanged.
    pub fn effective_media(&self) -> EffectiveMedia {
        let true_diameter_m = self.media.ball_diameter_m;
        let density_kg_m3 = self.media.density_kg_m3;
        let n_true = self.true_ball_count();

        let max_balls = self.simulation.max_balls as f32;
        if n_true > max_balls && max_balls > 0.0 {
            let scale_factor = (n_true / max_balls).sqrt();
            let ball_count = (n_true / (scale_factor * scale_factor)).round().max(1.0) as u32;
            EffectiveMedia {
                true_diameter_m,
                diameter_m: true_diameter_m * scale_factor,
                density_kg_m3,
                ball_count,
                scale_factor,
            }
        } else {
            EffectiveMedia {
                true_diameter_m,
                diameter_m: true_diameter_m,
                density_kg_m3,
                ball_count: n_true.round().max(0.0) as u32,
                scale_factor: 1.0,
            }
        }
    }
}

/// The media population actually handed to the DEM solver, after [`Params::effective_media`]
/// applies coarse-graining (if any). See that method's doc comment for the scaling law.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EffectiveMedia {
    /// True (UI) ball diameter, unchanged (m).
    pub true_diameter_m: f32,
    /// Diameter actually simulated (m). Equals `true_diameter_m` when `scale_factor == 1.0`.
    pub diameter_m: f32,
    /// Density actually simulated (kg/m^3). Always the true media density: with unit-depth disc
    /// masses, preserving solid area preserves total charge mass without rescaling density.
    pub density_kg_m3: f32,
    /// Number of DEM ball particles actually simulated.
    pub ball_count: u32,
    /// Coarse-graining scale factor `k = diameter_m / true_diameter_m` (>= 1.0; `1.0` = no
    /// coarse-graining applied).
    pub scale_factor: f32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_params_are_valid() {
        Params::default().validate().unwrap();
    }

    #[test]
    fn critical_speed_matches_known_value() {
        // A 0.6 m mill: Nc = 42.3 / sqrt(0.6) ~= 54.6 rpm.
        let mill = MillParams {
            diameter_m: 0.6,
            ..MillParams::default()
        };
        assert!((mill.critical_speed_rpm() - 54.61).abs() < 0.1);
    }

    #[test]
    fn percent_critical_round_trips_through_rpm() {
        let mill = MillParams {
            speed_mode: SpeedMode::PercentCritical,
            speed_value: 70.0,
            ..MillParams::default()
        };
        let rpm = mill.rpm();
        let mill_rpm = MillParams {
            speed_mode: SpeedMode::Rpm,
            speed_value: rpm,
            ..mill
        };
        assert!((mill_rpm.percent_critical() - 70.0).abs() < 1e-3);
    }

    #[test]
    fn clockwise_direction_negates_omega() {
        let ccw = MillParams {
            direction: Direction::CounterClockwise,
            ..MillParams::default()
        };
        let cw = MillParams {
            direction: Direction::Clockwise,
            ..ccw
        };
        assert!((ccw.omega() + cw.omega()).abs() < 1e-6);
        assert!(ccw.omega() > 0.0);
    }

    #[test]
    fn invalid_diameter_rejected() {
        let mill = MillParams {
            diameter_m: 0.0,
            ..MillParams::default()
        };
        assert!(mill.validate().is_err());
    }

    #[test]
    fn invalid_fill_fraction_rejected() {
        let media = MediaParams {
            fill_fraction: 1.5,
            ..MediaParams::default()
        };
        assert!(media.validate().is_err());
    }

    #[test]
    fn zero_lifters_needs_no_geometry_bounds() {
        // count == 0 (smooth wall) is valid even with degenerate geometry fields.
        let lifters = LiftersParams {
            count: 0,
            height_m: 0.0,
            base_width_m: 0.0,
            top_width_m: 0.0,
            phase_deg: 0.0,
        };
        assert!(lifters.validate().is_ok());
    }

    #[test]
    fn effective_media_is_unscaled_when_true_count_is_small() {
        // Large media in a small mill: well under max_balls, so no coarse-graining is applied.
        let params = Params {
            mill: MillParams {
                diameter_m: 1.0,
                ..MillParams::default()
            },
            media: MediaParams {
                ball_diameter_m: 0.3,
                fill_fraction: 0.30,
                ..MediaParams::default()
            },
            simulation: SimulationParams {
                max_balls: 2000,
                ..SimulationParams::default()
            },
            ..Params::default()
        };
        let eff = params.effective_media();
        assert!((eff.scale_factor - 1.0).abs() < 1e-6);
        assert!((eff.diameter_m - eff.true_diameter_m).abs() < 1e-6);
        assert!((eff.density_kg_m3 - params.media.density_kg_m3).abs() < 1e-6);
        assert!(eff.ball_count < params.simulation.max_balls);
    }

    #[test]
    fn effective_media_coarse_grains_when_true_count_is_large() {
        // Default scenario: D = 1 m, d = 10 mm, J = 0.30 => several thousand true discs, above
        // max_balls (600), so coarse-graining must kick in.
        let params = Params::default();
        let n_true = params.true_ball_count();
        assert!(n_true > params.simulation.max_balls as f32);

        let eff = params.effective_media();
        assert!(eff.scale_factor > 1.0);
        assert!(eff.diameter_m > eff.true_diameter_m);
        // Unit-depth disc masses: density is left unchanged, area preservation alone preserves mass.
        assert!((eff.density_kg_m3 - params.media.density_kg_m3).abs() < 1e-6);
        assert!(eff.ball_count <= params.simulation.max_balls);
        // Should land close to the target, not just "under" it.
        assert!(eff.ball_count as f32 > 0.5 * params.simulation.max_balls as f32);
    }

    #[test]
    fn effective_media_preserves_total_footprint_area() {
        let params = Params::default();
        let r_true = params.media.ball_diameter_m * 0.5;
        let n_true = params.true_ball_count();
        let area_true = n_true * std::f32::consts::PI * r_true * r_true;

        let eff = params.effective_media();
        let r_eff = eff.diameter_m * 0.5;
        let area_sim = eff.ball_count as f32 * std::f32::consts::PI * r_eff * r_eff;

        let rel_err = (area_sim - area_true).abs() / area_true;
        assert!(rel_err < 0.01, "relative area error too large: {rel_err}");
    }

    #[test]
    fn effective_media_preserves_total_charge_mass() {
        // With unit-depth disc masses (`m = rho * pi * r^2`, see `crate::dem::ball_mass`),
        // preserving total footprint area (checked separately above) at unchanged density
        // preserves total charge mass automatically -- this test confirms that end-to-end.
        let params = Params::default();
        let r_true = params.media.ball_diameter_m * 0.5;
        let n_true = params.true_ball_count();
        let mass_true = n_true * params.media.density_kg_m3 * std::f32::consts::PI * r_true.powi(2);

        let eff = params.effective_media();
        let r_eff = eff.diameter_m * 0.5;
        let mass_sim =
            eff.ball_count as f32 * eff.density_kg_m3 * std::f32::consts::PI * r_eff.powi(2);

        let rel_err = (mass_sim - mass_true).abs() / mass_true;
        assert!(rel_err < 0.01, "relative mass error too large: {rel_err}");
    }

    #[test]
    fn true_ball_count_scales_with_packing_fraction() {
        let base = Params::default();
        let denser = Params {
            media: MediaParams {
                packing_fraction_2d: 0.907,
                ..base.media
            },
            ..base
        };
        let ratio = denser.true_ball_count() / base.true_ball_count();
        assert!((ratio - 0.907 / base.media.packing_fraction_2d).abs() < 1e-3);
    }
}
