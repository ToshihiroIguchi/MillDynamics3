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
            speed_mode: SpeedMode::PercentCritical,
            speed_value: 70.0,
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
    /// Fraction of the drum's cross-sectional area occupied by media, including voids
    /// (i.e. the conventional "J" ball-filling fraction).
    pub fill_fraction: f32,
    /// Media (steel) density (kg/m^3).
    pub density_kg_m3: f32,
    pub restitution_ball_ball: f32,
    pub restitution_ball_wall: f32,
    pub friction_ball_ball: f32,
    pub friction_ball_wall: f32,
    /// Rolling-resistance coefficient (dimensionless torque coefficient).
    pub rolling_friction: f32,
}

impl Default for MediaParams {
    fn default() -> Self {
        Self {
            ball_diameter_m: 0.002,
            fill_fraction: 0.30,
            density_kg_m3: 7800.0,
            restitution_ball_ball: 0.5,
            restitution_ball_wall: 0.3,
            friction_ball_ball: 0.4,
            friction_ball_wall: 0.5,
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
        if !(self.density_kg_m3 > 0.0) {
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
    pub dye_pattern: DyePattern,
}

impl Default for SlurryParams {
    fn default() -> Self {
        Self {
            enabled: true,
            fill_fraction: 0.15,
            density_kg_m3: 1800.0,
            viscosity_pa_s: 0.5,
            rheology: Rheology::Newtonian,
            yield_stress_pa: 0.0,
            wall_no_slip: 1.0,
            ball_no_slip: 1.0,
            dye_pattern: DyePattern::LeftRight,
        }
    }
}

impl SlurryParams {
    pub fn validate(&self) -> Result<(), String> {
        if !(0.0..=0.9).contains(&self.fill_fraction) {
            return Err("slurry.fill_fraction must be in [0, 0.9]".into());
        }
        if !(self.density_kg_m3 > 0.0) {
            return Err("slurry.density_kg_m3 must be > 0".into());
        }
        if !(self.viscosity_pa_s >= 0.0 && self.viscosity_pa_s.is_finite()) {
            return Err("slurry.viscosity_pa_s must be >= 0".into());
        }
        if !(self.yield_stress_pa >= 0.0 && self.yield_stress_pa.is_finite()) {
            return Err("slurry.yield_stress_pa must be >= 0".into());
        }
        for (name, b) in [
            ("wall_no_slip", self.wall_no_slip),
            ("ball_no_slip", self.ball_no_slip),
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
    /// PBF density-constraint solver iterations per sub-step.
    pub pbf_iterations: u32,
    /// Wall-clock-to-simulation-time multiplier requested by the user (achieved rate may be
    /// lower and is reported back to the UI, see docs/PLAN.md ss4.1).
    pub time_scale: f32,
    /// Maximum wall-clock milliseconds the worker may spend simulating per animation frame.
    pub frame_budget_ms: f32,
    pub seed: u64,
}

impl Default for SimulationParams {
    fn default() -> Self {
        Self {
            resolution: 40,
            substeps: 4,
            pbf_iterations: 3,
            time_scale: 1.0,
            frame_budget_ms: 12.0,
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
        if !(self.time_scale > 0.0 && self.time_scale.is_finite()) {
            return Err("simulation.time_scale must be > 0".into());
        }
        if !(self.frame_budget_ms > 0.0 && self.frame_budget_ms.is_finite()) {
            return Err("simulation.frame_budget_ms must be > 0".into());
        }
        Ok(())
    }
}

/// The complete set of parameters needed to (re)construct a [`crate::Simulation`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Params {
    pub mill: MillParams,
    pub media: MediaParams,
    pub slurry: SlurryParams,
    pub lifters: LiftersParams,
    pub simulation: SimulationParams,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            mill: MillParams::default(),
            media: MediaParams::default(),
            slurry: SlurryParams::default(),
            lifters: LiftersParams::default(),
            simulation: SimulationParams::default(),
        }
    }
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
}
