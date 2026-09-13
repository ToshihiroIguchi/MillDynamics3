//! Drum geometry: signed-distance field (SDF), surface normal, and wall velocity.
//!
//! The drum is represented in its own rotating frame. To query geometry at a world-space point,
//! callers rotate the point into the drum frame by `-drum_angle` first (see [`Drum::sdf_world`]).
//!
//! Through milestone M1 the wall is a plain circle (lifters are added in M2, see docs/PLAN.md
//! ss3.1). `Drum::sdf` already takes the full [`LiftersParams`] so the M2 change is additive: it
//! only has to replace the body of `sdf_lifters`, not the calling convention.

use glam::Vec2;

use crate::params::LiftersParams;

/// Drum wall geometry and kinematics, derived once per parameter set (or per frame, since it is
/// cheap: no allocation, a handful of floats).
#[derive(Debug, Clone, Copy)]
pub struct Drum {
    pub radius_m: f32,
    pub omega: f32,
    pub lifters: LiftersParams,
}

impl Drum {
    pub fn new(radius_m: f32, omega: f32, lifters: LiftersParams) -> Self {
        Self {
            radius_m,
            omega,
            lifters,
        }
    }

    /// Signed distance to the nearest solid boundary, evaluated in the **drum-local frame**
    /// (i.e. `p` must already be rotated by `-drum_angle` relative to world space).
    ///
    /// Positive = inside free space (distance to the nearest wall); negative = inside solid.
    pub fn sdf(&self, p_local: Vec2) -> f32 {
        let dist_to_wall = self.radius_m - p_local.length();
        if self.lifters.count == 0 {
            dist_to_wall
        } else {
            dist_to_wall.min(self.sdf_lifters(p_local))
        }
    }

    /// Placeholder for M2: currently returns `f32::INFINITY` (no effect) since `count == 0` by
    /// default and this path is otherwise unreachable until lifter geometry lands.
    fn sdf_lifters(&self, _p_local: Vec2) -> f32 {
        debug_assert!(self.lifters.count > 0);
        f32::INFINITY
    }

    /// Signed distance and outward-pointing normal at a **world-space** point, given the current
    /// drum rotation angle (radians). The normal points from the wall towards free space.
    pub fn sdf_world(&self, p_world: Vec2, drum_angle: f32) -> (f32, Vec2) {
        let to_local = Vec2::from_angle(-drum_angle);
        let p_local = rotate(p_world, to_local);

        const EPS: f32 = 1e-4;
        let d = self.sdf(p_local);
        let dx = self.sdf(p_local + Vec2::new(EPS, 0.0)) - self.sdf(p_local - Vec2::new(EPS, 0.0));
        let dy = self.sdf(p_local + Vec2::new(0.0, EPS)) - self.sdf(p_local - Vec2::new(0.0, EPS));
        let grad_local = Vec2::new(dx, dy) / (2.0 * EPS);
        let normal_local = if grad_local.length_squared() > 1e-12 {
            grad_local.normalize()
        } else {
            -p_local.normalize_or_zero()
        };

        let to_world = Vec2::from_angle(drum_angle);
        let normal_world = rotate(normal_local, to_world);
        (d, normal_world)
    }

    /// Rigid-body velocity of the wall material point that is currently at `p_world`
    /// (`v = omega x r`, i.e. `omega * perp(p)` in 2D).
    pub fn wall_velocity(&self, p_world: Vec2) -> Vec2 {
        self.omega * Vec2::new(-p_world.y, p_world.x)
    }
}

/// Rotates a 2D vector by a unit vector representing `(cos, sin)` of the rotation angle.
fn rotate(v: Vec2, unit: Vec2) -> Vec2 {
    Vec2::new(
        v.x * unit.x - v.y * unit.y,
        v.x * unit.y + v.y * unit.x,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    fn smooth_drum(radius_m: f32, omega: f32) -> Drum {
        Drum::new(radius_m, omega, LiftersParams { count: 0, ..LiftersParams::default() })
    }

    #[test]
    fn sdf_zero_at_wall_positive_inside() {
        let drum = smooth_drum(0.3, 0.0);
        assert!((drum.sdf(Vec2::new(0.3, 0.0)) - 0.0).abs() < 1e-6);
        assert!(drum.sdf(Vec2::ZERO) > 0.0);
        assert!(drum.sdf(Vec2::new(0.35, 0.0)) < 0.0);
    }

    #[test]
    fn normal_points_inward_toward_center() {
        let drum = smooth_drum(0.3, 0.0);
        let (_, n) = drum.sdf_world(Vec2::new(0.3, 0.0), 0.0);
        // At the +x wall point, free space is toward -x, so the normal should point that way.
        assert!(n.x < -0.9, "normal={n:?}");
    }

    #[test]
    fn sdf_world_is_invariant_under_rotation() {
        let drum = smooth_drum(0.3, 0.0);
        let p = Vec2::new(0.25, 0.1);
        let (d0, _) = drum.sdf_world(p, 0.0);
        // Rotating both the query point and the drum angle by the same amount must not change
        // the signed distance (only the reported world-space normal direction should rotate).
        let angle = 0.7;
        let p_rotated = rotate(p, Vec2::from_angle(angle));
        let (d1, _) = drum.sdf_world(p_rotated, angle);
        assert!((d0 - d1).abs() < 1e-4);
    }

    #[test]
    fn wall_velocity_matches_omega_cross_r() {
        let drum = smooth_drum(0.3, 2.0);
        let p = Vec2::new(0.3, 0.0);
        let v = drum.wall_velocity(p);
        // omega x r for r=(R,0), omega=2 (ccw, +z) => v = (0, omega*R)
        assert!((v - Vec2::new(0.0, 0.6)).length() < 1e-5);
    }

    #[test]
    fn full_turn_angle_returns_same_normal() {
        let drum = smooth_drum(0.3, 0.0);
        let p = Vec2::new(0.3, 0.0);
        let (_, n0) = drum.sdf_world(p, 0.0);
        let (_, n1) = drum.sdf_world(p, 2.0 * PI);
        assert!((n0 - n1).length() < 1e-3);
    }
}
