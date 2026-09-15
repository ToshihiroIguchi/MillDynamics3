//! Drum geometry: signed-distance field (SDF), surface normal, and wall velocity.
//!
//! The drum is represented in its own rotating frame. To query geometry at a world-space point,
//! callers rotate the point into the drum frame by `-drum_angle` first (see [`Drum::sdf_world`]).
//!
//! The wall is a circle, optionally with `lifters.count` evenly-spaced trapezoidal bars protruding
//! inward (docs/PLAN.md ss3.1); `count == 0` (the default) is a perfectly smooth circle. The
//! lifter contribution and the circular wall are combined with `min` (nearest solid boundary
//! wins), and [`Drum::sdf_world`]'s numeric-gradient normal works unchanged for either case.

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

    /// Signed distance to the nearest lifter, positive = free space (i.e. already negated from
    /// the "positive = inside solid" convention [`lifter_cross_section_sdf`] uses internally, to
    /// match [`Drum::sdf`]'s overall convention). Each of the `lifters.count` bars is evaluated in
    /// its own (radial, tangential) frame about the drum centre and the nearest one wins.
    fn sdf_lifters(&self, p_local: Vec2) -> f32 {
        debug_assert!(self.lifters.count > 0);
        let n = self.lifters.count;
        let phase = self.lifters.phase_deg.to_radians();
        let mut max_inside = f32::NEG_INFINITY;
        for k in 0..n {
            let theta_i = phase + (k as f32) * std::f32::consts::TAU / (n as f32);
            let (sin_t, cos_t) = theta_i.sin_cos();
            // Radial/tangential coordinates of p_local about the drum centre, in lifter i's frame.
            let pr = p_local.x * cos_t + p_local.y * sin_t;
            let pt = -p_local.x * sin_t + p_local.y * cos_t;
            let inside = lifter_cross_section_sdf(
                pr,
                pt,
                self.radius_m,
                self.lifters.height_m,
                self.lifters.base_width_m * 0.5,
                self.lifters.top_width_m * 0.5,
            );
            max_inside = max_inside.max(inside);
        }
        -max_inside
    }

    /// Signed distance and outward-pointing normal at a **world-space** point, given the current
    /// drum rotation angle (radians). The normal points from the wall towards free space.
    ///
    /// A non-finite `p_world` returns `(0.0, Vec2::ZERO)` rather than computing anything: for an
    /// infinite input, `self.sdf` would return `-inf` at every probe, the central-difference
    /// gradient `(-inf) - (-inf)` would be `NaN`, the `length_squared() > 1e-12` fallback guard
    /// below is false for a NaN gradient so the `normalize_or_zero` branch runs, and
    /// `normalize_or_zero` on an infinite vector returns `Vec2::ZERO` -- so a caller's typical
    /// `p += -d * normal` correction would become `inf * 0 = NaN`, manufacturing a NaN position
    /// out of a merely-infinite input rather than just propagating one. Reporting "on the
    /// boundary, no usable normal" instead means every caller's `d < 0` (or similar) test
    /// declines to apply any correction, leaving cleanup to the caller's own state guard (e.g.
    /// [`crate::pbf::FluidParticles::step_coupled`]'s sanitize step).
    pub fn sdf_world(&self, p_world: Vec2, drum_angle: f32) -> (f32, Vec2) {
        if !p_world.is_finite() {
            return (0.0, Vec2::ZERO);
        }
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
    Vec2::new(v.x * unit.x - v.y * unit.y, v.x * unit.y + v.y * unit.x)
}

/// Signed distance to one lifter's trapezoidal cross-section, **positive when `(pr, pt)` is
/// inside the solid bar** (the opposite convention from [`Drum::sdf`]/[`Drum::sdf_lifters`], which
/// negate this before combining lifters with the circular wall -- see docs/PLAN.md ss3.1).
///
/// `pr`/`pt` are the point's radial/tangential coordinates about the drum centre, in the lifter's
/// own frame (its centreline along `pt = 0`). The bar spans `pr` in `[r - height, r]` (base at the
/// wall, top pointing inward) with tangential half-width `base_half` at the base and `top_half` at
/// the top -- a possibly-degenerate trapezoid (`top_half == base_half` is a rectangle, `top_half ==
/// 0` a triangle). The four vertices are wound counter-clockwise so each edge's left-hand normal
/// (`(-dy, dx)` of the edge direction) points inward; distance to each edge's *line* (not
/// segment) is exact on that edge and an approximation near corners, standard for a convex-SDF
/// combined via `min`/`max`.
fn lifter_cross_section_sdf(
    pr: f32,
    pt: f32,
    r: f32,
    height: f32,
    base_half: f32,
    top_half: f32,
) -> f32 {
    let r_top = r - height;
    let verts = [
        Vec2::new(r, -base_half),
        Vec2::new(r, base_half),
        Vec2::new(r_top, top_half),
        Vec2::new(r_top, -top_half),
    ];
    let p = Vec2::new(pr, pt);
    let mut d = f32::MAX;
    for i in 0..4 {
        let a = verts[i];
        let b = verts[(i + 1) % 4];
        let edge = b - a;
        let len = edge.length();
        if len < 1e-9 {
            continue;
        }
        let inward_normal = Vec2::new(-edge.y, edge.x) / len;
        d = d.min((p - a).dot(inward_normal));
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    fn smooth_drum(radius_m: f32, omega: f32) -> Drum {
        Drum::new(
            radius_m,
            omega,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        )
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

    #[test]
    fn sdf_world_never_manufactures_a_non_finite_result() {
        let drum = smooth_drum(0.3, 0.0);
        for p in [
            Vec2::new(f32::NAN, 0.0),
            Vec2::splat(f32::INFINITY),
            Vec2::new(0.0, f32::NEG_INFINITY),
        ] {
            let (d, n) = drum.sdf_world(p, 0.3);
            assert!(
                d.is_finite() && n.is_finite(),
                "sdf_world({p:?}) = ({d}, {n:?})"
            );
            // A caller applying the typical `p += -d * normal` correction must get a no-op,
            // not `inf * 0 = NaN`.
            assert!((-d * n).is_finite());
        }
    }

    fn lifter_drum(radius_m: f32, count: u32) -> Drum {
        Drum::new(
            radius_m,
            0.0,
            LiftersParams {
                count,
                height_m: 0.02,
                base_width_m: 0.03,
                top_width_m: 0.02,
                phase_deg: 0.0,
            },
        )
    }

    #[test]
    fn zero_lifters_matches_plain_circle_everywhere() {
        // A LiftersParams with count == 0 must produce exactly the plain-circle SDF regardless of
        // whatever geometry fields it carries (regression: adding real lifter geometry in M2 must
        // not change the M1 default, count == 0, behaviour).
        let with_lifters_off = Drum::new(
            0.4,
            0.0,
            LiftersParams {
                count: 0,
                height_m: 0.05,
                base_width_m: 0.05,
                top_width_m: 0.03,
                phase_deg: 12.0,
            },
        );
        let plain = smooth_drum(0.4, 0.0);
        for i in 0..64 {
            let angle = i as f32 * std::f32::consts::TAU / 64.0;
            let p = Vec2::new(angle.cos(), angle.sin()) * 0.4;
            assert!((with_lifters_off.sdf(p) - plain.sdf(p)).abs() < 1e-6);
        }
    }

    #[test]
    fn point_inside_a_lifter_is_solid() {
        // With 8 lifters and phase 0, lifter 0 is centred on +x at radius R. A point just inside
        // the wall along +x, which would be free space for a plain circle, must be solid (inside
        // the lifter) here.
        let drum = lifter_drum(0.3, 8);
        let p = Vec2::new(0.29, 0.0); // 10 mm inside the wall along the lifter's centreline
        assert!(
            drum.sdf(p) < 0.0,
            "expected solid (negative sdf) inside a lifter, got {}",
            drum.sdf(p)
        );

        // The plain circle at the same point would be free space, confirming the lifter (not the
        // wall itself) is what makes it solid.
        assert!(smooth_drum(0.3, 0.0).sdf(p) > 0.0);
    }

    #[test]
    fn point_between_lifters_is_free_space() {
        // Halfway between two of 8 evenly-spaced lifters (22.5 degrees off the nearest one's
        // centreline), a point close to the wall should still be free space: lifters are narrow
        // (30 mm base) compared to the gap between them at this radius.
        let drum = lifter_drum(0.3, 8);
        let angle = std::f32::consts::TAU / 16.0; // halfway between lifter 0 (0 deg) and lifter 1 (45 deg)
        let p = Vec2::new(angle.cos(), angle.sin()) * 0.29;
        assert!(
            drum.sdf(p) > 0.0,
            "expected free space between lifters, got {}",
            drum.sdf(p)
        );
    }

    #[test]
    fn lifter_normal_points_away_from_lifter_material() {
        let drum = lifter_drum(0.3, 8);
        // Just outside the lifter's tip (radially inward of its top face), the outward normal
        // should point further inward (away from the lifter, toward the drum centre / +x is the
        // lifter's centreline here, so "away" is -x).
        let tip_r = drum.radius_m - drum.lifters.height_m;
        let p = Vec2::new(tip_r - 0.005, 0.0);
        let (d, n) = drum.sdf_world(p, 0.0);
        assert!(
            d > 0.0,
            "point should be free space just inside the lifter's tip, got d={d}"
        );
        assert!(
            n.x < -0.5,
            "normal should point away from the lifter (toward -x), got n={n:?}"
        );
    }
}
