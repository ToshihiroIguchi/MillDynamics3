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

    /// Time of impact (a fraction `t` in `[0, 1]` of the sub-step) at which a disc of the given
    /// `radius`, whose centre moves in a straight **world-space** line from `x0_world` to
    /// `x1_world` while the drum itself rotates from `angle0` to `angle1`, first touches the
    /// wall/lifter solid. `None` if it never touches within the sub-step.
    ///
    /// Sphere-traces the drum's own static, rotation-free local-frame [`Drum::sdf`] along the
    /// ball's path re-expressed in the drum's *instantaneous* local frame at each traced `t` --
    /// i.e. this accounts for the drum's own rotation during the sub-step (a lifter sweeping into
    /// a resting ball is visible to this, unlike [`Drum::sdf_world`], which is a point query at
    /// one fixed angle and cannot see the wall move). Advancing `t` by `d / local_speed_bound` is
    /// a safe (never-overshooting) step because `sdf` is 1-Lipschitz in space and
    /// `local_speed_bound` is a rigorous upper bound on `d(local(t))/dt`: writing
    /// `local(t) = R(-angle(t)) * (x0 + t*(x1-x0))`, the product rule gives
    /// `d(local)/dt = R(-angle(t))*(x1-x0) - angle'(t) * perp(local(t))`, whose magnitude is at
    /// most `|x1-x0| + |angle1-angle0| * |x0 + t*(x1-x0)|` (rotation preserves length; `angle'(t)`
    /// is exactly the constant `angle1-angle0` since [`Drum`] assumes constant `omega` within a
    /// sub-step already) -- and `|x0 + t*(x1-x0)| <= |x0| + |x1-x0|` for any `t` in `[0, 1]`,
    /// with no assumption that the ball stays within the drum.
    ///
    /// Returns `None` (not `Some(0.0)`) when the ball already touches or overlaps the wall at
    /// `t = 0`: recovering an *existing* overlap -- including one this function itself produced
    /// by clamping the previous sub-step's advance to first contact -- is the discrete
    /// non-penetration solve's job (`dem.rs` step 3, with its own bounded-recovery policy), not
    /// this continuous-collision guard's. Treating "already touching" as a crossing here would
    /// report `t = 0` regardless of which way the ball is actually moving, freezing a ball that
    /// is separating from the wall just as readily as one still closing on it. The `t = 0` guard
    /// uses an exact `<= 0.0` test rather than the loop's own convergence tolerance below on
    /// purpose: a looser tolerance at `t = 0` would let the loop's first iteration re-trigger the
    /// exact freeze this guard exists to prevent, for any ball starting within that tolerance of
    /// the surface (routine immediately after a clamped contact, not an edge case).
    pub fn toi_swept(
        &self,
        x0_world: Vec2,
        x1_world: Vec2,
        angle0: f32,
        angle1: f32,
        radius: f32,
    ) -> Option<f32> {
        let local0 = rotate(x0_world, Vec2::from_angle(-angle0));
        if self.sdf(local0) - radius <= 0.0 {
            return None;
        }

        let path = x1_world - x0_world;
        let path_len = path.length();
        let dangle = angle1 - angle0;
        let pos_bound = x0_world.length() + path_len;
        let local_speed_bound = path_len + dangle.abs() * pos_bound;
        if local_speed_bound < 1e-9 {
            return None; // no relative motion, and the check above already ruled out t = 0 contact
        }

        // Exact `<= 0.0` throughout (not a loosened tolerance like `1e-6`): the guard above
        // already established `d(0) > 0` strictly, so the very first iteration below cannot
        // spuriously re-trigger at `t = 0` the way a looser per-iteration tolerance would.
        let mut t = 0.0f32;
        for _ in 0..32 {
            let angle_t = angle0 + dangle * t;
            let x_t = x0_world + path * t;
            let local = rotate(x_t, Vec2::from_angle(-angle_t));
            let d = self.sdf(local) - radius;
            if d <= 0.0 {
                return Some(t);
            }
            // `.max(1e-6)`: a minimum step-size floor so a `d` that is small-but-still-positive
            // cannot stall the loop just short of the 32-iteration cap -- unrelated to (and not a
            // reintroduction of) the stopping tolerance removed above.
            let step_t = (d / local_speed_bound).max(1e-6);
            t += step_t;
            if t >= 1.0 {
                return None;
            }
        }
        // Ran out of iterations without a clean root: treat the last sample as contact rather
        // than silently reporting "no contact" -- a false negative is the dangerous direction for
        // a tunnelling guard, a false positive only costs one extra (harmless) depenetration pass
        // in the discrete solve that follows.
        Some(t.min(1.0))
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
/// 0` a triangle). The four vertices are wound counter-clockwise.
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
    convex_polygon_sdf(Vec2::new(pr, pt), &verts)
}

/// Exact signed distance from `p` to a convex polygon with counter-clockwise-wound vertices,
/// positive when `p` is inside.
///
/// A naive convex-SDF -- `min` over each edge's distance to that edge's *infinite* supporting
/// line -- is exact for interior points (the nearest boundary point of an interior point is
/// always a perpendicular foot on some edge) but only an *approximation* for exterior points
/// whose nearest feature is a vertex, not an edge: the unclamped line distance for whichever edge
/// happens to be nearest always comes out smaller in magnitude than the true Euclidean distance
/// to that vertex (the line distance is one leg of the right triangle whose hypotenuse is the
/// true distance). That under-estimate made the lifter's convex corners register contact
/// *earlier* than geometrically correct -- not later -- which is also what inflated
/// `max_ball_wall_overlap_fraction` under `lifters.count > 0` (docs/PHYSICS.md §9): the SDF itself
/// was reporting less clearance than the ball truly had.
///
/// The fix: clamp each edge's closest-point projection to the segment itself (so a point beyond
/// an edge's endpoint measures to that endpoint, not the infinite line through it), and determine
/// the sign separately via a point-in-convex-polygon half-plane test. This is exact everywhere,
/// interior or exterior, at the same O(vertex count) cost as the line-based version.
fn convex_polygon_sdf(p: Vec2, verts: &[Vec2]) -> f32 {
    let n = verts.len();
    let mut min_dist_sq = f32::MAX;
    let mut inside = true;
    for i in 0..n {
        let a = verts[i];
        let b = verts[(i + 1) % n];
        let edge = b - a;
        let to_p = p - a;
        let len_sq = edge.length_squared();
        if len_sq < 1e-12 {
            continue;
        }
        let t = (to_p.dot(edge) / len_sq).clamp(0.0, 1.0);
        let closest = a + edge * t;
        min_dist_sq = min_dist_sq.min((p - closest).length_squared());
        // CCW winding => the inward half-plane is where the point is left of the edge.
        let cross = edge.x * to_p.y - edge.y * to_p.x;
        if cross < 0.0 {
            inside = false;
        }
    }
    let dist = min_dist_sq.sqrt();
    if inside {
        dist
    } else {
        -dist
    }
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
    fn toi_swept_finds_the_wall_crossing_of_a_fast_straight_shot() {
        // A ball fired straight at the (stationary, no-lifter) wall from well inside the drum,
        // moving far enough in one sub-step that it would tunnel through the wall entirely
        // without a swept test (docs/PHYSICS.md §9's former "no actual swept/continuous-
        // collision correction" gap). The true crossing (ball surface first touches the wall) is
        // where `|x(t)| = radius_m - ball_radius`.
        let drum = smooth_drum(0.3, 0.0);
        let ball_radius = 0.01;
        let x0 = Vec2::new(0.0, 0.0);
        let x1 = Vec2::new(0.5, 0.0); // would end up far outside the drum without clamping
        let true_t = (drum.radius_m - ball_radius) / (x1.x - x0.x);
        let t = drum
            .toi_swept(x0, x1, 0.0, 0.0, ball_radius)
            .expect("must find a wall crossing");
        assert!((t - true_t).abs() < 1e-3, "expected toi={true_t}, got {t}");
    }

    #[test]
    fn toi_swept_does_not_refreeze_a_ball_already_touching_the_wall() {
        // A ball resting exactly at the wall at t=0 must not be reported as a "crossing",
        // regardless of which way it is about to move: recovering an existing overlap is the
        // discrete solver's job, not this guard's. This is the exact regression an earlier
        // version of this function had -- treating "already touching" as `Some(0.0)` froze a
        // ball at its previous position every sub-step it stayed near the wall (including a
        // separating ball, not just an approaching one), because clamping to `t = 0` means "do
        // not move at all this sub-step" and CCD re-ran that same freeze on the very next
        // sub-step too, since the ball was still touching there as well.
        let drum = smooth_drum(0.3, 0.0);
        let ball_radius = 0.01;
        let x0 = Vec2::new(drum.radius_m - ball_radius, 0.0); // exactly touching
        let x1 = x0 + Vec2::new(0.001, 0.0); // moving slightly further into the wall
        let t = drum.toi_swept(x0, x1, 0.0, 0.0, ball_radius);
        assert!(
            t.is_none(),
            "expected no CCD crossing for an already-touching ball, got {t:?}"
        );
    }

    #[test]
    fn toi_swept_is_none_for_a_path_that_stays_inside() {
        let drum = smooth_drum(0.3, 0.0);
        let t = drum.toi_swept(Vec2::new(0.0, 0.0), Vec2::new(0.01, 0.0), 0.0, 0.0, 0.01);
        assert!(t.is_none(), "expected no wall crossing, got {t:?}");
    }

    #[test]
    fn toi_swept_catches_a_lifter_sweeping_into_a_resting_ball() {
        // The ball does not move at all this sub-step (x0 == x1); only the drum rotates a lifter
        // into it. `p = (0.29, 0.0)` is the same point `point_inside_a_lifter_is_solid` uses --
        // solid when queried at drum_angle = 0. A check using only the entry angle (as the
        // existing non-swept `drum.sdf_world(x, drum_angle)` narrow-phase does, evaluated once at
        // the sub-step's starting angle) would see the lifter far away at `angle0 = -0.3` and
        // report free space; the swept test must still catch the contact that exists by `t = 1`.
        let drum = lifter_drum(0.3, 1);
        let ball_radius = 0.002;
        let p = Vec2::new(0.29, 0.0);
        let angle0 = -0.3;
        let angle1 = 0.0;

        let (d_entry, _) = drum.sdf_world(p, angle0);
        assert!(
            d_entry - ball_radius > 0.0,
            "test setup: must be free space at the entry angle, got clearance {}",
            d_entry - ball_radius
        );

        let t = drum.toi_swept(p, p, angle0, angle1, ball_radius);
        assert!(
            t.is_some(),
            "a lifter sweeping onto a resting ball's position must register a crossing"
        );
    }

    #[test]
    fn lifter_corner_sdf_is_exact_not_a_vertex_underestimate() {
        // Regression for the exterior-corner SDF bug (second external review, finding "⑦"): the
        // old "min over each edge's infinite line" formula under-estimated the true distance to
        // a convex vertex. With this drum's default lifter geometry (height=0.02, base=0.03,
        // top=0.02, so r_top = 0.28 and the tip's near-+y vertex is at (0.28, 0.01)), the old
        // formula read 0.0121 m of clearance at the point checked below instead of the true
        // 0.0141 m (Euclidean distance to that vertex) -- a ball there was treated as already
        // closer to the lifter than it geometrically was. `lifter_drum` uses `phase_deg = 0.0`
        // and this drum has a single lifter, so at `drum_angle = 0.0` the lifter's own
        // (radial, tangential) frame coincides exactly with world (x, y).
        let drum = lifter_drum(0.3, 1);
        let p = Vec2::new(0.27, 0.02);
        let vertex = Vec2::new(0.28, 0.01);
        let true_dist = (p - vertex).length();
        let (d, _) = drum.sdf_world(p, 0.0);
        assert!(
            (d - true_dist).abs() < 1e-4,
            "expected exact corner distance {true_dist}, got {d} (old formula read ~0.0121)"
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
