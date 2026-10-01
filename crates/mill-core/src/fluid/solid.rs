//! Smooth solid boundaries (drum wall, lifters, balls) as precomputed kernel integrals.
//!
//! A solid region `S` acts on a fluid particle at `x` through two kernel integrals, evaluated once
//! into 1-D tables (Koschier & Bender's density-map idea, analytic for discs):
//!
//! - `Gamma(x) = integral_S W(|x - x'|) dA'`: the solid volume fraction the kernel sees, so the
//!   boundary's share of the density is `rho0 * Gamma`, its pressure acceleration is `-kappa rho0
//!   grad Gamma` (normal to the surface, hence the reaction on a ball passes through its centre:
//!   no pressure torque, exact buoyancy), and its `d rho / dt` term is `rho0 (v - v_S) . grad Gamma`.
//! - `T(x) = integral_S K(|x - x'|) e e^T dA'` with `K = |W'(r)| r / (r^2 + eta^2)` and `e` the pair
//!   unit vector: the viscous coupling tensor of the same pairwise-central form as fluid-fluid
//!   pairs. `T = A n n^T + B t t^T` with `n` the direction to the solid; the viscous force per unit
//!   mass is `-(8 mu / rho_i) T (v_i - v_S)`.
//!
//! These replace sums over discrete boundary particles, so the boundary has no surface roughness,
//! no gap to the first fluid row and no staggering: the solid half-space couples like fluid and
//! the no-slip plane is the interface itself. For a rigid body `v_S` is the body's velocity
//! *at the fluid particle* - identical to integrating over the body because central forces see
//! `(omega x delta) . delta = 0`.
//!
//! Tables are built for a disc of radius `a` (a ball; the drum wall is the complement of a disc),
//! as functions of the distance `r` to its centre, and for a flat half-space (lifters) as a
//! function of the signed distance `d` into the fluid.

use glam::Vec2;

use super::kernel::Kernel;
use crate::geometry::Drum;
use crate::params::LiftersParams;

/// `eta^2 = ETA_FACTOR * H^2` regularises `1 / r^2` for coincident particles (viscosity pairs).
pub(super) const ETA_FACTOR: f32 = 0.01;
/// Radial quadrature samples over the kernel support.
const RHO_SAMPLES: usize = 128;
/// Table entries over `[s_min, s_max]`.
const TABLE_ENTRIES: usize = 256;
/// `owner` of the drum wall (otherwise a ball index).
pub(super) const WALL: u32 = u32::MAX;

/// The solid's effect on one fluid particle.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SolidSample {
    pub owner: u32,
    pub gamma: f32,
    /// `grad Gamma` (points toward the solid).
    pub grad: Vec2,
    /// Viscous coupling tensor `T = b I + (a - b) n n^T`, `n` the unit direction to the solid.
    pub a: f32,
    pub b: f32,
    pub n: Vec2,
    /// Rigid velocity of the solid at the fluid particle.
    pub vel: Vec2,
    /// Offset of the fluid particle from the body's reference point (centre) for torque.
    pub lever: Vec2,
}

impl SolidSample {
    /// `T * v`.
    #[inline]
    pub fn apply_t(&self, v: Vec2) -> Vec2 {
        self.b * v + (self.a - self.b) * self.n.dot(v) * self.n
    }
}

/// `(Gamma, A, B)` of a solid disc of radius `a` seen from distance `r` of its centre.
fn disc_values(a: f64, r: f64, kernel: &Kernel, eta2: f64) -> (f64, f64, f64) {
    let h = kernel.support as f64;
    let (mut gamma, mut av, mut bv) = (0.0f64, 0.0f64, 0.0f64);
    // The same quadrature over a full circle normalises Gamma so that deep inside it is exactly 1.
    let mut gamma_full = 0.0f64;
    let drho = h / RHO_SAMPLES as f64;
    for k in 0..RHO_SAMPLES {
        let rho = (k as f64 + 0.5) * drho;
        // Half-angle of the part of the circle of radius `rho` around x that lies inside the
        // disc. The clamp makes the two limits fall out: pi when the circle is entirely inside
        // (rho < a - r), 0 when it misses the disc (rho < r - a or rho > a + r).
        let c = if r > 1e-12 {
            ((rho * rho + r * r - a * a) / (2.0 * rho * r)).clamp(-1.0, 1.0)
        } else if rho < a {
            -1.0
        } else {
            1.0
        };
        let phi = c.acos();
        let w = kernel.w(rho as f32) as f64;
        let k_pair = (kernel.dw_dr(rho as f32) as f64).abs() * rho / (rho * rho + eta2);
        gamma += w * rho * 2.0 * phi * drho;
        gamma_full += w * rho * 2.0 * std::f64::consts::PI * drho;
        let sc = phi.sin() * c;
        av += k_pair * rho * (phi + sc) * drho;
        bv += k_pair * rho * (phi - sc) * drho;
    }
    (gamma / gamma_full, av, bv)
}

/// Piecewise-linear table of `(Gamma, A, B)` against an abscissa `s`.
#[derive(Clone)]
struct Table {
    s_min: f32,
    ds: f32,
    gamma: Vec<f32>,
    a: Vec<f32>,
    b: Vec<f32>,
    /// Values for `s < s_min` and `s > s_max`.
    below: [f32; 3],
    above: [f32; 3],
}

impl Table {
    /// `(Gamma, dGamma/ds, A, B)`.
    #[inline]
    fn eval(&self, s: f32) -> (f32, f32, f32, f32) {
        let n = self.gamma.len();
        let t = (s - self.s_min) / self.ds;
        if t < 0.0 {
            return (self.below[0], 0.0, self.below[1], self.below[2]);
        }
        if t >= (n - 1) as f32 {
            return (self.above[0], 0.0, self.above[1], self.above[2]);
        }
        let i = t as usize;
        let f = t - i as f32;
        let lerp = |v: &[f32]| v[i] + f * (v[i + 1] - v[i]);
        (
            lerp(&self.gamma),
            (self.gamma[i + 1] - self.gamma[i]) / self.ds,
            lerp(&self.a),
            lerp(&self.b),
        )
    }
}

fn full_pair_integral(kernel: &Kernel, eta2: f64) -> f64 {
    let h = kernel.support as f64;
    let drho = h / RHO_SAMPLES as f64;
    let mut sum = 0.0;
    for k in 0..RHO_SAMPLES {
        let rho = (k as f64 + 0.5) * drho;
        let k_pair = (kernel.dw_dr(rho as f32) as f64).abs() * rho / (rho * rho + eta2);
        sum += k_pair * rho * drho;
    }
    std::f64::consts::PI * sum
}

/// Table of a disc of radius `a` against `r` (convex solid, e.g. a ball).
fn disc_table(a: f32, kernel: &Kernel) -> Table {
    let h = kernel.support;
    let eta2 = (ETA_FACTOR * h * h) as f64;
    let s_min = (a - h).max(0.0);
    let s_max = a + h;
    let ds = (s_max - s_min) / (TABLE_ENTRIES - 1) as f32;
    let (mut gamma, mut av, mut bv) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..TABLE_ENTRIES {
        let r = (s_min + i as f32 * ds) as f64;
        let (g, a_, b_) = disc_values(a as f64, r, kernel, eta2);
        gamma.push(g as f32);
        av.push(a_ as f32);
        bv.push(b_ as f32);
    }
    let full = full_pair_integral(kernel, eta2) as f32;
    Table {
        s_min,
        ds,
        gamma,
        a: av,
        b: bv,
        below: [1.0, full, full],
        above: [0.0, 0.0, 0.0],
    }
}

/// Table of the complement of a disc of radius `radius` (the drum wall: solid outside) against `r`.
fn wall_circle_table(radius: f32, kernel: &Kernel) -> Table {
    let disc = disc_table(radius, kernel);
    let full = disc.below[1];
    Table {
        s_min: disc.s_min,
        ds: disc.ds,
        gamma: disc.gamma.iter().map(|g| 1.0 - g).collect(),
        a: disc.a.iter().map(|a| full - a).collect(),
        b: disc.b.iter().map(|b| full - b).collect(),
        below: [0.0, 0.0, 0.0],
        above: [1.0, full, full],
    }
}

/// Table of a flat solid half-space against the signed distance `d` into the fluid (`d < 0` is
/// solid), built from a disc so large its curvature is negligible.
fn flat_table(kernel: &Kernel) -> Table {
    let h = kernel.support;
    let eta2 = (ETA_FACTOR * h * h) as f64;
    let big = 1.0e3 * h as f64;
    let ds = 2.0 * h / (TABLE_ENTRIES - 1) as f32;
    let (mut gamma, mut av, mut bv) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..TABLE_ENTRIES {
        let d = (-h + i as f32 * ds) as f64;
        let (g, a_, b_) = disc_values(big, big + d, kernel, eta2);
        gamma.push(g as f32);
        av.push(a_ as f32);
        bv.push(b_ as f32);
    }
    let full = full_pair_integral(kernel, eta2) as f32;
    // `d` increases into the fluid: deep solid is `s < s_min`.
    Table {
        s_min: -h,
        ds,
        gamma,
        a: av,
        b: bv,
        below: [1.0, full, full],
        above: [0.0, 0.0, 0.0],
    }
}

/// The drum wall (and lifters) as a smooth solid.
pub(super) struct WallModel {
    radius: f32,
    lifters: LiftersParams,
    circle: Table,
    flat: Table,
    support: f32,
}

impl WallModel {
    pub fn build(drum: &Drum, kernel: &Kernel) -> Self {
        Self {
            radius: drum.radius_m,
            lifters: drum.lifters,
            circle: wall_circle_table(drum.radius_m, kernel),
            flat: flat_table(kernel),
            support: kernel.support,
        }
    }

    pub fn matches(&self, drum: &Drum) -> bool {
        self.radius == drum.radius_m && self.lifters == drum.lifters
    }

    /// The wall's effect at `x` (`None` if farther than the kernel support from it).
    pub fn sample(&self, x: Vec2, drum: &Drum, angle: f32) -> Option<SolidSample> {
        let vel = drum.wall_velocity(x);
        if self.lifters.count == 0 {
            let r = x.length();
            if r < self.radius - self.support {
                return None;
            }
            let dir = if r > 1e-9 { x / r } else { Vec2::X };
            let (gamma, dg, a, b) = self.circle.eval(r);
            Some(SolidSample {
                owner: WALL,
                gamma,
                grad: dg * dir,
                a,
                b,
                n: dir,
                vel,
                lever: x,
            })
        } else {
            let (d, n_free) = drum.sdf_world(x, angle);
            if d > self.support {
                return None;
            }
            let (gamma, dg, a, b) = self.flat.eval(d);
            Some(SolidSample {
                owner: WALL,
                gamma,
                grad: dg * n_free,
                a,
                b,
                n: -n_free,
                vel,
                lever: x,
            })
        }
    }
}

/// Balls (all of one radius) as smooth discs.
pub(super) struct BallModel {
    pub radius: f32,
    table: Table,
    support: f32,
}

impl BallModel {
    pub fn empty() -> Self {
        Self {
            radius: 0.0,
            table: Table {
                s_min: 0.0,
                ds: 1.0,
                gamma: vec![0.0, 0.0],
                a: vec![0.0, 0.0],
                b: vec![0.0, 0.0],
                below: [0.0; 3],
                above: [0.0; 3],
            },
            support: 0.0,
        }
    }

    pub fn build(radius: f32, kernel: &Kernel) -> Self {
        Self {
            radius,
            table: disc_table(radius, kernel),
            support: kernel.support,
        }
    }

    /// Reach (from a ball centre) beyond which a ball does not affect a fluid particle.
    pub fn reach(&self) -> f32 {
        self.radius + self.support
    }

    pub fn sample(
        &self,
        owner: u32,
        x: Vec2,
        centre: Vec2,
        v: Vec2,
        omega: f32,
    ) -> Option<SolidSample> {
        let delta = x - centre;
        let r = delta.length();
        if r >= self.reach() {
            return None;
        }
        let dir = if r > 1e-9 { delta / r } else { Vec2::X };
        let (gamma, dg, a, b) = self.table.eval(r);
        Some(SolidSample {
            owner,
            gamma,
            grad: dg * dir,
            a,
            b,
            n: -dir,
            vel: v + omega * Vec2::new(-delta.y, delta.x),
            lever: delta,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kernel() -> Kernel {
        Kernel::new(2.0e-3)
    }

    #[test]
    fn flat_gamma_runs_from_one_inside_to_zero_outside_through_one_half() {
        let k = kernel();
        let t = flat_table(&k);
        let h = k.support;
        assert!((t.eval(-h * 0.999).0 - 1.0).abs() < 2e-3);
        assert!(t.eval(h * 0.999).0 < 2e-3);
        assert!((t.eval(0.0).0 - 0.5).abs() < 5e-3, "{}", t.eval(0.0).0);
        let mut prev = 2.0;
        for i in 0..200 {
            let d = -h + 2.0 * h * i as f32 / 199.0;
            let g = t.eval(d).0;
            assert!(g <= prev + 1e-6, "not monotone at d = {d}");
            prev = g;
        }
    }

    #[test]
    fn deep_inside_the_coupling_tensor_is_isotropic_with_the_full_value() {
        let k = kernel();
        let eta2 = (ETA_FACTOR * k.support * k.support) as f64;
        let full = full_pair_integral(&k, eta2) as f32;
        let t = flat_table(&k);
        let (_, _, a, b) = t.eval(-k.support * 0.999);
        assert!((a - full).abs() < 1e-2 * full && (b - full).abs() < 1e-2 * full);
        // at the interface the half-space gives about half of the full value
        let (_, _, a0, b0) = t.eval(0.0);
        assert!((a0 + b0 - full).abs() < 0.03 * full, "{} {} {full}", a0, b0);
        assert!(a0 > 0.0 && b0 > 0.0);
    }

    #[test]
    fn small_disc_table_matches_the_direct_kernel_integral_at_its_centre() {
        let k = kernel();
        let a = 0.5 * k.support;
        let t = disc_table(a, &k);
        // Gamma(r = 0) = integral_0^a W(rho) 2 pi rho drho
        let n = 20000;
        let mut direct = 0.0f64;
        for i in 0..n {
            let rho = (i as f64 + 0.5) * a as f64 / n as f64;
            direct +=
                k.w(rho as f32) as f64 * 2.0 * std::f64::consts::PI * rho * a as f64 / n as f64;
        }
        let g0 = t.eval(0.0).0 as f64;
        assert!((g0 - direct).abs() < 0.01 * direct, "{g0} vs {direct}");
    }

    #[test]
    fn gradient_matches_a_finite_difference_of_gamma() {
        let k = kernel();
        let wall = WallModel::build(
            &Drum::new(
                0.0315,
                0.0,
                LiftersParams {
                    count: 0,
                    ..LiftersParams::default()
                },
            ),
            &k,
        );
        let drum = Drum::new(
            0.0315,
            0.0,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        );
        let x = Vec2::new(0.0315 - 0.7e-3, 0.0);
        let g = |p: Vec2| wall.sample(p, &drum, 0.0).unwrap().gamma;
        let eps = 2.0e-5;
        let fd = (g(x + Vec2::new(eps, 0.0)) - g(x - Vec2::new(eps, 0.0))) / (2.0 * eps);
        let s = wall.sample(x, &drum, 0.0).unwrap();
        assert!(s.grad.x > 0.0, "grad must point toward the wall");
        assert!(
            (fd - s.grad.x).abs() < 0.05 * fd.abs(),
            "{fd} vs {}",
            s.grad.x
        );
    }

    #[test]
    fn a_hex_fluid_row_at_half_a_row_height_sees_rest_density() {
        // Fluid lattice (hex, spacing dx) in the half-space y > 0 with its first row at row_h/2
        // from a flat wall at y = 0 (solid below): fluid kernel sum + rho0 Gamma = rho0.
        let dx = 1.0e-3f32;
        let row_h = dx * 3.0f32.sqrt() * 0.5;
        let k = Kernel::new(2.0 * dx);
        let rho0 = 1800.0f32;
        let mut rho_unit = 0.0;
        for row in -4i32..=4 {
            let off = if row.rem_euclid(2) == 0 { 0.0 } else { 0.5 };
            for col in -4i32..=4 {
                rho_unit += k.w(Vec2::new((col as f32 + off) * dx, row as f32 * row_h).length());
            }
        }
        let m = rho0 / rho_unit;
        let flat = flat_table(&k);
        let y0 = 0.5 * row_h;
        let mut rho = rho0 * flat.eval(y0).0;
        for row in 0..6i32 {
            let off = if row.rem_euclid(2) == 0 { 0.0 } else { 0.5 };
            for col in -6i32..=6 {
                let p = Vec2::new((col as f32 + off) * dx, row as f32 * row_h);
                rho += m * k.w(p.length());
            }
        }
        assert!(
            (rho / rho0 - 1.0).abs() < 0.03,
            "first-row density ratio {}",
            rho / rho0
        );
    }
}
