//! Static solid boundary particles for the drum wall (Akinci et al. 2012), stored in the drum's
//! own rotating frame.
//!
//! The smooth wall is `RINGS` concentric rings outside the wall, ring `l` at `R + (l + 1/2) row_h`
//! (`row_h = dx sqrt(3)/2`, the fluid lattice's row spacing) so the wall surface (`sdf = 0`) lies
//! halfway between the outermost fluid row and the first solid row. Each ring has an integer
//! number of particles at ~`dx` arc spacing, alternate rings staggered by half a spacing (hex
//! pattern). Lifters, if any, are filled with a plain hex lattice of solid particles over their
//! cross-section (sites at least `row_h / 2` inside the surface). Every boundary particle carries
//! its own geometric volume `psi_b` (as mass: `psi_b = m * cell_area_b / (dx * row_h)`).

use glam::Vec2;

use crate::geometry::Drum;
use crate::params::LiftersParams;

/// Number of wall rings; 3 rows span `2.6 dx`, beyond the kernel support `2 dx`.
const RINGS: usize = 3;

pub struct Boundary {
    /// Positions in the drum-local frame.
    pub local: Vec<Vec2>,
    /// Per-particle boundary volume (mass-like, kg per metre of depth).
    pub psi: Vec<f32>,
    radius_m: f32,
    lifters: LiftersParams,
}

impl Boundary {
    /// `cell_mass` is the mass of one fluid lattice cell (`m`), `dx`/`row_h` the lattice spacings
    /// and `support` the kernel support radius.
    pub fn build(drum: &Drum, dx: f32, row_h: f32, cell_mass: f32) -> Self {
        let r = drum.radius_m;
        let mut local: Vec<Vec2> = Vec::new();
        let mut psi: Vec<f32> = Vec::new();

        for l in 0..RINGS {
            let radius = r + (l as f32 + 0.5) * row_h;
            let count = ((std::f32::consts::TAU * radius / dx).round() as usize).max(8);
            let arc = std::f32::consts::TAU * radius / count as f32;
            let offset = if l % 2 == 1 { 0.5 } else { 0.0 };
            for k in 0..count {
                let a = std::f32::consts::TAU * (k as f32 + offset) / count as f32;
                local.push(Vec2::new(a.cos(), a.sin()) * radius);
                psi.push(cell_mass * arc / dx);
            }
        }

        if drum.lifters.count > 0 {
            let n_rows = (2.0 * r / row_h).ceil() as i32 + 2;
            for row in -n_rows / 2..=n_rows / 2 {
                let y = row as f32 * row_h;
                let off = if row.rem_euclid(2) == 0 { 0.0 } else { 0.5 };
                let n_cols = (r / dx).ceil() as i32 + 2;
                for col in -n_cols..=n_cols {
                    let p = Vec2::new((col as f32 + off) * dx, y);
                    if p.length() <= r && drum.sdf(p) <= -0.5 * row_h {
                        local.push(p);
                        psi.push(cell_mass);
                    }
                }
            }
        }

        Self {
            local,
            psi,
            radius_m: r,
            lifters: drum.lifters,
        }
    }

    /// Whether this boundary was built for `drum`'s geometry (angular velocity is irrelevant).
    pub fn matches(&self, drum: &Drum) -> bool {
        self.radius_m == drum.radius_m && self.lifters == drum.lifters
    }
}

/// Boundary particles of one ball in the ball's own frame: concentric rings inside the ball at
/// `r - (l + 1/2) row_h` (the mirror image of the wall construction: the ball surface lies halfway
/// between the last fluid row and the first solid row), at most [`RINGS`] of them (deeper
/// particles are beyond the kernel support), plus a centre particle for a core the rings leave
/// uncovered. Each carries its geometric volume `psi = m * cell_area / (dx * row_h)`.
pub struct BallTemplate {
    pub radius: f32,
    pub local: Vec<Vec2>,
    pub psi: Vec<f32>,
}

impl BallTemplate {
    pub fn empty() -> Self {
        Self {
            radius: 0.0,
            local: Vec::new(),
            psi: Vec::new(),
        }
    }

    pub fn build(radius: f32, dx: f32, row_h: f32, cell_mass: f32) -> Self {
        let mut local = Vec::new();
        let mut psi = Vec::new();
        let mut last_ring = radius;
        for l in 0..RINGS {
            let rho = radius - (l as f32 + 0.5) * row_h;
            if rho <= 0.25 * dx {
                break;
            }
            let count = ((std::f32::consts::TAU * rho / dx).round() as usize).max(1);
            let arc = std::f32::consts::TAU * rho / count as f32;
            let offset = if l % 2 == 1 { 0.5 } else { 0.0 };
            for k in 0..count {
                let a = std::f32::consts::TAU * (k as f32 + offset) / count as f32;
                local.push(Vec2::new(a.cos(), a.sin()) * rho);
                psi.push(cell_mass * arc / dx);
            }
            last_ring = rho;
        }
        // Core left inside the innermost ring (a one-particle ball for very small radii).
        if last_ring - 0.5 * row_h > 0.5 * dx || local.is_empty() {
            local.push(Vec2::ZERO);
            psi.push(cell_mass);
        }
        Self { radius, local, psi }
    }
}
