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
use crate::grid::UniformGrid;
use crate::params::LiftersParams;

/// Number of wall rings; 3 rows span `2.6 dx`, beyond the kernel support `2 dx`.
const RINGS: usize = 3;

pub struct Boundary {
    /// Positions in the drum-local frame.
    pub local: Vec<Vec2>,
    /// Per-particle boundary volume (mass-like, kg per metre of depth).
    pub psi: Vec<f32>,
    grid: UniformGrid,
    radius_m: f32,
    lifters: LiftersParams,
}

impl Boundary {
    /// `cell_mass` is the mass of one fluid lattice cell (`m`), `dx`/`row_h` the lattice spacings
    /// and `support` the kernel support radius.
    pub fn build(drum: &Drum, dx: f32, row_h: f32, support: f32, cell_mass: f32) -> Self {
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

        let grid = UniformGrid::build(&local, support);
        Self {
            local,
            psi,
            grid,
            radius_m: r,
            lifters: drum.lifters,
        }
    }

    /// Whether this boundary was built for `drum`'s geometry (angular velocity is irrelevant).
    pub fn matches(&self, drum: &Drum) -> bool {
        self.radius_m == drum.radius_m && self.lifters == drum.lifters
    }

    pub fn for_each_near_local<F: FnMut(u32)>(&self, p_local: Vec2, f: F) {
        self.grid.for_each_near(p_local, f);
    }
}
