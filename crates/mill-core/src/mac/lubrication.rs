//! Sub-grid lubrication for near-contacting discs and for a disc near the drum wall.
//!
//! The resolved flow represents a squeeze film correctly down to a gap of about 4 cells and
//! saturates below (the force plateaus at roughly the value of a 2-cell gap). Below that, the
//! normal force along the line of centres is blended towards the 2D Reynolds squeeze-film model
//! with its first corrections: `F = c(h) u_n`,
//! `c = 3 sqrt(2) pi nu R^1.5 h^-1.5 (1 + c1 eps + c2 eps^2)`, `eps = h / R`, where `R` is the
//! reduced curvature radius, `1 / R = k1 + k2` with signed curvatures (`k = 1/a` for a disc, `k = -1/b`
//! for the concave drum wall). The corrections follow the quartic term of the gap profile,
//! `Q = R^3 (k1^3 + k2^3)`: `c1 = 0.788 + 0.265 Q`, `c2 = 0.006 + 0.057 (Q - 0.25)` (>= 0). They were
//! calibrated on the exact bipolar-coordinate solutions (disc against the drum wall,
//! `verify::eccentric_squeeze_force`, and two equal discs in unbounded fluid,
//! `verify::disc_pair_squeeze_force`): both are reproduced to 1 % for `eps <= 0.4` and 2.5 % at 0.8.
//! Units: fluid density 1, `nu` kinematic viscosity.

use std::f64::consts::PI;

/// Gap, in cells, above which the resolved force is used alone, and below `H_FULL` the model alone.
pub const H_START: f64 = 8.0;
pub const H_FULL: f64 = 3.0;

/// Signed curvature of a disc of radius `a`.
pub fn disc_curvature(a: f64) -> f64 {
    1.0 / a
}

/// Signed curvature of the concave drum wall of radius `b` (seen from inside).
pub fn wall_curvature(b: f64) -> f64 {
    -1.0 / b
}

/// Squeeze-film damping coefficient: force per unit approach speed at gap `h` between surfaces of
/// signed curvatures `k1` and `k2` (sum positive).
pub fn squeeze_coefficient(k1: f64, k2: f64, h: f64, nu: f64) -> f64 {
    let r = 1.0 / (k1 + k2);
    let eps = h / r;
    let q = r.powi(3) * (k1.powi(3) + k2.powi(3));
    let c1 = 0.788 + 0.265 * q;
    let c2 = (0.006 + 0.057 * (q - 0.25)).max(0.0);
    3.0 * 2f64.sqrt() * PI * nu * r.powf(1.5) / h.powf(1.5) * (1.0 + c1 * eps + c2 * eps * eps)
}

/// Shear (Couette) damping of a sliding film: force per unit relative surface speed,
/// `nu pi sqrt(2 r / h)` with `r = 1 / (k1 + k2)` (leading order for small gaps).
pub fn tangential_coefficient(k1: f64, k2: f64, h: f64, nu: f64) -> f64 {
    let r = 1.0 / (k1 + k2);
    nu * PI * (2.0 * r / h).sqrt()
}

/// Weight of the model in the blend, `0` for gaps of at least `H_START` cells and `1` below
/// `H_FULL` cells (smoothstep in between).
pub fn model_weight(h: f64, dx: f64) -> f64 {
    let t = ((H_START * dx - h) / ((H_START - H_FULL) * dx)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Fraction of the exact squeeze-film coefficient that the resolved flow delivers at a gap of
/// `h_cells` cells (measured on the wall and pair references, n = 128 and 256, E5a/E5e): about 1
/// down to 4 cells, 0.9 at 2, 0.25 at 1 and zero for touching surfaces (the values at 1 - 2 cells
/// scatter by +-0.05 with the disc size and position relative to the grid).
pub fn resolved_ratio(h_cells: f64) -> f64 {
    const TABLE: [(f64, f64); 6] = [
        (0.0, 0.0),
        (0.5, 0.06),
        (1.0, 0.25),
        (2.0, 0.90),
        (4.0, 0.99),
        (8.0, 1.0),
    ];
    if h_cells >= 8.0 {
        return 1.0;
    }
    let h = h_cells.max(0.0);
    for w in TABLE.windows(2) {
        if h <= w[1].0 {
            let t = (h - w[0].0) / (w[1].0 - w[0].0);
            return w[0].1 + t * (w[1].1 - w[0].1);
        }
    }
    1.0
}

/// Normal force opposing the approach after the blend: `grid_normal` is the resolved component,
/// `u_n` the approach speed (relative, positive when the gap closes).
#[allow(clippy::too_many_arguments)]
pub fn blended_normal_force(
    grid_normal: f64,
    k1: f64,
    k2: f64,
    h: f64,
    dx: f64,
    nu: f64,
    u_n: f64,
) -> f64 {
    let s = model_weight(h, dx);
    if s == 0.0 {
        return grid_normal;
    }
    let h = h.max(1e-6 / (k1 + k2));
    (1.0 - s) * grid_normal + s * squeeze_coefficient(k1, k2, h, nu) * u_n
}

/// Applies the wall squeeze-film blend to the resolved load `(fx, fy)` on a disc whose centre is at
/// `(cx, cy)` with velocity `(ux, uy)` inside the drum of radius `drum`: only the component along
/// the line of centres is blended. Returns the corrected `(fx, fy)`.
#[allow(clippy::too_many_arguments)]
pub fn wall_adjusted_load(
    load: (f64, f64),
    centre: (f64, f64),
    velocity: (f64, f64),
    a: f64,
    drum: f64,
    dx: f64,
    nu: f64,
) -> (f64, f64) {
    let dist = (centre.0 * centre.0 + centre.1 * centre.1).sqrt();
    let h = drum - a - dist;
    if h >= H_START * dx || dist < 1e-12 {
        return load;
    }
    // Unit vector from the drum centre towards the disc (the wall lies beyond).
    let (nx, ny) = (centre.0 / dist, centre.1 / dist);
    let approach = velocity.0 * nx + velocity.1 * ny;
    let grid_normal = -(load.0 * nx + load.1 * ny);
    let total = blended_normal_force(
        grid_normal,
        disc_curvature(a),
        wall_curvature(drum),
        h,
        dx,
        nu,
        approach,
    );
    let delta = grid_normal - total;
    (load.0 + delta * nx, load.1 + delta * ny)
}

/// One lubricated contact: disc `i` against disc `j`, or against the drum wall (`j = None`).
/// `n` is the unit vector from `i` towards the partner (outward for the wall), `weight` the model
/// weight of the correction and `c` the squeeze coefficient (force per unit approach speed).
/// `weight = 1 - resolved_ratio(gap / dx)`: the sub-grid force `-weight c u_rel,n` is *added* to
/// the resolved load, which keeps the common-mode (non-squeezing) part of the resolved force.
#[derive(Clone, Copy, Debug)]
pub struct Link {
    pub i: usize,
    pub j: Option<usize>,
    pub n: (f64, f64),
    pub weight: f64,
    pub c: f64,
    /// Shear coefficient (force per unit tangential slip speed, before `weight`).
    pub ct: f64,
    pub gap: f64,
}

/// All contacts with a gap below `H_START` cells for discs with the given centres and radii
/// inside the drum of radius `drum`.
pub fn links(centres: &[(f64, f64)], radii: &[f64], drum: f64, dx: f64, nu: f64) -> Vec<Link> {
    let mut out = Vec::new();
    for i in 0..centres.len() {
        for j in i + 1..centres.len() {
            let (ex, ey) = (centres[j].0 - centres[i].0, centres[j].1 - centres[i].1);
            let dist = (ex * ex + ey * ey).sqrt();
            let gap = dist - radii[i] - radii[j];
            if gap >= H_START * dx || dist < 1e-12 {
                continue;
            }
            let h = gap.max(1e-6 * radii[i].min(radii[j]));
            out.push(Link {
                i,
                j: Some(j),
                n: (ex / dist, ey / dist),
                weight: 1.0 - resolved_ratio(gap / dx),
                c: squeeze_coefficient(disc_curvature(radii[i]), disc_curvature(radii[j]), h, nu),
                gap,
                ct: tangential_coefficient(
                    disc_curvature(radii[i]),
                    disc_curvature(radii[j]),
                    h,
                    nu,
                ),
            });
        }
        let (cx, cy) = centres[i];
        let dist = (cx * cx + cy * cy).sqrt();
        let gap = drum - radii[i] - dist;
        if gap < H_START * dx && dist > 1e-12 {
            let h = gap.max(1e-6 * radii[i]);
            out.push(Link {
                i,
                j: None,
                n: (cx / dist, cy / dist),
                weight: 1.0 - resolved_ratio(gap / dx),
                c: squeeze_coefficient(disc_curvature(radii[i]), wall_curvature(drum), h, nu),
                gap,
                ct: tangential_coefficient(disc_curvature(radii[i]), wall_curvature(drum), h, nu),
            });
        }
    }
    out
}

/// Removes the model's share of the resolved normal force on every linked disc.
pub fn remove_grid_normal(forces: &mut [(f64, f64)], links: &[Link]) {
    for l in links {
        for body in [Some(l.i), l.j].into_iter().flatten() {
            let f = forces[body];
            let fn_ = f.0 * l.n.0 + f.1 * l.n.1;
            forces[body] = (f.0 - l.weight * fn_ * l.n.0, f.1 - l.weight * fn_ * l.n.1);
        }
    }
}

/// Model forces `-weight c (u_i - u_j) . n n` on the linked discs for the given velocities.
pub fn model_forces(links: &[Link], velocities: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut out = vec![(0.0, 0.0); velocities.len()];
    for l in links {
        let ui = velocities[l.i];
        let uj = l.j.map_or((0.0, 0.0), |j| velocities[j]);
        let approach = (ui.0 - uj.0) * l.n.0 + (ui.1 - uj.1) * l.n.1;
        let f = (
            -l.weight * l.c * approach * l.n.0,
            -l.weight * l.c * approach * l.n.1,
        );
        out[l.i].0 += f.0;
        out[l.i].1 += f.1;
        if let Some(j) = l.j {
            out[j].0 -= f.0;
            out[j].1 -= f.1;
        }
    }
    out
}
