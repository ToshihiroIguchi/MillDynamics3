//! Sub-grid lubrication for near-contacting discs and for a disc near the drum wall.
//!
//! The resolved flow represents a squeeze film correctly down to a gap of about 4 cells and
//! saturates below (the force plateaus at roughly the value of a 2-cell gap). Below that, the
//! normal force along the line of centres is blended towards the 2D Reynolds squeeze-film model
//! with the exact first corrections, calibrated against the bipolar-coordinate solution
//! (`verify::eccentric_squeeze_force`): `F = c(h) * u_n`,
//! `c = 3 sqrt(2) pi nu R^1.5 h^-1.5 (1 + 1.1 eps + 0.06 eps^2)`, `eps = h / R`, where `R` is the
//! reduced curvature radius (`1/R = 1/a1 + 1/a2` for two discs, `1/a - 1/b` against the concave
//! drum wall). Units: fluid density 1, `nu` kinematic viscosity.

use std::f64::consts::PI;

/// Gap, in cells, above which the resolved force is used alone, and below `H_FULL` the model alone.
pub const H_START: f64 = 8.0;
pub const H_FULL: f64 = 3.0;

/// Reduced curvature radius of two discs.
pub fn reduced_radius_discs(a1: f64, a2: f64) -> f64 {
    a1 * a2 / (a1 + a2)
}

/// Reduced curvature radius of a disc of radius `a` inside the drum of radius `b`.
pub fn reduced_radius_wall(a: f64, b: f64) -> f64 {
    a * b / (b - a)
}

/// Squeeze-film damping coefficient: force per unit approach speed at gap `h`.
pub fn squeeze_coefficient(r_eff: f64, h: f64, nu: f64) -> f64 {
    let eps = h / r_eff;
    3.0 * 2f64.sqrt() * PI * nu * r_eff.powf(1.5) / h.powf(1.5)
        * (1.0 + 1.1 * eps + 0.06 * eps * eps)
}

/// Weight of the model in the blend, `0` for gaps of at least `H_START` cells and `1` below
/// `H_FULL` cells (smoothstep in between).
pub fn model_weight(h: f64, dx: f64) -> f64 {
    let t = ((H_START * dx - h) / ((H_START - H_FULL) * dx)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Normal force on a disc (positive along the approach direction, i.e. opposing it) after the
/// blend: `grid_normal` is the resolved force component opposing the approach, `u_n` the
/// approach speed (relative, positive when the gap closes).
pub fn blended_normal_force(
    grid_normal: f64,
    r_eff: f64,
    h: f64,
    dx: f64,
    nu: f64,
    u_n: f64,
) -> f64 {
    let s = model_weight(h, dx);
    if s == 0.0 {
        return grid_normal;
    }
    let h = h.max(1e-6 * r_eff);
    (1.0 - s) * grid_normal + s * squeeze_coefficient(r_eff, h, nu) * u_n
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
        reduced_radius_wall(a, drum),
        h,
        dx,
        nu,
        approach,
    );
    let delta = grid_normal - total;
    (load.0 + delta * nx, load.1 + delta * ny)
}
