//! SPH kernels for the DFSPH fluid: the 2D cubic spline (density, pressure gradient, viscosity,
//! boundary volume) and the Akinci cohesion kernel (surface tension).

use glam::Vec2;

/// 2D cubic-spline kernel with compact support radius `support` (`H`).
///
/// `W(q) = sigma * (6 (q^3 - q^2) + 1)` for `q <= 1/2`, `sigma * 2 (1 - q)^3` for `1/2 < q <= 1`,
/// `q = r / H`, `sigma = 40 / (7 pi H^2)` (integrates to 1 over the disc of radius `H`).
#[derive(Debug, Clone, Copy)]
pub struct Kernel {
    pub support: f32,
    sigma: f32,
}

impl Kernel {
    pub fn new(support: f32) -> Self {
        Self {
            support,
            sigma: 40.0 / (7.0 * std::f32::consts::PI * support * support),
        }
    }

    #[inline]
    pub fn w(&self, r: f32) -> f32 {
        let q = r / self.support;
        if q >= 1.0 {
            0.0
        } else if q <= 0.5 {
            self.sigma * (6.0 * (q * q * q - q * q) + 1.0)
        } else {
            let t = 1.0 - q;
            self.sigma * 2.0 * t * t * t
        }
    }

    /// Radial derivative `dW/dr` (never positive).
    #[inline]
    pub fn dw_dr(&self, r: f32) -> f32 {
        let q = r / self.support;
        if q >= 1.0 {
            0.0
        } else if q <= 0.5 {
            self.sigma / self.support * 6.0 * (3.0 * q * q - 2.0 * q)
        } else {
            let t = 1.0 - q;
            -self.sigma / self.support * 6.0 * t * t
        }
    }

    /// Gradient with respect to the first argument; `delta = x_i - x_j`, `r = |delta|`.
    #[inline]
    pub fn grad(&self, delta: Vec2, r: f32) -> Vec2 {
        if r <= 1e-9 || r >= self.support {
            Vec2::ZERO
        } else {
            delta * (self.dw_dr(r) / r)
        }
    }
}

/// Akinci, Akinci & Teschner (2013) cohesion kernel, normalised for 2D (identical to the PBF
/// solver's, so `slurry.surface_tension_n_m` keeps its calibrated meaning). Zero at `r = 0` and
/// `r = h`, peaks near `0.6 h`, slightly negative below `h / 2`.
pub fn cohesion_kernel(r: f32, h: f32) -> f32 {
    if r <= 0.0 || r >= h {
        return 0.0;
    }
    let coeff = 35840.0 / (209.0 * std::f32::consts::PI * h.powi(8));
    let term = (h - r).powi(3) * r.powi(3);
    if r <= 0.5 * h {
        coeff * (2.0 * term - h.powi(6) / 64.0)
    } else {
        coeff * term
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_spline_integrates_to_one_over_the_disc() {
        let k = Kernel::new(0.01);
        let n = 20_000;
        let dr = k.support / n as f32;
        let mut sum = 0.0f64;
        for i in 0..n {
            let r = (i as f32 + 0.5) * dr;
            sum += (2.0 * std::f32::consts::PI * r * k.w(r) * dr) as f64;
        }
        assert!((sum - 1.0).abs() < 1e-3, "integral {sum}");
    }

    #[test]
    fn gradient_matches_finite_difference_and_points_inward() {
        let k = Kernel::new(0.01);
        for &r in &[0.001f32, 0.004, 0.005, 0.006, 0.009] {
            let eps = 1e-6;
            let fd = (k.w(r + eps) - k.w(r - eps)) / (2.0 * eps);
            let an = k.dw_dr(r);
            assert!(
                (fd - an).abs() <= 1e-2 * an.abs() + 1.0,
                "r {r}: fd {fd} analytic {an}"
            );
            assert!(an <= 0.0);
        }
        let g = k.grad(Vec2::new(0.004, 0.0), 0.004);
        assert!(g.x < 0.0 && g.y == 0.0);
    }
}
