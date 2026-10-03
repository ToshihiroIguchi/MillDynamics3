//! Generalised-Newtonian (Herschel-Bulkley) constitutive law for the augmented-Lagrangian
//! (ALG2) treatment of yield-stress fluids.
//!
//! Everything is kinematic (stress divided by density). The stress of a shear flow is
//! `tau = tau_y + k * rate^n` for `rate > 0`; in tensor form `tau = 2 eta(rate) D` with the
//! shear rate `rate = sqrt(2 D:D)`. Newtonian: `tau_y = 0, n = 1, k = nu`; Bingham: `n = 1`;
//! power law: `tau_y = 0`.

use std::f64::consts::SQRT_2;

#[derive(Clone, Copy, Debug)]
pub struct HerschelBulkley {
    /// Yield stress over density (m^2/s^2).
    pub tau_y: f64,
    /// Consistency over density (m^2 s^(n-2)).
    pub k: f64,
    /// Flow index.
    pub n: f64,
}

impl HerschelBulkley {
    pub fn newtonian(nu: f64) -> Self {
        Self {
            tau_y: 0.0,
            k: nu,
            n: 1.0,
        }
    }

    pub fn bingham(tau_y: f64, nu_p: f64) -> Self {
        Self {
            tau_y,
            k: nu_p,
            n: 1.0,
        }
    }

    /// Shear stress of a simple shear at `rate`.
    pub fn stress(&self, rate: f64) -> f64 {
        if rate <= 0.0 {
            return 0.0;
        }
        self.tau_y + self.k * rate.powf(self.n)
    }

    /// Shear rate that the d-step of ALG2 assigns to a trial stress tensor of Frobenius norm `z`
    /// with augmentation `r`: the root of `sqrt(2) (k g^n + tau_y) + r g / sqrt(2) = z`, or zero
    /// when `z <= sqrt(2) tau_y` (unyielded).
    pub fn trial_rate(&self, z: f64, r: f64) -> f64 {
        let yield_level = SQRT_2 * self.tau_y;
        if z <= yield_level {
            return 0.0;
        }
        if self.n == 1.0 {
            return (z - yield_level) / (SQRT_2 * self.k + r / SQRT_2);
        }
        let f = |g: f64| SQRT_2 * (self.k * g.powf(self.n) + self.tau_y) + r * g / SQRT_2 - z;
        let (mut lo, mut hi) = (0.0, SQRT_2 * z / r);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if f(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

/// Symmetric 2D tensor `[xx, xy, yy]`.
pub type Sym = [f64; 3];

/// Frobenius norm of a symmetric 2D tensor.
pub fn norm(t: &Sym) -> f64 {
    (t[0] * t[0] + 2.0 * t[1] * t[1] + t[2] * t[2]).sqrt()
}

/// One ALG2 update at a node: given the multiplier `lam` and the strain-rate tensor `d_u`
/// (`D(u)`), returns the new multiplier and the stress excess `S = lam_new - r d` whose divergence
/// is the explicit right-hand side of the next velocity solve. `omega` under-relaxes the
/// multiplier step (1 = ALG2); the fixed point does not depend on it.
pub fn alg2_update(law: &HerschelBulkley, r: f64, omega: f64, lam: &Sym, d_u: &Sym) -> (Sym, Sym) {
    let z = [
        lam[0] + r * d_u[0],
        lam[1] + r * d_u[1],
        lam[2] + r * d_u[2],
    ];
    let zn = norm(&z);
    let rate = law.trial_rate(zn, r);
    let scale = if zn > 0.0 { rate / SQRT_2 / zn } else { 0.0 };
    let d = [z[0] * scale, z[1] * scale, z[2] * scale];
    let new = [
        lam[0] + omega * r * (d_u[0] - d[0]),
        lam[1] + omega * r * (d_u[1] - d[1]),
        lam[2] + omega * r * (d_u[2] - d[2]),
    ];
    let s = [new[0] - r * d[0], new[1] - r * d[1], new[2] - r * d[2]];
    (new, s)
}
