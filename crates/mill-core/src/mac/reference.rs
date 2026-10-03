//! Independent reference solutions for the free-surface experiments.
//!
//! * `half_disc_sloshing`: linear potential-flow sloshing eigenvalues `K = omega^2 / g` of a
//!   half-full circular container (2D), by a Rayleigh-Ritz method on the harmonic basis
//!   `Re z^k, Im z^k` (the Neumann wall condition is the natural condition of the weak form).
//! * `rectangular_sloshing_omega`: Lamb's exact dispersion relation for a rectangular tank.

use std::f64::consts::PI;

/// `omega` of mode `mode` (1 = lowest) of a rectangular tank of width `width` and depth `depth`.
pub fn rectangular_sloshing_omega(g: f64, width: f64, depth: f64, mode: usize) -> f64 {
    let k = mode as f64 * PI / width;
    (g * k * (k * depth).tanh()).sqrt()
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`.
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p1, mut p2) = (1.0, 0.0);
            for j in 0..n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * j + 1) as f64 * z * p2 - j as f64 * p3) / (j + 1) as f64;
            }
            dp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let dz = p1 / dp;
            z -= dz;
            if dz.abs() < 1e-15 {
                break;
            }
        }
        x[i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
    (x, w)
}

/// Cholesky factor `L` (lower, row-major) of a symmetric positive definite matrix.
fn cholesky(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i * n + i] = s.sqrt();
            } else {
                l[i * n + j] = s / l[j * n + j];
            }
        }
    }
    Some(l)
}

/// Eigenvalues of a symmetric matrix (cyclic Jacobi), ascending.
fn jacobi_eigenvalues(mut a: Vec<f64>, n: usize) -> Vec<f64> {
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..i).map(move |j| (i, j)))
            .map(|(i, j)| a[i * n + j] * a[i * n + j])
            .sum();
        if off < 1e-28 {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let (c, s) = (1.0 / (t * t + 1.0).sqrt(), t / (t * t + 1.0).sqrt());
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
            }
        }
    }
    let mut ev: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    ev.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ev
}

/// The lowest `count` sloshing eigenvalues `K R = omega^2 R / g` (ascending) of a half-full
/// circular container of radius `R` (free surface along the diameter), from `terms` harmonic
/// polynomial pairs. Fluid occupies `y < 0`.
pub fn half_disc_sloshing(terms: usize, count: usize) -> Vec<f64> {
    // Work with R = 1: coordinates scaled by R, so K R is dimensionless.
    let nb = 2 * terms;
    let (gx, gw) = gauss_legendre(80);
    // Basis index b: k = b / 2 + 1, even b = Re z^k, odd b = Im z^k. Gradient of Re f = (Re f', -Im f'),
    // of Im f = (Im f', Re f') with f' = k z^(k-1).
    let grad = |b: usize, x: f64, y: f64| -> (f64, f64) {
        let k = b / 2 + 1;
        let (r, th) = ((x * x + y * y).sqrt(), y.atan2(x));
        let m = k as f64 * r.powi(k as i32 - 1);
        let (re, im) = (
            m * (th * (k as f64 - 1.0)).cos(),
            m * (th * (k as f64 - 1.0)).sin(),
        );
        if b.is_multiple_of(2) {
            (re, -im)
        } else {
            (im, re)
        }
    };
    let value = |b: usize, x: f64, y: f64| -> f64 {
        let k = b / 2 + 1;
        let (r, th) = ((x * x + y * y).sqrt(), y.atan2(x));
        let m = r.powi(k as i32);
        if b.is_multiple_of(2) {
            m * (th * k as f64).cos()
        } else {
            m * (th * k as f64).sin()
        }
    };
    let mut a = vec![0.0; nb * nb];
    for (ir, &xr) in gx.iter().enumerate() {
        let r = 0.5 * (xr + 1.0);
        for (it, &xt) in gx.iter().enumerate() {
            let th = -PI * 0.5 * (xt + 1.0);
            let wgt = gw[ir] * 0.5 * gw[it] * 0.5 * PI * r;
            let (x, y) = (r * th.cos(), r * th.sin());
            let g: Vec<(f64, f64)> = (0..nb).map(|b| grad(b, x, y)).collect();
            for i in 0..nb {
                for j in 0..=i {
                    a[i * nb + j] += wgt * (g[i].0 * g[j].0 + g[i].1 * g[j].1);
                }
            }
        }
    }
    for i in 0..nb {
        for j in 0..i {
            a[j * nb + i] = a[i * nb + j];
        }
    }
    // Surface (y = 0, x in [-1, 1]): zero-mean shifted values.
    let surf: Vec<Vec<f64>> = (0..nb)
        .map(|b| {
            let vals: Vec<f64> = gx.iter().map(|&x| value(b, x, 0.0)).collect();
            let mean: f64 = vals.iter().zip(&gw).map(|(v, w)| v * w).sum::<f64>() / 2.0;
            vals.iter().map(|v| v - mean).collect()
        })
        .collect();
    let mut bm = vec![0.0; nb * nb];
    for i in 0..nb {
        for j in 0..nb {
            bm[i * nb + j] = (0..gx.len()).map(|q| gw[q] * surf[i][q] * surf[j][q]).sum();
        }
    }
    // Generalised problem B x = mu A x, mu = 1 / K: with A = L L^T, solve (L^-1 B L^-T) y = mu y.
    let l = cholesky(&a, nb).expect("Ritz stiffness matrix not positive definite");
    // C = L^-1 B L^-T
    let mut tmp = vec![0.0; nb * nb]; // L^-1 B
    for col in 0..nb {
        for i in 0..nb {
            let mut s = bm[i * nb + col];
            for k in 0..i {
                s -= l[i * nb + k] * tmp[k * nb + col];
            }
            tmp[i * nb + col] = s / l[i * nb + i];
        }
    }
    let mut c = vec![0.0; nb * nb]; // (L^-1 B) L^-T: solve rows
    for row in 0..nb {
        for i in 0..nb {
            let mut s = tmp[row * nb + i];
            for k in 0..i {
                s -= l[i * nb + k] * c[row * nb + k];
            }
            c[row * nb + i] = s / l[i * nb + i];
        }
    }
    for i in 0..nb {
        for j in 0..i {
            let m = 0.5 * (c[i * nb + j] + c[j * nb + i]);
            c[i * nb + j] = m;
            c[j * nb + i] = m;
        }
    }
    let mut mu = jacobi_eigenvalues(c, nb);
    mu.reverse();
    mu.iter()
        .filter(|&&m| m > 1e-9)
        .take(count)
        .map(|&m| 1.0 / m)
        .collect()
}
