//! Dense linear algebra and the contact/lubrication/friction assembly of the Newton coupling of
//! many bodies (`BodyFlow::step_newton`).
//!
//! Unknowns per body are `(ux, uy, w)` with `w = r * omega`, so that the rotational mass is
//! `inertia / r^2` and a tangential contact speed is `(u_i - u_j) . t + w_i + w_j`.

use super::contact::{Contact, ContactParams};
use super::lubrication::Link;

/// Solves `a x = b` (row-major `n x n`) by Gaussian elimination with partial pivoting; the result
/// replaces `b`. Returns false for a singular matrix.
pub fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> bool {
    for c in 0..n {
        let mut p = c;
        for r in c + 1..n {
            if a[r * n + c].abs() > a[p * n + c].abs() {
                p = r;
            }
        }
        if a[p * n + c].abs() < 1e-300 {
            return false;
        }
        if p != c {
            for k in 0..n {
                a.swap(c * n + k, p * n + k);
            }
            b.swap(c, p);
        }
        let piv = a[c * n + c];
        for r in c + 1..n {
            let f = a[r * n + c] / piv;
            if f != 0.0 {
                for k in c..n {
                    a[r * n + k] -= f * a[c * n + k];
                }
                b[r] -= f * b[c];
            }
        }
    }
    for c in (0..n).rev() {
        let mut s = b[c];
        for k in c + 1..n {
            s -= a[c * n + k] * b[k];
        }
        b[c] = s / a[c * n + c];
    }
    true
}

/// Adds `w * h h^T` to `a` and `rhs_w * h` to `b`, where `h` has the entries `(index, value)`.
fn add_rank_one(a: &mut [f64], b: &mut [f64], n: usize, h: &[(usize, f64)], w: f64, rhs_w: f64) {
    for &(p, hp) in h {
        for &(q, hq) in h {
            a[p * n + q] += w * hp * hq;
        }
        b[p] += rhs_w * hp;
    }
}

/// Everything the dense solve needs for one step.
pub struct StepData<'a> {
    pub nb: usize,
    /// Mass per unknown (`m, m, I / r^2` for every body).
    pub mass: &'a [f64],
    pub dt: f64,
    /// Velocities at the start of the step (unknowns).
    pub start: &'a [f64],
    /// Fluid load at the linearisation point and its derivative `-dF/dV` (dense).
    pub fluid: &'a [f64],
    pub stiffness: &'a [f64],
    /// Unknowns at the linearisation point.
    pub point: &'a [f64],
    /// External (gravity) load per unknown.
    pub external: &'a [f64],
    pub links: &'a [Link],
    pub contacts: &'a [Contact],
    pub params: Option<ContactParams>,
    pub wall_speed: f64,
}

/// Implicit velocities at the end of the step for a linearised fluid: momentum balance with the
/// squeeze-film links, one-sided normal contacts (active set) and Coulomb friction (stick/slip
/// set). Returns `None` for a singular system.
pub fn solve_step(d: &StepData, guess: &[f64]) -> Option<Vec<f64>> {
    let n = 3 * d.nb;
    let mut a0 = d.stiffness.to_vec();
    let mut b0 = vec![0.0; n];
    for i in 0..n {
        a0[i * n + i] += d.mass[i] / d.dt;
        let mut s = d.mass[i] * d.start[i] / d.dt + d.fluid[i] + d.external[i];
        for j in 0..n {
            s += d.stiffness[i * n + j] * d.point[j];
        }
        b0[i] = s;
    }
    // Squeeze-film links: force on i is -w c n n^T (u_i - u_j).
    for l in d.links {
        let mut h = vec![(3 * l.i, l.n.0), (3 * l.i + 1, l.n.1)];
        if let Some(j) = l.j {
            h.push((3 * j, -l.n.0));
            h.push((3 * j + 1, -l.n.1));
        }
        add_rank_one(&mut a0, &mut b0, n, &h, l.weight * l.c, 0.0);
        // Shear of the film: load -ct (g . v - ws) g with g = (t, 1) on i and (-t, 1) on j.
        let t = (-l.n.1, l.n.0);
        let mut g = vec![(3 * l.i, t.0), (3 * l.i + 1, t.1), (3 * l.i + 2, 1.0)];
        let ws = match l.j {
            Some(j) => {
                g.push((3 * j, -t.0));
                g.push((3 * j + 1, -t.1));
                g.push((3 * j + 2, 1.0));
                0.0
            }
            None => d.wall_speed,
        };
        add_rank_one(
            &mut a0,
            &mut b0,
            n,
            &g,
            l.weight * l.ct,
            l.weight * l.ct * ws,
        );
    }
    let Some(p) = d.params else {
        let mut a = a0;
        return solve_dense(&mut a, &mut b0, n).then_some(b0);
    };
    let gamma = |c: &Contact| c.dashpot(d.dt);
    let push = |c: &Contact| c.push_out(p.beta, d.dt);
    let closing = |c: &Contact, v: &[f64]| {
        let vj = c.j.map_or((0.0, 0.0), |j| (v[3 * j], v[3 * j + 1]));
        (v[3 * c.i] - vj.0) * c.n.0 + (v[3 * c.i + 1] - vj.1) * c.n.1 + push(c)
    };
    let mut v = guess.to_vec();
    let mut normal: Vec<bool> = d.contacts.iter().map(|c| closing(c, &v) > 0.0).collect();
    let mut stick: Vec<bool> = vec![true; d.contacts.len()];
    let mut out = None;
    for _ in 0..10 {
        let mut a = a0.clone();
        let mut b = b0.clone();
        for (k, c) in d.contacts.iter().enumerate() {
            if !normal[k] {
                continue;
            }
            let mut h = vec![(3 * c.i, c.n.0), (3 * c.i + 1, c.n.1)];
            if let Some(j) = c.j {
                h.push((3 * j, -c.n.0));
                h.push((3 * j + 1, -c.n.1));
            }
            // Load on i: -gamma n (n.(ui-uj) + push): matrix gamma h h^T, rhs -gamma push h.
            add_rank_one(&mut a, &mut b, n, &h, gamma(c), -gamma(c) * push(c));
            if p.mu > 0.0 {
                let fnorm = (gamma(c) * closing(c, &v)).max(0.0);
                let t = (-c.n.1, c.n.0);
                let mut g = vec![(3 * c.i, t.0), (3 * c.i + 1, t.1), (3 * c.i + 2, 1.0)];
                let ws = if let Some(j) = c.j {
                    g.push((3 * j, -t.0));
                    g.push((3 * j + 1, -t.1));
                    g.push((3 * j + 2, 1.0));
                    0.0
                } else {
                    d.wall_speed
                };
                let kappa = 2.0 * c.m_eff / d.dt;
                let slip: f64 = g.iter().map(|&(q, gq)| gq * v[q]).sum::<f64>() - ws;
                if stick[k] && (kappa * slip).abs() <= p.mu * fnorm {
                    // load = -kappa (g.v - ws) g
                    add_rank_one(&mut a, &mut b, n, &g, kappa, kappa * ws);
                } else {
                    let ft = (p.mu * fnorm).copysign(slip);
                    for &(q, gq) in &g {
                        b[q] -= ft * gq;
                    }
                    stick[k] = false;
                }
            }
        }
        let mut x = b.clone();
        if !solve_dense(&mut a, &mut x, n) {
            return None;
        }
        let next: Vec<bool> = d.contacts.iter().map(|c| closing(c, &x) > 0.0).collect();
        // Return to stick whenever the normal set changes.
        let changed = next != normal;
        v = x.clone();
        if !changed {
            out = Some(x);
            break;
        }
        for s in stick.iter_mut() {
            *s = true;
        }
        normal = next;
        out = Some(x);
    }
    out
}
