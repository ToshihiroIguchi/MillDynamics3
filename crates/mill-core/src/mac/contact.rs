//! Solid contact between discs and with the drum wall, below the lubrication film.
//!
//! The squeeze film (`lubrication.rs`) stops an approach only asymptotically; real surfaces touch
//! at a roughness height. Contacts are active when the gap drops below `gap0` and act on the
//! velocities at the end of the step: a one-sided, implicit normal constraint (a stiff dashpot
//! with a Baumgarte-type push-out of the penetration, solved with the lubrication links) and
//! Coulomb friction with a regularised stick (explicit in the coupling iteration).

/// Material parameters of a contact.
#[derive(Clone, Copy, Debug)]
pub struct ContactParams {
    /// Roughness gap (m) below which solid contact starts.
    pub gap0: f64,
    /// Fraction of the penetration removed per step.
    pub beta: f64,
    /// Coulomb friction coefficient (0 disables friction).
    pub mu: f64,
}

/// One contact: disc `i` against disc `j`, or against the drum wall (`j = None`).
#[derive(Clone, Copy, Debug)]
pub struct Contact {
    pub i: usize,
    pub j: Option<usize>,
    /// Unit vector from `i` towards the partner (outward for the wall).
    pub n: (f64, f64),
    /// Penetration below `gap0` (m, >= 0).
    pub depth: f64,
    /// Effective mass (per unit density) of the pair.
    pub m_eff: f64,
    pub ri: f64,
    pub rj: f64,
}

/// Contacts for discs at the given predicted centres inside the drum.
pub fn contacts(
    centres: &[(f64, f64)],
    radii: &[f64],
    masses: &[f64],
    drum: f64,
    p: &ContactParams,
) -> Vec<Contact> {
    let mut out = Vec::new();
    for i in 0..centres.len() {
        for j in i + 1..centres.len() {
            let (ex, ey) = (centres[j].0 - centres[i].0, centres[j].1 - centres[i].1);
            let dist = (ex * ex + ey * ey).sqrt();
            let gap = dist - radii[i] - radii[j];
            if gap >= p.gap0 || dist < 1e-12 {
                continue;
            }
            out.push(Contact {
                i,
                j: Some(j),
                n: (ex / dist, ey / dist),
                depth: p.gap0 - gap,
                m_eff: masses[i] * masses[j] / (masses[i] + masses[j]),
                ri: radii[i],
                rj: radii[j],
            });
        }
        let (cx, cy) = centres[i];
        let dist = (cx * cx + cy * cy).sqrt();
        let gap = drum - radii[i] - dist;
        if gap < p.gap0 && dist > 1e-12 {
            out.push(Contact {
                i,
                j: None,
                n: (cx / dist, cy / dist),
                depth: p.gap0 - gap,
                m_eff: masses[i],
                ri: radii[i],
                rj: 0.0,
            });
        }
    }
    out
}

impl Contact {
    /// Normal approach speed (positive when closing).
    pub fn approach(&self, v: &[(f64, f64, f64)]) -> f64 {
        let vi = v[self.i];
        let vj = self.j.map_or((0.0, 0.0, 0.0), |j| v[j]);
        (vi.0 - vj.0) * self.n.0 + (vi.1 - vj.1) * self.n.1
    }

    /// Normal push-out speed target (m/s) and dashpot coefficient for a step of `dt`.
    pub fn push_out(&self, beta: f64, dt: f64) -> f64 {
        beta * self.depth / dt
    }

    pub fn dashpot(&self, dt: f64) -> f64 {
        self.m_eff / dt
    }

    /// Tangential slip speed at the contact point; `wall_speed` is the drum surface speed
    /// (counter-clockwise positive) used for a wall contact.
    pub fn slip(&self, v: &[(f64, f64, f64)], wall_speed: f64) -> f64 {
        let t = (-self.n.1, self.n.0);
        let vi = v[self.i];
        match self.j {
            Some(j) => {
                let vj = v[j];
                (vi.0 - vj.0) * t.0 + (vi.1 - vj.1) * t.1 + self.ri * vi.2 + self.rj * vj.2
            }
            None => vi.0 * t.0 + vi.1 * t.1 + self.ri * vi.2 - wall_speed,
        }
    }
}

/// Contact forces and torques on every disc for the velocities `v`: Coulomb friction with a
/// regularised stick (`kappa = m_eff / dt` per unit slip), limited by `mu` times the normal force
/// of the one-sided dashpot. Returns `(fx, fy, torque)` per disc.
pub fn friction_loads(
    list: &[Contact],
    v: &[(f64, f64, f64)],
    p: &ContactParams,
    dt: f64,
    wall_speed: f64,
) -> Vec<(f64, f64, f64)> {
    let mut out = vec![(0.0, 0.0, 0.0); v.len()];
    if p.mu <= 0.0 {
        return out;
    }
    for c in list {
        let fnorm = (c.dashpot(dt) * (c.approach(v) + c.push_out(p.beta, dt))).max(0.0);
        if fnorm <= 0.0 {
            continue;
        }
        let kappa = 0.5 * c.m_eff / dt;
        let slip = c.slip(v, wall_speed);
        let ft = (kappa * slip).clamp(-p.mu * fnorm, p.mu * fnorm);
        let t = (-c.n.1, c.n.0);
        out[c.i].0 -= ft * t.0;
        out[c.i].1 -= ft * t.1;
        out[c.i].2 -= c.ri * ft;
        if let Some(j) = c.j {
            out[j].0 += ft * t.0;
            out[j].1 += ft * t.1;
            out[j].2 -= c.rj * ft;
        }
    }
    out
}
