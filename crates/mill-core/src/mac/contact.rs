//! Solid contact between discs and with the drum wall, below the lubrication film.
//!
//! The squeeze film (`lubrication.rs`) stops an approach only asymptotically; real surfaces touch
//! at a roughness height. Contacts are active when the gap drops below `gap0` and act on the
//! velocities at the end of the step: a one-sided, implicit normal constraint (a stiff dashpot
//! with a Baumgarte-type push-out of the penetration, solved with the lubrication links) and
//! Coulomb friction with a regularised stick (explicit in the coupling iteration).

/// Penalty stiffness of the normal constraint relative to `m_eff / dt`.
pub const STIFFNESS: f64 = 100.0;

/// Material parameters of a contact.
#[derive(Clone, Copy, Debug)]
pub struct ContactParams {
    /// Gap (m) below which a contact is considered (the constraint itself acts at gap 0).
    pub gap0: f64,
    /// Fraction of an existing overlap removed per step.
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
    /// Surface gap (m); negative when overlapping.
    pub gap: f64,
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
                gap,
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
                gap,
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

    /// Push-out speed (m/s): the contact allows an approach of at most `gap / dt` (no overlap at
    /// the end of the step) and removes a fraction `beta` of an existing overlap per step.
    pub fn push_out(&self, beta: f64, dt: f64) -> f64 {
        if self.gap >= 0.0 {
            -self.gap / dt
        } else {
            -beta * self.gap / dt
        }
    }

    /// Stiffness of the normal constraint (force per unit speed beyond the allowed approach):
    /// `STIFFNESS` times the mass that stops in one step.
    pub fn dashpot(&self, dt: f64) -> f64 {
        STIFFNESS * self.m_eff / dt
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

/// Stiffness of the regularised friction (force per unit slip speed).
pub fn friction_stiffness(c: &Contact, dt: f64) -> f64 {
    0.5 * c.m_eff / dt
}

/// Per-disc linear stiffness `(translation, rotation)` of the friction of all contacts; it is added
/// to the coupling's preconditioner so that the explicit friction does not make the fixed point
/// overshoot.
pub fn friction_preconditioner(list: &[Contact], count: usize, dt: f64) -> Vec<(f64, f64)> {
    let mut out = vec![(0.0, 0.0); count];
    for c in list {
        let k = friction_stiffness(c, dt);
        out[c.i].0 += k;
        out[c.i].1 += c.ri * c.ri * k;
        if let Some(j) = c.j {
            out[j].0 += k;
            out[j].1 += c.rj * c.rj * k;
        }
    }
    out
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
        let kappa = friction_stiffness(c, dt);
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
