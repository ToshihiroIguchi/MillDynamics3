//! Rigid discs in the fluid: loads on a disc, and a fluid-disc coupling that stays stable for
//! light bodies (density ratio near 1, where the added mass of the fluid exceeds the body mass).
//!
//! The mesh is rebuilt around the disc every step ([`StaggeredMesh`] with a body flux), the flow
//! is carried over with [`FlowState`], and loads are the pressure on the cut cells plus the wall
//! traction of the flow relative to the body's rigid-body velocity. Units: fluid density 1.

use super::staggered::{Disc, FlowState, StaggeredFlow, StaggeredMesh};
/// Force and torque of the fluid on a disc (per unit density): pressure on the cut cells moved to
/// the wall point along the local gradient, plus the viscous traction of the flow relative to the
/// disc's rigid-body velocity (zero Dirichlet value at its surface, where the component-wise and
/// the full-stress tractions coincide). Torque is about the disc centre.
pub fn disc_load(flow: &StaggeredFlow, disc: &Disc) -> (f64, f64, f64) {
    let sm = flow.sm;
    let n = sm.n();
    let dx = sm.dx();
    let half = sm.mesh.grid.half;
    let (mut fx, mut fy, mut tq) = (0.0, 0.0, 0.0);
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            let c = i + n * j;
            if !flow.liq.cell[c] {
                continue;
            }
            let (x, y) = (-half + (i as f64 + 0.5) * dx, -half + (j as f64 + 0.5) * dx);
            let d = ((x - disc.cx).powi(2) + (y - disc.cy).powi(2)).sqrt();
            if (d - disc.r).abs() > 1.5 * dx {
                continue;
            }
            let (ax, ay) = (
                (sm.au[c + 1] - sm.au[c]) * dx,
                (sm.av[c + n] - sm.av[c]) * dx,
            );
            let grad = |lo: usize, hi: usize| -> f64 {
                match (flow.liq.cell[lo], flow.liq.cell[hi]) {
                    (true, true) => (flow.p[hi] - flow.p[lo]) / (2.0 * dx),
                    (true, false) => (flow.p[c] - flow.p[lo]) / dx,
                    (false, true) => (flow.p[hi] - flow.p[c]) / dx,
                    _ => 0.0,
                }
            };
            let (gx, gy) = (grad(c - 1, c + 1), grad(c - n, c + n));
            let (sx, sy) = (
                disc.cx + disc.r * (x - disc.cx) / d,
                disc.cy + disc.r * (y - disc.cy) / d,
            );
            let pw = flow.p[c] + gx * (sx - x) + gy * (sy - y);
            let (px, py) = (-pw * ax, -pw * ay);
            fx += px;
            fy += py;
            tq += (sx - disc.cx) * py - (sy - disc.cy) * px;
        }
    }
    // Viscous part from the relative velocity.
    let mut ur = flow.u.clone();
    let mut vr = flow.v.clone();
    for j in 0..n {
        for i in 0..n {
            let c = i + n * j;
            let (xu, yu) = sm.gu.centre(i, j);
            ur[c] = flow.u_star[c] - disc.velocity(xu, yu).0;
            let (xv, yv) = sm.gv.centre(i, j);
            vr[c] = flow.v_star[c] - disc.velocity(xv, yv).1;
        }
    }
    let zero = vec![0.0; n * n];
    let tag = |x: f64, y: f64| usize::from(disc.sdf(x, y).abs() < 0.75 * dx);
    let lu = sm.gu.wall_load_masked(
        &ur,
        &zero,
        flow.nu,
        |_, _| (0.0, 0.0),
        tag,
        2,
        &flow.liq.active_u,
    );
    let lv = sm.gv.wall_load_masked(
        &zero,
        &vr,
        flow.nu,
        |_, _| (0.0, 0.0),
        tag,
        2,
        &flow.liq.active_v,
    );
    // Loads are on the fluid; the disc feels the opposite. Torque about the origin -> centre.
    let (lfx, lfy) = (lu[1].fx + lv[1].fx, lu[1].fy + lv[1].fy);
    let lt_origin = lu[1].torque + lv[1].torque;
    let lt_centre = lt_origin - (disc.cx * lfy - disc.cy * lfx);
    (fx - lfx, fy - lfy, tq - lt_centre)
}

/// Result of advancing the fluid one step around a disc.
#[derive(Clone, Copy, Debug)]
pub struct BodyLoad {
    /// Force of the fluid on the disc and its torque about the centre (per unit density).
    pub fx: f64,
    pub fy: f64,
    pub torque: f64,
    /// Largest cell divergence right after the projection.
    pub divergence: f64,
}

/// A rigid disc with inertia. `mass` and `inertia` are per unit fluid density (a disc of the
/// fluid's own density has `mass = pi r^2`).
#[derive(Clone, Copy, Debug)]
pub struct Body {
    pub disc: Disc,
    pub mass: f64,
    pub inertia: f64,
    /// Acceleration of the previous step, used to predict the position of the next one.
    pub accel: (f64, f64),
}

/// Outcome of a coupled step.
#[derive(Clone, Copy, Debug)]
pub struct CoupledInfo {
    pub iterations: usize,
    pub residual: f64,
    pub load: BodyLoad,
}

/// Fluid inside the drum of radius `drum_radius` around one moving disc: owns the flow state and
/// the previous mesh (needed to initialise nodes that the disc uncovers).
pub struct BodyFlow {
    pub n: usize,
    pub half: f64,
    pub drum_radius: f64,
    pub nu: f64,
    state: Option<FlowState>,
    prev_mesh: Option<StaggeredMesh>,
    /// Advection scheme of the trial steps.
    pub weno: bool,
    /// Blend the force along the line of centres with the sub-grid squeeze film near the wall.
    pub lubrication: bool,
}

impl BodyFlow {
    pub fn new(n: usize, half: f64, drum_radius: f64, nu: f64) -> Self {
        Self {
            n,
            half,
            drum_radius,
            nu,
            state: None,
            prev_mesh: None,
            weno: true,
            lubrication: true,
        }
    }

    pub fn time(&self) -> f64 {
        self.state.as_ref().map_or(0.0, |s| s.time)
    }

    /// Advances a copy of the flow by `dt` with `disc` (position and velocity at the *end* of the
    /// step); returns the load on the disc and the new state, which `commit` makes current.
    pub fn trial(&self, disc: &Disc, dt: f64) -> (BodyLoad, FlowState, StaggeredMesh) {
        let (n, half, radius) = (self.n, self.half, self.drum_radius);
        let d = *disc;
        let mut sm = StaggeredMesh::new(n, half, move |x, y| {
            ((x * x + y * y).sqrt() - radius).max(-d.sdf(x, y))
        });
        sm.set_body_flux(&[d]);
        let dx = sm.dx();
        let wall = move |x: f64, y: f64| {
            if d.sdf(x, y).abs() < 0.6 * dx {
                d.velocity(x, y)
            } else {
                (0.0, 0.0)
            }
        };
        let mut flow = match self.state.clone() {
            Some(s) => StaggeredFlow::from_state(&sm, s, wall),
            None => {
                let mut f = StaggeredFlow::new(&sm, self.nu, wall);
                f.weno = self.weno;
                f.upwind = self.weno;
                f
            }
        };
        if let Some(old) = &self.prev_mesh {
            flow.initialise_fresh(old, move |x, y| d.velocity(x, y));
        }
        flow.step(dt);
        let (fx, fy, torque) = disc_load(&flow, disc);
        let (fx, fy) = if self.lubrication {
            super::lubrication::wall_adjusted_load(
                (fx, fy),
                (disc.cx, disc.cy),
                (disc.ux, disc.uy),
                disc.r,
                radius,
                dx,
                self.nu,
            )
        } else {
            (fx, fy)
        };
        let load = BodyLoad {
            fx,
            fy,
            torque,
            divergence: flow.projected_divergence,
        };
        let state = flow.into_state();
        (load, state, sm)
    }

    /// Makes the result of a `trial` the current state.
    pub fn commit(&mut self, state: FlowState, mesh: StaggeredMesh) {
        self.state = Some(state);
        self.prev_mesh = Some(mesh);
    }

    /// Advances fluid and disc together over `dt`. The disc equation is integrated with backward
    /// Euler, `m (u' - u) / dt = F_fluid(u') + F_ext`, and solved by fixed-point iteration on the
    /// new velocity; the iteration is preconditioned with the linear part of the fluid force,
    /// `-m_a du/dt - k u` (added mass `added_mass`, viscous drag per unit speed `drag_stiffness`),
    /// without which plain iteration diverges for light bodies or strong viscosity. Rotation is
    /// preconditioned with the torque stiffness `spin_stiffness` (`T ~ -kappa omega`).
    #[allow(clippy::too_many_arguments)]
    pub fn step_coupled(
        &mut self,
        body: &mut Body,
        external: (f64, f64, f64),
        dt: f64,
        added_mass: f64,
        drag_stiffness: f64,
        spin_stiffness: f64,
        tol: f64,
        max_iter: usize,
    ) -> CoupledInfo {
        let (cx, cy) = (body.disc.cx, body.disc.cy);
        let (u0, v0, w0) = (body.disc.ux, body.disc.uy, body.disc.omega);
        let (mut uk, mut vk, mut wk) = (u0, v0, w0);
        let mut info = None;
        let k0 = added_mass / dt + drag_stiffness;
        // Stiffness of the fluid force per component, refined by secant estimates between
        // successive iterates (the response of the cut-cell force to the body velocity changes with
        // the mesh and the history, so a fixed preconditioner can lose contraction).
        let (mut kx, mut ky) = (k0, k0);
        let mut kw = spin_stiffness;
        let mut last: Option<(f64, f64, f64, f64, f64, f64)> = None;
        for it in 1..=max_iter {
            let mut d = body.disc;
            d.ux = uk;
            d.uy = vk;
            d.omega = wk;
            // The position is predicted once (constant acceleration) and kept fixed during the
            // iteration: with it moving, the force would jump whenever the surface crosses a
            // node and the fixed-point map would no longer be contractive.
            d.cx = cx + dt * (u0 + 0.5 * dt * body.accel.0);
            d.cy = cy + dt * (v0 + 0.5 * dt * body.accel.1);
            let (load, state, mesh) = self.trial(&d, dt);
            let m = body.mass;
            let (fxk, fyk) = (load.fx + external.0, load.fy + external.1);
            let tk = load.torque + external.2;
            if let Some((pu, pv, pfx, pfy, pw, pt)) = last {
                let secant = |df: f64, du: f64, fallback: f64| {
                    let k = -df / du;
                    if du.abs() > 1e-9 * (1.0 + uk.abs()) && k.is_finite() && k > 0.0 {
                        k.clamp(0.1 * k0, 100.0 * k0)
                    } else {
                        fallback
                    }
                };
                kx = secant(fxk - pfx, uk - pu, kx);
                ky = secant(fyk - pfy, vk - pv, ky);
                let dw = wk - pw;
                let kk = -(tk - pt) / dw;
                if dw.abs() > 1e-9 * (1.0 + wk.abs()) && kk.is_finite() && kk > 0.0 {
                    kw = kk.clamp(0.1 * spin_stiffness, 100.0 * spin_stiffness.max(1e-12));
                }
            }
            last = Some((uk, vk, fxk, fyk, wk, tk));
            let un = (m * u0 / dt + fxk + kx * uk) / (m / dt + kx);
            let vn = (m * v0 / dt + fyk + ky * vk) / (m / dt + ky);
            let wn = (body.inertia * w0 / dt + tk + kw * wk) / (body.inertia / dt + kw);
            let rr = body.disc.r;
            let scale = ((un - u0).powi(2) + (vn - v0).powi(2) + (rr * (wn - w0)).powi(2))
                .sqrt()
                .max(1e-12);
            let residual =
                ((un - uk).powi(2) + (vn - vk).powi(2) + (rr * (wn - wk)).powi(2)).sqrt() / scale;
            uk = un;
            vk = vn;
            wk = wn;
            info = Some((it, residual, load, state, mesh, d));
            if residual < tol {
                break;
            }
        }
        let (iterations, residual, load, state, mesh, d) = info.expect("at least one iteration");
        // The committed trial used the previous iterate of the velocity; the position it moved
        // to is the one the disc takes.
        self.commit(state, mesh);
        body.accel = ((uk - u0) / dt, (vk - v0) / dt);
        body.disc = Disc {
            ux: uk,
            uy: vk,
            omega: wk,
            ..d
        };
        CoupledInfo {
            iterations,
            residual,
            load,
        }
    }
}

/// Per-body linear preconditioner data of the many-body coupling.
#[derive(Clone, Copy, Debug)]
pub struct BodyStiffness {
    pub added_mass: f64,
    pub drag_stiffness: f64,
    pub spin_stiffness: f64,
}

/// Outcome of a many-body coupled step.
#[derive(Clone, Debug)]
pub struct CoupledManyInfo {
    pub iterations: usize,
    pub residual: f64,
    pub loads: Vec<BodyLoad>,
    /// Number of lubricated contacts (disc-disc and disc-wall) in the step.
    pub links: usize,
}

/// Nearest-surface velocity of a set of discs at a point (used for wall values and fresh nodes):
/// distance to the nearest surface and the rigid-body velocity there.
fn nearest_velocity(discs: &[Disc], x: f64, y: f64) -> (f64, f64, f64) {
    let mut best = (f64::INFINITY, (0.0, 0.0));
    for d in discs {
        let s = d.sdf(x, y).abs();
        if s < best.0 {
            best = (s, d.velocity(x, y));
        }
    }
    (best.0, best.1 .0, best.1 .1)
}

/// Solves `(D_i + sum_links w c n n^T) u_i - sum_links w c n n^T u_j = base_i` for all discs by
/// Gauss-Seidel (the system is symmetric and diagonally dominant).
fn solve_lubricated(
    base: &[(f64, f64)],
    diag: &[(f64, f64)],
    links: &[super::lubrication::Link],
    guess: &[(f64, f64)],
) -> Vec<(f64, f64)> {
    let n = base.len();
    let mut a: Vec<[f64; 4]> = diag.iter().map(|d| [d.0, 0.0, 0.0, d.1]).collect();
    let mut nbr: Vec<Vec<(usize, [f64; 4])>> = vec![Vec::new(); n];
    for l in links {
        let w = l.weight * l.c;
        let nn = [
            w * l.n.0 * l.n.0,
            w * l.n.0 * l.n.1,
            w * l.n.0 * l.n.1,
            w * l.n.1 * l.n.1,
        ];
        for (k, v) in nn.iter().enumerate() {
            a[l.i][k] += v;
        }
        if let Some(j) = l.j {
            for (k, v) in nn.iter().enumerate() {
                a[j][k] += v;
            }
            nbr[l.i].push((j, nn));
            nbr[j].push((l.i, nn));
        }
    }
    let mut u = guess.to_vec();
    for _ in 0..2000 {
        let mut change = 0.0f64;
        let mut size = 1e-30f64;
        for i in 0..n {
            let mut r = base[i];
            for (j, nn) in &nbr[i] {
                r.0 += nn[0] * u[*j].0 + nn[1] * u[*j].1;
                r.1 += nn[2] * u[*j].0 + nn[3] * u[*j].1;
            }
            let m = a[i];
            let det = m[0] * m[3] - m[1] * m[2];
            let new = (
                (m[3] * r.0 - m[1] * r.1) / det,
                (m[0] * r.1 - m[2] * r.0) / det,
            );
            change = change.max((new.0 - u[i].0).abs().max((new.1 - u[i].1).abs()));
            size = size.max(new.0.abs().max(new.1.abs()));
            u[i] = new;
        }
        if change < 1e-13 * size {
            break;
        }
    }
    u
}

impl BodyFlow {
    /// Advances a copy of the flow by `dt` with several discs (positions and velocities at the end
    /// of the step); returns the raw resolved loads (force, torque about the centre), the largest
    /// post-projection divergence and the new state. No lubrication is applied here.
    pub fn trial_many(
        &self,
        discs: &[Disc],
        dt: f64,
    ) -> (Vec<(f64, f64, f64)>, f64, FlowState, StaggeredMesh) {
        let (n, half, radius) = (self.n, self.half, self.drum_radius);
        let ds: Vec<Disc> = discs.to_vec();
        let ds_sdf = ds.clone();
        let mut sm = StaggeredMesh::new(n, half, move |x, y| {
            let inner = ds_sdf
                .iter()
                .map(|d| d.sdf(x, y))
                .fold(f64::INFINITY, f64::min);
            ((x * x + y * y).sqrt() - radius).max(-inner)
        });
        sm.set_body_flux(&ds);
        let dx = sm.dx();
        let ds_wall = ds.clone();
        let wall = move |x: f64, y: f64| {
            let (s, ux, uy) = nearest_velocity(&ds_wall, x, y);
            if s < 0.6 * dx {
                (ux, uy)
            } else {
                (0.0, 0.0)
            }
        };
        let mut flow = match self.state.clone() {
            Some(s) => StaggeredFlow::from_state(&sm, s, wall),
            None => {
                let mut f = StaggeredFlow::new(&sm, self.nu, wall);
                f.weno = self.weno;
                f.upwind = self.weno;
                f
            }
        };
        if let Some(old) = &self.prev_mesh {
            let ds_fresh = ds.clone();
            flow.initialise_fresh(old, move |x, y| {
                let (_, ux, uy) = nearest_velocity(&ds_fresh, x, y);
                (ux, uy)
            });
        }
        flow.step(dt);
        let loads: Vec<(f64, f64, f64)> = discs.iter().map(|d| disc_load(&flow, d)).collect();
        let divergence = flow.projected_divergence;
        let state = flow.into_state();
        (loads, divergence, state, sm)
    }

    /// Advances fluid and several discs together over `dt`: like `step_coupled`, a per-body
    /// secant-preconditioned fixed point on the new velocities, plus the sub-grid squeeze film of
    /// every close pair and of every disc near the wall, which is stiff and coupled between
    /// bodies and therefore solved implicitly (Gauss-Seidel) inside each iteration.
    pub fn step_coupled_many(
        &mut self,
        bodies: &mut [Body],
        external: &[(f64, f64, f64)],
        setup: &[BodyStiffness],
        dt: f64,
        tol: f64,
        max_iter: usize,
    ) -> CoupledManyInfo {
        let nb = bodies.len();
        let dx = 2.0 * self.half / self.n as f64;
        let start: Vec<(f64, f64, f64)> = bodies
            .iter()
            .map(|b| (b.disc.ux, b.disc.uy, b.disc.omega))
            .collect();
        let pos: Vec<(f64, f64)> = bodies
            .iter()
            .map(|b| {
                (
                    b.disc.cx + dt * (b.disc.ux + 0.5 * dt * b.accel.0),
                    b.disc.cy + dt * (b.disc.uy + 0.5 * dt * b.accel.1),
                )
            })
            .collect();
        let radii: Vec<f64> = bodies.iter().map(|b| b.disc.r).collect();
        let lk = if self.lubrication {
            super::lubrication::links(&pos, &radii, self.drum_radius, dx, self.nu)
        } else {
            Vec::new()
        };
        let mut vel: Vec<(f64, f64, f64)> = start.clone();
        let mut kx: Vec<f64> = setup
            .iter()
            .map(|s| s.added_mass / dt + s.drag_stiffness)
            .collect();
        let mut ky = kx.clone();
        let k0 = kx.clone();
        let mut kw: Vec<f64> = setup.iter().map(|s| s.spin_stiffness).collect();
        type Iterate = (Vec<(f64, f64, f64)>, Vec<(f64, f64, f64)>);
        let mut last: Option<Iterate> = None;
        let mut result = None;
        for it in 1..=max_iter {
            let discs: Vec<Disc> = (0..nb)
                .map(|i| {
                    let mut d = bodies[i].disc;
                    d.ux = vel[i].0;
                    d.uy = vel[i].1;
                    d.omega = vel[i].2;
                    d.cx = pos[i].0;
                    d.cy = pos[i].1;
                    d
                })
                .collect();
            let (raw, div, state, mesh) = self.trial_many(&discs, dt);
            let forces: Vec<(f64, f64)> = raw.iter().map(|l| (l.0, l.1)).collect();
            let total: Vec<(f64, f64, f64)> = (0..nb)
                .map(|i| {
                    (
                        forces[i].0 + external[i].0,
                        forces[i].1 + external[i].1,
                        raw[i].2 + external[i].2,
                    )
                })
                .collect();
            if let Some((pvel, ptot)) = &last {
                let sec = |df: f64, du: f64, fb: f64, k0: f64, u: f64| {
                    let k = -df / du;
                    if du.abs() > 1e-9 * (1.0 + u.abs()) && k.is_finite() && k > 0.0 {
                        k.clamp(0.1 * k0, 100.0 * k0)
                    } else {
                        fb
                    }
                };
                for i in 0..nb {
                    kx[i] = sec(
                        total[i].0 - ptot[i].0,
                        vel[i].0 - pvel[i].0,
                        kx[i],
                        k0[i],
                        vel[i].0,
                    );
                    ky[i] = sec(
                        total[i].1 - ptot[i].1,
                        vel[i].1 - pvel[i].1,
                        ky[i],
                        k0[i],
                        vel[i].1,
                    );
                    let dw = vel[i].2 - pvel[i].2;
                    let kk = -(total[i].2 - ptot[i].2) / dw;
                    if dw.abs() > 1e-9 * (1.0 + vel[i].2.abs()) && kk.is_finite() && kk > 0.0 {
                        kw[i] = kk.clamp(
                            0.1 * setup[i].spin_stiffness,
                            100.0 * setup[i].spin_stiffness.max(1e-12),
                        );
                    }
                }
            }
            last = Some((vel.clone(), total.clone()));
            let base: Vec<(f64, f64)> = (0..nb)
                .map(|i| {
                    let m = bodies[i].mass;
                    (
                        m * start[i].0 / dt + total[i].0 + kx[i] * vel[i].0,
                        m * start[i].1 / dt + total[i].1 + ky[i] * vel[i].1,
                    )
                })
                .collect();
            let diag: Vec<(f64, f64)> = (0..nb)
                .map(|i| (bodies[i].mass / dt + kx[i], bodies[i].mass / dt + ky[i]))
                .collect();
            let guess: Vec<(f64, f64)> = vel.iter().map(|v| (v.0, v.1)).collect();
            let new_t = solve_lubricated(&base, &diag, &lk, &guess);
            let mut residual = 0.0f64;
            let mut new_vel = vel.clone();
            for i in 0..nb {
                let wn = (bodies[i].inertia * start[i].2 / dt + total[i].2 + kw[i] * vel[i].2)
                    / (bodies[i].inertia / dt + kw[i]);
                let rr = radii[i];
                let scale = ((new_t[i].0 - start[i].0).powi(2)
                    + (new_t[i].1 - start[i].1).powi(2)
                    + (rr * (wn - start[i].2)).powi(2))
                .sqrt()
                .max(1e-12);
                let r = ((new_t[i].0 - vel[i].0).powi(2)
                    + (new_t[i].1 - vel[i].1).powi(2)
                    + (rr * (wn - vel[i].2)).powi(2))
                .sqrt()
                    / scale;
                residual = residual.max(r);
                new_vel[i] = (new_t[i].0, new_t[i].1, wn);
            }
            vel = new_vel;
            let loads: Vec<BodyLoad> = (0..nb)
                .map(|i| BodyLoad {
                    fx: total[i].0 - external[i].0,
                    fy: total[i].1 - external[i].1,
                    torque: raw[i].2,
                    divergence: div,
                })
                .collect();
            result = Some((it, residual, loads, state, mesh));
            if residual < tol {
                break;
            }
        }
        let (iterations, residual, loads, state, mesh) = result.expect("at least one iteration");
        self.commit(state, mesh);
        for i in 0..nb {
            bodies[i].accel = ((vel[i].0 - start[i].0) / dt, (vel[i].1 - start[i].1) / dt);
            bodies[i].disc = Disc {
                cx: pos[i].0,
                cy: pos[i].1,
                ux: vel[i].0,
                uy: vel[i].1,
                omega: vel[i].2,
                ..bodies[i].disc
            };
        }
        CoupledManyInfo {
            iterations,
            residual,
            loads,
            links: lk.len(),
        }
    }
}
