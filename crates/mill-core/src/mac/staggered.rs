//! Staggered (MAC) variational cut-cell Navier-Stokes: the unknowns are the face-normal
//! velocities themselves (`u` on vertical faces, `v` on horizontal faces), so there is no
//! cell-velocity / face-flux pairing to be inconsistent at cut walls.
//!
//! * viscosity: Gibou symmetric Dirichlet Helmholtz on each face lattice (`viscous.rs`),
//! * pressure: aperture-weighted Neumann Poisson (E0) with an exact face projection,
//! * advection: central, other component by 4-face averaging, AB2 (accuracy first; a robust
//!   upwind/semi-Lagrangian scheme is a separate concern for the free-surface experiments),
//! * velocities at face nodes whose centre lies in the solid but whose face is partly open
//!   ("ghost" faces) and two further layers are filled by polynomial extrapolation after every
//!   stage.
//!
//! Storage: all fields are `n x n`, index `i + n j`. `u[i + n j]` sits at
//! `(-half + i dx, -half + (j + 1/2) dx)`, `v[i + n j]` at `(-half + (i + 1/2) dx, -half + j dx)`,
//! `p` at cell centres. Faces on the box border (`i = n` / `j = n`) are solid by contract.
//!
//! Units: density 1 (pressure is `p / rho`), kinematic viscosity `nu`.

use super::flow::Mesh;
use super::multigrid::PoissonSolver;
use super::viscous::{FluidGrid, Link};

/// Geometry of the staggered discretisation.
pub struct StaggeredMesh {
    pub mesh: Mesh,
    /// `u`-face lattice (offset `(0, 1/2)`) and `v`-face lattice (offset `(1/2, 0)`).
    pub gu: FluidGrid,
    pub gv: FluidGrid,
    /// 0 fluid node, 1..2 ghost layers (extrapolated), 255 unused.
    pub layer_u: Vec<u8>,
    pub layer_v: Vec<u8>,
    /// Aperture of the `u` face / `v` face at the same index.
    pub au: Vec<f64>,
    pub av: Vec<f64>,
    /// Offset of the open part's centroid from the face node along the face (metres); the face
    /// flux is sampled there: `u_node + off * du/dt` (second order at cut faces).
    pub off_u: Vec<f64>,
    pub off_v: Vec<f64>,
}

/// Centroid (as a fraction `t` of the segment `p0 -> p1`) of the part with `phi < 0`; 0.5 if none.
fn open_centroid(phi: &impl Fn(f64, f64) -> f64, p0: (f64, f64), p1: (f64, f64)) -> f64 {
    const M: usize = 16;
    let at = |t: f64| (p0.0 + t * (p1.0 - p0.0), p0.1 + t * (p1.1 - p0.1));
    let inside = |t: f64| {
        let (x, y) = at(t);
        phi(x, y) < 0.0
    };
    let (mut len, mut mom) = (0.0, 0.0);
    for k in 0..M {
        let (ta, tb) = (k as f64 / M as f64, (k + 1) as f64 / M as f64);
        let (ia, ib) = (inside(ta), inside(tb));
        let (a, b) = match (ia, ib) {
            (true, true) => (ta, tb),
            (false, false) => continue,
            _ => {
                let (mut lo, mut hi) = (ta, tb);
                for _ in 0..50 {
                    let mid = 0.5 * (lo + hi);
                    if inside(mid) == ia {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let t = 0.5 * (lo + hi);
                if ia {
                    (ta, t)
                } else {
                    (t, tb)
                }
            }
        };
        len += b - a;
        mom += 0.5 * (b * b - a * a);
    }
    if len > 0.0 {
        mom / len
    } else {
        0.5
    }
}

fn ghost_layers(fluid: &[bool], open: &[bool], n: usize) -> Vec<u8> {
    let mut layer: Vec<u8> = fluid.iter().map(|&f| if f { 0 } else { 255 }).collect();
    for depth in 1..=2u8 {
        let prev = layer.clone();
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if prev[c] != 255 {
                    continue;
                }
                let near = (i > 0 && prev[c - 1] == depth - 1)
                    || (i + 1 < n && prev[c + 1] == depth - 1)
                    || (j > 0 && prev[c - n] == depth - 1)
                    || (j + 1 < n && prev[c + n] == depth - 1);
                if near {
                    layer[c] = depth;
                }
            }
        }
    }
    for c in 0..n * n {
        if open[c] && layer[c] == 255 {
            layer[c] = 2;
        }
    }
    layer
}

impl StaggeredMesh {
    pub fn new(n: usize, half: f64, phi: impl Fn(f64, f64) -> f64) -> Self {
        let mesh = Mesh::new(n, half, &phi);
        let gu = FluidGrid::with_offset(n, half, 0.0, 0.5, &phi);
        let gv = FluidGrid::with_offset(n, half, 0.5, 0.0, &phi);
        let mut au = vec![0.0; n * n];
        let mut av = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                au[i + n * j] = mesh.wx[i + (n + 1) * j];
                av[i + n * j] = mesh.wy[i + n * j];
            }
        }
        let dx = mesh.grid.dx();
        let mut off_u = vec![0.0; n * n];
        let mut off_v = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                let (x0, y0) = (-half + i as f64 * dx, -half + j as f64 * dx);
                if au[c] > 0.0 {
                    off_u[c] = (open_centroid(&phi, (x0, y0), (x0, y0 + dx)) - 0.5) * dx;
                }
                if av[c] > 0.0 {
                    off_v[c] = (open_centroid(&phi, (x0, y0), (x0 + dx, y0)) - 0.5) * dx;
                }
            }
        }
        let open_u: Vec<bool> = au.iter().map(|&a| a > 0.0).collect();
        let open_v: Vec<bool> = av.iter().map(|&a| a > 0.0).collect();
        let layer_u = ghost_layers(&gu.fluid, &open_u, n);
        let layer_v = ghost_layers(&gv.fluid, &open_v, n);
        Self {
            mesh,
            gu,
            gv,
            layer_u,
            layer_v,
            au,
            av,
            off_u,
            off_v,
        }
    }

    pub fn n(&self) -> usize {
        self.mesh.n()
    }

    pub fn dx(&self) -> f64 {
        self.mesh.grid.dx()
    }
}

/// Fills the ghost layers of `q` by polynomial extrapolation from the interior along each grid
/// direction: quadratic `3 q1 - 3 q2 + q3` where three interior values exist, else linear
/// `2 q1 - q2`, else constant; directions of the best available order are averaged.
pub fn extend(q: &mut [f64], layer: &[u8], n: usize) {
    let ok = |a: isize, b: isize| a >= 0 && b >= 0 && (a as usize) < n && (b as usize) < n;
    for depth in 1..=2u8 {
        let snapshot = q.to_vec();
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if layer[c] != depth {
                    continue;
                }
                let at = |k: isize, di: isize, dj: isize| -> Option<f64> {
                    let (ii, jj) = (i as isize + k * di, j as isize + k * dj);
                    (ok(ii, jj) && layer[ii as usize + n * jj as usize] < depth)
                        .then(|| snapshot[ii as usize + n * jj as usize])
                };
                let mut best = (0u8, 0.0f64, 0.0f64);
                for (di, dj) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                    let cand = match (at(1, di, dj), at(2, di, dj), at(3, di, dj)) {
                        (Some(a), Some(b), Some(c3)) => Some((3u8, 3.0 * a - 3.0 * b + c3)),
                        (Some(a), Some(b), None) => Some((2, 2.0 * a - b)),
                        (Some(a), None, _) => Some((1, a)),
                        _ => None,
                    };
                    if let Some((order, val)) = cand {
                        if order > best.0 {
                            best = (order, val, 1.0);
                        } else if order == best.0 {
                            best.1 += val;
                            best.2 += 1.0;
                        }
                    }
                }
                if best.0 > 0 {
                    q[c] = best.1 / best.2;
                }
            }
        }
    }
}

type WallVelocity<'a> = Box<dyn Fn(f64, f64) -> (f64, f64) + 'a>;

struct Prev {
    u: Vec<f64>,
    v: Vec<f64>,
    au: Vec<f64>,
    av: Vec<f64>,
    dt: f64,
}

/// Flow state and integrator.
pub struct StaggeredFlow<'a> {
    pub sm: &'a StaggeredMesh,
    pub nu: f64,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub p: Vec<f64>,
    wall_velocity: WallVelocity<'a>,
    prev: Option<Prev>,
    poisson: PoissonSolver,
    pub time: f64,
    /// Largest `|div|` (per cell, aperture weighted) right after the last projection.
    pub projected_divergence: f64,
    /// Flux sampling corrections `off * d(velocity)/d(face direction)`, held fixed during a
    /// projection so that it stays exact.
    cu: Vec<f64>,
    cv: Vec<f64>,
}

impl<'a> StaggeredFlow<'a> {
    pub fn new(
        sm: &'a StaggeredMesh,
        nu: f64,
        wall_velocity: impl Fn(f64, f64) -> (f64, f64) + 'a,
    ) -> Self {
        let n = sm.n();
        let poisson = PoissonSolver::from_weights(
            n,
            sm.mesh.wx.clone(),
            sm.mesh.wy.clone(),
            vec![0.0; n * n],
        );
        Self {
            sm,
            nu,
            u: vec![0.0; n * n],
            v: vec![0.0; n * n],
            p: vec![0.0; n * n],
            wall_velocity: Box::new(wall_velocity),
            prev: None,
            poisson,
            time: 0.0,
            projected_divergence: 0.0,
            cu: vec![0.0; n * n],
            cv: vec![0.0; n * n],
        }
    }

    fn update_flux_corrections(&mut self) {
        let sm = self.sm;
        let n = sm.n();
        let inv = 1.0 / (2.0 * sm.dx());
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                self.cu[c] = sm.off_u[c] * (self.u[c + n] - self.u[c - n]) * inv;
                self.cv[c] = sm.off_v[c] * (self.v[c + 1] - self.v[c - 1]) * inv;
            }
        }
    }

    /// Sets the face velocities from `f(x, y)` on every open face and projects them.
    pub fn set_velocity(&mut self, f: impl Fn(f64, f64) -> (f64, f64)) {
        let sm = self.sm;
        let n = sm.n();
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if sm.layer_u[c] < 255 {
                    let (x, y) = sm.gu.centre(i, j);
                    self.u[c] = f(x, y).0;
                }
                if sm.layer_v[c] < 255 {
                    let (x, y) = sm.gv.centre(i, j);
                    self.v[c] = f(x, y).1;
                }
            }
        }
        self.project(1.0);
        self.p.fill(0.0);
        self.prev = None;
        self.init_pressure();
    }

    /// Aperture-weighted divergence of cell `(i, j)` (times `dx`).
    fn divergence(&self, i: usize, j: usize) -> f64 {
        let sm = self.sm;
        let n = sm.n();
        let fu = |c: usize| sm.au[c] * (self.u[c] + self.cu[c]);
        let fv = |c: usize| sm.av[c] * (self.v[c] + self.cv[c]);
        let right = if i + 1 < n { fu(i + 1 + n * j) } else { 0.0 };
        let top = if j + 1 < n { fv(i + n * (j + 1)) } else { 0.0 };
        right - fu(i + n * j) + top - fv(i + n * j)
    }

    pub fn max_divergence(&self) -> f64 {
        let n = self.sm.n();
        (0..n * n)
            .filter(|&c| self.sm.mesh.pressure_active[c])
            .map(|c| self.divergence(c % n, c / n).abs())
            .fold(0.0, f64::max)
    }

    /// Exact projection with the weighted operator; `dt_eff` is the pressure-update time scale.
    /// Returns the pressure increment; updates `u, v` on every open face (then re-extends).
    fn project(&mut self, dt_eff: f64) -> Vec<f64> {
        let sm = self.sm;
        let n = sm.n();
        let dx = sm.dx();
        self.update_flux_corrections();
        let mut rhs = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if sm.mesh.pressure_active[c] {
                    rhs[c] = -(dx / dt_eff) * self.divergence(i, j);
                }
            }
        }
        let mut phi = vec![0.0; n * n];
        self.poisson.solve(&rhs, &mut phi, 1e-12, 300);
        let k = dt_eff / dx;
        for j in 0..n {
            for i in 1..n {
                let c = i + n * j;
                if sm.au[c] > 0.0 {
                    self.u[c] -= k * (phi[c] - phi[c - 1]);
                }
            }
        }
        for j in 1..n {
            for i in 0..n {
                let c = i + n * j;
                if sm.av[c] > 0.0 {
                    self.v[c] -= k * (phi[c] - phi[c - n]);
                }
            }
        }
        self.projected_divergence = self.max_divergence();
        extend(&mut self.u, &sm.layer_u, n);
        extend(&mut self.v, &sm.layer_v, n);
        phi
    }

    /// Face gradient of a cell field, `(dp/dx at u nodes, dp/dy at v nodes)`.
    fn pressure_gradient(&self, p: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let sm = self.sm;
        let n = sm.n();
        let dx = sm.dx();
        let mut gx = vec![0.0; n * n];
        let mut gy = vec![0.0; n * n];
        for j in 1..n {
            for i in 1..n {
                let c = i + n * j;
                if sm.au[c] > 0.0 {
                    gx[c] = (p[c] - p[c - 1]) / dx;
                }
                if sm.av[c] > 0.0 {
                    gy[c] = (p[c] - p[c - n]) / dx;
                }
            }
        }
        (gx, gy)
    }

    /// Convective terms `(U . grad) u` at `u` nodes and `(U . grad) v` at `v` nodes, for fluid
    /// nodes (zero elsewhere); central differences on the extended fields.
    fn advection(&self, u: &[f64], v: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let sm = self.sm;
        let n = sm.n();
        let inv = 1.0 / (2.0 * sm.dx());
        let mut au = vec![0.0; n * n];
        let mut av = vec![0.0; n * n];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                if sm.gu.fluid[c] {
                    let vb = 0.25 * (v[c - 1] + v[c] + v[c - 1 + n] + v[c + n]);
                    au[c] = u[c] * (u[c + 1] - u[c - 1]) * inv + vb * (u[c + n] - u[c - n]) * inv;
                }
                if sm.gv.fluid[c] {
                    let ub = 0.25 * (u[c - n] + u[c + 1 - n] + u[c] + u[c + 1]);
                    av[c] = ub * (v[c + 1] - v[c - 1]) * inv + v[c] * (v[c + n] - v[c - n]) * inv;
                }
            }
        }
        (au, av)
    }

    /// Lattice Laplacian with Dirichlet wall values (Gibou links) at fluid nodes of `g`.
    fn laplacian(g: &FluidGrid, q: &[f64], bval: impl Fn(f64, f64) -> f64) -> Vec<f64> {
        let dx2 = g.dx() * g.dx();
        let mut out = vec![0.0; q.len()];
        for c in 0..q.len() {
            if !g.fluid[c] {
                continue;
            }
            let mut s = 0.0;
            for link in &g.links[c] {
                match *link {
                    Link::Fluid(nb) => s += q[nb] - q[c],
                    Link::Boundary { theta, bx, by } => s += (bval(bx, by) - q[c]) / theta,
                }
            }
            out[c] = s / dx2;
        }
        out
    }

    /// Pressure that balances the current acceleration `-(U.grad)u + nu lap u` (zero for a
    /// stationary state): the gradient part of the acceleration.
    pub fn init_pressure(&mut self) {
        let sm = self.sm;
        let n = sm.n();
        let wv = &self.wall_velocity;
        let lu = Self::laplacian(&sm.gu, &self.u, |x, y| wv(x, y).0);
        let lv = Self::laplacian(&sm.gv, &self.v, |x, y| wv(x, y).1);
        let (au, av) = self.advection(&self.u, &self.v);
        let (mut ax, mut ay) = (vec![0.0; n * n], vec![0.0; n * n]);
        for c in 0..n * n {
            if sm.gu.fluid[c] {
                ax[c] = -au[c] + self.nu * lu[c];
            }
            if sm.gv.fluid[c] {
                ay[c] = -av[c] + self.nu * lv[c];
            }
        }
        extend(&mut ax, &sm.layer_u, n);
        extend(&mut ay, &sm.layer_v, n);
        let u0 = std::mem::replace(&mut self.u, ax);
        let v0 = std::mem::replace(&mut self.v, ay);
        self.p = self.project(1.0);
        self.u = u0;
        self.v = v0;
    }

    /// One time step (BDF2 + AB2, incremental pressure correction; first step BDF1 + AB1).
    pub fn step(&mut self, dt: f64) {
        let sm = self.sm;
        let n = sm.n();
        let nn = n * n;
        let (au, av) = self.advection(&self.u, &self.v);
        let (gpx, gpy) = self.pressure_gradient(&self.p);
        let bdf2 = matches!(&self.prev, Some(pv) if (pv.dt - dt).abs() < 1e-14 * dt);
        let sigma = if bdf2 { 1.5 } else { 1.0 } / (self.nu * dt);
        let (mut su, mut sv) = (vec![0.0; nn], vec![0.0; nn]);
        let (mut xu, mut xv) = (vec![0.0; nn], vec![0.0; nn]);
        for c in 0..nn {
            let (eu, ev, pu, pv) = match (&self.prev, bdf2) {
                (Some(pv), true) => (
                    2.0 * au[c] - pv.au[c],
                    2.0 * av[c] - pv.av[c],
                    (4.0 * self.u[c] - pv.u[c]) / 3.0,
                    (4.0 * self.v[c] - pv.v[c]) / 3.0,
                ),
                _ => (au[c], av[c], self.u[c], self.v[c]),
            };
            su[c] = (-eu - gpx[c]) / self.nu;
            sv[c] = (-ev - gpy[c]) / self.nu;
            xu[c] = pu;
            xv[c] = pv;
        }
        let old = Prev {
            u: self.u.clone(),
            v: self.v.clone(),
            au,
            av,
            dt,
        };
        let wv = &self.wall_velocity;
        sm.gu
            .helmholtz(sigma)
            .solve(&mut xu, Some(&su), |x, y| wv(x, y).0, 1e-12);
        sm.gv
            .helmholtz(sigma)
            .solve(&mut xv, Some(&sv), |x, y| wv(x, y).1, 1e-12);
        extend(&mut xu, &sm.layer_u, n);
        extend(&mut xv, &sm.layer_v, n);
        self.u = xu;
        self.v = xv;
        let dt_eff = if bdf2 { dt * 2.0 / 3.0 } else { dt };
        let phi = self.project(dt_eff);
        for (p, f) in self.p.iter_mut().zip(&phi) {
            *p += f;
        }
        self.prev = Some(old);
        self.time += dt;
    }

    /// Wall load on the fluid summed over links by tag, from the two face lattices.
    pub fn wall_load(
        &self,
        tag: impl Fn(f64, f64) -> usize + Copy,
        tags: usize,
    ) -> Vec<super::viscous::WallLoad> {
        let sm = self.sm;
        let n = sm.n();
        let wv = &self.wall_velocity;
        let zero = vec![0.0; n * n];
        let lu = sm
            .gu
            .wall_load(&self.u, &zero, self.nu, |x, y| (wv(x, y).0, 0.0), tag, tags);
        let lv = sm
            .gv
            .wall_load(&zero, &self.v, self.nu, |x, y| (0.0, wv(x, y).1), tag, tags);
        lu.iter()
            .zip(&lv)
            .map(|(a, b)| super::viscous::WallLoad {
                fx: a.fx,
                fy: b.fy,
                torque: a.torque + b.torque,
            })
            .collect()
    }
}
