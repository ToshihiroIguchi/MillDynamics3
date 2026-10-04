//! Incompressible Navier-Stokes on the cut-cell grid: cell-centred velocity, face fluxes for
//! advection and an exact face projection with aperture weights (variational pressure), implicit
//! viscosity (`viscous.rs`), IMEX-BDF2 / AB2 time stepping with incremental pressure correction.
//!
//! Units: density 1 (pressure is `p / rho`), kinematic viscosity `nu`.

use super::multigrid::PoissonSolver;
use super::viscous::FluidGrid;

/// Aperture (fraction of the face segment inside the fluid) along the segment `p0 -> p1`, by
/// sampling `phi < 0` with bisection refinement of the crossings.
pub(super) fn aperture(phi: &impl Fn(f64, f64) -> f64, p0: (f64, f64), p1: (f64, f64)) -> f64 {
    const M: usize = 16;
    let at = |t: f64| (p0.0 + t * (p1.0 - p0.0), p0.1 + t * (p1.1 - p0.1));
    let inside = |t: f64| {
        let (x, y) = at(t);
        phi(x, y) < 0.0
    };
    // A segment farther from the boundary than its own length lies entirely on one side
    // (phi is a distance-like, 1-Lipschitz field).
    let seg = (p1.0 - p0.0).hypot(p1.1 - p0.1);
    let (mx, my) = at(0.5);
    let pm = phi(mx, my);
    if pm.abs() > seg {
        return if pm < 0.0 { 1.0 } else { 0.0 };
    }
    let mut total = 0.0;
    for k in 0..M {
        let (ta, tb) = (k as f64 / M as f64, (k + 1) as f64 / M as f64);
        let (ia, ib) = (inside(ta), inside(tb));
        total += match (ia, ib) {
            (true, true) => tb - ta,
            (false, false) => 0.0,
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
                    t - ta
                } else {
                    tb - t
                }
            }
        };
    }
    total
}

pub struct Mesh {
    pub grid: FluidGrid,
    /// Apertures of vertical faces `(n+1) x n` and horizontal faces `n x (n+1)`.
    pub wx: Vec<f64>,
    pub wy: Vec<f64>,
    /// Cells that carry a pressure unknown (some open face).
    pub pressure_active: Vec<bool>,
    /// Grid distance (4-neighbour steps) from the fluid-centre cells: 0 fluid, 1 and 2 ghost
    /// layers filled by extrapolation, 255 beyond.
    pub layer: Vec<u8>,
}

impl Mesh {
    pub fn new(n: usize, half: f64, phi: impl Fn(f64, f64) -> f64) -> Self {
        let grid = FluidGrid::new(n, half, &phi);
        let dx = grid.dx();
        let mut wx = vec![0.0; (n + 1) * n];
        let mut wy = vec![0.0; n * (n + 1)];
        for j in 0..n {
            let (y0, y1) = (-half + j as f64 * dx, -half + (j + 1) as f64 * dx);
            for i in 1..n {
                let x = -half + i as f64 * dx;
                wx[i + (n + 1) * j] = aperture(&phi, (x, y0), (x, y1));
            }
        }
        for j in 1..n {
            let y = -half + j as f64 * dx;
            for i in 0..n {
                let (x0, x1) = (-half + i as f64 * dx, -half + (i + 1) as f64 * dx);
                wy[i + n * j] = aperture(&phi, (x0, y), (x1, y));
            }
        }
        let mut pressure_active = vec![false; n * n];
        for j in 0..n {
            for i in 0..n {
                let s = wx[i + (n + 1) * j]
                    + wx[i + 1 + (n + 1) * j]
                    + wy[i + n * j]
                    + wy[i + n * (j + 1)];
                pressure_active[i + n * j] = s > 1e-12;
            }
        }
        let mut layer: Vec<u8> = grid
            .fluid
            .iter()
            .map(|&f| if f { 0 } else { 255 })
            .collect();
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
            if pressure_active[c] && layer[c] == 255 {
                layer[c] = 2;
            }
        }
        Self {
            grid,
            wx,
            wy,
            pressure_active,
            layer,
        }
    }

    pub fn n(&self) -> usize {
        self.grid.n
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
pub struct Flow<'a> {
    pub mesh: &'a Mesh,
    pub nu: f64,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    /// Face fluxes (normal velocity): vertical faces (x-velocity), horizontal faces (y-velocity).
    pub fu: Vec<f64>,
    pub fv: Vec<f64>,
    pub p: Vec<f64>,
    wall_velocity: WallVelocity<'a>,
    prev: Option<Prev>,
    poisson: PoissonSolver,
    pub time: f64,
    /// Diagnostics only: skip the pressure projection in `step`.
    pub skip_projection: bool,
    /// Diagnostics only: drop the convective term in `step`.
    pub skip_advection: bool,
}

impl<'a> Flow<'a> {
    pub fn new(
        mesh: &'a Mesh,
        nu: f64,
        wall_velocity: impl Fn(f64, f64) -> (f64, f64) + 'a,
    ) -> Self {
        let n = mesh.n();
        let poisson =
            PoissonSolver::from_weights(n, mesh.wx.clone(), mesh.wy.clone(), vec![0.0; n * n]);
        Self {
            mesh,
            nu,
            u: vec![0.0; n * n],
            v: vec![0.0; n * n],
            fu: vec![0.0; (n + 1) * n],
            fv: vec![0.0; n * (n + 1)],
            p: vec![0.0; n * n],
            wall_velocity: Box::new(wall_velocity),
            prev: None,
            poisson,
            time: 0.0,
            skip_projection: false,
            skip_advection: false,
        }
    }

    /// Sets the cell velocity from `f(x, y)` and makes the face fluxes divergence free.
    pub fn set_velocity(&mut self, f: impl Fn(f64, f64) -> (f64, f64)) {
        let g = &self.mesh.grid;
        for j in 0..g.n {
            for i in 0..g.n {
                let c = i + g.n * j;
                if self.mesh.pressure_active[c] {
                    let (x, y) = g.centre(i, j);
                    (self.u[c], self.v[c]) = f(x, y);
                }
            }
        }
        self.faces_from_cells();
        self.project(1.0);
        self.p.fill(0.0);
        self.prev = None;
        self.init_pressure();
    }

    pub fn faces_from_cells(&mut self) {
        let m = self.mesh;
        let n = m.n();
        for j in 0..n {
            for i in 1..n {
                let f = i + (n + 1) * j;
                self.fu[f] = if m.wx[f] > 0.0 {
                    0.5 * (self.u[i - 1 + n * j] + self.u[i + n * j])
                } else {
                    0.0
                };
            }
        }
        for j in 1..n {
            for i in 0..n {
                let f = i + n * j;
                self.fv[f] = if m.wy[f] > 0.0 {
                    0.5 * (self.v[i + n * (j - 1)] + self.v[i + n * j])
                } else {
                    0.0
                };
            }
        }
    }

    /// Fills the two ghost layers around the fluid-centre cells by polynomial extrapolation from
    /// the interior along each grid direction: quadratic `3 q1 - 3 q2 + q3` where three interior
    /// values exist, else linear `2 q1 - q2`, else constant. Directions of the best available
    /// order are averaged. (A quadratic ghost keeps the central first derivative second order.)
    pub fn extend(&self, q: &mut [f64]) {
        let m = self.mesh;
        let n = m.n();
        let ok = |a: isize, b: isize| a >= 0 && b >= 0 && (a as usize) < n && (b as usize) < n;
        for depth in 1..=2u8 {
            let snapshot = q.to_vec();
            for j in 0..n {
                for i in 0..n {
                    let c = i + n * j;
                    if m.layer[c] != depth {
                        continue;
                    }
                    let at = |k: isize, di: isize, dj: isize| -> Option<f64> {
                        let (ii, jj) = (i as isize + k * di, j as isize + k * dj);
                        (ok(ii, jj) && m.layer[ii as usize + n * jj as usize] < depth)
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

    pub fn divergence(&self, i: usize, j: usize) -> f64 {
        let m = self.mesh;
        let n = m.n();
        m.wx[i + 1 + (n + 1) * j] * self.fu[i + 1 + (n + 1) * j]
            - m.wx[i + (n + 1) * j] * self.fu[i + (n + 1) * j]
            + m.wy[i + n * (j + 1)] * self.fv[i + n * (j + 1)]
            - m.wy[i + n * j] * self.fv[i + n * j]
    }

    /// Exact face projection with the weighted operator; `dt_eff` is the pressure-update time
    /// scale. Returns `phi` (pressure increment); updates `fu, fv, u, v`.
    pub fn project(&mut self, dt_eff: f64) -> Vec<f64> {
        let m = self.mesh;
        let n = m.n();
        let dx = m.grid.dx();
        let mut rhs = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if m.pressure_active[c] {
                    rhs[c] = -(dx / dt_eff) * self.divergence(i, j);
                }
            }
        }
        let mut phi = vec![0.0; n * n];
        self.poisson.solve(&rhs, &mut phi, 1e-12, 300);
        let k = dt_eff / dx;
        for j in 0..n {
            for i in 1..n {
                let f = i + (n + 1) * j;
                if m.wx[f] > 0.0 {
                    self.fu[f] -= k * (phi[i + n * j] - phi[i - 1 + n * j]);
                }
            }
        }
        for j in 1..n {
            for i in 0..n {
                let f = i + n * j;
                if m.wy[f] > 0.0 {
                    self.fv[f] -= k * (phi[i + n * j] - phi[i + n * (j - 1)]);
                }
            }
        }
        let (gx, gy) = self.cell_gradient(&phi);
        for c in 0..n * n {
            if m.grid.fluid[c] {
                self.u[c] -= dt_eff * gx[c];
                self.v[c] -= dt_eff * gy[c];
            }
        }
        let (mut u, mut v) = (std::mem::take(&mut self.u), std::mem::take(&mut self.v));
        self.extend(&mut u);
        self.extend(&mut v);
        self.u = u;
        self.v = v;
        phi
    }

    /// Cell-centre gradient of a cell field: central differences on a ghost-extended copy.
    fn cell_gradient(&self, q: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let m = self.mesh;
        let n = m.n();
        let dx = m.grid.dx();
        let mut e = q.to_vec();
        self.extend(&mut e);
        let mut gx = vec![0.0; n * n];
        let mut gy = vec![0.0; n * n];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                if m.grid.fluid[c] {
                    gx[c] = (e[c + 1] - e[c - 1]) / (2.0 * dx);
                    gy[c] = (e[c + n] - e[c - n]) / (2.0 * dx);
                }
            }
        }
        (gx, gy)
    }

    /// Convective term `(u . grad) q` with the cell velocity and central differences on a
    /// ghost-extended copy of `q`.
    pub fn advection(&self, q: &[f64]) -> Vec<f64> {
        let m = self.mesh;
        let n = m.n();
        let dx = m.grid.dx();
        let mut e = q.to_vec();
        self.extend(&mut e);
        let mut out = vec![0.0; n * n];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                if m.grid.fluid[c] {
                    out[c] = (self.u[c] * (e[c + 1] - e[c - 1])
                        + self.v[c] * (e[c + n] - e[c - n]))
                        / (2.0 * dx);
                }
            }
        }
        out
    }

    /// Cell-centre Laplacian of a cell field with Dirichlet wall values (Gibou links).
    fn laplacian(&self, q: &[f64], bval: impl Fn(f64, f64) -> f64) -> Vec<f64> {
        use super::viscous::Link;
        let g = &self.mesh.grid;
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
        let m = self.mesh;
        let n = m.n();
        let wv = &self.wall_velocity;
        let lu = self.laplacian(&self.u, |x, y| wv(x, y).0);
        let lv = self.laplacian(&self.v, |x, y| wv(x, y).1);
        let au = self.advection(&self.u);
        let av = self.advection(&self.v);
        let (mut ax, mut ay) = (vec![0.0; n * n], vec![0.0; n * n]);
        for c in 0..n * n {
            if m.grid.fluid[c] {
                ax[c] = -au[c] + self.nu * lu[c];
                ay[c] = -av[c] + self.nu * lv[c];
            }
        }
        let (u0, v0, fu0, fv0) = (
            std::mem::replace(&mut self.u, ax),
            std::mem::replace(&mut self.v, ay),
            self.fu.clone(),
            self.fv.clone(),
        );
        let mut a = std::mem::take(&mut self.u);
        let mut b = std::mem::take(&mut self.v);
        self.extend(&mut a);
        self.extend(&mut b);
        self.u = a;
        self.v = b;
        self.faces_from_cells();
        self.p = self.project(1.0);
        self.u = u0;
        self.v = v0;
        self.fu = fu0;
        self.fv = fv0;
    }

    /// One time step (BDF2 + AB2, incremental pressure correction; first step BDF1 + AB1).
    pub fn step(&mut self, dt: f64) {
        let m = self.mesh;
        let nn = m.n() * m.n();
        let (mut au, mut av) = (self.advection(&self.u), self.advection(&self.v));
        if self.skip_advection {
            au.fill(0.0);
            av.fill(0.0);
        }
        let (gpx, gpy) = self.cell_gradient(&self.p);
        let bdf2 = matches!(&self.prev, Some(pv) if (pv.dt - dt).abs() < 1e-14 * dt);
        let sigma = if bdf2 { 1.5 } else { 1.0 } / (self.nu * dt);
        let (mut su, mut sv) = (vec![0.0; nn], vec![0.0; nn]);
        let (mut xu, mut xv) = (vec![0.0; nn], vec![0.0; nn]);
        for c in 0..nn {
            if !m.grid.fluid[c] {
                continue;
            }
            let (ex_u, ex_v, pu, pv) = match (&self.prev, bdf2) {
                (Some(pv), true) => (
                    2.0 * au[c] - pv.au[c],
                    2.0 * av[c] - pv.av[c],
                    (4.0 * self.u[c] - pv.u[c]) / 3.0,
                    (4.0 * self.v[c] - pv.v[c]) / 3.0,
                ),
                _ => (au[c], av[c], self.u[c], self.v[c]),
            };
            su[c] = (-ex_u - gpx[c]) / self.nu;
            sv[c] = (-ex_v - gpy[c]) / self.nu;
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
        let helm = m.grid.helmholtz(sigma);
        let wv = &self.wall_velocity;
        helm.solve(&mut xu, Some(&su), |x, y| wv(x, y).0, 1e-12);
        helm.solve(&mut xv, Some(&sv), |x, y| wv(x, y).1, 1e-12);
        self.extend(&mut xu);
        self.extend(&mut xv);
        self.u = xu;
        self.v = xv;
        self.faces_from_cells();
        if !self.skip_projection {
            let dt_eff = if bdf2 { dt * 2.0 / 3.0 } else { dt };
            let phi = self.project(dt_eff);
            for (p, f) in self.p.iter_mut().zip(&phi) {
                *p += f;
            }
        }
        self.prev = Some(old);
        self.time += dt;
    }

    /// Kinetic energy `sum 0.5 |u|^2 dx^2` over fluid cell centres.
    pub fn kinetic_energy(&self) -> f64 {
        let g = &self.mesh.grid;
        let dx2 = g.dx() * g.dx();
        (0..g.n * g.n)
            .filter(|&c| g.fluid[c])
            .map(|c| 0.5 * (self.u[c] * self.u[c] + self.v[c] * self.v[c]) * dx2)
            .sum()
    }

    pub fn max_face_divergence(&self) -> f64 {
        let n = self.mesh.n();
        (0..n * n)
            .filter(|&c| self.mesh.pressure_active[c])
            .map(|c| self.divergence(c % n, c / n).abs())
            .fold(0.0, f64::max)
    }
}
