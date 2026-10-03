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
use super::levelset::LevelSet;
use super::multigrid::PoissonSolver;
use super::viscous::{FluidGrid, Helmholtz, Link};

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
    /// Fraction of each cell inside the domain (8 x 8 sub-sampling).
    pub cell_fraction: Vec<f64>,
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

fn ghost_layers(fluid: &[bool], open: &[bool], n: usize, depth_max: u8) -> Vec<u8> {
    let mut layer: Vec<u8> = fluid.iter().map(|&f| if f { 0 } else { 255 }).collect();
    for depth in 1..=depth_max {
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
            layer[c] = depth_max;
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
        let mut cell_fraction = vec![0.0; n * n];
        const SUB: usize = 8;
        for j in 0..n {
            for i in 0..n {
                let mut inside = 0usize;
                for b in 0..SUB {
                    for a in 0..SUB {
                        let x = -half + (i as f64 + (a as f64 + 0.5) / SUB as f64) * dx;
                        let y = -half + (j as f64 + (b as f64 + 0.5) / SUB as f64) * dx;
                        if phi(x, y) < 0.0 {
                            inside += 1;
                        }
                    }
                }
                cell_fraction[i + n * j] = inside as f64 / (SUB * SUB) as f64;
            }
        }
        let open_u: Vec<bool> = au.iter().map(|&a| a > 0.0).collect();
        let open_v: Vec<bool> = av.iter().map(|&a| a > 0.0).collect();
        let layer_u = ghost_layers(&gu.fluid, &open_u, n, 2);
        let layer_v = ghost_layers(&gv.fluid, &open_v, n, 2);
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
            cell_fraction,
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
    let depth_max = layer
        .iter()
        .copied()
        .filter(|&l| l < 255)
        .max()
        .unwrap_or(0);
    for depth in 1..=depth_max {
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

/// Relative residual targets of the pressure and viscous solves.
const TOL_PROJECT: f64 = 1e-9;
const TOL_VISCOUS: f64 = 1e-9;

type WallVelocity<'a> = Box<dyn Fn(f64, f64) -> (f64, f64) + 'a>;

struct Prev {
    u: Vec<f64>,
    v: Vec<f64>,
    au: Vec<f64>,
    av: Vec<f64>,
    dt: f64,
}

/// Free-surface state: the level set (liquid where `psi < 0`) and body acceleration.
pub struct Surface {
    pub ls: LevelSet,
    /// Body acceleration `(gx, gy)` in m/s^2.
    pub gravity: (f64, f64),
    /// Liquid area the volume correction restores every step.
    pub volume0: f64,
    /// Reinitialisation sweeps per step.
    pub reinit_iterations: usize,
}

/// Smallest interface fraction used by the ghost-fluid pressure condition.
const THETA_MIN_FS: f64 = 1e-2;

/// Per-step liquid topology: which cells carry pressure, how faces connect them and which
/// face nodes are viscous unknowns.
pub struct Liquid {
    /// Cell carries a pressure unknown.
    pub cell: Vec<bool>,
    /// Face kinds (`u`, `v`): 0 inactive, 1 liquid-liquid, 2 liquid on the low side and air on
    /// the high side, 3 the reverse.
    pub ku: Vec<u8>,
    pub kv: Vec<u8>,
    /// Interface fraction of the centre-to-centre distance measured from the liquid centre.
    pub thu: Vec<f64>,
    pub thv: Vec<f64>,
    pub active_u: Vec<bool>,
    pub active_v: Vec<bool>,
    pub layer_u: Vec<u8>,
    pub layer_v: Vec<u8>,
}

impl Liquid {
    /// Everything inside the drum is liquid (no free surface).
    fn full(sm: &StaggeredMesh) -> Self {
        let n = sm.n();
        let kind = |a: &[f64]| a.iter().map(|&w| u8::from(w > 0.0)).collect::<Vec<u8>>();
        Self {
            cell: sm.mesh.pressure_active.clone(),
            ku: kind(&sm.au),
            kv: kind(&sm.av),
            thu: vec![1.0; n * n],
            thv: vec![1.0; n * n],
            active_u: sm.gu.fluid.clone(),
            active_v: sm.gv.fluid.clone(),
            layer_u: sm.layer_u.clone(),
            layer_v: sm.layer_v.clone(),
        }
    }

    fn from_level_set(sm: &StaggeredMesh, psi: &[f64]) -> Self {
        let n = sm.n();
        let cell: Vec<bool> = (0..n * n)
            .map(|c| psi[c] < 0.0 && sm.mesh.pressure_active[c])
            .collect();
        let mut out = Self {
            cell,
            ku: vec![0; n * n],
            kv: vec![0; n * n],
            thu: vec![1.0; n * n],
            thv: vec![1.0; n * n],
            active_u: vec![false; n * n],
            active_v: vec![false; n * n],
            layer_u: Vec::new(),
            layer_v: Vec::new(),
        };
        let classify = |lo: usize, hi: usize, aperture: f64, cell: &[bool]| -> (u8, f64) {
            if aperture <= 0.0 {
                return (0, 1.0);
            }
            match (cell[lo], cell[hi]) {
                (true, true) => (1, 1.0),
                (true, false) => (2, (psi[lo] / (psi[lo] - psi[hi])).clamp(THETA_MIN_FS, 1.0)),
                (false, true) => (3, (psi[hi] / (psi[hi] - psi[lo])).clamp(THETA_MIN_FS, 1.0)),
                (false, false) => (0, 1.0),
            }
        };
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if i >= 1 {
                    (out.ku[c], out.thu[c]) = classify(c - 1, c, sm.au[c], &out.cell);
                }
                if j >= 1 {
                    (out.kv[c], out.thv[c]) = classify(c - n, c, sm.av[c], &out.cell);
                }
                out.active_u[c] = out.ku[c] != 0 && sm.gu.fluid[c];
                out.active_v[c] = out.kv[c] != 0 && sm.gv.fluid[c];
            }
        }
        let open_u: Vec<bool> = out.ku.iter().map(|&k| k != 0).collect();
        let open_v: Vec<bool> = out.kv.iter().map(|&k| k != 0).collect();
        out.layer_u = ghost_layers(&out.active_u, &open_u, n, 3);
        out.layer_v = ghost_layers(&out.active_v, &open_v, n, 3);
        out
    }

    /// Aperture-weighted pressure Laplacian for this topology.
    fn poisson(&self, sm: &StaggeredMesh) -> PoissonSolver {
        let n = sm.n();
        let mut wx = vec![0.0; (n + 1) * n];
        let mut wy = vec![0.0; n * (n + 1)];
        let mut extra = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                match self.ku[c] {
                    1 => wx[i + (n + 1) * j] = sm.au[c],
                    2 => extra[c - 1] += sm.au[c] / self.thu[c],
                    3 => extra[c] += sm.au[c] / self.thu[c],
                    _ => {}
                }
                match self.kv[c] {
                    1 => wy[c] = sm.av[c],
                    2 => extra[c - n] += sm.av[c] / self.thv[c],
                    3 => extra[c] += sm.av[c] / self.thv[c],
                    _ => {}
                }
            }
        }
        PoissonSolver::from_weights(n, wx, wy, extra)
    }
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
    /// Helmholtz solvers for the last `sigma` (rebuilt only when `sigma` changes).
    helm: Option<(f64, Helmholtz<'a>, Helmholtz<'a>)>,
    pub liq: Liquid,
    pub surface: Option<Surface>,
}

impl<'a> StaggeredFlow<'a> {
    pub fn new(
        sm: &'a StaggeredMesh,
        nu: f64,
        wall_velocity: impl Fn(f64, f64) -> (f64, f64) + 'a,
    ) -> Self {
        let n = sm.n();
        let liq = Liquid::full(sm);
        let poisson = liq.poisson(sm);
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
            helm: None,
            liq,
            surface: None,
        }
    }

    /// Switches to a free-surface liquid described by `ls` (liquid where `psi < 0`) under the body
    /// acceleration `gravity`. The level set is given the domain cell fractions; its liquid area is
    /// the volume that is conserved from here on.
    pub fn enable_free_surface(&mut self, mut ls: LevelSet, gravity: (f64, f64)) {
        ls.weight = self.sm.cell_fraction.clone();
        let volume0 = ls.volume();
        self.liq = Liquid::from_level_set(self.sm, &ls.psi);
        self.poisson = self.liq.poisson(self.sm);
        self.surface = Some(Surface {
            ls,
            gravity,
            volume0,
            reinit_iterations: 2,
        });
        self.prev = None;
        self.helm = None;
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
                if self.liq.layer_u[c] < 255 {
                    let (x, y) = sm.gu.centre(i, j);
                    self.u[c] = f(x, y).0;
                }
                if self.liq.layer_v[c] < 255 {
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
        let fu = |c: usize| {
            if self.liq.ku[c] != 0 {
                sm.au[c] * (self.u[c] + self.cu[c])
            } else {
                0.0
            }
        };
        let fv = |c: usize| {
            if self.liq.kv[c] != 0 {
                sm.av[c] * (self.v[c] + self.cv[c])
            } else {
                0.0
            }
        };
        let right = if i + 1 < n { fu(i + 1 + n * j) } else { 0.0 };
        let top = if j + 1 < n { fv(i + n * (j + 1)) } else { 0.0 };
        right - fu(i + n * j) + top - fv(i + n * j)
    }

    pub fn max_divergence(&self) -> f64 {
        let n = self.sm.n();
        (0..n * n)
            .filter(|&c| self.liq.cell[c])
            .map(|c| self.divergence(c % n, c / n).abs())
            .fold(0.0, f64::max)
    }

    /// Exact projection with the weighted operator; `dt_eff` is the pressure-update time scale.
    /// Returns the pressure increment; updates `u, v` on every active face (then re-extends).
    fn project(&mut self, dt_eff: f64) -> Vec<f64> {
        let sm = self.sm;
        let n = sm.n();
        let dx = sm.dx();
        self.update_flux_corrections();
        let mut rhs = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if self.liq.cell[c] {
                    rhs[c] = -(dx / dt_eff) * self.divergence(i, j);
                }
            }
        }
        let mut phi = vec![0.0; n * n];
        self.poisson.solve(&rhs, &mut phi, TOL_PROJECT, 300);
        let k = dt_eff / dx;
        for j in 0..n {
            for i in 1..n {
                let c = i + n * j;
                match self.liq.ku[c] {
                    1 => self.u[c] -= k * (phi[c] - phi[c - 1]),
                    2 => self.u[c] -= k * (0.0 - phi[c - 1]) / self.liq.thu[c],
                    3 => self.u[c] -= k * (phi[c] - 0.0) / self.liq.thu[c],
                    _ => {}
                }
            }
        }
        for j in 1..n {
            for i in 0..n {
                let c = i + n * j;
                match self.liq.kv[c] {
                    1 => self.v[c] -= k * (phi[c] - phi[c - n]),
                    2 => self.v[c] -= k * (0.0 - phi[c - n]) / self.liq.thv[c],
                    3 => self.v[c] -= k * (phi[c] - 0.0) / self.liq.thv[c],
                    _ => {}
                }
            }
        }
        self.projected_divergence = self.max_divergence();
        extend(&mut self.u, &self.liq.layer_u, n);
        extend(&mut self.v, &self.liq.layer_v, n);
        phi
    }

    /// Face gradient of a cell field, `(dp/dx at u nodes, dp/dy at v nodes)`; air has `p = 0` at
    /// the interface fraction.
    fn pressure_gradient(&self, p: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let sm = self.sm;
        let n = sm.n();
        let dx = sm.dx();
        let mut gx = vec![0.0; n * n];
        let mut gy = vec![0.0; n * n];
        for j in 1..n {
            for i in 1..n {
                let c = i + n * j;
                gx[c] = match self.liq.ku[c] {
                    1 => (p[c] - p[c - 1]) / dx,
                    2 => -p[c - 1] / (self.liq.thu[c] * dx),
                    3 => p[c] / (self.liq.thu[c] * dx),
                    _ => 0.0,
                };
                gy[c] = match self.liq.kv[c] {
                    1 => (p[c] - p[c - n]) / dx,
                    2 => -p[c - n] / (self.liq.thv[c] * dx),
                    3 => p[c] / (self.liq.thv[c] * dx),
                    _ => 0.0,
                };
            }
        }
        (gx, gy)
    }

    /// Convective terms `(U . grad) u` at `u` nodes and `(U . grad) v` at `v` nodes, for active
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
                if self.liq.active_u[c] {
                    let vb = 0.25 * (v[c - 1] + v[c] + v[c - 1 + n] + v[c + n]);
                    au[c] = u[c] * (u[c + 1] - u[c - 1]) * inv + vb * (u[c + n] - u[c - n]) * inv;
                }
                if self.liq.active_v[c] {
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
            if self.liq.active_u[c] {
                ax[c] = -au[c] + self.nu * lu[c];
            }
            if self.liq.active_v[c] {
                ay[c] = -av[c] + self.nu * lv[c];
            }
        }
        extend(&mut ax, &self.liq.layer_u, n);
        extend(&mut ay, &self.liq.layer_v, n);
        let u0 = std::mem::replace(&mut self.u, ax);
        let v0 = std::mem::replace(&mut self.v, ay);
        self.p = self.project(1.0);
        self.u = u0;
        self.v = v0;
    }

    /// One time step (BDF2 + AB2, incremental pressure correction; first step BDF1 + AB1). With a
    /// free surface the liquid topology is rebuilt from the level set first, and the level set is
    /// advected with the new velocity at the end.
    pub fn step(&mut self, dt: f64) {
        let sm = self.sm;
        let n = sm.n();
        let nn = n * n;
        let gravity = self.surface.as_ref().map(|s| s.gravity);
        if let Some(surf) = &self.surface {
            self.liq = Liquid::from_level_set(sm, &surf.ls.psi);
            self.poisson = self.liq.poisson(sm);
            self.helm = None;
            for c in 0..nn {
                if !self.liq.cell[c] {
                    self.p[c] = 0.0;
                }
            }
            extend(&mut self.u, &self.liq.layer_u, n);
            extend(&mut self.v, &self.liq.layer_v, n);
        }
        let (au, av) = self.advection(&self.u, &self.v);
        let (gpx, gpy) = self.pressure_gradient(&self.p);
        let bdf2 = matches!(&self.prev, Some(pv) if (pv.dt - dt).abs() < 1e-14 * dt);
        let sigma = if bdf2 { 1.5 } else { 1.0 } / (self.nu * dt);
        let (gx, gy) = gravity.unwrap_or((0.0, 0.0));
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
            su[c] = (-eu - gpx[c] + gx) / self.nu;
            sv[c] = (-ev - gpy[c] + gy) / self.nu;
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
        if !matches!(&self.helm, Some((s, _, _)) if (s - sigma).abs() <= 1e-12 * sigma) {
            self.helm = Some((
                sigma,
                sm.gu.helmholtz_masked(sigma, &self.liq.active_u),
                sm.gv.helmholtz_masked(sigma, &self.liq.active_v),
            ));
        }
        let wv = &self.wall_velocity;
        if let Some((_, hu, hv)) = &self.helm {
            hu.solve(&mut xu, Some(&su), |x, y| wv(x, y).0, TOL_VISCOUS);
            hv.solve(&mut xv, Some(&sv), |x, y| wv(x, y).1, TOL_VISCOUS);
        }
        extend(&mut xu, &self.liq.layer_u, n);
        extend(&mut xv, &self.liq.layer_v, n);
        self.u = xu;
        self.v = xv;
        let dt_eff = if bdf2 { dt * 2.0 / 3.0 } else { dt };
        let phi = self.project(dt_eff);
        for (p, f) in self.p.iter_mut().zip(&phi) {
            *p += f;
        }
        self.prev = Some(old);
        self.time += dt;
        if let Some(surf) = &mut self.surface {
            let (mut uc, mut vc) = surf.ls.centre_velocity(&self.u, &self.v);
            let band = 4.0 * surf.ls.dx;
            for c in 0..nn {
                if surf.ls.psi[c].abs() > band {
                    uc[c] = 0.0;
                    vc[c] = 0.0;
                }
            }
            surf.ls.advect(&uc, &vc, dt);
            surf.ls.reinitialize(surf.reinit_iterations);
            surf.ls.correct_volume(surf.volume0);
        }
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
