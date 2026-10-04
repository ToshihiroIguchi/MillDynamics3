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
use super::levelset::{weno5, LevelSet};
use super::multigrid::PoissonSolver;
use super::rheology::{alg2_update, HerschelBulkley, Sym};
use super::variational::{pcg, ViscousOperator};
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
    /// Net flux (per unit length) of the moving solid bodies' surface velocity out of each cell's
    /// fluid part, see [`StaggeredMesh::set_body_flux`]; zero for static solids.
    pub body_flux: Vec<f64>,
}

/// A rigid disc moving through the domain: centre, radius, translation velocity and spin.
#[derive(Clone, Copy, Debug)]
pub struct Disc {
    pub cx: f64,
    pub cy: f64,
    pub r: f64,
    pub ux: f64,
    pub uy: f64,
    pub omega: f64,
}

impl Disc {
    /// Signed distance to the disc surface (negative inside).
    pub fn sdf(&self, x: f64, y: f64) -> f64 {
        ((x - self.cx).powi(2) + (y - self.cy).powi(2)).sqrt() - self.r
    }

    /// Velocity of the material point of the disc at `(x, y)`.
    pub fn velocity(&self, x: f64, y: f64) -> (f64, f64) {
        (
            self.ux - self.omega * (y - self.cy),
            self.uy + self.omega * (x - self.cx),
        )
    }
}

/// Length and first moment (both as fractions of the segment `p0 -> p1`) of the part with
/// `phi < 0`; crossings are located by bisection.
fn inside_moments(phi: &impl Fn(f64, f64) -> f64, p0: (f64, f64), p1: (f64, f64)) -> (f64, f64) {
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
    (len, mom)
}

/// Centroid (as a fraction `t` of the segment `p0 -> p1`) of the part with `phi < 0`; 0.5 if none.
fn open_centroid(phi: &impl Fn(f64, f64) -> f64, p0: (f64, f64), p1: (f64, f64)) -> f64 {
    const M: usize = 16;
    let at = |t: f64| (p0.0 + t * (p1.0 - p0.0), p0.1 + t * (p1.1 - p0.1));
    let inside = |t: f64| {
        let (x, y) = at(t);
        phi(x, y) < 0.0
    };
    // A segment farther from the boundary than its own length is entirely on one side
    // (phi is a distance-like, 1-Lipschitz field): the centroid is the midpoint.
    let seg = (p1.0 - p0.0).hypot(p1.1 - p0.1);
    let (mx, my) = at(0.5);
    if phi(mx, my).abs() > seg {
        return 0.5;
    }
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
                // Cells far from the boundary are entirely in or out (phi is 1-Lipschitz for
                // the signed-distance-like fields used here, so |phi| > half a diagonal decides).
                let phi_c = phi(-half + (i as f64 + 0.5) * dx, -half + (j as f64 + 0.5) * dx);
                if phi_c.abs() > 0.75 * dx {
                    cell_fraction[i + n * j] = if phi_c < 0.0 { 1.0 } else { 0.0 };
                    continue;
                }
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
        let layer_u = ghost_layers(&gu.fluid, &open_u, n, 3);
        let layer_v = ghost_layers(&gv.fluid, &open_v, n, 3);
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
            body_flux: vec![0.0; n * n],
        }
    }

    /// Sets `body_flux` for moving discs: the surface velocity of a rigid body is divergence free,
    /// so the flux it drives through the part of a cell boundary that the body covers equals,
    /// with the opposite sign, the flux through the body's own surface inside the cell; the
    /// pressure projection must add it to the divergence of the open faces.
    pub fn set_body_flux(&mut self, discs: &[Disc]) {
        let n = self.n();
        let dx = self.dx();
        let half = self.mesh.grid.half;
        let mut covered_u = vec![0.0; n * n];
        let mut covered_v = vec![0.0; n * n];
        for d in discs {
            let phi = |x: f64, y: f64| d.sdf(x, y);
            for j in 0..n {
                for i in 0..n {
                    let c = i + n * j;
                    let (x0, y0) = (-half + i as f64 * dx, -half + j as f64 * dx);
                    if (x0 - d.cx).abs() > d.r + 2.0 * dx || (y0 - d.cy).abs() > d.r + 2.0 * dx {
                        continue;
                    }
                    // u-face: vertical segment at x0 from y0 to y0 + dx, normal velocity ub_x.
                    let (len, mom) = inside_moments(&phi, (x0, y0), (x0, y0 + dx));
                    if len > 0.0 {
                        let alpha = d.ux - d.omega * (y0 - d.cy);
                        covered_u[c] += alpha * len - d.omega * dx * mom;
                    }
                    let (len, mom) = inside_moments(&phi, (x0, y0), (x0 + dx, y0));
                    if len > 0.0 {
                        let alpha = d.uy + d.omega * (x0 - d.cx);
                        covered_v[c] += alpha * len + d.omega * dx * mom;
                    }
                }
            }
        }
        self.body_flux.fill(0.0);
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                let right = if i + 1 < n { covered_u[c + 1] } else { 0.0 };
                let top = if j + 1 < n { covered_v[c + n] } else { 0.0 };
                self.body_flux[c] = right - covered_u[c] + top - covered_v[c];
            }
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
pub fn extend(q: &mut [f64], layer: &[u8], n: usize, max_order: u8) {
    extend_ordered(q, layer, n, &|_| max_order, &|_| None);
}

/// Same with the highest extrapolation order chosen per node.
pub fn extend_ordered(
    q: &mut [f64],
    layer: &[u8],
    n: usize,
    order_of: &dyn Fn(usize) -> u8,
    bounded: &dyn Fn(usize) -> Option<f64>,
) {
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
                let max_order = order_of(c);
                let at = |k: isize, di: isize, dj: isize| -> Option<f64> {
                    let (ii, jj) = (i as isize + k * di, j as isize + k * dj);
                    (ok(ii, jj) && layer[ii as usize + n * jj as usize] < depth)
                        .then(|| snapshot[ii as usize + n * jj as usize])
                };
                let mut best = (0u8, 0.0f64, 0.0f64);
                for (di, dj) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                    let (q1, q2, q3) = (at(1, di, dj), at(2, di, dj), at(3, di, dj));
                    let mut cand = match (q1, q2, q3) {
                        (Some(a), Some(b), Some(c3)) if max_order >= 3 => {
                            Some((3u8, 3.0 * a - 3.0 * b + c3))
                        }
                        (Some(a), Some(b), _) if max_order >= 2 => Some((2, 2.0 * a - b)),
                        (Some(a), _, _) => Some((1, a)),
                        _ => None,
                    };
                    if let Some(slack) = bounded(c) {
                        // Monotone extrapolation: stay within the range of the data it comes
                        // from, widened by `slack` times that range on both sides.
                        let vals: Vec<f64> = [q1, q2, q3].iter().flatten().copied().collect();
                        if let Some((o, val)) = cand {
                            let lo = vals.iter().cloned().fold(f64::INFINITY, f64::min);
                            let hi = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                            let w = slack * (hi - lo);
                            cand = Some((o, val.max(lo - w).min(hi + w)));
                        }
                    }
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

#[derive(Clone)]
struct Prev {
    u: Vec<f64>,
    v: Vec<f64>,
    au: Vec<f64>,
    av: Vec<f64>,
    dt: f64,
}

/// Free-surface state: the level set (liquid where `psi < 0`) and body acceleration.
#[derive(Clone)]
pub struct Surface {
    pub ls: LevelSet,
    /// Body acceleration `(gx, gy)` in m/s^2.
    pub gravity: (f64, f64),
    /// Liquid area the volume correction restores every step.
    pub volume0: f64,
    /// Reinitialisation sweeps per step.
    pub reinit_iterations: usize,
    /// Smallest interface fraction of the ghost-fluid pressure condition.
    pub theta_min: f64,
    /// Use the interface position predicted at the end of the step for the pressure condition.
    pub predictor: bool,
    /// Restore the liquid area every step by shifting the level set.
    pub correct_volume: bool,
}

impl Surface {
    /// Cell-centre velocity for advecting the level set: face averages within a band of the
    /// interface, zero elsewhere.
    fn cell_velocity(&self, u: &[f64], v: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let (mut uc, mut vc) = self.ls.centre_velocity(u, v);
        let band = 4.0 * self.ls.dx;
        for c in 0..uc.len() {
            if self.ls.psi[c].abs() > band {
                uc[c] = 0.0;
                vc[c] = 0.0;
            }
        }
        (uc, vc)
    }
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
    /// Unit normals of the level set (pointing into the air) at the `u` and `v` nodes.
    pub normal_u: Vec<(f64, f64)>,
    pub normal_v: Vec<(f64, f64)>,
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
            normal_u: vec![(0.0, 0.0); n * n],
            normal_v: vec![(0.0, 0.0); n * n],
        }
    }

    fn from_level_set(sm: &StaggeredMesh, psi: &[f64], theta_min: f64) -> Self {
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
            normal_u: vec![(0.0, 0.0); n * n],
            normal_v: vec![(0.0, 0.0); n * n],
        };
        let classify = |lo: usize, hi: usize, aperture: f64, cell: &[bool]| -> (u8, f64) {
            if aperture <= 0.0 {
                return (0, 1.0);
            }
            match (cell[lo], cell[hi]) {
                (true, true) => (1, 1.0),
                (true, false) => (2, (psi[lo] / (psi[lo] - psi[hi])).clamp(theta_min, 1.0)),
                (false, true) => (3, (psi[hi] / (psi[hi] - psi[lo])).clamp(theta_min, 1.0)),
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
        let at = |i: isize, j: isize| {
            let m = n as isize - 1;
            psi[(i.clamp(0, m) + n as isize * j.clamp(0, m)) as usize]
        };
        let dx = sm.dx();
        let unit = |gx: f64, gy: f64| {
            let g = (gx * gx + gy * gy).sqrt();
            if g > 1e-12 {
                (gx / g, gy / g)
            } else {
                (0.0, 0.0)
            }
        };
        for j in 0..n as isize {
            for i in 0..n as isize {
                let c = i as usize + n * j as usize;
                // u node between cells (i-1, j) and (i, j).
                let gx = (at(i, j) - at(i - 1, j)) / dx;
                let gy = (at(i - 1, j + 1) + at(i, j + 1) - at(i - 1, j - 1) - at(i, j - 1))
                    / (4.0 * dx);
                out.normal_u[c] = unit(gx, gy);
                // v node between cells (i, j-1) and (i, j).
                let gy = (at(i, j) - at(i, j - 1)) / dx;
                let gx = (at(i + 1, j - 1) + at(i + 1, j) - at(i - 1, j - 1) - at(i - 1, j))
                    / (4.0 * dx);
                out.normal_v[c] = unit(gx, gy);
            }
        }
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
    /// Multiplier of the linear-solver tolerances (1e-9 relative at 1).
    pub tol_scale: f64,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub p: Vec<f64>,
    /// Velocity right after the viscous solve of the last step (before the projection): the
    /// projection moves nodes that sit almost on a wall off their wall value, which a wall
    /// traction `(u_b - u) / theta` amplifies; loads on bodies are read from this field.
    pub u_star: Vec<f64>,
    pub v_star: Vec<f64>,
    wall_velocity: WallVelocity<'a>,
    prev: Option<Prev>,
    poisson: PoissonSolver,
    pub time: f64,
    /// Largest `|div|` (per cell, aperture weighted) right after the last projection.
    pub projected_divergence: f64,
    /// Iterations and residual of the last pressure solve.
    pub last_solve: (u32, f64),
    /// Flux sampling corrections `off * d(velocity)/d(face direction)`, held fixed during a
    /// projection so that it stays exact.
    cu: Vec<f64>,
    cv: Vec<f64>,
    /// Helmholtz solvers for the last `sigma` (rebuilt only when `sigma` changes).
    helm: Option<(f64, Helmholtz<'a>, Helmholtz<'a>)>,
    pub liq: Liquid,
    pub surface: Option<Surface>,
    /// Third-order upwind-biased advection instead of central differences.
    pub upwind: bool,
    /// Fifth-order WENO upwind advection (needs three ghost layers; overrides `upwind`).
    pub weno: bool,
    /// Drop the convective term (creeping flow; diagnostic).
    pub stokes: bool,
    /// Stress-free free-surface condition for the viscous solve (else zero normal derivative of
    /// each Cartesian component).
    pub stress_free_flux: bool,
    /// Stability over consistency at a free surface: clamp air-side ghost values to the range of
    /// the liquid data and use a zero advective gradient where the upstream node is air. Splashing
    /// flows survive (wall jets, thin sheets), but rotating flow along the surface is then
    /// reproduced only to a few percent that does not shrink with the grid (rigid-ring test).
    pub robust_surface: bool,
    /// Smagorinsky constant `C_s` of the explicit eddy viscosity (0 = off).
    pub smagorinsky: f64,
    /// Yield-stress / non-Newtonian fluid (augmented Lagrangian), see [`Self::set_rheology`].
    pub rheo: Option<Rheo>,
    /// Most ALG2 passes per step (1 = the multiplier lags one step) and the relative multiplier
    /// change below which the passes stop.
    pub alg_iterations: usize,
    pub alg_tolerance: f64,
    /// Under-relaxation of the multiplier step.
    pub alg_relax: f64,
    /// Passes the last step used.
    pub last_alg_passes: usize,
    /// Symmetric (variational) stress discretisation, see `mac/variational.rs`.
    pub variational: bool,
    visc_op: Option<ViscousOperator>,
    /// Diagnostic swirl body force `swirl * (-y, x)` (balanced only by stress).
    pub swirl: f64,
}

/// Augmented-Lagrangian state of a generalised-Newtonian fluid: the multiplier (the stress) and
/// the stress excess `S = lam - r d`. Normal components live at the cell centres (`*_c`, entries
/// `xx` and `yy`), the shear component at the cell corners (`*_k`, entry `xy`; corner `c` is the
/// lower-left corner of cell `c`), like the compact MAC discretisation of the strain rate.
#[derive(Clone)]
pub struct Rheo {
    pub law: HerschelBulkley,
    /// Augmentation parameter `r` (the implicit viscosity is `nu = r / 2`).
    pub r: f64,
    /// Newtonian viscosity solved with the wall-exact Gibou operator in the variational mode
    /// (the multiplier law then carries the remaining stress only).
    pub nu_base: f64,
    pub lam_c: Vec<Sym>,
    pub lam_k: Vec<Sym>,
    pub s_c: Vec<Sym>,
    pub s_k: Vec<Sym>,
}

/// Everything of a [`StaggeredFlow`] except the mesh and the wall-velocity closure: moving
/// bodies change the mesh every step, so the flow is carried from one mesh to the next.
#[derive(Clone)]
pub struct FlowState {
    pub nu: f64,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub p: Vec<f64>,
    prev: Option<Prev>,
    pub time: f64,
    pub surface: Option<Surface>,
    pub upwind: bool,
    pub weno: bool,
    pub stokes: bool,
    pub stress_free_flux: bool,
    pub robust_surface: bool,
    pub smagorinsky: f64,
    pub rheo: Option<Rheo>,
    pub alg_iterations: usize,
    pub alg_tolerance: f64,
    pub alg_relax: f64,
    pub variational: bool,
}

impl<'a> StaggeredFlow<'a> {
    /// Detaches the flow state from its mesh.
    pub fn into_state(self) -> FlowState {
        FlowState {
            nu: self.nu,
            u: self.u,
            v: self.v,
            p: self.p,
            prev: self.prev,
            time: self.time,
            surface: self.surface,
            upwind: self.upwind,
            weno: self.weno,
            stokes: self.stokes,
            stress_free_flux: self.stress_free_flux,
            robust_surface: self.robust_surface,
            smagorinsky: self.smagorinsky,
            rheo: self.rheo,
            alg_iterations: self.alg_iterations,
            alg_tolerance: self.alg_tolerance,
            alg_relax: self.alg_relax,
            variational: self.variational,
        }
    }

    /// Continues `state` on a new mesh `sm` (same grid, moved bodies). Velocities at nodes that
    /// the bodies uncovered are the ghost values extrapolated into the solid on the old mesh.
    pub fn from_state(
        sm: &'a StaggeredMesh,
        state: FlowState,
        wall_velocity: impl Fn(f64, f64) -> (f64, f64) + 'a,
    ) -> Self {
        let mut flow = Self::new(sm, state.nu, wall_velocity);
        flow.u = state.u;
        flow.v = state.v;
        flow.p = state.p;
        flow.prev = state.prev;
        flow.time = state.time;
        flow.upwind = state.upwind;
        flow.weno = state.weno;
        flow.stokes = state.stokes;
        flow.stress_free_flux = state.stress_free_flux;
        flow.robust_surface = state.robust_surface;
        flow.smagorinsky = state.smagorinsky;
        flow.rheo = state.rheo;
        flow.alg_iterations = state.alg_iterations;
        flow.alg_tolerance = state.alg_tolerance;
        flow.alg_relax = state.alg_relax;
        flow.variational = state.variational;
        if let Some(mut surf) = state.surface {
            surf.ls.weight = sm.cell_fraction.clone();
            flow.liq = Liquid::from_level_set(sm, &surf.ls.psi, surf.theta_min);
            flow.poisson = flow.liq.poisson(sm);
            flow.surface = Some(surf);
        }
        // Fresh nodes need valid values: re-extend into the new ghost layers.
        let mut tmp = std::mem::take(&mut flow.u);
        flow.extend_u(&mut tmp);
        flow.u = tmp;
        let mut tmp = std::mem::take(&mut flow.v);
        flow.extend_v(&mut tmp);
        flow.v = tmp;
        flow
    }

    /// Initialises what a moving body uncovered since `old` (the previous mesh): face nodes that
    /// were solid get the body's rigid velocity `rigid(x, y)` (they lie within a step's travel of
    /// its no-slip surface), cells that carried no pressure get the mean of their neighbours.
    pub fn initialise_fresh(
        &mut self,
        old: &StaggeredMesh,
        rigid: impl Fn(f64, f64) -> (f64, f64),
    ) {
        let sm = self.sm;
        let n = sm.n();
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if sm.gu.fluid[c] && !old.gu.fluid[c] {
                    let (x, y) = sm.gu.centre(i, j);
                    self.u[c] = rigid(x, y).0;
                    if let Some(pv) = &mut self.prev {
                        pv.u[c] = self.u[c];
                    }
                }
                if sm.gv.fluid[c] && !old.gv.fluid[c] {
                    let (x, y) = sm.gv.centre(i, j);
                    self.v[c] = rigid(x, y).1;
                    if let Some(pv) = &mut self.prev {
                        pv.v[c] = self.v[c];
                    }
                }
            }
        }
        let mut known: Vec<bool> = old.mesh.pressure_active.clone();
        for _ in 0..3 {
            let snapshot = known.clone();
            for j in 1..n - 1 {
                for i in 1..n - 1 {
                    let c = i + n * j;
                    if snapshot[c] || !sm.mesh.pressure_active[c] {
                        continue;
                    }
                    let (mut sum, mut cnt) = (0.0, 0.0);
                    for nb in [c - 1, c + 1, c - n, c + n] {
                        if snapshot[nb] {
                            sum += self.p[nb];
                            cnt += 1.0;
                        }
                    }
                    if cnt > 0.0 {
                        self.p[c] = sum / cnt;
                        known[c] = true;
                    }
                }
            }
        }
        let mut tmp = std::mem::take(&mut self.u);
        self.extend_u(&mut tmp);
        self.u = tmp;
        let mut tmp = std::mem::take(&mut self.v);
        self.extend_v(&mut tmp);
        self.v = tmp;
    }

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
            tol_scale: 1.0,
            u: vec![0.0; n * n],
            v: vec![0.0; n * n],
            p: vec![0.0; n * n],
            u_star: vec![0.0; n * n],
            v_star: vec![0.0; n * n],
            wall_velocity: Box::new(wall_velocity),
            prev: None,
            poisson,
            time: 0.0,
            projected_divergence: 0.0,
            last_solve: (0, 0.0),
            cu: vec![0.0; n * n],
            cv: vec![0.0; n * n],
            helm: None,
            liq,
            surface: None,
            upwind: false,
            weno: false,
            stokes: false,
            stress_free_flux: true,
            robust_surface: false,
            smagorinsky: 0.0,
            rheo: None,
            alg_iterations: 1,
            alg_tolerance: 1e-5,
            alg_relax: 1.0,
            last_alg_passes: 0,
            variational: false,
            visc_op: None,
            swirl: 0.0,
        }
    }

    /// Switches to a free-surface liquid described by `ls` (liquid where `psi < 0`) under the body
    /// acceleration `gravity`. The level set is given the domain cell fractions; its liquid area is
    /// the volume that is conserved from here on.
    pub fn enable_free_surface(&mut self, mut ls: LevelSet, gravity: (f64, f64)) {
        ls.weight = self.sm.cell_fraction.clone();
        let volume0 = ls.volume();
        self.liq = Liquid::from_level_set(self.sm, &ls.psi, THETA_MIN_FS);
        self.poisson = self.liq.poisson(self.sm);
        self.surface = Some(Surface {
            ls,
            gravity,
            volume0,
            reinit_iterations: 2,
            theta_min: THETA_MIN_FS,
            predictor: true,
            correct_volume: true,
        });
        self.prev = None;
        self.helm = None;
    }

    /// Ghost extension of a `u`-lattice field: quadratic extrapolation along the grid lines
    /// (walls and air side). With `robust_surface` the air-side values are clamped to the range of
    /// the liquid data they come from.
    fn extend_u(&self, q: &mut [f64]) {
        extend_ordered(q, &self.liq.layer_u, self.sm.n(), &|_| 3, &|c| {
            (self.robust_surface && self.sm.gu.fluid[c] && !self.liq.active_u[c]).then_some(0.0)
        });
    }

    fn extend_v(&self, q: &mut [f64]) {
        extend_ordered(q, &self.liq.layer_v, self.sm.n(), &|_| 3, &|c| {
            (self.robust_surface && self.sm.gv.fluid[c] && !self.liq.active_v[c]).then_some(0.0)
        });
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
        right - fu(i + n * j) + top - fv(i + n * j) + sm.body_flux[i + n * j]
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
        let st = self
            .poisson
            .solve(&rhs, &mut phi, TOL_PROJECT * self.tol_scale, 300);
        self.last_solve = (st.iterations, st.residual);
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
        let mut tmp = std::mem::take(&mut self.u);
        self.extend_u(&mut tmp);
        self.u = tmp;
        let mut tmp = std::mem::take(&mut self.v);
        self.extend_v(&mut tmp);
        self.v = tmp;
        phi
    }

    /// Right-hand-side data (times `dx`) for the dropped free-surface links of the two Helmholtz
    /// solves. A zero normal derivative of each Cartesian component is not the stress-free
    /// condition (`(grad u + grad u^T) n = 0`): it damps the rigid rotation of a liquid ring at the
    /// free surface by an error that does not shrink with the grid. Here the derivative along the
    /// link direction `e` is built from the stress-free normal derivative
    /// `n . grad u_i = -(grad u)^T n` and the tangential derivative `t . grad u_i` of the current
    /// (lagged) velocity: `g = (e.n) (-sum_j n_j A_ji) + (e.t) (sum_j t_j A_ij)` with
    /// `A_ij = d u_i / d x_j`.
    fn free_surface_flux(&self) -> (Vec<f64>, Vec<f64>) {
        let sm = self.sm;
        let n = sm.n();
        let dx = sm.dx();
        let mut fu = vec![0.0; n * n];
        let mut fv = vec![0.0; n * n];
        let dirs: [(f64, f64); 4] = [(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                if self.liq.active_u[c] {
                    let a = self.gradient_at_u(c);
                    fu[c] = self.link_flux(
                        &sm.gu,
                        &self.liq.active_u,
                        c,
                        0,
                        &a,
                        self.liq.normal_u[c],
                        &dirs,
                        dx,
                    );
                }
                if self.liq.active_v[c] {
                    let a = self.gradient_at_v(c);
                    fv[c] = self.link_flux(
                        &sm.gv,
                        &self.liq.active_v,
                        c,
                        1,
                        &a,
                        self.liq.normal_v[c],
                        &dirs,
                        dx,
                    );
                }
            }
        }
        (fu, fv)
    }

    /// Velocity gradient `A[component][direction] = d u_i / d x_j` at the `u` node `c`.
    fn gradient_at_u(&self, c: usize) -> [[f64; 2]; 2] {
        let n = self.sm.n();
        let dx = self.sm.dx();
        let (u, v) = (&self.u, &self.v);
        [
            [
                (u[c + 1] - u[c - 1]) / (2.0 * dx),
                (u[c + n] - u[c - n]) / (2.0 * dx),
            ],
            [
                (v[c] + v[c + n] - v[c - 1] - v[c - 1 + n]) / (2.0 * dx),
                0.5 * (v[c - 1 + n] - v[c - 1] + v[c + n] - v[c]) / dx,
            ],
        ]
    }

    /// Same at the `v` node `c`.
    fn gradient_at_v(&self, c: usize) -> [[f64; 2]; 2] {
        let n = self.sm.n();
        let dx = self.sm.dx();
        let (u, v) = (&self.u, &self.v);
        [
            [
                0.5 * (u[c + 1] - u[c] + u[c + 1 - n] - u[c - n]) / dx,
                (u[c] + u[c + 1] - u[c - n] - u[c + 1 - n]) / (2.0 * dx),
            ],
            [
                (v[c + 1] - v[c - 1]) / (2.0 * dx),
                (v[c + n] - v[c - n]) / (2.0 * dx),
            ],
        ]
    }

    /// Smagorinsky eddy-viscosity acceleration `div(nu_t grad u)` (explicit, component-wise,
    /// `nu_t = (C_s dx)^2 |S|`) at the active nodes; free-surface and wall links carry no flux.
    fn eddy_acceleration(&self) -> (Vec<f64>, Vec<f64>) {
        let sm = self.sm;
        let n = sm.n();
        let dx = sm.dx();
        let cs2 = (self.smagorinsky * dx).powi(2);
        let strain = |a: [[f64; 2]; 2]| {
            (2.0 * (a[0][0] * a[0][0] + a[1][1] * a[1][1]) + (a[0][1] + a[1][0]).powi(2)).sqrt()
        };
        let mut nut_u = vec![0.0; n * n];
        let mut nut_v = vec![0.0; n * n];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                if self.liq.active_u[c] {
                    nut_u[c] = cs2 * strain(self.gradient_at_u(c));
                }
                if self.liq.active_v[c] {
                    nut_v[c] = cs2 * strain(self.gradient_at_v(c));
                }
            }
        }
        let lap = |q: &[f64], nut: &[f64], active: &[bool], c: usize| -> f64 {
            let mut s = 0.0;
            for k in [c - 1, c + 1, c - n, c + n] {
                if active[k] {
                    s += 0.5 * (nut[c] + nut[k]) * (q[k] - q[c]);
                }
            }
            s / (dx * dx)
        };
        let mut au = vec![0.0; n * n];
        let mut av = vec![0.0; n * n];
        for j in 2..n - 2 {
            for i in 2..n - 2 {
                let c = i + n * j;
                if self.liq.active_u[c] {
                    au[c] = lap(&self.u, &nut_u, &self.liq.active_u, c);
                }
                if self.liq.active_v[c] {
                    av[c] = lap(&self.v, &nut_v, &self.liq.active_v, c);
                }
            }
        }
        (au, av)
    }

    /// Sum over the dropped links of node `c` of `dx * g`, see [`Self::free_surface_flux`].
    #[allow(clippy::too_many_arguments)]
    fn link_flux(
        &self,
        grid: &FluidGrid,
        active: &[bool],
        c: usize,
        comp: usize,
        a: &[[f64; 2]; 2],
        normal: (f64, f64),
        dirs: &[(f64, f64); 4],
        dx: f64,
    ) -> f64 {
        let (nx, ny) = normal;
        if nx == 0.0 && ny == 0.0 {
            return 0.0;
        }
        let (tx, ty) = (-ny, nx);
        let mut total = 0.0;
        for (k, link) in grid.links[c].iter().enumerate() {
            if let Link::Fluid(nb) = *link {
                if active[nb] {
                    continue;
                }
                let (ex, ey) = dirs[k];
                let normal_part = -(nx * a[0][comp] + ny * a[1][comp]);
                let tangential_part = tx * a[comp][0] + ty * a[comp][1];
                let g = (ex * nx + ey * ny) * normal_part + (ex * tx + ey * ty) * tangential_part;
                total += dx * g;
            }
        }
        total
    }

    /// Kinetic energy with the face masses the projection is orthogonal for: aperture for
    /// liquid-liquid faces, aperture times the interface fraction for liquid-air faces.
    pub fn weighted_energy(&self) -> f64 {
        let sm = self.sm;
        let n = sm.n();
        let dx2 = sm.dx() * sm.dx();
        let mut e = 0.0;
        for c in 0..n * n {
            let wu = match self.liq.ku[c] {
                1 => sm.au[c],
                2 | 3 => sm.au[c] * self.liq.thu[c],
                _ => 0.0,
            };
            let wv = match self.liq.kv[c] {
                1 => sm.av[c],
                2 | 3 => sm.av[c] * self.liq.thv[c],
                _ => 0.0,
            };
            e += 0.5 * dx2 * (wu * self.u[c] * self.u[c] + wv * self.v[c] * self.v[c]);
        }
        e
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
    /// nodes (zero elsewhere): central differences, or third-order upwind-biased ones when
    /// `upwind` is set; the other velocity component is the 4-face average. A node whose upstream
    /// neighbour is air gets a zero gradient (no information comes from air).
    fn advection(&self, u: &[f64], v: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let sm = self.sm;
        let n = sm.n();
        if self.stokes {
            return (vec![0.0; n * n], vec![0.0; n * n]);
        }
        let dx = sm.dx();
        let (inv2, inv6) = (0.5 / dx, 1.0 / (6.0 * dx));
        let air_u: Vec<bool> = (0..n * n)
            .map(|c| sm.gu.fluid[c] && !self.liq.active_u[c])
            .collect();
        let air_v: Vec<bool> = (0..n * n)
            .map(|c| sm.gv.fluid[c] && !self.liq.active_v[c])
            .collect();
        let d = |q: &[f64], air: &[bool], c: usize, stride: usize, vel: f64| -> f64 {
            if self.weno {
                let dd = |k: usize| (q[k] - q[k - stride]) / dx;
                return if vel > 0.0 {
                    weno5(
                        dd(c - 2 * stride),
                        dd(c - stride),
                        dd(c),
                        dd(c + stride),
                        dd(c + 2 * stride),
                    )
                } else {
                    weno5(
                        dd(c + 3 * stride),
                        dd(c + 2 * stride),
                        dd(c + stride),
                        dd(c),
                        dd(c - stride),
                    )
                };
            }
            if !self.upwind {
                return (q[c + stride] - q[c - stride]) * inv2;
            }
            let (up1, up2, dn, sign) = if vel > 0.0 {
                (c - stride, c - 2 * stride, c + stride, 1.0)
            } else {
                (c + stride, c + 2 * stride, c - stride, -1.0)
            };
            if air[up1] && self.robust_surface {
                return 0.0;
            }
            sign * (q[up2] - 6.0 * q[up1] + 3.0 * q[c] + 2.0 * q[dn]) * inv6
        };
        let mut au = vec![0.0; n * n];
        let mut av = vec![0.0; n * n];
        let margin = if self.weno { 3 } else { 2 };
        for j in margin..n - margin {
            for i in margin..n - margin {
                let c = i + n * j;
                if self.liq.active_u[c] {
                    let vb = 0.25 * (v[c - 1] + v[c] + v[c - 1 + n] + v[c + n]);
                    au[c] = u[c] * d(u, &air_u, c, 1, u[c]) + vb * d(u, &air_u, c, n, vb);
                }
                if self.liq.active_v[c] {
                    let ub = 0.25 * (u[c - n] + u[c + 1 - n] + u[c] + u[c + 1]);
                    av[c] = ub * d(v, &air_v, c, 1, ub) + v[c] * d(v, &air_v, c, n, v[c]);
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

    /// Copy of `q` whose solid-side neighbours of the fluid nodes carry wall-aware ghost values:
    /// the linear extrapolation through the wall point (value `bval`) and the node, or, for a node
    /// closer than a quarter spacing to the wall (where that amplifies by `1 / theta`), through the
    /// wall point and the next node away from it. The strain rate computed from these values sees
    /// the no-slip wall; extrapolating the interior alone does not.
    fn wall_aware(g: &FluidGrid, q: &[f64], bval: impl Fn(f64, f64) -> f64) -> Vec<f64> {
        let n = g.n;
        let mut out = q.to_vec();
        for c in 0..n * n {
            if !g.fluid[c] {
                continue;
            }
            for (k, link) in g.links[c].iter().enumerate() {
                if let Link::Boundary { theta, bx, by } = *link {
                    let nb = match k {
                        0 => c.wrapping_sub(1),
                        1 => c + 1,
                        2 => c.wrapping_sub(n),
                        _ => c + n,
                    };
                    if nb >= n * n || g.fluid[nb] {
                        continue;
                    }
                    let wall = bval(bx, by);
                    let opposite = match &g.links[c][k ^ 1] {
                        Link::Fluid(o) => Some(*o),
                        _ => None,
                    };
                    out[nb] = match opposite {
                        Some(o) if theta < 0.25 => {
                            wall + (wall - q[o]) * (1.0 - theta) / (1.0 + theta)
                        }
                        _ => wall + (wall - q[c]) * (1.0 - theta) / theta,
                    };
                }
            }
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
        self.extend_u(&mut ax);
        self.extend_v(&mut ay);
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
        let prof = std::env::var_os("MILL_PROFILE").is_some();
        let t0 = std::time::Instant::now();
        let sm = self.sm;
        let n = sm.n();
        let nn = n * n;
        let gravity = self.surface.as_ref().map(|s| s.gravity);
        // Predictor: the interface position at the end of the step, advected with the old velocity,
        // defines the topology used by the pressure condition (keeps the coupling second order).
        let mut start_velocity = None;
        if let Some(surf) = &self.surface {
            let vel = surf.cell_velocity(&self.u, &self.v);
            let mut pred = surf.ls.clone();
            if surf.predictor {
                pred.advect(&vel.0, &vel.1, dt);
            }
            self.liq = Liquid::from_level_set(sm, &pred.psi, surf.theta_min);
            start_velocity = Some(vel);
            self.poisson = self.liq.poisson(sm);
            self.helm = None;
            for c in 0..nn {
                if !self.liq.cell[c] {
                    self.p[c] = 0.0;
                }
            }
            let mut tmp = std::mem::take(&mut self.u);
            self.extend_u(&mut tmp);
            self.u = tmp;
            let mut tmp = std::mem::take(&mut self.v);
            self.extend_v(&mut tmp);
            self.v = tmp;
        }
        // ALG2: re-solve the step with the refreshed stress excess until the multiplier settles
        // (one pass for a Newtonian fluid, or in the lagged mode `alg_iterations = 1`).
        let passes = if self.rheo.is_some() {
            self.alg_iterations.max(1)
        } else {
            1
        };
        let saved = (
            self.u.clone(),
            self.v.clone(),
            self.p.clone(),
            self.prev.clone(),
            self.time,
        );
        self.last_alg_passes = 0;
        if self.variational && self.rheo.is_some() {
            let wv = &self.wall_velocity;
            self.visc_op = Some(ViscousOperator::build(
                sm,
                &self.liq.active_u,
                &self.liq.active_v,
                &self.liq.cell,
                &sm.cell_fraction,
                &|x, y| wv(x, y),
            ));
        }
        for pass in 0..passes {
            if pass > 0 {
                (self.u, self.v, self.p, self.prev, self.time) = saved.clone();
            }
            let (au, av) = self.advection(&self.u, &self.v);
            let t_adv = t0.elapsed();
            let (eddy_u, eddy_v) = if self.smagorinsky > 0.0 {
                self.eddy_acceleration()
            } else {
                (vec![0.0; nn], vec![0.0; nn])
            };
            let (gpx, gpy) = self.pressure_gradient(&self.p);
            let (sdu, sdv) = self.stress_divergence();
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
                let (fxs, fys) = if self.swirl != 0.0 {
                    let (i, j) = (c % n, c / n);
                    (
                        -self.swirl * sm.gu.centre(i, j).1,
                        self.swirl * sm.gv.centre(i, j).0,
                    )
                } else {
                    (0.0, 0.0)
                };
                su[c] = (-eu - gpx[c] + gx + eddy_u[c] + sdu[c] + fxs) / self.nu;
                sv[c] = (-ev - gpy[c] + gy + eddy_v[c] + sdv[c] + fys) / self.nu;
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
            let t_helm = t0.elapsed();
            let flux =
                (self.surface.is_some() && self.stress_free_flux).then(|| self.free_surface_flux());
            let wv = &self.wall_velocity;
            if let (true, Some(op), Some(rh), Some((_, hu, hv))) =
                (self.variational, &self.visc_op, &self.rheo, &self.helm)
            {
                let r = rh.r;
                let nu_base = rh.nu_base;
                let sigma_t = sigma * self.nu;
                let mut b = vec![0.0; 2 * nn];
                for c in 0..nn {
                    b[c] = self.nu * (su[c] + sigma * xu[c]);
                    b[nn + c] = self.nu * (sv[c] + sigma * xv[c]);
                }
                let s: Vec<[f64; 3]> = rh.s_c.iter().chain(rh.s_k.iter()).cloned().collect();
                let mut bt = vec![0.0; 2 * nn];
                op.transpose_with_wall(&s, r, &mut bt);
                for (bi, ti) in b.iter_mut().zip(&bt) {
                    *bi -= ti;
                }
                let mask: Vec<bool> = self
                    .liq
                    .active_u
                    .iter()
                    .chain(self.liq.active_v.iter())
                    .copied()
                    .collect();
                let mut x: Vec<f64> = xu.iter().chain(xv.iter()).copied().collect();
                // Warm start from the current velocity.
                x[..nn].copy_from_slice(&self.u);
                x[nn..].copy_from_slice(&self.v);
                let nu = self.nu;
                let zero_wall = |_: f64, _: f64| 0.0;
                let apply = |x: &[f64], out: &mut [f64]| {
                    op.apply(x, out);
                    for i in 0..2 * nn {
                        out[i] = sigma_t * x[i] + r * out[i];
                    }
                    if nu_base > 0.0 {
                        let (lu, _) =
                            Self::gibou_parts(&sm.gu, &x[..nn], &self.liq.active_u, &zero_wall);
                        let (lv, _) =
                            Self::gibou_parts(&sm.gv, &x[nn..], &self.liq.active_v, &zero_wall);
                        for c in 0..nn {
                            out[c] += nu_base * lu[c];
                            out[nn + c] += nu_base * lv[c];
                        }
                    }
                };
                let precond = |res: &[f64], z: &mut [f64]| {
                    let mut zu = vec![0.0; nn];
                    let mut zv = vec![0.0; nn];
                    let ru: Vec<f64> = res[..nn].iter().map(|v| v / nu).collect();
                    let rv: Vec<f64> = res[nn..].iter().map(|v| v / nu).collect();
                    hu.solve(&mut zu, Some(&ru), |_, _| 0.0, 1e-2);
                    hv.solve(&mut zv, Some(&rv), |_, _| 0.0, 1e-2);
                    z[..nn].copy_from_slice(&zu);
                    z[nn..].copy_from_slice(&zv);
                };
                if nu_base > 0.0 {
                    let (_, wu) =
                        Self::gibou_parts(&sm.gu, &self.u, &self.liq.active_u, &|x, y| wv(x, y).0);
                    let (_, wvv) =
                        Self::gibou_parts(&sm.gv, &self.v, &self.liq.active_v, &|x, y| wv(x, y).1);
                    for c in 0..nn {
                        b[c] += nu_base * wu[c];
                        b[nn + c] += nu_base * wvv[c];
                    }
                }
                let (it, res) = pcg(&apply, &precond, &mask, &b, &mut x, 1e-9, 300);
                if std::env::var("ALG_DEBUG").is_ok() {
                    eprintln!("pcg {it} its, residual {res:.2e}");
                }
                xu.copy_from_slice(&x[..nn]);
                xv.copy_from_slice(&x[nn..]);
            } else if let Some((_, hu, hv)) = &self.helm {
                let (fu, fv) = match &flux {
                    Some((a, b)) => (Some(&a[..]), Some(&b[..])),
                    None => (None, None),
                };
                hu.solve_with_flux(
                    &mut xu,
                    Some(&su),
                    fu,
                    |x, y| wv(x, y).0,
                    TOL_VISCOUS * self.tol_scale,
                );
                hv.solve_with_flux(
                    &mut xv,
                    Some(&sv),
                    fv,
                    |x, y| wv(x, y).1,
                    TOL_VISCOUS * self.tol_scale,
                );
            }
            let t_solve = t0.elapsed();
            self.extend_u(&mut xu);
            self.extend_v(&mut xv);
            self.u_star = xu.clone();
            self.v_star = xv.clone();
            self.u = xu;
            self.v = xv;
            let dt_eff = if bdf2 { dt * 2.0 / 3.0 } else { dt };
            let phi = self.project(dt_eff);
            for (p, f) in self.p.iter_mut().zip(&phi) {
                *p += f;
            }
            if prof {
                eprintln!(
                    "STEP adv {:.1} helm-build {:.1} helm-solve {:.1} project {:.1} ms",
                    t_adv.as_secs_f64() * 1e3,
                    (t_helm - t_adv).as_secs_f64() * 1e3,
                    (t_solve - t_helm).as_secs_f64() * 1e3,
                    (t0.elapsed() - t_solve).as_secs_f64() * 1e3
                );
            }
            self.prev = Some(old);
            self.time += dt;
            self.last_alg_passes = pass + 1;
            let change = self.update_rheology();
            if std::env::var("ALG_DEBUG").is_ok() {
                let mut vm = 0.0f64;
                for c in 0..nn {
                    if self.liq.active_u[c] {
                        vm = vm.max(self.u[c].abs());
                    }
                }
                eprintln!("pass {pass} change {change:.3e} vmax {vm:.3e}");
            }
            if pass > 0 && change < self.alg_tolerance {
                break;
            }
        }
        if let (Some(surf), Some(v0)) = (&self.surface, start_velocity) {
            let v1 = surf.cell_velocity(&self.u, &self.v);
            let uc: Vec<f64> = v0.0.iter().zip(&v1.0).map(|(a, b)| 0.5 * (a + b)).collect();
            let vc: Vec<f64> = v0.1.iter().zip(&v1.1).map(|(a, b)| 0.5 * (a + b)).collect();
            let surf = self.surface.as_mut().expect("free surface enabled");
            surf.ls.advect(&uc, &vc, dt);
            surf.ls.reinitialize(surf.reinit_iterations);
            if surf.correct_volume {
                surf.ls.correct_volume(surf.volume0);
            }
        }
    }

    /// Switches to the generalised-Newtonian fluid `law`, treated by the augmented-Lagrangian
    /// method with augmentation `r` (kinematic, 1/s * m^2/s... i.e. twice the implicit viscosity):
    /// the velocity solves use the constant viscosity `r / 2` and the explicit divergence of the
    /// stress excess `S = lam - r d`; `lam` and `d` are updated once per step from the new
    /// velocity (time-lagged ALG2, exact at a steady state).
    pub fn set_rheology(&mut self, law: HerschelBulkley, r: f64) {
        let nn = self.sm.n() * self.sm.n();
        self.nu = 0.5 * r;
        self.helm = None;
        self.rheo = Some(Rheo {
            law,
            r,
            nu_base: 0.0,
            lam_c: vec![[0.0; 3]; nn],
            lam_k: vec![[0.0; 3]; nn],
            s_c: vec![[0.0; 3]; nn],
            s_k: vec![[0.0; 3]; nn],
        });
    }

    /// Variational hybrid: the Newtonian viscosity `nu_base` stays in the wall-exact Gibou solve,
    /// `law` (with `k` = 0 for a Bingham fluid whose plastic viscosity is `nu_base`) is handled by
    /// the augmented Lagrangian with the symmetric operator of `mac/variational.rs`.
    pub fn set_rheology_hybrid(&mut self, law: HerschelBulkley, r: f64, nu_base: f64) {
        self.set_rheology(law, r);
        self.variational = true;
        if let Some(rh) = &mut self.rheo {
            rh.nu_base = nu_base;
        }
        self.nu = nu_base + 0.5 * r;
        self.helm = None;
    }

    /// `-lap x` of the Gibou operator on the active nodes with homogeneous walls (links to
    /// inactive nodes dropped), and the wall-data part `sum bval / theta` of its right-hand side.
    fn gibou_parts(
        g: &FluidGrid,
        x: &[f64],
        active: &[bool],
        bval: &dyn Fn(f64, f64) -> f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let dx2 = g.dx() * g.dx();
        let mut neg_lap = vec![0.0; x.len()];
        let mut wall = vec![0.0; x.len()];
        for c in 0..x.len() {
            if !g.fluid[c] || !active[c] {
                continue;
            }
            let mut s = 0.0;
            let mut w = 0.0;
            for link in &g.links[c] {
                match *link {
                    Link::Fluid(nb) => {
                        if active[nb] {
                            s += x[c] - x[nb];
                        }
                    }
                    Link::Boundary { theta, bx, by } => {
                        s += x[c] / theta;
                        w += bval(bx, by) / theta;
                    }
                }
            }
            neg_lap[c] = s / dx2;
            wall[c] = w / dx2;
        }
        (neg_lap, wall)
    }

    /// ALG2 d- and lambda-steps from the current velocity. Strain rates are compact MAC
    /// differences: `Dxx`, `Dyy` at cell centres, `Dxy` at corners; the missing components of the
    /// tensor at each location are the averages of the four surrounding values.
    #[allow(clippy::needless_range_loop)]
    fn update_rheology(&mut self) -> f64 {
        let Some(mut rh) = self.rheo.take() else {
            return 0.0;
        };
        let n = self.sm.n();
        if let (true, Some(op)) = (self.variational, &self.visc_op) {
            let nn = n * n;
            let x: Vec<f64> = self.u.iter().chain(self.v.iter()).copied().collect();
            let strain = op.strain(&x);
            let (law, r, omega) = (rh.law, rh.r, self.alg_relax);
            let (mut delta, mut scale) = (0.0f64, 0.0f64);
            for p in 0..2 * nn {
                let (lam, s) = if p < nn {
                    (&mut rh.lam_c[p], &mut rh.s_c[p])
                } else {
                    (&mut rh.lam_k[p - nn], &mut rh.s_k[p - nn])
                };
                if !op.valid(p) {
                    *lam = [0.0; 3];
                    *s = [0.0; 3];
                    continue;
                }
                let old = *lam;
                (*lam, *s) = alg2_update(&law, r, omega, &old, &strain[p]);
                for k in 0..3 {
                    delta = delta.max((lam[k] - old[k]).abs());
                    scale = scale.max(lam[k].abs());
                }
            }
            self.rheo = Some(rh);
            return delta / scale.max(1e-300);
        }
        let inv = 1.0 / self.sm.dx();
        let wv = &self.wall_velocity;
        let uw = Self::wall_aware(&self.sm.gu, &self.u, |x, y| wv(x, y).0);
        let vw = Self::wall_aware(&self.sm.gv, &self.v, |x, y| wv(x, y).1);
        let (u, v) = (&uw, &vw);
        let (lu, lv) = (&self.liq.layer_u, &self.liq.layer_v);
        let known = |c: usize, cu: &[usize], cv: &[usize]| {
            cu.iter().all(|&k| lu[c + k] < 255) && cv.iter().all(|&k| lv[c + k] < 255)
        };
        let mut dxx = vec![f64::NAN; n * n];
        let mut dyy = vec![f64::NAN; n * n];
        let mut dxy = vec![f64::NAN; n * n];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let c = i + n * j;
                if known(c, &[0, 1], &[0, n]) {
                    // Deviatoric part: exactly cancels the implicit Laplacian of the velocity
                    // solve for any (also not discretely divergence-free) iterate.
                    let half_diff = 0.5 * ((u[c + 1] - u[c]) - (v[c + n] - v[c])) * inv;
                    dxx[c] = half_diff;
                    dyy[c] = -half_diff;
                }
                if known(c, &[0], &[0]) && lu[c - n] < 255 && lv[c - 1] < 255 {
                    dxy[c] = 0.5 * ((u[c] - u[c - n]) + (v[c] - v[c - 1])) * inv;
                }
            }
        }
        let (law, r) = (rh.law, rh.r);
        let omega = self.alg_relax;
        let (mut delta, mut scale) = (0.0f64, 0.0f64);
        let (mut done_c, mut done_k) = (vec![false; n * n], vec![false; n * n]);
        for j in 2..n - 2 {
            for i in 2..n - 2 {
                let c = i + n * j;
                // Air has no stress: no multiplier in air cells and at corners that touch one.
                let air = |k: usize| self.sm.mesh.pressure_active[k] && !self.liq.cell[k];
                let cell_air = air(c);
                let corner_air = cell_air || air(c - 1) || air(c - n) || air(c - n - 1);
                if cell_air {
                    (rh.lam_c[c], rh.s_c[c]) = ([0.0; 3], [0.0; 3]);
                }
                if corner_air {
                    (rh.lam_k[c], rh.s_k[c]) = ([0.0; 3], [0.0; 3]);
                }
                // Cell centre: normal components of the cell, shear from the four corners.
                let cxy = [dxy[c], dxy[c + 1], dxy[c + n], dxy[c + n + 1]];
                if !cell_air && dxx[c].is_finite() && cxy.iter().all(|x| x.is_finite()) {
                    let d = [dxx[c], 0.25 * cxy.iter().sum::<f64>(), dyy[c]];
                    done_c[c] = self.sm.mesh.pressure_active[c];
                    let old = rh.lam_c[c];
                    (rh.lam_c[c], rh.s_c[c]) = alg2_update(&law, r, omega, &old, &d);
                    for k in [0, 2] {
                        delta = delta.max((rh.lam_c[c][k] - old[k]).abs());
                        scale = scale.max(rh.lam_c[c][k].abs());
                    }
                }
                // Corner: shear of the corner, normal components from the four cells around it.
                let cells = [c - n - 1, c - n, c - 1, c];
                if !corner_air
                    && dxy[c].is_finite()
                    && cells
                        .iter()
                        .all(|&k| dxx[k].is_finite() && dyy[k].is_finite())
                {
                    let mx = 0.25 * cells.iter().map(|&k| dxx[k]).sum::<f64>();
                    let my = 0.25 * cells.iter().map(|&k| dyy[k]).sum::<f64>();
                    let d = [mx, dxy[c], my];
                    done_k[c] = cells.iter().all(|&k| self.sm.mesh.pressure_active[k]);
                    let old = rh.lam_k[c];
                    (rh.lam_k[c], rh.s_k[c]) = alg2_update(&law, r, omega, &old, &d);
                    delta = delta.max((rh.lam_k[c][1] - old[1]).abs());
                    scale = scale.max(rh.lam_k[c][1].abs());
                }
            }
        }
        // Beyond the walls the stress excess is extrapolated from the fluid side (computing it
        // from extrapolated velocities, which carry no wall information, is unreliable).
        let air = |k: usize| self.sm.mesh.pressure_active[k] && !self.liq.cell[k];
        let none = vec![false; n * n];
        let layer_c = ghost_layers(&done_c, &none, n, 2);
        let layer_k = ghost_layers(&done_k, &none, n, 2);
        for comp in 0..3 {
            let mut qc: Vec<f64> = rh.s_c.iter().map(|t| t[comp]).collect();
            extend(&mut qc, &layer_c, n, 2);
            let mut qk: Vec<f64> = rh.s_k.iter().map(|t| t[comp]).collect();
            extend(&mut qk, &layer_k, n, 2);
            for c in n + 1..n * n {
                if !done_c[c] && !air(c) {
                    rh.s_c[c][comp] = qc[c];
                }
                if !done_k[c] && !(air(c) || air(c - 1) || air(c - n) || air(c - n - 1)) {
                    rh.s_k[c][comp] = qk[c];
                }
            }
        }
        self.rheo = Some(rh);
        delta / scale.max(1e-300)
    }

    /// Divergence of the stress excess `S` at the active `u` and `v` nodes (zero without a
    /// rheology): compact differences of the cell-centre normal and corner shear values.
    fn stress_divergence(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.sm.n();
        let mut au = vec![0.0; n * n];
        let mut av = vec![0.0; n * n];
        let Some(rh) = &self.rheo else {
            return (au, av);
        };
        if self.variational {
            return (au, av);
        }
        let inv = 1.0 / self.sm.dx();
        for j in 2..n - 2 {
            for i in 2..n - 2 {
                let c = i + n * j;
                if self.liq.active_u[c] {
                    au[c] =
                        (rh.s_c[c][0] - rh.s_c[c - 1][0] + rh.s_k[c + n][1] - rh.s_k[c][1]) * inv;
                }
                if self.liq.active_v[c] {
                    av[c] =
                        (rh.s_k[c + 1][1] - rh.s_k[c][1] + rh.s_c[c][2] - rh.s_c[c - n][2]) * inv;
                }
            }
        }
        (au, av)
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

    /// Physical torque of the fluid on a drum wall that rotates rigidly at `omega` about the
    /// origin, per unit density: the wall shear of the flow *relative* to the rigid rotation
    /// (zero Dirichlet value at the wall, where the component-wise and the full-stress tractions
    /// coincide). The rigid-rotation part of the component-wise traction, `-2 nu omega A`, is an
    /// artefact of the vector Laplacian and would swamp a partially wetted wall. Only liquid nodes
    /// count. Pressure exerts no torque on a circle.
    pub fn rotating_wall_torque(&self, omega: f64) -> f64 {
        self.rotating_wall_torque_of(omega, f64::NEG_INFINITY)
    }

    /// As `rotating_wall_torque`, counting only boundary points at a radius of at least
    /// `min_radius` (the drum wall, excluding bodies inside it).
    pub fn rotating_wall_torque_of(&self, omega: f64, min_radius: f64) -> f64 {
        let drum = move |x: f64, y: f64| usize::from((x * x + y * y).sqrt() < min_radius);
        let sm = self.sm;
        let n = sm.n();
        let mut ur = self.u.clone();
        let mut vr = self.v.clone();
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                let (_, y) = sm.gu.centre(i, j);
                ur[c] = self.u[c] + omega * y;
                let (x, _) = sm.gv.centre(i, j);
                vr[c] = self.v[c] - omega * x;
            }
        }
        let zero = vec![0.0; n * n];
        let lu = sm.gu.wall_load_masked(
            &ur,
            &zero,
            self.nu,
            |_, _| (0.0, 0.0),
            drum,
            2,
            &self.liq.active_u,
        );
        let lv = sm.gv.wall_load_masked(
            &zero,
            &vr,
            self.nu,
            |_, _| (0.0, 0.0),
            drum,
            2,
            &self.liq.active_v,
        );
        lu[0].torque + lv[0].torque
    }
}
