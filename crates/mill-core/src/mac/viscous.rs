//! Cell-centred Dirichlet Helmholtz/Poisson on an arbitrary (level-set) fluid region with
//! second-order accurate solution (Gibou et al. 2002, symmetric `1/theta` form) and the discrete
//! wall reaction used for force/torque.
//!
//! For constant viscosity the incompressible viscous term is `mu * Laplacian` acting on each
//! Cartesian velocity component, so one scalar solver serves `u` and `v`.

use super::multigrid::{PoissonSolver, SolveStats};

const THETA_MIN: f64 = 1e-3;
const LEFT: usize = 0;
const RIGHT: usize = 1;
const DOWN: usize = 2;
const UP: usize = 3;

/// Link from a fluid cell to one neighbour cell centre.
#[derive(Clone, Copy, Debug)]
pub enum Link {
    Fluid(usize),
    /// Boundary crossed at fraction `theta` of the centre-to-centre distance, at `(bx, by)`.
    Boundary {
        theta: f64,
        bx: f64,
        by: f64,
    },
}

#[derive(Clone, Debug)]
pub struct FluidGrid {
    pub n: usize,
    pub half: f64,
    /// Lattice offset in units of `dx`: node `(i, j)` sits at `(-half + (i + ox) dx, -half + (j + oy) dx)`.
    pub ox: f64,
    pub oy: f64,
    pub fluid: Vec<bool>,
    pub links: Vec<[Link; 4]>,
}

/// Wall reaction on the fluid, summed over the links whose boundary point has a given tag.
#[derive(Clone, Copy, Debug, Default)]
pub struct WallLoad {
    pub fx: f64,
    pub fy: f64,
    /// Torque about the origin exerted on the fluid by the wall.
    pub torque: f64,
}

impl FluidGrid {
    /// Fluid is where `phi < 0`. The region must stay away from the grid border.
    pub fn new(n: usize, half: f64, phi: impl Fn(f64, f64) -> f64) -> Self {
        Self::with_offset(n, half, 0.5, 0.5, phi)
    }

    /// Same on a shifted lattice (`ox = oy = 0.5` are the cell centres; `(0, 0.5)` and `(0.5, 0)`
    /// are the staggered `u` and `v` face lattices).
    pub fn with_offset(
        n: usize,
        half: f64,
        ox: f64,
        oy: f64,
        phi: impl Fn(f64, f64) -> f64,
    ) -> Self {
        let dx = 2.0 * half / n as f64;
        let centre =
            |i: usize, j: usize| (-half + (i as f64 + ox) * dx, -half + (j as f64 + oy) * dx);
        let mut fluid = vec![false; n * n];
        for j in 0..n {
            for i in 0..n {
                let (x, y) = centre(i, j);
                fluid[i + n * j] = phi(x, y) < 0.0;
            }
        }
        let dummy = Link::Boundary {
            theta: 1.0,
            bx: 0.0,
            by: 0.0,
        };
        let mut links = vec![[dummy; 4]; n * n];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if !fluid[c] {
                    continue;
                }
                let (x, y) = centre(i, j);
                let nbrs: [(isize, isize); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];
                for (k, (di, dj)) in nbrs.iter().enumerate() {
                    let (ni, nj) = (i as isize + di, j as isize + dj);
                    let inside = ni >= 0 && nj >= 0 && (ni as usize) < n && (nj as usize) < n;
                    if inside && fluid[ni as usize + n * nj as usize] {
                        links[c][k] = Link::Fluid(ni as usize + n * nj as usize);
                        continue;
                    }
                    let (xn, yn) = (x + *di as f64 * dx, y + *dj as f64 * dx);
                    let (mut lo, mut hi) = (0.0, 1.0);
                    for _ in 0..60 {
                        let mid = 0.5 * (lo + hi);
                        if phi(x + mid * (xn - x), y + mid * (yn - y)) < 0.0 {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    let t = (0.5 * (lo + hi)).max(THETA_MIN);
                    links[c][k] = Link::Boundary {
                        theta: t,
                        bx: x + t * (xn - x),
                        by: y + t * (yn - y),
                    };
                }
            }
        }
        Self {
            n,
            half,
            ox,
            oy,
            fluid,
            links,
        }
    }

    pub fn dx(&self) -> f64 {
        2.0 * self.half / self.n as f64
    }

    pub fn centre(&self, i: usize, j: usize) -> (f64, f64) {
        let dx = self.dx();
        (
            -self.half + (i as f64 + self.ox) * dx,
            -self.half + (j as f64 + self.oy) * dx,
        )
    }

    /// Operator `(sigma - Laplacian) x`, `sigma >= 0` in 1/m^2 (`1 / (nu dt)`); `sigma = 0` is the
    /// Poisson problem.
    pub fn helmholtz(&self, sigma: f64) -> Helmholtz<'_> {
        self.helmholtz_masked(sigma, &self.fluid)
    }

    /// Same on the nodes with `active[c]` only: links to wall boundaries stay Dirichlet, links to
    /// inactive neighbour nodes are dropped (zero normal derivative, e.g. a free surface).
    pub fn helmholtz_masked(&self, sigma: f64, active: &[bool]) -> Helmholtz<'_> {
        let n = self.n;
        let dx2 = self.dx() * self.dx();
        let mut wx = vec![0.0; (n + 1) * n];
        let mut wy = vec![0.0; n * (n + 1)];
        let mut extra = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if !self.fluid[c] || !active[c] {
                    continue;
                }
                extra[c] = sigma * dx2;
                for (k, link) in self.links[c].iter().enumerate() {
                    match link {
                        Link::Fluid(nb) => {
                            if !active[*nb] {
                                continue;
                            }
                            match k {
                                RIGHT => wx[i + 1 + (n + 1) * j] = 1.0,
                                UP => wy[i + n * (j + 1)] = 1.0,
                                LEFT => wx[i + (n + 1) * j] = 1.0,
                                _ => wy[i + n * j] = 1.0,
                            }
                        }
                        Link::Boundary { theta, .. } => extra[c] += 1.0 / theta,
                    }
                }
            }
        }
        let _ = DOWN;
        Helmholtz {
            grid: self,
            sigma,
            active: active
                .iter()
                .zip(&self.fluid)
                .map(|(&a, &f)| a && f)
                .collect(),
            solver: PoissonSolver::from_weights(n, wx, wy, extra),
        }
    }

    /// Reaction of the walls on the fluid for the field `(u, v)` with boundary velocity `bv`;
    /// `tag(bx, by)` selects which boundary a link belongs to. Force = `mu (u_b - u_c) / theta`.
    ///
    /// This is the traction of the *component-wise* operator, `mu grad(u) . n`. The physical
    /// traction `mu (grad u + grad u^T) . n` differs by the `grad u^T . n` term, which for a rigid
    /// wall rotating at `Omega` integrates to a torque only: see [`rigid_wall_torque_correction`].
    /// Forces need no correction.
    pub fn wall_load(
        &self,
        u: &[f64],
        v: &[f64],
        mu: f64,
        bv: impl Fn(f64, f64) -> (f64, f64),
        tag: impl Fn(f64, f64) -> usize,
        tags: usize,
    ) -> Vec<WallLoad> {
        let mut out = vec![WallLoad::default(); tags];
        for c in 0..self.n * self.n {
            if !self.fluid[c] {
                continue;
            }
            for link in &self.links[c] {
                if let Link::Boundary { theta, bx, by } = *link {
                    let (ub, vb) = bv(bx, by);
                    let (fx, fy) = (mu * (ub - u[c]) / theta, mu * (vb - v[c]) / theta);
                    let w = &mut out[tag(bx, by)];
                    w.fx += fx;
                    w.fy += fy;
                    w.torque += bx * fy - by * fx;
                }
            }
        }
        out
    }
}

pub struct Helmholtz<'a> {
    grid: &'a FluidGrid,
    sigma: f64,
    active: Vec<bool>,
    solver: PoissonSolver,
}

impl Helmholtz<'_> {
    /// Solves `(sigma - Laplacian) x = sigma x_old` (or `source` added) with Dirichlet values
    /// `bval(bx, by)`; `x` is both the initial guess slot and the result (`x_old` is read first).
    pub fn solve(
        &self,
        x: &mut [f64],
        source: Option<&[f64]>,
        bval: impl Fn(f64, f64) -> f64,
        tol: f64,
    ) -> SolveStats {
        let g = self.grid;
        let dx2 = g.dx() * g.dx();
        let mut b = vec![0.0; x.len()];
        for c in 0..x.len() {
            if !self.active[c] {
                continue;
            }
            b[c] = self.sigma * dx2 * x[c] + source.map_or(0.0, |s| dx2 * s[c]);
            for link in &g.links[c] {
                if let Link::Boundary { theta, bx, by } = *link {
                    b[c] += bval(bx, by) / theta;
                }
            }
        }
        let mut sol = vec![0.0; x.len()];
        let stats = self.solver.solve(&b, &mut sol, tol, 500);
        for c in 0..x.len() {
            if self.active[c] {
                x[c] = sol[c];
            }
        }
        stats
    }
}

/// Torque to add to the component-wise wall torque (on the fluid) to obtain the physical torque
/// of the full viscous stress, for a rigid wall rotating at `omega` (rad/s) and enclosing `area`
/// (m^2): `+2 mu omega area` for a body surrounded by fluid, `-2 mu omega area` for a drum wall
/// enclosing the fluid. (Check: fluid in rigid rotation inside its drum then exerts no torque.)
pub fn rigid_wall_torque_correction(mu: f64, omega: f64, area: f64, fluid_inside: bool) -> f64 {
    let sign = if fluid_inside { -1.0 } else { 1.0 };
    sign * 2.0 * mu * omega * area
}
