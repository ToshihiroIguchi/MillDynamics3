//! Reference problems shared by the unit tests and `examples/grid_probe.rs`.

use super::{CircleDomain, PoissonSolver, SolveStats};
use std::f64::consts::PI;

/// Manufactured Neumann Poisson problem `p = cos(pi r^2 / R^2)` on the unit disc; returns the
/// volume-weighted L2 error, the CG statistics and the multigrid level count.
pub fn verify_manufactured(n: usize) -> (f64, SolveStats, usize) {
    let r0 = 1.0;
    let dom = CircleDomain::new(r0, n, 0.05);
    let solver = PoissonSolver::new(&dom);
    let k = PI / (r0 * r0);
    let exact = |x: f64, y: f64| (k * (x * x + y * y)).cos();
    let source = |x: f64, y: f64| {
        let rho = x * x + y * y;
        4.0 * rho * k * k * (k * rho).cos() + 4.0 * k * (k * rho).sin()
    };
    let mut b = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            b[i + n * j] = dom.integrate_cell(i, j, source);
        }
    }
    let mut p = vec![0.0; n * n];
    let stats = solver.solve(&b, &mut p, 1e-12, 200);
    let lv = solver.fine();
    let (mut mean_e, mut cnt) = (0.0, 0.0);
    let mut ex = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let c = i + n * j;
            if lv.diag[c] > 1e-12 {
                let (x, y) = dom.cell_center(i, j);
                ex[c] = exact(x, y);
                mean_e += ex[c];
                cnt += 1.0;
            }
        }
    }
    mean_e /= cnt;
    let (mut e2, mut v) = (0.0, 0.0);
    for j in 0..n {
        for i in 0..n {
            let c = i + n * j;
            if lv.diag[c] > 1e-12 {
                let vol = dom.integrate_cell(i, j, |_, _| 1.0);
                let d = p[c] - (ex[c] - mean_e);
                e2 += d * d * vol;
                v += vol;
            }
        }
    }
    ((e2 / v).sqrt(), stats, solver.level_count())
}

use super::viscous::{rigid_wall_torque_correction, FluidGrid};

/// Steady Stokes Taylor-Couette flow (inner disc `r1` rotating at `omega`, outer wall `r2`
/// fixed) solved componentwise with Dirichlet walls. Returns (max velocity error relative to
/// `omega r1`, inner-wall torque relative error, outer-wall torque relative error).
pub fn verify_couette_stokes(n: usize) -> (f64, f64, f64) {
    let (r1, r2, omega, mu) = (0.25, 0.5, 1.0, 1.0);
    let grid = FluidGrid::new(n, 0.55, |x, y| {
        let r = (x * x + y * y).sqrt();
        (r - r2).max(r1 - r)
    });
    let ut = |r: f64| omega * r1 * r1 * (r2 * r2 / r - r) / (r2 * r2 - r1 * r1);
    let bv = |x: f64, y: f64| {
        let r = (x * x + y * y).sqrt();
        let s = ut(r.clamp(r1, r2)) / r;
        if r < 0.5 * (r1 + r2) {
            (-omega * y, omega * x)
        } else {
            let _ = s;
            (0.0, 0.0)
        }
    };
    let h = grid.helmholtz(0.0);
    let (mut u, mut v) = (vec![0.0; n * n], vec![0.0; n * n]);
    h.solve(&mut u, None, |x, y| bv(x, y).0, 1e-13);
    h.solve(&mut v, None, |x, y| bv(x, y).1, 1e-13);
    let mut emax = 0.0f64;
    for j in 0..n {
        for i in 0..n {
            let c = i + n * j;
            if grid.fluid[c] {
                let (x, y) = grid.centre(i, j);
                let r = (x * x + y * y).sqrt();
                let s = ut(r) / r;
                emax = emax.max((u[c] + s * y).abs()).max((v[c] - s * x).abs());
            }
        }
    }
    let tag = |x: f64, y: f64| usize::from((x * x + y * y).sqrt() > 0.5 * (r1 + r2));
    let loads = grid.wall_load(&u, &v, mu, bv, tag, 2);
    let t_exact = 4.0 * PI * mu * omega * r1 * r1 * r2 * r2 / (r2 * r2 - r1 * r1);
    // Torque on the fluid from the inner wall is positive; on the outer wall negative.
    (
        emax / (omega * r1),
        (loads[0].torque + rigid_wall_torque_correction(mu, omega, PI * r1 * r1, false) - t_exact)
            / t_exact,
        (-loads[1].torque - t_exact) / t_exact,
    )
}

/// Positive zeros of `J1` (McMahon expansion, relative error < 1e-4 from the first zero on).
fn j1_zeros(count: usize) -> Vec<f64> {
    (1..=count)
        .map(|n| {
            let b = (n as f64 + 0.25) * PI;
            b - 0.375 / b + 12.0 / (512.0 * b * b * b)
        })
        .collect()
}

const START_SUBSTEPS: usize = 32;

/// One checkpoint of the spin-up comparison.
#[derive(Clone, Copy, Debug)]
pub struct SpinUpPoint {
    /// `nu t / R^2`.
    pub t_nu: f64,
    /// Relative error of `L(t) / L_rigid`.
    pub l_err: f64,
    /// Relative error of the physical wall torque against `dL/dt`.
    pub torque_err: f64,
}

/// Impulsively started rotation of a circular wall (fluid initially at rest), Stokes regime,
/// BDF2 in time (first steps backward Euler on a refined step). Compares the angular momentum
/// and the wall torque with the Bessel-series solution at `t_nu` in `checkpoints`.
/// `dt_nu` is the time step in `nu / R^2` units.
pub fn verify_spinup(n: usize, dt_nu: f64, checkpoints: &[f64]) -> Vec<SpinUpPoint> {
    let (radius, omega, nu) = (0.5, 1.0, 1.0);
    let grid = FluidGrid::new(n, radius * 1.1, |x, y| (x * x + y * y).sqrt() - radius);
    let dt = dt_nu * radius * radius / nu;
    let bv = |x: f64, y: f64| (-omega * y, omega * x);
    let zeros = j1_zeros(400);
    let m = n * n;
    let mut u = vec![0.0; m];
    let mut v = vec![0.0; m];
    let (mut u_prev, mut v_prev) = (u.clone(), v.clone());
    let l_rigid = {
        let mut s = 0.0;
        for j in 0..n {
            for i in 0..n {
                if grid.fluid[i + n * j] {
                    let (x, y) = grid.centre(i, j);
                    s += omega * (x * x + y * y);
                }
            }
        }
        s * grid.dx() * grid.dx()
    };
    let l_of = |u: &[f64], v: &[f64]| {
        let mut s = 0.0;
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if grid.fluid[c] {
                    let (x, y) = grid.centre(i, j);
                    s += x * v[c] - y * u[c];
                }
            }
        }
        s * grid.dx() * grid.dx()
    };
    let mut out = Vec::new();
    let mut t = 0.0;
    let mut step = 0usize;
    let t_end = checkpoints.iter().cloned().fold(0.0, f64::max) * radius * radius / nu;
    let area = PI * radius * radius;
    while t < t_end - 1e-12 {
        // Start-up: backward-Euler substeps for the first step (impulsive start is singular).
        let (steps, h, bdf2) = if step == 0 {
            (START_SUBSTEPS, dt / START_SUBSTEPS as f64, false)
        } else {
            (1, dt, true)
        };
        for _ in 0..steps {
            let sigma = if bdf2 { 1.5 / (nu * h) } else { 1.0 / (nu * h) };
            let helm = grid.helmholtz(sigma);
            let (un, vn) = (u.clone(), v.clone());
            if bdf2 {
                for c in 0..m {
                    u[c] = (4.0 * un[c] - u_prev[c]) / 3.0;
                    v[c] = (4.0 * vn[c] - v_prev[c]) / 3.0;
                }
            }
            helm.solve(&mut u, None, |x, y| bv(x, y).0, 1e-13);
            helm.solve(&mut v, None, |x, y| bv(x, y).1, 1e-13);
            if bdf2 {
                u_prev = un;
                v_prev = vn;
            }
            t += h;
        }
        // After the start-up interval the previous level is the initial state (at rest).
        step += 1;
        let tn = nu * t / (radius * radius);
        if let Some(&cp) = checkpoints
            .iter()
            .find(|&&c| (c - tn).abs() < 0.5 * dt_nu * 0.999)
        {
            let l_exact = 1.0
                - 8.0
                    * zeros
                        .iter()
                        .map(|l| (-l * l * tn).exp() / (l * l))
                        .sum::<f64>();
            let dl_exact = l_rigid
                * (nu / (radius * radius))
                * 8.0
                * zeros.iter().map(|l| (-l * l * tn).exp()).sum::<f64>();
            let load = grid.wall_load(&u, &v, nu, bv, |_, _| 0, 1)[0];
            let tq = load.torque + rigid_wall_torque_correction(nu, omega, area, true);
            out.push(SpinUpPoint {
                t_nu: cp,
                l_err: l_of(&u, &v) / l_rigid / l_exact - 1.0,
                torque_err: tq / dl_exact - 1.0,
            });
        }
    }
    out
}

use super::flow::{Flow, Mesh};

/// Result of the Navier-Stokes Taylor-Couette hold test.
#[derive(Clone, Copy, Debug)]
pub struct CouetteNs {
    pub steps: usize,
    pub u_err: f64,
    pub inner_torque_err: f64,
    pub outer_torque_err: f64,
    pub divergence: f64,
}

/// Starts from the exact Taylor-Couette state (inner disc `r1` rotating at 1 rad/s, outer wall
/// `r2` fixed, kinematic viscosity `nu`, density 1) and integrates the full Navier-Stokes
/// equations to `t_end`; the exact state is stationary, so any change is discretisation error.
pub fn verify_couette_ns(n: usize, nu: f64, t_end: f64, cfl: f64) -> CouetteNs {
    let (r1, r2, omega) = (0.25, 0.5, 1.0);
    let mesh = Mesh::new(n, 0.55, |x, y| {
        let r = (x * x + y * y).sqrt();
        (r - r2).max(r1 - r)
    });
    let ut = |r: f64| omega * r1 * r1 * (r2 * r2 / r - r) / (r2 * r2 - r1 * r1);
    let exact = |x: f64, y: f64| {
        let r = (x * x + y * y).sqrt();
        let s = ut(r) / r;
        (-s * y, s * x)
    };
    let wall = move |x: f64, y: f64| {
        if (x * x + y * y).sqrt() < 0.5 * (r1 + r2) {
            (-omega * y, omega * x)
        } else {
            (0.0, 0.0)
        }
    };
    let mut flow = Flow::new(&mesh, nu, wall);
    flow.set_velocity(exact);
    let dx = mesh.grid.dx();
    let dt = cfl * dx / (omega * r1);
    let steps = (t_end / dt).ceil() as usize;
    let dt = t_end / steps as f64;
    for _ in 0..steps {
        flow.step(dt);
    }
    let g = &mesh.grid;
    let mut emax = 0.0f64;
    for j in 0..n {
        for i in 0..n {
            let c = i + n * j;
            if g.fluid[c] {
                let (x, y) = g.centre(i, j);
                let (eu, ev) = exact(x, y);
                emax = emax.max((flow.u[c] - eu).abs()).max((flow.v[c] - ev).abs());
            }
        }
    }
    let tag = |x: f64, y: f64| usize::from((x * x + y * y).sqrt() > 0.5 * (r1 + r2));
    let loads = g.wall_load(&flow.u, &flow.v, nu, wall, tag, 2);
    let t_exact = 4.0 * PI * nu * omega * r1 * r1 * r2 * r2 / (r2 * r2 - r1 * r1);
    CouetteNs {
        steps,
        u_err: emax / (omega * r1),
        inner_torque_err: (loads[0].torque
            + rigid_wall_torque_correction(nu, omega, PI * r1 * r1, false)
            - t_exact)
            / t_exact,
        outer_torque_err: (-loads[1].torque - t_exact) / t_exact,
        divergence: flow.max_face_divergence(),
    }
}

use super::staggered::{StaggeredFlow, StaggeredMesh};

/// Same hold test as [`verify_couette_ns`] on the staggered face-velocity scheme.
pub fn verify_couette_mac(n: usize, nu: f64, t_end: f64, cfl: f64) -> CouetteNs {
    let (r1, r2, omega) = (0.25, 0.5, 1.0);
    let sm = StaggeredMesh::new(n, 0.55, |x, y| {
        let r = (x * x + y * y).sqrt();
        (r - r2).max(r1 - r)
    });
    let ut = |r: f64| omega * r1 * r1 * (r2 * r2 / r - r) / (r2 * r2 - r1 * r1);
    let exact = |x: f64, y: f64| {
        let r = (x * x + y * y).sqrt();
        let s = ut(r) / r;
        (-s * y, s * x)
    };
    let wall = move |x: f64, y: f64| {
        if (x * x + y * y).sqrt() < 0.5 * (r1 + r2) {
            (-omega * y, omega * x)
        } else {
            (0.0, 0.0)
        }
    };
    let mut flow = StaggeredFlow::new(&sm, nu, wall);
    flow.set_velocity(exact);
    let dx = sm.dx();
    let dt = cfl * dx / (omega * r1);
    let steps = (t_end / dt).ceil() as usize;
    let dt = t_end / steps as f64;
    let mut divergence = 0.0f64;
    for _ in 0..steps {
        flow.step(dt);
        divergence = divergence.max(flow.projected_divergence);
    }
    let mut emax = 0.0f64;
    for j in 0..n {
        for i in 0..n {
            let c = i + n * j;
            if sm.gu.fluid[c] {
                let (x, y) = sm.gu.centre(i, j);
                emax = emax.max((flow.u[c] - exact(x, y).0).abs());
            }
            if sm.gv.fluid[c] {
                let (x, y) = sm.gv.centre(i, j);
                emax = emax.max((flow.v[c] - exact(x, y).1).abs());
            }
        }
    }
    let tag = |x: f64, y: f64| usize::from((x * x + y * y).sqrt() > 0.5 * (r1 + r2));
    let loads = flow.wall_load(tag, 2);
    let t_exact = 4.0 * PI * nu * omega * r1 * r1 * r2 * r2 / (r2 * r2 - r1 * r1);
    CouetteNs {
        steps,
        u_err: emax / (omega * r1),
        inner_torque_err: (loads[0].torque
            + rigid_wall_torque_correction(nu, omega, PI * r1 * r1, false)
            - t_exact)
            / t_exact,
        outer_torque_err: (-loads[1].torque - t_exact) / t_exact,
        divergence,
    }
}

/// Diagnostic: error history of the staggered hold test (`(step, max u err, max u err at
/// ghost faces)` at the listed step counts, step 0 = right after `set_velocity`).
pub fn diagnose_couette_mac(n: usize, nu: f64, cfl: f64, at: &[usize]) -> Vec<(usize, f64, f64)> {
    let (r1, r2, omega) = (0.25, 0.5, 1.0);
    let sm = StaggeredMesh::new(n, 0.55, |x, y| {
        let r = (x * x + y * y).sqrt();
        (r - r2).max(r1 - r)
    });
    let ut = |r: f64| omega * r1 * r1 * (r2 * r2 / r - r) / (r2 * r2 - r1 * r1);
    let exact = |x: f64, y: f64| {
        let r = (x * x + y * y).sqrt();
        let s = ut(r) / r;
        (-s * y, s * x)
    };
    let wall = move |x: f64, y: f64| {
        if (x * x + y * y).sqrt() < 0.5 * (r1 + r2) {
            (-omega * y, omega * x)
        } else {
            (0.0, 0.0)
        }
    };
    let mut flow = StaggeredFlow::new(&sm, nu, wall);
    flow.set_velocity(exact);
    let dt = cfl * sm.dx() / (omega * r1);
    let measure = |flow: &StaggeredFlow| {
        let (mut fl, mut gh) = (0.0f64, 0.0f64);
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if sm.au[c] > 0.0 {
                    let (x, y) = sm.gu.centre(i, j);
                    let e = (flow.u[c] - exact(x, y).0).abs();
                    if sm.gu.fluid[c] {
                        fl = fl.max(e);
                    } else {
                        gh = gh.max(e);
                    }
                }
                if sm.av[c] > 0.0 {
                    let (x, y) = sm.gv.centre(i, j);
                    let e = (flow.v[c] - exact(x, y).1).abs();
                    if sm.gv.fluid[c] {
                        fl = fl.max(e);
                    } else {
                        gh = gh.max(e);
                    }
                }
            }
        }
        (fl / (omega * r1), gh / (omega * r1))
    };
    let mut out = Vec::new();
    let last = *at.last().unwrap();
    for s in 0..=last {
        if at.contains(&s) {
            let (a, b) = measure(&flow);
            out.push((s, a, b));
        }
        flow.step(dt);
    }
    out
}
