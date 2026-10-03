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

use super::bodies::{disc_load, Body, BodyFlow, BodyLoad};
use super::staggered::{Disc, FlowState, StaggeredFlow, StaggeredMesh};

/// Same hold test as [`verify_couette_ns`] on the staggered face-velocity scheme.
pub fn verify_couette_mac(n: usize, nu: f64, t_end: f64, cfl: f64) -> CouetteNs {
    verify_couette_mac_scheme(n, nu, t_end, cfl, false)
}

/// Hold test with a choice of advection scheme (`upwind`: third-order upwind-biased).
pub fn verify_couette_mac_scheme(
    n: usize,
    nu: f64,
    t_end: f64,
    cfl: f64,
    upwind: bool,
) -> CouetteNs {
    verify_couette_mac_flow(n, nu, t_end, cfl, &|f| {
        f.upwind = upwind;
        f.weno = upwind;
    })
}

/// Hold test with the flow configured by `configure`.
pub fn verify_couette_mac_flow(
    n: usize,
    nu: f64,
    t_end: f64,
    cfl: f64,
    configure: &dyn Fn(&mut StaggeredFlow),
) -> CouetteNs {
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
    configure(&mut flow);
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

/// Impulsive spin-up of a circular wall with the full staggered Navier-Stokes stepper (the flow
/// is azimuthal, so the Stokes Bessel solution is exact for NS too). Same checkpoints and
/// metrics as [`verify_spinup`]; the first step is covered by `START_SUBSTEPS` equal substeps.
pub fn verify_spinup_mac(n: usize, dt_nu: f64, checkpoints: &[f64]) -> Vec<SpinUpPoint> {
    let (radius, omega, nu) = (0.5, 1.0, 1.0);
    let sm = StaggeredMesh::new(n, radius * 1.1, |x, y| (x * x + y * y).sqrt() - radius);
    let dt = dt_nu * radius * radius / nu;
    let mut flow = StaggeredFlow::new(&sm, nu, |x, y| (-omega * y, omega * x));
    let zeros = j1_zeros(400);
    let dx2 = sm.dx() * sm.dx();
    let moment = |f: &dyn Fn(f64, f64, usize) -> f64| {
        let mut s = 0.0;
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if sm.gu.fluid[c] {
                    let (x, y) = sm.gu.centre(i, j);
                    s += -y * f(x, y, c);
                }
                if sm.gv.fluid[c] {
                    let (x, y) = sm.gv.centre(i, j);
                    s += x * f(x, y, c + n * n);
                }
            }
        }
        s * dx2
    };
    // L_rigid: u = -omega y, v = omega x.
    let l_rigid = moment(&|x, y, c| if c < n * n { -omega * y } else { omega * x });
    let area = PI * radius * radius;
    let t_end = checkpoints.iter().cloned().fold(0.0, f64::max) * radius * radius / nu;
    let mut out = Vec::new();
    let mut t = 0.0;
    let mut first = true;
    while t < t_end - 1e-12 {
        let (steps, h) = if first {
            (START_SUBSTEPS, dt / START_SUBSTEPS as f64)
        } else {
            (1, dt)
        };
        first = false;
        for _ in 0..steps {
            flow.step(h);
            t += h;
        }
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
            let l = moment(&|_, _, c| {
                if c < n * n {
                    flow.u[c]
                } else {
                    flow.v[c - n * n]
                }
            });
            let load = flow.wall_load(|_, _| 0, 1)[0];
            let tq = load.torque + rigid_wall_torque_correction(nu, omega, area, true);
            out.push(SpinUpPoint {
                t_nu: cp,
                l_err: l / l_rigid / l_exact - 1.0,
                torque_err: tq / dl_exact - 1.0,
            });
        }
    }
    out
}

use super::levelset::LevelSet;

/// Result of the rigid-rotation level-set test.
#[derive(Clone, Copy, Debug)]
pub struct LevelSetRotation {
    pub steps: usize,
    /// `sum |indicator - exact indicator| dx^2 / area` after the rotation.
    pub shape_l1: f64,
    /// `V_end / V_0 - 1` without the volume correction.
    pub volume_drift: f64,
}

/// Rotates a disc (or the slotted Zalesak disc) rigidly about the origin for `revs` revolutions
/// (period 1) with WENO5/TVD-RK3 advection and reinitialisation every 2 steps; compares the final
/// shape with the initial one.
pub fn verify_levelset_rotation(
    n: usize,
    slotted: bool,
    revs: f64,
    reinit_every: usize,
) -> LevelSetRotation {
    let half = 0.55;
    let disc = |x: f64, y: f64| (x * x + (y - 0.25) * (y - 0.25)).sqrt() - 0.15;
    let shape = move |x: f64, y: f64| {
        if slotted {
            let slot = (x.abs() - 0.025).max((y - 0.225).abs() - 0.125);
            disc(x, y).max(-slot)
        } else {
            disc(x, y)
        }
    };
    let mut ls = LevelSet::new(n, half, shape);
    ls.reinitialize(20);
    let v0 = ls.volume();
    let omega = 2.0 * PI;
    let dx = ls.dx;
    let mut uc = vec![0.0; n * n];
    let mut vc = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = ls.centre(i, j);
            uc[i + n * j] = -omega * y;
            vc[i + n * j] = omega * x;
        }
    }
    let dt0 = 0.5 * dx / (omega * 0.5);
    let steps = (revs / dt0).ceil() as usize;
    let dt = revs / steps as f64;
    for s in 0..steps {
        ls.advect(&uc, &vc, dt);
        if reinit_every > 0 && s % reinit_every == reinit_every - 1 {
            ls.reinitialize(2);
        }
    }
    let mut diff = 0.0;
    let mut area = 0.0;
    for j in 0..n {
        for i in 0..n {
            let (x, y) = ls.centre(i, j);
            let exact = shape(x, y) < 0.0;
            diff += f64::from(u8::from(exact != (ls.psi[i + n * j] < 0.0)));
            area += f64::from(u8::from(exact));
        }
    }
    LevelSetRotation {
        steps,
        shape_l1: diff / area,
        volume_drift: ls.volume() / v0 - 1.0,
    }
}

use super::levelset::LevelSet as FsLevelSet;

/// Result of the still-pool test.
#[derive(Clone, Copy, Debug)]
pub struct StillPool {
    pub steps: usize,
    /// `max |velocity| / sqrt(g D)` over active nodes at the end.
    pub spurious: f64,
    /// `max |p - p_hydrostatic| / (g depth)` over liquid cells at the end.
    pub pressure_err: f64,
    /// Liquid area drift of the level set (`V / V0 - 1`).
    pub volume_drift: f64,
    /// Highest interface position error `max |y_interface - level|`, in cells.
    pub level_err_cells: f64,
}

/// A circular drum (radius 0.5) half filled with still liquid (`psi = y - level`), gravity
/// downward; integrates `t_end` seconds and checks that nothing moves and the pressure stays
/// hydrostatic.
pub fn verify_still_pool(n: usize, level: f64, nu: f64, t_end: f64) -> StillPool {
    let radius = 0.5;
    let g = 9.81;
    let sm = StaggeredMesh::new(n, 0.55, |x, y| (x * x + y * y).sqrt() - radius);
    let mut flow = StaggeredFlow::new(&sm, nu, |_, _| (0.0, 0.0));
    let ls = FsLevelSet::new(n, 0.55, |_, y| y - level);
    flow.enable_free_surface(ls, (0.0, -g));
    let dx = sm.dx();
    let dt0 = 0.3 * dx / (g * 2.0 * radius).sqrt();
    let steps = (t_end / dt0).ceil() as usize;
    let dt = t_end / steps as f64;
    for _ in 0..steps {
        flow.step(dt);
    }
    let mut vmax = 0.0f64;
    for c in 0..n * n {
        if flow.liq.active_u[c] {
            vmax = vmax.max(flow.u[c].abs());
        }
        if flow.liq.active_v[c] {
            vmax = vmax.max(flow.v[c].abs());
        }
    }
    let surf = flow.surface.as_ref().expect("free surface enabled");
    let depth = level + radius;
    let mut perr = 0.0f64;
    for j in 0..n {
        for i in 0..n {
            let c = i + n * j;
            if flow.liq.cell[c] {
                let (_, y) = surf.ls.centre(i, j);
                perr = perr.max((flow.p[c] - g * (level - y)).abs());
            }
        }
    }
    // Interface height in the middle column: linear interpolation of psi.
    let mut lerr = 0.0f64;
    for i in n / 2 - 2..n / 2 + 2 {
        for j in 0..n - 1 {
            let (a, b) = (surf.ls.psi[i + n * j], surf.ls.psi[i + n * (j + 1)]);
            if a < 0.0 && b >= 0.0 {
                let (_, ya) = surf.ls.centre(i, j);
                let yi = ya + dx * a / (a - b);
                lerr = lerr.max((yi - level).abs() / dx);
            }
        }
    }
    StillPool {
        steps,
        spurious: vmax / (g * 2.0 * radius).sqrt(),
        pressure_err: perr / (g * depth),
        volume_drift: surf.ls.volume() / surf.volume0 - 1.0,
        level_err_cells: lerr,
    }
}

/// Result of a sloshing run.
#[derive(Clone, Copy, Debug)]
pub struct Sloshing {
    pub steps: usize,
    /// Measured angular frequency (rad/s) from zero crossings of the first-mode amplitude.
    pub omega: f64,
    /// Amplitude decay rate (1/s) from successive extrema (0 if fewer than two found).
    pub damping: f64,
    pub volume_drift: f64,
}

/// Interface height of the column `i` by linear interpolation of `psi` (highest crossing).
fn column_height(ls: &FsLevelSet, i: usize) -> Option<f64> {
    let n = ls.n;
    let mut best = None;
    for j in 0..n - 1 {
        let (a, b) = (ls.psi[i + n * j], ls.psi[i + n * (j + 1)]);
        if a < 0.0 && b >= 0.0 {
            let (_, ya) = ls.centre(i, j);
            best = Some(ya + ls.dx * a / (a - b));
        }
    }
    best
}

/// Extracts `(omega, damping)` from a time series `(t, a)` by zero crossings and extrema.
fn oscillation(series: &[(f64, f64)]) -> (f64, f64) {
    let mut crossings = Vec::new();
    for w in series.windows(2) {
        if w[0].1 * w[1].1 < 0.0 {
            crossings.push(w[0].0 + (w[1].0 - w[0].0) * w[0].1 / (w[0].1 - w[1].1));
        }
    }
    let omega = if crossings.len() >= 2 {
        PI * (crossings.len() - 1) as f64 / (crossings[crossings.len() - 1] - crossings[0])
    } else {
        f64::NAN
    };
    // Peak amplitudes between crossings.
    let mut peaks: Vec<(f64, f64)> = Vec::new();
    let mut best = (0.0, 0.0f64);
    let mut next = 0;
    for &(t, a) in series {
        while next < crossings.len() && t > crossings[next] {
            if best.1 > 0.0 {
                peaks.push(best);
            }
            best = (0.0, 0.0);
            next += 1;
        }
        if a.abs() > best.1 {
            best = (t, a.abs());
        }
    }
    let damping = if peaks.len() >= 2 {
        let (a, b) = (peaks[0], peaks[peaks.len() - 1]);
        (a.1 / b.1).ln() / (b.0 - a.0)
    } else {
        0.0
    };
    (omega, damping)
}

/// Small-amplitude first-mode sloshing in a rectangular tank `[-width/2, width/2] x [-width/2,
/// -width/2 + height]` filled to `depth` above the floor, started from a cosine interface of
/// amplitude `amp` at rest. Returns the measured frequency and decay.
pub fn verify_sloshing_rect(
    n: usize,
    width: f64,
    depth: f64,
    amp: f64,
    nu: f64,
    t_end: f64,
    upwind: bool,
) -> Sloshing {
    verify_sloshing_rect_with(n, width, depth, amp, nu, t_end, &|f| f.upwind = upwind)
}

/// [`verify_sloshing_rect`] with the flow configured by `configure`.
pub fn verify_sloshing_rect_with(
    n: usize,
    width: f64,
    depth: f64,
    amp: f64,
    nu: f64,
    t_end: f64,
    configure: &dyn Fn(&mut StaggeredFlow),
) -> Sloshing {
    let g = 9.81;
    let floor = -0.5 * width;
    let top = floor + 0.95;
    let sm = StaggeredMesh::new(n, 0.55, |x, y| {
        (x.abs() - 0.5 * width).max(floor - y).max(y - top)
    });
    let mut flow = StaggeredFlow::new(&sm, nu, |_, _| (0.0, 0.0));
    configure(&mut flow);
    let k = PI / width;
    let level = floor + depth;
    let ls = FsLevelSet::new(n, 0.55, |x, y| {
        y - level - amp * (k * (x + 0.5 * width)).cos()
    });
    flow.enable_free_surface(ls, (0.0, -g));
    let dx = sm.dx();
    let dt0 = 0.3 * dx / (g * depth).sqrt();
    let steps = (t_end / dt0).ceil() as usize;
    let dt = t_end / steps as f64;
    let mut series = Vec::with_capacity(steps);
    let columns: Vec<usize> = (0..n)
        .filter(|&i| {
            let (x, _) = flow.surface.as_ref().unwrap().ls.centre(i, 0);
            x.abs() < 0.5 * width - 2.0 * dx
        })
        .collect();
    for s in 0..steps {
        flow.step(dt);
        let ls = &flow.surface.as_ref().unwrap().ls;
        // First-mode amplitude: projection of the interface elevation on cos(k (x + w/2)).
        let mut a = 0.0;
        for &i in &columns {
            if let Some(h) = column_height(ls, i) {
                let (x, _) = ls.centre(i, 0);
                a += (h - level) * (k * (x + 0.5 * width)).cos() * dx;
            }
        }
        series.push(((s + 1) as f64 * dt, a * 2.0 / width));
    }
    let (omega, damping) = oscillation(&series);
    let surf = flow.surface.as_ref().unwrap();
    Sloshing {
        steps,
        omega,
        damping,
        volume_drift: surf.ls.volume() / surf.volume0 - 1.0,
    }
}

/// Small-amplitude sloshing of a half-full circular drum (radius 0.5) started from a tilted flat
/// interface `y = amp x / R` at rest; the tilt amplitude `sum eta x / sum x^2 * R` is followed in
/// time. The reference frequency is `sqrt(K g / R)` with `K R` from
/// [`super::reference::half_disc_sloshing`].
pub fn verify_sloshing_circle(n: usize, amp: f64, nu: f64, t_end: f64, upwind: bool) -> Sloshing {
    verify_sloshing_circle_with(n, amp, nu, t_end, &|f| f.upwind = upwind)
}

/// [`verify_sloshing_circle`] with the flow configured by `configure`.
pub fn verify_sloshing_circle_with(
    n: usize,
    amp: f64,
    nu: f64,
    t_end: f64,
    configure: &dyn Fn(&mut StaggeredFlow),
) -> Sloshing {
    let radius = 0.5;
    let g = 9.81;
    let sm = StaggeredMesh::new(n, 0.55, |x, y| (x * x + y * y).sqrt() - radius);
    let mut flow = StaggeredFlow::new(&sm, nu, |_, _| (0.0, 0.0));
    configure(&mut flow);
    let ls = FsLevelSet::new(n, 0.55, |x, y| y - amp * x / radius);
    flow.enable_free_surface(ls, (0.0, -g));
    let dx = sm.dx();
    let dt0 = 0.3 * dx / (g * radius).sqrt();
    let steps = (t_end / dt0).ceil() as usize;
    let dt = t_end / steps as f64;
    let columns: Vec<usize> = (0..n)
        .filter(|&i| {
            let (x, _) = flow.surface.as_ref().unwrap().ls.centre(i, 0);
            x.abs() < radius - 2.0 * dx
        })
        .collect();
    let mut series = Vec::with_capacity(steps);
    for s in 0..steps {
        flow.step(dt);
        let ls = &flow.surface.as_ref().unwrap().ls;
        let (mut num, mut den) = (0.0, 0.0);
        for &i in &columns {
            if let Some(h) = column_height(ls, i) {
                let (x, _) = ls.centre(i, 0);
                num += h * x;
                den += x * x;
            }
        }
        series.push(((s + 1) as f64 * dt, num / den * radius));
    }
    let (omega, damping) = oscillation(&series);
    let surf = flow.surface.as_ref().unwrap();
    Sloshing {
        steps,
        omega,
        damping,
        volume_drift: surf.ls.volume() / surf.volume0 - 1.0,
    }
}

/// Result of the dam-break run: front position history and conservation.
#[derive(Clone, Debug)]
pub struct DamBreak {
    pub steps: usize,
    /// `(T, Z)` samples with `T = t sqrt(g / a)`, `Z = (x_front - x_wall) / a`.
    pub front: Vec<(f64, f64)>,
    pub volume_drift: f64,
    pub max_speed: f64,
}

/// Collapse of a square liquid column of side `a` against the left wall of a rectangular tank
/// `1.0 x 0.95` (liquid released at `t = 0`); samples the surge front (rightmost liquid cell in
/// the bottom row) at the dimensionless times `ts`.
pub fn verify_dam_break(n: usize, a: f64, nu: f64, ts: &[f64], upwind: bool) -> DamBreak {
    let g = 9.81;
    let (left, floor) = (-0.5, -0.5);
    let top = floor + 0.95;
    let sm = StaggeredMesh::new(n, 0.55, |x, y| {
        (x - 0.5).max(left - x).max(floor - y).max(y - top)
    });
    let mut flow = StaggeredFlow::new(&sm, nu, |_, _| (0.0, 0.0));
    flow.upwind = upwind;
    let ls = FsLevelSet::new(n, 0.55, |x, y| (x - (left + a)).max(y - (floor + a)));
    flow.enable_free_surface(ls, (0.0, -g));
    flow.surface.as_mut().expect("surface").ls.reinitialize(20);
    let dx = sm.dx();
    let dt = 0.2 * dx / (2.0 * (g * a).sqrt());
    let t_end = ts.iter().cloned().fold(0.0, f64::max) * (a / g).sqrt();
    let steps = (t_end / dt).ceil() as usize;
    let mut front = Vec::new();
    let mut next = 0;
    let mut max_speed = 0.0f64;
    // Bottom row: first row whose centre is above the floor.
    let jrow = ((floor - (-0.55)) / dx).floor() as usize;
    for s in 0..steps {
        flow.step(dt);
        let t = (s + 1) as f64 * dt;
        let tt = t * (g / a).sqrt();
        while next < ts.len() && tt >= ts[next] {
            let ls = &flow.surface.as_ref().expect("surface").ls;
            let mut xf = left;
            for i in 0..n {
                if ls.psi[i + n * jrow] < 0.0 && sm.mesh.pressure_active[i + n * jrow] {
                    let (x, _) = ls.centre(i, jrow);
                    // Interface position by linear interpolation towards the next cell.
                    let (p0, p1) = (ls.psi[i + n * jrow], ls.psi[i + 1 + n * jrow]);
                    xf = if p1 >= 0.0 {
                        x + dx * p0 / (p0 - p1)
                    } else {
                        x
                    };
                }
            }
            front.push((tt, (xf - left) / a));
            next += 1;
        }
        for c in 0..n * n {
            if flow.liq.active_u[c] {
                max_speed = max_speed.max(flow.u[c].abs());
            }
        }
    }
    let surf = flow.surface.as_ref().expect("surface");
    DamBreak {
        steps,
        front,
        volume_drift: surf.ls.volume() / surf.volume0 - 1.0,
        max_speed,
    }
}

/// Outcome of the wall-impact reproduction.
#[derive(Clone, Copy, Debug)]
pub struct WallImpact {
    pub steps: usize,
    /// Time at which the speed exceeded 50 m/s, if it did.
    pub blew_up_at: Option<f64>,
    pub final_speed: f64,
    pub volume_drift: f64,
}

/// A liquid layer of `h_cells` cells thickness on the floor of a tank (`x in [-0.5, 0]` at first),
/// moving to the right with `u0` and hitting the wall at `x = wall_x`. `configure` tweaks the flow
/// (diagnostic switches) before the run.
pub fn verify_wall_impact(
    n: usize,
    h_cells: f64,
    u0: f64,
    nu: f64,
    t_end: f64,
    wall_x: f64,
    configure: &dyn Fn(&mut StaggeredFlow),
) -> WallImpact {
    let g = 9.81;
    let left = -0.5;
    let dx0 = 1.1 / n as f64;
    let floor = if std::env::var("ALIGN").is_ok() {
        -0.55 + (0.05 / dx0).round() * dx0
    } else {
        -0.5
    };
    let top = floor + 0.95;
    let sm = StaggeredMesh::new(n, 0.55, |x, y| {
        (x - wall_x).max(left - x).max(floor - y).max(y - top)
    });
    let dx = sm.dx();
    let h = h_cells * dx;
    let mut flow = StaggeredFlow::new(&sm, nu, |_, _| (0.0, 0.0));
    flow.upwind = true;
    flow.weno = true;
    let ls = FsLevelSet::new(n, 0.55, |x, y| (x - 0.0).max(-0.4 - x).max(y - (floor + h)));
    flow.enable_free_surface(ls, (0.0, -g));
    flow.surface.as_mut().expect("surface").ls.reinitialize(20);
    configure(&mut flow);
    flow.set_velocity(|x, y| {
        if x < 0.02 && x > -0.42 && y < floor + h + 0.02 {
            (u0, 0.0)
        } else {
            (0.0, 0.0)
        }
    });
    let dt = 0.2 * dx / (u0 + 2.0 * (g * h).sqrt()).max(1.0);
    let steps = (t_end / dt).ceil() as usize;
    let mut blew = None;
    let mut speed = 0.0f64;
    for s in 0..steps {
        flow.step(dt);
        speed = 0.0;
        for c in 0..n * n {
            if flow.liq.active_u[c] {
                speed = speed.max(flow.u[c].abs());
            }
            if flow.liq.active_v[c] {
                speed = speed.max(flow.v[c].abs());
            }
        }
        if speed.is_nan() || speed >= 50.0 {
            blew = Some((s + 1) as f64 * dt);
            break;
        }
    }
    let surf = flow.surface.as_ref().expect("surface");
    WallImpact {
        steps,
        blew_up_at: blew,
        final_speed: speed,
        volume_drift: surf.ls.volume() / surf.volume0 - 1.0,
    }
}

/// Rigid translation of a liquid block (side `side`) at `u0` in zero gravity inside a large
/// tank: the exact solution keeps the velocity uniform and the pressure zero. Returns the largest
/// deviation of the speed from `u0` over active nodes after `t_end`, or `None` on blow-up.
pub fn verify_translation(
    n: usize,
    side: f64,
    u0: f64,
    nu: f64,
    t_end: f64,
    configure: &dyn Fn(&mut StaggeredFlow),
) -> Option<(f64, f64)> {
    let sm = StaggeredMesh::new(n, 0.55, |x, y| x.abs().max(y.abs()) - 0.5);
    let mut flow = StaggeredFlow::new(&sm, nu, |_, _| (0.0, 0.0));
    flow.upwind = true;
    flow.weno = true;
    let ls = FsLevelSet::new(n, 0.55, |x, y| (x + 0.25).abs().max(y.abs()) - 0.5 * side);
    flow.enable_free_surface(ls, (0.0, 0.0));
    flow.surface.as_mut().expect("surface").ls.reinitialize(20);
    configure(&mut flow);
    flow.set_velocity(|x, y| {
        if (x + 0.25).abs().max(y.abs()) < 0.5 * side + 0.02 {
            (u0, 0.0)
        } else {
            (0.0, 0.0)
        }
    });
    let dt = 0.2 * sm.dx() / u0;
    let steps = (t_end / dt).ceil() as usize;
    for _ in 0..steps {
        flow.step(dt);
    }
    let mut dev = 0.0f64;
    let mut vmax = 0.0f64;
    for c in 0..n * n {
        if flow.liq.active_u[c] {
            dev = dev.max((flow.u[c] - u0).abs());
        }
        if flow.liq.active_v[c] {
            vmax = vmax.max(flow.v[c].abs());
        }
    }
    (dev.is_finite() && dev < 50.0).then_some((dev, vmax))
}

/// Like [`verify_wall_impact`] but reports the speed history only: time-averaged largest
/// `|v|` over active nodes in windows of `t_end / 10`, to see whether the transverse velocity
/// grows exponentially (instability) from the exact uniform-translation start.
pub fn track_sheet_growth(
    n: usize,
    h_cells: f64,
    u0: f64,
    t_end: f64,
    align: bool,
    configure: &dyn Fn(&mut StaggeredFlow),
) -> Vec<(f64, f64, f64, (usize, usize))> {
    let left = -0.5;
    let dx0 = 1.1 / n as f64;
    let floor = if align {
        -0.55 + (0.05 / dx0).round() * dx0
    } else {
        -0.5
    };
    let top = floor + 0.95;
    let sm = StaggeredMesh::new(n, 0.55, |x, y| {
        (x - 0.5).max(left - x).max(floor - y).max(y - top)
    });
    let dx = sm.dx();
    let h = h_cells * dx;
    let mut flow = StaggeredFlow::new(&sm, 1e-6, |_, _| (0.0, 0.0));
    flow.upwind = true;
    flow.weno = true;
    let ls = FsLevelSet::new(n, 0.55, |x, y| (x - 0.0).max(-0.4 - x).max(y - (floor + h)));
    flow.enable_free_surface(ls, (0.0, 0.0));
    flow.surface.as_mut().expect("surface").ls.reinitialize(20);
    configure(&mut flow);
    flow.set_velocity(|x, y| {
        if x < 0.02 && x > -0.42 && y < floor + h + 0.02 {
            (u0, 0.0)
        } else {
            (0.0, 0.0)
        }
    });
    let dt = 0.2 * dx / u0;
    let steps = (t_end / dt).ceil() as usize;
    let mut out = Vec::new();
    let every = (steps / 10).max(1);
    for s in 0..steps {
        flow.step(dt);
        if s % every == every - 1 {
            let (mut vmax, mut udev) = (0.0f64, 0.0f64);
            let mut at = (0, 0);
            for c in 0..n * n {
                if flow.liq.active_v[c] {
                    if flow.v[c].abs() > vmax {
                        at = (c % n, c / n);
                    }
                    vmax = vmax.max(flow.v[c].abs());
                }
                if flow.liq.active_u[c] {
                    udev = udev.max((flow.u[c] - u0).abs());
                }
            }
            out.push(((s + 1) as f64 * dt, vmax, udev, at));
            if vmax.is_nan() || udev.is_nan() || vmax >= 50.0 || udev >= 50.0 {
                break;
            }
        }
    }
    out
}

/// Result of the rimming-flow test.
#[derive(Clone, Debug)]
pub struct Rimming {
    pub steps: usize,
    /// Polar angle from the bottom (counter-clockwise), simulated and thin-film film thickness.
    pub phi: Vec<f64>,
    pub h_sim: Vec<f64>,
    pub h_theory: Vec<f64>,
    pub max_rel_err: f64,
    pub rms_rel_err: f64,
    pub volume_drift: f64,
    /// Thin-film flux `q` fitted to the liquid volume, over its maximum.
    pub q_over_qmax: f64,
    /// Largest film thickness over the drum radius.
    pub h_over_r: f64,
}

/// Thin-film (Moffatt 1977) steady thickness `h(phi)` of a viscous film on the inside of a
/// cylinder of radius `radius` rotating at `omega` (wall speed `u = omega radius`), for the flux
/// `q`; lubrication limit `h << R`, `u^2 / (g R) << 1`. `phi` is measured from the bottom, in the
/// direction of rotation.
fn moffatt_thickness(g: f64, nu: f64, u: f64, q: f64, phi: f64) -> f64 {
    let a = g * phi.sin() / (3.0 * nu);
    let f = |h: f64| u * h - a * h * h * h - q;
    let (mut lo, mut hi) = (
        0.0,
        if a > 1e-14 {
            (u / (3.0 * a)).sqrt()
        } else {
            q / u
        },
    );
    if f(hi) < 0.0 {
        return f64::NAN; // no film solution at this flux
    }
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// A viscous liquid film in a horizontally rotating circular drum (radius 0.5), started from a
/// uniform ring of thickness `h0` and integrated for `revs` revolutions; compares the final film
/// thickness with the thin-film solution of the same liquid volume.
pub fn verify_rimming(n: usize, nu: f64, omega: f64, h0: f64, revs: f64) -> Rimming {
    let (radius, g) = (0.5, 9.81);
    let sm = StaggeredMesh::new(n, 0.55, |x, y| (x * x + y * y).sqrt() - radius);
    let mut flow = StaggeredFlow::new(&sm, nu, move |x, y| (-omega * y, omega * x));
    flow.upwind = true;
    flow.weno = true;
    let ls = FsLevelSet::new(n, 0.55, |x, y| (radius - h0) - (x * x + y * y).sqrt());
    flow.enable_free_surface(ls, (0.0, -g));
    flow.set_velocity(move |x, y| {
        if (x * x + y * y).sqrt() > radius - h0 - 0.02 {
            (-omega * y, omega * x)
        } else {
            (0.0, 0.0)
        }
    });
    let dx = sm.dx();
    let u_wall = omega * radius;
    let dt = 0.2 * dx / u_wall;
    let t_end = revs * 2.0 * PI / omega;
    let steps = (t_end / dt).ceil() as usize;
    let dt = t_end / steps as f64;
    for _ in 0..steps {
        flow.step(dt);
    }
    let surf = flow.surface.as_ref().expect("free surface enabled");
    let ls = &surf.ls;
    // Bilinear psi at a point (cell centres).
    let psi_at = |x: f64, y: f64| -> f64 {
        let (fx, fy) = ((x + 0.55) / dx - 0.5, (y + 0.55) / dx - 0.5);
        let (i, j) = (fx.floor() as usize, fy.floor() as usize);
        let (sx, sy) = (fx - i as f64, fy - j as f64);
        let at = |a: usize, b: usize| ls.psi[a + n * b];
        (1.0 - sx) * (1.0 - sy) * at(i, j)
            + sx * (1.0 - sy) * at(i + 1, j)
            + (1.0 - sx) * sy * at(i, j + 1)
            + sx * sy * at(i + 1, j + 1)
    };
    let samples = 72;
    let mut phi = Vec::new();
    let mut h_sim = Vec::new();
    for k in 0..samples {
        let p = 2.0 * PI * (k as f64 + 0.5) / samples as f64;
        let dir = (p.sin(), -p.cos());
        // March inward from the wall until psi turns positive (air).
        let mut r = radius - 0.25 * dx;
        let mut h = f64::NAN;
        let mut prev = psi_at(dir.0 * r, dir.1 * r);
        while r > radius - 0.3 {
            let rn = r - 0.05 * dx;
            let cur = psi_at(dir.0 * rn, dir.1 * rn);
            if prev < 0.0 && cur >= 0.0 {
                let t = prev / (prev - cur);
                h = radius - (r - t * 0.05 * dx);
                break;
            }
            prev = cur;
            r = rn;
        }
        phi.push(p);
        h_sim.push(h);
    }
    // Fit the flux q to the liquid volume (annulus area = integral of (R h - h^2 / 2) dphi).
    let volume = surf.volume0;
    let u = u_wall;
    let h_c = (nu * u / g).sqrt();
    let q_max = 2.0 / 3.0 * u * h_c;
    let area = |q: f64| -> f64 {
        let m = 720;
        (0..m)
            .map(|k| {
                let p = 2.0 * PI * (k as f64 + 0.5) / m as f64;
                let h = moffatt_thickness(g, nu, u, q, p);
                (radius * h - 0.5 * h * h) * 2.0 * PI / m as f64
            })
            .sum::<f64>()
    };
    let (mut lo, mut hi) = (0.0, q_max * (1.0 - 1e-9));
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if area(mid) < volume {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let q = 0.5 * (lo + hi);
    let h_theory: Vec<f64> = phi
        .iter()
        .map(|&p| moffatt_thickness(g, nu, u, q, p))
        .collect();
    let rel: Vec<f64> = h_sim
        .iter()
        .zip(&h_theory)
        .map(|(a, b)| (a - b).abs() / b)
        .collect();
    Rimming {
        steps,
        max_rel_err: rel.iter().cloned().fold(0.0, f64::max),
        rms_rel_err: (rel.iter().map(|e| e * e).sum::<f64>() / rel.len() as f64).sqrt(),
        volume_drift: surf.ls.volume() / surf.volume0 - 1.0,
        q_over_qmax: q / q_max,
        h_over_r: h_theory.iter().cloned().fold(0.0, f64::max) / radius,
        phi,
        h_sim,
        h_theory,
    }
}

/// Result of the rigid-rotation liquid ring test.
#[derive(Clone, Copy, Debug)]
pub struct RigidRing {
    pub steps: usize,
    /// Largest deviation of the free surface from the exact circle, in cells.
    pub surface_err_cells: f64,
    /// Largest velocity deviation from the rigid rotation over `omega R`.
    pub velocity_err: f64,
    /// Where it occurs: polar radius and angle (degrees).
    pub velocity_err_at: (f64, f64),
    /// Largest `|p - p_exact| / (rho omega^2 R^2 / 2)` over liquid cells.
    pub pressure_err: f64,
    pub volume_drift: f64,
}

/// Exact steady solution with a free surface: a liquid ring in rigid rotation at `omega` inside a
/// drum of radius 0.5 rotating at `omega`, gravity downward. The pressure `p = omega^2 r^2 / 2 -
/// g y + C` is constant on a circle of radius `r_i` centred at `(0, g / omega^2)`, so that circle
/// is an exact free surface for any viscosity. The run starts from this state and checks that it
/// is preserved.
pub fn verify_rigid_ring(
    n: usize,
    omega: f64,
    r_i: f64,
    nu: f64,
    t_end: f64,
    g: f64,
    configure: &dyn Fn(&mut StaggeredFlow),
) -> RigidRing {
    let radius = 0.5;
    let e = g / (omega * omega);
    let sm = StaggeredMesh::new(n, 0.55, |x, y| (x * x + y * y).sqrt() - radius);
    let mut flow = StaggeredFlow::new(&sm, nu, move |x, y| (-omega * y, omega * x));
    flow.upwind = true;
    flow.weno = true;
    configure(&mut flow);
    let ls = FsLevelSet::new(n, 0.55, move |x, y| {
        r_i - (x * x + (y - e) * (y - e)).sqrt()
    });
    flow.enable_free_surface(ls, (0.0, -g));
    flow.surface.as_mut().expect("surface").ls.reinitialize(20);
    flow.set_velocity(move |x, y| (-omega * y, omega * x));
    let dx = sm.dx();
    let dt = 0.2 * dx / (omega * radius);
    let steps = (t_end / dt).ceil() as usize;
    let dt = t_end / steps as f64;
    for _ in 0..steps {
        flow.step(dt);
    }
    let surf = flow.surface.as_ref().expect("surface");
    let ls = &surf.ls;
    // Surface error: interface radius from the circle centre along rays.
    let mut surface_err = 0.0f64;
    let psi_at = |x: f64, y: f64| -> f64 {
        let (fx, fy) = ((x + 0.55) / dx - 0.5, (y + 0.55) / dx - 0.5);
        let (i, j) = (fx.floor() as usize, fy.floor() as usize);
        let (sx, sy) = (fx - i as f64, fy - j as f64);
        let at = |a: usize, b: usize| ls.psi[a + n * b];
        (1.0 - sx) * (1.0 - sy) * at(i, j)
            + sx * (1.0 - sy) * at(i + 1, j)
            + (1.0 - sx) * sy * at(i, j + 1)
            + sx * sy * at(i + 1, j + 1)
    };
    for k in 0..90 {
        let a = 2.0 * PI * (k as f64 + 0.5) / 90.0;
        let dir = (a.cos(), a.sin());
        // March outward from inside the circle (air, psi > 0) until psi turns negative.
        let (mut r, step) = (r_i - 0.2, 0.05 * dx);
        let mut prev = psi_at(dir.0 * r, e + dir.1 * r);
        while r < r_i + 0.2 {
            let rn = r + step;
            let (x, y) = (dir.0 * rn, e + dir.1 * rn);
            if x.abs() > 0.5 || (x * x + y * y).sqrt() > radius {
                break;
            }
            let cur = psi_at(x, y);
            if prev >= 0.0 && cur < 0.0 {
                let t = prev / (prev - cur);
                surface_err = surface_err.max(((r + t * step) - r_i).abs() / dx);
                break;
            }
            prev = cur;
            r = rn;
        }
    }
    let mut verr = 0.0f64;
    let mut verr_at = (0.0, 0.0);
    let mut perr = 0.0f64;
    let top = e + r_i;
    for j in 0..n {
        for i in 0..n {
            let c = i + n * j;
            if flow.liq.active_u[c] {
                let (x, y) = sm.gu.centre(i, j);
                let d = (flow.u[c] + omega * y).abs();
                if d > verr {
                    verr = d;
                    verr_at = ((x * x + y * y).sqrt(), y.atan2(x).to_degrees());
                }
            }
            if flow.liq.active_v[c] {
                let (x, y) = sm.gv.centre(i, j);
                let d = (flow.v[c] - omega * x).abs();
                if d > verr {
                    verr = d;
                    verr_at = ((x * x + y * y).sqrt(), y.atan2(x).to_degrees());
                }
            }
            if flow.liq.cell[c] {
                let (x, y) = ls.centre(i, j);
                let exact = 0.5 * omega * omega * (x * x + y * y - top * top) - g * (y - top);
                perr = perr.max((flow.p[c] - exact).abs());
            }
        }
    }
    RigidRing {
        steps,
        surface_err_cells: surface_err,
        velocity_err: verr / (omega * radius),
        velocity_err_at: verr_at,
        pressure_err: perr / (0.5 * omega * omega * radius * radius),
        volume_drift: surf.ls.volume() / surf.volume0 - 1.0,
    }
}

/// Result of the viscous-slurry drum run.
#[derive(Clone, Debug)]
pub struct DrumSlurry {
    pub steps: usize,
    /// Mean wall torque on the liquid per unit density over the measured revolutions.
    pub torque: f64,
    /// Mean gravity balance `A g x_cm` over the same window (equals the wall torque in steady
    /// state because the pressure exerts no torque on a circular wall).
    pub torque_gravity: f64,
    /// RMS deviation of the instantaneous wall torque from its mean, relative to the mean.
    pub torque_ripple: f64,
    /// Mean wall torque of each measured revolution.
    pub per_rev: Vec<f64>,
    pub volume_drift: f64,
    /// Film thickness at 72 angles from the bottom (direction of rotation), at the end.
    pub film: Vec<f64>,
    /// Mean gravity torque from the liquid faces (`sum -g x dV` over liquid v-faces) as a check
    /// of the level-set centroid.
    pub torque_gravity_faces: f64,
    /// Mean torque of the pressure on the cut wall cells (zero for an exact circle).
    pub torque_pressure: f64,
    /// Rate of change of the liquid angular momentum over the measured window.
    pub l_rate: f64,
}

/// Viscous liquid in a horizontally rotating drum of radius 0.5 filled to the area fraction
/// `fill`, started from rest with a flat surface. `omega` is the wall rotation rate; the run
/// settles for `settle` revolutions and then averages the wall torque over `measure` revolutions.
pub fn verify_drum_slurry(
    n: usize,
    nu: f64,
    omega: f64,
    fill: f64,
    settle: f64,
    measure: f64,
) -> DrumSlurry {
    let (radius, g) = (0.5, 9.81);
    let sm = StaggeredMesh::new(n, 0.55, |x, y| (x * x + y * y).sqrt() - radius);
    let mut flow = StaggeredFlow::new(&sm, nu, move |x, y| (-omega * y, omega * x));
    flow.upwind = true;
    flow.weno = true;
    flow.stokes = std::env::var("E3_STOKES").is_ok();
    flow.robust_surface = std::env::var("E3_ROBUST").is_ok();
    // Uniform ring of the requested area fraction, started in rigid rotation (a flat pool meeting
    // the moving wall would create a contact-line singularity that stretches a sub-grid sheet).
    let r_i = radius * (1.0 - fill).sqrt();
    let ls = FsLevelSet::new(n, 0.55, move |x, y| r_i - (x * x + y * y).sqrt());
    flow.enable_free_surface(ls, (0.0, -g));
    flow.set_velocity(move |x, y| {
        if (x * x + y * y).sqrt() > r_i - 0.02 {
            (-omega * y, omega * x)
        } else {
            (0.0, 0.0)
        }
    });
    let dx = sm.dx();
    let dt = 0.2 * dx / (omega * radius);
    let rev = 2.0 * PI / omega;
    let per_rev_steps = (rev / dt).ceil() as usize;
    let dt = rev / per_rev_steps as f64;
    let settle_steps = (settle * per_rev_steps as f64).round() as usize;
    for k in 0..settle_steps {
        flow.step(dt);
        if std::env::var("E3_DUMP").is_ok() {
            let mut best = (0.0f64, 0usize, 0usize);
            for j in 0..n {
                for i in 0..n {
                    let c = i + n * j;
                    if flow.liq.active_u[c] && flow.u[c].abs() > best.0 {
                        best = (flow.u[c].abs(), i, j);
                    }
                }
            }
            if best.0 > 1.45 {
                let (bi, bj) = (best.1, best.2);
                println!("EXCEED at step {k}: |u|={:.4} node ({bi},{bj})", best.0);
                let ls = &flow.surface.as_ref().unwrap().ls;
                for j in (bj.saturating_sub(3)..(bj + 4).min(n)).rev() {
                    let mut row = String::new();
                    for i in bi.saturating_sub(3)..(bi + 4).min(n) {
                        let c = i + n * j;
                        row += &format!(
                            "{}{:>7.3}({:>5.2}) ",
                            if flow.liq.active_u[c] { '*' } else { ' ' },
                            flow.u[c],
                            ls.psi[c] / sm.dx()
                        );
                    }
                    println!("{row}");
                }
                break;
            }
        }
        if std::env::var("E3_TRACE").is_ok() && k % 25 == 0 {
            let mut best = (0.0f64, 0.0f64, 0.0f64);
            for j in 0..n {
                for i in 0..n {
                    let c = i + n * j;
                    if flow.liq.active_u[c] && flow.u[c].abs() > best.0 {
                        let (x, y) = sm.gu.centre(i, j);
                        best = (flow.u[c].abs(), x, y);
                    }
                    if flow.liq.active_v[c] && flow.v[c].abs() > best.0 {
                        let (x, y) = sm.gv.centre(i, j);
                        best = (flow.v[c].abs(), x, y);
                    }
                }
            }
            let vol = flow.surface.as_ref().map(|s| s.ls.volume()).unwrap_or(0.0);
            println!(
                "step {k} liquid max {:.4} at ({:.3},{:.3}) vol {vol:.6}",
                best.0, best.1, best.2
            );
        }
    }
    let ang_mom = |flow: &StaggeredFlow| -> f64 {
        let mut l = 0.0;
        for jj in 0..n {
            for ii in 0..n {
                let c = ii + n * jj;
                if flow.liq.active_u[c] {
                    let (_, y) = sm.gu.centre(ii, jj);
                    l -= y * flow.u[c] * dx * dx * sm.au[c];
                }
                if flow.liq.active_v[c] {
                    let (x, _) = sm.gv.centre(ii, jj);
                    l += x * flow.v[c] * dx * dx * sm.av[c];
                }
            }
        }
        l
    };
    let l_start = ang_mom(&flow);
    let measure_revs = measure.round().max(1.0) as usize;
    let mut series = Vec::new();
    let mut gravity = Vec::new();
    let mut gravity_faces = Vec::new();
    let mut pressure_torque = Vec::new();
    let mut per_rev = Vec::new();
    for _ in 0..measure_revs {
        let mut sum = 0.0;
        for _ in 0..per_rev_steps {
            flow.step(dt);
            let t = flow.rotating_wall_torque(omega);
            sum += t;
            series.push(t);
            let surf = flow.surface.as_ref().expect("free surface enabled");
            gravity.push(surf.volume0 * g * surf.ls.centroid_x());
            let mut gf = 0.0;
            for jj in 0..n {
                for ii in 0..n {
                    let c = ii + n * jj;
                    if flow.liq.active_v[c] {
                        let (x, _) = sm.gv.centre(ii, jj);
                        gf += g * x * dx * dx * sm.av[c];
                    }
                }
            }
            gravity_faces.push(gf);
            let mut tp = 0.0;
            for jj in 0..n {
                for ii in 0..n {
                    let c = ii + n * jj;
                    if !flow.liq.cell[c] {
                        continue;
                    }
                    let fx = (sm.au[c + 1] - sm.au[c]) * dx;
                    let fy = (sm.av[c + n] - sm.av[c]) * dx;
                    let (x, y) = (
                        -0.55 + (ii as f64 + 0.5) * dx,
                        -0.55 + (jj as f64 + 0.5) * dx,
                    );
                    tp += flow.p[c] * (x * fy - y * fx);
                }
            }
            pressure_torque.push(tp);
        }
        per_rev.push(sum / per_rev_steps as f64);
    }
    let count = series.len() as f64;
    let mean = series.iter().sum::<f64>() / count;
    let ripple = (series.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / count).sqrt()
        / mean.abs().max(1e-30);
    let surf = flow.surface.as_ref().expect("free surface enabled");
    let ls = &surf.ls;
    let psi_at = |x: f64, y: f64| -> f64 {
        let (fx, fy) = ((x + 0.55) / dx - 0.5, (y + 0.55) / dx - 0.5);
        let (i, j) = (fx.floor() as usize, fy.floor() as usize);
        let (sx, sy) = (fx - i as f64, fy - j as f64);
        let at = |a: usize, b: usize| ls.psi[a + n * b];
        (1.0 - sx) * (1.0 - sy) * at(i, j)
            + sx * (1.0 - sy) * at(i + 1, j)
            + (1.0 - sx) * sy * at(i, j + 1)
            + sx * sy * at(i + 1, j + 1)
    };
    let film: Vec<f64> = (0..72)
        .map(|k| {
            let p = 2.0 * PI * (k as f64 + 0.5) / 72.0;
            let dir = (p.sin(), -p.cos());
            let mut r = radius - 0.25 * dx;
            let mut prev = psi_at(dir.0 * r, dir.1 * r);
            while r > 0.05 {
                let rn = r - 0.05 * dx;
                let cur = psi_at(dir.0 * rn, dir.1 * rn);
                if prev < 0.0 && cur >= 0.0 {
                    return radius - (r - prev / (prev - cur) * 0.05 * dx);
                }
                prev = cur;
                r = rn;
            }
            f64::NAN
        })
        .collect();
    DrumSlurry {
        film,
        torque_gravity_faces: gravity_faces.iter().sum::<f64>() / count,
        torque_pressure: pressure_torque.iter().sum::<f64>() / count,
        l_rate: (ang_mom(&flow) - l_start) / (measure_revs as f64 * rev),
        steps: settle_steps + series.len(),
        torque: mean,
        torque_gravity: gravity.iter().sum::<f64>() / count,
        torque_ripple: ripple,
        per_rev,
        volume_drift: surf.ls.volume() / surf.volume0 - 1.0,
    }
}

/// Result of the fixed-disc buoyancy test.
#[derive(Clone, Copy, Debug)]
pub struct Buoyancy {
    pub steps: usize,
    /// Vertical pressure force on the disc over the displaced weight `g pi r^2`, minus one, with
    /// the cell-centre pressure used as is.
    pub force_err_plain: f64,
    /// Same with the pressure moved to the wall point along the local pressure gradient.
    pub force_err: f64,
    /// Horizontal force over the displaced weight.
    pub side_force: f64,
    /// Torque about the disc centre over `g r^3`.
    pub torque: f64,
    pub spurious: f64,
}

/// Fixed disc (radius `r`, centre `(0.05, -0.12)`) submerged in a still pool of depth up to
/// `y = 0.2` inside the drum of radius 0.5: the pressure force on the disc must equal the weight
/// of the displaced liquid (per unit density, `g pi r^2`) and its torque must vanish.
pub fn verify_buoyancy(n: usize, r: f64, t_end: f64) -> Buoyancy {
    let (radius, g, level) = (0.5, 9.81, 0.2);
    let (cx, cy) = (0.05, -0.12);
    let sm = StaggeredMesh::new(n, 0.55, move |x, y| {
        let drum = (x * x + y * y).sqrt() - radius;
        let ball = r - ((x - cx).powi(2) + (y - cy).powi(2)).sqrt();
        drum.max(ball)
    });
    let mut flow = StaggeredFlow::new(&sm, 1e-3, |_, _| (0.0, 0.0));
    flow.upwind = true;
    flow.weno = true;
    let ls = FsLevelSet::new(n, 0.55, move |_, y| y - level);
    flow.enable_free_surface(ls, (0.0, -g));
    let dx = sm.dx();
    let dt0 = 0.3 * dx / (g * 2.0 * radius).sqrt();
    let steps = (t_end / dt0).ceil() as usize;
    let dt = t_end / steps as f64;
    for _ in 0..steps {
        flow.step(dt);
    }
    let mut vmax = 0.0f64;
    for c in 0..n * n {
        if flow.liq.active_u[c] {
            vmax = vmax.max(flow.u[c].abs());
        }
        if flow.liq.active_v[c] {
            vmax = vmax.max(flow.v[c].abs());
        }
    }
    let centre =
        |i: usize, j: usize| (-0.55 + (i as f64 + 0.5) * dx, -0.55 + (j as f64 + 0.5) * dx);
    let (mut fx, mut fy, mut tq) = (0.0, 0.0, 0.0);
    let mut fy0 = 0.0;
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            let c = i + n * j;
            if !flow.liq.cell[c] {
                continue;
            }
            let (x, y) = centre(i, j);
            let d = ((x - cx).powi(2) + (y - cy).powi(2)).sqrt();
            if (d - r).abs() > 1.5 * dx {
                continue;
            }
            // Force on the fluid from the solid is p (a_r - a_l, a_t - a_b) dx; the disc feels the
            // opposite.
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
            let (sx, sy) = (cx + r * (x - cx) / d, cy + r * (y - cy) / d);
            let pw = flow.p[c] + gx * (sx - x) + gy * (sy - y);
            fy0 -= flow.p[c] * ay;
            let (px, py) = (-pw * ax, -pw * ay);
            fx += px;
            fy += py;
            tq += (sx - cx) * py - (sy - cy) * px;
        }
    }
    let weight = g * PI * r * r;
    Buoyancy {
        steps,
        force_err_plain: fy0 / weight - 1.0,
        force_err: fy / weight - 1.0,
        side_force: fx / weight,
        torque: tq / (g * r * r * r),
        spurious: vmax / (g * 2.0 * radius).sqrt(),
    }
}

/// Steady Stokes drag (force per unit density, i.e. `4 pi nu |D|`) on a disc of radius `a`
/// translating at `u` inside a fixed concentric cylinder of radius `b`: the stream function
/// `psi = sin(theta) (A r + B / r + C r^3 + D r ln r)` with `f(a) = u a`, `f'(a) = u`,
/// `f(b) = f'(b) = 0`; the net force is carried by the Stokeslet coefficient `D`.
pub fn stokes_annulus_drag(a: f64, b: f64, u: f64, nu: f64) -> f64 {
    let rows = |r: f64| {
        (
            [r, 1.0 / r, r * r * r, r * r.ln()],
            [1.0, -1.0 / (r * r), 3.0 * r * r, r.ln() + 1.0],
        )
    };
    let mut m = [[0.0; 5]; 4];
    m[0][..4].copy_from_slice(&rows(a).0);
    m[0][4] = u * a;
    m[1][..4].copy_from_slice(&rows(a).1);
    m[1][4] = u;
    m[2][..4].copy_from_slice(&rows(b).0);
    m[3][..4].copy_from_slice(&rows(b).1);
    for k in 0..4 {
        let piv = (k..4)
            .max_by(|&i, &j| m[i][k].abs().total_cmp(&m[j][k].abs()))
            .unwrap_or(k);
        m.swap(k, piv);
        for i in k + 1..4 {
            let f = m[i][k] / m[k][k];
            let pivot_row = m[k];
            for (mij, pkj) in m[i].iter_mut().zip(pivot_row).skip(k) {
                *mij -= f * pkj;
            }
        }
    }
    let mut x = [0.0; 4];
    for k in (0..4).rev() {
        let s: f64 = (k + 1..4).map(|j| m[k][j] * x[j]).sum();
        x[k] = (m[k][4] - s) / m[k][k];
    }
    4.0 * PI * nu * x[3].abs()
}

/// Steady Stokes force (per unit density, `nu` the kinematic viscosity) on a disc of radius `a`
/// translating at `u` along the line of centres inside a fixed cylinder of radius `b`, with the
/// centres a distance `e` apart. Exact bipolar-coordinate solution: with `x + i y = i k cot(...)`
/// the stream function is `psi = f(xi) sin(eta) / (cosh xi - cos eta)`, where only the first
/// harmonic survives, `f = A (sinh 2 s - 2 s) + B (cosh 2 s - 1)` in `s = xi - xi_outer`; the
/// force is the Stokeslet strength `4 pi nu |d| / k` of the term `d xi` in `f`, which gives
/// `F = 8 pi nu u sinh(2 D) / (C^2 - S sinh(2 D))` with `D = xi_inner - xi_outer`,
/// `C = cosh 2 D - 1`, `S = sinh 2 D - 2 D`. Depends on the bipolar distance `D` only.
pub fn eccentric_squeeze_force(a: f64, b: f64, e: f64, u: f64, nu: f64) -> f64 {
    // k from e = sqrt(k^2 + b^2) - sqrt(k^2 + a^2) (decreasing in k).
    let (mut lo, mut hi) = (0.0f64, 1e9f64);
    for _ in 0..200 {
        let k = 0.5 * (lo + hi);
        let ek = (k * k + b * b).sqrt() - (k * k + a * a).sqrt();
        if ek > e {
            lo = k;
        } else {
            hi = k;
        }
    }
    let k = 0.5 * (lo + hi);
    let delta = (k / a).asinh() - (k / b).asinh();
    let t = 2.0 * delta;
    let s2 = t.sinh();
    let c = 2.0 * delta.sinh().powi(2);
    // sinh t - t, by its series for small t (cancellation).
    let s = if t < 0.2 {
        let t2 = t * t;
        t * t2 / 6.0
            * (1.0 + t2 / 20.0 * (1.0 + t2 / 42.0 * (1.0 + t2 / 72.0 * (1.0 + t2 / 110.0))))
    } else {
        s2 - t
    };
    8.0 * PI * nu * u * s2 / (c * c - s * s2).abs()
}

/// Steps over which the force of a translating-disc run is averaged.
const WINDOW: usize = 20;

/// Result of a translating-disc run.
#[derive(Clone, Copy, Debug)]
pub struct MovingDisc {
    pub steps: usize,
    /// Force on the disc along its direction of motion (per unit density).
    pub drag: f64,
    /// Mean and relative RMS deviation of the drag over the last `WINDOW` steps.
    pub drag_mean: f64,
    pub drag_noise: f64,
    pub side: f64,
    pub torque: f64,
    pub divergence: f64,
    pub centre_x: f64,
}

/// Disc of radius `a` translating at `u` (+x) inside the drum of radius 0.5, viscosity `nu`, fluid
/// started at rest, integrated to `t_end` with step `dt`. With `moving` the disc advances (the mesh
/// is rebuilt every step), else its surface velocity is imposed at a fixed position `x0`.
#[allow(clippy::too_many_arguments)]
pub fn verify_moving_disc(
    n: usize,
    a: f64,
    x0: f64,
    u: f64,
    omega: f64,
    nu: f64,
    t_end: f64,
    dt: f64,
    moving: bool,
) -> MovingDisc {
    let radius = 0.5;
    let steps = (t_end / dt).round() as usize;
    let mut disc = Disc {
        cx: x0,
        cy: 0.0,
        r: a,
        ux: u,
        uy: 0.0,
        omega,
    };
    let build = |d: &Disc| {
        let dd = *d;
        let mut sm = StaggeredMesh::new(n, 0.55, move |x, y| {
            ((x * x + y * y).sqrt() - radius).max(-dd.sdf(x, y))
        });
        sm.set_body_flux(&[*d]);
        sm
    };
    let mut state: Option<FlowState> = None;
    let mut prev_sm: Option<StaggeredMesh> = None;
    let mut last = (0.0, 0.0, 0.0, 0.0);
    let mut window: Vec<f64> = Vec::new();
    for step in 0..steps {
        if moving {
            disc.cx = x0 + u * dt * (step + 1) as f64;
        }
        let sm = build(&disc);
        let d = disc;
        let dx = sm.dx();
        let wall = move |x: f64, y: f64| {
            if d.sdf(x, y).abs() < 0.6 * dx {
                d.velocity(x, y)
            } else {
                (0.0, 0.0)
            }
        };
        let mut flow = match state.take() {
            Some(s) => StaggeredFlow::from_state(&sm, s, wall),
            None => StaggeredFlow::new(&sm, nu, wall),
        };
        if let Some(old) = &prev_sm {
            flow.initialise_fresh(old, move |x, y| d.velocity(x, y));
        }
        flow.step(dt);
        if step + WINDOW >= steps {
            let (fx, _, _) = disc_load(&flow, &disc);
            window.push(-fx);
        }
        if step + 1 == steps {
            let (fx, fy, tq) = disc_load(&flow, &disc);
            last = (fx, fy, tq, flow.projected_divergence);
        }
        state = Some(flow.into_state());
        prev_sm = Some(sm);
    }
    let mean = window.iter().sum::<f64>() / window.len() as f64;
    let rms = (window.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / window.len() as f64).sqrt();
    MovingDisc {
        steps,
        drag: -last.0,
        drag_mean: mean,
        drag_noise: rms / mean,
        side: last.1,
        torque: last.2,
        divergence: last.3,
        centre_x: disc.cx,
    }
}

/// Viscous torque per unit spin rate of a disc of radius `a` in a concentric drum of radius `b`
/// (steady Couette flow): `4 pi nu a^2 b^2 / (b^2 - a^2)`.
pub fn confined_spin_stiffness(a: f64, b: f64, nu: f64) -> f64 {
    4.0 * PI * nu * a * a * b * b / (b * b - a * a)
}

/// Added mass (per unit density) of a disc of radius `a` moving along a diameter of a concentric
/// cylinder of radius `b` in inviscid fluid: `pi a^2 (b^2 + a^2) / (b^2 - a^2)`.
pub fn confined_added_mass(a: f64, b: f64) -> f64 {
    PI * a * a * (b * b + a * a) / (b * b - a * a)
}

/// Result of the added-mass tests.
#[derive(Clone, Copy, Debug)]
pub struct AddedMass {
    pub steps: usize,
    /// Measured effective added mass over the confined inviscid value, minus one.
    pub ratio_err: f64,
    /// Mean fixed-point iterations per step (coupled run), else 0.
    pub iterations: f64,
    /// Disc velocity at the end over `F t / (m + m_a)`, minus one (coupled run), else 0.
    pub speed_err: f64,
    pub divergence: f64,
}

/// Prescribed acceleration: the disc (radius `a`, concentric) is accelerated from rest at `a0`
/// along +x in a drum of radius 0.5 filled with fluid of viscosity `nu`; the fluid reaction at
/// `t_end`, `-F = m_a a0` for an inviscid fluid, gives the effective added mass.
pub fn verify_added_mass_prescribed(
    n: usize,
    a: f64,
    a0: f64,
    nu: f64,
    t_end: f64,
    dt: f64,
) -> AddedMass {
    let radius = 0.5;
    let steps = (t_end / dt).round() as usize;
    let mut bf = BodyFlow::new(n, 0.55, radius, nu);
    let mut load = BodyLoad {
        fx: 0.0,
        fy: 0.0,
        torque: 0.0,
        divergence: 0.0,
    };
    for step in 0..steps {
        let t = dt * (step + 1) as f64;
        let disc = Disc {
            cx: 0.5 * a0 * t * t,
            cy: 0.0,
            r: a,
            ux: a0 * t,
            uy: 0.0,
            omega: 0.0,
        };
        let (l, state, mesh) = bf.trial(&disc, dt);
        bf.commit(state, mesh);
        load = l;
    }
    let ma = confined_added_mass(a, radius);
    AddedMass {
        steps,
        ratio_err: -load.fx / (ma * a0) - 1.0,
        iterations: 0.0,
        speed_err: 0.0,
        divergence: load.divergence,
    }
}

/// Free disc of mass `ratio * pi a^2` driven by a constant external force `force` along +x in
/// the same drum (coupled fluid-disc steps); compares the velocity at `t_end` with the inviscid
/// short-time value `F t / (m + m_a)`.
pub fn verify_added_mass_coupled(
    n: usize,
    a: f64,
    ratio: f64,
    force: f64,
    nu: f64,
    t_end: f64,
    dt: f64,
) -> AddedMass {
    let radius = 0.5;
    let steps = (t_end / dt).round() as usize;
    let mass = ratio * PI * a * a;
    let ma = confined_added_mass(a, radius);
    let mut bf = BodyFlow::new(n, 0.55, radius, nu);
    let mut body = Body {
        disc: Disc {
            cx: 0.0,
            cy: 0.0,
            r: a,
            ux: 0.0,
            uy: 0.0,
            omega: 0.0,
        },
        mass,
        inertia: 0.5 * mass * a * a,
        accel: (0.0, 0.0),
    };
    let (mut iters, mut div) = (0.0, 0.0f64);
    for _ in 0..steps {
        let info = bf.step_coupled(
            &mut body,
            (force, 0.0, 0.0),
            dt,
            ma,
            0.0,
            confined_spin_stiffness(a, radius, nu),
            1e-6,
            30,
        );
        iters += info.iterations as f64;
        div = div.max(info.load.divergence);
    }
    let expected = force * t_end / (mass + ma);
    AddedMass {
        steps,
        ratio_err: 0.0,
        iterations: iters / steps as f64,
        speed_err: body.disc.ux / expected - 1.0,
        divergence: div,
    }
}

/// Free disc pushed by a constant force in a viscous (Stokes) fluid: after the transient
/// (rate `drag / (m + m_a)`) it moves at the terminal speed `F / k` with `k` the steady
/// concentric Stokes drag coefficient. Returns the relative error of the speed at `t_end` and the
/// mean number of coupling iterations.
pub fn verify_terminal_stokes(
    n: usize,
    a: f64,
    ratio: f64,
    force: f64,
    nu: f64,
    t_end: f64,
    dt: f64,
) -> (f64, f64) {
    let radius = 0.5;
    let steps = (t_end / dt).round() as usize;
    let mass = ratio * PI * a * a;
    let ma = confined_added_mass(a, radius);
    let mut bf = BodyFlow::new(n, 0.55, radius, nu);
    let mut body = Body {
        disc: Disc {
            cx: 0.0,
            cy: 0.0,
            r: a,
            ux: 0.0,
            uy: 0.0,
            omega: 0.0,
        },
        mass,
        inertia: 0.5 * mass * a * a,
        accel: (0.0, 0.0),
    };
    let mut iters = 0.0;
    for _ in 0..steps {
        let info = bf.step_coupled(
            &mut body,
            (force, 0.0, 0.0),
            dt,
            ma,
            stokes_annulus_drag(a, radius, 1.0, nu),
            confined_spin_stiffness(a, radius, nu),
            1e-6,
            40,
        );
        iters += info.iterations as f64;
    }
    let k = stokes_annulus_drag(a, radius, 1.0, nu);
    (body.disc.ux / (force / k) - 1.0, iters / steps as f64)
}

/// Result of the falling-disc run.
#[derive(Clone, Debug)]
pub struct FallingDisc {
    pub steps: usize,
    /// Vertical position of the disc at the sampled times.
    pub y: Vec<f64>,
    pub vy_max: f64,
    pub iterations: f64,
    pub divergence: f64,
    /// Largest `|u|` of the fluid at the end (blow-up check).
    pub umax: f64,
    pub ok: bool,
}

/// Disc of radius `a` and density ratio `ratio` released at rest at `(0, y0)` in the drum of
/// radius 0.5 under gravity 9.81 (reduced gravity on the body, buoyancy through the fluid
/// pressure is replaced by the net weight `(ratio - 1) pi a^2 g`), viscosity `nu`; `dt` is
/// `cfl dx / sqrt(g a)`. Samples `y` every `t_end / 5`.
pub fn verify_falling_disc(
    n: usize,
    a: f64,
    ratio: f64,
    y0: f64,
    nu: f64,
    t_end: f64,
    cfl: f64,
) -> FallingDisc {
    let (radius, g) = (0.5, 9.81);
    let dx = 1.1 / n as f64;
    let dt0 = cfl * dx / (g * a).sqrt();
    let steps = (t_end / dt0).ceil() as usize;
    let dt = t_end / steps as f64;
    let mass = ratio * PI * a * a;
    let ma = confined_added_mass(a, radius);
    let mut bf = BodyFlow::new(n, 0.55, radius, nu);
    let mut body = Body {
        disc: Disc {
            cx: 0.0,
            cy: y0,
            r: a,
            ux: 0.0,
            uy: 0.0,
            omega: 0.0,
        },
        mass,
        inertia: 0.5 * mass * a * a,
        accel: (0.0, 0.0),
    };
    let weight = (ratio - 1.0) * PI * a * a * g;
    let (mut iters, mut div, mut vmax) = (0.0, 0.0f64, 0.0f64);
    let mut y = Vec::new();
    let mut ok = true;
    for step in 0..steps {
        let drag_est = stokes_annulus_drag(a, radius, 1.0, nu).min(1e3);
        let info = bf.step_coupled(
            &mut body,
            (0.0, -weight, 0.0),
            dt,
            ma,
            drag_est,
            confined_spin_stiffness(a, radius, nu),
            1e-5,
            30,
        );
        iters += info.iterations as f64;
        div = div.max(info.load.divergence);
        vmax = vmax.max(body.disc.uy.abs());
        if !(body.disc.cy.is_finite() && body.disc.cy.abs() < radius) || info.residual > 1e-2 {
            ok = false;
            break;
        }
        if (step + 1) % (steps / 5).max(1) == 0 {
            y.push(body.disc.cy);
        }
    }
    FallingDisc {
        steps,
        y,
        vy_max: vmax,
        iterations: iters / steps as f64,
        divergence: div,
        umax: 0.0,
        ok,
    }
}

/// Free disc driven by a constant external torque in a drum at rest: the steady spin rate is the
/// torque over the confined Couette stiffness. Returns the relative error of the spin rate at
/// `t_end`, the largest translation speed (must stay at round-off by symmetry) and the mean
/// iterations per step.
pub fn verify_free_spin(
    n: usize,
    a: f64,
    ratio: f64,
    torque: f64,
    nu: f64,
    t_end: f64,
    dt: f64,
) -> (f64, f64, f64) {
    let radius = 0.5;
    let steps = (t_end / dt).round() as usize;
    let mass = ratio * PI * a * a;
    let ma = confined_added_mass(a, radius);
    let kappa = confined_spin_stiffness(a, radius, nu);
    let mut bf = BodyFlow::new(n, 0.55, radius, nu);
    let mut body = Body {
        disc: Disc {
            cx: 0.0,
            cy: 0.0,
            r: a,
            ux: 0.0,
            uy: 0.0,
            omega: 0.0,
        },
        mass,
        inertia: 0.5 * mass * a * a,
        accel: (0.0, 0.0),
    };
    let (mut iters, mut vmax) = (0.0, 0.0f64);
    for _ in 0..steps {
        let info = bf.step_coupled(
            &mut body,
            (0.0, 0.0, torque),
            dt,
            ma,
            stokes_annulus_drag(a, radius, 1.0, nu),
            kappa,
            1e-8,
            30,
        );
        iters += info.iterations as f64;
        vmax = vmax.max(body.disc.ux.hypot(body.disc.uy));
    }
    (
        body.disc.omega / (torque / kappa) - 1.0,
        vmax,
        iters / steps as f64,
    )
}
