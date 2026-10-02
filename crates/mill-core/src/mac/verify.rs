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
