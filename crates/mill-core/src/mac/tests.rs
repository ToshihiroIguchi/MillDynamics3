use super::verify::verify_manufactured;
use super::*;

#[test]
fn poisson_converges_second_order() {
    let (e1, s1, _) = verify_manufactured(32);
    let (e2, s2, _) = verify_manufactured(64);
    let (e3, s3, _) = verify_manufactured(128);
    let o12 = (e1 / e2).log2();
    let o23 = (e2 / e3).log2();
    eprintln!("poisson L2: {e1:.3e} {e2:.3e} {e3:.3e} orders {o12:.2} {o23:.2}");
    eprintln!(
        "iterations: {} {} {}",
        s1.iterations, s2.iterations, s3.iterations
    );
    assert!(o23 > 1.7, "order {o23}");
    assert!(
        s3.iterations <= s1.iterations + 6,
        "iteration count grows with grid size"
    );
}

#[test]
fn projection_is_divergence_free_to_roundoff() {
    let n = 96;
    let dom = CircleDomain::new(0.5, n, 0.05);
    let solver = PoissonSolver::new(&dom);
    let mut seed = 12345u64;
    let mut rnd = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    };
    let mut u: Vec<f64> = dom
        .wx
        .iter()
        .map(|&w| if w > 0.0 { rnd() } else { 0.0 })
        .collect();
    let mut v: Vec<f64> = dom
        .wy
        .iter()
        .map(|&w| if w > 0.0 { rnd() } else { 0.0 })
        .collect();
    let div = |u: &[f64], v: &[f64]| -> Vec<f64> {
        let mut d = vec![0.0; n * n];
        for j in 0..n {
            for i in 0..n {
                d[i + n * j] = dom.wx[i + 1 + (n + 1) * j] * u[i + 1 + (n + 1) * j]
                    - dom.wx[i + (n + 1) * j] * u[i + (n + 1) * j]
                    + dom.wy[i + n * (j + 1)] * v[i + n * (j + 1)]
                    - dom.wy[i + n * j] * v[i + n * j];
            }
        }
        d
    };
    let d0 = div(&u, &v);
    let rhs: Vec<f64> = d0.iter().map(|x| -x).collect();
    let mut p = vec![0.0; n * n];
    let st = solver.solve(&rhs, &mut p, 1e-14, 300);
    for j in 0..n {
        for i in 0..=n {
            if dom.wx[i + (n + 1) * j] > 0.0 {
                u[i + (n + 1) * j] -= p[i + n * j] - p[i - 1 + n * j];
            }
        }
    }
    for j in 0..=n {
        for i in 0..n {
            if dom.wy[i + n * j] > 0.0 {
                v[i + n * j] -= p[i + n * j] - p[i + n * (j - 1)];
            }
        }
    }
    let d1 = div(&u, &v);
    let m0 = d0.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let m1 = d1.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    eprintln!(
        "div before {m0:.3e} after {m1:.3e} ({} it, res {:.1e})",
        st.iterations, st.residual
    );
    assert!(m1 <= 1e-10 * m0.max(1.0));
}

#[test]
fn solver_is_deterministic() {
    let (a, ..) = verify_manufactured(48);
    let (b, ..) = verify_manufactured(48);
    assert_eq!(a.to_bits(), b.to_bits());
}
