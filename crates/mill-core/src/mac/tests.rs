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

#[test]
fn half_disc_sloshing_reference_is_converged() {
    let a = super::reference::half_disc_sloshing(12, 1)[0];
    let b = super::reference::half_disc_sloshing(16, 1)[0];
    assert!((a - b).abs() < 1e-5, "{a} vs {b}");
    assert!((b - 1.355727).abs() < 1e-5, "{b}");
}

#[test]
fn level_set_disc_rotates_with_small_volume_drift() {
    let r = super::verify::verify_levelset_rotation(96, false, 1.0, 2);
    assert!(r.volume_drift.abs() < 2e-3, "{r:?}");
    assert!(r.shape_l1 < 0.05, "{r:?}");
}

#[test]
fn still_pool_stays_at_rest() {
    let r = super::verify::verify_still_pool(48, -0.1, 1e-6, 0.3);
    assert!(r.spurious < 1e-4, "{r:?}");
    assert!(r.pressure_err < 0.05, "{r:?}");
    assert!(r.volume_drift.abs() < 1e-4, "{r:?}");
}

#[test]
fn submerged_disc_feels_the_displaced_weight() {
    let r = super::verify::verify_buoyancy(64, 0.1, 0.1);
    assert!(r.force_err.abs() < 5e-3, "{r:?}");
    assert!(r.torque.abs() < 1e-2, "{r:?}");
}

#[test]
fn translating_disc_drag_matches_stokes_annulus() {
    use super::verify::{stokes_annulus_drag, verify_moving_disc};
    let reference = stokes_annulus_drag(0.2, 0.5, 0.01, 1.0);
    let frozen = verify_moving_disc(64, 0.2, 0.0, 0.01, 0.0, 1.0, 0.6, 0.02, false);
    assert!((frozen.drag / reference - 1.0).abs() < 1e-2, "{frozen:?}");
    // The moving disc (mesh rebuilt every step) must reproduce the frozen-position drag.
    let still = verify_moving_disc(64, 0.2, 0.03, 0.1, 0.0, 1.0, 0.3, 0.01, false);
    let moved = verify_moving_disc(64, 0.2, 0.0, 0.1, 0.0, 1.0, 0.3, 0.01, true);
    assert!(moved.drag_noise < 0.05, "{moved:?}");
    assert!(
        (moved.drag_mean / still.drag_mean - 1.0).abs() < 2e-2,
        "{moved:?} vs {still:?}"
    );
}

#[test]
fn light_free_disc_reaches_the_stokes_terminal_speed() {
    // Density ratio 0.5 (added mass 3x the body mass): plain fixed-point coupling diverges.
    let (err, iters) = super::verify::verify_terminal_stokes(48, 0.2, 0.5, 1.0, 1.0, 0.05, 0.005);
    assert!(err.abs() < 1e-2, "{err}");
    assert!(iters < 6.0, "{iters}");
}

#[test]
fn squeeze_references_and_model_agree() {
    use super::lubrication::{disc_curvature, squeeze_coefficient, wall_curvature};
    use super::verify::{disc_pair_squeeze_force, eccentric_squeeze_force, stokes_annulus_drag};
    // Exact eccentric solution reduces to the concentric Stokes drag.
    let conc = stokes_annulus_drag(0.2, 0.5, 1.0, 1.0);
    let ecc = eccentric_squeeze_force(0.2, 0.5, 1e-6, 1.0, 1.0);
    assert!((ecc / conc - 1.0).abs() < 1e-4, "{ecc} vs {conc}");
    // Disc against the drum wall: model over exact, eps = h / R_eff up to 0.8.
    for a in [0.03, 0.05] {
        let b = 0.5;
        let r_eff = a * b / (b - a);
        for eps in [0.01, 0.1, 0.4, 0.8] {
            let h = eps * r_eff;
            let exact = eccentric_squeeze_force(a, b, b - a - h, 1.0, 1.0);
            let model = squeeze_coefficient(disc_curvature(a), wall_curvature(b), h, 1.0);
            assert!(
                (model / exact - 1.0).abs() < 0.01,
                "wall a={a} eps={eps}: {}",
                model / exact
            );
        }
    }
    // Two equal discs (relative approach speed is twice the speed of each).
    let a = 1.0;
    for eps in [0.01, 0.1, 0.4, 0.8] {
        let h = eps * a / 2.0;
        let exact = disc_pair_squeeze_force(a, h, 1.0, 1.0);
        let model = 2.0 * squeeze_coefficient(disc_curvature(a), disc_curvature(a), h, 1.0);
        assert!(
            (model / exact - 1.0).abs() < 0.01,
            "pair eps={eps}: {}",
            model / exact
        );
    }
}

#[test]
fn bingham_annular_couette_holds_analytic_profile() {
    use super::rheology::HerschelBulkley;
    use super::verify::verify_couette_hb;
    let law = HerschelBulkley::bingham(0.03, 0.01);
    let res = verify_couette_hb(64, law, 0.02, 3.0, 0.5, false, 1);
    assert!(res.r_yield_exact < 0.5, "case must contain a plug");
    assert!(res.u_err < 0.03, "velocity error {}", res.u_err);
    assert!(
        res.c_mean_err.abs() < 0.01,
        "torque error {}",
        res.c_mean_err
    );
}
