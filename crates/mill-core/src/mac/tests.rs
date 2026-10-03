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

#[test]
fn variational_operator_annulus_residual() {
    use super::staggered::StaggeredMesh;
    use super::variational::ViscousOperator;
    let (r1, r2, omega) = (0.25, 0.5, 1.0);
    let mut previous: Option<f64> = None;
    for n in [64usize, 128] {
        let sm = StaggeredMesh::new(n, 0.55, |x, y| {
            let r = (x * x + y * y).sqrt();
            (r - r2).max(r1 - r)
        });
        let nn = n * n;
        let wall = move |x: f64, y: f64| {
            if (x * x + y * y).sqrt() < 0.5 * (r1 + r2) {
                (-omega * y, omega * x)
            } else {
                (0.0, 0.0)
            }
        };
        let ut = |r: f64| omega * r1 * r1 * (r2 * r2 / r - r) / (r2 * r2 - r1 * r1);
        let exact = |x: f64, y: f64| {
            let r = (x * x + y * y).sqrt();
            let s = ut(r) / r;
            (-s * y, s * x)
        };
        let frac: Vec<f64> = if std::env::var("FRAC_ONE").is_ok() {
            sm.cell_fraction
                .iter()
                .map(|&f| f64::from(f > 0.0))
                .collect()
        } else {
            sm.cell_fraction.clone()
        };
        let op = ViscousOperator::build(
            &sm,
            &sm.gu.fluid,
            &sm.gv.fluid,
            &sm.mesh.pressure_active,
            &frac,
            &wall,
        );
        let mut x = vec![0.0; 2 * nn];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                let (px, py) = sm.gu.centre(i, j);
                x[c] = exact(px, py).0;
                let (px, py) = sm.gv.centre(i, j);
                x[nn + c] = exact(px, py).1;
            }
        }
        let strain = op.strain(&x);
        let mut f = vec![0.0; 2 * nn];
        op.transpose(&strain, &mut f);
        // Residual force per node by distance to the nearest wall (in cells).
        let dx = sm.dx();
        let mut worst = [0.0f64; 4];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                if sm.gu.fluid[c] {
                    let (px, py) = sm.gu.centre(i, j);
                    let r = (px * px + py * py).sqrt();
                    let d = ((r - r1).min(r2 - r) / dx) as usize;
                    // Exact solution has zero stress force: compare with the viscous scale 1/dx^2 * u.
                    worst[d.min(3)] = worst[d.min(3)].max(f[c].abs() / (omega * r1));
                }
            }
        }
        // Active nodes that no valid point refers to carry no viscous coupling at all.
        let mut seen = vec![false; 2 * nn];
        let ones: Vec<[f64; 3]> = vec![[1.0; 3]; op.points()];
        let mut probe = vec![0.0; 2 * nn];
        op.transpose(&ones, &mut probe);
        for (k, v) in probe.iter().enumerate() {
            seen[k] = *v != 0.0;
        }
        let mut orphans = [0usize; 4];
        for j in 0..n {
            for i in 0..n {
                let c = i + n * j;
                for (lat, g) in [(0, &sm.gu), (1, &sm.gv)] {
                    if g.fluid[c] && !seen[c + lat * nn] {
                        let (px, py) = g.centre(i, j);
                        let r = (px * px + py * py).sqrt();
                        let d = ((r - r1).min(r2 - r) / dx) as usize;
                        orphans[d.min(3)] += 1;
                    }
                }
            }
        }
        // Stokes solve with the operator: A x = -B^T W g b0 (wall data), compare with the exact field.
        {
            use super::variational::pcg;
            let zero = vec![[0.0; 3]; op.points()];
            let mut rhs = vec![0.0; 2 * nn];
            op.transpose_with_wall(&zero, -1.0, &mut rhs);
            let mask: Vec<bool> = sm
                .gu
                .fluid
                .iter()
                .chain(sm.gv.fluid.iter())
                .copied()
                .collect();
            let mut sol = x.clone();
            for (k, m) in mask.iter().enumerate() {
                if *m {
                    sol[k] = 0.0;
                }
            }
            let apply = |v: &[f64], out: &mut [f64]| op.apply(v, out);
            let ident = |r: &[f64], z: &mut [f64]| z.copy_from_slice(r);
            let (it, res) = pcg(&apply, &ident, &mask, &rhs, &mut sol, 1e-10, 20000);
            let mut err = 0.0f64;
            let mut at = (0usize, 0.0, 0.0);
            for k in 0..2 * nn {
                if mask[k] && (sol[k] - x[k]).abs() > err {
                    err = (sol[k] - x[k]).abs();
                    let c = k % nn;
                    let g = if k < nn { &sm.gu } else { &sm.gv };
                    let (px, py) = g.centre(c % n, c / n);
                    at = (k / nn, px, py);
                }
            }
            let rr = (at.1 * at.1 + at.2 * at.2).sqrt();
            eprintln!(
                "worst node lattice {} at ({:.3},{:.3}) r={rr:.4} (wall at 0.25/0.5)",
                at.0, at.1, at.2
            );
            eprintln!(
                "n={n} Stokes solve: {it} its, residual {res:.1e}, max velocity error {err:.3e}"
            );
            assert!(res < 1e-8, "CG did not converge");
            let limit = if n == 64 { 0.06 } else { 0.02 };
            assert!(err < limit * omega * r1 * 4.0, "n={n}: error {err}");
            if let Some(p) = previous {
                assert!(err < 0.6 * p, "error must fall with the grid: {p} -> {err}");
            }
            previous = Some(err);
        }
        eprintln!("orphan active nodes by distance 0,1,2,>=3: {orphans:?}");
        eprintln!("n={n} residual force/(omega r1) by distance 0,1,2,>=3 cells: {worst:?}");
    }
}
