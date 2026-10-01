//! Solver-agnostic verification suite: checks the fluid solver and ball-fluid coupling against
//! analytic solutions. Only `harness` (in `common/verification_common.rs`, shared with the
//! `verification_probe` example) touches solver APIs, so a solver rewrite edits that module only.
//!
//! Tolerances live in `gates` and must not be loosened to make a solver pass.
//! Tests that fail on the current PBF solver are `#[ignore]`d with the measured numbers; run them
//! with `cargo test --release -p mill-core --test verification -- --ignored`.

#[path = "common/verification_common.rs"]
mod common;

use common::{analytic, gates, harness};
use std::sync::OnceLock;

/// Panics listing every `(label, value, limit)` with `value > limit` (or NaN).
fn assert_all_le(items: &[(String, f64, f64)]) {
    let bad: Vec<String> = items
        .iter()
        .filter(|(_, v, lim)| v.is_nan() || *v > *lim)
        .map(|(l, v, lim)| format!("{l}: {v:.5} > {lim}"))
        .collect();
    assert!(bad.is_empty(), "gates failed:\n  {}", bad.join("\n  "));
}

fn hydro(res: u32) -> &'static harness::Hydrostatic {
    static A: OnceLock<harness::Hydrostatic> = OnceLock::new();
    static B: OnceLock<harness::Hydrostatic> = OnceLock::new();
    match res {
        25 => A.get_or_init(|| harness::hydrostatic(25)),
        50 => B.get_or_init(|| harness::hydrostatic(50)),
        _ => unreachable!(),
    }
}

fn spin(res: u32) -> &'static harness::SpinUp {
    static A: OnceLock<harness::SpinUp> = OnceLock::new();
    static B: OnceLock<harness::SpinUp> = OnceLock::new();
    match res {
        25 => A.get_or_init(|| harness::spin_up(25)),
        50 => B.get_or_init(|| harness::spin_up(50)),
        _ => unreachable!(),
    }
}

fn couette(res: u32) -> &'static harness::TaylorCouette {
    static A: OnceLock<harness::TaylorCouette> = OnceLock::new();
    static B: OnceLock<harness::TaylorCouette> = OnceLock::new();
    match res {
        25 => A.get_or_init(|| harness::taylor_couette(25)),
        50 => B.get_or_init(|| harness::taylor_couette(50)),
        _ => unreachable!(),
    }
}

// ------------------------------------------------------------------ analytic helpers

#[test]
fn analytic_helpers_are_self_consistent() {
    use analytic::*;
    assert!(j0(2.404_825_557_695_773).abs() < 1e-7);
    assert!(j1(3.831_705_970_207_512_5).abs() < 1e-7);
    assert!((j0(30.0) - (-0.086_367_983_581_040_2)).abs() < 1e-7);
    let z = j1_zeros(5000);
    assert!((z[0] - 3.831_705_970_2).abs() < 1e-8);
    // sum 1/lambda_n^2 = 1/8 (tail ~ 1/(pi^2 N)).
    let s: f64 = z.iter().map(|l| 1.0 / (l * l)).sum();
    assert!((s - 0.125).abs() < 5e-5, "sum 1/lambda^2 = {s}");
    // L(0) = 0 and L(inf) = L_inf.
    assert!(spinup_angular_momentum_fraction_closed(0.0, 5000).abs() < 1e-3);
    assert!((spinup_angular_momentum_fraction_closed(5.0, 50) - 1.0).abs() < 1e-12);
    // Quadrature of the velocity series agrees with the closed form.
    for t in [0.005, 0.02, 0.05, 0.1, 0.2, 0.5] {
        let q = spinup_angular_momentum_fraction(t);
        let c = spinup_angular_momentum_fraction_closed(t, 2000);
        assert!((q - c).abs() < 2e-3, "t={t}: quadrature {q} vs closed {c}");
    }
    // The fluid velocity equals the wall velocity (s = 1) at the wall and relaxes to rigid
    // rotation (u = s) at late times.
    let zs = j1_zeros(400);
    assert!((spinup_velocity(1.0, 0.1, &zs) - 1.0).abs() < 1e-6);
    assert!((spinup_velocity(0.5, 10.0, &zs) - 0.5).abs() < 1e-12);
}

// ------------------------------------------------------------------ 1. hydrostatic

fn hydro_compression(res: u32) {
    let h = hydro(res);
    assert_all_le(&[
        (
            format!("res {res} worst-5% compression"),
            h.worst5_c,
            gates::HYDRO_WORST5_COMPRESSION,
        ),
        (
            format!("res {res} mean compression"),
            h.mean_c,
            gates::HYDRO_MEAN_COMPRESSION,
        ),
    ]);
}

fn hydro_ke(res: u32) {
    let h = hydro(res);
    assert_all_le(&[(
        format!("res {res} KE/(M g dx)"),
        h.ke_over_mgdx,
        gates::HYDRO_KE_OVER_MGDX,
    )]);
}

fn hydro_wall(res: u32) {
    let h = hydro(res);
    assert!(h.bins_used >= 4, "too few wall bins ({})", h.bins_used);
    assert_all_le(&[(
        format!("res {res} wall pressure L2 (bins {})", h.bins_used),
        h.wall_l2_rel,
        gates::HYDRO_WALL_PRESSURE_L2,
    )]);
}

#[test]
#[ignore = "PBF: worst-5% compression 0.270, mean 0.047 (gates 0.02 / 0.001)"]
fn hydrostatic_compression_res25() {
    hydro_compression(25);
}
#[test]
#[ignore = "PBF: worst-5% compression 0.914, mean 0.138 (gates 0.02 / 0.001)"]
fn hydrostatic_compression_res50() {
    hydro_compression(50);
}
#[test]
fn hydrostatic_kinetic_energy_res25() {
    hydro_ke(25);
}
#[test]
fn hydrostatic_kinetic_energy_res50() {
    hydro_ke(50);
}
#[test]
#[ignore = "PBF: wall pressure L2 error 0.229 (gate 0.05)"]
fn hydrostatic_wall_pressure_res25() {
    hydro_wall(25);
}
#[test]
#[ignore = "PBF: wall pressure L2 error 0.197 (gate 0.05)"]
fn hydrostatic_wall_pressure_res50() {
    hydro_wall(50);
}

// ------------------------------------------------------------------ 2. spin-up

fn spin_l(res: u32) {
    let s = spin(res);
    let items: Vec<_> = s
        .samples
        .iter()
        .map(|(t, ls, le)| {
            (
                format!(
                    "res {res} |L_sim-L_exact|/Linf @ t/tau={t:.3} (sim {ls:.4}, exact {le:.4})"
                ),
                (ls - le).abs(),
                gates::SPINUP_L_ERR,
            )
        })
        .collect::<Vec<_>>();
    assert_all_le(&items);
}

fn spin_torque(res: u32) {
    let s = spin(res);
    assert_all_le(&[(
        format!(
            "res {res} torque mismatch (no-gravity variant {:.4})",
            s.torque_mismatch_no_gravity
        ),
        s.torque_mismatch,
        gates::SPINUP_TORQUE_MISMATCH,
    )]);
}

#[test]
#[ignore = "PBF: |L-Lexact|/Linf 0.52/0.63/0.44/0.20 at t/tau 0.02/0.05/0.1/0.2 (gate 0.03)"]
fn spin_up_angular_momentum_res25() {
    spin_l(25);
}
#[test]
#[ignore = "PBF: |L-Lexact|/Linf 0.52/0.60/0.47/0.23 at t/tau 0.02/0.05/0.1/0.2 (gate 0.03)"]
fn spin_up_angular_momentum_res50() {
    spin_l(50);
}
#[test]
#[ignore = "PBF: wall torque vs dL/dt mismatch 0.468 (gate 0.03)"]
fn spin_up_wall_torque_balance_res25() {
    spin_torque(25);
}
#[test]
#[ignore = "PBF: wall torque vs dL/dt mismatch 0.451 (gate 0.03)"]
fn spin_up_wall_torque_balance_res50() {
    spin_torque(50);
}

// ------------------------------------------------------------------ 3. Taylor-Couette

fn tc_ball(res: u32) {
    let t = couette(res);
    assert!(t.torque_ball < 0.0, "torque on ball should oppose its spin");
    assert_all_le(&[(
        format!(
            "res {res} ball torque {:.5e} vs T {:.5e} (clamp hits {})",
            t.torque_ball, t.t_analytic, t.clamp_hits
        ),
        (t.torque_ball.abs() - t.t_analytic).abs() / t.t_analytic,
        gates::TC_TORQUE_REL,
    )]);
}

fn tc_wall(res: u32) {
    let t = couette(res);
    assert_all_le(&[(
        format!(
            "res {res} wall torque {:.5e} vs T {:.5e}",
            t.torque_wall, t.t_analytic
        ),
        (t.torque_wall.abs() - t.t_analytic).abs() / t.t_analytic,
        gates::TC_TORQUE_REL,
    )]);
}

#[test]
#[ignore = "PBF: ball torque ~1e-3 of analytic (rel err 0.999, gate 0.03)"]
fn taylor_couette_ball_torque_res25() {
    tc_ball(25);
}
#[test]
#[ignore = "PBF: ball torque ~1e-4 of analytic (rel err 1.000, gate 0.03)"]
fn taylor_couette_ball_torque_res50() {
    tc_ball(50);
}
#[test]
#[ignore = "PBF: wall torque 3% of analytic (rel err 0.966, gate 0.03)"]
fn taylor_couette_wall_torque_res25() {
    tc_wall(25);
}
#[test]
#[ignore = "PBF: wall torque 10% of analytic (rel err 0.895, gate 0.03)"]
fn taylor_couette_wall_torque_res50() {
    tc_wall(50);
}

// ------------------------------------------------------------------ 4. energy closure

#[test]
#[ignore = "PBF: |unattributed|/shaft_work 0.174 (gate 0.02); takes ~80 s"]
fn energy_closure_unattributed() {
    let b = harness::energy_closure(25);
    assert!(b.shaft_work_j > 0.0, "no shaft work");
    assert_all_le(&[(
        format!(
            "unattributed {:.4e} J vs shaft work {:.4e} J over {:.3} s",
            b.unattributed_j, b.shaft_work_j, b.elapsed_s
        ),
        (b.unattributed_j / b.shaft_work_j).abs(),
        gates::ENERGY_UNATTRIBUTED_FRAC,
    )]);
}

// ------------------------------------------------------------------ 5. dry substep convergence

#[test]
fn dry_substep_convergence() {
    let seeds = [1u64, 2];
    let mean = |s: u32| {
        seeds
            .iter()
            .map(|&seed| harness::dry_mean_power(s, seed).0)
            .sum::<f64>()
            / seeds.len() as f64
    };
    // `simulation.substeps` is validated to [1, 16], so 4/8/16 stand in for 8/16/32.
    let p16 = mean(16);
    let items: Vec<_> = [4u32, 8]
        .iter()
        .map(|&s| {
            let p = mean(s);
            (
                format!("substeps {s} power {p:.5e} W vs 16 -> {p16:.5e} W"),
                (p - p16).abs() / p16,
                gates::DRY_SUBSTEP_REL,
            )
        })
        .collect();
    assert_all_le(&items);
}
