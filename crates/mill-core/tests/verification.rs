//! Solver-agnostic verification suite: checks the fluid solver and ball-fluid coupling against
//! analytic solutions. Only `harness` (in `common/verification_common.rs`, shared with the
//! `verification_probe` example) touches solver APIs, so a solver rewrite edits that module only.
//!
//! Tolerances live in `gates` and must not be loosened to make a solver pass.
//! Tests that fail on the current PBF solver are `#[ignore]`d with the measured numbers; run them
//! with `cargo test --release -p mill-core --test verification -- --ignored`.

#[path = "common/verification_common.rs"]
mod common;

use common::harness::Solver;
use common::{analytic, gates, harness};
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

/// Panics listing every `(label, value, limit)` with `value > limit` (or NaN).
fn assert_all_le(items: &[(String, f64, f64)]) {
    let bad: Vec<String> = items
        .iter()
        .filter(|(_, v, lim)| v.is_nan() || *v > *lim)
        .map(|(l, v, lim)| format!("{l}: {v:.5} > {lim}"))
        .collect();
    assert!(bad.is_empty(), "gates failed:\n  {}", bad.join("\n  "));
}

/// Memoises one expensive harness run per key so every gate test on the same run shares it.
type Cells<T> = Mutex<HashMap<(Solver, u32), &'static OnceLock<T>>>;
struct Cache<T: 'static>(OnceLock<Cells<T>>);

impl<T: 'static> Cache<T> {
    const fn new() -> Self {
        Self(OnceLock::new())
    }
    fn get(&'static self, key: (Solver, u32), run: impl FnOnce() -> T) -> &'static T {
        let cell: &'static OnceLock<T> = {
            let mut map = self.0.get_or_init(Default::default).lock().unwrap();
            map.entry(key)
                .or_insert_with(|| Box::leak(Box::new(OnceLock::new())))
        };
        cell.get_or_init(run)
    }
}

fn hydro(solver: Solver, res: u32) -> &'static harness::Hydrostatic {
    static C: Cache<harness::Hydrostatic> = Cache::new();
    C.get((solver, res), || harness::hydrostatic(solver, res))
}

fn spin(solver: Solver, res: u32) -> &'static harness::SpinUp {
    static C: Cache<harness::SpinUp> = Cache::new();
    C.get((solver, res), || harness::spin_up(solver, res))
}

fn rot(solver: Solver, res: u32) -> &'static harness::RotatingDrum {
    static C: Cache<harness::RotatingDrum> = Cache::new();
    C.get((solver, res), || harness::rotating_drum(solver, res))
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

fn hydro_compression(solver: Solver, res: u32) {
    let h = hydro(solver, res);
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

fn hydro_ke(solver: Solver, res: u32) {
    let h = hydro(solver, res);
    assert_all_le(&[(
        format!("res {res} KE/(M g dx)"),
        h.ke_over_mgdx,
        gates::HYDRO_KE_OVER_MGDX,
    )]);
}

fn hydro_wall(solver: Solver, res: u32) {
    let h = hydro(solver, res);
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
    hydro_compression(Solver::Pbf, 25);
}
#[test]
#[ignore = "PBF: worst-5% compression 0.914, mean 0.138 (gates 0.02 / 0.001)"]
fn hydrostatic_compression_res50() {
    hydro_compression(Solver::Pbf, 50);
}
#[test]
fn hydrostatic_kinetic_energy_res25() {
    hydro_ke(Solver::Pbf, 25);
}
#[test]
fn hydrostatic_kinetic_energy_res50() {
    hydro_ke(Solver::Pbf, 50);
}
#[test]
#[ignore = "PBF: wall pressure L2 error 0.229 (gate 0.05)"]
fn hydrostatic_wall_pressure_res25() {
    hydro_wall(Solver::Pbf, 25);
}
#[test]
#[ignore = "PBF: wall pressure L2 error 0.197 (gate 0.05)"]
fn hydrostatic_wall_pressure_res50() {
    hydro_wall(Solver::Pbf, 50);
}

// ------------------------------------------------------------------ 2. spin-up

fn spin_l_items(solver: Solver, res: u32) -> Vec<(String, f64, f64)> {
    let s = spin(solver, res);
    s.samples
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
        .collect::<Vec<_>>()
}

fn spin_l(solver: Solver, res: u32) {
    assert_all_le(&spin_l_items(solver, res));
}

fn spin_torque_items(solver: Solver, res: u32) -> Vec<(String, f64, f64)> {
    let s = spin(solver, res);
    vec![(
        format!(
            "res {res} torque mismatch (no-gravity variant {:.4})",
            s.torque_mismatch_no_gravity
        ),
        s.torque_mismatch,
        gates::SPINUP_TORQUE_MISMATCH,
    )]
}

fn spin_torque(solver: Solver, res: u32) {
    assert_all_le(&spin_torque_items(solver, res));
}

#[test]
#[ignore = "PBF: |L-Lexact|/Linf 0.52/0.63/0.44/0.20 at t/tau 0.02/0.05/0.1/0.2 (gate 0.03)"]
fn spin_up_angular_momentum_res25() {
    spin_l(Solver::Pbf, 25);
}
#[test]
#[ignore = "PBF: |L-Lexact|/Linf 0.52/0.60/0.47/0.23 at t/tau 0.02/0.05/0.1/0.2 (gate 0.03)"]
fn spin_up_angular_momentum_res50() {
    spin_l(Solver::Pbf, 50);
}
#[test]
#[ignore = "PBF: wall torque vs dL/dt mismatch 0.468 (gate 0.03)"]
fn spin_up_wall_torque_balance_res25() {
    spin_torque(Solver::Pbf, 25);
}
#[test]
#[ignore = "PBF: wall torque vs dL/dt mismatch 0.451 (gate 0.03)"]
fn spin_up_wall_torque_balance_res50() {
    spin_torque(Solver::Pbf, 50);
}

// ------------------------------------------------------------------ DFSPH: cases 1 and 2

fn hydro_all_gates(res: u32) {
    let h = hydro(Solver::Dfsph, res);
    assert!(h.bins_used >= 4, "too few wall bins ({})", h.bins_used);
    // Iteration-cap and wall-backstop counts are reported by the probe but not gated: the outcome
    // they would guard (density error, penetration energy) is gated directly below / in the
    // ledger, and thin films pressed against the wall legitimately trip the backstop.
    let mut items = Vec::new();
    items.extend([
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
        (
            format!("res {res} KE/(M g dx)"),
            h.ke_over_mgdx,
            gates::HYDRO_KE_OVER_MGDX,
        ),
        (
            format!("res {res} wall pressure L2 (bins {})", h.bins_used),
            h.wall_l2_rel,
            gates::HYDRO_WALL_PRESSURE_L2,
        ),
    ]);
    assert_all_le(&items);
}

fn spin_all_gates(res: u32) {
    let mut items = Vec::new();
    items.extend(spin_l_items(Solver::Dfsph, res));
    items.extend(spin_torque_items(Solver::Dfsph, res));
    assert_all_le(&items);
}

#[test]
fn dfsph_hydrostatic_res25() {
    hydro_all_gates(25);
}
#[test]
fn dfsph_hydrostatic_res50() {
    hydro_all_gates(50);
}
#[test]
fn dfsph_spin_up_res25() {
    spin_all_gates(25);
}
#[test]
fn dfsph_spin_up_res50() {
    spin_all_gates(50);
}

// ------------------------------------------------------------------ 6. rotating drum, fluid only

fn rot_convergence_items(solver: Solver, coarse: u32, fine: u32) -> Vec<(String, f64, f64)> {
    let (a, b) = (rot(solver, coarse), rot(solver, fine));
    let rel = |x: f64, y: f64| (y - x).abs() / x.abs();
    vec![
        (
            format!(
                "power {:.5} W (res {coarse}) -> {:.5} W (res {fine})",
                a.power_w, b.power_w
            ),
            rel(a.power_w, b.power_w),
            gates::ROT_CONVERGENCE_REL,
        ),
        (
            format!(
                "torque {:.5e} N m (res {coarse}) -> {:.5e} N m (res {fine})",
                a.torque_nm, b.torque_nm
            ),
            rel(a.torque_nm, b.torque_nm),
            gates::ROT_CONVERGENCE_REL,
        ),
    ]
}

fn rot_convergence(solver: Solver, coarse: u32, fine: u32) {
    assert_all_le(&rot_convergence_items(solver, coarse, fine));
}

fn rot_power_impulse(solver: Solver, res: u32) {
    let d = rot(solver, res);
    assert_all_le(&[(
        format!(
            "res {res} power {:.5} W vs omega*T_impulse {:.5} W",
            d.power_w, d.power_from_impulse_w
        ),
        d.power_mismatch,
        gates::ROT_POWER_IMPULSE_REL,
    )]);
}

fn rot_dfsph_convergence(coarse: u32, fine: u32) {
    assert_all_le(&rot_convergence_items(Solver::Dfsph, coarse, fine));
}

/// Slow in release (res 80 has ~4x the particles of res 40 and the settle + 1 s window is
/// ~700 outer steps with several internal steps each).
#[test]
#[ignore = "DFSPH gate not yet met (first-order convergence): power 0.1742 -> 0.2408 W, torque 1.645e-2 -> 2.275e-2 N m (change 0.383, gate 0.03); 8 rev settle + 3 rev window, about 5 min"]
fn dfsph_rotating_drum_convergence_25_50() {
    rot_dfsph_convergence(25, 50);
}
/// Slow in release (see `dfsph_rotating_drum_convergence_25_50`).
#[test]
#[ignore = "DFSPH gate not yet met (first-order convergence): power 0.2210 -> 0.2892 W, torque 2.088e-2 -> 2.731e-2 N m (change 0.308, gate 0.03); about 20 min"]
fn dfsph_rotating_drum_convergence_40_80() {
    rot_dfsph_convergence(40, 80);
}
#[test]
fn dfsph_rotating_drum_power_matches_impulse() {
    for res in [25, 40] {
        rot_power_impulse(Solver::Dfsph, res);
    }
}

#[test]
#[ignore = "PBF (expected not to converge): power 6.592 -> 5.398 W, torque 0.6225 -> 0.5098 N m, change 0.181 (gate 0.03); res 40 -> 80 power 5.785 -> 4.788 W"]
fn pbf_rotating_drum_convergence_25_50() {
    rot_convergence(Solver::Pbf, 25, 50);
}

// ------------------------------------------------------------------ DFSPH: fluid-only energy closure

fn ledger_gate(label: &str, l: &harness::FluidLedger) {
    assert!(l.wall_work_j > 0.0, "{label}: no wall work");
    assert_all_le(&[(
        format!(
            "{label}: unattributed {:.4e} J vs wall work {:.4e} J (dE_mech {:.4e}, viscous dissipation {:.4e})",
            l.unattributed_j, l.wall_work_j, l.d_mech_j, l.viscous_dissipation_j
        ),
        l.unattributed_frac,
        gates::FLUID_ENERGY_UNATTRIBUTED_FRAC,
    )]);
}

#[test]
#[ignore = "DFSPH gate not yet met: |unattributed|/wall_work 0.346 (gate 0.02): dE_mech 7.3e-3 J, wall work 3.93e-1 J, viscous dissipation 5.22e-1 J (res 80: 0.259). The second density solve re-adds KE the implicit viscosity removed"]
fn dfsph_energy_closure_rotating_drum_res40() {
    let d = rot(Solver::Dfsph, 40);
    ledger_gate(
        "rotating drum res 40",
        d.ledger.as_ref().expect("DFSPH ledger"),
    );
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
