use super::*;
use crate::params::LiftersParams;

const DT: f32 = 1.0 / 960.0;

fn drum(radius: f32, omega: f32, lifters: u32) -> Drum {
    Drum::new(
        radius,
        omega,
        LiftersParams {
            count: lifters,
            ..LiftersParams::default()
        },
    )
}

fn slurry(fill: f32) -> SlurryParams {
    SlurryParams {
        fill_fraction: fill,
        ..SlurryParams::default()
    }
}

fn total_momentum(f: &Fluid) -> Vec2 {
    f.v.iter().fold(Vec2::ZERO, |a, &v| a + v) * f.particle_mass
}

fn angular_momentum_about(f: &Fluid, c: Vec2, vc: Vec2) -> f64 {
    f.x.iter()
        .zip(&f.v)
        .map(|(&x, &v)| {
            let r = x - c;
            let u = v - vc;
            (r.x * u.y - r.y * u.x) as f64
        })
        .sum::<f64>()
        * f.particle_mass as f64
}

/// A free blob of lattice particles well away from the wall (R = 0.1 m, blob radius 0.03 m).
fn free_blob() -> (Fluid, Drum, SlurryParams) {
    let d = drum(0.1, 0.0, 0);
    let s = slurry(0.5);
    let mut f = Fluid::new(&s, &d, 0.0, 16);
    let row_h = f.dx * 3.0f32.sqrt() * 0.5;
    let mut x = Vec::new();
    for row in -8i32..=8 {
        let off = if row.rem_euclid(2) == 0 { 0.0 } else { 0.5 };
        for col in -8i32..=8 {
            let p = Vec2::new((col as f32 + off) * f.dx, row as f32 * row_h);
            if p.length() < 0.03 {
                x.push(p);
            }
        }
    }
    f.v = vec![Vec2::ZERO; x.len()];
    f.dye = vec![0.0; x.len()];
    f.kappa = vec![0.0; x.len()];
    f.x = x;
    (f, d, s)
}

#[test]
fn lattice_interior_density_equals_rest_density() {
    let d = drum(0.0315, 0.0, 0);
    let s = slurry(0.9);
    let f = Fluid::new(&s, &d, 0.0, 30);
    let rho = f.densities();
    let mut checked = 0;
    for (p, r) in f.x.iter().zip(&rho) {
        if p.length() < 0.0315 - 3.0 * f.dx {
            // Away from the wall and from the (partial) free-surface row.
            if p.y < 0.0 {
                assert!(
                    (r / f.rest_density - 1.0).abs() < 5e-3,
                    "interior density ratio {}",
                    r / f.rest_density
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 100, "only {checked} interior particles checked");
}

#[test]
fn wall_adjacent_density_is_near_rest_density_once_settled() {
    // The seeded lattice does not follow the curved wall, so the wall layer has to settle
    // against the boundary particles first; then it must read bulk density, not a deficit or
    // pile-up.
    let d = drum(0.0315, 0.0, 0);
    let s = slurry(0.6);
    let mut f = Fluid::new(&s, &d, 0.0, 20);
    for _ in 0..240 {
        f.step(&d, 0.0, &s, DT);
    }
    let rho = f.densities();
    let (mut sum, mut count, mut worst) = (0.0f32, 0, 0.0f32);
    for (p, r) in f.x.iter().zip(&rho) {
        if 0.0315 - p.length() < 1.5 * f.dx && p.y < -0.01 {
            let e = r / f.rest_density - 1.0;
            sum += e;
            worst = worst.max(e.abs());
            count += 1;
        }
    }
    assert!(count > 20);
    assert!(
        (sum / count as f32).abs() < 0.03 && worst < 0.15,
        "wall layer density error mean {} worst {worst}",
        sum / count as f32
    );
}

#[test]
fn free_fall_blob_conserves_momentum_under_internal_forces() {
    let (mut f, d, s) = free_blob();
    // Internal shear so viscosity and pressure are active.
    for (v, x) in f.v.iter_mut().zip(&f.x) {
        *v = Vec2::new(-x.y, x.x) * 30.0 + Vec2::new(0.05 * x.x / 0.03, 0.0);
    }
    let p0 = total_momentum(&f);
    let steps = 24;
    for _ in 0..steps {
        f.step(&d, 0.0, &s, DT);
    }
    let expected_dp_y = f.particle_mass * f.len() as f32 * GRAVITY * DT * steps as f32;
    let dp = total_momentum(&f) - p0;
    let scale = f.particle_mass * f.len() as f32 * 0.05;
    assert!(
        (dp.x).abs() < 1e-3 * scale
            && (dp.y - expected_dp_y).abs() < 1e-3 * expected_dp_y.abs().max(scale),
        "dp {dp:?} expected y {expected_dp_y}"
    );
}

#[test]
fn free_fall_blob_conserves_angular_momentum_about_its_centre() {
    let (mut f, d, s) = free_blob();
    for (v, x) in f.v.iter_mut().zip(&f.x) {
        *v = Vec2::new(-x.y, x.x) * 40.0 + Vec2::new(0.02 * x.x / 0.03, 0.0);
    }
    let l0 = angular_momentum_about(&f, Vec2::ZERO, Vec2::ZERO);
    let steps = 24;
    for _ in 0..steps {
        f.step(&d, 0.0, &s, DT);
    }
    let n = f.len() as f32;
    let c = f.x.iter().fold(Vec2::ZERO, |a, &x| a + x) / n;
    let vc = f.v.iter().fold(Vec2::ZERO, |a, &v| a + v) / n;
    let l1 = angular_momentum_about(&f, c, vc);
    assert!(l0.abs() > 0.0);
    assert!(
        (l1 - l0).abs() < 2e-2 * l0.abs(),
        "L0 {l0} L1 {l1} (viscous decay allowed, but central forces cannot create or shift it)"
    );
}

#[test]
fn rigidly_rotating_blob_feels_no_viscous_force() {
    let (mut f, d, s) = free_blob();
    for (v, x) in f.v.iter_mut().zip(&f.x) {
        *v = Vec2::new(-x.y, x.x) * 20.0;
    }
    let w = f.prepare(&d, 0.0);
    let ke0 = kinetic_energy_f64(&f.v, f.particle_mass);
    let mut stats = FluidStats::default();
    let mut slurry = s;
    slurry.surface_tension_n_m = 0.0;
    f.viscosity_solve(&w, &slurry, DT, &mut stats);
    let ke1 = kinetic_energy_f64(&f.v, f.particle_mass);
    let _ = d;
    assert!(
        (ke1 - ke0).abs() < 1e-4 * ke0,
        "rigid rotation lost {} of {ke0} J to viscosity",
        ke0 - ke1
    );
}

#[test]
fn viscosity_is_a_no_op_at_zero_viscosity() {
    let (mut f, d, mut s) = free_blob();
    for (v, x) in f.v.iter_mut().zip(&f.x) {
        *v = Vec2::new(x.y, 0.0) * 10.0;
    }
    s.viscosity_pa_s = 0.0;
    let before = f.v.clone();
    let w = f.prepare(&d, 0.0);
    let mut stats = FluidStats::default();
    f.viscosity_solve(&w, &s, DT, &mut stats);
    assert_eq!(before, f.v);
    assert_eq!(stats.viscosity_iterations, 0);
}

#[test]
fn still_drum_stays_nearly_still_and_inside() {
    let d = drum(0.0315, 0.0, 0);
    let s = slurry(0.35);
    let mut f = Fluid::new(&s, &d, 0.0, 15);
    let mut hits = 0;
    for _ in 0..240 {
        let st = f.step(&d, 0.0, &s, DT);
        hits += st.wall_backstop_hits;
    }
    assert_eq!(hits, 0, "wall backstop fired");
    for p in &f.x {
        assert!(p.is_finite() && p.length() < 0.0315);
    }
    let vmax = f.v.iter().map(|v| v.length()).fold(0.0f32, f32::max);
    assert!(vmax < 0.05, "residual speed {vmax} m/s in a still pool");
}

#[test]
fn two_runs_are_bit_identical() {
    let d = drum(0.0315, 2.0, 0);
    let s = slurry(0.35);
    let run = || {
        let mut f = Fluid::new(&s, &d, 0.0, 12);
        let mut angle = 0.0f32;
        for _ in 0..60 {
            f.step(&d, angle, &s, DT);
            angle += d.omega * DT;
        }
        f.x.iter()
            .chain(&f.v)
            .flat_map(|p| [p.x.to_bits(), p.y.to_bits()])
            .collect::<Vec<u32>>()
    };
    assert_eq!(run(), run());
}

#[test]
fn a_non_finite_particle_is_reset() {
    let d = drum(0.0315, 0.0, 0);
    let s = slurry(0.35);
    let mut f = Fluid::new(&s, &d, 0.0, 12);
    f.x[3] = Vec2::new(f32::NAN, 0.0);
    f.step(&d, 0.0, &s, DT);
    assert!(f.x.iter().all(|p| p.is_finite()));
    assert!(f.v.iter().all(|v| v.is_finite()));
}

#[test]
fn a_lifter_drum_keeps_the_fluid_inside() {
    let d = drum(0.0315, 3.0, 4);
    let s = slurry(0.3);
    let mut f = Fluid::new(&s, &d, 0.0, 15);
    let mut angle = 0.0f32;
    for _ in 0..240 {
        f.step(&d, angle, &s, DT);
        angle += d.omega * DT;
    }
    for p in &f.x {
        let (sdf, _) = d.sdf_world(*p, angle);
        assert!(p.is_finite() && sdf > -1e-4, "particle outside: sdf {sdf}");
    }
}
