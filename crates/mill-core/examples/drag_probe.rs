//! Drag probe: hydrodynamic drag force on ONE ball pushed at a prescribed constant velocity U
//! through a quiescent, fully filled slurry pool in the shipped PBF solver, compared with the
//! 2D-cylinder reference per unit depth
//!   F_ref = 0.5 * rho * U^2 * d * C_D(Re),   C_D = 1 + 10 * Re^(-2/3)   (White's fit),
//!   Re = U d / nu,  nu = mu / rho.
//! The fit is only meaningful for Re > ~10; below Re = 5 the ratio is still printed but marked
//! 'ref n/a (Re<5)'.
//!
//! Setup: drum of radius `--rdrum` (not rotating, no lifters, slurry fill fraction 1.0), one ball
//! of diameter `--d`. Every step the ball's position/velocity are overwritten (it travels a straight
//! line of 1.2 R_drum through the centre, at most 4000 steps) once with +U and once with -U.
//! Drag is odd in U, so everything even in U (buoyancy, hydrostatics, lattice bias, pool noise)
//! cancels in the difference:
//!   F_meas = -(mean impulse_x(+U) - mean impulse_x(-U)) / (2 dt)   (force per unit depth, N/m).
//! The first max(2 d / U, 20 steps) are discarded as start-up; the rest is averaged, and the
//! standard error of F_meas comes from 4 block means of each run. Keep U >~ 0.1 m/s: a still
//! pool has particle noise of order 0.01 m/s that swamps smaller signals, so reach low Re by
//! raising `--nu` (U = Re nu / d), e.g. `--nu 1e-3 --re 10,100`.
//!
//! Resolution: the fluid lattice spacing is dx = rdrum / res, so `--d`, `--rdrum` and `--res`
//! together set d/dx (default 0.01 / (0.125/25) = 2). Results are reported as a function of d/dx.
//!
//! Per-mechanism split: `CouplingImpulses` / `FluidStepStats` have no per-mechanism accumulators
//! (only the total ball impulse, clamp count and total fluid momentum change), so only the total
//! is printed.
//!
//! Measured 2026-10-19 (docs/VERIFICATION.md "Ball drag measurement"): at the shipped resolution
//! (d/dx 2-4) the total drag is 5-17x the cylinder value for Re 1e3-3e3 whatever the form-drag
//! closure does, and within ~+-25 % of it at Re 1e4.
//!
//! `cargo run --release -p mill-core --example drag_probe -- --re 1000,10000 --res 25 --substeps 16`
//! Options: `--nu a,b` (m2/s, default 1e-6), `--re a,b` (default 10,100,1000,10000; U = Re nu / d),
//! `--res a,b` (25), `--substeps 8`, `--d 0.01`, `--rdrum 0.125`,
//! `--coarse K` (simulate radius K d/2 at the same density, true radius d/2), `--wet w` and
//! `--st sigma` (override slurry wettability / surface tension, to ablate adhesion and cohesion), `--csv`.

use glam::Vec2;
use mill_core::dem::DemState;
use mill_core::geometry::Drum;
use mill_core::params::{EffectiveMedia, LiftersParams, Params};
use mill_core::pbf::FluidParticles;

/// Comma-separated f32 list for `key`, or `default` when absent.
fn list(args: &[String], key: &str, default: &[f32]) -> Vec<f32> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1))
        .map(|v| v.split(',').filter_map(|s| s.trim().parse().ok()).collect())
        .unwrap_or_else(|| default.to_vec())
}

struct Cfg {
    nu: f32,
    wet: Option<f32>,
    st: Option<f32>,
    res: u32,
    substeps: u32,
    d: f32,
    rdrum: f32,
    coarse: f32,
}

/// Mean x impulse on the ball (kg m/s per metre depth) over the measured window, and its standard
/// error from 4 block means. `u` is signed: the ball travels in the direction of its sign.
fn mean_impulse_x(c: &Cfg, u: f32) -> (f64, f64) {
    let params = Params::default();
    let mut slurry = params.slurry;
    slurry.fill_fraction = 1.0;
    slurry.density_kg_m3 = 1000.0;
    slurry.viscosity_pa_s = 1000.0 * c.nu;
    if let Some(w) = c.wet {
        slurry.wettability = w;
    }
    if let Some(t) = c.st {
        slurry.surface_tension_n_m = t;
    }
    let iters = params.simulation.pbf_iterations;
    let dt = 1.0 / (60.0 * c.substeps as f32);
    let drum = Drum::new(
        c.rdrum,
        0.0,
        LiftersParams {
            count: 0,
            ..LiftersParams::default()
        },
    );
    let kc = c.coarse.max(1.0);
    let eff = EffectiveMedia {
        true_diameter_m: c.d,
        diameter_m: kc * c.d,
        density_kg_m3: 7800.0,
        ball_count: 1,
        scale_factor: kc,
    };
    let mut dem = DemState::new(&eff, c.rdrum, 1);
    assert_eq!(dem.balls.len(), 1, "expected exactly one ball");
    let r_sim = 0.5 * kc * c.d;
    let mut f = FluidParticles::seed_lattice(&slurry, c.rdrum, c.res, &[Vec2::ZERO], r_sim);

    // The ball crosses the pool in a straight line (no wrap-around: re-entering its own wake
    // would pollute the average). Total path 1.2 R_drum (at least 5 d of clearance to the wall at
    // both ends), at most MAX_STEPS steps; start-up = the first 2 d of travel (at least 20 steps).
    // The U = 0 baseline holds the ball at the start point for the same number of steps.
    const MAX_STEPS: usize = 4000;
    let travel = 1.2 * c.rdrum;
    let speed = u.abs();
    let total = ((travel / (speed * dt)).ceil() as usize).clamp(60, MAX_STEPS);
    let settle = ((2.0 * c.d / (speed * dt)).ceil() as usize)
        .max(20)
        .min(total / 2);
    let start = -0.5 * travel * u.signum();
    let mut samples: Vec<f64> = Vec::with_capacity(total);
    for s in 0..total {
        let t = s as f32 * dt;
        dem.balls.x[0] = Vec2::new(start + u * t, 0.0);
        dem.balls.v[0] = Vec2::new(u, 0.0);
        dem.balls.omega[0] = 0.0;
        let (imp, _st) = f.step_coupled(&drum, 0.0, &slurry, iters, dt, &dem.balls);
        if s >= settle {
            samples.push(imp.impulses[0].x as f64);
        }
    }
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len().max(1) as f64;
    let blocks: Vec<f64> = samples
        .chunks((samples.len() / 4).max(1))
        .map(mean)
        .collect();
    let m = mean(&samples);
    let var = blocks
        .iter()
        .map(|b| (b - mean(&blocks)).powi(2))
        .sum::<f64>()
        / (blocks.len().max(2) - 1) as f64;
    (m, (var / blocks.len().max(1) as f64).sqrt())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv = args.iter().any(|a| a == "--csv");
    let one = |key: &str, default: f32| list(&args, key, &[default])[0];
    let nus = list(&args, "--nu", &[1e-6]);
    let re_given = args.iter().any(|a| a == "--re");
    let res_list = list(&args, "--re", &[10.0, 100.0, 1000.0, 10000.0]);
    let _ = re_given; // U is always derived from Re (the default Re list applies when absent)
    let ress = list(&args, "--res", &[25.0]);
    let (d, rdrum) = (one("--d", 0.01), one("--rdrum", 0.125));
    let substeps = one("--substeps", 8.0) as u32;
    let coarse = one("--coarse", 1.0);
    let opt = |key: &str| args.iter().position(|a| a == key).map(|_| one(key, 0.0));
    let (wet, st) = (opt("--wet"), opt("--st"));
    let dt = 1.0 / (60.0 * substeps as f32);

    let header = "nu Re U res d/dx F_meas(+-se) F_ref R=F_meas/F_ref(+-se)";
    if csv {
        println!("nu,Re,U,res,d_over_dx,F_meas,F_se,F_ref,R,R_se,note");
    } else {
        println!("drag_probe d {d} rdrum {rdrum} substeps {substeps} coarse {coarse} dt {dt:.3e}");
        println!("{header}");
    }
    for &nu in &nus {
        for &re in &res_list {
            let u = re * nu / d;
            for &res in &ress {
                let dx = rdrum / res;
                {
                    let ddx = d / dx;
                    if u * dt / dx > 0.5 {
                        if csv {
                            println!("{nu:e},{re},{u:.4e},{res},{ddx:.2},,,,skipped (CFL)");
                        } else {
                            println!(
                                "{nu:.2e} {re:>7.1} {u:.3e} {res:>3} {ddx:>5.2} skipped (CFL)"
                            );
                        }
                        continue;
                    }
                    let c = Cfg {
                        nu,
                        wet,
                        st,
                        res: res as u32,
                        substeps,
                        d,
                        rdrum,
                        coarse,
                    };
                    let (fwd, se_f) = mean_impulse_x(&c, u);
                    let (bwd, se_b) = mean_impulse_x(&c, -u);
                    let f_meas = -(fwd - bwd) / (2.0 * dt as f64);
                    let f_se = (se_f.powi(2) + se_b.powi(2)).sqrt() / (2.0 * dt as f64);
                    let cd = 1.0 + 10.0 * (re as f64).powf(-2.0 / 3.0);
                    let f_ref = 0.5 * 1000.0 * (u as f64).powi(2) * d as f64 * cd;
                    let ratio = f_meas / f_ref;
                    let r_se = f_se / f_ref;
                    let note = if re < 5.0 { "ref n/a (Re<5)" } else { "" };
                    if csv {
                        println!(
                            "{nu:e},{re},{u:.4e},{res},{ddx:.2},{f_meas:.5e},{f_se:.2e},{f_ref:.5e},{ratio:.4},{r_se:.3},{note}"
                        );
                    } else {
                        println!(
                            "{nu:.2e} {re:>7.1} {u:.3e} {res:>3} {ddx:>5.2} {f_meas:>11.4e}(+-{f_se:.1e}) {f_ref:>11.4e} {ratio:>8.4}(+-{r_se:.3}) {note}"
                        );
                    }
                }
            }
        }
    }
}
