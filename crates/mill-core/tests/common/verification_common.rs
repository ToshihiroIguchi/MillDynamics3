//! Shared code for the solver-agnostic verification suite (`tests/verification.rs`) and the
//! `verification_probe` example (which includes this file via `#[path]`).
//!
//! * `harness` is the ONLY module that touches solver APIs (`FluidParticles`, `Simulation`, ...).
//!   Swapping the fluid solver means editing `harness` only.
//! * `analytic` holds the closed-form references (pure maths, no solver types).
//! * `gates` holds every tolerance, shared by the tests and the probe so they cannot drift.
//! * `report` turns harness results into flat rows for the probe.
//!
//! All quantities are SI and per metre of mill depth (2D).
#![allow(dead_code)]

/// Tolerances. Never loosen these to make the current solver pass.
pub mod gates {
    pub const HYDRO_WORST5_COMPRESSION: f64 = 0.02;
    pub const HYDRO_MEAN_COMPRESSION: f64 = 0.001;
    pub const HYDRO_KE_OVER_MGDX: f64 = 0.05;
    pub const HYDRO_WALL_PRESSURE_L2: f64 = 0.05;
    pub const SPINUP_L_ERR: f64 = 0.03;
    pub const SPINUP_TORQUE_MISMATCH: f64 = 0.03;
    pub const TC_TORQUE_REL: f64 = 0.03;
    pub const ENERGY_UNATTRIBUTED_FRAC: f64 = 0.02;
    pub const DRY_SUBSTEP_REL: f64 = 0.03;
    /// Fractions of `R^2/nu` at which the spin-up angular momentum is gated.
    pub const SPINUP_FRACTIONS: [f64; 4] = [0.02, 0.05, 0.1, 0.2];
}

pub mod analytic {
    use std::f64::consts::PI;

    fn bessel_series(nu: u32, x: f64) -> f64 {
        // sum_k (-1)^k / (k! (k+nu)!) (x/2)^(2k+nu)
        let half = 0.5 * x;
        let mut term = half.powi(nu as i32);
        for i in 1..=nu {
            term /= i as f64;
        }
        let mut sum = term;
        for k in 1..200u32 {
            term *= -half * half / (k as f64 * (k + nu) as f64);
            sum += term;
            if term.abs() < 1e-18 * sum.abs().max(1e-300) {
                break;
            }
        }
        sum
    }

    fn bessel_hankel(nu: u32, x: f64) -> f64 {
        // Hankel asymptotic expansion, 9 terms; accurate to ~1e-8 for x >= 12.
        let mu = 4.0 * (nu * nu) as f64;
        let (mut p, mut q) = (0.0, 0.0);
        let mut a = 1.0;
        for k in 0..9u32 {
            if k > 0 {
                let j = k as f64;
                a *= (mu - (2.0 * j - 1.0) * (2.0 * j - 1.0)) / (j * 8.0 * x);
            }
            match k % 4 {
                0 => p += a,
                1 => q += a,
                2 => p -= a,
                _ => q -= a,
            }
        }
        let chi = x - (0.5 * nu as f64 + 0.25) * PI;
        (2.0 / (PI * x)).sqrt() * (p * chi.cos() - q * chi.sin())
    }

    pub fn j0(x: f64) -> f64 {
        let x = x.abs();
        if x < 12.0 {
            bessel_series(0, x)
        } else {
            bessel_hankel(0, x)
        }
    }

    pub fn j1(x: f64) -> f64 {
        let s = x.signum();
        let x = x.abs();
        s * if x < 12.0 {
            bessel_series(1, x)
        } else {
            bessel_hankel(1, x)
        }
    }

    /// First `n` positive zeros of `J1` (Newton from McMahon's guess).
    pub fn j1_zeros(n: usize) -> Vec<f64> {
        (1..=n)
            .map(|k| {
                let beta = (k as f64 + 0.25) * PI;
                let mut x = beta - 3.0 / (8.0 * beta);
                for _ in 0..30 {
                    let d = j0(x) - j1(x) / x;
                    let dx = j1(x) / d;
                    x -= dx;
                    if dx.abs() < 1e-13 {
                        break;
                    }
                }
                x
            })
            .collect()
    }

    /// Impulsively started rigid rotation of the cylinder wall (no-slip, radius `R`), fluid
    /// initially at rest. With `s = r/R`, `T = nu t / R^2` and `lambda_n` the positive zeros of
    /// `J1`:
    ///
    ///   u_theta / (omega0 R) = s + 2 sum_n J1(lambda_n s) / (lambda_n J0(lambda_n)) exp(-lambda_n^2 T)
    ///
    /// (derivation: `u = omega0 r - v`, `v` solves the axisymmetric diffusion equation with
    /// `v(R)=0`, `v(r,0)=omega0 r`; the Fourier-Bessel coefficient of `s` in `J1(lambda_n s)` is
    /// `2 / (lambda_n J2(lambda_n)) = -2 / (lambda_n J0(lambda_n))` because `J1(lambda_n)=0`
    /// gives `J2 = -J0`).
    pub fn spinup_velocity(s: f64, t_nu: f64, zeros: &[f64]) -> f64 {
        let mut u = s;
        for &l in zeros {
            let e = (-l * l * t_nu).exp();
            if e < 1e-18 {
                break;
            }
            u += 2.0 * j1(l * s) / (l * j0(l)) * e;
        }
        u
    }

    /// `L(t) / L_inf` for a full disc, by midpoint quadrature in `r` of
    /// `integral rho u r dA` normalised by rigid rotation (`L_inf/(...) = 1/4` in these units).
    pub fn spinup_angular_momentum_fraction(t_nu: f64) -> f64 {
        let zeros = j1_zeros(400);
        let m = 2000usize;
        let mut sum = 0.0;
        for i in 0..m {
            let s = (i as f64 + 0.5) / m as f64;
            sum += spinup_velocity(s, t_nu, &zeros) * s * s;
        }
        4.0 * sum / m as f64
    }

    /// Closed form of the same quantity: `1 - 8 sum_n exp(-lambda_n^2 T) / lambda_n^2`.
    pub fn spinup_angular_momentum_fraction_closed(t_nu: f64, n_terms: usize) -> f64 {
        1.0 - 8.0
            * j1_zeros(n_terms)
                .iter()
                .map(|l| (-l * l * t_nu).exp() / (l * l))
                .sum::<f64>()
    }

    /// Taylor-Couette torque per metre of depth, inner cylinder `r1` rotating at `omega1`, outer
    /// cylinder `r2` at rest: `T = 4 pi mu omega1 r1^2 r2^2 / (r2^2 - r1^2)`.
    pub fn taylor_couette_torque(mu: f64, omega1: f64, r1: f64, r2: f64) -> f64 {
        4.0 * PI * mu * omega1 * r1 * r1 * r2 * r2 / (r2 * r2 - r1 * r1)
    }
}

/// The only module that touches solver types.
pub mod harness {
    use glam::Vec2;
    use mill_core::dem::DemState;
    use mill_core::geometry::Drum;
    use mill_core::params::{EffectiveMedia, LiftersParams, Params, SpeedMode};
    use mill_core::pbf::{FluidParticles, WALL_IMPULSE_BINS};
    use mill_core::{EnergyBudget, Simulation};

    use super::analytic;
    use super::gates;

    pub const DT: f32 = 1.0 / 480.0;
    const G: f64 = 9.81;
    /// `R^2 / nu` used by the spin-up and Taylor-Couette cases (s).
    pub const TAU_S: f64 = 0.25;
    /// The impulsive start is resolved with a 4x finer step than the production sub-step.
    pub const SPINUP_DT: f32 = DT / 4.0;

    fn make_drum(radius: f32, omega: f32) -> Drum {
        Drum::new(
            radius,
            omega,
            LiftersParams {
                count: 0,
                ..LiftersParams::default()
            },
        )
    }

    // ---------------------------------------------------------------- 1. hydrostatic

    #[derive(Debug, Clone)]
    pub struct Hydrostatic {
        pub res: u32,
        pub n: usize,
        pub mean_c: f64,
        pub worst5_c: f64,
        pub ke_over_mgdx: f64,
        pub wall_l2_rel: f64,
        pub bins_used: usize,
        pub y_surface: f64,
    }

    pub fn hydrostatic(res: u32) -> Hydrostatic {
        let params = Params::default();
        let r = params.mill.radius_m();
        let slurry = params.slurry; // fill 0.35, default viscosity
        let iters = params.simulation.pbf_iterations;
        let drum = make_drum(r, 0.0);
        let mut f = FluidParticles::seed_lattice(&slurry, r, res, &[], 0.0);
        let dx = (r / res as f32) as f64;
        for _ in 0..(1.0 / DT) as usize {
            f.step(&drum, 0.0, &slurry, iters, DT);
        }
        let window = (0.5 / DT) as usize;
        let mass_total = f.particle_mass as f64 * f.len() as f64;
        let (mut mean_acc, mut worst_acc, mut samples) = (0.0, 0.0, 0usize);
        let mut ke_acc = 0.0;
        let mut bins = [0.0f64; WALL_IMPULSE_BINS];
        for k in 0..window {
            let st = f.step(&drum, 0.0, &slurry, iters, DT);
            for (b, v) in bins.iter_mut().zip(st.wall_normal_impulse_by_angle) {
                *b += v as f64;
            }
            let ke: f64 = 0.5
                * f.particle_mass as f64
                * f.v.iter().map(|v| v.length_squared() as f64).sum::<f64>();
            ke_acc += ke;
            if k % 8 == 0 {
                let rho0 = f.rest_density;
                let mut c: Vec<f64> = f
                    .densities()
                    .iter()
                    .map(|r| (r / rho0 - 1.0).max(0.0) as f64)
                    .collect();
                mean_acc += c.iter().sum::<f64>() / c.len() as f64;
                c.sort_by(|a, b| b.partial_cmp(a).unwrap());
                let top = &c[..(c.len() / 20).max(1)];
                worst_acc += top.iter().sum::<f64>() / top.len() as f64;
                samples += 1;
            }
        }
        let ke_mean = ke_acc / window as f64;
        // Free-surface height: 98th percentile of particle y.
        let mut ys: Vec<f32> = f.x.iter().map(|p| p.y).collect();
        ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let y_s = ys[((ys.len() as f32 * 0.98) as usize).min(ys.len() - 1)] as f64;

        let t_win = window as f64 * DT as f64;
        let rho0 = f.rest_density as f64;
        let dth = 2.0 * std::f64::consts::PI / WALL_IMPULSE_BINS as f64;
        let (mut num, mut den, mut used) = (0.0, 0.0, 0usize);
        for (k, b) in bins.iter().enumerate() {
            let th0 = -std::f64::consts::PI + k as f64 * dth;
            let th1 = th0 + dth;
            let mut y_max = (r as f64 * th0.sin()).max(r as f64 * th1.sin());
            if th0 < std::f64::consts::FRAC_PI_2 && th1 > std::f64::consts::FRAC_PI_2 {
                y_max = r as f64;
            }
            if y_max >= y_s - 2.0 * f.h as f64 {
                continue;
            }
            let th_c = th0 + 0.5 * dth;
            let depth = y_s - r as f64 * th_c.sin();
            let f_an = rho0 * G * depth.max(0.0) * r as f64 * dth;
            let f_sim = b / t_win;
            num += (f_sim - f_an).powi(2);
            den += f_an * f_an;
            used += 1;
        }
        Hydrostatic {
            res,
            n: f.len(),
            mean_c: mean_acc / samples as f64,
            worst5_c: worst_acc / samples as f64,
            ke_over_mgdx: ke_mean / (mass_total * G * dx),
            wall_l2_rel: if den > 0.0 {
                (num / den).sqrt()
            } else {
                f64::NAN
            },
            bins_used: used,
            y_surface: y_s,
        }
    }

    // ---------------------------------------------------------------- 2. spin-up

    #[derive(Debug, Clone)]
    pub struct SpinUp {
        pub res: u32,
        pub n: usize,
        /// Particle-sum `L_inf = omega0 * sum m r^2` (actual lattice).
        pub l_inf: f64,
        /// `(t / tau, L_sim / L_inf, L_exact / L_inf)`.
        pub samples: Vec<(f64, f64, f64)>,
        /// `|dL - sum(J_wall) - sum(J_gravity)| / |dL|` over [0.02, 0.2] tau.
        pub torque_mismatch: f64,
        /// Same without the gravity-torque correction.
        pub torque_mismatch_no_gravity: f64,
    }

    pub fn spin_up(res: u32) -> SpinUp {
        let params = Params::default();
        let r = params.mill.radius_m();
        let mut slurry = params.slurry;
        // `validate()` caps the fill at 0.9, but `seed_lattice` accepts any value: use 1.0 so the
        // lattice is a full disc (the analytic solution needs one).
        slurry.fill_fraction = 1.0;
        slurry.viscosity_pa_s = slurry.density_kg_m3 * r * r / TAU_S as f32;
        let iters = params.simulation.pbf_iterations;
        let omega0 = 2.0f32;
        let drum = make_drum(r, omega0);
        let mut f = FluidParticles::seed_lattice(&slurry, r, res, &[], 0.0);
        let m = f.particle_mass as f64;
        let l_inf = omega0 as f64 * m * f.x.iter().map(|p| p.length_squared() as f64).sum::<f64>();
        let dt = SPINUP_DT;
        let n_steps = (0.2 * TAU_S / dt as f64).round() as usize;
        let mut l = vec![0.0f64; n_steps + 1];
        let mut cum_wall = vec![0.0f64; n_steps + 1];
        let mut cum_grav = vec![0.0f64; n_steps + 1];
        let mut angle = 0.0f32;
        for k in 1..=n_steps {
            angle += omega0 * dt;
            let st = f.step(&drum, angle, &slurry, iters, dt);
            cum_wall[k] = cum_wall[k - 1] + st.wall_angular_impulse as f64;
            let sum_x: f64 = f.x.iter().map(|p| p.x as f64).sum();
            cum_grav[k] = cum_grav[k - 1] - m * G * sum_x * dt as f64;
            l[k] =
                f.x.iter()
                    .zip(&f.v)
                    .map(|(x, v)| (x.x * v.y - x.y * v.x) as f64)
                    .sum::<f64>()
                    * m;
        }
        let step_of = |frac: f64| ((frac * TAU_S / dt as f64).round() as usize).clamp(1, n_steps);
        let samples = gates::SPINUP_FRACTIONS
            .iter()
            .map(|&fr| {
                let k = step_of(fr);
                let t_nu = k as f64 * dt as f64 / TAU_S;
                (
                    t_nu,
                    l[k] / l_inf,
                    analytic::spinup_angular_momentum_fraction(t_nu),
                )
            })
            .collect();
        let (k1, k2) = (step_of(0.02), n_steps);
        let dl = l[k2] - l[k1];
        let jw = cum_wall[k2] - cum_wall[k1];
        let jg = cum_grav[k2] - cum_grav[k1];
        SpinUp {
            res,
            n: f.len(),
            l_inf,
            samples,
            torque_mismatch: (dl - jw - jg).abs() / dl.abs(),
            torque_mismatch_no_gravity: (dl - jw).abs() / dl.abs(),
        }
    }

    // ---------------------------------------------------------------- 3. Taylor-Couette

    #[derive(Debug, Clone)]
    pub struct TaylorCouette {
        pub res: u32,
        pub n: usize,
        pub t_analytic: f64,
        /// Mean signed torque (per metre) delivered to the ball (expected `-T`).
        pub torque_ball: f64,
        /// Mean signed torque the wall delivers to the fluid (expected `-T`).
        pub torque_wall: f64,
        pub clamp_hits: u64,
        pub steps: usize,
    }

    pub fn taylor_couette(res: u32) -> TaylorCouette {
        let params = Params::default();
        let r2 = params.mill.radius_m();
        let r1 = 0.5 * r2;
        let mut slurry = params.slurry;
        slurry.fill_fraction = 1.0; // full annulus (ball sites are skipped by seed_lattice)
        slurry.viscosity_pa_s = slurry.density_kg_m3 * r2 * r2 / TAU_S as f32;
        let iters = params.simulation.pbf_iterations;
        let omega1 = 1.0f32;
        let drum = make_drum(r2, 0.0);
        let eff = EffectiveMedia {
            true_diameter_m: 2.0 * r1,
            diameter_m: 2.0 * r1,
            density_kg_m3: 7800.0,
            ball_count: 1,
            scale_factor: 1.0,
        };
        let mut dem = DemState::new(&eff, r2, 1);
        assert_eq!(dem.balls.len(), 1, "expected exactly one ball");
        let mut f = FluidParticles::seed_lattice(&slurry, r2, res, &[Vec2::ZERO], r1);
        let (mut tb, mut tw, mut clamp, mut cnt) = (0.0f64, 0.0f64, 0u64, 0usize);
        let settle = (3.0 * TAU_S / DT as f64).round() as usize;
        let measure = (1.0 * TAU_S / DT as f64).round() as usize;
        for k in 0..settle + measure {
            dem.balls.x[0] = Vec2::ZERO;
            dem.balls.v[0] = Vec2::ZERO;
            dem.balls.omega[0] = omega1;
            dem.balls.theta[0] += omega1 * DT;
            let (imp, st) = f.step_coupled(&drum, 0.0, &slurry, iters, DT, &dem.balls);
            if k >= settle {
                tb += imp.angular_impulses[0] as f64;
                tw += st.wall_angular_impulse as f64;
                clamp += imp.clamp_hits as u64;
                cnt += 1;
            }
        }
        let t_win = cnt as f64 * DT as f64;
        TaylorCouette {
            res,
            n: f.len(),
            t_analytic: analytic::taylor_couette_torque(
                slurry.viscosity_pa_s as f64,
                omega1 as f64,
                r1 as f64,
                r2 as f64,
            ),
            torque_ball: tb / t_win,
            torque_wall: tw / t_win,
            clamp_hits: clamp,
            steps: cnt,
        }
    }

    // ---------------------------------------------------------------- 4/5. whole simulation

    fn percent_critical(params: &mut Params, percent: f32) {
        params.mill.speed_mode = SpeedMode::PercentCritical;
        params.mill.speed_value = percent;
    }

    fn advance(sim: &mut Simulation, seconds: f64) {
        let n = (seconds * 60.0).round() as usize;
        for _ in 0..n {
            sim.step(1.0 / 60.0);
        }
    }

    /// Default params, wet, 60 % critical, given fluid resolution: settle 8 s, then 2 revolutions.
    pub fn energy_closure(res: u32) -> EnergyBudget {
        let mut p = Params::default();
        p.simulation.resolution = res;
        percent_critical(&mut p, 60.0);
        let rev_s = 60.0 / p.mill.rpm() as f64;
        let mut sim = Simulation::new(p).expect("valid params");
        advance(&mut sim, 8.0);
        sim.reset_energy_budget();
        advance(&mut sim, 2.0 * rev_s);
        sim.energy_budget()
    }

    /// Dry mean shaft power (W per metre) over a 4 s window after a 6 s settle, and the
    /// `effective_substeps` actually used. `simulation.substeps` is validated to [1, 16].
    pub fn dry_mean_power(substeps: u32, seed: u64) -> (f64, u32) {
        let mut p = Params::default();
        p.slurry.enabled = false;
        p.simulation.substeps = substeps;
        p.simulation.seed = seed;
        // Large media so `effective_substeps()` (which auto-raises the rate for small, coarse-
        // grained balls -- at the default 2 mm media it is 16 for every request) equals the
        // requested `substeps`; otherwise the comparison would be vacuous.
        p.media.ball_diameter_m = 0.01;
        p.media.fill_fraction = 0.5;
        percent_critical(&mut p, 60.0);
        let eff = p.effective_substeps();
        assert_eq!(eff, substeps, "effective_substeps must equal the request");
        let mut sim = Simulation::new(p).expect("valid params");
        advance(&mut sim, 6.0);
        sim.reset_energy_budget();
        advance(&mut sim, 4.0);
        let b = sim.energy_budget();
        (b.shaft_work_j / b.elapsed_s, eff)
    }
}

/// Flat rows for the probe.
pub mod report {
    use super::{gates, harness};

    pub struct Row {
        pub case: &'static str,
        pub res: u32,
        pub metric: String,
        pub value: f64,
        pub gate: String,
    }

    fn row(case: &'static str, res: u32, metric: &str, value: f64, gate: String) -> Row {
        Row {
            case,
            res,
            metric: metric.to_string(),
            value,
            gate,
        }
    }

    pub fn hydrostatic(res: u32) -> Vec<Row> {
        let h = harness::hydrostatic(res);
        vec![
            row("hydrostatic", res, "n_particles", h.n as f64, "-".into()),
            row(
                "hydrostatic",
                res,
                "worst5_compression",
                h.worst5_c,
                format!("<= {}", gates::HYDRO_WORST5_COMPRESSION),
            ),
            row(
                "hydrostatic",
                res,
                "mean_compression",
                h.mean_c,
                format!("<= {}", gates::HYDRO_MEAN_COMPRESSION),
            ),
            row(
                "hydrostatic",
                res,
                "ke_over_mgdx",
                h.ke_over_mgdx,
                format!("<= {}", gates::HYDRO_KE_OVER_MGDX),
            ),
            row(
                "hydrostatic",
                res,
                "wall_pressure_l2_rel",
                h.wall_l2_rel,
                format!("<= {}", gates::HYDRO_WALL_PRESSURE_L2),
            ),
            row(
                "hydrostatic",
                res,
                "wall_bins_used",
                h.bins_used as f64,
                "-".into(),
            ),
            row("hydrostatic", res, "y_surface_m", h.y_surface, "-".into()),
        ]
    }

    pub fn spin_up(res: u32) -> Vec<Row> {
        let s = harness::spin_up(res);
        let mut rows = vec![row("spin_up", res, "n_particles", s.n as f64, "-".into())];
        for (t, ls, le) in &s.samples {
            rows.push(row(
                "spin_up",
                res,
                &format!("L_sim/Linf@t/tau={t:.3}"),
                *ls,
                "-".into(),
            ));
            rows.push(row(
                "spin_up",
                res,
                &format!("L_exact/Linf@t/tau={t:.3}"),
                *le,
                "-".into(),
            ));
            rows.push(row(
                "spin_up",
                res,
                &format!("|dL|/Linf@t/tau={t:.3}"),
                (ls - le).abs(),
                format!("<= {}", gates::SPINUP_L_ERR),
            ));
        }
        rows.push(row(
            "spin_up",
            res,
            "torque_mismatch",
            s.torque_mismatch,
            format!("<= {}", gates::SPINUP_TORQUE_MISMATCH),
        ));
        rows.push(row(
            "spin_up",
            res,
            "torque_mismatch_no_gravity",
            s.torque_mismatch_no_gravity,
            "-".into(),
        ));
        rows
    }

    pub fn taylor_couette(res: u32) -> Vec<Row> {
        let t = harness::taylor_couette(res);
        let g = format!("<= {}", gates::TC_TORQUE_REL);
        vec![
            row("taylor_couette", res, "n_particles", t.n as f64, "-".into()),
            row(
                "taylor_couette",
                res,
                "T_analytic",
                t.t_analytic,
                "-".into(),
            ),
            row(
                "taylor_couette",
                res,
                "torque_ball_signed",
                t.torque_ball,
                "-".into(),
            ),
            row(
                "taylor_couette",
                res,
                "ball_rel_err",
                (t.torque_ball.abs() - t.t_analytic).abs() / t.t_analytic,
                g.clone(),
            ),
            row(
                "taylor_couette",
                res,
                "torque_wall_signed",
                t.torque_wall,
                "-".into(),
            ),
            row(
                "taylor_couette",
                res,
                "wall_rel_err",
                (t.torque_wall.abs() - t.t_analytic).abs() / t.t_analytic,
                g,
            ),
            row(
                "taylor_couette",
                res,
                "clamp_hits",
                t.clamp_hits as f64,
                "-".into(),
            ),
        ]
    }

    pub fn energy(res: u32) -> Vec<Row> {
        let b = harness::energy_closure(res);
        let w = |j: f64| j / b.elapsed_s;
        let pct = |j: f64| 100.0 * j / b.shaft_work_j;
        let mut rows = Vec::new();
        let terms: [(&str, f64); 10] = [
            ("shaft_work", b.shaft_work_j),
            ("dem_wall_work", b.dem_wall_work_j),
            ("fluid_wall_work", b.fluid_wall_work_j),
            ("ball_contact_dissipation", b.ball_contact_dissipation_j),
            ("fluid_wall_slip", b.fluid_wall_slip_j),
            ("fluid_viscous", b.fluid_viscous_j),
            ("fluid_clamp_removed", b.fluid_clamp_removed_j),
            ("interface_created", b.interface_created_j),
            ("delta_mechanical", b.delta_mechanical_j),
            ("unattributed", b.unattributed_j),
        ];
        rows.push(row("energy", res, "elapsed_s", b.elapsed_s, "-".into()));
        for (name, j) in terms {
            rows.push(row("energy", res, &format!("{name}_W"), w(j), "-".into()));
            rows.push(row(
                "energy",
                res,
                &format!("{name}_pct_shaft"),
                pct(j),
                "-".into(),
            ));
        }
        rows.push(row(
            "energy",
            res,
            "|unattributed|/shaft",
            (b.unattributed_j / b.shaft_work_j).abs(),
            format!("<= {}", gates::ENERGY_UNATTRIBUTED_FRAC),
        ));
        rows
    }

    pub fn dry() -> Vec<Row> {
        let seeds = [1u64, 2];
        let mean = |s: u32| {
            let v: Vec<(f64, u32)> = seeds
                .iter()
                .map(|&seed| harness::dry_mean_power(s, seed))
                .collect();
            (v.iter().map(|x| x.0).sum::<f64>() / v.len() as f64, v[0].1)
        };
        let (p16, e16) = mean(16);
        let mut rows = vec![row(
            "dry_substeps",
            0,
            "power_W_substeps16",
            p16,
            format!("(effective {e16})"),
        )];
        for s in [4u32, 8] {
            let (p, e) = mean(s);
            rows.push(row(
                "dry_substeps",
                0,
                &format!("power_W_substeps{s}"),
                p,
                format!("(effective {e})"),
            ));
            rows.push(row(
                "dry_substeps",
                0,
                &format!("rel_diff_vs_16_substeps{s}"),
                (p - p16).abs() / p16,
                format!("<= {}", gates::DRY_SUBSTEP_REL),
            ));
        }
        rows
    }
}
