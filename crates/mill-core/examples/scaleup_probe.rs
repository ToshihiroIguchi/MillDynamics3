//! Validation probe for the scale-up runner (`mill_core::scaleup`).
//!
//! `cargo run --release -p mill-core --example scaleup_probe -- [--cases conv,sim,sweep,sens,wet]
//!  [--seeds 3]`
//!
//! Prints markdown tables. Cases:
//! - `conv`: sensitivity of the stress-energy summary to the sub-step (`dt_ratio`) and the contact
//!   iterations, relative to the finest setting.
//! - `sim`: dry similarity check. The same `D/d`, same `% critical`, ball size and drum scaled
//!   together. A rigid-contact model should give the same dimensionless numbers; the residual is the
//!   absolute-speed constants left in the solver (restitution threshold, friction speed scale).
//! - `sweep`: `D/d` sweep with power-law fits against the drop (+1) and shear (-1) expectations.
//! - `sens`: +-20 % contact coefficients.
//! - `wet`: lubrication-cushion viscosity sweep at fixed geometry.

use mill_core::params::Params;
use mill_core::scaleup::{self, ScaleupConfig};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1).cloned())
}

fn base() -> Params {
    let mut p = Params::default();
    p.slurry.enabled = false;
    p
}

fn show(names: &[&str], mean: &[f64], err: &[f64], all: &[String]) -> String {
    names
        .iter()
        .map(|n| {
            let k = all.iter().position(|a| a == n).unwrap();
            format!("{:.4e} ± {:.1e}", mean[k], err[k])
        })
        .collect::<Vec<_>>()
        .join(" | ")
}

const KEY: [&str; 5] = [
    "specific_power_w_per_kg",
    "impacts_total_per_ball_per_s",
    "mean_se_scaled_drop",
    "mean_se_scaled_shear",
    "stress_intensity_scaled",
];

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cases = arg(&args, "--cases").unwrap_or_else(|| "conv,sim,sweep,sens,wet".into());
    let n_seeds: u64 = arg(&args, "--seeds")
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let seeds: Vec<u64> = (1..=n_seeds).collect();
    let cfg = ScaleupConfig {
        seeds: seeds.clone(),
        ..Default::default()
    };
    println!("seeds {seeds:?}; key = {KEY:?}\n");

    if cases.contains("conv") {
        println!("## conv: sub-step and iterations (63 mm drum, 2 mm ZrO2, dry)\n");
        println!("| dt_ratio | iters | {} |", KEY.join(" | "));
        println!("|---|---|---|---|---|---|---|");
        for (dt, it) in [
            (0.1, 6),
            (0.05, 6),
            (0.025, 6),
            (0.05, 12),
            (0.025, 12),
            (0.0125, 12),
            (0.025, 20),
        ] {
            let c = ScaleupConfig {
                dt_ratio: dt,
                dem_iterations: it,
                ..cfg.clone()
            };
            let e = scaleup::run_ensemble(&base(), &c).unwrap();
            println!(
                "| {dt} | {it} | {} |",
                show(&KEY, &e.mean, &e.std_err, &e.scalar_names)
            );
        }
        println!();
    }

    if cases.contains("sim") {
        println!("## sim: dry similarity at D/d = 31.5, 59 %Nc (drum and ball scaled together)\n");
        println!("| D (m) | d (m) | balls | {} |", KEY.join(" | "));
        println!("|---|---|---|---|---|---|---|---|");
        for scale in [1.0f32, 2.0, 4.0] {
            let mut p = base();
            let pc = p.mill.percent_critical();
            p.mill.speed_mode = mill_core::params::SpeedMode::PercentCritical;
            p.mill.speed_value = pc;
            p.mill.diameter_m *= scale;
            p.media.ball_diameter_m *= scale;
            let e = scaleup::run_ensemble(&p, &cfg).unwrap();
            println!(
                "| {} | {} | {} | {} |",
                p.mill.diameter_m,
                p.media.ball_diameter_m,
                e.first.summary.n_balls,
                show(&KEY, &e.mean, &e.std_err, &e.scalar_names)
            );
        }
        println!();
    }

    if cases.contains("sweep") {
        println!("## sweep: D/d (ball 2 mm fixed, 59 %Nc)\n");
        let ratios = [20.0, 30.0, 45.0, 60.0, 90.0];
        let s = scaleup::sweep_dd(&base(), &cfg, &ratios).unwrap();
        println!("| D/d | balls | {} |", KEY.join(" | "));
        println!("|---|---|---|---|---|---|---|");
        for p in &s.points {
            println!(
                "| {} | {} | {} |",
                p.dd_ratio,
                p.n_balls,
                show(
                    &KEY,
                    &p.ensemble.mean,
                    &p.ensemble.std_err,
                    &p.ensemble.scalar_names
                )
            );
        }
        println!("\n| quantity | fitted exponent in D/d | expected | R^2 |\n|---|---|---|---|");
        for f in &s.fits {
            println!(
                "| {} | {:.3} | {} | {:.3} |",
                f.name,
                f.exponent,
                f.expected_exponent
                    .map(|e| format!("{e:+.1}"))
                    .unwrap_or("-".into()),
                f.r2
            );
        }
        println!();
    }

    if cases.contains("sens") {
        println!("## sens: +-20 % contact coefficients\n");
        println!("| case | {} |", KEY.join(" | "));
        println!("|---|---|---|---|---|---|");
        let cases = scaleup::sensitivity(&base(), &cfg, 0.2).unwrap();
        let names = cases[0].scalar_names.clone();
        let zero = vec![0.0; names.len()];
        for c in &cases {
            println!("| {} | {} |", c.label, show(&KEY, &c.mean, &zero, &names));
        }
        println!();
    }

    if cases.contains("wet") {
        println!("## wet: lubrication-cushion viscosity (Pa s)\n");
        println!("| mu | St_ref | {} |", KEY.join(" | "));
        println!("|---|---|---|---|---|---|---|");
        for mu in [0.0f32, 0.001, 0.01, 0.1, 1.0] {
            let c = ScaleupConfig {
                wet_viscosity_pa_s: mu,
                ..cfg.clone()
            };
            let e = scaleup::run_ensemble(&base(), &c).unwrap();
            println!(
                "| {mu} | {} | {} |",
                e.first
                    .summary
                    .reference_stokes_number
                    .map(|s| format!("{s:.3e}"))
                    .unwrap_or("dry".into()),
                show(&KEY, &e.mean, &e.std_err, &e.scalar_names)
            );
        }
    }
}
