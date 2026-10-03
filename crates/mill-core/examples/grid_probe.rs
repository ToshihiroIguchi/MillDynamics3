//! Experiment tables for the fixed-grid solver track. `cargo run --release -p mill-core --example
//! grid_probe -- --exp E0`.

use mill_core::mac::verify::{
    diagnose_couette_mac, eccentric_squeeze_force, stokes_annulus_drag, track_sheet_growth,
    verify_buoyancy, verify_couette_mac_scheme, verify_couette_ns, verify_couette_stokes,
    verify_dam_break, verify_drum_slurry, verify_levelset_rotation, verify_manufactured,
    verify_moving_disc, verify_rigid_ring, verify_rimming, verify_sloshing_circle,
    verify_sloshing_rect, verify_spinup, verify_spinup_mac, verify_still_pool, verify_translation,
    verify_wall_impact,
};
use std::time::Instant;

fn main() {
    let exp = std::env::args()
        .skip_while(|a| a != "--exp")
        .nth(1)
        .unwrap_or_else(|| "E0".into());
    match exp.as_str() {
        "E0" => e0(),
        "E1a" => e1a(),
        "E1s" => e1s(),
        "E1b" => e1b(),
        "R0" => r0(),
        "R2" => r2(),
        "R2d" => r2d(),
        "R2s" => r2s(),
        "E2ls" => e2ls(),
        "E2a" => e2a(),
        "E2ref" => e2ref(),
        "E2b" => e2b(),
        "E2bc" => e2bc(),
        "E2c" => e2c(),
        "E2w" => e2w(),
        "E2t" => e2t(),
        "E2g" => e2g(),
        "E2d" => e2d(),
        "E2r" => e2r(),
        "E3s" => e3s(),
        "E4a" => e4a(),
        "E4b" => e4b(),
        "E4s" => e4s(),
        "E4m" => e4m(),
        "E4d" => e4d(),
        "E4e" => e4e(),
        "E4k" => e4k(),
        "E4h" => e4h(),
        "E4p" => e4p(),
        "E4y" => e4y(),
        "E4r" => e4rot(),
        "E5a" => e5a(),
        "E5b" => e5b(),
        "E5t" => e5t(),
        "E5c" => e5c(),
        "E5d" => e5d(),
        "E5e" => e5e(),
        "E5f" => e5f(),
        "E5g" => e5g(),
        other => eprintln!("unknown experiment {other}"),
    }
}

fn e0() {
    println!("E0 manufactured Neumann Poisson on the unit disc (cut cells)");
    println!(
        "{:>6} {:>12} {:>7} {:>6} {:>10} {:>9}",
        "n", "L2 error", "order", "iters", "residual", "ms"
    );
    let mut prev: Option<f64> = None;
    for n in [32usize, 64, 128, 256, 512] {
        let t = Instant::now();
        let (e, st, _) = verify_manufactured(n);
        let order = prev.map_or(f64::NAN, |p| (p / e).log2());
        println!(
            "{n:>6} {e:>12.3e} {order:>7.2} {:>6} {:>10.1e} {:>9.1}",
            st.iterations,
            st.residual,
            t.elapsed().as_secs_f64() * 1e3
        );
        prev = Some(e);
    }
}

fn e1a() {
    println!("E1a Stokes Taylor-Couette (R1=0.25, R2=0.5), componentwise Dirichlet solve");
    println!(
        "{:>6} {:>12} {:>14} {:>14}",
        "n", "max u err", "inner torque", "outer torque"
    );
    for n in [32usize, 64, 128, 256, 512] {
        let (e, ti, to) = verify_couette_stokes(n);
        println!(
            "{n:>6} {e:>12.3e} {:>13.3}% {:>13.3}%",
            100.0 * ti,
            100.0 * to
        );
    }
}

fn e1s() {
    println!("E1 spin-up of a disc (impulsive wall), BDF2; errors of L(t) and wall torque");
    let cps = [0.02, 0.05, 0.1, 0.3, 1.0];
    for (n, dt) in [
        (64usize, 0.01),
        (128, 0.01),
        (128, 0.005),
        (256, 0.005),
        (256, 0.0025),
    ] {
        let pts = verify_spinup(n, dt, &cps);
        print!("n={n:<4} dt={dt:<7}");
        for p in &pts {
            print!(
                " | T={:<4} L{:>7.3}% Q{:>7.3}%",
                p.t_nu,
                100.0 * p.l_err,
                100.0 * p.torque_err
            );
        }
        println!();
    }
}

fn e1b() {
    println!("E1b Navier-Stokes Taylor-Couette hold test from the exact state (t_end = 10)");
    println!(
        "{:>5} {:>7} {:>6} {:>10} {:>12} {:>12} {:>10}",
        "n", "Re", "steps", "u err", "inner T", "outer T", "div"
    );
    for re in [10.0f64, 1000.0] {
        let nu = 0.25 * 0.25 / re;
        for n in [64usize, 128, 256] {
            let r = verify_couette_ns(n, nu, 10.0, 0.25);
            println!(
                "{n:>5} {re:>7} {:>6} {:>10.2e} {:>11.3}% {:>11.3}% {:>10.1e}",
                r.steps,
                r.u_err,
                100.0 * r.inner_torque_err,
                100.0 * r.outer_torque_err,
                r.divergence
            );
        }
    }
}

fn r0() {
    println!("R0.1 hold test, time-step dependence (t_end = 10)");
    println!(
        "{:>5} {:>7} {:>6} {:>7} {:>10} {:>12} {:>12}",
        "n", "Re", "cfl", "steps", "u err", "inner T", "outer T"
    );
    for re in [10.0f64, 1000.0] {
        let nu = 0.25 * 0.25 / re;
        for n in [64usize, 128] {
            for cfl in [0.25f64, 0.0625] {
                let r = verify_couette_ns(n, nu, 10.0, cfl);
                println!(
                    "{n:>5} {re:>7} {cfl:>6} {:>7} {:>10.2e} {:>11.3}% {:>11.3}%",
                    r.steps,
                    r.u_err,
                    100.0 * r.inner_torque_err,
                    100.0 * r.outer_torque_err
                );
            }
        }
    }
}

fn r2() {
    println!("R2 staggered hold test from the exact state (t_end = 10, cfl 0.25)");
    println!(
        "{:>5} {:>7} {:>6} {:>10} {:>12} {:>12} {:>10} {:>9}",
        "n", "Re", "steps", "u err", "inner T", "outer T", "div", "ms/step"
    );
    let ns: Vec<usize> = std::env::args()
        .skip_while(|a| a != "--n")
        .skip(1)
        .map(|a| a.parse().unwrap())
        .collect();
    let upwind = std::env::args().any(|a| a == "--upwind");
    println!("advection: {}", if upwind { "upwind3" } else { "central" });
    let ns = if ns.is_empty() {
        vec![64, 128, 256]
    } else {
        ns
    };
    for re in [10.0f64, 1000.0] {
        let nu = 0.25 * 0.25 / re;
        for &n in &ns {
            let t = Instant::now();
            let r = verify_couette_mac_scheme(n, nu, 10.0, 0.25, upwind);
            println!(
                "{n:>5} {re:>7} {:>6} {:>10.2e} {:>11.3}% {:>11.3}% {:>10.1e} {:>9.2}",
                r.steps,
                r.u_err,
                100.0 * r.inner_torque_err,
                100.0 * r.outer_torque_err,
                r.divergence,
                t.elapsed().as_secs_f64() * 1e3 / r.steps as f64
            );
        }
    }
}

fn r2d() {
    println!("R2d error history (fluid-node / ghost-face max error, units of omega r1)");
    for re in [10.0f64] {
        let nu = 0.25 * 0.25 / re;
        for n in [64usize, 128] {
            for (s, a, b) in diagnose_couette_mac(n, nu, 0.25, &[0, 1, 2, 5, 20, 100, 400]) {
                println!("n={n:<4} Re={re:<5} step {s:<4} fluid {a:.3e} ghost {b:.3e}");
            }
        }
    }
}

fn r2s() {
    println!("R2s spin-up of a disc, staggered NS stepper; errors of L(t) and wall torque");
    let cps = [0.02, 0.05, 0.1, 0.3, 1.0];
    for (n, dt) in [(64usize, 0.01), (128, 0.01), (128, 0.005), (256, 0.005)] {
        let pts = verify_spinup_mac(n, dt, &cps);
        print!("n={n:<4} dt={dt:<7}");
        for p in &pts {
            print!(
                " | T={:<4} L{:>7.3}% Q{:>7.3}%",
                p.t_nu,
                100.0 * p.l_err,
                100.0 * p.torque_err
            );
        }
        println!();
    }
}

fn e2ls() {
    println!("E2 level set: rigid rotation, 1 revolution (WENO5 + RK3, reinit every 2 steps)");
    println!(
        "{:>9} {:>5} {:>6} {:>12} {:>14}",
        "shape", "n", "steps", "L1 shape", "volume drift"
    );
    for (slotted, every) in [(false, 0usize), (false, 2), (true, 0), (true, 2)] {
        for n in [64usize, 128, 256] {
            let r = verify_levelset_rotation(n, slotted, 1.0, every);
            println!(
                "{:>9} {n:>5} {:>6} {:>12.3e} {:>13.4}%  reinit_every={every}",
                if slotted { "zalesak" } else { "disc" },
                r.steps,
                r.shape_l1,
                100.0 * r.volume_drift
            );
        }
    }
}

fn e2a() {
    println!("E2a still pool in a drum (R = 0.5), t_end = 1 s, nu = 1e-6");
    println!(
        "{:>5} {:>7} {:>6} {:>12} {:>12} {:>12} {:>10} {:>9}",
        "n", "level", "steps", "|u|/sqrt(gD)", "p err", "volume", "y err/dx", "ms/step"
    );
    for level in [0.0f64, -0.13] {
        for n in [64usize, 128] {
            let t = Instant::now();
            let r = verify_still_pool(n, level, 1e-6, 1.0);
            println!(
                "{n:>5} {level:>7} {:>6} {:>12.2e} {:>11.3}% {:>11.4}% {:>10.3} {:>9.2}",
                r.steps,
                r.spurious,
                100.0 * r.pressure_err,
                100.0 * r.volume_drift,
                r.level_err_cells,
                t.elapsed().as_secs_f64() * 1e3 / r.steps as f64
            );
        }
    }
}

fn e2ref() {
    println!("Half-full circular container: K R = omega^2 R / g, Rayleigh-Ritz (harmonic basis)");
    for terms in [8usize, 12, 16, 20] {
        let k = mill_core::mac::reference::half_disc_sloshing(terms, 4);
        println!("terms={terms:<3} {k:.6?}");
    }
}

fn e2b() {
    use mill_core::mac::reference::rectangular_sloshing_omega;
    println!("E2b small-amplitude sloshing, rectangular tank 1.0 x depth 0.5, amp 2 mm, nu = 1e-6");
    let omega_ref = rectangular_sloshing_omega(9.81, 1.0, 0.5, 1);
    println!(
        "reference omega = {omega_ref:.5} rad/s (T = {:.4} s)",
        2.0 * std::f64::consts::PI / omega_ref
    );
    println!(
        "{:>5} {:>8} {:>6} {:>10} {:>10} {:>10} {:>10}",
        "n", "scheme", "steps", "omega", "err", "damping", "volume"
    );
    let upwind_only = std::env::args().any(|a| a == "--upwind");
    for (n, nu) in [
        (64usize, 1e-6),
        (128, 1e-6),
        (128, 1e-4),
        (128, 1e-3),
        (128, 1e-2),
    ] {
        for upwind in [false, true] {
            if upwind_only && !upwind {
                continue;
            }
            let r = verify_sloshing_rect(n, 1.0, 0.5, 0.002, nu, 3.0, upwind);
            print!("nu={nu:<6}");
            println!(
                "{n:>5} {:>8} {:>6} {:>10.5} {:>9.3}% {:>10.4} {:>9.4}%",
                if upwind { "upwind3" } else { "central" },
                r.steps,
                r.omega,
                100.0 * (r.omega / omega_ref - 1.0),
                r.damping,
                100.0 * r.volume_drift
            );
        }
    }
}

fn e2bc() {
    let k = mill_core::mac::reference::half_disc_sloshing(16, 1)[0];
    let omega_ref = (k * 9.81 / 0.5).sqrt();
    println!(
        "E2b half-full circular drum R = 0.5, tilt amp 25 mm, Ritz K R = {k:.6}, omega = {omega_ref:.5} rad/s"
    );
    println!(
        "{:>5} {:>8} {:>8} {:>6} {:>10} {:>10} {:>10}",
        "n", "nu", "steps", "scheme", "omega", "err", "damping"
    );
    for (n, nu) in [(64usize, 1e-6), (128, 1e-6), (256, 1e-6), (128, 1e-3)] {
        let r = verify_sloshing_circle(n, 0.025, nu, 3.0, true);
        println!(
            "{n:>5} {nu:>8} {:>8} {:>6} {:>10.5} {:>9.3}% {:>10.4}",
            r.steps,
            "upwind3",
            r.omega,
            100.0 * (r.omega / omega_ref - 1.0),
            r.damping
        );
    }
}

fn e2c() {
    let ts = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
    println!("E2c dam break of a square column a = 0.2 in a 1.0 tank, nu = 1e-6 (Z = x_front / a)");
    print!("{:>5} {:>7} {:>8}", "n", "scheme", "volume");
    for t in ts {
        print!("  T={t:<4}");
    }
    println!("  max|u|");
    for (n, upwind) in [(64usize, true), (128, true), (256, true), (128, false)] {
        let r = verify_dam_break(n, 0.2, 1e-6, &ts, upwind);
        print!(
            "{n:>5} {:>7} {:>7.4}%",
            if upwind { "upwind3" } else { "central" },
            100.0 * r.volume_drift
        );
        for (_, z) in &r.front {
            print!("  {z:>6.3}");
        }
        println!("  {:.2}", r.max_speed);
    }
    println!("Ritter shallow-water upper bound: Z = 1 + 2 T");
}

fn e2w() {
    println!("E2w liquid layer (u0 = 2.5 m/s, gravity) hitting a wall at x = 0.5 (t_end = 0.8 s): blow-up time or final speed");
    println!("{:>5} {:>8} {:>12} {:>10}", "n", "h/dx", "blow-up", "speed");
    for n in [64usize, 128, 256] {
        for h in [0.6, 1.0, 1.7, 3.0, 6.0] {
            let r = verify_wall_impact(n, h, 2.5, 1e-6, 0.8, 0.5, &|f| f.upwind = true);
            println!(
                "{n:>5} {h:>8} {:>12} {:>10.2}",
                r.blew_up_at.map_or("-".to_string(), |t| format!("{t:.4}")),
                r.final_speed
            );
        }
    }
}

fn e2t() {
    println!("E2t rigid translation of a liquid block (side 0.2, u0 = 2.5, no gravity), t = 0.1 s");
    for n in [64usize, 128] {
        match verify_translation(n, 0.2, 2.5, 1e-6, 0.1, &|_| {}) {
            Some((dev, vmax)) => println!("n={n:<4} max|u-u0| = {dev:.3e}  max|v| = {vmax:.3e}"),
            None => println!("n={n:<4} BLEW UP"),
        }
    }
}

fn e2g() {
    println!("E2g thin sheet sliding on a floor, u0 = 2.5, no gravity: max|v| at t = 0.1");
    for n in [64usize, 128] {
        for h in [0.6, 1.0, 1.7, 3.0] {
            let hist = track_sheet_growth(n, h, 2.5, 0.1, true, &|f| f.upwind = true);
            let last = hist.last().unwrap();
            println!(
                "n={n:<4} h/dx={h:<4} t={:.4} max|v|={:.3e} at {:?}",
                last.0, last.1, last.3
            );
        }
    }
}

fn e2d() {
    println!(
        "E2d rimming flow (R = 0.5, Omega = 1, nu = 0.02, ring h0 = 16 mm) vs Moffatt thin film"
    );
    let revs: f64 = std::env::args()
        .skip_while(|a| a != "--revs")
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(2.0);
    for n in [128usize, 256] {
        let r = verify_rimming(n, 0.02, 1.0, 0.016, revs);
        println!(
            "n={n:<4} steps={} q/qmax={:.3} h/R(max)={:.3} max rel err={:.2}% rms={:.2}% volume={:.4}%",
            r.steps,
            r.q_over_qmax,
            r.h_over_r,
            100.0 * r.max_rel_err,
            100.0 * r.rms_rel_err,
            100.0 * r.volume_drift
        );
        for k in (0..r.phi.len()).step_by(6) {
            println!(
                "   phi={:6.1} deg  h_sim={:.4}  h_theory={:.4}",
                r.phi[k].to_degrees(),
                r.h_sim[k],
                r.h_theory[k]
            );
        }
    }
}

fn e2r() {
    use mill_core::mac::staggered::StaggeredFlow;
    println!("E2r rigidly rotating liquid ring (exact steady free surface), R = 0.5, r_i = 0.3, nu = 1e-3, t = 1 s");
    println!(
        "{:>16} {:>9} {:>4} {:>11} {:>10} {:>10}",
        "variant", "Omega,g", "n", "surface/dx", "velocity", "pressure"
    );
    type Cfg = Box<dyn Fn(&mut StaggeredFlow)>;
    let variants: Vec<(&str, Cfg)> = vec![
        ("default (weno)", Box::new(|f| f.weno = true)),
        (
            "robust surface",
            Box::new(|f| {
                f.weno = true;
                f.robust_surface = true;
            }),
        ),
    ];
    for (name, cfg) in &variants {
        for (omega, g) in [(8.0f64, 0.0f64), (8.0, 9.81)] {
            for n in [64usize, 128, 256] {
                let r = verify_rigid_ring(n, omega, 0.3, 1e-3, 1.0, g, cfg.as_ref());
                println!(
                    "{name:>16} {omega:>4},{g:<4} {n:>4} {:>11.3} {:>10.2e} {:>10.2e}",
                    r.surface_err_cells, r.velocity_err, r.pressure_err
                );
            }
        }
    }
}

fn e3s() {
    println!("E3s viscous slurry in a rotating drum, R = 0.5, Fr = omega^2 R / g = 0.36, fill 0.3");
    println!(
        "{:>6} {:>5} {:>6} {:>11} {:>11} {:>8} {:>8} {:>9}  per-rev torque",
        "Re", "n", "steps", "torque", "gravity", "g-bal", "ripple", "volume"
    );
    let args: Vec<String> = std::env::args().collect();
    let arg = |name: &str, default: f64| -> f64 {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    };
    let settle = arg("--settle", 2.0);
    let fill = arg("--fill", 0.3);
    let measure = arg("--measure", 2.0);
    let omega = (0.36 * 9.81 / 0.5f64).sqrt();
    let res: Vec<f64> = match args
        .iter()
        .position(|a| a == "--re")
        .and_then(|i| args.get(i + 1))
    {
        Some(v) => v.split(',').filter_map(|x| x.parse().ok()).collect(),
        None => vec![0.2],
    };
    let ns: Vec<usize> = match args
        .iter()
        .position(|a| a == "--n")
        .and_then(|i| args.get(i + 1))
    {
        Some(v) => v.split(',').filter_map(|x| x.parse().ok()).collect(),
        None => vec![64, 128],
    };
    for re in res {
        let nu = omega * 0.25 / re;
        for &n in &ns {
            let t0 = Instant::now();
            let r = verify_drum_slurry(n, nu, omega, fill, settle, measure);
            println!(
                "{re:>6} {n:>5} {:>6} {:>11.6} {:>11.6} {:>7.2}% {:>7.2}% {:>8.4}%  {:.5?}  ({:.0?})",
                r.steps,
                r.torque,
                r.torque_gravity,
                100.0 * (r.torque / r.torque_gravity - 1.0),
                100.0 * r.torque_ripple,
                100.0 * r.volume_drift,
                r.per_rev,
                t0.elapsed()
            );
            println!(
                "       gravity from faces {:.6}  pressure torque {:.6}  wall-gravity {:.6}  dL/dt {:.6}",
                r.torque_gravity_faces,
                r.torque_pressure,
                r.torque - r.torque_gravity_faces,
                r.l_rate
            );
            let line: Vec<String> = r
                .film
                .iter()
                .step_by(6)
                .map(|h| format!("{h:.4}"))
                .collect();
            println!("       film(phi=2.5deg step 30deg) {}", line.join(" "));
        }
    }
}

fn e4a() {
    println!(
        "E4a fixed disc r = 0.1 in a still pool (level 0.2): pressure force vs displaced weight"
    );
    println!(
        "{:>5} {:>6} {:>6} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "n", "d/dx", "steps", "F_y plain", "F_y", "F_x", "torque", "spurious"
    );
    for n in [64usize, 128, 256] {
        let r = verify_buoyancy(n, 0.1, 0.2);
        println!(
            "{n:>5} {:>6.1} {:>6} {:>9.3}% {:>9.3}% {:>10.2e} {:>10.2e} {:>10.2e}",
            0.2 / (1.1 / n as f64),
            r.steps,
            100.0 * r.force_err_plain,
            100.0 * r.force_err,
            r.side_force,
            r.torque,
            r.spurious
        );
    }
}

fn e4b() {
    println!("E4b translating disc a = 0.2 in the drum (R = 0.5), nu = 1; fluid from rest");
    println!(
        "{:>5} {:>8} {:>8} {:>6} {:>11} {:>11} {:>9} {:>10} {:>10}",
        "n", "mode", "x_end", "steps", "drag", "reference", "err", "side", "diverg"
    );
    let (a, nu) = (0.2, 1.0);
    let u0 = 0.01;
    let reference = stokes_annulus_drag(a, 0.5, u0, nu);
    for n in [64usize, 128, 256] {
        let r = verify_moving_disc(n, a, 0.0, u0, 0.0, nu, 0.6, 0.01, false);
        println!(
            "{n:>5} {:>8} {:>8.4} {:>6} {:>11.6e} {:>11.6e} {:>8.3}% {:>10.2e} {:>10.2e}",
            "frozen",
            r.centre_x,
            r.steps,
            r.drag,
            reference,
            100.0 * (r.drag / reference - 1.0),
            r.side,
            r.divergence
        );
    }
    println!(
        "moving vs frozen (U = 0.1, window mean of the last 20 steps, frozen at the mean position)"
    );
    for n in [64usize, 128, 256] {
        for dt in [0.01, 0.005] {
            let t_end = 0.4;
            let mid = 0.1 * dt * 0.5 * 20.0;
            let frozen = verify_moving_disc(n, a, 0.04 - mid, 0.1, 0.0, nu, t_end, dt, false);
            let moved = verify_moving_disc(n, a, 0.0, 0.1, 0.0, nu, t_end, dt, true);
            println!(
                "{n:>5} dt={dt:<6} frozen {:>9.5} moving {:>9.5} diff {:>7.3}% noise {:>6.2}% (frozen {:>5.2}%)",
                frozen.drag_mean,
                moved.drag_mean,
                100.0 * (moved.drag_mean / frozen.drag_mean - 1.0),
                100.0 * moved.drag_noise,
                100.0 * frozen.drag_noise
            );
        }
    }
}

fn e4s() {
    println!("E4s frozen disc a = 0.2, U = 0.1, nu = 1: steady drag vs centre position");
    for n in [64usize, 128] {
        let row: Vec<String> = (0..=10)
            .map(|k| {
                let x0 = 0.004 * k as f64;
                let r = verify_moving_disc(n, 0.2, x0, 0.1, 0.0, 1.0, 0.4, 0.01, false);
                format!("{:.4}", r.drag_mean)
            })
            .collect();
        println!("n={n:<4} x0=0..0.04 step 0.004: {}", row.join(" "));
    }
}

fn e4rot() {
    println!("E4r disc a = 0.2 spinning at omega = 0.1 in the drum (R = 0.5), nu = 1: torque vs exact Couette");
    let (a, b, omega, nu) = (0.2, 0.5, 0.1, 1.0);
    let exact = 4.0 * std::f64::consts::PI * nu * omega * a * a * b * b / (b * b - a * a);
    for n in [64usize, 128, 256] {
        let r = verify_moving_disc(n, a, 0.0, 0.0, omega, nu, 0.6, 0.01, false);
        println!(
            "n={n:<4} torque {:.6} exact {:.6} err {:+.3}% (side force {:.1e})",
            r.torque,
            exact,
            100.0 * (r.torque / exact - 1.0),
            r.side
        );
    }
}

fn e4m() {
    use mill_core::mac::verify::{
        confined_added_mass, verify_added_mass_coupled, verify_added_mass_prescribed,
    };
    let a = 0.2;
    let only_terminal = std::env::args().any(|x| x == "--terminal");
    println!(
        "E4m added mass: disc a = {a} in the drum R = 0.5, m_a (confined, inviscid) = {:.5}",
        confined_added_mass(a, 0.5)
    );
    if !only_terminal {
        println!("prescribed acceleration a0 = 3 m/s^2 from rest, t = 0.05 s");
        println!(
            "{:>5} {:>8} {:>8} {:>7} {:>10} {:>10}",
            "n", "nu", "dt", "steps", "m_a err", "div"
        );
        for n in [64usize, 128] {
            for nu in [1e-2, 1e-3, 1e-4] {
                for dt in [0.005, 0.0025] {
                    let r = verify_added_mass_prescribed(n, a, 3.0, nu, 0.05, dt);
                    println!(
                        "{n:>5} {nu:>8.0e} {dt:>8} {:>7} {:>9.3}% {:>10.1e}",
                        r.steps,
                        100.0 * r.ratio_err,
                        r.divergence
                    );
                }
            }
        }
    }
    println!("terminal speed in Stokes flow (nu = 1, F = 1, t = 0.05 s) vs F / k, k = concentric Stokes drag per unit speed");
    for n in [64usize, 128] {
        for ratio in [0.5, 1.0, 3.0] {
            for dt in [0.005, 0.0025] {
                let (err, iters) =
                    mill_core::mac::verify::verify_terminal_stokes(n, a, ratio, 1.0, 1.0, 0.05, dt);
                println!(
                    "n={n:<4} rho_s={ratio:<4} dt={dt:<7} iters {iters:.1} speed err {:+.3}%",
                    100.0 * err
                );
            }
        }
    }
    if only_terminal {
        return;
    }
    println!("coupled free disc, constant force 1.0, t = 0.05 s: u(t) vs F t / (m + m_a)");
    println!(
        "{:>5} {:>6} {:>8} {:>8} {:>6} {:>10} {:>10}",
        "n", "rho_s", "nu", "dt", "iters", "u err", "div"
    );
    for n in [64usize, 128] {
        for ratio in [0.5, 1.0, 1.25, 3.0] {
            for dt in [0.005, 0.0025] {
                let r = verify_added_mass_coupled(n, a, ratio, 1.0, 1e-3, 0.05, dt);
                println!(
                    "{n:>5} {ratio:>6} {:>8.0e} {dt:>8} {:>6.1} {:>9.3}% {:>10.1e}",
                    1e-3,
                    r.iterations,
                    100.0 * r.speed_err,
                    r.divergence
                );
            }
        }
    }
}

fn e4k() {
    use mill_core::mac::bodies::BodyFlow;
    use mill_core::mac::staggered::Disc;
    println!(
        "E4k fluid force response dF/du of the first and later steps (nu = 1, a = 0.2, n = 64)"
    );
    for dt in [0.005, 0.0025] {
        let mut bf = BodyFlow::new(64, 0.55, 0.5, 1.0);
        for step in 0..6 {
            let mk = |u: f64| Disc {
                cx: 0.0,
                cy: 0.0,
                r: 0.2,
                ux: u,
                uy: 0.0,
                omega: 0.0,
            };
            let mky = |u: f64| Disc {
                cx: 0.0,
                cy: 0.0,
                r: 0.2,
                ux: 0.0,
                uy: u,
                omega: 0.0,
            };
            let (m1, _, _) = bf.trial(&mky(0.0), dt);
            let (m2, _, _) = bf.trial(&mky(0.01), dt);
            println!(
                "   lateral: Fy(0)={:.5} Fy(0.01)={:.5} K_y={:.2}",
                m1.fy,
                m2.fy,
                -(m2.fy - m1.fy) / 0.01
            );
            let (l1, _, _) = bf.trial(&mk(0.0), dt);
            let (l2, st, mesh) = bf.trial(&mk(0.01), dt);
            let (l3, _, _) = bf.trial(&mk(0.02), dt);
            println!(
                "dt={dt} step {step}: F(0)={:.4} F(0.01)={:.4} F(0.02)={:.4}  K = {:.2} / {:.2}",
                l1.fx,
                l2.fx,
                l3.fx,
                -(l2.fx - l1.fx) / 0.01,
                -(l3.fx - l2.fx) / 0.01
            );
            bf.commit(st, mesh);
        }
        println!("  expected K = m_a/dt + k = {:.2}", 0.17354 / dt + 65.4);
    }
}

fn e4h() {
    use mill_core::mac::bodies::{Body, BodyFlow};
    use mill_core::mac::staggered::Disc;
    use mill_core::mac::verify::{confined_added_mass, stokes_annulus_drag};
    let (a, nu, dt) = (0.2, 1.0, 0.005);
    let mass = std::f64::consts::PI * a * a;
    let mut bf = BodyFlow::new(64, 0.55, 0.5, nu);
    let mut body = Body {
        disc: Disc {
            cx: 0.0,
            cy: 0.0,
            r: a,
            ux: 0.0,
            uy: 0.0,
            omega: 0.0,
        },
        mass,
        inertia: 0.5 * mass * a * a,
        accel: (0.0, 0.0),
    };
    let (ma, k) = (
        confined_added_mass(a, 0.5),
        stokes_annulus_drag(a, 0.5, 1.0, nu),
    );
    for step in 0..14 {
        let info = bf.step_coupled(
            &mut body,
            (1.0, 0.0, 0.0),
            dt,
            ma,
            k,
            mill_core::mac::verify::confined_spin_stiffness(a, 0.5, nu),
            1e-6,
            40,
        );
        println!(
            "step {step}: u {:.5} v {:.2e} iters {} residual {:.2e} F_fluid {:.4}",
            body.disc.ux, body.disc.uy, info.iterations, info.residual, info.load.fx
        );
    }
}

fn e4p() {
    use mill_core::mac::bodies::BodyFlow;
    use mill_core::mac::staggered::Disc;
    let (a, nu, dt) = (0.2, 1.0, 0.005);
    let mut bf = BodyFlow::new(64, 0.55, 0.5, nu);
    let mut x = 0.0;
    for step in 0..30 {
        let u = 0.015;
        x += u * dt;
        let d = Disc {
            cx: x,
            cy: 0.0,
            r: a,
            ux: u,
            uy: 0.0,
            omega: 0.0,
        };
        let (l, st, mesh) = bf.trial(&d, dt);
        let vmax = st.v.iter().fold(0.0f64, |m, q| m.max(q.abs()));
        println!(
            "step {step}: Fx {:.5} Fy {:.2e} max|v| {:.2e} div {:.1e}",
            l.fx, l.fy, vmax, l.divergence
        );
        bf.commit(st, mesh);
    }
}

fn e4y() {
    use mill_core::mac::bodies::BodyFlow;
    use mill_core::mac::staggered::Disc;
    println!("E4y lateral force on a disc moving along x at lateral offset y (nu = 1, U = 0.015)");
    let (a, nu, dt) = (0.2, 1.0, 0.005);
    for y in [0.0, 1e-4, 1e-3, 5e-3, 2e-2] {
        let mut bf = BodyFlow::new(64, 0.55, 0.5, nu);
        let mut x = 0.0;
        let mut last = (0.0, 0.0);
        for _ in 0..30 {
            x += 0.015 * dt;
            let d = Disc {
                cx: x,
                cy: y,
                r: a,
                ux: 0.015,
                uy: 0.0,
                omega: 0.0,
            };
            let (l, st, mesh) = bf.trial(&d, dt);
            bf.commit(st, mesh);
            last = (l.fx, l.fy);
        }
        println!(
            "y = {y:.0e}: Fx {:.5} Fy {:.4e}  Fy/y = {:.3e}",
            last.0,
            last.1,
            last.1 / y.max(1e-30)
        );
    }
}

fn e4d() {
    use mill_core::mac::verify::verify_falling_disc;
    println!(
        "E4d falling disc a = 0.1 released at y = 0.15 in the drum, t = 0.6 s: y at 5 sample times"
    );
    for nu in [1e-3, 1e-4] {
        for ratio in [1.0, 1.05, 1.2] {
            for n in [64usize, 128] {
                let t = Instant::now();
                let r = verify_falling_disc(n, 0.1, ratio, 0.15, nu, 0.6, 0.3);
                let ys: Vec<String> = r.y.iter().map(|y| format!("{y:.4}")).collect();
                println!(
                    "nu={nu:.0e} rho={ratio:<4} n={n:<4} steps {:>5} iters {:.1} ok={} vymax {:.3} y: {} ({:.0?})",
                    r.steps, r.iterations, r.ok, r.vy_max, ys.join(" "), t.elapsed()
                );
            }
        }
    }
}

fn e4e() {
    println!("E4e free disc a = 0.2 under torque T = 0.05 in the drum at rest (nu = 1): steady spin vs T / kappa");
    for n in [64usize, 128, 256] {
        for ratio in [0.5, 1.0, 3.0] {
            let (err, v, it) =
                mill_core::mac::verify::verify_free_spin(n, 0.2, ratio, 0.05, 1.0, 0.8, 0.01);
            println!(
                "n={n:<4} rho={ratio:<4} spin err {:+.3}%  max translation {v:.1e}  iters {it:.1}",
                100.0 * err
            );
        }
    }
}

/// Leading 2D squeeze-film force on a disc of radius `a` approaching the concave drum wall `b` at
/// speed `u` across a gap `h` (Reynolds equation, parabolic gap): `3 sqrt(2) pi nu u R^1.5 / h^1.5`.
fn squeeze_lub(a: f64, b: f64, u: f64, nu: f64, h: f64) -> f64 {
    let r_eff = a * b / (b - a);
    3.0 * 2f64.sqrt() * std::f64::consts::PI * nu * u * r_eff.powf(1.5) / h.powf(1.5)
}

fn e5a() {
    let (a, nu, u) = (0.05, 1.0, 0.01);
    let far = stokes_annulus_drag(a, 0.5, u, nu);
    println!("E5a squeeze: disc a = {a} toward the drum wall (R = 0.5), nu = 1, U = {u}");
    println!("concentric drag (far) = {far:.4e}");
    for e in [1e-4, 0.01, 0.1, 0.3, 0.44, 0.449, 0.4499] {
        println!(
            "exact e = {e:<7} -> {:.5e}  (lub {:.5e}, h = {:.2e})",
            eccentric_squeeze_force(a, 0.5, e, u, nu),
            squeeze_lub(a, 0.5, u, nu, 0.45 - e),
            0.45 - e
        );
    }
    println!(
        "{:>5} {:>6} {:>9} {:>11} {:>11} {:>9} {:>9}",
        "n", "h/dx", "h/a", "drag", "exact", "ratio", "blended"
    );
    for n in [64usize, 128, 256] {
        let dx = 1.1 / n as f64;
        for f in [16.0, 8.0, 4.0, 2.0, 1.0, 0.5, 0.25] {
            let h = f * dx;
            let x0 = 0.5 - a - h;
            let r = verify_moving_disc(n, a, x0, u, 0.0, nu, 2.0, 0.005, false);
            let model = eccentric_squeeze_force(a, 0.5, x0, u, nu);
            let blended = mill_core::mac::lubrication::blended_normal_force(
                r.drag_mean,
                mill_core::mac::lubrication::disc_curvature(a),
                mill_core::mac::lubrication::wall_curvature(0.5),
                h,
                dx,
                nu,
                u,
            );
            println!(
                "{n:>5} {f:>6} {:>9.4} {:>11.4e} {:>11.4e} {:>9.3} {:>9.3}",
                h / a,
                r.drag_mean,
                model,
                r.drag_mean / model,
                blended / model
            );
        }
    }
}

fn e5b() {
    let (a, nu, u) = (0.05, 1.0, 0.01);
    println!("E5b squeeze vs eccentricity, a = {a}");
    println!(
        "{:>5} {:>6} {:>6} {:>11} {:>11} {:>9}",
        "n", "e", "h/dx", "drag", "exact", "ratio"
    );
    for n in [128usize, 256] {
        for e in [0.0, 0.1, 0.2, 0.3, 0.35, 0.4] {
            let r = verify_moving_disc(n, a, e, u, 0.0, nu, 0.8, 0.01, false);
            let ex = eccentric_squeeze_force(a, 0.5, e.max(1e-6), u, nu);
            println!(
                "{n:>5} {e:>6} {:>6.1} {:>11.4e} {:>11.4e} {:>9.4}",
                (0.45 - e) / (1.1 / n as f64),
                r.drag_mean,
                ex,
                r.drag_mean / ex
            );
        }
    }
}

fn e5t() {
    let (a, nu, u) = (0.05, 1.0, 0.01);
    println!("E5t steadiness: squeeze, a = {a}, e = 0.4, n = 128");
    for t_end in [0.4, 0.8, 1.6, 3.2] {
        for dt in [0.01, 0.005] {
            let r = verify_moving_disc(128, a, 0.4, u, 0.0, nu, t_end, dt, false);
            let ex = eccentric_squeeze_force(a, 0.5, 0.4, u, nu);
            println!(
                "t_end {t_end:>4} dt {dt:<6} ratio {:.4} noise {:.2e}",
                r.drag_mean / ex,
                r.drag_noise
            );
        }
    }
}

fn e5c() {
    let (nu, u, b) = (1.0, 1.0, 0.5);
    println!("E5c exact squeeze force over the leading lubrication term");
    println!(
        "{:>6} {:>9} {:>11} {:>11} {:>9} {:>9}",
        "a", "h/Reff", "exact", "lead", "ex/lead", "(ex-far)/lead"
    );
    for a in [0.03, 0.05, 0.1] {
        let r_eff = a * b / (b - a);
        let far = stokes_annulus_drag(a, b, u, nu);
        for eps in [0.003, 0.01, 0.03, 0.1, 0.2, 0.4, 0.8] {
            let h = eps * r_eff;
            let ex = eccentric_squeeze_force(a, b, b - a - h, u, nu);
            let lead = squeeze_lub(a, b, u, nu, h);
            println!(
                "{a:>6} {eps:>9} {ex:>11.4e} {lead:>11.4e} {:>9.4} {:>9.4}",
                ex / lead,
                (ex - far) / lead
            );
        }
    }
}

/// Time for a disc pushed by a constant force to close the gap from `h0` to `h` (quasi-steady
/// Stokes motion with the exact squeeze resistance).
fn squeeze_travel_time(a: f64, b: f64, nu: f64, force: f64, h0: f64, h: f64) -> f64 {
    let steps = 4000;
    let (l0, l1) = (h.ln(), h0.ln());
    let mut t = 0.0;
    for k in 0..steps {
        let lx = l0 + (l1 - l0) * (k as f64 + 0.5) / steps as f64;
        let hh = lx.exp();
        // dt = c(h) dh / F with dh = h dlnh
        t += eccentric_squeeze_force(a, b, b - a - hh, 1.0, nu) / force * hh * (l1 - l0)
            / steps as f64;
    }
    t
}

fn e5d() {
    use mill_core::mac::bodies::{Body, BodyFlow};
    use mill_core::mac::staggered::Disc;
    use mill_core::mac::verify::confined_added_mass;
    let (a, b, nu, force, h0) = (0.05, 0.5, 1.0, 5.0, 0.06);
    let levels = [
        0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002,
    ];
    println!("E5d disc a = {a} pushed to the drum wall by F = {force}, nu = {nu}, rho_s = 1.2, h0 = {h0}");
    for n in [64usize, 128] {
        let dx = 1.1 / n as f64;
        for lub in [true, false] {
            let mut bf = BodyFlow::new(n, 0.55, b, nu);
            bf.lubrication = lub;
            let mass = 1.2 * std::f64::consts::PI * a * a;
            let c0 = eccentric_squeeze_force(a, b, b - a - h0, 1.0, nu);
            let mut body = Body {
                disc: Disc {
                    cx: b - a - h0,
                    cy: 0.0,
                    r: a,
                    ux: force / c0,
                    uy: 0.0,
                    omega: 0.0,
                },
                mass,
                inertia: 0.5 * mass * a * a,
                accel: (0.0, 0.0),
            };
            let ma = confined_added_mass(a, b);
            let dt = 0.005;
            let mut t = 0.0;
            let mut next = 0;
            let mut iters = 0usize;
            let mut steps = 0usize;
            let start = Instant::now();
            println!("n={n} (dx={dx:.4}) lubrication={lub}");
            println!(
                "{:>9} {:>8} {:>9} {:>9} {:>9}",
                "gap", "gap/dx", "t_sim", "t_ref", "ratio"
            );
            while next < levels.len() && t < 12.0 {
                let h_before = b - a - body.disc.cx;
                let k = eccentric_squeeze_force(a, b, body.disc.cx, 1.0, nu);
                let info = bf.step_coupled(
                    &mut body,
                    (force, 0.0, 0.0),
                    dt,
                    ma,
                    k,
                    4.0 * std::f64::consts::PI * nu * a * a,
                    1e-3,
                    40,
                );
                t += dt;
                steps += 1;
                iters += info.iterations;
                let h = b - a - body.disc.cx;
                while next < levels.len() && h <= levels[next] {
                    let frac = (h_before - levels[next]) / (h_before - h).max(1e-15);
                    let ts = t - dt + dt * frac;
                    let tr = squeeze_travel_time(a, b, nu, force, h0, levels[next]);
                    println!(
                        "{:>9.4} {:>8.2} {:>9.4} {:>9.4} {:>9.3}",
                        levels[next],
                        levels[next] / dx,
                        ts,
                        tr,
                        ts / tr
                    );
                    next += 1;
                }
                if h < 1e-4 || !h.is_finite() {
                    break;
                }
            }
            println!(
                "  steps {steps}, mean iterations {:.1}, {:.1} s",
                iters as f64 / steps.max(1) as f64,
                start.elapsed().as_secs_f64()
            );
        }
    }
}

fn e5e() {
    use mill_core::mac::bodies::BodyFlow;
    use mill_core::mac::lubrication::{links, model_forces, remove_grid_normal};
    use mill_core::mac::staggered::Disc;
    use mill_core::mac::verify::disc_pair_squeeze_force;
    let (a, nu, u) = (0.05, 1.0, 0.01);
    println!(
        "E5e two equal discs a = {a} squeezed towards each other at +-{u} (drum R = 0.5), nu = 1"
    );
    println!(
        "{:>5} {:>6} {:>9} {:>11} {:>11} {:>9} {:>9}",
        "n", "h/dx", "h/a", "grid", "exact", "ratio", "blended"
    );
    for n in [128usize, 256] {
        let dx = 1.1 / n as f64;
        for f in [16.0, 8.0, 4.0, 2.0, 1.0, 0.5, 0.25] {
            let h = f * dx;
            let c = a + 0.5 * h;
            let discs = [
                Disc {
                    cx: -c,
                    cy: 0.0,
                    r: a,
                    ux: u,
                    uy: 0.0,
                    omega: 0.0,
                },
                Disc {
                    cx: c,
                    cy: 0.0,
                    r: a,
                    ux: -u,
                    uy: 0.0,
                    omega: 0.0,
                },
            ];
            let mut bf = BodyFlow::new(n, 0.55, 0.5, nu);
            let dt = 0.005;
            let mut fx = 0.0;
            for _ in 0..400 {
                let (loads, _, state, mesh) = bf.trial_many(&discs, dt);
                bf.commit(state, mesh);
                fx = loads[0].0;
            }
            let exact = disc_pair_squeeze_force(a, h, u, nu);
            let lk = links(&[(-c, 0.0), (c, 0.0)], &[a, a], 0.5, dx, nu);
            let mut forces = vec![(fx, 0.0), (-fx, 0.0)];
            remove_grid_normal(&mut forces, &lk);
            let model = model_forces(&lk, &[(u, 0.0), (-u, 0.0)]);
            let blended = -(forces[0].0 + model[0].0);
            println!(
                "{n:>5} {f:>6} {:>9.4} {:>11.4e} {:>11.4e} {:>9.3} {:>9.3}",
                h / a,
                -fx,
                exact,
                -fx / exact,
                blended / exact
            );
        }
    }
}

fn pair_travel_time(a: f64, nu: f64, force: f64, h0: f64, h: f64) -> f64 {
    use mill_core::mac::verify::disc_pair_squeeze_force;
    let steps = 4000;
    let (l0, l1) = (h.ln(), h0.ln());
    let mut t = 0.0;
    for k in 0..steps {
        let lx = l0 + (l1 - l0) * (k as f64 + 0.5) / steps as f64;
        let hh = lx.exp();
        // Each disc moves at U = F / c, the gap closes at 2U: dt = c dh / (2 F).
        t +=
            disc_pair_squeeze_force(a, hh, 1.0, nu) / (2.0 * force) * hh * (l1 - l0) / steps as f64;
    }
    t
}

fn e5f() {
    use mill_core::mac::bodies::{Body, BodyFlow, BodyStiffness};
    use mill_core::mac::staggered::Disc;
    use mill_core::mac::verify::disc_pair_squeeze_force;
    let (a, nu, force, h0) = (0.05, 1.0, 5.0, 0.06);
    let levels = [
        0.03, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002,
    ];
    println!("E5f two discs a = {a} pushed together by F = {force} each, nu = {nu}, rho_s = 1.2, h0 = {h0}");
    for n in [64usize, 128] {
        let dx = 1.1 / n as f64;
        for lub in [true, false] {
            let mut bf = BodyFlow::new(n, 0.55, 0.5, nu);
            bf.lubrication = lub;
            let mass = 1.2 * std::f64::consts::PI * a * a;
            let c0 = disc_pair_squeeze_force(a, h0, 1.0, nu);
            let u0 = force / c0;
            let cc = a + 0.5 * h0;
            let mk = |cx: f64, ux: f64| Body {
                disc: Disc {
                    cx,
                    cy: 0.0,
                    r: a,
                    ux,
                    uy: 0.0,
                    omega: 0.0,
                },
                mass,
                inertia: 0.5 * mass * a * a,
                accel: (0.0, 0.0),
            };
            let mut bodies = [mk(-cc, u0), mk(cc, -u0)];
            let setup = [BodyStiffness {
                added_mass: std::f64::consts::PI * a * a,
                drag_stiffness: c0,
                spin_stiffness: 4.0 * std::f64::consts::PI * nu * a * a,
            }; 2];
            let ext = [(force, 0.0, 0.0), (-force, 0.0, 0.0)];
            let dt = 0.005;
            let mut t = 0.0;
            let mut next = 0;
            let (mut iters, mut steps) = (0usize, 0usize);
            let start = Instant::now();
            println!("n={n} (dx={dx:.4}) lubrication={lub}");
            println!(
                "{:>9} {:>8} {:>9} {:>9} {:>9}",
                "gap", "gap/dx", "t_sim", "t_ref", "ratio"
            );
            while next < levels.len() && t < 12.0 {
                let h_before = bodies[1].disc.cx - bodies[0].disc.cx - 2.0 * a;
                let info = bf.step_coupled_many(&mut bodies, &ext, &setup, dt, 1e-3, 40);
                t += dt;
                steps += 1;
                iters += info.iterations;
                let h = bodies[1].disc.cx - bodies[0].disc.cx - 2.0 * a;
                while next < levels.len() && h <= levels[next] {
                    let frac = (h_before - levels[next]) / (h_before - h).max(1e-15);
                    let ts = t - dt + dt * frac;
                    let tr = pair_travel_time(a, nu, force, h0, levels[next]);
                    println!(
                        "{:>9.4} {:>8.2} {:>9.4} {:>9.4} {:>9.3}",
                        levels[next],
                        levels[next] / dx,
                        ts,
                        tr,
                        ts / tr
                    );
                    next += 1;
                }
                if h < 1e-4 || !h.is_finite() {
                    break;
                }
            }
            println!(
                "  steps {steps}, mean iterations {:.1}, {:.1} s",
                iters as f64 / steps.max(1) as f64,
                start.elapsed().as_secs_f64()
            );
        }
    }
}

fn e5g() {
    use mill_core::mac::bodies::BodyFlow;
    use mill_core::mac::staggered::Disc;
    let (a, nu, u) = (0.05, 1.0, 0.01);
    let r_eff = a * 0.5 / (0.5 - a);
    println!("E5g disc a = {a} sliding along the drum wall at U = {u} (no spin), nu = 1");
    println!(
        "{:>5} {:>6} {:>9} {:>11} {:>11} {:>9}",
        "n", "h/dx", "h/a", "F_slide", "lead", "F/lead"
    );
    for n in [128usize, 256] {
        let dx = 1.1 / n as f64;
        for f in [32.0, 16.0, 8.0, 4.0, 2.0, 1.0, 0.5, 0.25] {
            let h = f * dx;
            let e = 0.5 - a - h;
            let disc = Disc {
                cx: e,
                cy: 0.0,
                r: a,
                ux: 0.0,
                uy: u,
                omega: 0.0,
            };
            let mut bf = BodyFlow::new(n, 0.55, 0.5, nu);
            let mut fy = 0.0;
            for _ in 0..400 {
                let (loads, _, state, mesh) = bf.trial_many(&[disc], 0.005);
                bf.commit(state, mesh);
                fy = loads[0].1;
            }
            let lead = 2.0 * 2f64.sqrt() * std::f64::consts::PI * nu * u * (r_eff / h).sqrt();
            println!(
                "{n:>5} {f:>6} {:>9.4} {:>11.4e} {:>11.4e} {:>9.3}",
                h / a,
                -fy,
                lead,
                -fy / lead
            );
        }
    }
}
