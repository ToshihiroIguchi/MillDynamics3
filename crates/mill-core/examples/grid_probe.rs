//! Experiment tables for the fixed-grid solver track. `cargo run --release -p mill-core --example
//! grid_probe -- --exp E0`.

use mill_core::mac::verify::{
    diagnose_couette_mac, verify_couette_mac, verify_couette_ns, verify_couette_stokes,
    verify_manufactured, verify_spinup, verify_spinup_mac,
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
    let ns = if ns.is_empty() {
        vec![64, 128, 256]
    } else {
        ns
    };
    for re in [10.0f64, 1000.0] {
        let nu = 0.25 * 0.25 / re;
        for &n in &ns {
            let t = Instant::now();
            let r = verify_couette_mac(n, nu, 10.0, 0.25);
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
