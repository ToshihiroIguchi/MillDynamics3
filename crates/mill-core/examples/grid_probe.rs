//! Experiment tables for the fixed-grid solver track. `cargo run --release -p mill-core --example
//! grid_probe -- --exp E0`.

use mill_core::mac::verify::{
    verify_couette_ns, verify_couette_stokes, verify_manufactured, verify_spinup,
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
