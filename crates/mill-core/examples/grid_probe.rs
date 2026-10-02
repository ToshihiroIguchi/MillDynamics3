//! Experiment tables for the fixed-grid solver track. `cargo run --release -p mill-core --example
//! grid_probe -- --exp E0`.

use mill_core::mac::verify::verify_manufactured;
use std::time::Instant;

fn main() {
    let exp = std::env::args()
        .skip_while(|a| a != "--exp")
        .nth(1)
        .unwrap_or_else(|| "E0".into());
    match exp.as_str() {
        "E0" => e0(),
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
