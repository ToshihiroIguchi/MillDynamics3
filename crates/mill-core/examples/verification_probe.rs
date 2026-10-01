//! Prints the verification-suite metrics (tests/verification.rs) at several fluid resolutions.
//! `cargo run --release -p mill-core --example verification_probe -- [--csv] [--res 15,25,...]
//!  [--cases hydro,spin,couette,energy,dry]`
//! The setup lives in `tests/common/verification_common.rs`, shared with the test file.

#[path = "../tests/common/verification_common.rs"]
mod common;

use common::report::{self, Row};

fn arg_value(flag: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    let csv = std::env::args().any(|a| a == "--csv");
    let res_list: Vec<u32> = arg_value("--res")
        .map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect())
        .unwrap_or_else(|| vec![15, 25, 40, 60, 100]);
    let cases: Vec<String> = arg_value("--cases")
        .map(|v| v.split(',').map(str::to_string).collect())
        .unwrap_or_else(|| {
            ["hydro", "spin", "couette", "energy", "dry"]
                .map(String::from)
                .to_vec()
        });
    let has = |c: &str| cases.iter().any(|x| x == c);

    if csv {
        println!("case,res,metric,value,gate");
    } else {
        println!(
            "{:<15} {:>4}  {:<38} {:>14}  gate",
            "case", "res", "metric", "value"
        );
    }
    let emit = |rows: Vec<Row>| {
        for r in rows {
            if csv {
                println!(
                    "{},{},{},{:.8e},{}",
                    r.case, r.res, r.metric, r.value, r.gate
                );
            } else {
                println!(
                    "{:<15} {:>4}  {:<38} {:>14.6e}  {}",
                    r.case, r.res, r.metric, r.value, r.gate
                );
            }
        }
    };
    for &res in &res_list {
        if has("hydro") {
            emit(report::hydrostatic(res));
        }
        if has("spin") {
            emit(report::spin_up(res));
        }
        if has("couette") {
            emit(report::taylor_couette(res));
        }
        if has("energy") {
            emit(report::energy(res));
        }
    }
    if has("dry") {
        emit(report::dry());
    }
}
