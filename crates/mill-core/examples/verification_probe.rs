//! Prints the verification-suite metrics (tests/verification.rs) at several fluid resolutions.
//! `cargo run --release -p mill-core --example verification_probe -- [--csv] [--res 15,25,...]
//!  [--cases hydro,spin,rot,couette,energy,dry] [--solver pbf|dfsph|both]`
//! `--solver` (default `both`) applies to the fluid-only cases `hydro`, `spin` and `rot`;
//! Taylor-Couette, `energy` and `dry` always run the PBF solver (balls are not in DFSPH yet) and
//! are skipped for `--solver dfsph`. Each fluid-only case prints a `ms_per_step` row (wall-clock
//! milliseconds per outer `step` call over the whole run).
//! The setup lives in `tests/common/verification_common.rs`, shared with the test file.

#[path = "../tests/common/verification_common.rs"]
mod common;

use common::harness::Solver;
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
            ["hydro", "spin", "rot", "couette", "energy", "dry"]
                .map(String::from)
                .to_vec()
        });
    let has = |c: &str| cases.iter().any(|x| x == c);
    let solver_arg = arg_value("--solver").unwrap_or_else(|| "both".to_string());
    let solvers: Vec<Solver> = match solver_arg.as_str() {
        "pbf" => vec![Solver::Pbf],
        "dfsph" => vec![Solver::Dfsph],
        "both" => vec![Solver::Pbf, Solver::Dfsph],
        other => panic!("--solver must be pbf, dfsph or both (got {other})"),
    };
    let run_pbf = solvers.contains(&Solver::Pbf);

    if csv {
        println!("case,solver,res,metric,value,gate");
    } else {
        println!(
            "{:<15} {:<5} {:>4}  {:<38} {:>14}  gate",
            "case", "solver", "res", "metric", "value"
        );
    }
    let emit = |rows: Vec<Row>| {
        for r in rows {
            if csv {
                println!(
                    "{},{},{},{},{:.8e},{}",
                    r.case, r.solver, r.res, r.metric, r.value, r.gate
                );
            } else {
                println!(
                    "{:<15} {:<5} {:>4}  {:<38} {:>14.6e}  {}",
                    r.case, r.solver, r.res, r.metric, r.value, r.gate
                );
            }
        }
    };
    for &res in &res_list {
        for &solver in &solvers {
            if has("hydro") {
                emit(report::hydrostatic(solver, res));
            }
            if has("spin") {
                emit(report::spin_up(solver, res));
            }
            if has("rot") {
                emit(report::rotating_drum(solver, res));
            }
        }
        if run_pbf && has("couette") {
            emit(report::taylor_couette(res));
        }
        if run_pbf && has("energy") {
            emit(report::energy(res));
        }
    }
    if run_pbf && has("dry") {
        emit(report::dry());
    }
}
