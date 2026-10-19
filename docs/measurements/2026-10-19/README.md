# Measurements of 2026-10-19 (saved before a shutdown; work interrupted)

Raw logs of the strict-gate re-measurement and the ball-drag measurement. Nothing here is
committed yet; the VERIFICATION/PHYSICS/TECHNOLOGY docs are NOT updated with these numbers.

## Files
- `gates/step1_summary.md` - strict gates, 2026-10-01 vs 2026-10-19, shipped PBF at commit f8c2d16
  (6 pass / 16 fail of 22 run). `gates/verification_all.log` is the raw output.
- `drag/signed_res25.log`, `signed_res25_s16.log`, `signed_res50_s16.log` - VALID drag
  measurements (+U / -U difference, standard errors). `ablate_wet0_st0.log` - adhesion and
  cohesion switched off (VALID). `smoke.log`, `scan_a.log` - INVALID (first probe version, U = 0
  baseline subtraction was noise dominated); keep only as history.
  The `k` column in the valid logs is the multiplier on the form-drag rate (0, 0.5, 1) from a
  temporary hook that has since been removed from `pbf.rs`.
- `ab/base.csv`, `ab/base.log` - `tolerance_probe` baseline, 8 seeds, configs
  ship_fast / ui_fast / ship_v0.01, recorded with the code of f8c2d16 (hook = 1.0, bit-identical).
- `refs_and_survey.md` - the verified literature list used for docs/TECHNOLOGY.md.
- `plan.md` - the approved plan for this task.

## Key results
- Ball drag in water (d = 10 mm, d/dx 2-4, quiescent full pool): total drag is 5-17x the cylinder
  reference at Re 1e3-3e3 (F ~ 0.7-0.9 N/m at U = 0.1 m/s), nearly independent of the form-drag
  closure strength k, of resolution (d/dx 2 -> 4) and of sub-steps (8 -> 16); adhesion/cohesion off
  only lowers it by ~25 %. At Re 1e4 it is ~0.9-1.9x and depends on sub-steps (4.6 vs 8.0 N/m,
  res 25). Closure k = 0 / 0.5 / 1 changes the total by <= ~20 % (within noise). Cause of the excess
  drag is NOT identified (candidates: overlap push, density-constraint noise).
- Reference: White's fit `C_D = 1 + 10 Re^(-2/3)` could not be confirmed in a literature search;
  the subcritical plateau C_D ~ 1.2 (Roshko 1961, Zdravkovich, Achenbach) was. Treat the
  reference as C_D 1.0-1.2 (about +-20 %).

## State of the working tree when work stopped
- `crates/mill-core/src/pbf.rs`: form-drag rate now has the textbook factor 1/2
  (`0.5 * rest_density * BALL_FORM_DRAG_COEFFICIENT * 2 r_t * |v_rel| / m_t`), temporary hook
  `FORM_DRAG_SCALE` removed, comments updated. UNVERIFIED, UNCOMMITTED.
- `crates/mill-core/examples/drag_probe.rs`: new probe, UNCOMMITTED.

## To resume (not yet run)
1. `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`,
   `cargo test --workspace`.
2. `cargo run --release -p mill-core --example tolerance_probe -- --compare docs/measurements/2026-10-19/ab/base.csv --seeds 8 --only ship_fast,ui_fast,ship_v0.01`
3. `cargo test --release -p mill-core --test verification -- --include-ignored --skip dfsph --nocapture --test-threads=3`
   (compare Taylor-Couette and energy closure with `gates/step1_summary.md`).
4. `./target/release/examples/drag_probe --re 1000,3000,10000 --res 25 --substeps 8 --csv`
   (and `--substeps 16`, `--res 50 --substeps 16`) and compare with the `k = 1` rows above.
5. If green: update docs/VERIFICATION.md (new dated section), PHYSICS.md section 6.2/9 and
   TECHNOLOGY.md 4.6 / 5.3 with these numbers; fix the stale test pointer at `dem.rs:307` /
   `PHYSICS.md:788`; commit (no push).
