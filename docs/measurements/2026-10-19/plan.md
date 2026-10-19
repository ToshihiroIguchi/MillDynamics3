# Plan: re-measure the strict analytic gates, then decide the ball form-drag closure

## Context
`docs/TECHNOLOGY.md` (just written) §4.6 and §5.3 flag two open points:
1. The shipped PBF solver's strict analytic gates (hydrostatic, spin-up, Taylor-Couette, energy
   closure, rotating-drum convergence) are only known from the **2026-10-01** baseline. They have
   not been re-measured since the ball no-slip boundary (10-05), central-force viscosity and
   squeeze films (10-06), and the water default (10-08). Some harness inputs also changed after the
   numbers were recorded: spin-up and Taylor-Couette inherit `Params::default()` (density now 1000
   instead of 1800), energy closure now runs on water (0.001 Pa s).
2. Ball form drag (`crates/mill-core/src/pbf.rs:1865-1874`) is
   `1/tau_form = rho0 * C_D * (2 r_t) * |v_rel| / m_t` with `C_D = 1` (`pbf.rs:115`). The textbook
   2D-cylinder drag is `1/2 rho C_D d |v| v`, so the code effectively uses `C_D = 2`. No
   measurement of PBF drag against any reference exists in the repo (only Stokes-regime settling in
   `examples/pbf_bench.rs`). Also the same `beta` drives the ball-spin relaxation (`a_rot`,
   `pbf.rs:1907`), and several comments/tests are stale (Stokes-era).

Goal: (A) report honestly which strict gates the shipped solver meets today; (B) decide the
form-drag closure from a measurement against a literature reference, not by taste; (C) keep every
existing invariant and test green and document the result.

## Recommendation for form drag (hypothesis, to be confirmed by Step 2)
- Viscous (Stokes/low-Re) drag is already supplied by the ball-surface no-slip boundary in the
  implicit viscosity solve (`BALL_BC_GAIN`, `solve_viscosity_with_balls`), and the overlap push
  (step 3.5) already delivers some resolved pressure drag. The closure should therefore add only the
  **inertial/form part, high-Re limit**: `1/tau_form = (1/2) rho0 C_D (2 r_t) |v_rel| / m_t`, with
  `C_D = 1.0` (the Re -> infinity limit of White's cylinder fit `C_D ~ 1 + 10 Re^(-2/3)`, Re < 2e5;
  also Roshko-type subcritical value ~1.0-1.2). Do NOT add the `10/Re^(2/3)` term: it would double
  count the no-slip boundary's viscous drag.
- Net effect vs now: form drag halves. Because the resolved overlap push may already contribute
  and the factor 2 might have been compensating for it, **the decision is made by measurement**
  (Step 2), with a pre-registered rule.

## Steps

### Step 1 - Re-measure the strict gates (no code changes; shipped solver as-is)
Run (release, background; `cargo test` builds with LTO so the first build is slow):
```
cargo test --release -p mill-core --test verification -- --ignored --exact <name> --nocapture
```
for every PBF/`Simulation` test: `hydrostatic_compression_res{25,50}`,
`hydrostatic_wall_pressure_res{25,50}`, `spin_up_angular_momentum_res{25,50}`,
`spin_up_wall_torque_balance_res{25,50}`, `taylor_couette_ball_torque_res{25,50}`,
`taylor_couette_wall_torque_res{25,50}`, `energy_closure_unattributed` (~80 s),
`pbf_rotating_drum_convergence_25_50`, plus the default-run `hydrostatic_kinetic_energy_res*` and
`dry_substep_convergence`. DFSPH tests are not the shipped solver: run only the quick
`dfsph_hydrostatic_*`/`dfsph_spin_up_*` for completeness, skip `dfsph_rotating_drum_*` (5 and 20
min) and `dfsph_energy_closure_*` (state "not re-run").
Also `cargo run --release -p mill-core --example verification_probe -- --res 25,50 --solver pbf`
for the table form. Capture the **measured numbers**, not only pass/fail (the tests assert on gate
constants in `tests/common/verification_common.rs:14-38`).

Report: a table "gate | gate value | 2026-10-01 value | today's value | pass?", with the harness
inputs actually used (density, viscosity, dt, resolution) because they changed. Harness is pinned
to 1800 kg/m3 / 50 Pa s for hydrostatic and rotating drum (`verification_common.rs:455-456`), so
those do NOT describe the water default: state that, and add one water data point only if cheap
(the still-pool probe `examples/still_pool_probe.rs`).
Expected outcome (stated in advance, to be tested): most strict gates still fail.

### Step 2 - Measure ball drag against a reference (new probe, no solver change yet)
New example `crates/mill-core/examples/drag_probe.rs` (pattern: the Taylor-Couette adapter,
`verification_common.rs:768-819`, which overwrites a ball's state every step and calls
`FluidParticles::step_coupled`; and `pbf_bench.rs` for an isolated ball):
- A ball of d = 2 mm (also one run at coarse-graining k = 2 to check independence) at the centre of
  a quiescent pool in a drum large relative to the ball (blockage <= ~5 %), pushed at a prescribed
  constant velocity U over several diameters (ball state overwritten each step). Drag force =
  -(sum of coupling impulses on the ball)/dt averaged after the start-up transient; record the
  split by mechanism (overlap push, form-drag closure, implicit no-slip solve) from
  `CouplingImpulses`/`FluidStepStats`.
- Re = U d / nu over {1, 10, 100, 1e3, 1e4} by varying U and nu (nu = 1e-6 water; also nu >= 0.1
  where the grid reference `mac/verify.rs:1739` `stokes_annulus_drag` exists and Re is low).
- Resolutions 25, 40, 50 and `substeps` 8/16.
- Reference: `F_ref = 1/2 rho U^2 d C_D(Re)` per unit depth with White's fit
  `C_D = 1 + 10 Re^(-2/3)` (Re > ~10) and, for Re <= 1, the grid/Stokes-annulus value. Verify the
  literature values (White, Fluid Mechanics; Tritton 1959; Roshko 1961) with WebSearch before
  citing, and record exact source.
- Variants of the closure strength `k` multiplying `inv_tau_form`: **k = 0 (closure off), 0.5, 1.0
  (today)**, switched by a temporary `#[cfg(test)]`/env-free parameter in the probe path (not
  shipped). Output ratio `R = F_meas / F_ref` per (Re, resolution, k).

Pre-registered decision rule (Re >= 1e3 is where form drag matters; Re ~ 1-10 is no-slip territory):
1. If some k has `0.75 <= R <= 1.33` for Re in {1e3, 1e4} at resolutions 25-50 -> choose that k,
   and prefer the physically meaningful one (k = 0.5 with `C_D = 1` == textbook).
2. If k = 0 already satisfies it -> the closure is redundant: remove step 6.5 (and the `c_rot`
   machinery it alone needs) in a separate, ablation-tested commit.
3. If none satisfies it -> keep the best k, document the residual error as a model limit; do not
   tune further.
Also check R at Re ~ 1-100 does not worsen vs k = 1 (double-counting guard).

### Step 3 - Apply the decision (numerical core: main session, per CLAUDE.md; not delegated)
Only if Step 2 selects a change in `crates/mill-core/src/pbf.rs`:
- Edit `inv_tau_form` (`:1865-1873`) to the chosen form, `BALL_FORM_DRAG_COEFFICIENT` doc comment
  (`:108-115`, stale "default 50 Pa*s"), the step 6.5 comment block (`:1702-1769`, still describes
  `tau_stokes`; `inv_tau_stokes = 0.0` stays or is deleted).
- Fix the stale test pointer `dem.rs:307` and `docs/PHYSICS.md:788` (test
  `ball_drag_relaxation_fraction_is_independent_of_coarse_graining` no longer exists).
- Re-check tests that depend on drag magnitude and adjust ONLY if the new value is justified by
  Step 2 data, with the reason in the comment:
  `pbf.rs:2670 ball_drag_matches_two_dimensional_stokes_scaling` (50 % margin, estimated ~15 % of
  its prediction is form drag), `:2736 ball_drag_gains_form_drag_at_high_relative_speed_and_low_
  viscosity` (> 100 x Stokes), `:2778 ball_drag_never_overshoots_the_local_fluid_velocity`,
  `coupling.rs:151, :578, :788, :804`, `metrics.rs:1012` (bound 0.20), and the energy invariants
  `coupling.rs:930, :1149`, `lib.rs:759`.
- Add a regression test from the probe: drag ratio R within the accepted band at one (Re, res)
  point, fast (<~10 s).

### Step 4 - Verify the change (same-session A/B, build both from `git stash`/`pop`)
1. `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`,
   `cargo test --workspace` (default suite) - all green.
2. Full-system energy invariants (the three above) with the new drag; record `worst_margin_j` before
   vs after (single-seed metric is noisy: also run 2 seeds).
3. `examples/tolerance_probe.rs` A/B, **>= 8 seeds** (3 seeds underestimated noise, VERIFICATION.md
   :1706-1708), power/torque/impact rate/toe-shoulder/compression/energy-closure gates. A physics
   change is expected to shift some metrics: report the shifts, do not hide them.
4. `examples/perf_probe.rs`: hashes WILL change (physics change; expected). Confirm speed unchanged.
5. Re-run Step 1's Taylor-Couette and energy-closure tests with the new code (form drag touches both
   ball torque via `a_rot` and the energy ledger).
6. Spot-check in the browser (Playwright or run skill) that a default run looks sane (no floating
   balls, speed badge unchanged).

### Step 5 - Documentation (delegate to sonnet/haiku per CLAUDE.md; English only)
- `docs/VERIFICATION.md`: new dated section with the Step 1 table, the Step 2 drag table (R vs Re,
  resolution, k), the decision and the rule used.
- `docs/PHYSICS.md` §6.2 / §9: current closure and the fixed stale statements.
- `docs/TECHNOLOGY.md`: update §4.6 (gate status with date) and the §5.3 form-drag bullet
  (replace "not settled" with the measured outcome); fix the Notes-level inconsistency it records.
- One commit per milestone on `main`: (1) probe + gate re-measurement docs, (2) closure change +
  tests, (3) docs. No push.

## Critical files
- `crates/mill-core/src/pbf.rs` (form-drag closure `:1702-2006`, constant `:115`)
- `crates/mill-core/examples/drag_probe.rs` (new); reuse `tests/common/verification_common.rs`
  adapter pattern, `examples/pbf_bench.rs`, `mac/verify.rs:1739`
- `crates/mill-core/tests/verification.rs`, `examples/verification_probe.rs`,
  `examples/tolerance_probe.rs`, `examples/perf_probe.rs` (existing, reused)
- docs: `VERIFICATION.md`, `PHYSICS.md`, `TECHNOLOGY.md`

## Model delegation (CLAUDE.md cost policy)
- Main session: Step 2 interpretation, Step 3 (pbf.rs numerics), review of Step 4 results.
- sonnet: write `drag_probe.rs`, regression test, doc edits.
- haiku: run the long cargo commands in the background and return the raw numbers.

## Verification of the whole task
Done when: (a) Step 1 table exists with today's numbers and honest pass/fail; (b) the drag decision
is traceable to measured R(Re) and the pre-registered rule; (c) fmt/clippy/test green, energy
invariants hold, tolerance A/B and shifts reported; (d) docs updated with dates and numbers.
Caveats to state in the report: gates run at the harness's pinned inputs (not water) for
hydrostatic/rotating drum; DFSPH rotating-drum tests skipped; timings are same-session only.
