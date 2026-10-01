# Verification suite and PBF baseline

Purpose: gate the fluid solver and the ball-fluid coupling against analytic solutions and an
energy budget, independently of the solver implementation. Added before the planned DFSPH +
boundary-particle rewrite so the rewrite is judged by numbers, not by feel. Goal it serves: power
and torque must not change when resolution, coarse-graining `k`, or sub-steps change.

Tolerance for power/torque convergence (per resolution doubling and per `k` change): **+-3 %**.

## Where things live

| What | Where |
|---|---|
| Gated tests (26) | `crates/mill-core/tests/verification.rs` |
| Solver adapter (`harness`), analytic helpers, all tolerances (`gates`) | `crates/mill-core/tests/common/verification_common.rs` |
| Table / CSV probe | `cargo run --release -p mill-core --example verification_probe -- [--csv] [--res 15,25] [--cases hydro,spin,couette,energy,dry]` |
| Full-system convergence + energy columns | `coarse_graining_probe` (`power_mean_w`, `torque_mean_nm`, `budget_*_w`, `unattributed_pct_of_shaft`, `torque_from_impulse_nm`) |
| Energy ledger (read-only) | `Simulation::energy_budget()` / `reset_energy_budget()`, `EnergyBudget` in `lib.rs` |

Tests that fail on the current PBF solver are `#[ignore]`d with the measured numbers in the reason:
`cargo test --release -p mill-core --test verification -- --ignored`. A solver swap edits only
`harness`.

## Cases and gates

1. **Hydrostatic** (63 mm drum, fill 0.35): worst-5 % compression <= 0.02, mean <= 0.001,
   KE/(M g dx) <= 0.05, wall-pressure profile (wall normal impulse per polar-angle bin vs
   rho g (y_surface - y_wall) x arc) L2 <= 0.05.
2. **Impulsive spin-up** (full disc, `R^2/nu = 0.25 s`, dt 1/1920): angular momentum vs the
   Bessel-series solution within 0.03 of `L_inf` at `t nu/R^2 = 0.02, 0.05, 0.1, 0.2`; wall angular
   impulse vs `dL/dt` (exact gravity torque of the lattice included) within 0.03.
3. **Taylor-Couette** (ball `R1 = R2/2`, prescribed rotation): ball torque and wall torque within
   3 % of `T = 4 pi mu w R1^2 R2^2 / (R2^2 - R1^2)`.
4. **Energy closure** (browser defaults, 60 %Nc, wet): `|unattributed| / shaft_work <= 0.02` over >= 1 rev.
5. **Dry sub-step convergence** (10 mm media, 2 seeds; `substeps` is capped at 16): 4 vs 16 and
   8 vs 16 within 3 %.

Energy ledger terms (J/m depth; `unattributed` = dE - (shaft - ball contact dissipation - fluid
wall slip - fluid viscous - clamp removed + interface created)): the net energy change no
instrumented term explains. A closed solver has `unattributed ~ 0`.

## Baseline: current PBF solver (2026-10-01)

Only the dry sub-step test and the hydrostatic KE gate pass; everything else fails.

| Metric (gate) | res 15 | 25 | 40 | 60 | 100 |
|---|---|---|---|---|---|
| Hydrostatic worst-5 % compression (0.02) | 0.098 | 0.27 | 0.68 | 1.45 | 4.49 |
| Hydrostatic mean compression (0.001) | 0.024 | 0.047 | 0.099 | 0.21 | 1.10 |
| Wall-pressure L2 (0.05) | 0.35 | 0.23 | 0.19 | 0.22 | 0.48 |
| Spin-up `|dL|/L_inf` at 0.02 (0.03) | 0.54 | 0.52 | 0.51 | 0.52 | 0.52 |
| Spin-up wall-torque mismatch (0.03) | 0.47 | 0.47 | 0.48 | 0.45 | 0.45 |
| Taylor-Couette ball torque error (0.03) | 0.996 | 0.999 | 1.000 | 1.000 | 1.000 |
| Taylor-Couette wall torque error (0.03) | 0.86 | 0.97 | 1.00 | 0.98 | 0.99 |
| Energy `|unattributed|/shaft` (0.02) | 0.20 | 0.17 | 0.80 | 4.35 | 8.40 |

Full-system, browser defaults, 60 %Nc, k=1, 3 seeds, 12 rev window (W per metre depth):

| | res 15 | 25 | 40 | 60 |
|---|---|---|---|---|
| Wet power | 15.86 | 15.03 | 13.01 | 10.25 |
| Wet torque (N m) | 1.498 | 1.419 | 1.229 | 0.968 |
| Wet `unattributed` (% of shaft) | +21 | +18 | +79 | +429 |
| Dry power (any res) | 3.169 | | | |

`k` sweep at res 40: dry power 3.169 / 3.641 / 3.893 W at k = 1 / 2 / 4 (+15 %, +23 %); wet power
13.01 / 10.55 / 8.92 W (-19 %, -31 %). Wet ball-contact dissipation rises 1.33 -> 2.52 W with k.
Raw outputs: `fullsys_res`, `fullsys_k`, `verification_probe` (txt/csv; scratchpad, regenerate
with the commands above).

## Reading the numbers

- **Compression error grows with resolution** even without balls: the wall step has no solid
  density, so particles stack at the wall (see PHYSICS.md section 9).
- **The fluid spins up far too slowly**: at `t nu/R^2 = 0.02` the exact angular momentum is 53 %
  of rigid rotation, PBF gives ~0-1 %. Resolution-independent, so it is a model deficiency
  (velocity-blend no-slip on wall-touching particles only), not discretisation error.
- **The ball barely drives the fluid** in Taylor-Couette: torque is ~1e-5 N m against 3e-2
  analytic, with no coupling clamp hits. The heuristic Stokes-relaxation coupling is not a
  no-slip boundary.
- **Positive `unattributed` that grows with resolution**: PBF position projection reconstructs
  velocities from over-compression at the wall (energy creation); the implicit viscosity step then
  removes it (`fluid_viscous` 11 -> 48 W, res 15 -> 60). The wet viscous "dissipation" at high
  resolution is largely numerical, which is why wet power falls with resolution. (Interpretation,
  consistent with the ledger but not separately proven.)
- **Dry ledger residual -7 % (k=1) to -12 % (k=4)**: DEM predict-step energy not captured by
  `dissipated_energy_j`'s definition; small, systematic, not investigated.
- Dry power is resolution-independent but rises with `k`; wet power falls with `k`. The old claim
  that bulk quantities are `k`-invariant (PHYSICS.md section 9) is disproven.

## DFSPH baseline (Phase B, 2026-10-01)

New fluid-only solver `crates/mill-core/src/fluid/` (DFSPH + Akinci wall boundary particles; see
its module doc for the algorithm). The browser still runs PBF until balls are coupled (Phase C).
`verification_probe --solver both` runs both solvers; the DFSPH tests in `verification.rs` run in
the default suite (hydrostatic and spin-up at res 25 and 50, power = omega x torque identity).

| Metric (gate) | res 15 | 25 | 40 | 60 | 100 |
|---|---|---|---|---|---|
| Hydrostatic worst-5 % compression (0.02) | 2.6e-4 | 3.7e-4 | 3.8e-4 | 4.5e-4 | 1.7e-4 |
| Hydrostatic mean compression (0.001) | 5.8e-5 | 7.4e-5 | 6.1e-5 | 7.7e-5 | 4.7e-5 |
| Hydrostatic wall-pressure L2 (0.05) | 0.043 | 0.033 | 0.034 | 0.034 | 0.028 |
| Spin-up `|dL|/L_inf` at 0.02 (0.03) | 0.036 | 0.012 | 0.019 | 0.022 | 0.022 |
| Spin-up at 0.05 | 0.017 | 0.001 | 0.006 | 0.006 | 0.005 |
| Spin-up at 0.1 | 0.004 | 0.017 | 0.014 | 0.013 | 0.010 |
| Spin-up at 0.2 | 0.026 | 0.026 | 0.028 | 0.031 | 0.024 |
| Spin-up wall-torque mismatch (0.03) | 6e-4 | 5e-5 | 3e-4 | 9e-4 | 2e-3 |

PBF numbers for the same rows are in the table above (compression 0.098 -> 4.49, spin-up error
~0.52). Marginal points: spin-up res 15 at 0.02 (0.036) and res 60 at 0.2 (0.031).

**Case 6, rotating partially filled drum, fluid only** (default slurry, 60 % of critical speed,
8 revolutions settle + 3 revolutions measured; the charge sloshes once per revolution with wall
power swinging between about +2.4 and -1.8 W, so the window must be a whole number of revolutions):

| | res 25 | 40 | 50 | 80 |
|---|---|---|---|---|
| Power (W per metre) | 0.174 | 0.221 | 0.241 | 0.289 |
| Torque (N m per metre) | 0.0165 | 0.0209 | 0.0227 | 0.0273 |
| `|unattributed|` / wall work | 0.52 | 0.35 | 0.32 | 0.26 |
| ms / step (native release) | 9.8 | 38.6 | 37.7 | 165 |

PBF for the same case: power 6.59 -> 5.40 W (25 -> 50), also unconverged and 20-30x larger.
**DFSPH does not meet the +-3 % convergence gate here** (+38 % for 25 -> 50, +31 % for 40 -> 80;
roughly first order in `dx`, extrapolating to ~0.31 W).

### What was found while getting here
- **Operator-splitting error dominates at this viscosity.** At 50 Pa s one step's viscous
  diffusion length (`sqrt(nu dt)` ~ 7 mm) is a quarter of the drum radius. Viscosity before the
  pressure solve lets wall friction carry the pool's weight (11.1 of 19.4 N in a still pool);
  viscosity after it undoes the pressure's decompression (persistent 13 % over-density). The
  scheme kept is gravity -> density solve -> viscosity -> second density solve: quiet hydrostatic
  state and rotating-drum power stable under halving and quartering `dt` (at 5 Pa s: 2.45 / 2.49 /
  2.45 W). Its price: the second solve re-adds kinetic energy the implicit viscosity removed, so the
  fluid energy ledger does not close (26-52 % of wall work).
- **The wall is a Navier-slip wall unless corrected.** A central pair force across the 0.87 dx gap
  between the first fluid and solid rows couples tangential velocity weakly (slip length ~2.5 dx).
  `WALL_COUPLING_FACTOR = 10` (viscosity.rs) restores no-slip against the Bessel spin-up solution;
  it is a calibrated constant, identical at every resolution. Late-time spin-up error saturates at
  ~0.026-0.031 for any larger factor.
- Wall backstop hits (a thin film dragged down the wall rests at the wall plane) and density-solve
  iteration-cap hits (deep columns, Jacobi at relaxation 0.5) are reported by the probe but not
  gated; the outcomes they would guard (density error, backstop kinetic energy ~1e-5 J) are.
- The old spin-up energy-closure gate had no settle phase and compared a transient against a tiny
  wall work; it was removed. Closure is judged on case 6 (ignored test, numbers above).
- Cost: DFSPH is 2-6x slower per outer step than PBF at res >= 40 (CG ~100-180 iterations).

## Phase C0: ball boundary particles, Taylor-Couette (2026-10-01) - gate NOT met

`Fluid::step_with_balls` / `Fluid::new_with_balls` (ball boundary rings, prescribed ball motion,
per-ball impulse and angular impulse out) are in the tree; `Simulation` still runs PBF. Taylor-Couette
(ball `R1 = R2/2` prescribed at 1 rad/s, full annulus, `R2^2/nu = 0.25 s`, settle 3 tau, measure 1 tau,
analytic torque 2.969e-2 N m per metre) on the DFSPH solver, `|torque| / analytic - 1`:

| | ball torque | wall torque |
|---|---|---|
| res 25 | +6.2 % | +14.1 % |
| res 50 | +14.9 % | +17.5 % |
| res 50, dt/2 and dt/4 | +13.2 %, +17.2 % | +15.2 %, +20.2 % |
| res 25 / 50, wall factor 3 | +2.0 % / +11.2 % | +7.0 % / +13.1 % |
| res 25 / 50, wall factor 20 | +9.3 % / +11.3 % | +18.4 % / +16.3 % |
| res 40, gravity switched off | -16.6 % | -17.1 % |

Gate (both within 3 % at res 25 and 50): **failed**. Observations:
- The error is not the wall-coupling factor (3 to 20 moves it by a few percent only) and does not
  converge with resolution; ball and wall torque disagree by 2-8 % where they must be equal in steady
  state.
- With gravity the bulk velocity profile is 20-25 % below the analytic one at both resolutions; with
  gravity off it is ~5 % above and the torques are 17 % low, so gravity changes the sign of the error
  even though a full annulus has no free surface. The seeded lattice does not fill the annulus
  exactly (voids against the Cartesian lattice, 7 % of the area at res 25 between the lattice and
  the boundary planes), which introduces a small free surface and non-axisymmetric flow.
- Applying one viscosity step to the analytic profile changes the bulk by 1 % rms per step (3 % with
  a 3 dx kernel support), i.e. the discrete viscous operator is not consistent on this curved,
  disordered configuration at the resolutions used.
Next step needs a decision (see the Phase C report): the ball-surface viscous layer is not accurate
enough to build the two-way coupling on without further wall/operator work.
