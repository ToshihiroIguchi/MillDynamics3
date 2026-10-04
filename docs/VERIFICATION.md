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

## Smooth solid boundary (2026-10-01, replaces the boundary-particle layer) - gate still NOT met

Diagnosis behind the change (all measured): the bulk viscous operator is accurate (pure-bulk residual
0.7-1.4 % on the analytic Taylor-Couette profile; effective viscosity on the hex lattice 0.965 nu,
isotropic to 4 digits over orientation); the pressure part of the ball torque is 1e-4..1e-3 N m
against 3e-2 total (pressure roughness is not the problem); but the viscous torque depends on
how the boundary layer is represented (discrete rows, 0.87 dx gap, density at the boundary) far more
than 3 % allows, and a mass that is ~1 % too large for the available space pushes particles into the
boundary and inflates the torque (+23 % at res 25, +69 % at res 50).
The wall and balls are now smooth solids (`fluid/solid.rs`): tabulated kernel integrals `Gamma`
(volume fraction) and `T` (viscous coupling tensor) per body radius, `grad Gamma` normal to the
surface (no pressure torque on a ball, exact buoyancy), no calibration factor. The particle mass is
`rho0 * fill * area / N` exactly.

Taylor-Couette (ball R1 = R2/2 at 1 rad/s, fill 0.75 so the fluid mass equals the annulus, settle
3 tau, measure 1 tau, `|torque| / analytic - 1`, gravity off, dt = 1/480 unless noted):

| res | 15 | 25 | 28 | 35 | 50 |
|---|---|---|---|---|---|
| ball torque | -16.6 % | -13.7 % | -20.2 % | -18.3 % | -14.7 % |
| wall torque | -16.6 % | -14.6 % | -20.8 % | -19.0 % | -15.3 % |

res 35 vs time step: dt = 1/480 -18.3 %, dt/4 -3.1 %, dt/16 +6.9 % (first order in dt: the
pressure-viscosity operator splitting, since `dt nu / dx^2 ~ 6-20`). Other orderings at res 35:
without the second density solve -8.9 %, without the first -24 %, incremental pressure correction
(previous kappa before viscosity) -13.4 % (dt) / -11.6 % (dt/4). With gravity (res 25 / 50) the
torque is +1.6 / +66 % before the exact-mass rule and ~-14 % after. Radial velocity noise is 5 % of
the tangential velocity and does not explain it (eddy diffusion ~1e-7 against nu 4e-3).

Phase B gates on the smooth boundary without any calibration factor (previously met with the
calibrated particle boundary): hydrostatic wall-pressure L2 0.104 / 0.060 (res 25 / 50, gate 0.05);
spin-up `|dL|/L_inf` at t/tau 0.02 / 0.05 / 0.1 / 0.2: 0.057 / 0.037 / 0.014 / 0.012 (res 25),
0.042 / 0.021 / 0.003 / 0.024 (res 50) with the angular momentum normalised by the *current* second
moment (the lattice relaxes after seeding; the old seeded-moment normalisation was biased by up to
3 %). The error at early times falls with resolution (0.080 / 0.057 / 0.047 / 0.042 at res 15 / 25 /
40 / 60) but is above 0.03; the DFSPH hydrostatic and spin-up tests are `#[ignore]`d with these
numbers.

Conclusion: at the production time step the torque error is dominated by a first-order
pressure-viscosity splitting error; it can be bought down only by a smaller step (x4 for -3 % at
res 35, and `dt nu / dx^2 <= ~3` would mean ~25 internal steps at the defaults) or by a monolithic
pressure-viscosity solve.

## Correction: steady-state Taylor-Couette (2026-10-02) - earlier TC numbers were transients

The Taylor-Couette runs above settled for 3 tau and measured 1 tau. The torque keeps relaxing for
about 10 tau (slow drift, plus stick-slip bursts of +-10..20 % at small time steps), so those numbers
were not steady-state values, and the statement "first order in dt" (-18 / -3 / +7 % at dt, dt/4,
dt/16) was a coincidence of unsteady windows and is **retracted**. Steady-state protocol: settle 8 tau,
measure 2-4 tau, gravity off, smooth boundary, exact-fill mass, `|torque| / analytic - 1` (ball and
wall torque agree to 0.1 %):

| | res 15 | 25 | 35 | 50 |
|---|---|---|---|---|
| dt nu / dx^2 <= 3, solid coupling x1 | -23.6 % | -25.7 % | -24.0 % | -17.8 % |
| dt nu / dx^2 <= 3, solid coupling x6 | -12.3 % | - | -14.3 % | -11.2 % |
| dt nu / dx^2 <= 3, solid coupling x10 | -7.0 % | -1.7 % | -9.0 % | -11.8 % |

res 25, solid coupling x1, versus the viscous number `C = dt nu / dx^2` per internal step: C = 12 / 3 /
1 / 0.3 gives -26 / -26 / -18 / -20 % (at 0.3 the torque oscillates +-5 % in bursts). Sensitivities at
res 25: bulk viscous coefficient x1.25 gives -10.6 % (x0.8 gives -38 %), i.e. torque ~ coefficient^0.8.

Conclusions: (1) shrinking the time step moves the steady error by only a few percent, so operator
splitting is not the main error; (2) the solid viscous coupling is too weak by a large factor (x10
brings res 25 to -1.7 %), but a single factor does not give a resolution-independent result
(-7 / -2 / -9 / -12 %), i.e. no calibration reaches +-3 % across resolutions; (3) the bulk
operator itself is worth +3 % (its 0.965 nu lattice factor). `Fluid::max_viscous_number` and
`max_internal_steps` stay as hooks (default: off).

## Phase C2: Simulation switched to DFSPH (2026-10-03) - coupled charge not yet stable

`Simulation` now runs `fluid::Fluid` (balls are unknowns of the implicit viscous solve, wall and balls
are smooth solids); `FluidView` lets metrics/surface read either solver; `coupling::step_dfsph` is the
new step. PBF stays in the tree for its own tests and `perf_probe`.

Behaviour of the coupled default scene (244 balls of 2 mm, fill 0.30, slurry 0.35, 50 Pa s) over the
first 150 DEM sub-steps:
- res 25 (ball radius 0.8 dx): without a pore handler, ball speed bursts to 4-19 m/s (dry: 0.2 m/s) the
  moment the charge settles on the wall. Cause: fluid particles trapped in pores narrower than the
  kernel (ball against wall / ball) are squeezed out by the density solve at ~4 m/s and the reaction
  goes to the balls. A position-level pore handler (`GAMMA_MAX = 0.85`, particles moved down the
  summed-volume-fraction gradient, counted as backstop hits) keeps res 25 bounded (ball speed <= 2.4
  m/s early, ~0.5 m/s settled).
- res 60 (ball diameter 3.8 dx): still bursts after ~130 sub-steps (ball and fluid up to 12 m/s,
  internal steps pinned at the cap of 32, 300+ backstop hits per step).
- 40 large balls (diameter ~6 dx) at res 40: stable (<= 0.8 m/s), but ~100 ball-backstop hits per
  sub-step.
Four tests written against the PBF behaviour are `#[ignore]`d with reasons (three full-system stability
tests, one compression-error test). Not yet done: C3 measurements (power vs resolution, energy closure,
ms/frame) because the resolved dense charge is not stable enough to measure.

## Grid-solver track — E0 (infrastructure), 2026-10-03

Decision: replace the particle fluid with a fixed-grid incompressible solver, built as gated
experiments E0..E9 (module `crates/mill-core/src/mac/`; `grid.rs` is the unrelated spatial hash).
Scope: water (1e-3 Pa·s) to yield-stress paste. Stop rules and the decisions D1 (water target)
and D2 (coarse-graining validity) are the user's.

E0 result (`grid_probe --exp E0`, Neumann Poisson `p = cos(pi r^2/R^2)` on a cut-cell disc, face
apertures exact, volume-weighted L2 error against the analytic field):

| n | L2 error | order | PCG iterations (1e-12) |
|---|---|---|---|
| 32 | 2.39e-3 | – | 12 |
| 64 | 4.52e-4 | 2.41 | 13 |
| 128 | 1.28e-4 | 1.82 | 15 |
| 256 | 3.65e-5 | 1.81 | 16 |
| 512 | 6.36e-6 | 2.52 | 18 |

Projection of a random face field: divergence 1.8 -> 1e-14 (round-off), bit-identical reruns.
Gate change versus the plan: the multigrid is a piecewise-constant aggregation (exact Galerkin
coarse operators, over-correction 1.8) used as a CG preconditioner, so "residual reduction per
V-cycle" is replaced by "PCG iteration count nearly grid independent" (12 -> 18 over 16x cells,
slow log growth). Gate met; order alternates 1.8..2.5 with the cut-cell geometry.

### E1 — Newtonian viscous flow (2026-10-03)

**E1a, viscous operator alone: gate met.** Cell-centred componentwise Dirichlet Helmholtz
(Gibou symmetric `1/theta` form) on the level-set region, multigrid-preconditioned CG.
- Steady Stokes Taylor-Couette (R1 = 0.25, R2 = 0.5): max velocity error 4e-3 -> 2e-5 (n 32 -> 512,
  second order); wall torque error 0.03 % at n = 32 and <= 0.03 % up to n = 512, both walls.
- The wall load is the discrete link reaction `mu (u_b - u_c) / theta`. It is the traction of the
  component-wise operator `mu grad(u).n`; the physical traction `mu (grad u + grad u^T).n` differs by
  `-mu Omega t` on a rigid wall rotating at `Omega`, i.e. a torque `+-2 mu Omega A` (A = enclosed area,
  `+` body in fluid, `-` drum enclosing fluid). Without it the inner torque is off by exactly
  `-2 pi mu omega R1^2` (-37.5 % here). Check: fluid in rigid rotation exerts no torque.
- Impulsive spin-up of a disc (BDF2, start-up by 32 backward-Euler substeps, previous level = rest):
  `L(t)` error <= 0.2 % for `nu t/R^2 >= 0.05` and torque error <= 0.3 % for 0.05..0.1 at
  `dt = 0.0025 R^2/nu`; second order in dt; identical for n = 64, 128, 256 (no spatial error).
  (Torque at `nu t/R^2 >= 0.3` is a 1 % difference of large cancelling terms; not used as a gate.)

**E1b, full Navier-Stokes hold test (exact Taylor-Couette state, integrate to t = 10): gate NOT met.**
Cell-centred velocity, aperture-weighted exact face projection (E0 operator), IMEX BDF2/AB2,
incremental pressure correction, two-layer ghost extension of the velocity, central advection.

| n | Re | u err (of omega R1) | inner T | outer T |
|---|---|---|---|---|
| 64 / 128 / 256 | 10 | 9.1e-3 / 5.0e-3 / 2.8e-3 | -1.5 / +0.5 / -0.2 % | -3.9 / -2.7 / -1.8 % |
| 64 / 128 / 256 | 1000 | 1.0e-2 / 7.6e-3 / 4.2e-3 | -0.7 / +6.2 / +2.2 % | -8.8 / -3.1 / -1.7 % |

Velocity converges at first order (0.85), the outer torque at order ~0.5, the inner torque is noisy.
Gate was torque <= 1 % with order >= 1.5. Face divergence is round-off (exact projection holds).

Targeted fixes tried, in order (each judged on the same hold test):
1. Cells whose centre is solid but that are mostly fluid (cut slivers) had no velocity and were
   replaced by a neighbour value: max face divergence of the exact field 0.11. Filling them by
   extrapolation: 3e-3. (Large gain.)
2. Aperture-normalised / then central advection with ghost values: near-wall advection error of
   the exact field 3.1e-2 -> 1.2e-2; interior second order. (Large gain.)
3. Initial pressure from the projected acceleration (small gain, needed for a clean start).
4. Quadratic instead of linear ghost extrapolation: ghost error max 1.1e-2 -> 5.1e-3, no gain in torque.
5. Rotational incremental pressure form (`p += phi - nu div u*`): WORSE (u err 1.8e-1 at n = 256,
   Re 0.1): the cut-cell divergence is too noisy. Reverted.

Diagnosis: the viscous operator is superconvergent (E1a), the projection is exact in the face
sense, but the collocated cell-velocity / face-flux pairing is inconsistent at cut walls. The
centripetal balance `u.grad u = grad p` is violated by O(1) truncation error in cells whose centre is
within ~dx of the wall (adv error max 1.1e-2 not decreasing with n; rms decreasing only ~dx^0.3),
which feeds a first-order velocity error into the wall layer and hence the wall torque. Stopped here
per the stop rule; options are reported to the user.

### E1b remediation — R0 diagnostics and R2 staggered scheme (2026-10-03)

**R0.1 (time-step dependence).** The cell-centred scheme (`mac/flow.rs`) was rerun at CFL 0.25 and
0.0625 (n = 64/128, Re 10/1000): velocity error and wall torques are unchanged to within noise
(e.g. n = 128, Re 1000: u err 7.60e-3 and torque +6.236 % / -3.135 % at both CFL). The error is
therefore spatial, not fractional-step splitting, so R1 was dropped and R2 was built.

**R2 (`mac/staggered.rs`, `grid_probe --exp R2`).** Face-normal velocities are the unknowns, Gibou
Dirichlet Helmholtz on each face lattice, E0 weighted projection, quadratic ghost extension,
central advection (AB2) and BDF2 viscosity. Two fixes were needed:
1. ghost faces must carry the analytic/extrapolated value on the lattice before the first
   projection (step-0 error);
2. the face flux must be sampled at the centroid of the open part of a cut face:
   `flux = aperture * (u_node + off * du/dt)`, with the correction held fixed during the projection
   so it stays exact. Without it the steady error was O(1) in divergence at cut cells and did not
   converge (u err 8e-3 at n = 64 and 128); with it 8e-4 / 3.5e-4.

Hold test from the exact Taylor-Couette state, t_end = 10, CFL 0.25:

| n | Re | u err | inner T | outer T | div | ms/step |
|---|----|-------|---------|---------|-----|---------|
| 64 | 10 | 7.98e-4 | -0.080 % | -0.102 % | 1e-16 | 3.1 |
| 128 | 10 | 3.47e-4 | 0.055 % | -0.103 % | 8e-17 | 28 |
| 256 | 10 | 8.72e-5 | 0.015 % | -0.330 % | 8e-17 | 153 |
| 64 | 1000 | 4.68e-3 | 1.263 % | 0.013 % | 8e-17 | 4.4 |
| 128 | 1000 | 2.71e-4 | 0.021 % | -0.075 % | 8e-17 | 18 |
| 256 | 1000 | 7.79e-5 | 0.042 % | -0.345 % | 8e-17 | 88 |

Gate: torque <= 1 % at mid resolution **met** (<= 0.11 % at n = 128), velocity order 1.2 (64->128)
and 2.0 (128->256) at Re 10, 4.1 and 1.8 at Re 1000 (**met** on the finer pair), face divergence
round-off **met**. Caveats: (a) the outer torque at n = 256 (-0.33 %) is slightly worse than at
n = 128, still far inside 1 %; (b) the NS spin-up rerun on the staggered stepper is not done yet
(the Stokes operator case is unchanged); (c) cost is dominated by rebuilding the multigrid
hierarchies every step (88-153 ms/step at n = 256) and must be cached before E3; (d) advection is
plain central, adequate for this smooth test but a robust scheme is still needed for water.

**R2 follow-up: spin-up rerun on the staggered NS stepper** (`grid_probe --exp R2s`; azimuthal flow, so
the Stokes Bessel solution is exact for NS). Angular momentum `L(t)` error: n = 256, dt = 0.005 nu/R^2:
-1.0 % (T = 0.02, dt-limited start-up), -0.10 % (0.05), -0.02 % (0.1), 0.003 % (0.3), at least as
good as the cell-centred E1a stepper. Wall torque error: 2.5 / 0.22 / 1.07 % at T = 0.02 / 0.05 / 0.1
(E1a: -10.5 / -1.4 / -0.55 % at the same dt). Late-time torque ratios are meaningless
(dL/dt -> 0, relative error blows up); in absolute terms the rigid-rotation residual converges
fast (0.40 -> 6.6e-3 -> 6e-4 for n = 64/128/256, against the 2 mu Omega A = 1.57 scale), but the
transient torque at T = 0.3 has an absolute error ~6e-3 (0.4 % of that scale) that does not
shrink between n = 128 and 256 (cut-cell jitter). Judged acceptable for the gate; revisit in E4.

Solver tolerances for the staggered stepper were relaxed 1e-12 -> 1e-9 (relative residual) with
identical errors and projected divergence <= 8e-15; Helmholtz hierarchies are cached per sigma.
Cost at n = 128: 15-19 ms/step (was 20-30).

## Grid-solver track — E2 (free surface), status 2026-10-03

Code: `mac/levelset.rs` (WENO5 + TVD-RK3 level set, Newton-contour reinitialisation, volume
correction), free-surface layer in `mac/staggered.rs` (`Liquid` topology rebuilt per step from a
predicted interface, ghost-fluid Neumann/Dirichlet pressure with `1/theta`, masked viscous
Helmholtz, gravity as a face force, predictor/corrector level-set coupling, optional third-order
upwind advection), reference solutions in `mac/reference.rs`. Probes: `grid_probe --exp E2ls | E2a |
E2b | E2bc | E2c | E2ref`.

* **Level set** (rigid rotation, 1 revolution, n = 64/128/256): disc volume drift -0.19 / 0.009 /
  0.031 % (uncorrected, reinit every 2 steps), advection alone is second order. The Zalesak slotted
  disc loses its slot corners to reinitialisation (volume drift up to 0.45 %), covered by the
  global volume correction that the solver applies every step.
* **E2a still pool** (R = 0.5, 1 s, water): spurious velocity 1.05e-6 / 7.6e-7 sqrt(gD) (gate 1e-4),
  hydrostatic pressure error 1.75 / 0.42 % (n = 64 / 128, level 0) and 2.4 / 0.57 % (level -0.13)
  (gate 1 % at n = 128), volume exact. **Met.**
* **E2b sloshing.** Rectangular tank 1.0 x 0.5 (exact Lamb dispersion, omega = 5.31655 rad/s): frequency
  error -0.036 / -0.019 % at n = 64 / 128 (nu = 1e-6); viscous damping rate rises monotonically with
  nu (0.0002 / 0.014 / 0.075 / 0.351 1/s for nu = 1e-6 .. 1e-2). Half-full circular drum R = 0.5:
  independent Rayleigh-Ritz reference K R = 1.355727 (converged to 6 digits, harmonic basis;
  omega = 5.15746 rad/s): frequency error -0.05 / +0.02 / -0.28 % at n = 64 / 128 / 256 (tilt amp
  25 mm). **Met** (<= 1 %). A real defect was found and fixed on the way: with the interface advanced
  after the velocity update the free surface injected energy (amplitude +1 % per period, first
  order in dx); a predictor (interface advected with the old velocity defines the pressure
  topology) and a trapezoidal corrector (final level-set step with the mean of old and new
  velocity) removed it (growth rate 0.017 -> 0.0005 1/s).
* **E2c dam break** (square column a = 0.2 in a 1.0 tank, water). Front position Z = x/a at
  T = t sqrt(g/a) = 0.5, 1, 1.5, 2, 2.5: n = 64: 1.380 2.015 2.759 3.590 4.508; n = 128: 1.339 1.929
  2.710 3.576 4.488; n = 256: 1.315 1.912 2.704 3.623 4.480. 128 -> 256 changes by <= 1.8 %; all values
  below the shallow-water bound 1 + 2T. The Martin & Moyce (1952) data table could not be obtained
  from the web sources available to this session (the paper is paywalled), so no comparison with
  the experiment was made and no numbers were invented. **NOT MET: at wall impact (T ~ 2.8) the
  jet that climbs the wall (liquid sheet ~1.7 cells thick) blows up for n >= 128 (n = 64 survives).**
  Diagnosis: growth factor 2-3 per step independent of dt (CFL 0.2 -> 0.05), advection scheme
  (central / third-order / first-order upwind), time scheme (BDF2 / BDF1), incremental pressure,
  predictor, volume correction, theta_min (0.01 .. 1.0), ghost extrapolation order, wall alignment
  and Poisson convergence (9-13 iterations) are all excluded; larger viscosity delays it (nu =
  1e-4: T = 3.06, 1e-3: 3.9) and nu = 1e-2 is stable. It is a grid-scale inviscid instability of
  sheets thinner than ~2 cells, not yet understood. Open: (d) rimming flow was not run.

### E2 — free-surface investigation, final state of 2026-10-03 (gate NOT met)

Root causes found, in the order they were found:

1. **Interface coupling energy injection** (fixed): predictor/corrector level-set coupling (see above).
2. **Trailing-edge instability of liquid sheets** (fixed): a sheet translating uniformly at 2.5 m/s
   (`grid_probe --exp E2g`, exact solution, no gravity) grew a transverse velocity exponentially from
   round-off at its rear end (1e-8 -> 0.24 m/s in 0.1 s). Cause: the air-side ghost values are
   extrapolated quadratically along the grid lines from the liquid; an upwind stencil whose upstream
   node is such a ghost feeds the extrapolation back with a gain above one. Fifth-order WENO
   upwind advection (`StaggeredFlow::weno`) removes it consistently (max |v| 4e-7 at all sheet
   thicknesses 0.6 .. 3 dx and n = 64 / 128) and keeps the hold test (Re 1000, n = 128: torque
   -0.36 %).
3. **Free-surface viscous condition** (fixed, now consistent): a zero normal derivative of each
   Cartesian velocity component is not the stress-free condition; it gave an error that does not
   shrink with the grid (7.5 % of omega R for a rigidly rotating ring at n = 64 and 128, nu = 1e-3).
   The dropped links now carry the stress-free normal derivative built from the lagged velocity
   gradient (`StaggeredFlow::free_surface_flux`).
4. **Stability tricks that break consistency** (found, not used by default): clamping the air-side
   ghost values to the range of the liquid data and zeroing the advective gradient when the
   upstream node is air (`robust_surface`) keeps splashing flows alive (wall impact stable for 14 of
   15 cases, dam break stable to T = 4 at n = 64 .. 256) but leaves a 4-6 % velocity error along the
   surface of a rotating liquid ring that does not decrease with n (rigid ring, g = 0: 6.5e-2 /
   5.2e-2 / 3.9e-2 at n = 64 / 128 / 256); without them the same test gives 2.5e-3 / 5.0e-4 at
   n = 64 / 128 (order 2.3).

Validation test added: `grid_probe --exp E2r`, liquid ring in rigid rotation (exact for g = 0 only:
with gravity the circle centred at (0, g/omega^2) is not invariant under the rotation, so it is not a
steady solution; an earlier version of the test used it and its failures were a test error).

State of the four E2 items with the default (consistent) configuration, n = 64 / 128 / 256:
* (a) still pool: met (unchanged).
* (b) sloshing: rectangular tank frequency error -0.036 / -0.019 %; half-full circle
  -0.10 / -0.04 % at n = 64 / 128 but **63 % at n = 256** (the run goes wrong: the interface jitters
  at the free-surface nodes; the ring test shows the same: velocity error 1.05e-2 at n = 256).
* (c) dam break: **fails at wall impact** for all n (blow-up between T = 3 and 3.5) with the default
  scheme; before the impact the front positions agree to 1-2 % across n. With `robust_surface` it
  survives.
* (d) rimming flow (Moffatt): the film shows Rayleigh-Taylor-type ripples on the upper half that grow
  with resolution (max error 16 % at n = 128, 45 % at n = 256); the thin-film solution is not a
  stable reference without surface tension.
* wall impact of a liquid layer (`--exp E2w`, 15 cases): fails in all but one case with the default
  scheme: the jet that climbs the wall gains energy in the convective step (advective power
  +0.0006 .. +0.006 per step, exponential, `explicit_power` diagnostic during development); smaller
  time steps make it worse (it is a per-step instability, not a CFL one); eddy viscosity
  (Smagorinsky C_s up to 0.17), first-order time stepping, no predictor, other theta limits,
  bounding of newly wetted nodes and a normal (Aslam) continuation of the velocity into the air
  only delay it.

Conclusion: the free-surface solver is accurate and stable for smooth flows (hold, still pool,
sloshing at n <= 128, rigid ring), but not for jets/splashes, and there is a conflict between
consistency and robustness at the free surface that two targeted fixes did not resolve. Stop rule:
reported to the user; no building on top.

## E3s — viscous slurry in a rotating drum (first trial, 2026-10-03, gate NOT met)

Probe: `grid_probe --exp E3s --re <list> --n <list> --fill f`. Drum R = 0.5, Fr = 0.36, fill 0.3, ring start in rigid rotation.

- A flat pool started against the moving wall blows up at every n, Re and fill (also without convection, also with
  `robust_surface`): the no-slip contact line drags a sub-grid sheet. A uniform ring start is stable.
- `wall_load` summed air nodes and carried the vector-Laplacian artefact `2 nu Omega A`; replaced by
  `StaggeredFlow::rotating_wall_torque` (shear of the flow relative to the rigid wall rotation, liquid nodes only).
- Re = 20: wall torque 0.2257 / 0.2119 / 0.2171 at n = 64 / 128 / 256 (spread +-3 %, not monotone); the exact steady balance
  `torque = A g x_cm` is missed by 9 / 17 / 11 %. Re = 0.2: the rigid-rotation torque is ~0 and the traction noise
  (about nu * velocity error / dx) exceeds it.

### E3s follow-up: torque diagnostics (Re = 20, fill 0.3 unless noted)

| n | wall torque | gravity (centroid) | gravity (faces) | pressure torque | dL/dt |
|---|---|---|---|---|---|
| 64 | 0.2256 | 0.2073 | 0.2045 | 1.2e-3 | 2e-4 |
| 128 | 0.2118 | 0.1814 | 0.1796 | -4e-5 | -1e-4 |
| 256 | 0.2171 | 0.1948 | 0.1940 | -3e-5 | -4e-4 |

- The run is steady (dL/dt ~ 0) and the pressure torque is negligible, so `wall torque = A g x_cm` must hold; the wall
  torque misses it by +0.021 / +0.032 / +0.023 (fill 0.6: +0.034/+0.038 at n = 64/128; fill 0.8: +0.041/+0.051).
- Cause: the whole torque is carried by a velocity jump of ~0.2 % of the wall speed over the wall-adjacent cell
  (`tau / (nu R perimeter) * dx`), so a relative velocity error of 4e-4 already gives ~20 % torque error. The traction
  measurement is ill-conditioned for a nearly rigidly co-rotating liquid.
- Steady film thickness h(phi) at phi = 2.5, 32.5, ... deg (n = 64 / 128 / 256), e.g. phi = 242.5 deg:
  0.0389 / 0.0465 / 0.0427; 212.5 deg: 0.0387 / 0.0449 / 0.0405. Differences of +-8 % that are not monotone in n
  (n = 128 is the outlier), i.e. the film (5-10 cells thick) is not converged and the centroid torque inherits +-7 %.

## E4a — fixed disc in a still pool (2026-10-03, PASSED)

Disc r = 0.1 at (0.05, -0.12) in the drum (R = 0.5), flat pool at y = 0.2, ghost-fluid free surface, 0.2 s. Pressure force on
the disc = `-sum p (a_r - a_l, a_t - a_b) dx` over the cut cells, with the cell pressure moved to the wall point along the local
gradient (`grid_probe --exp E4a`):

| n | d/dx | F_y/(g pi r^2) - 1 (plain / corrected) | F_x | torque/(g r^3) | spurious |
|---|---|---|---|---|---|
| 64 | 11.6 | -1.10 % / -0.12 % | -5e-5 | 2e-4 | 1.8e-5 |
| 128 | 23.3 | 0.03 % / -0.03 % | -4e-5 | -3e-4 | 4e-6 |
| 256 | 46.5 | 0.15 % / 0.00 % | -8e-6 | 1e-6 | 2e-6 |

Gate (<= 0.5 %, zero torque): met. The cell-centre pressure alone is not enough at d/dx ~ 12 (1.1 %).

## E4b — prescribed translating / spinning disc (2026-10-03, PASSED)

Machinery added to `mac/staggered.rs`: `Disc` (centre, radius, velocity, spin), `StaggeredMesh::set_body_flux` (flux of the body's
surface velocity through the part of each cell boundary the body covers, exact by the divergence theorem; added to the divergence
of the projection), `FlowState` / `into_state` / `from_state` (flow carried to a rebuilt mesh every step),
`initialise_fresh` (nodes that the body uncovered get its rigid velocity, new pressure cells the neighbour mean) and
`u_star`/`v_star` (velocity before the projection, used for body loads: the projection moves nodes that sit almost on the
wall off their wall value and a wall traction `(u_b - u)/theta` amplifies that by 1/theta, giving 10 x spikes in the force when
a node happens to lie on the surface). Loads: pressure on the cut cells moved to the wall point along the local gradient, plus
the traction of the flow *relative to the rigid body velocity* (zero Dirichlet value, so the component-wise and full-stress
tractions coincide); near-wall links (theta < 0.3) use a two-node gradient.

Disc a = 0.2 in the drum R = 0.5, nu = 1 (`grid_probe --exp E4b | E4r | E4s`):

| test | n = 64 | 128 | 256 |
|---|---|---|---|
| translating disc, steady drag vs 2D Stokes cylinder-in-cylinder (Stokeslet coefficient of `psi = sin(theta)(Ar + B/r + Cr^3 + D r ln r)`), U = 0.01 | -0.03 % | +0.10 % | +0.02 % |
| spinning disc, torque vs exact Couette `4 pi nu omega a^2 b^2/(b^2-a^2)` | +0.07 % | +0.06 % | +0.04 % |
| moving (mesh rebuilt each step) vs frozen at the mean position, U = 0.1, dt 0.01 / 0.005 | -0.12 / -0.51 % | -0.20 / -0.08 % | -0.11 / -0.04 % |
| drag noise of the moving disc (rms over 20 steps) | 0.6 / 0.2 % | 0.4 / 0.3 % | 0.4 / 0.3 % |

Steady drag vs eccentric position is smooth (n = 128: 6.543 ... 6.681 for x0 = 0 ... 0.04). Gate (<= 2 % at d/dx >= 8): met
(here d/dx = 23 ... 93). Before the `u_star` fix the same test jittered by up to +-25 % depending on where the surface fell
between nodes. E1 hold tests (R2) are unchanged (torque <= 0.1 %).

## E4c — free disc coupled to the fluid (2026-10-03, PASSED for the Stokes/added-mass checks)

`mac/bodies.rs`: `BodyFlow::trial/commit/step_coupled`. The disc equation is integrated with backward Euler and solved with
a fixed-point iteration per step on the new body velocity; the position is predicted once (constant acceleration) and frozen
during the iteration (with a moving position the force jumps whenever the surface crosses a node and the map stops being
contractive). The iteration is preconditioned with the linear part of the fluid force, `K = m_a/dt + k` (added mass, Stokes
drag per unit speed), refined by secant estimates per component after the second iterate (measured K varied 82 -> 420 with
the mesh history; a fixed preconditioner then diverged exponentially), and a torque stiffness for the spin (without it the
spin iteration diverged because `dt kappa / I > 1`). 3 iterations per step.

Disc a = 0.2 in the drum (R = 0.5), force 1.0 along +x, `grid_probe --exp E4m | E4k | E4h | E4p | E4y`:

- terminal speed in Stokes flow (nu = 1, concentric drag coefficient `k = 65.4`), after 0.05 s: error
  -0.01 / +0.02 % (rho_s = 0.5 / 1, n = 64, dt 0.005), -0.2 / -0.25 % (dt 0.0025), +0.07 / +0.09 % (n = 128, dt 0.005);
  rho_s = 3: -0.5 ... -1.2 % (not yet fully relaxed).
- prescribed acceleration a0 = 3 m/s^2: the effective added mass error is `-2 % + 650 % sqrt(nu)`: 62.8 / 18.4 / 4.9 % at nu = 1e-2 / 1e-3 / 1e-4
  (n = 128), i.e. the excess is the (physical) viscous Stokes-layer force and the inviscid added mass `pi a^2 (b^2 + a^2)/(b^2 - a^2)`
  is met to ~2 % (an extrapolation, not a direct measurement).
- coupled run at nu = 1e-3: stable for rho_s = 0.5 ... 3 (added mass up to 3 x the body mass), 3 iterations per step; speed
  deficit vs the inviscid value -3.5 ... -7.8 % (viscous), independent of dt and n to 0.4 %.

## E4 summary (2026-10-03) — rigid bodies in fluid: gate MET

| item | result |
|---|---|
| (a) fixed disc buoyancy | <= 0.12 % at d/dx = 12, torque ~ 0 |
| (b) prescribed disc, Stokes drag / Couette torque | <= 0.1 % (n = 64 ... 256); moving mesh vs frozen <= 0.5 %, force noise <= 0.6 % |
| (c) free disc, Stokes terminal speed | <= 0.25 % (rho_s = 0.5, 1); added mass consistent with the inviscid value to ~2 % (extrapolated) |
| (d) neutrally buoyant / light disc at low viscosity (`--exp E4d`: nu = 1e-3, 1e-4, rho_s = 1, 1.05, 1.2, dt = 0.3 dx / sqrt(g a)) | stable, 3 iterations per step; rho_s = 1 stays at rest exactly; trajectories n = 64 vs 128 within 2 % of the drop distance |
| (e) free disc under torque, steady spin vs `T / kappa` (confined Couette) | -0.07 / -0.06 / -0.04 % at n = 64 / 128 / 256, all rho_s |

Not done: a published 2D sedimentation benchmark (Glowinski 2001 / Wan & Turek 2006): their tables are not available here and are not
reconstructed from memory. The Stokes terminal speed, Couette torque and added mass are exact references instead. Open items: several
discs and disc-disc / disc-wall closeness (E5 lubrication), the disc inside a *free surface* flow (E8; the machinery is
independent: `FlowState` carries the level set), and the cost of the mesh rebuild per step (about 5 ms at n = 128).

## E5a — disc against the drum wall: squeeze film (2026-10-03, wall normal force PASSED)

Reference: exact Stokes force on a disc translating along the line of centres inside an eccentric
cylinder (`verify::eccentric_squeeze_force`, bipolar coordinates, only the first harmonic survives,
force = Stokeslet strength `4 pi nu |d| / k`, `F = 8 pi nu u sinh 2D / (C^2 - S sinh 2D)`). Checked
against the concentric solution (e -> 0, 5 digits) and the Reynolds asymptote `3 sqrt2 pi nu u
R^1.5 h^-1.5` (ratio 1.002 at h = 1e-4). The exact/leading ratio is a function of eps = h/R_eff only:
`1 + 1.1 eps + 0.06 eps^2` to 0.3 % for eps <= 0.4 (a = 0.03 ... 0.1, R = 0.5).

Resolved solver (a = 0.05, nu = 1, U = 0.01, fixed position, dt = 0.005, t = 2; dt = 0.01 and
t < 1.6 leave a 1-5 % time-splitting / start-up error, see E5t), force over exact:

| h/dx | n=128 (a/dx 5.8) | n=256 (a/dx 11.6) |
|---|---|---|
| 16 | 0.987 | 1.001 |
| 8 | 0.982 | 0.997 |
| 4 | 0.965 | 0.972 |
| 2 | 0.895 | 0.824 |
| 1 | 0.34 | 0.29 |
| 0.5, 0.25 | 0.1, 0.01 | 0.08, 0.03 |

The resolved film is accurate to about 3 % down to a gap of 4 cells and saturates below (force
plateau near the value of a 2-cell gap; at n = 64, a/dx = 2.9, it even goes negative).
Sub-grid model (`mac/lubrication.rs`): the force along the line of centres is blended,
`(1 - s) F_grid + s c(h) u_n`, `c = 3 sqrt2 pi nu R^1.5 h^-1.5 (1 + 1.1 eps + 0.06 eps^2)`, `s` a
smoothstep from 0 at h = 8 dx to 1 at 3 dx. Blended force over exact: 0.978 ... 0.999 for all
h/dx = 16 ... 0.25 at n = 128 and 256.

E5d dynamic test: a disc (rho_s = 1.2) pushed to the wall by a constant force F = 5 (nu = 1),
coupled solver, 3 iterations per step. Without the sub-grid film the disc reaches the wall in
0.4 - 0.55 s (gap -> 0). With it the gap closes as the exact quasi-steady law: times to pass
h = 0.03, 0.02, ..., 2e-4 differ from the exact ODE by a constant 0.006 - 0.03 s (the far field of
a fluid started at rest needs about R^2/nu to build up, E5t) and the *intervals* agree to
0.984 - 0.9997 (n = 64 and 128; the first interval, 3.5 -> 2.3 dx, is the worst).

## E5 summary (2026-10-03) — near-contact lubrication: normal (squeeze) gate MET

Model generalised to signed curvatures (`mac/lubrication.rs`): `c = 3 sqrt2 pi nu R^1.5 h^-1.5
(1 + c1 eps + c2 eps^2)`, `c1 = 0.788 + 0.265 Q`, `Q = R^3 (k1^3 + k2^3)` (quartic term of the gap
profile). Second exact reference: two equal discs in unbounded fluid,
`disc_pair_squeeze_force` (`F = 8 pi nu u cosh 2xi / (2 xi cosh 2xi - sinh 2xi)`, `cosh xi = 1 + h/2a`;
reduces to `3 pi nu u (a/h)^1.5`). The model matches the wall solution to 0.3 % (a <= 0.05, eps <=
0.8; 1.8 % for a = 0.1) and the pair solution to 0.8 % (eps <= 1.6) (unit test
`squeeze_references_and_model_agree`).

Many-body coupling (`BodyFlow::trial_many`, `step_coupled_many`): per-body secant-preconditioned
fixed point; the stiff pair/wall squeeze terms are solved implicitly by Gauss-Seidel inside each
iteration (a Jacobi update has spectral radius c / (c + m/dt) -> 1). 3 - 4 iterations per step.

| test | result |
|---|---|
| E5a static disc against the wall, blended force over exact, h/dx 16 ... 0.25, n = 128, 256 | 0.978 ... 0.999 |
| E5e static two discs, blended over exact (unbounded), h <= 4 dx | 0.994 ... 1.000 |
| E5e resolved only, h = 4 dx / 2 dx / 1 dx (n = 256) | 0.98 / 0.84 / 0.22 (plateau below) |
| E5d disc pushed to the wall, intervals of the closing time 0.03 -> 2e-4 vs exact ODE | 0.984 ... 0.9997 (n = 64, 128) |
| E5f two discs pushed together, same | 0.989 ... 0.9999 |
| same without the sub-grid film | discs overlap after 0.15 - 0.55 s (gap -> 0) |

At h/a >= 1.4 the resolved two-disc force exceeds the unbounded exact value by 2 - 17 %: this is
the confinement of the drum (it does not shrink with n), not a grid error.

Open (not part of this gate): tangential and rotational lubrication. A disc sliding along the wall
(E5g, a = 0.05, n = 256, U = 0.01) is followed by the resolved flow down to a gap of about 1 cell
(`F / (2 sqrt2 pi nu U sqrt(R/h))` = 1.23 at h = 1 dx, then 0.88 at 0.5 dx and 0.57 at 0.25 dx,
where the leading term would need to be added); the missing force is `O(h/R)` of the squeeze term
and has no exact reference here. Also not done: lubrication between a disc and a free surface (E8),
a local (non-O(n^2) per disc) load evaluation for many bodies (E7/E8), non-equal-radius pair
references.

## E6 summary (non-Newtonian: Herschel-Bulkley by ALG2)

Implementation: `mac/rheology.rs` (law, d-step) and `StaggeredFlow::set_rheology`. The velocity
solves keep the constant-viscosity Gibou Helmholtz (`nu = r / 2`); the stress excess
`S = lambda - r d` (normal parts at cell centres, shear at corners, deviatoric strain, zero in air,
extrapolated beyond walls) enters as an explicit divergence; `lambda`, `d` are updated from the new
velocity, once per step or in `alg_iterations` passes per step.

E6a annular Couette (r1 = 0.25 at 1 rad/s, r2 = 0.5 fixed, K = 0.01, tau_y = 0.03 gives a plug beyond
r = 0.41), `grid_probe --exp E6a`, started from the analytic profile, n = 128, r = 0.02..0.08:

| law | ring-mean torque error | max velocity error / (Omega r1) |
|---|---|---|
| Newtonian | 0.01 % | 0.03 % |
| power law n = 0.6 | 0.02 % | 0.08 % |
| Bingham | 0.1 % | 0.3 - 0.4 % |
| Herschel-Bulkley n = 0.6 | 0.1 % | 0.3 - 0.4 % |

The result changes by < 0.5 % between r = 0.02 and 0.08. The stress in a plug is not unique (any
self-equilibrated field is admissible), so the torque constant is compared as an azimuthal ring mean.
From rest the approach to steady state is slow for large r (r = 0.32: 10 % at t = 12).

E6d tilted surface (20 degrees) in a stationary drum, `grid_probe --exp E6d` (n = 64): the surface
is held (slope 0.364 -> 0.361 over 2 s at tau_y/rho = 8) but a creep remains in the plug that falls
only as 1/r, not with the number of passes: 5e-3 m/s at r = 4, 3e-4 at r = 40 (tau_y/rho = 2.. 8),
6e-5 at r = 400; the gate is 1e-5 m/s. For partial yielding (tau_y/rho = 0.5) r = 40..400 diverges
within a few steps (the explicit stress and the implicit Laplacian are different discretisations
near cut cells, an error that scales with r), r = 4 is stable and creeps at 4e-4 m/s.
Gate (d) is NOT met. Tried: more passes, air masking, wall extrapolation of the stress, deviatoric
strain, an exact same-operator cancellation (worse).

### E6 follow-up: variational stress operator (2026-10-04)

Diagnosis of the plug creep. A static test without a free surface (`grid_probe --exp E6s`: a drum
filled with a Bingham fluid, driven only by the swirl force `eps (-y, x)`, exact solution `u = 0`)
showed that the creep (3e-3 m/s at r = 1, falling only as 1/r) is not a free-surface effect and
that more ALG2 passes made it worse: the strain operator `D(u)` and the stress divergence of the
original scheme are not adjoint, so the multiplier winds up in the null space of the divergence.

Fix: `mac/variational.rs`. The deviatoric strain at cell centres and corners (including
no-slip ghost velocities through the wall point, corners that touch solid cells keep the cells
that exist) is one sparse linear map `B`; the stress force is exactly `-B^T W g tau`, the implicit
operator `B^T W g B` is solved by CG preconditioned with the scalar Gibou multigrid. Air needs no
condition (points that touch air carry no stress). Used alone it is first order at walls (Stokes
annulus velocity error 3.8 % / 1.2 % of `omega r1` at n = 64 / 128, the Gibou scheme is ~0.1 %), so
the default is the **hybrid** (`set_rheology_hybrid`): the Newtonian (plastic) viscosity stays in
the Gibou solve, only the yield stress goes through the variational augmented Lagrangian.

Static swirl hold (n = 64, tau_y = 1, eps = 2, max speed after 1 s):

| scheme | r | passes | max speed |
|---|---|---|---|
| original | 10 | 1 | 8e-7 (lagged, no hold of a real load) |
| original | 10 | 20 | 1e-4..3e-4 |
| variational | 100 | 8 | 4.7e-5 |
| variational | 1000 | 8 | 5e-6 |
| hybrid | 1000 | 6 | 2e-6 |

Tilted surface (20 degrees, n = 64, 1 s, hybrid, r = 1000, 5 passes): tau_y/rho = 8: slope 0.3629
(initial 0.3640), bulk speed 3.3e-5 m/s (was 5e-3) with interface-layer spikes of 1e-2 m/s for a few
tenths of a second; tau_y/rho = 0.5 (partial yield) slumps steadily (slope 0.359 -> 0.295) and is
stable at r = 1000 (the original scheme diverged for r >= 40). Bingham annulus (hybrid, n = 128):
velocity error 0.7 %, speed in the plug 1.5e-6 m/s (was 2e-4).

Gate (d) (<= 1e-5 m/s with a free surface) is still not met: 3.3e-5 bulk, interface spikes. Cost is
the open issue: 5-20 passes of 200-400 CG iterations per step (the preconditioner ignores the yield
operator), minutes per simulated second at n = 64.

## E7 — ball size versus grid (2026-10-04): decision D2 input

E7a, one disc (a = 0.05) translating in the drum, Stokes, 8 sub-cell positions per resolution, drag over the exact
eccentric solution (`grid_probe --exp E7a`):

| d/dx | 32 | 16 | 8 | 4 | 2 | 1 |
|---|---|---|---|---|---|---|
| mean error | -0.5 % | -0.6 % | -1.2 % | -3.6 % | -8 % | -26 % |
| worst case over positions | -0.6 % | -1.2 % | -1.8 % | -6.1 % | -15 % | -36 % |

A single resolved disc stays within 3 % down to d/dx = 8.

E7b, 20 discs (hexagonal cluster, gap 0.4 a, rho_s = 2) settling in Stokes flow, same geometry at every grid
(`grid_probe --exp E7b --n N`), mean speed at t = 0.1 relative to the d/dx = 32 run (0.07095):

| d/dx | 32 | 16 | 8 | 4 | 2 |
|---|---|---|---|---|---|
| speed error | reference | -2.0 % | -9.9 % | -24.9 % | diverged |

Finding: the E5 normal-lubrication blend removed the whole resolved normal force on every linked disc, including
the common-mode drag, so a dense cluster fell 7 x too fast. Fixed: the sub-grid force is now additive,
`-(1 - resolved_ratio(gap/dx)) c u_rel,n`, with the resolved fraction tabulated from E5 (`lubrication::resolved_ratio`).
The static pair/wall E5 accuracy drops from 0.98-1.00 to about 0.90-1.10 at gaps of 1-2 cells (the resolved fraction
there scatters with disc size and position). Not done: tangential lubrication, the unresolved (porosity + drag) model.

D2 so far: resolved coupling is within 3 % only for d/dx >= 8 (isolated) and d/dx >= 16 (dense, gap < a/2);
below d/dx ~ 4 it fails. n = 440 with 20 discs costs about 700 s for 0.1 s.

### E7c — k-consistency of the resolved coupling (2026-10-04)

Same area, same total force, same cluster footprint, different coarse-graining: 28 discs of
radius a = 0.04 (d/dx = 16, n = 220) against 7 discs of radius 2a = 0.08 (d/dx = 32, n = 220),
hex pitch 2.4 radius, nu = 1, rho_s = 2, mean settling speed at t = 0.1:

| case | mean vy | note |
|---|---|---|
| 28 x a   | -0.0558 | d/dx = 16: about -2 % from the d/dx = 32 value (E7b) |
| 7 x 2a   | -0.0579 | |

The k = 2 coarse-grained cluster settles 3.7 % faster raw, about 1.5 % after correcting the
fine case for its own grid error. A dense cluster (gap 0.4 radius) behaves as a nearly
impermeable blob, so permeability (which scales with the radius squared) matters little here;
a dilute or fluid-dominated charge would show the radius-squared dependence and is not covered.
Consequence for D2: coarse-graining raises d/dx, so the resolved coupling needs no unresolved
(porosity + drag) model as long as k is chosen to keep d/dx >= 8 (isolated) or 16 (dense); k = 2 is
k-consistent to about 2-4 % in the dense case. Open: k = 4, dilute charges, rotation/tangential lubrication.

### E8a — discs in a filled rotating drum (2026-10-04)

D2 fixed by the user: resolved coupling only, k keeps d/dx >= 8 (isolated) .. 16 (dense).
E8a is the first DEM-like system on the grid: 12 discs (a = 0.04, rho_s = 2) in a fully filled
drum (R = 0.5, omega = 2 rad/s, nu = 0.1, g = 50, lubricated, no contact springs, no friction),
`grid_probe --exp E8a`. Power is read as the gravity moment of the discs (steady state: equals
omega x torque); the wall torque of the fluid is too noisy (0.06 .. 0.64 within 0.1 s, the
E3 conditioning problem again).

Infrastructure: rotating drum wall and drum-only wall torque in `BodyFlow`; many-body coupling
with Anderson acceleration (depth 5, scaled velocities) and a frozen secant preconditioner after 3
iterations (a preconditioner that changes between iterations defeats the acceleration; with it
the fixed point needs 9-26 iterations instead of hitting the limit of 40 and diverging at
t = 0.33 s). The fluid force response to a body velocity is smooth and linear (perturbations 1e-1
.. 1e-5 agree within 1 %). Cost: mesh built once per coupling step (positions are fixed during the
iteration), distance early-outs in the aperture/centroid/cell-fraction construction, pruned disc
sdf, linear tolerances 1e-6 instead of 1e-9 inside the coupling (moment changes 0.02-0.2 %):
n = 128, 12 discs, 0.5 s simulated: ~390 s -> 26 s (15x). The grid size matters: n = 110 has only
two multigrid levels and a 3x slower projection than n = 128.

Result: the system does not reach a steady state. The cluster rides around the drum (angle > 2.9 rad,
period ~3.5 s) and the gravity moment oscillates between -1.04 and +0.16. Up to t = 1.2 s the
moment agrees to 1 % for n = 128/192/256 and dt, 0.005/0.0025 (-1.035 .. -1.043); later the phase of
the orbit differs between runs (t = 2.4: -0.35, -0.38, -0.16, -0.15) and three of the four 6 s runs
(n = 192, 256, n = 128 with dt/2) went to NaN at t = 4.8-5.4, where the cluster reaches the wall and
close pair gaps. Only n = 128, dt = 0.005 finished.

Verdict: E8a does not give a +-3 % power statement. Missing for a meaningful full-system test:
contact (normal spring, Coulomb friction), tangential/rotational lubrication, a free surface
(the E2/E3 gates are open), and a configuration that settles to a steady state. Earlier 2 s runs
(n = 110/165/220, before the coupling fixes) gave a mean moment over t = 1..2 of
-0.968/-0.965/-0.950 (transient, not a steady state).
