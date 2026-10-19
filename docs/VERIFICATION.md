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


## E8b — solid contact and a jammed pile (2026-10-04)

Added `mac/contact.rs`: one-sided implicit normal constraint (allowed approach `gap/dt`, Baumgarte
fraction `beta` of an existing overlap, stiffness 100 `m_eff/dt`, solved with the lubrication links by
projected Gauss-Seidel) and regularised Coulomb friction (explicit in the coupling iteration, its
stiffness added to the preconditioner). `grid_probe --exp E8a --mu <friction> [--gap0 --beta --rho]`.

Results (n = 128, d/dx = 9.3, nu = 0.1, omega = 2, 12 discs, dt = 0.005):
- rho_s = 2 (light discs orbit with the fluid): with contact, 6 s runs finish without NaN (before: NaN
  at t ~ 5 when the cluster reaches the wall). Mean gravity moment over t = 3..6: -0.3736 (mu = 0) and
  -0.3739 (mu = 0.5); the cluster still orbits (period ~3.5 s), no steady state.
- rho_s = 8 (the pile sits at the bottom, angle ~0.55 rad): NaN at t ~ 0.8-1.2 s for mu = 0 and 0.5,
  for every contact variant tried (soft dashpot: wall penetration under the stack weight; stiff
  constraint: the many-body fixed point stops converging in one step, 40 iterations, residual 1.8,
  then overlap). With friction the discs also spun up to 10-40 rad/s (explicit friction gain > 1)
  until the friction stiffness went into the preconditioner; mu = 0 fails anyway, so the cause is
  the stiff jammed contact network inside the fixed point, not friction.

Verdict: stop rule reached (two targeted fixes). The contact model works for a free cluster, but the
fixed-point coupling cannot handle a jammed pile. No +-3 % statement yet.


## E8c — Newton coupling of all bodies and tangential lubrication (2026-10-04)

`BodyFlow::step_newton` (+ `mac/newton.rs`): one Newton iteration on all body velocities `(ux, uy,
r omega)`; the fluid load is linearised with a dense, symmetrised Jacobian (finite differences every
8 steps, Broyden updates between) and the dense system with the squeeze-film links, the one-sided
contacts (active set), Coulomb friction (stick/slip set) and the new tangential film shear
(`lubrication::tangential_coefficient`, `nu pi sqrt(2 r / h)`, weighted like the normal film) is solved
exactly; step cap 0.3 m/s per iteration. Probe: `grid_probe --exp E8a --newton`.

- Light discs (rho_s = 2, 1.2 s): Newton and fixed-point coupling agree (gravity moment -1.0347 vs
  -1.0360, 0.1 %), cost about equal.
- Heavy pile (rho_s = 8, mu = 0 and 0.5): still NaN. Time of failure: 1.0 s (no tangential film) ->
  1.14 s (with it) -> ~1.1 s (symmetrised Jacobian, step cap). A finite-difference Jacobian at every step
  (1 h run) fails at the same time, so staleness of K is not the cause; the iteration stops
  converging (residual 1-5 for 20 iterations) while the pile slides up the wall with gaps of
  2e-5..5e-5 m = 0.003 dx, i.e. far below the grid scale, where the model is a sum of lubrication
  laws (sub-grid) and one-sided contacts, a non-smooth system.
- Before the tangential film was added the discs spun at 10-15 rad/s even with mu = 0; this is gone.

Verdict: stop rule reached again. Open causes: non-smooth active-set cycling, accuracy of the
tangential/rotational film at gaps << dx (only the leading Couette term, not calibrated against the
exact bipolar solution as the normal film was), the unresolved fluid in the pores of a pile.


## E8d — contacts split off the coupling, hydrodynamic radius (2026-10-04)

Two changes, probe flags `--split` and `--shrink <cells>` (`grid_probe --exp E8a`):
- `BodyFlow::split_contact`: contacts are no longer solved inside the nonlinear fluid coupling; after
  the coupled step a Moreau-type velocity projection (`contact::project`, sequential impulses, normal
  one-sided + Coulomb friction, 200 sweeps) acts on the fluid-step velocities, and the centres are
  `x0 + dt v`. Alone this did NOT cure the heavy pile (NaN at 0.8 s): the first trial of the failing
  step already blew up inside the fluid solver (`max|u|` 4e5, body speeds 0.17 m/s, pair gap 5e-7 m),
  i.e. the cause was the fluid mesh, not the contact non-smoothness.
- `BodyFlow::hydro_shrink`: the discs given to the fluid mesh/flux have radius `r - shrink * dx`
  (hydrodynamic radius), so sub-grid gaps never form sliver cells; lubrication and contact keep the
  true radius.

Result (rho_s = 8, mu = 0.5, n = 110, d/dx = 8, nu = 0.1, omega = 2, 12 discs, shrink = 0.15): first
run of the heavy pile to t = 2 s without NaN; the pile slides up to ~0.93 rad, gravity moment
settles near -6.8 (mean over the 2nd half -6.42). No sensitivity study yet (shrink, grid, dt), so
no +-3 % statement; the shrink biases the displaced volume by ~2 shrink dx / d.

### E8d sensitivity (rho_s = 8, mu = 0.5, split contacts, t = 2 s)

| variant | gravity moment, mean t = 1..2 | at t = 2 | pile angle at t = 2 |
|---|---|---|---|
| n = 110, shrink 0.10 | -5.98 | -6.19 | 0.82 |
| n = 110, shrink 0.15 | -6.42 | -6.82 | 0.93 |
| n = 110, shrink 0.25 | -5.93 | -6.47 | 0.87 |
| n = 110, shrink 0.15, dt = 0.0025 | -6.30 | -6.82 | 0.95 |
| n = 128, shrink 0.15 | -5.39 | -5.64 | 0.74 |

All five runs finish without NaN, but the spread is about +-9 % (mean) and the pile is still
climbing the wall (angle not steady), so a window mean of a transient is not a +-3 % metric. The
largest shift comes from the grid (n 110 -> 128: -16 %), i.e. the sliding pile is not grid-converged
yet; shrink and dt matter less (<= 7 %). Gate NOT met.

### E8d static benchmark: 12 heavy discs settling in a still drum (omega = 0)

`grid_probe --exp E8a --rho 8 --mu 0.5 --split --omega 0 --t 2 --shrink 0.15`, n = 96 / 110 / 128:
no NaN, 8-17 coupling iterations typical (peaks 27-40), the pile settles to the drum bottom (angle
|0.05| rad or less). Residual gravity moment (window mean t = 1..2): -0.036 / +0.168 / +0.090, i.e.
|moment| <= 0.2 against ~6 in the rotating case (<= 3 %), but it does not converge to zero: the
pile keeps rocking by a few hundredths of a radian over 2 s (contact projection + sub-grid film
leave a slow creep). Verdict: stable, physically plausible, not a quantitative gate yet.

### Open after E8d
- Tangential/rotational film (step 2) is still the uncalibrated Couette leading term. A calibration
  needs a trusted reference (exact bipolar-coordinate solution of a cylinder translating/rotating
  next to a wall); it was not done unattended because the 2D lubrication coefficient is subtle
  (the pressure-driven part changes the net shear force) and must not be guessed.
- The dynamic pile is not grid-converged (n 110 -> 128: -16 %); find the source (sliding speed vs
  grid) before any +-3 % claim.

### E8d grid study at constant physical shrink (1.5e-3 m), rho_s = 8, mu = 0.5, t = 2 s

| n (d/dx) | shrink (cells) | gravity moment, mean t = 1..2 | at t = 2 | pile angle at t = 2 |
|---|---|---|---|---|
| 96 (7.0) | 0.131 | -4.73 | -4.71 | 0.60 |
| 110 (8.0) | 0.150 | -6.42 | -6.82 | 0.93 |
| 128 (9.3) | 0.175 | -5.72 | -6.05 | 0.80 |
| 160 (11.6) | 0.218 | -6.71 | -7.71 | 1.15 |

Not converged (spread -4.7 .. -6.7, not monotone at n = 110): the pile climbs the wall faster
the finer the grid. Making the shrink a fixed physical length did not remove the trend (n = 128:
-5.39 -> -5.72), so the shrink bias is not the main cause. What differs with n is how much of the
wall/pair gap is carried by the grid and how much by the sub-grid film models (blend near dx), and
the tangential film is the uncalibrated part: the sliding resistance of the pile, hence its climb
speed, is set there. Conclusion: the tangential film calibration (needs a trusted reference) is the
blocking item for any +-3 % statement on the dynamic pile; further grid or dt sweeps before it would
only re-measure the same scatter.


## E8e — calibrated tangential film, Newton path with split contacts (2026-10-04)

Finding: the tangential film existed only in the Newton path (`newton.rs`); every earlier
`step_coupled_many` pile run (E8b, E8d tables above) had NO tangential film, and its grid scatter
was the missing sliding/spin resistance.

Derivation (2D lubrication, film `h = h0 + x^2/(2 r)`, `r = 1/(k_i + k_j)`, unbounded fluid so
`int p' dx = 0`, steady in the frame of the contact point): flux `Q = (2/3)(V_i + V_j) h0`, and with
`K = 2 pi nu sqrt(2 r / h0)` the tangential force on surface `i` is `K (phi_i (V_i + V_j) - V_i)`
(`phi_i = k_i / (k_i + k_j)`), the force on `j` is minus that (no net film force), and the torque on
`i` is `-K r_i V_i` (pressure acts through the centre). `V` are surface speeds relative to the
contact point `c = phi_i u_i + phi_j u_j` (tangential). For a cylinder next to a flat wall this gives
`F = 2 pi mu sqrt(2a/h) U = 2^(3/2) pi mu U sqrt(a/h)`, the leading term of the exact bipolar
solution (Jeffrey and Onishi 1981), and a torque that depends only on the spin. The resistance matrix
is symmetric (checked by hand, e.g. `dF_i/dw_j = dT_j/du_i = -K phi_i r_j`). The previous code used
`nu pi sqrt(2r/h)` on the slip speed only: half the exact scale and wrong coupling.
The reciprocity of the 4x4 block was checked algebraically, not by a unit test; the 2D
exact-solution comparison is the asymptote above, not a numerical evaluation of the full series.

Implementation: `lubrication::tangential_coefficient` (2 pi), `Link::phi`, film block in
`newton::solve_step`; `step_newton` now also supports `split_contact` and `hydro_shrink`.

Heavy pile, `--newton --split --mu 0.5 --rho 8`, constant physical shrink 1.5e-3 m, t = 2 s,
gravity moment (iterations per step 4-10, no NaN):

| t | n = 96 | n = 110 | n = 128 |
|---|---|---|---|
| 0.6 | -5.48 | -5.47 | -5.30 |
| 1.0 | -7.77 | -7.78 | -7.77 |
| 1.2 | -8.01 | -8.08 | -8.21 |
| 1.4 | -7.81 | -7.84 | -8.10 |
| 2.0 | -4.82 | -4.86 | -5.12 |

Compared with E8d (spread -4.7..-6.7 on the window mean, trajectories not alike) the grid
dependence is now small: up to t = 1.2 the three grids agree within 3.3 % (peak gravity moment
-8.01 / -8.08 / -8.21, spread 2.4 %), later up to ~6 % as the pile goes over the top (angle > 2 rad,
a cascade). The window mean of E8d is not a valid metric here (the pile is a transient); the peak
gravity moment and the trajectory up to the peak are. Gate: MET for the pre-peak window and the peak
on n = 96..128 (d/dx = 7..9.3) for rho_s = 8, mu = 0.5, nu = 0.1; NOT yet shown for dt, shrink,
other mu, rho_s = 2, water / 50 Pa.s / Bingham, or n > 128.

### E8e sensitivity (2026-10-04)

Base case n = 110, dt = 0.005, shrink 1.5e-3 m, mu = 0.5, rho_s = 8, nu = 0.1, `--newton --split`;
one knob changed per run. Gravity moment:

| run | t=0.6 | t=1.0 | t=1.2 (peak) | t=1.6 | t=2.0 |
|---|---|---|---|---|---|
| base (n=110) | -5.47 | -7.78 | -8.08 | -7.31 | -4.86 |
| dt = 0.0025 | -5.42 | -7.81 | -8.15 | -7.10 | -4.45 |
| dt = 0.01 | -5.59 | -7.86 | -8.00 | -6.59 | -4.00 |
| shrink x0.5 | -5.38 | -7.76 | -8.09 | -7.31 | -5.19 |
| shrink x2 | -5.50 | -7.88 | -8.26 | -7.38 | -5.03 |
| mu = 0.3 | -5.47 | -7.77 | -8.06 | -7.35 | -5.06 |
| mu = 0.8 | -5.47 | -7.78 | -8.07 | -7.35 | -5.06 |
| n = 160 (d/dx 12.4) | -5.37 | -7.82 | -8.14 | -7.43 | -4.76 |

Peak gravity moment over all runs (incl. n = 96 and 128 from above): -7.997 .. -8.261, i.e. within
+-1.7 % of the mean (3.3 % full range); the largest single deviation from base is shrink x2 (+2.3 %).
Up to t = 1.0 all runs lie within 2 %. Friction coefficient barely matters in this scenario.
After the peak the trajectory is sensitive to dt: dt = 0.01 is 10 % off at t = 1.6 (and 18 % at
t = 2.0), dt <= 0.005 is within ~5 % (t = 1.6) and the spread at t = 2.0 (-4.45 .. -5.19) is the
cascade regime. Gate on the peak metric: MET for dt in [0.0025, 0.01], shrink x0.5..x2, mu 0.3..0.8,
n = 96..160 at rho_s = 8, nu = 0.1. Still open: rho_s = 2, other viscosities, Bingham, post-peak
trajectory (use dt <= 0.005).

### E8e fluid sweep (2026-10-05)

`--newton --split`, dt = 0.005, shrink 1.5e-3 m, mu = 0.5; n = 96 (d/dx 7) vs n = 128 (d/dx 9.3).
Peak gravity moment (n=96 / n=128):

| rho_s | nu | peak | difference |
|---|---|---|---|
| 2 | 0.1 | -1.052 / -1.059 | 0.6 % |
| 8 | 1.0 (slurry-like) | -7.10 / -7.22 | 1.7 % |
| 8 | 0.01 (water-like) | NaN / NaN | both runs fail by t = 0.4 |

At rho_s = 2 and nu = 1 the peak is grid-converged within the +-3 % gate (trajectories agree to
~2 % up to t = 1.2; after the peak they diverge as in the nu = 0.1 case). At nu = 0.01 the run goes
NaN at both resolutions (vmax ~ 2.6 m/s, persistent pair overlap -1.9e-3 m with the split contact
list reporting 0 coupled contacts): the water-like regime is NOT met and joins E2 (high-Re /
free-surface jet robustness) as an open item. Bingham not yet run.

#### E8e low-viscosity diagnosis (2026-10-05)

Setup as the fluid sweep (`--newton --split`, rho_s = 8, mu = 0.5, n = 96, shrink 0.131 dx). Boundary and
switches at nu = 0.01 (water-like):

| variant | result |
|---|---|
| nu = 0.03 | runs, gravity moment -7.72 at t = 1.0 |
| nu = 0.01 | NaN at t ~ 0.3-0.4 |
| nu = 0.01, dt = 0.0025 | NaN earlier (t < 0.2) |
| nu = 0.01, mu = 0 or 0.1 | NaN at t ~ 0.3 (friction is not the cause) |
| nu = 0.01, central advection (`--central`) | NaN (the WENO scheme is not the cause) |
| nu = 0.01, Jacobian refreshed every step | NaN at t ~ 0.35 (Jacobian staleness is not the cause) |

Mechanism seen in the per-disc dump: the pile slides along the bottom of the drum at about the wall
speed (1 m/s) with wall gaps of order -1e-4 m, disc spins chatter in sign from step to step (|w| up to
6 rad/s, then growing to 1e2 rad/s), and the Newton solve stops converging (residual ~1 at 20
iterations) a few steps before the NaN. At nu = 0.1 the film damps the spins and the same setup is
stable. So the water-like regime is not met; the likely cause is the non-smooth (sub-grid gap, 1/sqrt(h)
film) response combined with weak viscous damping, which a linear fluid Jacobian cannot represent.

#### E8e Bingham paste (2026-10-05)

`BodyFlow::rheology` (ALG2 Bingham, plastic viscosity nu = 0.1, augmentation r = 2 nu; the lubrication film keeps
the plastic viscosity) and probe `--tau` (kinematic yield stress). Setup: `--newton --split`, rho_s = 8, mu = 0.5,
dt = 0.005, shrink 1.5e-3 m, t = 1.2 s, n = 96 vs 128. Gravity moment at the end / peak of the window:

| tau_y/rho (m^2/s^2) | n = 96 | n = 128 | difference |
|---|---|---|---|
| 2.5 (Bn ~ 1) | -7.558 (still rising) | -7.681 | 1.6 % |
| 10 (Bn ~ 4) | -6.491 (peak at t = 0.96) | -6.496 | 0.1 % |

The whole trajectory agrees within about 2 % (tau = 10: 1.8 % at t = 1.2). The gate is met for a paste in the
viscous range (nu >= 0.1); the water-like limit stays open (see above). Not tested: Herschel-Bulkley exponent
n != 1, hybrid/variational mode, yield stress with the free surface.

Further switches at nu = 0.01 (n = 96, probe `--nolub` disables the sub-grid film): without the film the run fails
earlier (NaN before t = 0.16, as expected: discs overlap the wall); with a 0.4 dx hydrodynamic shrink instead of 0.131
it still fails (wall torque 2e2 at t = 0.24, NaN by 0.32). So neither the film nor sliver cells explain the failure;
the water-like regime is parked together with E2 (stop rule: more than two attempts without meeting the gate).

### E3s resolution sweep (2026-10-05, Re = 20, fill 0.3, Fr = 0.36)

| n | wall torque | gravity moment (centroid) |
|---|---|---|
| 64 | 0.2257 | 0.2073 |
| 80 | 0.2220 | 0.2041 |
| 96 | 0.2181 | 0.1997 |
| 112 | 0.2090 | 0.1936 |
| 128 | 0.2120 | 0.1815 |
| 160 | 0.2107 | 0.1888 |
| 192 | 0.2146 | 0.1912 |
| 256 | 0.2171 | 0.1948 |

Settling 8 instead of 2 revolutions changes nothing (n = 64, 128: < 0.1 %), so the state is steady. The gravity moment
spans 0.1815-0.2073 (+-7 %), falling to n = 128 and rising again; the +-3 % gate is not met for the pure-slurry drum.
The steady film is only 2-5 cells thick at n <= 128 and not converged. Pure-slurry E3 is closed as not met; the
full-system E8 sweep (nu >= 0.1) is the accepted evidence for the power metric.

### Grid track status summary (2026-10-05)

| item | status |
|---|---|
| E0, E1, E4, E5, E7 | gates met |
| E6 (yield stress) | Couette gate met; plug creep 3e-5 against the 1e-5 gate (hybrid mode), closed as practically sufficient; fixing it needs a yield-aware preconditioner for the variational solve |
| E8 full system, nu >= 0.1 (slurry, Bingham paste) | peak gravity moment within 2 % between n = 96 and n = 128 (gate +-3 % met) |
| E8, nu = 0.01 (water-like) | NaN; joins E2 (high-Re / free-surface robustness) |
| E2, E3 (pure slurry) | gates not met (see sections above) |

### Simulation back on PBF (2026-10-05)

A screenshot of the web app (63 mm drum, 100 rpm, slurry 0.35, 244 balls) showed balls outside the drum and balls
jumping inside the charge. `examples/escape_probe` reproduces it headless on the DFSPH-backed `Simulation`
(c52d3e7): single balls at rest receive 6 m/s in one sub-step from the fluid impulse and leave the wall by up to
88 mm at t = 0.5-0.75 s; with the slurry off, or on the PBF-backed `Simulation` (c52d3e7~1), the worst wall overlap is
5e-5 m and no ball leaves the drum over 6 s. DFSPH coupling remains unfit for the dense charge ("bursts at res 60"), so
`Simulation` uses PBF again (c52d3e7 reverted); the DFSPH code stays in `fluid/` for the grid-track comparison.

## PBF vs grid benchmark, first pass (2026-10-05)

Case: 10 heavy discs (a = 0.04, rho_s/rho_f = 8) in a drum of radius 0.5 filled with slurry
(nu = 0.1, omega = 1, g = 9.81). Metric: net gravity moment of the discs about the axis (grid
units, rho_f = 1; the PBF value from `pbf_bench` is divided by rho_f = 1000).
Grid: `grid_probe --exp E8a --n 96 --nu 0.1 --omega 1 --rho 8 --m 10 --g 9.81 --mu 0.5 --split
--shrink 0.15 --newton --t 4`. PBF: `pbf_bench --t 4 [--res --substeps]` (slurry fill 0.9).

| solver / setting | peak gravity moment |
|---|---|
| grid n=96 | -1.30 (t = 2.8, angle 1.6 rad) |
| PBF res 40, substeps 8 (app default) | -0.926 |
| PBF res 60 / 80 | -0.782 / -0.741 |
| PBF substeps 4 / 16 | -0.758 / -1.000 |

Findings:
- PBF depends strongly on the cost knobs: -0.74 .. -1.00 (about +-15 %) against the 3 % gate.
- Absolute PBF-vs-grid difference (about 30 %) is not yet meaningful: initial conditions differ
  (grid starts in solid-body rotation with a hex block of discs, PBF spins up from rest with a
  random packing, and keeps a 10 % air gap).
- Next if continued: make both start identically (shared initial disc positions and spin-up),
  then attribute the PBF error to drag, lubrication or viscosity terms.

### PBF vs grid, identical start (2026-10-05)

`Simulation::with_initial_balls` + `pbf_bench --match 1 --m 10` start PBF from the grid's hex block of
10 discs with balls and slurry in solid-body rotation (the 10 % air gap remains). Same case as above
(nu 0.1, rho_s 8, omega 1). Peak gravity moment in grid units (grid n=96: -1.30):

| PBF setting | peak | 2nd-half mean |
|---|---|---|
| res 40 / substeps 8 | -0.637 | -0.615 |
| res 80 / substeps 8 | -0.457 | -0.367 |
| res 40 / substeps 16 | -0.966 | -0.837 |

With a shared start PBF is 26-65 % below the grid and the cost-knob spread grows to about +-35 %, so the
earlier +-15 % was optimistic. Next: attribute to drag / lubrication / viscosity terms.

### PBF drag attribution: isolated disc (2026-10-05)

`pbf_bench --match 1 --m 1 --y0 0 --omega 0.001` drops one disc (rho_s 8, a = 0.04, nu 0.1) from the
centre of the slurry-filled drum; the grid's exact concentric Stokes result (`stokes_annulus_drag`)
gives a terminal speed of 0.423 (grid units).

| PBF setting | settling speed before the wall |
|---|---|
| res 40 / substeps 8 | 0.8-0.9 |
| res 80 / substeps 8 | 0.9-1.4 |

PBF settles 2-3x faster than the exact Stokes drag allows and faster still at finer resolution, so the
**ball drag term** (step 6.5 of `pbf.rs`: relaxation to an entrained fluid mass with free-space
`tau_stokes = rho r^2 / 4 mu`) is the main source of the resolution dependence. A 2D disc has no
free-space Stokes drag (it depends on the confinement through a logarithm), so that constant has no
physical basis. Fix direction: replace it by a drag calibrated on the confined (annulus/pair) solutions
the grid track provides, then re-run the benchmark above. Lubrication and viscosity are not yet isolated.

#### Correction: the drag rate is not the lever (2026-10-05)

Replacing `tau_stokes` by a confined (Kuwabara) rate and scaling it x1 / x3 changed the isolated-disc
settling speed by < 5 %: the closure is already saturated (`beta` ~ 1), so the ball simply follows its
entrained fluid. A viscosity sweep shows the resistance is not viscous at all in PBF:

| nu | exact Stokes terminal speed | PBF plateau speed |
|---|---|---|
| 0.1 | 0.423 | 0.8-0.9 |
| 0.3 | 0.141 | 0.63-0.69 |
| 1 | 0.042 | 0.54-0.58 |
| 3 | 0.014 | 0.35-0.43 (no clamp hits) |

A 30x error at nu = 3 means the viscous momentum exchange between the entrained blob and the resolved
fluid around it is far too weak (the fluid below the disc is displaced almost freely). The drag
replacement was reverted. Fixing this needs the fluid response around a ball (ball-surface no-slip in the
viscosity solve / a resolved-flow boundary), not the relaxation constant.

### PBF ball no-slip boundary (2026-10-05)

The viscous exchange between a ball and the resolved fluid was far too weak (previous section), so the
Stokes part of step 6.5 is gone. Instead step 7 treats each ball's circumference (`dx = h/2` ghost
points carrying the rigid-body surface velocity) as a Dirichlet boundary of the Morris viscosity
Laplacian, with the **balls' velocity and spin as extra unknowns of the same implicit CG solve**
(`solve_viscosity_with_balls`: SPD, Jacobi-preconditioned, exact linear/angular momentum exchange).
Treating the ball explicitly was unstable at >= 10 Pa*s (v_rms 0.88 m/s in a still drum, 14 % clamp
hits); starting the solve from the ball velocity *after* this sub-step's earlier impulses avoids
double counting with form drag and overlap pushes. Form drag (inertial) stays in step 6.5.
`BALL_BC_GAIN = 0.5` compensates the single ghost layer; the isolated disc is insensitive to it
between 0.5 and 10, the dense pile varies least with resolution at 0.25-0.5.

Isolated disc (`--match 1 --m 1 --y0 0 --omega 0.001`), settling speed, exact concentric Stokes in
brackets: nu 0.1 0.37-0.39 (0.423), nu 3 0.02-0.05 (0.014); res 40 and 80 agree within 2 %
(before: 0.8-1.4 and 0.35-0.43, i.e. 2-3x and 30x too fast).

Dense pile (`--match 1 --m 10`, nu 0.1, grid n=96: -1.30), peak gravity moment:

| setting | before | after |
|---|---|---|
| res 40 / 60 / 80, substeps 8 | -0.64 / - / -0.46 | -1.04 / -1.03 / -1.04 |
| res 40, substeps 4 / 8 / 12 / 16 | - / -0.64 / - / -0.97 | -0.99 / -1.04 / -1.16 / -1.18 |
| nu 1, res 40 / 80 / substeps 16 | - | -1.19 / -1.15 / -1.16 |

Resolution dependence is gone (1 %), and at nu = 1 all knobs agree within 3.5 %. At nu = 0.1 the
sub-step count still moves the peak by +-9 % and the absolute level is ~20 % below the grid, so the
+-3 % goal is **not yet met** at low viscosity; candidate causes: no sub-grid lubrication/contact
film between discs, the 10 % air gap, and the position-based sub-step dependence. k (coarse-graining)
is not re-measured (the matched start does not scale with k); a coarse ball now sees its own radius
in the boundary, i.e. drag grows only logarithmically with k, consistent with E7c for dense charges.
Cost: +45 % per frame (realtime 6.0 -> 10.5, balanced 18 -> 28, accuracy 54 -> 77 ms, same session).
Bingham not re-run. The old test of k-independence of the relaxation fraction was removed with the
Stokes closure; all other tests pass.

#### Sub-step dependence is a splitting error, not iterations (2026-10-05)

Isolated disc settling speed vs sub-steps 4 / 8 / 16 (exact Stokes in brackets): nu 0.03 0.42 / 0.52 /
0.53; nu 0.1 0.47 / 0.40 / 0.36 (0.423); nu 3 0.075 / 0.040 / 0.031 (0.014). Raising `--iters` (3 -> 12)
or `--diters` (2 -> 8) changes nothing, so the density/contact solver is not the cause. The error shrinks
roughly quadratically with dt: the viscous step acts on the velocity after the position-based pressure
projection, so within a sub-step the displaced flow is not resisted (operator-splitting error). For the
dense pile at nu 0.1 the peak converges at 12-16 sub-steps (-1.16 / -1.18, i.e. ~9 % below the grid's
-1.30); 4-8 sub-steps are not converged. Consequence: +-3 % across sub-steps needs >= 12 sub-steps at
nu 0.1 (cost), or a splitting fix (Strang half-steps) that is not yet tried. Moving the solve to the start
of the step does not help (the scheme is cyclic). `pbf_bench` gained `--iters` and `--diters`.

Strang splitting tried (viscosity half-step before the predict and after the projection, env-gated
experiment, reverted): nu 3 isolated-disc speed at 4 / 8 / 16 sub-steps 0.065 / 0.035 / 0.029 versus
0.075 / 0.040 / 0.031 without it. The dependence is essentially unchanged, so the operator-splitting
explanation above is at most a small part; the dominant dt-dependent term is still unidentified
(not the iteration counts, not the splitting order). Practical rule stays: >= 12 sub-steps at nu 0.1.

### PBF: Bingham and coarse-graining k (2026-10-05)

**Bingham.** `pbf_bench --tau 2.5` gives bit-identical results to the Newtonian runs (-1.04 / -1.04 /
-1.18 at res 40 / 80 substeps 8 / res 40 substeps 16): `pbf.rs` has no yield-stress model at all (the
Bingham / Herschel-Bulkley effective viscosity is still marked "future work" next to `mean_shear_rate`;
`Params` carries the fields but nothing reads them). The grid track meets the +-3 % gate for paste
(tau 2.5 / 10, n96 vs n128 within 1.6 % / 0.1 %), PBF cannot represent it. Closing this needs a
shear-rate-dependent viscosity in step 7 (Papanastasiou-style regularisation is the obvious first try).

**k.** Unmatched start (`--k 1/2/3`, nu 0.1, res 40, substeps 12, 4 s), peak gravity moment:
k=1 (10 balls) -1.23, k=2 (2 balls) -0.88, k=3 (1 ball) -0.89. Ball counts round (2 x 4 = 8 and
1 x 9 = 9 of 10 particle-masses), so per unit mass -1.23 / -1.10 / -0.98: roughly -10 % and -20 %, with
random packings that differ per run, so this is indicative only, not a pass. Consistent with the
no-slip boundary seeing the coarse radius (drag grows only logarithmically with k). A clean k test
needs a matched start at the same total mass for each k.

#### PBF Bingham: first implementation (2026-10-05)

`pbf.rs` step 7 now uses a Papanastasiou-regularised effective viscosity per particle
(`mu + tau_y (1 - exp(-gamma/gamma0)) / gamma`, unyielded viscosity capped at 101 mu, `BINGHAM_PLASTIC_RATIO`),
from the particle's own SPH shear rate of the previous solve; `morris_weights` takes per-particle mu and the
ball boundary uses the fluid particle's mu. All tests pass.

Grid reference at the same case (`grid_probe --exp E8a --n 96 --nu 0.1 --omega 1 --rho 8 --m 10 --g 9.81 --mu 0.5
--split --shrink 0.15 --newton --t 4 --tau T`): peak gravity moment **-0.865 (tau 2.5) / -0.855 (tau 10)**, against
-1.30 for the Newtonian fluid. PBF (`--tau`): tau 2.5 res 40/80 substeps 8 and res 40 substeps 16
-1.15 / -1.02 / -1.26; tau 10 -1.15 / -1.08 / -1.16; versus -1.04 Newtonian. Varying the plastic ratio
10 / 100 / 1000 / 10000 (tau 2.5, res 40, substeps 12): -1.28 / - / -1.18 / -1.17.

So PBF responds to the yield stress only weakly (the grid peak drops by a third) and sits ~30 % above the grid
with ~+-10 % knob spread: **the Bingham target is not met**. The yield stress changes the PBF result by less than
the stiffness of the regularisation does, which points at the pile dynamics (what the near-ball
viscosity does when the fluid yields) rather than the plug itself; not yet diagnosed. The implementation
is kept (physically motivated, no regression, Newtonian path bit-identical) but unverified against the grid.

### PBF at water-like viscosity (2026-10-05)

Matched start (`--match 1 --m 10`, k = 1, rho_s 8, omega 1, 4 s), peak gravity moment, no grid reference exists
(the grid NaNs at nu 0.01). res 40 / 80 at substeps 8 and 16:

| nu | 40/8 | 80/8 | 40/16 | 80/16 | spread |
|---|---|---|---|---|---|
| 0.01 | -0.94 | -0.95 | -1.06 | -0.81 | about +-13 % |
| 0.003 | -0.80 | -0.74 | -0.83 | -0.65 | about +-12 % |

PBF runs without NaN or clamp trouble, but the result is not converged in either knob (res 80 / substeps 16 is
the lowest and the trend is not monotone), and this is with k = 1, so coarse-graining is not the main cause of the
water spread. Converging water needs a resolved wake/boundary layer, i.e. finer fluid and more sub-steps
(substeps is capped at 16), not a different k.

### PBF: central-force viscosity, squeeze films, converged settings (2026-10-06)

**Method.** The 10-disc pile is chaotic, so single runs mislead (the peak is a maximum statistic, with
sign flips for few balls). Results below are seeded ensembles (`pbf_bench --seed 1..n`, +-0.1 a jitter of
every disc, n = 6-24, mean and relative std of one member; standard error = std / sqrt(n)). The metric is
the 2nd-half mean of the **torque-balance power** `-omega * sum m g x` over balls *and fluid*, printed next to
`P draw`; the peak gravity moment of the discs is kept for the grid comparison. Grid units map to SI as before.

**1. The viscous operator conserved no angular momentum.** The Morris form `c (v_i - v_j)` gives a spurious
force on a rigidly rotating fluid wherever the kernel sum is one-sided (free surface, drum wall). Pure slurry
(no balls, rho_s irrelevant), `P draw` against the torque balance: nu 1: 1069 vs 218 W (5x), nu 0.1: 182 vs 69 W
(2.6x). The central-force form `4 c e (e . (v_i - v_j))` (same Laplacian for incompressible flow, exactly zero
for rigid rotation, angular momentum conserving): 164 vs 167 W and 63 vs 64 W. With balls (10 discs) P draw / balance:
nu 0.1 -3.5 %, nu 1 +12 % (was +73 %). `pbf.rs` uses the central form (`ViscWeights`, unit test
`viscous_operator_vanishes_for_rigid_rotation_even_at_a_free_edge`).

**2. Cost-knob convergence (central form, nu 0.1, 10 discs, ensemble n = 8, peak / torque-balance W):**

| res / substeps | peak | balance |
|---|---|---|
| 40 / 8 | -1150 | 1203 |
| 40 / 16 | -1192 | 1234 |
| 60 / 16 | -1061 | 1080 |
| 80 / 16 | -1053 | 1050 |
| 60 / 32 | -1206 | 1215 |
| 60 / 64 | -1209 | 1230 |
| 80 / 32 | -1188 | 1184 |

Standard error is about 1-1.5 %. Resolution and sub-steps are not independent: res 40 / 16 and res 60 / 32 agree,
res 60 / 16 does not (the sub-step needed shrinks with the particle spacing). The converged level is
about -1.2 (grid -1.30, -7 %). **Converged within +-3 % for res >= 60 and substeps >= 32**, roughly 9x the
cost of res 40 / substeps 8. nu 1: peak -1183 (40/8), -1155 / -1169 / -1174 (res 40 / 60 / 80 at 16), -1131 (60/32),
i.e. +-2 %; the torque balance there is NOT converged in the sub-step (1070, 824, 662 W at 8 / 16 / 32);
the cause is not found (pure slurry with the same wall is stable to +-2 %, so it comes from the pile-fluid coupling).
Ablations at nu 0.1: ball drag closure (6.5) off: no change; overlap-push reaction (3.5) off: spread does not
change once the momentum ledger is kept; squeeze links, gain 1/2/4 and range 1-3x: +-1 %; ball boundary gain 0.25-4:
+-3 %; pbf/dem iteration counts: no change. Step 6.6 (buoyancy) is inactive for balls larger than `h` (the fluid
density at the ball centre is zero; buoyancy then comes from the pressure through 3.5).

**3. Squeeze films** (`squeeze_links`, `mac::lubrication` model, 2D Reynolds film with the calibrated
corrections). Normal film only, approach side only (a symmetric film kept the settled charge jittering at
v_rms 0.3 m/s and resisted the contact solver's separation), stiffness capped at `dt c / m <= 8`. The gain is
small and does not move the level (nu 0.1: -1148 off, -1152 on); it matters for coarse balls (k = 2: 14 %
with a symmetric film).

**4. Coarse-graining k (matched mass: 36 balls at k=1, 9 at k=2, 4 at k=3 on a lattice scaled by k).** Torque
balance, ensemble: k=1 2003 W (n = 24, std 8 %), k=2 2148 W (+7 %, n = 24); k=3 2314 W (+15 %, n = 8, std 12 %:
four balls are too few for a stable mean). Scaling the film gain with k (geometry says one coarse face carries
~k pair films) has no effect (k^0 vs k^0.5: 2148 vs 2142 W). **k is not within +-3 %**: about +7 % at k = 2.

**5. Grid comparison (E8a, 10 discs).** Grid nu 0.1: -1.30; grid nu 1 (new): -1.025. PBF converged: nu 0.1 -1.2
(-7 %), nu 1 -1.13..-1.17 (+10..+14 %): the viscosity dependence of PBF (-1.2 to -1.15) is much weaker than the
grid's (-1.30 to -1.03). Pure slurry E3s (fill 0.3, Fr 0.36, grid `--exp E3s`, P = omega rho T): Re 20 grid 600 W,
PBF 631 W (+5 %); Re 40 grid 386 W, PBF 463 W; at Re 2 and 10 the grid holds a uniform rimming film (54 / 539 W)
while PBF holds a pool (1255 / 889 W): two different equilibria from different starts, not comparable.
Switching surface tension (0.072 N/m) and wettability (0.6) off does not change the pool.

**6. Negative result: drum wall as a Dirichlet boundary** (ghost points on the wall in the viscosity solve, like
the balls'): the pure-slurry sub-step spread drops from 3.6 % to 1.3 %, but P grows 1.2-1.6x against the legacy
per-sub-step blend and moves away from the grid (Re 20: 775 vs 600 W, blend 631 W); gain 0.1-0.5 does not fix it.
Reverted. Open: why the grid at Re 2-10 prefers the film (initial condition?) and the nu 1 torque-balance dependence
on the sub-step.

**7. Re-measure with surface tension and wettability off (`--st 0 --wet 0`, res 60 / substeps 32, n = 8).** No effect:
nu 0.1 peak -1191 (app defaults -1206), balance 1209 W (1215 W); nu 1 peak -1131, balance 664 W (same as before).
So the earlier tables stand and the grid comparison is not distorted by those two forces.

**8. Bingham with the central form (nu 0.1, same settings, peak / balance W).** Newtonian -1191 / 1209; tau 2.5
-1164 / 623; tau 10 -1040 / 146. Grid peak -0.865 (tau 2.5) / -0.855 (tau 10): PBF is +35 % / +22 % above it.
The yield stress cuts PBF's power strongly (1209 -> 623 -> 146 W) but hardly moves the pile's peak (-2 % / -13 %),
whereas the grid's peak drops by a third. The Bingham target is still not met; the cause is not found.

### PBF: water as the default (2026-10-08)

The slurry default is now water (1000 kg/m3, 0.001 Pa*s). Two things followed. (1) A still water pool in the
hydrostatic harness run has residual KE/(M g dx) = 3.85 against the 0.05 gate (the gate was set for 50 Pa*s):
the harness and the DFSPH unit tests stay pinned to 1800 kg/m3 / 50 Pa*s and the gate is not loosened; a still
pool at water viscosity is visibly noisier. (2) `escape_probe` at the water default: 0 balls outside over 6 s.

**App geometry (2 mm balls, 63 mm drum, fill 0.30, 60 % Nc, wet, no coarse-graining; `coarse_graining_probe
--browser-defaults --reference`, 3 seeds, 3 + 4 revolutions, substeps 16).** Mill power against the converged
reference (res 80, substeps 32, 3.026 W); seed sd is <= 1.5 %, so the scatter below is a real resolution effect:

| fluid res | 15 | 20 | 25 | 30 | 40 | 50 | 60 |
|---|---|---|---|---|---|---|---|
| power vs reference | +19.4 % | +10.5 % | +4.6 % | +2.5 % | -4.9 % | -1.0 % | +1.7 % |

Resolution >= 30 is within about +-5 % (not monotone: 40 is the worst), >= 50 within +-2 %; **the +-3 % gate is met
only at res 50 / 60**. Coarse-graining at the app geometry (Manual k = 2): res 25 +20.2 %, res 40 +11.7 %
(not met, consistent with the +7 % of the matched-mass pile test; no coarse-graining is needed at the
shipped defaults anyway: N_true ~ 244 < max_balls). The reference peak ball overlap is 13.4 % of the radius
(target < 10 %).

**Cost (native, `perf_probe`, 120 frames, ms/frame):** res 15 / 150 balls 11.9; res 25 / 300 balls 66.8; res 30 / 300
balls 72.2; res 40 / 600 balls 160.6; res 50 / 600 balls 201.9. Presets are now Realtime 15, Balanced 30, Accuracy 50
(`web/src/params/presets.ts`), i.e. the quality tiers cost +8 % (Balanced) and +26 % (Accuracy) against the
previous 25 / 40 and give +-3 % / +-2 % power instead of +5 % / -5 %.

**Still water pool (`examples/still_pool_probe`, stopped drum, no balls, res 50, 3174 particles).** The fluid never
comes to rest: RMS speed 0.015 m/s and peak 0.05-0.08 m/s, constant from 0.5 s to 10 s (sqrt(g dx) = 0.079 m/s), so
it is particle noise, not a decaying slosh. It is damped by viscosity (RMS 0.0054 at 0.1 Pa*s, 0.00035 at 50 Pa*s),
is independent of surface tension, wettability and sub-steps, and worse at coarser resolution (res 30: 0.017). More
density iterations (3 -> 5 / 8 / 12 / 20) lower it only to about 0.009. Under-relaxing the Jacobi correction by 0.7 / 0.5 / 0.3
lowers the RMS to 0.0065 / 0.0053 / 0.0050 but makes the pool worse compressed (res 50, 50 Pa*s: mean compression 0.138 -> 0.172,
wall-pressure L2 0.20 -> 0.25), so it is not a fix. Raising the adaptive extra iterations 6 -> 30 helps only partly
(mean compression 0.099, wall-pressure L2 0.16, noise 0.0107) at several times the cost. Cause: the pressure-only Jacobi
projection moves information one neighbour layer per iteration, so a deep pool stays compressed (mean compression 0.047 at
res 25, 0.138 at res 50) and the unconverged residue shows up as velocity noise. A real fix is a different pressure solver
(global pressure Poisson solve or multigrid), not a parameter.


### PBF: global pressure solve (PCG) tried and rejected (2026-10-08)

The density projection's Jacobi step is a diagonal approximation of `J J^T lambda = -C` (J: constraint
gradients of the active rows). I replaced it with a preconditioned CG solve of that system
(residual 1e-3 in about 20 iterations, positions `J^T lambda`).

- Still pool (res 50, water): v_rms rose from 0.015 to 0.037 m/s; mean compression only fell 3.0 % to
  2.2 %. So the persistent noise is NOT an unconverged Jacobi iteration (the earlier diagnosis was
  wrong): an exact one-sided projection plus position-based velocity reconstruction itself makes noise.
- Ensemble (nu 0.1, res 40, 6 seeds, k = 1/2/3): torque-balance reference 2048/1551/1651 W (Jacobi) vs
  2022/1565/1663 W (CG); P draw 1751/1292/1227 W vs 1737/1291/1232 W. No change in the k drift.
- Reverted. Noise would need a velocity-level treatment (divergence solve / damping), not a better
  position projection.

### PBF: coarse-graining k re-measured, no correction applied (2026-10-09)

Matched-mass ensembles (`--match 1 --nu 0.1 --t 7 --res 60 --substeps 16`, 16 seeds; k = 1: 36 balls, k = 2: 9,
k = 3: 4), torque-balance reference, mean over the 2nd half. A temporary gain `BALL_BC_GAIN * k^q` on the ball
no-slip boundary was scanned, then removed.

| case | reference (W) | vs k = 1 (2192 W, std 8 %) |
|---|---|---|
| k = 2, q = -1 | 2262 | +3.2 % |
| k = 2, q = -0.5 | 2221 | +1.3 % |
| k = 2, q = 0 (current) | 2126 | -3 % |
| k = 2, q = +0.5 | 2055 | -6.2 % |
| k = 2, q = +1 | 2043 | -6.8 % |
| k = 3, q = 0 (current) | 2395 | +9 % (std 9 %) |
| k = 3, q = -0.5 | 2481 | +13 % |

The earlier "+7 % at k = 2" was within the ensemble noise; k = 2 is about -3 % and k = 3 about +9 % with the
current gain. Scaling the gain with k helps k = 2 but moves k = 3 the wrong way, so it is not adopted. No k
correction is applied. The grid is only an offline reference; the app does not show it.

**Bingham yield-stress sweep (2026-10-09; nu 0.1, 10 discs, res 60, substeps 16, t = 4, 12 seeds; tau in grid units).**

| tau | peak gravity moment | torque-balance reference (W) |
|---|---|---|
| 0 | -1027 | 1045 |
| 2.5 | -1256 | 1315 |
| 10 | -1153 | 1097 |
| 40 | -1056 | 947 |

The pile's peak does not fall with the yield stress (it first rises 22 %), whereas the grid's peak drops by a third
(-0.865 against -1.30). Hypothesis, not tested: `shear_rates` uses fluid neighbours only, so the fluid next to a ball
(which shears against the ball's no-slip surface) looks unsheared and gets the capped viscosity `(1 + 100) mu`, i.e.
extra drag where the real fluid has yielded. The ball ghost velocities would have to enter the shear estimate.

**Ball ghost points in the Bingham shear rate: no effect (2026-10-09, reverted).** Adding the balls' no-slip ghost
points (translation + rotation velocity) as neighbours of the shear-rate estimate gave, for tau 2.5 / 10 / 40, peak
-1261 / -1155 / -1053 and reference 1321 / 1098 / 943 W against -1256 / -1153 / -1056 and 1315 / 1097 / 947 W
before: identical within the 1 % ensemble error. The hypothesis above is refuted; the weak, non-monotone yield-stress
response of the pile peak has another cause (not found).

## Scale-up runner (`mill_core::scaleup`, 2026-10-09; runner since removed)

Note: the separate runner and its probe were removed on 2026-10-10; the impact statistics are now
recorded in the main simulation (`stress.rs`, docs/SCALEUP.md). The measurements below were made
with the runner's own settings (no slurry, `dt = dt_ratio sqrt(d/g)`, 6-12 contact iterations,
restitution threshold and friction speed scale `0.1428 sqrt(g d)`, a lubrication-cushion
restitution `exp(-35/St)` in the wet case) and are not what the main simulation's defaults give.

Probe: `cargo run --release -p mill-core --example scaleup_probe -- --cases conv,sim,sweep,sens,wet`
(3 seeds, 63 mm drum, 2 mm ZrO2, 59 % Nc, dry, 6 recorded revolutions; `± ` is the standard error over
seeds). Key columns: specific power (W/kg, 3D ball mass), impacts per ball per second, mean
`SE/(m g d)` of drop-type and shear-type impacts, stress intensity (mean over all impacts).
Full method and conventions: `docs/SCALEUP.md`.

**Similarity check (met).** Drum and ball scaled together by 1 / 2 / 4 at constant % Nc, ball-relative
numerical settings (`dt = dt_ratio sqrt(d/g)`, restitution threshold and friction speed scale
`0.1428 sqrt(g d)`): mean `SE/(m g d)` drop 0.0767 / 0.0769 / 0.0757, shear 0.2835 / 0.2850 / 0.2838,
stress intensity 0.1735 / 0.1745 / 0.1729 (all within 1.5 %); specific power x1.41 / x2.00 (= sqrt of the
scale), impact rate per ball x0.70 / x0.50. (A first attempt that scaled the geometry but left the speed
at 100 rpm, i.e. a different % Nc, differed by up to +90 %: the probe, not the solver, was wrong.) Unit test
`scaleup::tests::geometrically_similar_mills_give_the_same_dimensionless_stress_energy`.

**Numerical convergence (NOT met at ±3 %).** Relative to nothing finer is known to be converged:

| dt_ratio | iters | spec. power | impacts/ball/s | SE drop | SE shear | intensity |
|---|---|---|---|---|---|---|
| 0.1 | 6 | 1.005 | 149.7 | 0.0734 | 0.268 | 0.162 |
| 0.05 | 6 | 1.064 | 164.8 | 0.0767 | 0.284 | 0.174 |
| 0.025 | 6 | 1.074 | 170.7 | 0.0782 | 0.303 | 0.186 |
| 0.05 | 12 | 1.069 | 148.3 | 0.0827 | 0.304 | 0.189 |
| 0.025 | 12 | 1.085 | 156.3 | 0.0841 | 0.318 | 0.198 |
| 0.0125 | 12 | 1.102 | 173.7 | 0.0865 | 0.313 | 0.198 |
| 0.025 | 20 | 1.099 | 152.1 | 0.0845 | 0.324 | 0.202 |

Specific power and the mean stress energies move by a few percent per halving of `dt` or doubling of the
iterations (intensity 0.162 -> 0.202 across the table); the **impact count does not converge in `dt`**
(148 -> 156 -> 174 at 12 iterations): finer steps resolve more small re-contacts. Treat the impact
frequency as a number that is only comparable between runs with identical settings. Defaults are now
`dt_ratio = 0.025`, `dem_iterations = 12`; comparisons must keep both fixed across conditions.

**D/d sweep and the drop/shear scaling (NOT confirmed).** 59 % Nc, `d` = 2 mm, D/d = 20 / 30 / 45 / 60 / 90
(98 ... 1993 balls): fitted exponent of mean drop `SE/(m g d)` in D/d is -0.31 (expected +1, R^2 0.69),
shear +0.31 (expected -1, R^2 0.80), stress intensity +0.06 (R^2 0.20); shear energy jumps between D/d 45 and 60
(0.288 -> 0.391). The simple free-fall / shear-rate scaling assumes a cataracting (drop-dominated) regime;
at 59 % Nc the charge cascades and the "drop" class is mostly gentle contacts, so the expectation does not
apply, and the power-law extrapolation of the stress energies is **not reliable here**. Only specific
power (exponent 0.41, R^2 0.998) and, less well, the impact rate (0.59, R^2 0.96) follow a clean power
law. The panel's R^2 < 0.95 warning exists for exactly this case. A test at a cataracting speed (e.g. 80-90 % Nc)
has not been done.

**Contact sensitivity (ratios not robust).** +-20 % on a coefficient changes: friction ball-ball: shear SE
+17 % / -8 %; friction ball-wall: shear SE -7 % / +5 %, drop SE -6 % / +3 %; restitution ball-ball: impact rate
-11 % / +14 %, drop SE -1 % / +7 %; specific power stays within +-1 %. So conclusions that rest on
differences below ~10-15 % in stress-energy or impact-rate ratios are within the unknown contact coefficients.

**Wet cushion.** mu = 0 / 0.001 / 0.01 / 0.1 / 1 Pa s: impact rate 165 / 135 / 96 / 94 / 94, drop SE 0.077 /
0.082 / 0.090 / 0.088 / 0.089; the effect saturates once the typical impact Stokes number is below ~100
(restitution already ~0). Specific power is unchanged (+-1 %), i.e. the cushion does not change the power;
wet effects other than impact damping (drag, buoyancy, rheology) are not in this runner. The constant 35 was
not calibrated against `mac/`.

## Impact statistics of the main simulation (`examples/stress_probe`, 2026-10-10)

`cargo run --release -p mill-core --example stress_probe -- --seeds 2 --cases dt,iters,deb`: default scene,
slurry off (63 mm drum, 2 mm ZrO2, 244 balls, 59 % Nc), `Simulation::step` at 60 Hz frames, mean of the
`Metrics.stress` EMA values over 4 revolutions after 4 revolutions of warm-up, 2 seeds. Impact episodes
are debounced (`stress::DEBOUNCE_FACTOR = 0.25`, i.e. re-contacts of a pair within 0.25 sqrt(d/g) = 3.6 ms
count once). Sub-steps below 16 are raised to 16 by `effective_substeps` at this ball size (identical rows).

| 12 iterations | impacts/ball/s | SE drop | SE shear | intensity | W/kg |
|---|---|---|---|---|---|
| 16 sub-steps, no debounce | 136.2 | 0.0820 | 0.3025 | 0.1881 | 1.026 |
| 32, no debounce | 145.2 | 0.0826 | 0.3144 | 0.1953 | 1.041 |
| 64, no debounce | 158.5 | 0.0877 | 0.3116 | 0.1971 | 1.064 |
| 16, debounce 0.25 | 96.2 | 0.1001 | 0.3478 | 0.2150 | 1.026 |
| 32, debounce 0.25 | 99.9 | 0.1032 | 0.3589 | 0.2242 | 1.041 |
| 64, debounce 0.25 | 107.0 | 0.1081 | 0.3457 | 0.2210 | 1.064 |

| 32 sub-steps, debounce 0.25 | impacts/ball/s | SE drop | SE shear | intensity | W/kg |
|---|---|---|---|---|---|
| 2 iterations | 115.8 | 0.0937 | 0.3140 | 0.1926 | 1.036 |
| 6 | 105.0 | 0.1015 | 0.3422 | 0.2135 | 1.027 |
| 12 | 99.9 | 0.1032 | 0.3589 | 0.2242 | 1.041 |
| 20 | 97.7 | 0.1062 | 0.3595 | 0.2266 | 1.037 |

Debounce window sweep (12 iterations; factor x sub-steps 16 / 32 / 64, impacts/ball/s): 0.25: 96 / 100 / 107;
0.5: 81 / 83 / 89; 1: 64 / 66 / 71; 2: 47 / 48 / 50.

Findings (honest): (1) the **impact count is not converged in the sub-step**: with debounce it still rises
+4 % (16 -> 32) and +7 % (32 -> 64); without debounce +7 % / +9 %. Lengthening the debounce window
lowers the absolute count (arbitrarily) but leaves a +7 ... +11 % rise over 16 -> 64 sub-steps, so the
dependence is not just micro-bounces of resting contacts; the debounce is a modelling choice
and the count is comparable only between runs with identical settings. (2) The mean stress energies move
less: intensity 0.215 / 0.224 / 0.221 over 16 / 32 / 64 sub-steps (within ~4 %), and from the default 2 to
12 / 20 contact iterations drop SE +10 % / +13 %, shear SE +14 % / +14.5 %, intensity +16 % / +18 %:
iterations 12 vs 20 differ by 1-3 %, so 12 iterations are enough. The app default (2 iterations) understates
the stress energies by roughly 15 %. (3) Specific power moves < 4 % over all settings (1.03-1.06 W/kg).
Recommended for impact statistics: 32 sub-steps, 12 iterations, coarse-graining off ("Impact statistics"
preset); expected accuracy about 5 % for stress energies, 10 % or worse for impact counts.

### D/d sweep of the main simulation, two speeds (`stress_probe --cases dd`, 2026-10-10)

2 mm ball, drum diameter varied (D/d = 20 / 30 / 45 / 60 = 98 / 221 / 498 / 886 balls), 16 sub-steps, 12
iterations, debounce on, dry, 1 seed (no error bars; seed scatter of earlier runs was about 1-3 %).

| %Nc, D/d | impacts/ball/s | SE drop | SE shear | intensity | W/kg |
|---|---|---|---|---|---|
| 59, 20 | 73.3 | 0.1283 | 0.3213 | 0.2235 | 0.813 |
| 59, 30 | 93.7 | 0.1044 | 0.3421 | 0.2155 | 0.986 |
| 59, 45 | 106.1 | 0.0900 | 0.3706 | 0.2154 | 1.233 |
| 59, 60 | 104.6 | 0.0908 | 0.3667 | 0.2151 | 1.456 |
| 85, 20 | 89.6 | 0.1681 | 0.5070 | 0.3400 | 1.191 |
| 85, 30 | 113.4 | 0.1461 | 0.5128 | 0.3234 | 1.451 |
| 85, 45 | 120.9 | 0.1200 | 0.5974 | 0.3351 | 1.741 |
| 85, 60 | 98.7 | 0.1088 | 0.8257 | 0.4180 | 1.967 |

Between D/d 20 and 60 the log-log slope of the dimensionless drop stress energy is -0.31 (59 % Nc) and
-0.40 (85 % Nc) against the free-fall expectation +1; the shear energy +0.12 and +0.44 against -1;
specific power +0.53 and +0.46 (a clean ~sqrt(D) law at both speeds); stress intensity is flat at 59 % Nc
and rises only at the last point at 85 %. The free-fall / shear-rate scaling laws are **not** supported
at either speed, i.e. stress energy at a larger D/d cannot be extrapolated from smaller D/d by those
laws. Plausible reason (not tested): D/d also changes the charge geometry (number of ball layers, ratio of
free path to ball size), which the velocity-only argument ignores. What does hold is Froude similarity at a
*fixed* D/d (earlier runner check: 1.5 % over a 4x scale): dimensionless quantities are exactly
transferable only between geometrically similar mills.

### Contact-onset impact counting: tried and rejected (2026-10-10)

Counting an impact only on the first sub-step a pair touches (instead of every sub-step whose
pre-solve approach speed exceeds the threshold) was implemented and measured (2 seeds, 12 iterations,
impacts/ball/s at 16 / 32 / 64 sub-steps): without debounce 55.9 / 75.2 / 99.7, with the 0.25 sqrt(d/g)
debounce 52.8 / 67.5 / 82.6, versus 96.2 / 99.9 / 107.0 for the shipped definition (debounce only).
The onset count rises +28 % / +22 % per doubling (finer steps split one physical contact into several
touch/separate cycles) and the mean stress energies fall about 10 % per doubling, so it is **worse**
and was removed. A sub-step-independent impact definition therefore remains open; the shipped definition
is the least dependent measured (+4 % / +7 % per doubling).

### D/d sweep with error bars (`stress_probe --seeds 3 --cases dd`, 2026-10-10)

Same setup as the 1-seed sweep above, 3 seeds, mean ± standard error:

| %Nc, D/d | impacts/ball/s | SE drop | SE shear | intensity | W/kg |
|---|---|---|---|---|---|
| 59, 20 | 73.3 ± 0.3 | 0.1266 ± 0.0012 | 0.329 ± 0.006 | 0.225 ± 0.002 | 0.815 ± 0.002 |
| 59, 30 | 93.7 ± 0.5 | 0.1049 ± 0.0005 | 0.344 ± 0.004 | 0.217 ± 0.002 | 0.991 ± 0.002 |
| 59, 45 | 107.4 ± 2.6 | 0.0900 ± 0.0028 | 0.368 ± 0.010 | 0.214 ± 0.002 | 1.237 ± 0.009 |
| 59, 60 | 99.6 ± 3.0 | 0.0865 ± 0.0022 | 0.398 ± 0.019 | 0.224 ± 0.006 | 1.441 ± 0.008 |
| 85, 20 | 89.4 ± 0.7 | 0.1647 ± 0.0022 | 0.511 ± 0.009 | 0.341 ± 0.004 | 1.183 ± 0.005 |
| 85, 30 | 113.4 ± 0.4 | 0.1444 ± 0.0013 | 0.514 ± 0.001 | 0.323 ± 0.001 | 1.458 ± 0.008 |
| 85, 45 | 122.3 ± 2.2 | 0.1216 ± 0.0019 | 0.586 ± 0.017 | 0.332 ± 0.005 | 1.743 ± 0.015 |
| 85, 60 | 116.8 ± 9.4 | 0.1211 ± 0.0068 | 0.638 ± 0.096 | 0.352 ± 0.034 | 2.054 ± 0.044 |

Slopes between D/d 20 and 60: drop SE -0.35 (59 %) / -0.28 (85 %) (free-fall expectation +1), shear SE
+0.17 / +0.20 (expectation -1), specific power +0.52 / +0.50. The single-seed 85 % / D/d 60 shear value
(0.826) was noise (3 seeds: 0.638 ± 0.096); scatter grows with D/d (more balls, longer revolutions) and
at 85 % Nc. Conclusion unchanged: the free-fall / shear-rate laws are not supported; do not extrapolate
stress energies in D/d.

### CORRECTION (2026-10-10): the D/d sweeps above for D/d >= 60 were coarse-grained

`stress_probe --cases dd` did not turn coarse-graining off, so with the default `max_balls = 600` every row
with more than 600 balls (D/d 60 and larger: 886 balls -> k = 1.2, D/d 100 -> k = 2.0, ...) was simulated
with enlarged balls. **The D/d 60 rows of the 1-seed and 3-seed sweeps and all D/d 100-350 rows measured on
2026-10-10 before this correction are invalid** (they flatly looked "intensity ~ 0.22 for every D/d up to
350"; that was an artefact). D/d 20-45 rows (<= 498 balls) are valid. The probe now forces
`CoarseGrainingMode::Off` and prints the ball count and k per row. It also shows how strongly a small
coarse-graining distorts the statistics: D/d 60 at k = 1.2 gave intensity 0.224 / shear 0.398, at k = 1.00
0.256 / 0.500.

Valid 59 % Nc, k = 1.00, 2 seeds (D/d 45 / 60 / 100; D/d 150 and 250 were interrupted):

| D/d (balls) | impacts/ball/s | SE drop | SE shear | intensity | W/kg |
|---|---|---|---|---|---|
| 45 (498) | 104.9 ± 1.2 | 0.0876 ± 0.0024 | 0.377 ± 0.007 | 0.2164 ± 0.0010 | 1.228 ± 0.005 |
| 60 (886) | 100.9 ± 0.8 | 0.0831 ± 0.0035 | 0.500 ± 0.003 | 0.2557 ± 0.0026 | 1.370 ± 0.001 |
| 100 (2460) | 150.6 ± 1.5 | 0.1080 ± 0.0004 | 0.548 ± 0.003 | 0.2769 ± 0.0020 | 1.666 ± 0.003 |

So with true-size balls the stress intensity does **not** stay flat in D/d: 0.214-0.225 for D/d 20-45,
then 0.256 (60) and 0.277 (100), and the impact frequency per ball jumps from about 100 to 151 between D/d
60 and 100. The trend for a 1 m drum with 2 mm media (D/d = 500, about 61 500 discs) is therefore
unknown, and not flat.

### D/d 150 and 250 with true-size balls (`stress_probe --cases dd --pcs 59`, 2026-10-12)

59 % Nc, k = 1.00, 16 sub-steps, 12 iterations, dry. D/d 150: 2 seeds; D/d 250: **1 seed** (no error bar; the
probe does not report the steady-state drift, so a not-yet-settled charge cannot be excluded for this row).

| D/d (balls) | impacts/ball/s | SE drop | SE shear | intensity | W/kg |
|---|---|---|---|---|---|
| 150 (5535) | 199.3 ± 10.5 | 0.1775 ± 0.0016 | 0.664 ± 0.047 | 0.361 ± 0.010 | 1.872 ± 0.048 |
| 250 (15375) | 208.9 | 0.330 | 1.139 | 0.625 | 2.195 |

Together with the valid rows above, stress intensity vs D/d (20-250): 0.225, 0.217, 0.216, 0.256, 0.277,
0.361, 0.625; impact rate per ball: 73, 94, 105, 101, 151, 199, 209; specific power: 0.81, 0.99, 1.23, 1.37,
1.67, 1.87, 2.20 (growth slows: exponent about 0.30 between D/d 100 and 250, against about 0.5 below D/d 100).
The stress energies grow steeply at large D/d, so no flat or simple power-law trend exists to extrapolate to
D/d = 500 (a 1 m drum with 2 mm media); the D/d 250 jump (x1.7 in intensity from D/d 150) may be a regime
change of the charge (not examined) or an unsettled run. Neither the value nor even the sign of the trend
beyond D/d 250 is known.

### What coarse-graining preserves (`stress_probe --cases cg`, 2026-10-12)

D/d = 45 (90 mm drum, 2 mm ZrO2), 59 % Nc, dry, 16 sub-steps, 12 iterations, 2 seeds, `Manual` coarse-graining
with k = 1 / 1.5 / 2 / 3 (498 / 221 / 125 / 55 balls; total charge mass and footprint preserved). Window means.

| k | power W/m | dissipated W/m | KE J/m | toe deg | shoulder deg | centroid y mm | intensity | shear SE | drop SE |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 7.686 | 7.239 | 0.0425 | 227.5 | 80.5 | -25.20 | 0.216 | 0.377 | 0.088 |
| 1.5 | 7.772 | 7.305 | 0.0697 | 227.4 | 75.8 | -24.50 | 0.221 | 0.351 | 0.105 |
| 2 | 8.018 | 7.549 | 0.0759 | 225.7 | 74.6 | -24.01 | 0.222 | 0.338 | 0.114 |
| 3 | 8.758 | 8.295 | 0.1013 | 221.8 | 70.7 | -23.04 | 0.229 | 0.299 | 0.160 |

Relative to k = 1: power draw per metre +1 % / +4 % / +14 %; dissipated power the same; toe angle -0.1 / -1.8 /
-5.7 deg; shoulder -4.7 / -5.9 / -9.8 deg; charge centroid height +3 / +5 / +9 %; **kinetic energy +64 / +79 /
+138 %** (not preserved); the per-impact statistics are not preserved either: drop SE +19 / +30 / +82 %, shear SE
-7 / -10 / -21 %. The all-impact mean (intensity) moved only +2 / +3 / +6 %, but through opposite-signed
changes of the two classes, and only because the intensity is flat in D/d at D/d_eff = D/(k d) = 15-45
(0.22 for D/d 20-45 above): the coarse-grained statistics match those of a true-size mill at **D/d_eff**
(k = 3 at D/d 45 gives 0.229, true D/d 20 gives 0.225; k = 1.2 at D/d 60 gave 0.224, true D/d 50 is
about 0.23, true D/d 60 is 0.256). So coarse-graining simulates a mill with larger balls, not the real one.

### Coarse-graining at D/d 150 and 250 (`stress_probe --cases cg --cgdd 150|250 --ks ...`, 2026-10-12)

Same protocol as above (59 % Nc, dry, 16 sub-steps, 12 iterations, 4 warm-up + 4 measured revolutions unless
noted). k = 1 at D/d 150 (5535 balls, 1 seed): power 133.5 W/m after 4 warm-up revolutions and 136.0 W/m after
12 (+2 %: settled), so the k = 1 reference is 136 W/m; the 2-seed 1.87 W/kg of the earlier sweep converts to
130 W/m, consistent. At D/d 250 only the 1-seed 2.195 W/kg of the k = 1 sweep exists (424 W/m by conversion).

| D/d | k (balls) | power W/m | vs k = 1 | toe deg | shoulder deg | centroid y mm | KE J/m | intensity | shear SE | drop SE |
|---|---|---|---|---|---|---|---|---|---|---|
| 150 | 1 (5535) | 136.0 | - | 231.8 | 91.6 | -85.7 | 4.04 | 0.368 | 0.633 | 0.195 |
| 150 | 2 (1384) | 159.0 | +17 % | 227.8 | 86.6 | -85.2 | 1.17 | 0.264 | 0.548 | 0.077 |
| 150 | 3 (615) | 162.7 | +20 % | 227.8 | 79.2 | -83.3 | 1.80 | 0.226 | 0.393 | 0.091 |
| 150 | 5 (221) | 160.1 | +18 % | 227.6 | 76.3 | -81.7 | 2.45 | 0.219 | 0.347 | 0.104 |
| 250 | 2 (3844) | 537.6 | +27 % | 229.8 | 90.8 | -141.9 | 13.2 | 0.295 | 0.601 | 0.104 |
| 250 | 3 (1708) | 581.8 | +37 % | 228.0 | 86.0 | -140.9 | 7.31 | 0.279 | 0.565 | 0.084 |
| 250 | 5 (615) | 591.1 | +39 % | 227.6 | 78.1 | -137.9 | 9.17 | 0.231 | 0.400 | 0.091 |

(D/d 150 k = 2 with 12 warm-up revolutions: 159.0 W/m, intensity 0.264; k = 2 at 4 warm-up: 158.6, 0.261.)

Findings: (1) the **power draw is not preserved** at large D/d: +17...20 % at D/d 150 and +27...39 % at D/d
250 (versus +1 / +4 / +14 % at D/d 45), but nearly independent of k for k = 2...5 (160 +- 2 W/m at D/d 150), so
the discrepancy looks like a step between k = 1 and any coarse-graining, not a gradual k effect; the cause is
not known. (2) The charge **centroid height** is preserved within 1 % (k = 2) and 5 % (k = 5) at D/d 150
and the **toe angle** within 4-5 deg, the shoulder angle within 5-15 deg. (3) Kinetic energy and every
per-impact statistic change strongly and not monotonically in k (e.g. D/d 150 KE 4.0 / 1.2 / 1.8 / 2.5 J/m for
k = 1 / 2 / 3 / 5). Stress intensity at k = 2 (0.264) is far from k = 1 (0.368).

## 2026-10-18: approximate-speedup tolerance gate (tolerance_probe)

`examples/tolerance_probe.rs` is the physical-metric A/B gate for solver changes that are not bit-exact. Seeds
1..N (default 3; 8 were used for the decisions below), 2 settle + 3 measured revolutions. Configs: `ui_fast` /
`ui_accuracy` (legacy 1 m drum, coarse-graining Auto), `ship_fast` / `ship_acc` (shipped mill), and a viscosity
sweep `ship_v0.01` / `ship_v1` / `ship_v50` / `ship_acc_v50` (Fast tier unless acc). `--baseline FILE` writes a
CSV, `--compare FILE` gates, `--only a,b` selects configs.

| Metric | Gate |
|---|---|
| power, torque, specific power, impacts | `|delta| <= max(2 sqrt(se_a^2 + se_b^2), 2 %)` |
| toe / shoulder angle | same, with a 2 degree floor |
| unattributed / shaft energy fraction | `|new| <= |baseline| x 1.10 + 0.005` |
| compression mean and worst-5 % | `new <= baseline x 1.10 + 1e-4` |
| shear rate, ms/frame | informational |

The energy rule is relative because the reference solver is already off by 0.17 / -0.04 / -0.09 / +0.38 depending
on config, so an absolute limit would always fail.

**Result for the M7 changes** (BoundaryGroup + simplified pair weight, then CG warm start), 8 seeds vs the
previous solver: 0 gated failures on ship_fast, ship_acc (3 seeds), ship_v1, ship_v50, ship_v0.01 and
ship_acc_v50 (3 seeds vs d7ced07).

| Config | Power delta |
|---|---|
| ship_fast | -0.02 % |
| ship_v1 | -0.05 % |
| ship_v50 | -0.17 % |
| ship_acc (3 seeds; within 2 x SE = 2.4 %) | +2.0 % |
| ship_v0.01 | -0.15 % |
| ship_acc_v50 | +0.19 % |

**Methodology lessons.**

- A 3-seed standard error underestimates the noise: a bit-different but mathematically equivalent solver showed
  +1.3..+2.4 % power shifts in water at 3 seeds that vanished at 8 seeds. Use >= 8 seeds for decisions near the
  2 % floor and compare the same seeds.
- The reference solver is deterministic: a `--compare` against its own baseline gives delta 0.
- The Accuracy tier at the 63 mm drum has mean compression 0.18 / worst-5 % 0.62 even at the reference
  (pre-existing).


## 2026-10-19: strict gates re-measured, ball drag measured (`examples/drag_probe`)

Raw logs: `docs/measurements/2026-10-19/`. Commit f8c2d16 code, then the form-drag factor change below.

**Strict gates (shipped PBF, release, 22 tests run; DFSPH rotating drum / energy not re-run).** 6 pass, 16 fail.
Harness inputs: hydrostatic and rotating drum pinned 1800 kg/m3 / 50 Pa s; spin-up and Taylor-Couette use
`Params::default()` (density 1000); energy closure is the full water default at res 25.

| gate | limit | 2026-10-01 | 2026-10-19 |
|---|---|---|---|
| hydrostatic worst-5 % / mean compression, res 25 | 0.02 / 0.001 | 0.270 / 0.047 | 0.266 / 0.047 |
| hydrostatic worst-5 % / mean compression, res 50 | 0.02 / 0.001 | 0.914 / 0.138 | 0.919 / 0.138 |
| wall-pressure L2, res 25 / 50 | 0.05 | 0.229 / 0.197 | 0.224 / 0.200 |
| spin-up L error at t/tau .02/.05/.1/.2, res 25 | 0.03 | .52/.63/.44/.20 | .514/.607/.419/.174 |
| spin-up wall torque vs dL/dt | 0.03 | 0.468 / 0.451 | **pass** / **pass** |
| Taylor-Couette ball torque rel. error, res 25 / 50 | 0.03 | 0.999 / 1.000 | 0.093 / 0.158 |
| Taylor-Couette wall torque rel. error, res 25 / 50 | 0.03 | 0.966 / 0.895 | 0.571 / 0.929 |
| energy closure, unattributed / shaft work | 0.02 | 0.174 | 0.092 |
| rotating drum 25 to 50, power change (1800 / 50) | 0.03 | 0.181 | 0.055 |
| hydrostatic KE, dry sub-step | | pass | pass |

After the factor change below, the same tests gave Taylor-Couette ball torque 0.106 / 0.153, wall torque 0.606 / 0.943,
energy closure 0.0915, rotating drum 0.0553: within noise of the table.

**Ball drag.** One 10 mm ball (d/dx 2 at res 25, 4 at res 50) is moved through a quiescent full water pool at
speed +U and -U; drag = -(impulse(+U) - impulse(-U)) / (2 dt) per metre depth. Reference
`F = 1/2 rho U^2 d C_D`, with $C_D = 1+10 Re^{-2/3}$ (White; the form could not be confirmed in a literature search,
the subcritical plateau ~1.0-1.2 was), so the reference is uncertain by about 20 %.

| Re | res / sub-steps | R = F_meas/F_ref, closure rate x0 / x0.5 / x1 |
|---|---|---|
| 1e3 | 25 / 8 | 16.1 / 18.4 / 16.3 |
| 1e3 | 25 / 16 | 13.0 / 14.7 / 15.3 |
| 1e3 | 50 / 16 | 17.1 / 16.9 / 16.1 |
| 3e3 | 25 / 8 | 4.7 / 5.2 / 5.0 |
| 3e3 | 50 / 16 | 4.6 / 3.9 / 3.0 |
| 1e4 | 25 / 8 | 0.90 / 0.99 / 1.12 |
| 1e4 | 25 / 16 | 1.57 / 1.77 / 1.89 |
| 1e4 | 50 / 16 | 1.25 / 1.30 / 1.27 |

Standard errors are 5-10 % of R (1e4) and up to 10 % (1e3). With adhesion and cohesion off (res 25 / 8): R = 12.3 (1e3),
4.6 (3e3), 0.95-1.20 (1e4). Result: the total drag in water is 5-17 times the cylinder value for Re 1e3-3e3 whatever the
closure does, and depends on sub-steps at Re 1e4. The closure is therefore not the lever. The pre-registered rule
(case 3: no closure strength satisfies 0.75 <= R <= 1.33 at Re 1e3-1e4) applies: no tuning. The cause of the
excess is not identified (candidates: overlap push, density-constraint noise; not isolated).

**Change made.** The form-drag rate carried no factor 1/2 (it acted as $C_D = 2$); it now reads
`0.5 * rho0 * C_D * (2 r_t) * |v_rel| / m_t`, `C_D = 1`. Same-session A/B, 8 seeds, ship_fast / ui_fast / ship_v0.01:
0 of 27 gated metrics fail (power +0.1 / -1.6 / +0.5 %, toe/shoulder within 2.4 deg, impacts within 2.4 %);
`cargo test --workspace` green (187 unit tests, 5 verification tests), fmt and clippy clean.

**Where the excess drag comes from (same day, temporary per-mechanism recorder, since removed;
log `docs/measurements/2026-10-19/drag/split_res25_s8.log`).** x-impulse on the ball per mechanism, res 25,
8 sub-steps, N/m (reference 0.055 / 0.47 / 5.1):

| Re | total | overlap push (step 3.5) | adhesion | form-drag closure | buoyancy | implicit no-slip solve |
|---|---|---|---|---|---|---|
| 1e3 | 1.01 | 0.62 | 0.39 | 0.004 | 0 | ~0 |
| 3e3 | 2.45 | 1.97 | 0.46 | 0.014 | 0 | 0 |
| 1e4 | 5.06 | 4.58 | -0.23 | 0.71 | 0 | 0 |

The overlap push, which arrests the normal velocity of every fluid particle that enters the ball's contact
radius relative to the ball surface (`pbf.rs` step 3.5), carries 60-95 % of the drag; the form-drag closure is
negligible below Re 1e4, and the viscous no-slip solve is negligible in water. Interpretation (not tested by
disabling the mechanism): the density constraint does not know the ball, so a ball moving through a full pool
displaces fluid only through this inelastic momentum exchange, which acts as a dissipative collision rather than a
pressure; drag is then set by the exchange rule, not by C_D. At Re 1e3 adhesion adds 0.39 N/m (38 %).
Open: replace the exchange by a boundary density/pressure for the ball (the DFSPH track's solid kernel integrals
are the existing design), or test a weaker exchange against this probe; neither is done.
