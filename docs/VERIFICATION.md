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
