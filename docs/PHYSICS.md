# MillDynamics3 — Physics Reference

Equations and algorithms actually implemented in `crates/mill-core/src/{params,geometry,grid,dem,
pbf,coupling,surface,metrics}.rs`, as of the `wip/physics-fixes-metrics-panel` branch. This
documents the real, current source — not the design narrative in `docs/PLAN.md` ss3, whose ss3.2
and ss3.4 text predates the coarse-graining/mass-model and coupling fixes described below (call
sites are cited so the two can be cross-checked; where they disagree, this document and the source
win).

See `docs/PLAN.md` ss3 for the original design rationale (solver-choice tables, milestone history)
and `docs/PARAMETERS.md` for the UI-facing parameter list.

---

## 1. Units & conventions

- SI throughout the core: metres, kilograms, seconds, Pascal-seconds, radians (`params.rs` module
  doc comment). The UI may display mm/rpm and converts in `web/src/params/schema.ts`; nothing in
  `mill-core` itself uses a non-SI unit.
- **2D cross-section, unit depth-of-1-metre convention.** Every ball and fluid particle is modelled
  as a unit-depth (1 m) slice, not a 3D sphere/voxel. Consequently every mass, energy and power
  value produced anywhere in this crate is *per metre of mill length* (`dem.rs` module doc comment
  and `DemStepStats` doc comment); a caller wanting an absolute value for a real mill multiplies by
  the mill's actual axial length.
- World frame: drum centre at the origin, gravity along `−y` (`GRAVITY = -9.81` m/s^2 in both
  `dem.rs` and `pbf.rs`), drum rotates about the origin.
- Angle convention used consistently by `geometry.rs`, `dem.rs` and `metrics.rs`: `atan2(y, x)` in
  `[0, 2*pi)`, measured counter-clockwise from `+x`; "down" (6 o'clock) is `3*pi/2`. `metrics.rs`
  additionally exposes `to_vertical_degrees` to convert to the mill-literature "0 deg at 12
  o'clock, clockwise-increasing for CCW rotation" convention used by acceptance-test literature
  values.
- Determinism: every run is reproducible from `Params` + `seed` (`rng.rs`'s small deterministic
  PRNG seeds initial lattice jitter; `grid.rs`'s `UniformGrid` uses a `BTreeMap`, not `HashMap`,
  specifically so broad-phase candidate-pair iteration order — and hence f32 summation order in the
  contact/kernel accumulations — is identical run-to-run for identical input; see that module's
  doc comment for the full argument). `crates/mill-core/src/dem.rs`'s `ContactBook` similarly uses
  a `HashMap` for the hot per-substep accumulation (measured ~2x faster than `BTreeMap` at default
  ball counts) but replays results through a sorted `Vec` (`sorted_ball_ball`/`sorted_ball_wall`)
  before any pass that iterates *across* contacts and mutates shared state sequentially (friction,
  restitution, rolling resistance), which is where non-deterministic order would actually change
  the physical result.
- Fixed sub-step: `crate::FIXED_DT = 1/240 s` (`lib.rs`) is the sub-step size at the *default*
  `substeps = 4`; `Simulation::fixed_sub_dt() = 1 / (60 * substeps)` is the general form for
  whatever `substeps` an instance's params actually specify. Two ways to advance a `Simulation`:
  `step(dt)` splits a caller-supplied `dt` into `substeps` equal sub-steps
  (`sub_dt = dt / substeps`) and is what the native test suite/benches use; `step_fixed()` advances
  exactly one `fixed_sub_dt()`-sized sub-step, for a caller (the `web/` worker) running its own
  fixed-sub-step accumulator loop against wall-clock time (docs/PLAN.md ss4.1) — the two are
  equivalent when `step(dt)` is called with `dt == substeps * fixed_sub_dt()` (i.e. one nominal
  frame at 1x time scale), see `tests::step_fixed_called_substeps_times_matches_one_step_call_at_1x`.
  A `step_fixed()` caller should call `reset_frame_stats()` once per frame instead of relying on
  `step`'s own internal reset, so `max_substep_displacement_over_diameter`/`coupling_clamp_hits`
  (both "since the last check" diagnostics, see ss8) still cover the whole frame. `simulation.
  dem_iterations` and `simulation.pbf_iterations` control the respective solvers' per-substep
  Gauss-Seidel/Jacobi iteration counts.

---

## 2. Drum geometry & SDF (`geometry.rs`)

The drum (`Drum { radius_m, omega, lifters }`) is represented in its own rotating frame; callers
rotate a world-space point into that frame by `-drum_angle` before querying it.

- **Base wall.** `dist_to_wall(p_local) = radius_m − |p_local|`; positive = free space, negative =
  inside solid (`Drum::sdf`).
- **Lifters.** `lifters.count == 0` (the default) means a perfectly smooth circle: `sdf = dist_to_wall`
  exactly, with no extra evaluation cost. For `count = n > 0`, each of the `n` evenly-spaced bars is
  a trapezoidal cross-section (`lifter_cross_section_sdf`): in its own local (radial `pr`,
  tangential `pt`) frame it spans `pr` in `[R − height, R]` with tangential half-width `base_width/2`
  at the base (`pr = R`) and `top_width/2` at the tip (`pr = R − height`); a possibly-degenerate
  trapezoid (`top_width == base_width` is a rectangle, `top_width == 0` a triangle). The signed
  distance to the convex quad is the minimum, over its four counter-clockwise-wound edges, of the
  point's distance to that edge's infinite line along the edge's inward normal — exact on-edge,
  an approximation near corners (a standard, documented simplification for a `min`/`max`-combined
  convex SDF). `Drum::sdf_lifters` takes the maximum "inside" value over all `n` bars (nearest
  lifter wins) and negates it to match the wall's positive-outside convention; `Drum::sdf` then
  takes `dist_to_wall.min(sdf_lifters(..))` (nearest solid boundary, of wall or any lifter, wins).
- **Normal.** `Drum::sdf_world(p_world, drum_angle)` rotates into the local frame, evaluates `sdf`,
  and takes a central-difference numeric gradient (`EPS = 1e-4`) to get the local normal, falling
  back to `-p_local.normalize_or_zero()` if the gradient is degenerate (e.g. exactly at the drum
  centre); the normal is then rotated back to world space. Points from the wall towards free space.
  A non-finite `p_world` short-circuits to `(0.0, Vec2::ZERO)` rather than letting the numeric
  gradient manufacture a NaN from `-inf − -inf` (documented in detail on `sdf_world`) — this keeps a
  caller's typical `p += -d * normal` correction a safe no-op instead of injecting NaN.
- **Wall velocity.** `wall_velocity(p) = omega * perp(p) = omega * (-p.y, p.x)`, i.e. `omega x r` in
  2D, evaluated at the world-space point currently coincident with that wall material point.

---

## 3. Media / coarse-graining (`params.rs`)

Balls are seeded from an *effective* media population, not the raw UI ball diameter, whenever the
true population would be too large to step in real time.

- **True disc count** (`Params::true_ball_count`):
  `N_true = fill_fraction * packing_fraction_2d * drum_area / (pi * r_true^2)`
  where `drum_area = pi * radius_m^2` and `r_true = ball_diameter_m / 2`. `fill_fraction` (the
  conventional "J" ball-filling fraction) is the charge footprint *including voids*;
  `packing_fraction_2d` (default `0.82`, valid range `[0.5, 0.907]`) converts that footprint into a
  solid disc count — random close packing of equal discs is ~0.82, `pi / (2*sqrt(3)) ~= 0.907` is
  the hexagonal upper bound.
- **Coarse-graining** (`Params::effective_media`): if `N_true > simulation.max_balls`, the true
  population is replaced by `N_sim ~= max_balls` larger, lighter discs with scale factor
  `k = sqrt(N_true / max_balls)`:
  - `diameter_m = true_diameter_m * k` (preserves total footprint area, `N_sim * pi * r_eff^2 ~=
    N_true * pi * r_true^2`, i.e. the fill fraction the user set).
  - `density_kg_m3` is **left unchanged** at the true media density. Because balls are modelled as
    unit-depth discs (`ball_mass = rho * pi * r^2`, ss4 below), preserving the total solid
    footprint area at constant density automatically preserves total charge mass (`mass ~ r^2`, and
    `N_sim * r_eff^2 == N_true * r_true^2` by construction) — no compensating `1/k` density scaling
    is needed or applied. (`Params::effective_media`'s doc comment states this explicitly, and
    `effective_media_preserves_total_charge_mass`/`effective_media_coarse_grains_when_true_count_is_large`
    in `params.rs`'s test module assert it directly.)
  - `ball_count = round(N_true / k^2)` (`~= max_balls`).
  - When `N_true <= max_balls`, `scale_factor = 1.0` and the true media parameters pass through
    unchanged.
- This is a standard coarse-grained DEM approximation: it preserves bulk charge mass and footprint
  area, not the true interstitial void structure or single-collision statistics at the real
  particle size. `EffectiveMedia { true_diameter_m, diameter_m, density_kg_m3, ball_count,
  scale_factor }` carries both the true and simulated values so the approximation is never silent
  (surfaced in the UI's derived-values panel and in `metrics::Metrics`'s
  `true_ball_count`/`simulated_ball_count`/`coarse_graining_factor`/`effective_ball_diameter_m`).

---

## 4. Ball solver — XPBD rigid discs (`dem.rs`)

Non-penetration (ball-ball and ball-wall/lifter) is solved as a rigid, zero-compliance geometric
position constraint, projected iteratively — unconditionally stable at any sub-step size, unlike an
explicit spring-dashpot contact (whose stable time step would need to be on the order of
microseconds at this project's wall speeds/ball stiffness; see `docs/PLAN.md`'s intro section for
the derivation of why that was rejected). Friction, restitution and rolling resistance are added as
position/velocity corrections on top of the converged contact solve.

### 4.1 Ball properties

- `ball_mass(diameter_m, density_kg_m3) = density_kg_m3 * pi * r^2` (`r = diameter_m / 2`) — a
  unit-depth disc's areal mass, kg per metre of mill length. This is the same 2D-slice convention
  the fluid uses (`particle_mass = rest_density * dx^2`, ss5.1), which is what makes ball<->fluid
  momentum exchange (ss6) dimensionally meaningful. `dem.rs`'s doc comment on `ball_mass` notes an
  earlier version used a 3D-sphere mass here instead, which made a fluid particle hundreds of times
  heavier than a coarse-grained ball.
- `ball_inertia(mass, radius_m) = 0.5 * mass * radius_m^2` — moment of inertia of a **uniform disc**
  about its centre (`I = 1/2 m r^2`), consistent with the unit-depth-disc mass model above. This is
  *not* the 3D-sphere value `2/5 m r^2`.
- All balls in a population currently share one radius/mass/inertia (`Balls` struct; a size
  distribution is future work), seeded on a hexagonal lattice filling the drum's cross-section with
  small random jitter (`Balls::seed_lattice`, deterministic via `seed`).

### 4.2 Per-substep algorithm (`DemState::step_with_external_forces`)

Given a fixed sub-step `dt`, drum pose, media parameters and iteration count:

0. **Sanitize.** Any incoming non-finite `x`/`v` (which can only have entered from outside this
   struct, e.g. a corrupted `external` impulse) is reset to `(Vec2::ZERO, Vec2::ZERO)` before
   anything else runs, mirroring the fluid solver's own step-0 guard (ss5.2).
1. **Predict.** For each ball: `v.y += GRAVITY * dt`; if external (fluid coupling) impulses are
   supplied (`Option<&CouplingImpulses>`), `v += impulse * inv_mass` and `omega += angular_impulse *
   inv_inertia` are applied here too — as a direct velocity change, with **no extra `* dt`** (see
   ss6 for why). Then `x += v * dt`, `theta += omega * dt`. Pre-predict velocity (`v_pre`) and
   pre-step position/orientation (`x0`, `theta0`) are saved for later steps.
2. **Broad-phase.** A `UniformGrid` is built over the predicted ball positions with cell size
   `max(2*r*1.05, 2*r + max_ball_speed*dt)` — a fixed 5% pad, or the largest predicted
   per-sub-step displacement if that's bigger. `for_each_candidate_pair` enumerates ball-ball
   candidate pairs once, reused across all solver iterations this sub-step. Ball-wall contact is
   checked directly per ball (no broad-phase needed: one `Drum::sdf_world` query per ball). The
   displacement margin doesn't make this a continuous/swept check (ss9 still applies) — it only
   keeps a fast-moving pair from being lost between sub-steps, so step 3's now-bounded recovery
   below gets a chance to act on it before the overlap becomes severe.
3. **Solve non-penetration** (`dem_iterations` Gauss-Seidel passes, zero compliance, with a
   **bounded recovery rate**):
   - Ball-ball: `C = |x_i - x_j| - (r_i + r_j)`; if `C < 0`, `C_eff = max(C, -MAX_RECOVERY_FRACTION
     * 2r)` (`MAX_RECOVERY_FRACTION = 0.2`), `d_lambda = -C_eff / (w_i + w_j)` (`w = 1/mass`), apply
     `x_i += w_i * d_lambda * n_hat`, `x_j -= w_j * d_lambda * n_hat`, accumulate `lambda_n`
     (computed from the *capped* `C_eff`, so a partially-recovered contact correctly carries less
     normal force this sub-step) for that `(i, j)` pair in `ContactBook`.
   - Ball-wall/lifter: `C = drum.sdf_world(x_i) - r_i`; same `C_eff` cap; if `C < 0`,
     `d_lambda = -C_eff / w_i` (wall has `w_wall = 0`, infinite mass), apply
     `x_i += w_i * d_lambda * n_hat`, accumulate `lambda_n` for ball `i`.
   - **Why the cap.** Recovering a very deep overlap in full, in one iteration, hands step 4 a
     position delta that becomes an unphysically large separation velocity (`Δx/dt`) — the DEM
     half of the fluidised-charge energy-injection bug (ss9/git history). At `dem_iterations = 4`
     the cap still allows recovering up to `0.8` diameters of overlap per sub-step, so ordinary
     small overlaps are unaffected; it only throttles the pathological case.
4. **Reconstruct velocities**: `v = (x - x0) / dt`, `omega = angle_diff(theta, theta0) / dt`.
5. **Friction** (Coulomb-clamped position correction, iterating contacts in a deterministic sorted
   order, with `v`/`omega` kept in sync *per contact* rather than reconstructed once at the end):
   - Ball-ball: `v_t = (v_i - v_j).dot(t_hat) - r*omega_i - r*omega_j`; the raw correction that
     would fully cancel `v_t` this sub-step is `raw = -v_t * dt / w_sum_t` (`w_sum_t = w_i + w_j +
     r^2*w_rot_i + r^2*w_rot_j`), clamped to `[-mu*lambda_n, +mu*lambda_n]`
     (`media.friction_ball_ball`) and applied as a position correction (linear on `x`, rotational on
     `theta`, both scaled by the respective inverse mass/inertia) — **and immediately as the
     matching velocity/`omega` change too** (`dv = w*d_lambda_t*t_hat/dt`,
     `domega = -r*w_rot*d_lambda_t/dt`), so a *later* contact touching an already-corrected ball
     computes its own `v_t` against the true current state, not a stale one. Each individual
     contact's correction is provably non-energy-increasing given the `v`/`omega` it's computed
     against (a standard "reduce relative velocity toward zero" impulse); without the per-contact
     sync, a ball touching more than one contact (the common case under gravity load) could have a
     later contact's correction computed against data that didn't yet reflect an earlier one from
     the same pass, breaking that guarantee by a few percent of the local kinetic energy per
     sub-step — found via the energy-balance regression test in ss9/git history.
   - Ball-wall: same form (including the immediate velocity sync), with the wall's own
     `wall_velocity` substituted for the "other body"'s velocity and no wall-side rotational term;
     `media.friction_ball_wall` is the clamp coefficient. The tangential *position* multiplier
     `d_lambda_t` (not itself an impulse — see ss6's impulse-vs-force distinction) divided by `dt`
     gives the actual tangential impulse the wall delivered; dotted with the wall's velocity there,
     that gives the work done against wall friction this contact (`wall_work_j`, summed over all
     ball-wall contacts) — the wall's *normal* impulse does no work since the wall's velocity is
     purely tangential. An earlier version omitted the `/ dt` here and under-reported power draw by
     a factor of `dt` (240x at this project's default sub-step) relative to `dissipated_energy_j`.
   - Velocities are reconstructed a second time from the total position/orientation change once
     more (redundant with the per-contact sync above in exact arithmetic; kept as a cheap
     self-healing resync against float summation-order drift, same as step 4).
6. **Restitution** (bidirectional, using the *pre-solve* approach velocity `v_pre` to choose the
   target, applied to *every* contact regardless of approach speed):
   - Ball-ball: `v_n_pre = (v_pre_i - v_pre_j).dot(n_hat)`. If `v_n_pre` is more negative than
     `-RESTITUTION_VELOCITY_THRESHOLD` (`0.02` m/s) it's a fresh impact — restitution `e =
     media.restitution_ball_ball` applies and this contact is counted in
     `DemStepStats::collision_count`, with impact energy `E = 0.5 * (mass/2) * v_n_pre^2` (reduced
     mass `mass/2` for two equal masses) binned into `impact_energy_histogram`. Otherwise it's a
     resting/sliding contact and `e = 0`. Either way the target relative normal velocity is
     `target = max(0, -e * v_n_pre)`, and the *actual* current relative normal velocity `v_n_now`
     is driven to that target by an impulse split by inverse mass — **in both directions**, not
     only when `v_n_now` falls short of `target`.
   - Ball-wall: same, with the wall's "infinite mass" (reduced mass is just the ball's own,
     `E = 0.5 * mass * v_n_pre^2` for a fresh impact), `e = media.restitution_ball_wall`.
   - **Why bidirectional.** Step 3's depenetration recovers overlap by moving positions; step 4
     turns that raw position change into velocity with no regard for how much separation speed is
     physically justified. An earlier version of this pass only ever *added* separation velocity
     when `v_n_now` fell short of `target` (bailing out otherwise), so it could never remove
     velocity step 3 had over-manufactured — letting a deep-overlap recovery inject essentially
     unbounded velocity into the charge instead of being capped by the very restitution law meant
     to bound it. This was the other half of the DEM side of the fluidised-charge energy-injection
     bug (the depenetration cap in step 3 is the first half); together they make the solver
     provably unable to increase the ball population's total mechanical energy beyond what the
     wall's own friction work supplies (`dem::tests::total_energy_never_increases_in_a_still_drum`,
     `steady_cascading_charge_never_gains_more_energy_than_the_wall_supplies` in `lib.rs`).
7. **Rolling resistance**: for each ball with nonzero `omega` and nonzero accumulated normal
   impulse (`total_lambda_n`, summed over its ball-ball and ball-wall contacts this sub-step), the
   angular deceleration is capped at `max_delta_omega = media.rolling_friction *
   (total_lambda_n/dt) * r / inertia * dt`, applied toward zero and clamped so it cannot overshoot
   past `omega = 0` (i.e. cannot reverse the sign of `omega` in one sub-step) — a dimensionless
   torque-coefficient model of rolling friction.

`DemStepStats` (per-substep diagnostics, all per-metre-of-mill-depth) returned by this function:
`wall_work_j`, `collision_count`, `impact_energy_histogram` (12 log-spaced bins from `1e-6` J to
`1.0` J, `impact_energy_bin_edges`), `dissipated_energy_j`, and
`max_substep_displacement_over_diameter` (see ss9).

`dissipated_energy_j = e_after_predict + wall_work_j - e_final`, where `e = mechanical_energy_j`
(translational + rotational KE + gravitational PE, `mechanical_energy_j`'s own doc comment) is
evaluated right after the predict step (i.e. including gravity/external forces but before any
contact correction) and again after the full contact solve (steps 3-7) — a net, mechanism-agnostic
dissipation estimate that does not separately attribute loss to friction vs. restitution vs.
rolling resistance vs. step 3's bounded depenetration recovery. Adding back `wall_work_j` is what
keeps this non-negative: it is exactly the residual of the invariant
`crate::tests::steady_cascading_charge_never_gains_more_energy_than_the_wall_supplies` proves holds
every sub-step (`e_final <= e_after_predict + wall_work_j`), so it is bounded below by the same
argument, not by clamping. An earlier version used translational KE only and omitted `wall_work_j`
entirely (`ke_after_predict - ke_final`), which read hundreds of watts negative whenever the
moving wall pumped energy into the charge faster than friction/restitution removed it — an
entirely ordinary situation under active cascading, not a thermodynamics violation, but one an
external review understandably flagged as looking like one. **Does not include fluid-side viscous
dissipation** (ss6.2's drag removes mechanical energy from the balls that never shows up here), so
`Simulation::dissipated_power_w < Simulation::power_draw_w` in a coupled steady state is expected,
not a leak. Regression tests:
`dem::tests::dissipated_energy_is_never_materially_negative_while_cascading`,
`tests::dissipated_power_approaches_power_draw_in_a_settled_dry_charge` (`crates/mill-core/src/lib.rs`).

---

## 5. Slurry solver — PBF + implicit viscosity (`pbf.rs`)

Position Based Fluids with 2D Poly6 (density) / Spiky-gradient (pressure-gradient) kernels, an
iterative density-constraint (incompressibility) solve, drum-wall boundary projection with
no-slip velocity blending, and an implicit (backward-Euler, conjugate-gradient) Newtonian viscosity
solve. `s_corr` artificial pressure is implemented but disabled.

### 5.1 Seeding & kernels

- `FluidParticles::seed_lattice`: spacing `dx = drum_radius_m / resolution`, kernel radius `h = 2 *
  dx`, `rest_density = slurry.density_kg_m3`. Particles are placed on a **hexagonal** lattice
  (columns `dx` apart, rows `row_h = dx * sqrt(3)/2` apart), so `particle_mass = rest_density * dx *
  row_h` (2D unit-depth convention, matching `ball_mass`) — the hex-cell area, **not** the
  square-lattice `dx^2` an earlier version of this document stated. Deriving both the mass and the
  seeded site count from that same hex-cell area keeps the total seeded mass, the seeded footprint
  (`fill_fraction`), and the lattice's own Poly6 density sum (which reads `rest_density` to within
  ~0.1% right at seeding) all consistent at once; using the square-lattice cell for the mass while
  seeding hexagonally left every run ~15% over-dense at `t = 0` in an earlier version. Fill:
  `slurry.fill_fraction` of the drum's cross-sectional area from the bottom upward; lattice sites
  overlapping an already-seeded ball are skipped.
- **Poly6** (density kernel): `W(r^2, h) = 4/(pi*h^8) * (h^2 - r^2)^3` for `r < h`, else `0`.
- **Spiky gradient**: `grad_W(delta, r, h) = -30/(pi*h^5) * (h-r)^2 * (delta/r)` for `0 < r < h`,
  else `0` — points from `x_i` toward its neighbour (Spiky decreases with `r`).
- `FluidParticles::densities()`: `rho_i = sum_j m_j * W_poly6(r_ij, h)`, including the particle's
  own self-contribution at `r = 0`.

### 5.2 Per-substep algorithm (`FluidParticles::step_coupled`; `step` is the ball-free special case)

0. **Sanitize.** Any incoming non-finite `x`/`v` (which can only have entered from outside this
   struct, e.g. a diverged ball position read in steps 3.5/6.5) is reset to `(Vec2::ZERO,
   Vec2::ZERO)` before anything else runs — every correction this solver itself performs is
   individually bounded, so this guard exists purely against externally-introduced NaN/inf.
1. **Predict**: `v.y += GRAVITY * dt`; `x += v * dt`.
2. **Neighbour lists**: a `UniformGrid` over the predicted positions with cell size `h`;
   `build_neighbor_lists` keeps only candidate pairs actually within `h` (mutual: `j` in
   `neighbors[i]` implies `i` in `neighbors[j]`), rebuilt once per sub-step and reused across every
   density-constraint iteration.
3. **Density-constraint solve** (`pbf_iterations` Jacobi-style passes — all deltas computed from
   current positions, then applied together):
   - `C_i = rho_i/rest_density - 1`.
   - `lambda_i = -C_i / (sum_j |grad_W_ij/rest_density|^2 + grad_self^2 + EPSILON_RELAX)`,
     `EPSILON_RELAX = 200.0` (CFM-style relaxation against divide-by-zero for sparse
     neighbourhoods).
   - `delta_p_i = 1/rest_density * sum_j (lambda_i + lambda_j + s_corr_ij) * grad_W_ij`, applied
     directly to `x`.
   - **`s_corr` is implemented but disabled** (`S_CORR_K = 0.0`, formula kept as `-k *
     (W(r)/W(delta_q))^n` with `delta_q = 0.2h`, `n = 4`): at this project's real-SI scale
     (`rest_density` ~1000-2000 kg/m^3, small `h`), the literature-default `k = 0.1` produced a
     correction 10-30x larger than the density-constraint terms it's meant to gently supplement,
     causing runaway dispersal instead of preventing clustering. Disabling it gives a stable
     settled puddle at ~1% mean density error over 2 s of simulated time; several smaller `k`
     values were tried without finding a stable working point. Documented in `pbf.rs` as revisit-if
     visual clustering artefacts are observed, not resolved.
4. **Ball overlap projection** (coupling step 1 — see ss6.1).
5. **Boundary projection** (position only): for any particle with `drum.sdf_world(x) = d < 0`,
   `x += -d * normal`; marks the particle `touched_wall` for the no-slip blend below.
6. **Velocity reconstruction**: `v = (x - x0) / dt` for every particle (`x0` = pre-predict
   position).
7. **No-slip wall blending**: for particles marked `touched_wall`, `v = (1-beta)*v + beta*v_wall`
   (`v_wall = drum.wall_velocity(x)`), `beta = slurry.wall_no_slip` clamped to `[0,1]`.
8. **Ball viscous drag** (coupling step 2 — see ss6.2).
9. **Buoyancy** (coupling step 3 — see ss6.3).
10. **Implicit Newtonian viscosity** (`mu = slurry.viscosity_pa_s`; skipped entirely, at zero
    iteration cost, when `mu == 0`): density is recomputed once more at the final,
    post-boundary-projection positions, `morris_weights` are built, and `solve_implicit_viscosity`
    solves the diffusion system in place. `mean_shear_rate` is computed from the same final
    positions/densities/neighbour lists.
11. **Fluid speed clamp** (stability backstop): any particle speed above `v_max =
    FLUID_SPEED_SAFETY_FACTOR * (|omega|*radius_m + sqrt(4*|GRAVITY|*radius_m))`,
    `FLUID_SPEED_SAFETY_FACTOR = 5.0`, is rescaled down to `v_max`. Deliberately not a CFL-style
    `h/dt` bound (that would be coupled to `simulation.resolution` and could clip real motion at
    high resolution instead of only blow-ups). Momentum removed here is *not* credited back to any
    ball.
12. **Coupling impulse clamp** (see ss6.4).

### 5.3 Implicit viscosity detail

The Newtonian viscous term `(mu/rho) * laplacian(v)` is discretised with the Morris (1997) SPH
viscous Laplacian (`morris_weights`):

```
c_ij = m_j * (mu_i + mu_j) / (rho_i * rho_j) * |x_ij . grad_W_ij| / (|x_ij|^2 + eta^2)
eta^2 = VISCOSITY_ETA_FACTOR * h^2,  VISCOSITY_ETA_FACTOR = 0.01
```

(single-phase fluid here, so `mu_i = mu_j = mu` and `mu_i + mu_j = 2*mu`). `(Lv)_i = sum_j c_ij *
(v_i - v_j)` is the resulting graph Laplacian (`laplacian_apply`) — symmetric and positive
semi-definite since every `c_ij >= 0`. The backward-Euler system `(I + dt*L) v_new = v` is solved
matrix-free by conjugate gradient (`solve_implicit_viscosity`), operating directly on the `Vec2`
field with the real inner product `<a,b> = sum_i a_i . b_i` (mathematically equivalent to running
scalar CG on the x- and y-components independently, without the bookkeeping of splitting them).
Relative-residual stopping tolerance `VISCOSITY_CG_TOLERANCE = 1e-3` (`|r| <= tol * |b|`), warm-started
from `v` itself, iteration ceiling `VISCOSITY_CG_MAX_ITERS = 50`. `A = I + dt*L` is SPD whenever `dt >
0`, so CG is guaranteed to converge in exact arithmetic; the 50-iteration ceiling is a generous
practical bound, not a tuned-tight one. This replaced an earlier explicit XSPH mixing scheme whose
effective kinematic viscosity saturated around `~h^2/dt` regardless of the input value, so it could
not represent a genuinely low-viscosity (near-water) slurry; the implicit solve is unconditionally
stable for any `viscosity_pa_s >= 0` at the fixed sub-step and uses the value directly in Pa*s.

`mean_shear_rate`: `gamma_dot = sqrt(2 * D:D)`, `D = sym(grad v)`, with `grad_v_i = sum_j (m_j/rho_j)
* (v_j - v_i) (x) grad_W_ij` — the standard SPH velocity-gradient estimator (same neighbour-weighting
form the density constraint and viscosity already use). Surfaced via `FluidStepStats` for
`metrics::Metrics::mean_shear_rate_per_s` and reserved for a future Bingham/Herschel-Bulkley
extension (not implemented; `Rheology::Bingham` exists in `params.rs` as a reserved enum variant but
only `Rheology::Newtonian` is wired into the solver).

`nu = mu/rho` uses each particle's own *measured* SPH density, not `slurry.density_kg_m3` directly
— these agree closely in a well-relaxed interior but diverge near a low-density boundary (free
surface, sparse neighbourhood), where the resulting locally elevated `nu` is itself physically
reasonable (a thin/rarefied film resists shearing more, relative to its own mass, than the bulk).

---

## 6. Two-way ball<->fluid coupling (`coupling.rs`, `pbf.rs` steps 3.5/6.5/6.6/7.5/8)

`coupling::step` orchestrates one shared sub-step: the fluid solves first (`step_coupled`), using
*last* sub-step's ball positions/velocities as a fixed boundary and producing this sub-step's
reaction impulses on the balls; the balls then advance using those impulses
(`DemState::step_with_external_forces`), so their new state becomes the fluid's boundary for the
*next* sub-step (a staggered/semi-implicit scheme).

**Impulses, not forces.** `CouplingImpulses { impulses: Vec<Vec2>, angular_impulses: Vec<f32>,
clamp_hits: u32, fluid_momentum_change: Vec2 }` accumulates a linear + angular impulse per ball. The
exchange is expressed and applied as impulses (`delta_v = impulse * inv_mass`, no extra `* dt`)
rather than as forces re-integrated with another `* dt` on the ball side; the doc comment on
`CouplingImpulses` notes an earlier force-based version double-round-tripped through `dt` and was
visibly unstable (balls flung out of the charge) at this project's default sub-step rate.

### 6.1 Overlap projection (pbf.rs step 3.5)

After the density-constraint solve, any fluid particle within `contact_radius = balls.radius +
ADHESION_CONTACT_MARGIN * h` (`ADHESION_CONTACT_MARGIN = 0.05`, previously `0.25` — see ss6.1a for
why it shrank) of a ball centre is pushed out along the connecting normal by `push_mag =
min(contact_radius - dist, 0.5 * balls.radius)` — capped at half the ball's radius so a particle
deeply embedded in the contact zone cannot produce an outsized single-substep *position* correction.
The grid used to find a ball's nearby fluid particles here is sized to cover the wider of this step's
`contact_radius` and ss6.1a's `adhesion_radius` (it used to default to `2 * balls.radius`, which was
smaller than `contact_radius` at this project's defaults and silently missed real neighbours —
regression test: `pbf::tests::coupling_finds_ball_neighbours_beyond_two_radii`).

**The momentum exchanged is bounded separately from the position correction.** The full push
implies a velocity of `push_mag / dt` once step 5 reconstructs velocity from position — unbounded
in momentum as `dt` shrinks or as particles pile up against a ball, since the position cap alone
doesn't limit it. Instead, the actual momentum exchanged is `impulse_mag = mass * min(push_mag/dt,
arrest_speed)`, where `arrest_speed = max(0, -(v_i - v_surface).dot(n_hat))` is how much of the
particle's velocity *toward* the ball's surface (`v_surface = balls.v + balls.omega x lever`) this
contact actually needs to arrest — zero for a particle already moving away, or at rest relative to
the ball. The ball's Newton's-third-law reaction is `-impulse_mag * n_hat` (angular component via
`lever = n_hat * balls.radius`); the part of the position push beyond `impulse_mag` (a purely
geometric correction, position without momentum) is tracked separately and subtracted back out of
the fluid particle's own velocity right after step 5's reconstruction, so it never shows up as
implicit fluid speed.

This split fixes one contributor to the coupling half of the fluidised-charge energy-injection bug:
before it, a purely resting/settled overlap (particle not actually approaching the ball,
`arrest_speed = 0`) still exchanged the full `mass * push / dt` every sub-step it persisted, behind
a stability clamp (ss6.4) permissive enough to let it through at this project's default parameters.
This mechanism alone was not the dominant one, though — see ss6.2 for the fix that turned out to
matter most for the reported symptom (a still-drum coupled charge failing to settle at all).

### 6.1a Ball<->slurry adhesion / media wettability (pbf.rs step 3.6)

ss6.1's overlap projection only ever pushes fluid *away* from a ball's surface — nothing pulled it
back, a geometric 180 degree (fully non-wetting) contact angle baked into the discretisation
regardless of any parameter. `slurry.wettability` (`[0, 1]`, default `0.6`) adds the missing
attraction: for a fluid particle in the shell `(shell_start, adhesion_radius]` just outside
`contact_radius`,

```
shell_start     = contact_radius + ADHESION_SHELL_DEAD_ZONE * h   (ADHESION_SHELL_DEAD_ZONE = 1e-3)
adhesion_radius = contact_radius
                  + min(ADHESION_RANGE_FACTOR * h, ADHESION_RANGE_BALL_RADII * balls.radius)
                                                    (ADHESION_RANGE_FACTOR = 1.0,
                                                     ADHESION_RANGE_BALL_RADII = 0.5)
accel_mag       = wettability * ADHESION_ACCEL_FACTOR * |GRAVITY| (ADHESION_ACCEL_FACTOR = 2.0)
```

`shell_start` (not `contact_radius` itself) is where the shell begins: ss6.1's push lands an
overlapping particle *exactly* at `contact_radius`, and without this small dead zone, ordinary
floating-point rounding in that push could leave the particle a hair past `contact_radius`, which
this step's smooth taper (zero in the limit, but not exactly zero a few ULPs off the boundary) would
then treat as genuine shell membership — caught by
`coupling::tests::an_approaching_overlapping_fluid_particle_gives_the_ball_the_opposite_reaction`'s
momentum-conservation check, which has no tolerance for a second mechanism sneaking in.

Each shell particle gets a smooth bump weight `taper = 1 - (2t - 1)^2`, `t = (dist - shell_start) /
(adhesion_radius - shell_start)` (zero at both edges), producing a taper-weighted mean pull direction
`pull_dir` (toward the weighted-mean fluid side) and a coverage fraction `coverage = min(1,
sum(taper))` — a fully-surrounded ball (symmetric coverage, `dir_bar ~= 0`) has no net direction to
pull and is skipped, same as a fully-entrained ball feeling no net drag in ss6.2.

**The ball's own velocity response is relaxed by the contributing fluid mass, exactly like ss6.2's
`a_lin`.** At this project's coupling resolution a ball's shell typically holds only one or two
sub-resolution fluid particles, so `m_contrib = particle_mass * sum(taper)` is routinely orders of
magnitude below `balls.mass`. Applying the target closing velocity `accel_mag * coverage * dt`
straight to the ball and letting the exact-conservation split recoil it back onto that tiny fluid
mass gives the fluid an unphysical velocity spike (the same mass-ratio amplification ss6.2 already
had to fix for viscous drag) — masked, not caught, by the fluid speed clamp (ss6.4) absorbing the
excess on the fluid side only, which broke momentum conservation. The fix mirrors ss6.2 exactly:

```
attach = m_contrib / (balls.mass + m_contrib)          <= 1, -> 0 as m_contrib -> 0
dv_b   = pull_dir * (accel_mag * coverage * dt * attach)
```

`coupling.impulses[b] += balls.mass * dv_b` (plus the matching angular term via `lever = pull_dir *
balls.radius`); each contributing particle `j` absorbs `dv_j = -(balls.mass * frac_j / mass) * dv_b`
(`frac_j = taper_j / sum(taper)`), so the fluid side's total momentum change is exactly `-balls.mass *
dv_b` regardless of `m_contrib`, `balls.mass`, or how many particles share the shell — regression
test: `coupling::tests::adhesion_pulls_a_ball_and_nearby_fluid_together_only_when_wettability_is_positive`.

**Balls are processed sequentially (Gauss-Seidel), not Jacobi**, for the identical reason as ss6.2: a
fluid particle can sit in more than one ball's shell at once at this project's coupling resolution,
and computing every ball's pull against the same stale fluid snapshot would let their individually-
bounded reactions stack.

**Applied as a deferred velocity delta, not directly to `self.v`.** Unlike ss6.1's push (a position
correction whose *excess* velocity is tracked and subtracted back out), step 3.6 moves no position at
all — it is a pure velocity kick. Applying it to `self.v` in place at step 3.6's point in the
pipeline would be silently discarded by ss5.2 step 5's velocity reconstruction (`v = (x - x0) / dt`,
driven only by position deltas, which runs immediately after). Step 3.6's contribution is instead
accumulated separately and added back into `self.v` right after that reconstruction (and after ss6.1's
excess subtraction), the same "defer past the reconstruction" pattern ss6.1 already uses for its own
excess term.

**Deliberately not a true Young's-equation contact angle.** That needs a matching fluid-fluid
cohesion term, which this crate does not (re-)implement (`S_CORR_K` stays `0`, ss5.2 step 3's
artificial-pressure term) — an explicit scope cut, not an oversight. `wettability` only fixes the
more basic defect: slurry previously could not cling to media at all, regardless of any parameter.

**Two defects fixed after an external review reported ball+slurry clumps floating indefinitely
through the air, unrelated to either mechanism above.**

1. **The shell was wider than the ball.** `ADHESION_RANGE_FACTOR * h` alone (no cap) reached
   roughly 2.6 ball radii out from the surface at this project's default coupling resolution
   (fluid `h` several times a coarse-grained ball's own radius) — fluid held at that range reads
   as a separate floating clump rather than a wetting film. `ADHESION_RANGE_BALL_RADII` caps the
   shell at half the ball's own radius regardless of how fine the fluid resolution is.
2. **Nothing tested whether there was any actual liquid there.** A single sub-resolution particle
   passing through the shell pulled at full strength, identical to a particle touching the edge of
   a real pool — there was no submergence, local-density, or neighbour-count gate. Combined with
   (1)'s wide shell, a ball with a couple of stray shell particles could pick up roughly
   `wettability * 2g` of pull, comparable to gravity itself, holding the clump together
   indefinitely rather than letting it disperse.

The fix (`adhesion_wetness`) gates each contributing particle's `taper` weight by how "wet" its
own already-computed SPH density (`pbf.rs` step 3's `density[]`, one iteration stale — cheap, and
adequate for a gate) reads, before it contributes to `weight_sum`/`dir_sum`/`coverage`/`m_contrib`/
`attach` and the final `frac`-weighted reaction split:

```
self_frac = particle_mass * poly6(0, h) / rest_density   (the density an isolated particle with
                                                            zero real neighbours reads, ~0.28 at
                                                            this project's defaults)
wet_j     = ((density_j / rest_density - self_frac) / (1 - self_frac)).clamp(0, 1)
taper_j  *= wet_j
```

`wet_j = 0` for an isolated particle (its measured density is exactly `self_frac * rest_density`,
having no real neighbours) and approaches `1` as `density_j` approaches `rest_density` (genuine
bulk liquid) — so a ball surrounded only by scattered stray droplets now receives essentially no
pull, while a ball at the edge of a real pool still wets normally. Because a small *mutually-close*
cluster of droplets locally reads an elevated (if still sub-`rest_density`) density among
themselves, the gate reduces rather than fully eliminates adhesion for a tight multi-droplet clump
that isn't connected to any bulk liquid — a known, accepted limitation of using purely local
density as the wetness signal, not requiring a matching fluid-fluid cohesion term (see the
scope-cut note above) to fix properly. Regression tests:
`coupling::tests::an_airborne_ball_does_not_carry_a_floating_slurry_clump` (several mutually-
isolated scattered droplets around an airborne ball produce ~zero net impulse) and
`coupling::tests::adhesion_pulls_a_ball_and_nearby_fluid_together_only_when_wettability_is_positive`
(a particle with enough nearby density still gets pulled).

### 6.2 Viscous no-slip drag (pbf.rs step 6.5)

Each ball relaxes toward its locally-entrained fluid mass via a centre-of-mass relaxation, not a
per-particle blend. Every fluid particle `i` within `h_c = balls.radius + h` of a ball contributes
a smooth poly6 taper weight `phi_i = poly6(|r_i|^2, h_c) / poly6(0, h_c)` (`r_i = x_i - x_b`), giving
an entrained mass `m_ent = sum_i particle_mass * phi_i` and a weighted mean fluid velocity `v_bar =
(sum_i particle_mass * phi_i * v_i) / m_ent`:

```
tau   = rho_ball * r^2 / (4 * mu)                     (Stokes-regime disc relaxation-time closure)
beta  = slurry.ball_no_slip * (1 - exp(-dt / tau))     (unchanged closure)
a_lin = beta * m_ent / (balls.mass + m_ent)            <= 1 for every mu, dt
dv_b  = a_lin * (v_bar - v_b)
```

`rho_ball = balls.mass / (pi * balls.radius^2)`; `mu = slurry.viscosity_pa_s` is the *physical*
slurry viscosity. The reaction is distributed back across the same weighted fluid neighbours
(`dv_i = -(balls.mass / m_ent) * dv_b * phi_i`), so linear momentum is conserved exactly
(`sum_i particle_mass * dv_i == -balls.mass * dv_b`). The angular part mirrors this using
`balls.inertia` and a tangential taper field in place of `balls.mass` and `v_bar`, plus a
`balls.mass * cross2(r_bar, dv_b)` correction to the ball's angular impulse -- the reaction torque
from spreading the linear reaction over an off-centre weighted mean position `r_bar`, needed for
*exact* (not just per-mechanism) angular momentum conservation once that reaction isn't applied
through the ball's own centre.

`a_lin <= 1` by construction, so the ball can never be driven past `v_bar` regardless of `mu`, `dt`,
or the entrained/ball mass ratio. An earlier version applied `beta` independently to *every* nearby
fluid particle and let the ball absorb the sum of all those reactions -- amplifying the ball's
actual response by the entrained/ball mass ratio (`~9x` at this project's defaults), which the
stability clamp (ss6.4) then hid as "occasional clamping is normal" instead of surfacing as the bug
it was (the ball<->fluid coupling clamp rate went from ~3.4% of ball-substeps at 0.5 Pa*s to ~50%
at 50/200 Pa*s before this fix). In the small-`beta`, entrained-mass-dominated limit (`m_ent >>
balls.mass`, `dt << tau`) the corrected formula reduces analytically to ordinary 2D Stokes drag,
`impulse ~= 4*pi*mu*(v_fluid - v_ball)*dt`, independent of the ball/fluid mass ratio -- verified
directly by `pbf::tests::ball_drag_matches_two_dimensional_stokes_scaling`.

**Balls are processed sequentially (Gauss-Seidel), not Jacobi.** `self.v` is updated immediately
after each ball's contribution, not accumulated into a separate array and applied once after every
ball has been visited. At this project's coupling resolution `h_c` easily spans several
neighbouring balls' worth of a packed charge, so a fluid particle commonly sits within more than
one ball's entrainment radius at once; computing every ball's `v_bar` against the same stale fluid
state and only summing the results afterwards let several individually-bounded (`a_lin <= 1`)
corrections stack on the same fluid particle well beyond what any single ball's own bound was meant
to allow. Sequential application means each ball after the first already sees the fluid's
up-to-date state, so its own relaxation is bounded against reality rather than against a snapshot
every other overlapping ball is also independently correcting.

**The angular reaction's exact-conservation normalizer needs a relative, not absolute, threshold.**
The tangential field distributing the ball's angular reaction back to the fluid with zero net
linear momentum is scaled by `c_rot = -i_b * dw_b / s`, where
`s = sum_i mass*phi_i^2*(|r_i|^2 - cross2(r_i, p_bar))` is a geometric normalizer that can land
anywhere from comparable to its own positive-definite part (`sum_i mass*phi_i^2*|r_i|^2`) down to
many orders of magnitude smaller, essentially at random, whenever the entrained neighbour count is
small (routine at this project's fine coupling resolution relative to a ball's size — 2-4 fluid
neighbours per ball is typical, not an edge case). Guarding this division with a fixed absolute
epsilon (`|s| > 1e-9`) let a modest torque get divided by an almost-cancelled denominator and
amplified by 2-4 orders of magnitude (`c_rot` observed in the hundreds to ~13000, in a *stationary*
drum with no rotation at all to drive it) -- this dominated the reported fluidised-charge symptom
for the coupled (slurry-enabled) case: with the overlap-push and Gauss-Seidel fixes above alone, a
ball+slurry charge released in a still drum still failed to settle (`v_rms ~2.5` m/s, fluid pinned
at its speed clamp, ss5.2 step 11). The fix compares `|s|` against a *fraction* of its own
positive-definite part instead (`|s| > 0.1 * sum_i mass*phi_i^2*|r_i|^2`), which correctly detects
near-total cancellation regardless of scale; when it fails, the rotational reaction (both `dw_b` and
its fluid-side redistribution) is skipped for that ball this sub-step, same as the pre-existing
degenerate case. Regression test: `coupling::tests::a_coupled_charge_in_a_still_drum_settles_with_slurry_on`.

### 6.3 Buoyancy (pbf.rs step 6.6)

Balls are typically sub-resolution relative to the fluid spacing and do not contribute to the
density-constraint sum, so buoyancy is modelled directly rather than emerging from the PBF pressure
field. Each ball samples the local fluid density by reusing the same Poly6 kernel and the fluid's
own step-2 neighbour grid at the ball's position:

```
rho_local = sum_j m_j * W_poly6(|x_ball - x_j|^2, h)   (over nearby fluid particles j)
rho_eff = min(rho_local, rest_density)
impulse_on_ball = (0, -rho_eff * (pi * r^2) * GRAVITY * dt)   (acts through the centroid: no torque)
```

Clamping to `rest_density` avoids over-buoyancy from a locally compacted pocket and tapers smoothly
to zero as a ball nears the free surface (lower sampled density there) instead of an on/off cutoff.
The reaction is applied immediately as a velocity change split across the contributing fluid
particles, weighted by each one's share `w_j / rho_local` of the sampled local density — this must
run after step 6's velocity reconstruction, not before, or the reconstruction would overwrite it.
This mechanism did not exist in the original coupling design (`docs/PLAN.md` ss3.4 describes
buoyancy as emergent from the density-constraint push alone, with an Akinci-style boundary-particle
extension noted as future work); it is now a distinct, directly-modelled step.

### 6.4 Stability clamp (pbf.rs step 8)

The accumulated per-ball impulse this sub-step is clamped to
`impulse_clamp(mass, dt) = F_CLAMP_G_MULTIPLE * mass * 9.81 * dt`, `F_CLAMP_G_MULTIPLE = 60.0`
(angular impulse clamped to that magnitude times `balls.radius`). `pbf.rs`'s doc comment on
`F_CLAMP_G_MULTIPLE` records that this constant used to be `3.0`, tuned down purely to survive an
earlier 3D-sphere-vs-2D-disc mass mismatch between balls and fluid particles, then `20.0` once both
were unit-depth discs with dimensionally consistent masses. It is `60.0` specifically because of
step 6.2's drag term, which is legitimately larger (though still `a_lin <= 1` bounded) than the
overlap-push mechanism (ss6.1) once `beta` saturates (`mu` above roughly 10 Pa*s at this project's
ball sizes) -- a velocity-matching relaxation's natural impulse scale doesn't shrink with `dt` the
way a position-capped correction does, so a clamp tuned only against the latter would clip
physically legitimate no-slip corrections at the former's scale. This is unaffected by ss6.1's
overlap-push fix: that fix bounds the *momentum* the overlap push can exchange by the same kind of
physical argument the clamp exists to backstop, so it doesn't change how large a legitimate drag
correction can be -- lowering the clamp to match ss6.1's now-much-smaller contribution was tried and
confirmed wrong (it broke `pbf::tests::ball_drag_never_overshoots_the_local_fluid_velocity`, a
scenario with no overlap-push contribution at all). Residual clamp rate at the project's default
sub-step rate is ~5.5-6.5% of ball-substeps at 50 and 200 Pa*s
(`coupling::tests::cascading_charge_keeps_coupling_clamp_hits_rare_once_settled_at_*`, tightened
from an earlier `<15%` bound to `<10%` now that ss6.1's contribution is negligible), now
attributable almost entirely to step 6.2's drag, not ss6.1's overlap push -- persistent clamping
still indicates a real problem, just a rarer one than before either fix. `clamp_hits` (count of
balls clamped this sub-step) and `fluid_momentum_change` (sum of fluid momentum change from every
mechanism above, for a Newton's-third-law cross-check against `sum(impulses)` when no clamp fired)
are both tracked in `CouplingImpulses` purely as diagnostics.

---

## 7. Free surface extraction (`surface.rs`)

Produces closed (rarely open) polylines approximating the slurry free surface for rendering, via:

1. **Splat.** Every fluid particle splats a smooth, compactly-supported cubic-falloff kernel
   (`splat_kernel(r2, radius) = (1 - r2/radius^2)^3` for `r2 < radius^2`, else 0; unnormalized —
   only its threshold-relative shape matters) of radius `SPLAT_RADIUS_FACTOR * fluid.h`
   (`SPLAT_RADIUS_FACTOR = 1.5`) onto a `GRID_SIZE x GRID_SIZE` (`GRID_SIZE = 128`) scalar occupancy
   field spanning `[-drum.radius_m, drum.radius_m]` on each axis (`build_field`).
2. **Mask.** Any grid point outside the wall or inside a lifter bar (per `drum.sdf`, evaluated in
   the drum-local frame) is zeroed — this is what stops the marching-squares contour from bulging
   through the wall in the first place.
3. **Threshold.** Marching squares runs at `THRESHOLD_FRACTION * phi_full`
   (`THRESHOLD_FRACTION = 0.5`), where `phi_full` is computed **analytically**
   (`analytic_phi_full`), not read off the grid's measured maximum:
   `phi_full = number_density * integral(K dA) = (rest_density / particle_mass) * (pi *
   splat_radius^2 / 4)`. `number_density = rest_density / particle_mass = 1 / (dx * row_h)` (the
   hex lattice's true number density, `row_h = dx * sqrt(3)/2` — see ss5.1's `particle_mass`
   correction) is exact at the PBF density constraint's equilibrium; the kernel's integral has
   the closed form `pi * splat_radius^2 / 4` (substitute `u = r^2/splat_radius^2`). By symmetry, the
   field value exactly at a flat bulk/free-surface boundary is exactly half the deep-bulk value, so
   the `0.5` contour lands precisely on the physical surface. Using the grid's own measured maximum
   instead is unsafe: at high rotation speed the slurry can centrifuge into a thin wall-hugging film
   whose field values are depressed everywhere, while a transient over-compacted PBF pocket
   elsewhere can inflate the measured maximum well past the true bulk value and starve the film's
   threshold, making the real free surface invisible to marching squares. If the analytic threshold
   yields no contour at all (a degenerate case), extraction retries once with the grid's measured
   peak as the threshold reference instead, as a graceful fallback (`extract_surface`).
4. **Contour.** Standard marching squares (`cell_edge_pairs`, with the two ambiguous 4-corner cases
   resolved by a fixed, non-adaptive choice) produces segments per cell, which are chained into
   polylines by exact `EdgeId` matching (`extend_chain`).
5. **Smooth & clamp.** One pass of Chaikin corner-cutting (`chaikin_smooth_closed`) softens the
   blocky contour; every resulting point still outside the wall (up to one grid cell of overshoot
   from linear interpolation/smoothing) is projected back onto the wall along its SDF normal
   (`project_inside_wall`).

`flatten_polygons` serialises the result as `[n_polys, len_0, x, y, ..., len_1, x, y, ...]` for the
wasm/frontend boundary.

---

## 8. Metrics & grinding instrumentation (`metrics.rs`, `lib.rs`)

`metrics::compute(balls, fluid, drum, drum_angle)` is a pure function of the current populations and
derives:

- **Toe/shoulder angles** (`charge_toe_shoulder`): balls within `WALL_MARGIN_BALL_RADII (2.5) *
  balls.radius` of the wall approximate the charge's outer layer; if there are at least
  `MIN_WALL_BALLS (6)` such balls and their angular gap exceeds `CENTRIFUGE_GAP_THRESHOLD_RAD (0.3
  * pi)` (i.e. not centrifuged), the two ends of the largest angular gap in that wall-adjacent
  population are the toe/shoulder (assignment to toe vs. shoulder depends on `drum.omega`'s sign).
  `None`/`None` otherwise.
- **Charge centroid** (`charge_centroid`): plain unweighted mean of all ball positions (all balls
  share one mass, so this is the correct centroid).
- **Slurry pool angular extent** (`slurry_pool_angular_extent`): fluid particles within
  `POOL_WALL_MARGIN_H (1.5) * fluid.h` of the wall (per `drum.sdf_world`, so lifters count),
  clustered the same wraparound-aware way as toe/shoulder.
- **Free-surface line fit** (`fluid_free_surface_line`): fluid particles *not* near the wall (the
  pool's own exposed top), fit by the standard 2x2 covariance-matrix principal-eigenvector method;
  returns `(angle_rad mod pi, offset_m)` where `offset_m` is the signed perpendicular distance of
  the fitted line from the drum centre. `None` if fewer than 3 such particles.
- **Pool depth at bottom** (`slurry_pool_depth_at_bottom`): among fluid particles within
  `POOL_DEPTH_BAND_RAD (+/- 5 deg)` of straight-down, `drum.radius_m - min(distance_from_centre)`.
- **Lacey mixing index** (`mixing_index`): particles binned into a `MIXING_GRID_SIZE (16) x 16`
  grid; `M = (S0^2 - S^2) / (S0^2 - Sr^2)` where `S^2` is the variance of per-cell mean dye value
  (unweighted across occupied cells), `S0^2 = p*(1-p)` is the fully-segregated variance (`p` =
  overall mean dye fraction), `Sr^2 = S0^2 / n_bar` is the fully-mixed variance using the average
  occupied-cell particle count as the per-cell sample size. Clamped to `[0,1]`; degenerate cases
  (near-zero denominator) return `1.0`.
- **Total ball kinetic energy** (`total_kinetic_energy_j`): `sum(0.5 * mass * |v|^2)`.
- **Max ball-ball overlap fraction** (`max_ball_overlap_fraction`): `(2r - dist)/r` over
  broad-phase candidate pairs (same `UniformGrid` cell size as the DEM contact solve) — note the
  denominator is the ball **radius**, not the diameter: two centres coincident reads `2.0` (200%).
- **Max ball-wall overlap fraction** (`max_ball_wall_overlap_fraction`): `(r - d)/r` (`d` from
  `Drum::sdf_world`, same per-radius normalisation as above) over every ball. Companion to
  `max_ball_overlap_fraction`, which is ball-ball only and leaves wall penetration invisible.
- **Max/mean fluid density error fraction** (`max_fluid_density_error_fraction`,
  `mean_fluid_density_error_fraction`): `|rho_i - rest_density| / rest_density`, reusing
  `FluidParticles::densities()`. **Includes free-surface/boundary neighbour deficiency**, which the
  one-sided PBF density constraint (ss5.2 step 3) never corrects by design -- this reading can look
  large (tens of percent) on an entirely healthy run, since it is structurally dominated by
  under-dense free-surface/near-wall particles in any SPH-family method. Not, by itself, evidence
  the solver has failed to converge; see the compression-error pair below for that.
- **Max/mean fluid compression error fraction** (`max_fluid_compression_error_fraction`,
  `mean_fluid_compression_error_fraction`): `(rho_i/rest_density - 1).max(0.0)`, term-for-term the
  one-sided quantity the PBF density constraint's own `C_i` (ss5.2 step 3) actually drives toward
  zero, so this is the honest convergence readout -- a healthy run keeps it small (a percent or two)
  regardless of how large the absolute-value reading above looks. Added after an external review
  read the absolute-value density error (72% max, 16-26% mean on an ordinary run) as evidence of a
  broken incompressibility constraint; the compression-only reading over the same population stays
  small. Regression: `metrics::tests::settled_puddle_compression_error_is_small`.

**Grinding/solver diagnostics are not derivable from `(balls, fluid, drum)` alone** — `compute`
zero/empty-defaults these fields; `Simulation::metrics()` fills them in from its own accumulated
state after calling `compute`:

- `power_draw_w`, `torque_nm`, `collision_rate_per_s`, `dissipated_power_w`,
  `impact_energy_histogram` (`ImpactEnergyHistogram { bin_edges_j, counts_per_s }`): EMA-smoothed
  (time constant `GRINDING_STATS_EMA_TAU_S = 1.0` s, i.e. `alpha = 1 - exp(-sub_dt/tau)` applied
  every sub-step in `Simulation::update_grinding_stats`) from `DemStepStats`. `power_draw_w =
  wall_work_j / sub_dt`; `torque_nm = power_draw_w / omega` (0 when `|omega| <= 1e-6`, since torque
  is undefined rather than infinite at zero rotation); `collision_rate_per_s =
  collision_count/sub_dt`; `dissipated_power_w = dissipated_energy_j/sub_dt`; each histogram bin's
  EMA is `count/sub_dt`.
- `coupling_clamp_hits`: sum of `CouplingImpulses::clamp_hits` over every sub-step of the most
  recent `Simulation::step` call — reset each call, **not** EMA-smoothed (meant to read as "did
  this happen just now", not a smoothed rate).
- `effective_ball_diameter_m`, `simulated_ball_count`, `true_ball_count`, `coarse_graining_factor`:
  from `Params::effective_media`/`Params::true_ball_count` (ss3).
- `max_substep_displacement_over_diameter`: the largest value of
  `DemStepStats::max_substep_displacement_over_diameter` seen across every sub-step of the most
  recent `Simulation::step` call — reset each call, **not** EMA-smoothed (a tunnelling-risk spike
  is exactly the kind of transient an average would hide). See ss9.
- `mean_shear_rate_per_s`, `viscosity_solver_iterations`: from the most recently completed
  sub-step's `FluidStepStats`.

---

## 9. Known limitations

- **Tunnelling / continuous-collision guard is diagnostic, with a partial mitigation.**
  `DemStepStats::max_substep_displacement_over_diameter = |v| * dt / (2*radius)` (the fraction of
  a ball's own diameter it moved in one sub-step) is computed every sub-step and surfaced through
  `Simulation::max_substep_displacement_over_diameter` and `metrics::Metrics::
  max_substep_displacement_over_diameter`. There is still **no actual swept/continuous-collision
  correction** — the broad-phase (`dem.rs` step 2) now includes a per-sub-step displacement margin
  (ss4.2) so a fast pair isn't *lost* between sub-steps, and step 3's recovery is now
  rate-limited (`MAX_RECOVERY_FRACTION`, ss4.2) so a deep overlap decays over several sub-steps
  instead of becoming a single unphysical velocity spike (ss6's restitution fix removes any excess
  regardless) — but a ball can still, in principle, pass fully through another within one discrete
  sub-step without either ever registering as a candidate pair at all. This was observed at
  ~0.85-0.86 at this project's *former* default parameters (`max_balls = 2000`, `substeps = 4`) —
  already at the point where tunnelling risk is a live, not merely theoretical, concern, and
  plausibly connected to the fluidised-charge energy-injection bug this crate's history documents
  (a solver already at its tunnelling limit gives the (now-fixed) unbounded depenetration/coupling
  mechanisms above the largest overlaps to react to). The current defaults (`max_balls = 600`,
  `substeps = 8`, see `docs/PARAMETERS.md`) bring the typical ratio to roughly `0.16` at the
  default mill speed **with no lifters** -- comfortably under 1, though this is a UI-configurable
  parameter, not a solver-enforced bound, so a user can still push it back into risky territory
  (larger `max_balls`, fewer `substeps`, faster rotation, smaller media). **With `lifters.count = 8`
  at the default speed, both `max_substep_displacement_over_diameter` and `max_ball_overlap_fraction`
  read measurably higher** than the no-lifter figure above -- cataracting balls launched off a lifter
  reach higher peak speeds than the no-lifter cascading case. Measured in-browser (default params,
  `lifters.count = 8`, steady cataracting after ~20 s of sim time):
  `max_substep_displacement_over_diameter` ~0.30 (vs. ~0.16-0.21 with no lifters),
  `max_ball_overlap_fraction` ~25-27% (vs. ~20-45%, similar order). An external review's screenshot
  of a lifters-enabled run reading `substep_displacement = 0.402` is consistent with this pattern
  (measurably higher than the no-lifter figure, not comparable to it), not evidence of a regression
  -- both stay well under the XPBD `< 1` criterion. See `docs/METRICS.md` for the full measured
  reading, including the metrics not discussed here.
  - **A separate, notable finding from the same measurement: `max_ball_wall_overlap_fraction`
    reads ~31-33% with `lifters.count = 8`, versus ~1-2% with no lifters**, and this was
    consistently elevated across repeated readings (not a one-off transient). Plausible
    contributors, not yet root-caused: the lifter SDF's near-corner approximation
    (`Drum::sdf_lifters`, ss2) and/or genuinely harder impacts against a lifter face during
    cataracting (a ball striking a lifter has less opportunity to be caught by a gradually
    increasing normal force the way a shallow approach to the smooth circular wall does). Worth
    a follow-up investigation; not fixed as part of this pass.
- **`s_corr` artificial pressure is implemented but disabled** (`S_CORR_K = 0.0` in `pbf.rs`) — see
  ss5.2. Revisit only if visual clustering artefacts are observed; no working non-zero `k` was found
  for this project's SI parameter scale.
- **Bingham/Herschel-Bulkley rheology is not implemented.** `params::Rheology::Bingham` exists as a
  reserved enum variant (with `yield_stress_pa` already a validated `SlurryParams` field) but only
  `Rheology::Newtonian` is wired into `pbf.rs`'s viscosity solve; `mean_shear_rate` is computed and
  surfaced specifically in anticipation of this future extension.
- **No ball size distribution.** `Balls` (and `EffectiveMedia`) model one shared radius/mass/inertia
  for the whole population; per-particle size variation is not implemented.
- **Buoyancy has no boundary-particle contribution to the fluid's own density field.** Balls do not
  contribute to the PBF density-constraint sum (`step_coupled` step 3), so a ball never "occupies"
  space from the fluid's perspective the way an Akinci-style boundary-particle scheme would; the
  directly-modelled buoyancy impulse (ss6.3) is a substitute for, not an emergent consequence of,
  incompressibility.
- **2D cross-section only.** Every value in this crate is per metre of unit mill depth; no axial
  (out-of-plane) transport, end-wall effects, or 3D particle shape are modelled. This is a
  deliberate, project-wide scope decision (see `docs/PLAN.md`), not a defect, but it means results
  are qualitative relative to a real 3D mill.
- **The broad-phase displacement margin (ss4.2) covers the largest single-ball speed, not the
  largest pair *closing* speed.** `cell_size = max(2r*1.05, 2r + max_speed*dt)` uses `max_speed`,
  the fastest any one ball moves this sub-step; a head-on pair (each moving toward the other at
  close to `max_speed`) closes at up to twice that rate, so the margin can in principle be half of
  what step 2's own rationale intends for that specific pair. Widening it to cover the worst-case
  pairwise closing speed would cost more candidate pairs per sub-step (and, at the broad-phase
  level, cannot distinguish a genuinely closing pair from a merely fast one); left as a known gap,
  not fixed, pending the M6 performance pass creating headroom for a broader margin.
- **`dissipated_power_w` excludes fluid-side viscous dissipation.** It is a purely DEM-side
  accounting (ss4.2's `dissipated_energy_j`, net of the wall's own work input); the slurry's
  viscosity (ss5.3) also removes real mechanical energy from the balls (ss6.2's drag) that never
  appears in this metric. A coupled steady-state run legitimately shows
  `dissipated_power_w < power_draw_w`; this is not a missing-energy bug.
- **No fluid-fluid cohesion/surface-tension term.** `S_CORR_K` stays `0` (ss5.2 step 3); the only
  attraction anywhere in this solver is the ball<->fluid adhesion term (ss6.1a, gated by
  `slurry.wettability`), not a true Young's-equation contact angle. The free surface has no surface
  tension.
