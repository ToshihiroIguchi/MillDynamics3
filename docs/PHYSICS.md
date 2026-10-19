# MillDynamics3 — Physics Reference

Equations and algorithms actually implemented in `crates/mill-core/src/{params,geometry,grid,dem,
pbf,coupling,surface,metrics,stress,lib}.rs` and `mac/lubrication.rs` (the squeeze-film model that
`pbf.rs` imports), as of the 2026-10-20 source on `main`. This documents the real, current source
— not the design narrative in `docs/PLAN.md` ss3, whose ss3.2 and ss3.4 text predates the
coarse-graining/mass-model and coupling fixes described below (call sites are cited as
`path:line` so the two can be cross-checked; where they disagree, this document and the source
win; line numbers drift when comments are edited, so search for the named symbol if a line is off).

**Reading convention.** Sections 0–8 describe the *current* algorithm. Dated passages (for
example "2026-09-26") are *history*: they record what was measured or fixed on that date, with the
defaults, presets and mechanisms of that date. Where such a passage names a default (`50 Pa·s`,
`max_balls = 150`, a 1 m drum, `1/240 s`), that was the value at the time, not today's; the current
defaults are listed in §3.1.

See `docs/PLAN.md` ss3 for the original design rationale (solver-choice tables, milestone history)
and `docs/PARAMETERS.md` for the UI-facing parameter list.

---

## 0. Overview: what is simulated and how

This section is the big picture. Sections 1–9 give the exact equations, constants and known
limitations of each piece.

### 0.1 The physical system

A horizontal tumbling ball mill, seen as a **2D cross-section** perpendicular to its axis. A circular
drum (optionally with lifter bars) rotates at a prescribed angular speed `omega`. Inside are
(a) many rigid circular grinding balls and (b) a viscous liquid slurry. Gravity acts along `-y`.
Everything is per metre of axial length (unit-depth discs). The drum drives everything: the wall
drags balls by friction and slurry by viscous no-slip, the charge is lifted, and then cascades or
cataracts back down. Balls and slurry push on each other (buoyancy, drag, displacement).

### 0.2 Models at a glance

| Part | Physics | Algorithm | Source | Details |
|---|---|---|---|---|
| Drum | Rigid rotating boundary, optional trapezoidal lifters | Signed-distance function (SDF) in the drum frame; wall velocity `omega x r` | `geometry.rs` | §2 |
| Balls | Rigid-body contact: non-penetration, Coulomb friction, restitution, rolling resistance, gravity | XPBD discs: predict, project overlaps (Gauss–Seidel), rebuild velocities, then friction and restitution | `dem.rs` | §4 |
| Slurry | Incompressible, viscous (Newtonian) liquid with surface tension and wetting | Position Based Fluids (SPH-style particles with a density constraint) plus implicit viscosity | `pbf.rs` | §5, §6.1a |
| Ball–slurry coupling | Buoyancy, form drag, viscous drag (ball-surface no-slip in the implicit solve), sub-grid squeeze films, displacement of fluid by balls | Impulse exchange between the two solvers each sub-step; joint implicit ball–fluid viscosity solve | `coupling.rs`, `pbf.rs`, `mac/lubrication.rs` | §5.3, §6 |
| Free surface | – (visualisation only) | Splat particles to a grid, marching squares | `surface.rs` | §7 |
| Outputs | Power, torque, energy budget, 3D-converted impact (stress) statistics, mixing | Accumulated from the contact/wall work and the recorded impacts of each sub-step | `metrics.rs`, `stress.rs`, `lib.rs` | §8 |

### 0.3 One simulation sub-step

`Simulation::step_fixed` advances one fixed sub-step `dt = 1/(60 * effective_substeps)` s
(`Simulation::fixed_sub_dt`, `lib.rs:426`). The requested `simulation.substeps` defaults to 8
(`params.rs:488`); `Params::effective_substeps` (§3.3) raises it when the media are small, to 16 at
the shipped 63 mm drum / 2 mm ball defaults, i.e. `dt = 1/960 s`. (`FIXED_DT = 1/240 s`,
`lib.rs:41`, is used only by the benches.) One sub-step runs `coupling::step_recorded`
(`coupling.rs:71`):

1. **Slurry first.** Using the balls' last positions and velocities as moving obstacles, advance the
   fluid (`FluidParticles::step_coupled`, `pbf.rs:1174`):
   1. apply gravity and predict positions;
   2. build neighbour lists (uniform grid with cell size = kernel radius `h`);
   3. iterate the **density constraint** so the particles stay at rest density, which makes the
      fluid incompressible (this is the "pressure" solve of PBF); the pass count is adaptive;
   4. push fluid out of balls (momentum-bounded), add **cohesion** (surface tension) and
      **adhesion** (wetting) forces, and project particles out of the wall;
   5. rebuild velocities from the position change, then blend in the wall velocity (no-slip);
   6. apply **form drag** and **buoyancy** as impulses, recording the equal and opposite reaction
      on each ball;
   7. solve **implicit viscosity** `(I + dt L) v = v*` by conjugate gradient, jointly with the
      ball velocities (ball-surface no-slip boundary and sub-grid squeeze-film links) when balls
      are near the fluid. Implicit means it is stable at any viscosity.
2. **Balls second.** Apply the fluid reaction impulses, then (`DemState::step_recorded`,
   `dem.rs:530`):
   1. predict (gravity, `x += v dt`);
   2. broad-phase pair search on a counting-sort uniform grid;
   3. clamp each ball's predicted advance to its swept time of impact with the wall or lifters
      (`Drum::toi_swept`, step 2.5; only for severe predicted overlaps);
   4. project ball–ball and ball–wall/lifter overlaps (`C = distance - (r_i + r_j) >= 0`) with a
      capped recovery rate;
   5. rebuild velocities from the position change;
   6. apply Coulomb friction (static/kinetic blend), restitution (with an energy-injection guard)
      and rolling resistance.
3. **Bookkeeping.** Rotate the drum by `omega * dt`. Accumulate wall work (power, torque), the
   energy budget and the recorded impacts (§8).

The staggered order (fluid with old ball state, then balls with new fluid impulses) is semi-implicit
and keeps the coupling stable.

### 0.4 Key ideas behind the algorithm choices

- **Why position-based (XPBD) contact instead of spring–dashpot DEM?** Real ball stiffness would
  require microsecond time steps at these speeds. Treating contact as a rigid constraint solved by
  projection is stable at any sub-step; the price is that contact forces are reconstructed, not
  integrated (§4).
- **Why PBF for the slurry?** Particle methods handle the violent free surface (splashing,
  cascading) with no mesh, and PBF's density constraint is unconditionally stable. Viscosity is
  solved implicitly so that one code path covers water-like slurry (the default, 1 mPa·s) up to
  paste-like slurry (tens of Pa·s).
- **Coarse-graining.** Optional (default Off). If the real ball count is too large for real-time
  stepping, a smaller number of larger balls can be used at the same total mass and footprint area
  (`Auto` or `Manual` mode, §3). The UI shows when this is active.
- **Determinism.** The same `Params` and `seed` always give the same trajectory: fixed time step,
  deterministic neighbour ordering and a seeded PRNG.

### 0.5 What the model does *not* capture

It is a qualitative 2D model: no axial flow and no 3D packing (2D void fraction is much lower
than in 3D), no particle breakage or grinding kinetics (only collision/stress statistics are
measured), no slurry-borne fine solids, only a regularised Bingham yield-stress option (§5.3, no
Herschel–Bulkley law) in the shipped solver, and slurry percolation into a settled bed is limited
by the fluid resolution. §9 lists the measured limitations.

### 0.6 Experimental solvers (not used by the shipped simulator)

`crates/mill-core/src/fluid/` (a DFSPH replacement for PBF) and `crates/mill-core/src/mac/` (a
fixed-grid Navier–Stokes track with cut cells, multigrid-preconditioned CG, resolved ball coupling
and Herschel–Bulkley rheology) are verification and rebuild experiments, documented in
`docs/VERIFICATION.md`. `Simulation` runs DEM + PBF only: `lib.rs` and `coupling.rs` call neither
`fluid/` nor `mac/`. The single cross-use is that `pbf.rs` imports the squeeze-film model
`crate::mac::lubrication` (`pbf.rs:552`, inside `squeeze_links`), which the shipped solver *does*
use (§6.5). Everything else in `mac/` and all of `fluid/` is unused by the simulator.

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
  PRNG seeds initial lattice jitter; `grid.rs`'s `UniformGrid` is a counting-sort CSR structure —
  no hash map at all since the 2026-09-27 M6 solver pass, see docs/PERF.md — built so broad-phase
  candidate-pair iteration order — and hence f32 summation order in the contact/kernel
  accumulations — is identical run-to-run for identical input; see that module's doc comment for
  the full argument). `crates/mill-core/src/dem.rs`'s `ContactBook` and `pbf.rs`'s neighbour lists
  are likewise plain arrays (parallel to the broad-phase's own pair list, or CSR-indexed) rather
  than a hash map, for the same reason plus raw speed in that hot per-substep accumulation loop;
  `ContactBook` still replays results through a sorted `Vec` (`sorted_ball_ball`/`sorted_ball_wall`)
  before any pass that iterates *across* contacts and mutates shared state sequentially (friction,
  restitution, rolling resistance), which is where non-deterministic order would actually change
  the physical result.
- Fixed sub-step: the solver steps by `Simulation::fixed_sub_dt() = 1 / (60 *
  effective_substeps)` (`lib.rs:426`; `effective_substeps` is the requested `simulation.substeps`,
  default 8, auto-raised for small media, §3.3 — 16 at the shipped defaults, so `1/960 s`).
  `crate::FIXED_DT = 1/240 s` (`lib.rs:41`) is a reference constant used only by
  `benches/step.rs`, never by the solver; it equals `fixed_sub_dt()` only at `effective_substeps
  = 4`. Two ways to advance a `Simulation`: `step(dt)` (`lib.rs:456`) splits a caller-supplied `dt`
  into `effective_substeps` equal sub-steps (`sub_dt = dt / n`) and is what the native test
  suite/benches use; `step_fixed()` (`lib.rs:438`) advances exactly one `fixed_sub_dt()`-sized
  sub-step, for a caller (the `web/` worker) running its own fixed-sub-step accumulator loop
  against wall-clock time (docs/PLAN.md ss4.1) — the two are equivalent when `step(dt)` is called
  with `dt == n * fixed_sub_dt()` (i.e. one nominal frame at 1x time scale), see
  `tests::step_fixed_called_substeps_times_matches_one_step_call_at_1x`. A `step_fixed()` caller
  should call `reset_frame_stats()` once per frame instead of relying on `step`'s own internal
  reset, so `max_substep_displacement_over_diameter`/`coupling_clamp_hits` (both "since the last
  check" diagnostics, see ss8) still cover the whole frame. `simulation.dem_iterations` (default 2)
  and `simulation.pbf_iterations` (default 3, a floor — §5.2) control the respective solvers'
  per-substep Gauss-Seidel/Jacobi iteration counts.

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
  trapezoid (`top_width == base_width` is a rectangle, `top_width == 0` a triangle). The four
  counter-clockwise-wound vertices are `(R, -b/2)`, `(R, b/2)`, `(R-H, t/2)`, `(R-H, -t/2)` in
  `(pr, pt)` (`lifter_cross_section_sdf`, `geometry.rs:211`). The signed distance to this convex
  polygon is **exact** (`convex_polygon_sdf`, `geometry.rs:247`): the magnitude is the minimum
  Euclidean distance to the four edges as closed *segments* (the projection parameter is clamped to
  `[0,1]`, so a point whose nearest feature is a vertex measures to that vertex), and the sign
  comes from a separate point-in-convex-polygon half-plane test (positive inside the solid, in this
  helper's own convention). An earlier version took the minimum distance to each edge's *infinite*
  line, which is exact in the interior but under-estimates the distance to an exterior vertex; that
  approximation no longer exists (history: §9). `Drum::sdf_lifters` (`geometry.rs:50`) takes the
  maximum "inside" value over all `n` bars (nearest lifter wins) and negates it to match the wall's
  positive-outside convention; `Drum::sdf` (`geometry.rs:37`) then takes
  `dist_to_wall.min(sdf_lifters(..))` (nearest solid boundary, of wall or any lifter, wins).
- **Normal.** `Drum::sdf_world(p_world, drum_angle)` rotates into the local frame, evaluates `sdf`,
  and takes a central-difference numeric gradient (`EPS = 1e-4`) to get the local normal, falling
  back to `-p_local.normalize_or_zero()` if the gradient is degenerate (e.g. exactly at the drum
  centre); the normal is then rotated back to world space. Points from the wall towards free space.
  A non-finite `p_world` short-circuits to `(0.0, Vec2::ZERO)` rather than letting the numeric
  gradient manufacture a NaN from `-inf − -inf` (documented in detail on `sdf_world`) — this keeps a
  caller's typical `p += -d * normal` correction a safe no-op instead of injecting NaN.
- **Wall velocity.** `wall_velocity(p) = omega * perp(p) = omega * (-p.y, p.x)`, i.e. `omega x r` in
  2D, evaluated at the world-space point currently coincident with that wall material point
  (`Drum::wall_velocity`, `geometry.rs:112`).
- **Swept time of impact** (`Drum::toi_swept`, `geometry.rs:146`). Returns the first fraction
  `tau` in `[0,1]` of a sub-step at which a disc of radius `r`, whose centre moves on the straight
  world-space line `x0 -> x1` while the drum turns from `theta0` to `theta1`, touches the solid —
  seeing the *moving* wall, which a point query at one fixed angle (`sdf_world`) cannot (a lifter
  sweeping into a resting ball). It sphere-traces the static drum-frame SDF along the path
  expressed in the instantaneous drum frame, advancing

  ```
  d_tau = max( (d(tau) - r) / ( |x1 - x0| + |theta1 - theta0| * (|x0| + |x1 - x0|) ),  1e-6 )
  ```

  where the denominator is a rigorous upper bound on the speed of the point in the drum frame
  (`d(local)/dt = R(-theta(t)) (x1-x0) - theta' perp(local)`, `|x0 + t (x1-x0)| <= |x0| + |x1-x0|`),
  so, the SDF being 1-Lipschitz, the step never overshoots. It stops at the first `tau` with
  `d - r <= 0` and returns it; it returns `None` if the ball already touches or overlaps the solid
  at `tau = 0` (recovering an existing overlap is the discrete solver's job, and treating it as a
  crossing would freeze a ball that is separating), if there is no relative motion, or if the path
  stays clear; after 32 iterations without a clean root it conservatively reports contact (a false
  negative is the dangerous direction for a tunnelling guard). How the DEM uses it: §4.2 step 2.5.

---

## 3. Media / coarse-graining (`params.rs`)

### 3.1 Current defaults

`Params::default()` is the 2 mm ball in a 63 mm drum, a bench-scale configuration whose
`N_true ~ 244` needs no coarse-graining (values from the `Default` impls, `params.rs`):

| Group | Value (source) |
|---|---|
| Drum | `mill.diameter_m = 0.063`, 100 rpm (~59 % of `Nc = 42.3/sqrt(D) = 168.5 rpm`), counter-clockwise (`params.rs:62-72`) |
| Media | `ball_diameter_m = 0.002`, `fill_fraction` (J) `0.30`, `packing_fraction_2d = 0.82`, `density_kg_m3 = 6000` (YSZ-like), restitution ball–ball / ball–wall `0.7 / 0.5`, friction ball–ball / ball–wall `0.25 / 0.35` (static = kinetic), `friction_velocity_scale_m_s = 0.02`, `rolling_friction = 0.01` (`params.rs:200-226`) |
| Slurry | water: `density_kg_m3 = 1000`, `viscosity_pa_s = 0.001`, `fill_fraction = 0.35`, Newtonian, `wall_no_slip = ball_no_slip = 1`, `wettability = 0.6`, `surface_tension_n_m = 0.072` (`params.rs:350-373`) |
| Lifters | `count = 0` (smooth wall); if enabled: height 0.020 m, base width 0.030 m, top width 0.020 m (`params.rs:417-427`) |
| Simulation | `resolution = 40`, `substeps = 8` (requested), `pbf_iterations = 3`, `dem_iterations = 2`, `max_balls = 600`, `coarse_graining_mode = Off`, `coarse_graining_k = 1.0`, `time_scale = 1`, `frame_budget_ms = 12`, `seed = 1` (`params.rs:484-498`) |

The web quality presets (`web/src/params/presets.ts`: "Fast (default)" 300/25, "Balanced" 300/30,
"Accuracy" 600/50, "Impact statistics" 600/25 as `max_balls`/`resolution`) override
`max_balls` / `resolution` of `Params::default()` (which has `resolution = 40`); the
"Impact statistics" preset also sets `substeps = 32`, `dem_iterations = 12` and coarse-graining
`off`. At the shipped media size only `resolution` matters because nothing coarse-grains. Older
text in the dated history of §9 quotes the defaults of its day (50 Pa·s slurry, a 1 m drum,
Realtime 150/15, automatic coarse-graining); those are not current.

### 3.2 True population and coarse-graining

Balls are seeded from an *effective* media population (`Params::effective_media`,
`params.rs:582`). By default (`CoarseGrainingMode::Off`) that is the true population.

- **True disc count** (`Params::true_ball_count`, `params.rs:556`):
  `N_true = fill_fraction * packing_fraction_2d * drum_area / (pi * r_true^2)`
  where `drum_area = pi * radius_m^2` and `r_true = ball_diameter_m / 2`. `fill_fraction` (the
  conventional "J" ball-filling fraction) is the charge footprint *including voids*;
  `packing_fraction_2d` (default `0.82`, valid range `[0.5, 0.907]`) converts that footprint into a
  solid disc count — random close packing of equal discs is ~0.82, `pi / (2*sqrt(3)) ~= 0.907` is
  the hexagonal upper bound. At the defaults `N_true = 0.30 * 0.82 * 0.0315^2 / 0.001^2 ~ 244`.
- **Mode** (`CoarseGrainingMode`, `params.rs:35`; `simulation.coarse_graining_mode`):
  - `Off` (**default**): `k = 1`, `N_sim = N_true` regardless of `max_balls`;
  - `Auto`: `k = sqrt(N_true / max_balls)` if `N_true > max_balls`, else `k = 1`;
  - `Manual`: `k = max(1, simulation.coarse_graining_k)`, ignoring `max_balls`.
- **Effective media** (for `k > 1`):
  - `diameter_m = true_diameter_m * k` (preserves total footprint area, `N_sim * pi * r_eff^2 ~=
    N_true * pi * r_true^2`, i.e. the fill fraction the user set).
  - `density_kg_m3` is **left unchanged** at the true media density. Because balls are modelled as
    unit-depth discs (`ball_mass = rho * pi * r^2`, ss4 below), preserving the total solid
    footprint area at constant density automatically preserves total charge mass (`mass ~ r^2`, and
    `N_sim * r_eff^2 == N_true * r_true^2` by construction) — no compensating `1/k` density scaling
    is needed or applied (`effective_media_preserves_total_charge_mass` /
    `effective_media_coarse_grains_when_true_count_is_large` in `params.rs`'s test module assert it).
  - `ball_count = round(N_true / k^2)` (at least 1; `~= max_balls` in `Auto` mode).
  - When the resolved `k <= 1`, `scale_factor = 1.0` and the true media parameters pass through
    unchanged (`ball_count = round(N_true)`).
- This is a standard coarse-grained DEM approximation: it preserves bulk charge mass and footprint
  area, not the true interstitial void structure or single-collision statistics at the real
  particle size (and per-impact statistics and the power draw change measurably with `k`:
  docs/VERIFICATION.md, "What coarse-graining preserves"). `EffectiveMedia { true_diameter_m,
  diameter_m, density_kg_m3, ball_count, scale_factor }` (`params.rs:709`) carries both the true and
  simulated values so the approximation is never silent (surfaced in the UI's derived-values panel
  and in `metrics::Metrics`'s `true_ball_count`/`simulated_ball_count`/`coarse_graining_factor`/
  `effective_ball_diameter_m`). The true radius survives in `Balls::true_radius` (`dem.rs:308`) and
  is used only by the form-drag rate (§6.2).

### 3.3 Automatic accuracy scaling

Two derived quantities are raised automatically (never silently: the parameters panel shows
"Effective slurry resolution" and "Effective sub-steps"). With `d_eff = effective_media().diameter_m`:

- **Fluid resolution** (`Params::effective_fluid_resolution`, `params.rs:648`). The fluid lattice
  spacing is `dx = R / resolution`. To keep one fluid particle from spanning several balls (which
  saturated the coupling clamps and injected energy) the resolution is raised so that
  `dx <= d_eff`:

  ```
  res_eff = min(200, max(resolution, ceil(R / d_eff)))
  ```

  At the defaults `R / d_eff = 15.75`, so `res_eff = max(40, 16) = 40` (the request stands).
- **Sub-step count** (`Params::effective_substeps`, `params.rs:690`). The XPBD tunnelling
  criterion `max|v| dt / (2r) < 1` (§4.2, §8) worsens as `d_eff` shrinks. With a deliberately
  generous reference impact speed `v_ref = sqrt(2 g D)` (free fall across the drum diameter `D`) and
  `dt = 1 / (60 n)`, requiring `v_ref dt <= TARGET_RATIO * d_eff` with `TARGET_RATIO = 0.5`
  (`params.rs:692`, half the `< 1` criterion) gives

  ```
  n_min = ceil( v_ref / (60 * 0.5 * d_eff) ),   n_s = max( substeps, min(n_min, 16) )
  ```

  The auto-raise stops at 16; an explicit request above 16 (the validated range is 1–64) is honoured
  as is (`params.rs:702`). At the defaults `v_ref = sqrt(2 * 9.81 * 0.063) = 1.11 m/s`,
  `n_min = 19`, so `n_s = 16` and `dt = 1/960 s`. The legacy 1 m drum with a 63 mm ball needs no
  raise (`effective_substeps_matches_requested_for_the_legacy_1m_drum`). Clamping at 16 means a
  pathological combination (a very small ball with `max_balls` far above `N_true`) cannot make a
  frame unboundedly expensive; in that residual regime the displacement ratio may still approach or
  exceed 1, a disclosed limit of the auto-raise (§9).

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
2. **Broad-phase** (`dem.rs:597`-`dem.rs:633`). A `UniformGrid` is built over the predicted ball
   positions with cell size

   ```
   s = max( 2r * 1.05,  2r + v_max * dt + n_it * MAX_RECOVERY_FRACTION * 2r )
   ```

   (`dem.rs:626`-`dem.rs:630`; `v_max` is the largest post-predict ball speed, `n_it =
   dem_iterations.max(1)`): a fixed 5% pad, or the largest predicted per-sub-step displacement
   *plus the motion step 3 itself can add* (each of the `n_it` passes can move a ball by up to
   `MAX_RECOVERY_FRACTION * 2r`), whichever is larger. `for_each_candidate_pair` enumerates
   ball-ball candidate pairs once, reused across all solver iterations and steps 5–7 of this
   sub-step. The solve-margin term is load-bearing: without it a third ball just outside the
   predict-only margin can be pushed into a genuine overlap during step 3 and never be recorded as
   a contact (ablation: a −22.3 J full-system energy violation, history in §9). Ball-wall contact
   is checked directly per ball (no broad-phase needed: one `Drum::sdf_world` query per ball). The
   margin makes the *discrete* ball–ball search robust but is not a swept ball–ball check (§9).
2.5. **Swept wall/lifter clamp** (`dem.rs:635`-`dem.rs:706`). For each ball, the end-of-step
   constraint `C_end = d(x_i; theta_D + omega dt) - r` is evaluated. Only if the predicted overlap
   is deeper than one bounded-recovery pass can absorb, `C_end < -MAX_RECOVERY_FRACTION * 2r`
   (`dem.rs:679`), `Drum::toi_swept` (§2) is called; if it reports a time of impact `tau`, the
   predicted advance is truncated to

   ```
   x_i <- x0_i + (x_i - x0_i) * min( tau + CCD_CONTACT_SKIN_FRACTION * r / |x_i - x0_i|,  1 )
   ```

   `CCD_CONTACT_SKIN_FRACTION = 0.01` (`dem.rs:72`) advances a small physical depth past the exact
   crossing so the resulting state is a genuine small overlap that steps 5–7 see (an exact touch
   has `C = 0`, would record no `lambda_n`, and would silently drop friction, restitution and
   rolling resistance). Swept detection is applied to the **wall and lifters only**: a ball–ball
   version froze motion in dense beds and made packed charges worse, and an *unconditional* wall
   version repeatedly truncated ordinary tangential sliding along the curved wall (the one-sub-step
   chord dips inside the true arc almost every sub-step) — hence the depth gate (comments at
   `dem.rs:635`-`dem.rs:673`).
3. **Solve non-penetration** (`dem_iterations` Gauss-Seidel passes, zero compliance, with a
   **bounded recovery rate**):
   - Ball-ball: `C = |x_i - x_j| - (r_i + r_j)`; if `C < 0`, `C_eff = max(C, -MAX_RECOVERY_FRACTION
     * 2r)` (`MAX_RECOVERY_FRACTION = 0.2`), `d_lambda = -C_eff / (w_i + w_j)` (`w = 1/mass`), apply
     `x_i += w_i * d_lambda * n_hat`, `x_j -= w_j * d_lambda * n_hat`, accumulate `lambda_n`
     (computed from the *capped* `C_eff`, so a partially-recovered contact correctly carries less
     normal force this sub-step) for that `(i, j)` pair in `ContactBook`.
   - Ball-wall/lifter: `C = drum.sdf_world(x_i, drum_angle_next) - r_i`; same `C_eff` cap; if
     `C < 0`, `d_lambda = -C_eff / w_i` (wall has `w_wall = 0`, infinite mass), apply
     `x_i += w_i * d_lambda * n_hat`, accumulate `lambda_n` for ball `i`. `drum_angle_next =
     drum_angle + drum.omega * dt` -- the drum's angle at the **end** of this sub-step, not its
     start: balls have already been predicted to `t + dt` (step 1), so the boundary they are
     tested/projected against must be where the (possibly rotating, lifter-bearing) wall actually
     is at `t + dt`. Steps 5/6's own `drum.sdf_world`/`wall_velocity` calls use the same
     `drum_angle_next` for the same reason; an earlier version evaluated all of these against the
     sub-step's *start* angle, a one-sub-step kinematic lag an external review flagged.
   - **Why the cap.** Recovering a very deep overlap in full, in one iteration, hands step 4 a
     position delta that becomes an unphysically large separation velocity (`Δx/dt`) — the DEM
     half of the fluidised-charge energy-injection bug (ss9/git history). At the current default
     `dem_iterations = 2` the cap still allows recovering up to `0.4` diameters of overlap per
     sub-step, so ordinary small overlaps are unaffected; it only throttles the pathological case.
4. **Reconstruct velocities**: `v = (x - x0) / dt`, `omega = angle_diff(theta, theta0) / dt`.
5. **Friction** (Coulomb-clamped position correction, iterating contacts in a deterministic sorted
   order, with `v`/`omega` kept in sync *per contact* rather than reconstructed once at the end):
   - Ball-ball: `v_t = (v_i - v_j).dot(t_hat) - r*omega_i - r*omega_j`; the raw correction that
     would fully cancel `v_t` this sub-step is `raw = -v_t * dt / w_sum_t` (`w_sum_t = w_i + w_j +
     r^2*w_rot_i + r^2*w_rot_j`), clamped to `[-mu_eff*lambda_n, +mu_eff*lambda_n]`
     (`mu_eff` below, from `media.friction_ball_ball` / `friction_ball_ball_static`) and applied
     as a position correction (linear on `x`, rotational on
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
     the clamp coefficient is `mu_eff` built from `media.friction_ball_wall` /
     `friction_ball_wall_static`. `wall_velocity` is sampled at the
     **contact point** (`x_i - r*n_hat`), not the ball's centre — the wall is rigid but not a point,
     so its velocity varies across the ball's own radius, and the contact point is the only
     location where "zero relative tangential slip" is meaningful; sampling at the centre instead
     (an earlier version's behaviour, an external review flagged) undercounted the wall's tangential
     speed by `omega_ball * r`, reading a perfectly rolling ball as slipping. The tangential
     *position* multiplier `d_lambda_t` (not itself an impulse — see ss6's impulse-vs-force
     distinction) divided by `dt` gives the actual tangential impulse the wall delivered; dotted
     with the wall's (contact-point) velocity, that gives the work done against wall friction this
     contact. An earlier version omitted the `/ dt` here and under-reported power draw by a factor
     of `dt` (a factor 240 at the then-default sub-step) relative to `dissipated_energy_j`.
     **`wall_work_j`** additionally sums the wall's *normal* impulse dotted with its own
     contact-point velocity (`(lambda_n/dt) * n_hat.dot(v_wall)`): zero for a smooth cylindrical
     wall, whose normal is always radial while its velocity is purely tangential there, but nonzero
     for a lifter face, whose normal has a large circumferential component and does real work
     lifting the charge — a contribution an earlier version omitted entirely, materially
     under-reporting power draw/torque on a lifters-enabled mill (an external review finding).
   - Velocities are reconstructed a second time from the total position/orientation change once
     more (redundant with the per-contact sync above in exact arithmetic; kept as a cheap
     self-healing resync against float summation-order drift, same as step 4).
   - **Static/kinetic blend** (`effective_friction`, `dem.rs:128`). The Coulomb coefficient used in
     both clamps depends smoothly on the contact's relative tangential speed `v_t`:

     ```
     mu_eff(v_t) = mu_k + (mu_s - mu_k) * exp( -(v_t / v0)^2 ),   v0 = friction_velocity_scale_m_s = 0.02 m/s
     ```

     an even, everywhere-differentiable function (`mu_s` as `v_t -> 0`, `mu_k` for `|v_t| >> v0`),
     chosen over a hard Karnopp-style stick/slip switch to avoid a velocity discontinuity that makes
     a contact chatter between "just stuck" and "just slipping" in this fixed-sub-step solve. At the
     shipped defaults `mu_s = mu_k` (ball–ball `0.25`, ball–wall `0.35`; `params.rs:217-222`), and
     the function returns `mu_k` exactly via a fast path when `mu_s <= mu_k`, so the default
     behaviour is the single-coefficient Coulomb friction (bit-identical to before the split).
6. **Restitution** (bidirectional, using the *pre-solve* approach velocity `v_pre` to choose the
   target, applied to *every* contact regardless of approach speed):
   - Ball-ball: `v_n_pre = (v_pre_i - v_pre_j).dot(n_hat)`. If `v_n_pre` is more negative than
     `-RESTITUTION_VELOCITY_THRESHOLD` (`0.02` m/s, an absolute constant, `dem.rs:32`, used by
     the solver in every build: `dem.rs:919`, `dem.rs:981`) it's a fresh impact — restitution `e =
     media.restitution_ball_ball` applies and this contact is counted in
     `DemStepStats::collision_count`, with 2D-disc impact energy `E = 0.5 * (mass/2) * v_n_pre^2`
     (reduced mass `mass/2` for two equal masses) binned into the internal 12-bin
     `impact_energy_histogram` (these two DEM-internal counters are *not* what the simulator
     reports; see the end of this section and §8). Otherwise it's a
     resting/sliding contact and `e = 0`. Either way the target relative normal velocity is
     `target = max(0, -e * v_n_pre, v_n_pre)`, and the *actual* current relative normal velocity
     `v_n_now` is driven to that target by an impulse split by inverse mass — **in both
     directions**, not only when `v_n_now` falls short of `target`. The third term in `target`
     (`max(..., v_n_pre)`) floors it at the pair's own pre-solve separation speed: without it, a
     pair that entered this step *already separating* faster than `-e * v_n_pre` (e.g. from step
     3's depenetration) would be decelerated below the speed it already had — an artificial
     attraction between dry rigid discs, violating the one-sided normal-force (Signorini)
     condition, which an external review flagged. This floor is inert for an approaching pair
     (`v_n_pre <= 0`, so `v_n_pre` itself is `<= 0` and never binds), so it does not reopen the
     energy-injection bug the bidirectional pass below exists to close.
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
   torque-coefficient model of rolling friction. **Known simplification** (ss9): this damps every
   ball's *world* `omega` toward zero rather than the *relative* spin at each contact, so it
   applies no reaction torque to a ball-ball contact's partner (angular momentum is not conserved)
   and, for a ball-wall contact, targets zero rather than `drum.omega` (fighting, rather than
   assisting, a ball correctly rolling without slipping on a rotating drum).
`DemStepStats` (per-substep diagnostics, all per-metre-of-mill-depth) returned by this function:
`wall_work_j`, `collision_count`, `impact_energy_histogram` (12 log-spaced bins from `1e-6` J to
`1.0` J, `impact_energy_bin_edges`, `dem.rs:75`-`dem.rs:93`), `dissipated_energy_j`,
`max_substep_displacement_over_diameter` (see ss9) and `coupling_work_j` (kinetic-energy change
caused by the applied fluid impulses, for the energy budget, §8).

**What is reported as collision rate and impact histogram.** The 12-bin 2D-disc histogram and
`collision_count` above are internal to the DEM step. The simulator's *reported* collision rate and
impact-energy histogram come from the impact log: `step_recorded` (`dem.rs:530`) appends one
`ImpactRecord { kind, ball, partner, normal_speed, slip_speed, restitution }` (`dem.rs:152`) per
impact whose approach speed exceeds the **size-relative** threshold
`restitution_velocity_threshold(r) = 0.1428 * sqrt(2 g r)` (`dem.rs:44`; `0.1428 sqrt(g d)` equals
the historical 0.02 m/s exactly at `d = 2 mm`), and `stress.rs` converts these to 3D statistics
(28 dimensionless log-bins, debounce of re-contacts, 3D sphere mass; §8). Recording only reads state
the solver already computed, so the trajectory is bit-identical with or without it. The solver
itself keeps the absolute `0.02 m/s` threshold for choosing the restitution target.

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
no-slip velocity blending, and an implicit (backward-Euler, conjugate-gradient) viscosity solve
with an optional regularised-Bingham yield stress, solved jointly with the ball velocities where
balls are near the fluid. There is no `s_corr` artificial pressure (removed; see §5.2 step 3).

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
3. **Density-constraint solve** (`pbf.rs:1241`-`pbf.rs:1334`; Jacobi-style passes — all deltas
   computed from current positions, then applied together; the pass count is adaptive, below):
   - `C_i = max(0, rho_i/rest_density - 1)` — **one-sided**: only over-dense particles are
     projected apart (`pbf.rs:1306`). Macklin & Müller's published PBF leaves `C_i` unclamped and
     suppresses the resulting tensile instability with an artificial-pressure term; here the clamp
     removes the instability at its source.
   - Gradients carry the particle mass: with `rho_i = m_f * sum_j W_ij`, `dC_i/dx_j = (m_f /
     rest_density) * grad_W_ij`, so `g_ij = (m_f/rest_density) * grad_W_ij` and `g_ii = sum_j g_ij`
     (`pbf.rs:1285`). Omitting the `m_f` factor (a bug until 2026-09-30) made the correction too weak
     by a factor of order `1/m_f`, negligible at a 1 m drum but ~1000x at the 63 mm default.
   - `lambda_i = -C_i / (sum_j |g_ij|^2 + |g_ii|^2 + EPSILON_RELAX)`,
     `EPSILON_RELAX = 200.0` (`pbf.rs:48`; CFM-style relaxation against divide-by-zero for sparse
     neighbourhoods; `pbf.rs:1309`).
   - `delta_p_i = (m_f / rest_density) * sum_j (lambda_i + lambda_j) * grad_W_ij` (`pbf.rs:1321`),
     applied directly to `x`.
   - **Adaptive pass count.** `simulation.pbf_iterations` (default 3) is a *floor*. After each pass
     the loop stops once at least that many passes have run and the largest compression error
     measured at the start of that pass, `max_i C_i`, is `<= COMPRESSION_CONVERGENCE_TOLERANCE =
     0.02` (`pbf.rs:61`); otherwise it continues up to `ADAPTIVE_MAX_EXTRA_ITERATIONS = 6` further
     passes (`pbf.rs:57`, `pbf.rs:1326`-`pbf.rs:1333`). A violent impact that leaves a large
     residual therefore costs extra passes only when needed.
   - History: an earlier version added an artificial-pressure term (`s_corr`) here, the standard
     PBF remedy for the tensile instability an *unclamped* density constraint produces; since this
     constraint is one-sided that instability does not arise, and
     `s_corr` was never found stable at this project's real-SI parameter scale (`rest_density`
     ~1000-2000 kg/m^3) in the time available, so it was removed entirely rather than kept as a
     disabled dead option. This solver's fluid-fluid attraction is now the explicit, separately
     tunable cohesion term (ss6.1a step 3.6a, `slurry.surface_tension_n_m`), not an artefact of
     this constraint.
4. **Ball overlap projection** (source step 3.5, `pbf.rs:1337`-`pbf.rs:1444`, coupling step 1 — see
   ss6.1) and **cohesion / adhesion** (source step 3.6, `pbf.rs:1447`-`pbf.rs:1631`, ss6.1a).
5. **Boundary projection** (source step 4, `pbf.rs:1634`; position only): for any particle with
   `drum.sdf_world(x, drum_angle_next) = d < 0`, `x += -d * normal`; marks the particle
   `touched_wall` for the no-slip blend below. Uses `drum_angle_next = drum_angle + drum.omega *
   dt` (the drum's angle at the **end** of this sub-step), not `drum_angle`, since particles have
   already been predicted to `t + dt` in step 1 — the same one-sub-step lag fix as `dem.rs`'s
   `DemState::step_with_external_forces` (ss4.2).
6. **Velocity reconstruction** (source step 5, `pbf.rs:1655`): `v = (x - x0) / dt` for every
   particle (`x0` = pre-predict position); the geometric excess of the overlap push is subtracted
   and the deferred cohesion/adhesion kicks are added.
7. **No-slip wall blending** (source step 6, `pbf.rs:1676`): for particles marked `touched_wall`,
   `v = (1-beta)*v + beta*v_wall` (`v_wall = drum.wall_velocity(x)`), `beta = slurry.wall_no_slip`
   clamped to `[0,1]`. The implied
   impulse this delivers to each such particle (`mass * (v_new - v_old)`), dotted with the wall's
   own velocity there, is summed into `FluidStepStats::wall_work_j` — the fluid-side counterpart of
   `dem::DemStepStats::wall_work_j`, folded into `Simulation::power_draw_w`/`torque_nm` (ss8) so a
   wet mill's motor load reflects viscous drag on the slurry, not only ball-wall friction.
8. **Ball form drag** (source step 6.5, `pbf.rs:1703`; coupling step 2 — see ss6.2).
9. **Buoyancy** (source step 6.6, `pbf.rs:2019`; coupling step 3 — see ss6.3).
10. **Implicit viscosity** (source step 7, `pbf.rs:2129`-`pbf.rs:2300`; `mu =
    slurry.viscosity_pa_s`; skipped entirely, at zero iteration cost, when `mu == 0`): density is
    recomputed once more at the final, post-boundary-projection positions; if the rheology is
    Bingham with a positive yield stress, the per-particle effective viscosity `mu_eff` replaces
    `mu` (ss5.3); `morris_weights` are built; then either the fluid-only
    `solve_implicit_viscosity` (no ball ghost points and no squeeze links this sub-step) or the
    joint ball–fluid `solve_viscosity_with_balls` solves the system in place (ss5.3, ss6.2, ss6.5).
    `mean_shear_rate` is computed from the same final positions/densities/neighbour lists.
11. **Fluid speed clamp** (source step 7.5, `pbf.rs:2310`; stability backstop): any particle speed
    above `v_max = FLUID_SPEED_SAFETY_FACTOR * (|omega|*radius_m + sqrt(4*|GRAVITY|*radius_m))`,
    `FLUID_SPEED_SAFETY_FACTOR = 5.0` (`pbf.rs:974`, `fluid_speed_clamp_v_max` at `pbf.rs:979`), is
    rescaled down to `v_max` — viscosity (step 10) can push a
    particle over the ceiling via a persistently-squeezed position correction (step 6 reconstructs
    velocity from position, so nothing else bounds the resulting speed until the wall's own
    projection does, far above the physically-attainable range). Deliberately not a CFL-style
    `h/dt` bound (that would be coupled to `simulation.resolution` and could clip real motion at
    high resolution instead of only blow-ups). Momentum removed here is *not* credited back to any
    ball.
12. **Coupling impulse clamp** (source step 8, `pbf.rs:2324`; see ss6.4).

### 5.3 Implicit viscosity detail

The viscous term `(mu/rho) * laplacian(v)` is discretised with the Morris (1997) SPH viscous
coefficients (`morris_weights`, `pbf.rs:312`):

```
c_ij = m_j * (mu_i + mu_j) / (rho_i * rho_j) * |x_ij . grad_W_ij| / (|x_ij|^2 + eta^2)
eta^2 = VISCOSITY_ETA_FACTOR * h^2,  VISCOSITY_ETA_FACTOR = 0.01      (pbf.rs:106)
```

(`mu_i` is the particle's effective viscosity: `mu` for a Newtonian slurry, `mu_eff` of the Bingham
law below otherwise; `x_ij = x_i - x_j`; `c_ij >= 0` and symmetric in `i, j`.)

**Central-force operator.** The coefficients are applied in a *central-force* form that acts along
the line of centres only (`laplacian_apply_into`, `pbf.rs:359`):

```
(L v)_i = sum_j CENTRAL_VISCOSITY_FACTOR * c_ij * e_ij * ( e_ij . (v_i - v_j) ),
          e_ij = (x_i - x_j)/|x_i - x_j|,    CENTRAL_VISCOSITY_FACTOR = 4      (pbf.rs:384)
```

The factor 4 is the 2D angular average (`<e_a e_b e_k e_l>`) that recovers the same continuum
Laplacian for incompressible flow as the plain Morris form `sum_j c_ij (v_i - v_j)`. The plain
form exerts a spurious force on a rigidly rotating fluid wherever the kernel sum is one-sided (free
surface, drum wall) — it damped solid-body rotation there (5x the torque-balance power at 1 Pa·s,
docs/VERIFICATION.md). The central form vanishes exactly for rigid rotation and conserves angular
momentum (`pbf::tests::viscous_operator_vanishes_for_rigid_rotation_even_at_a_free_edge`). `L` is
symmetric positive semi-definite since every `c_ij >= 0`.

**Implicit solve.** The backward-Euler system `(I + dt*L) v_new = v` is solved matrix-free by
conjugate gradient, operating directly on the `Vec2` field with the real inner product `<a,b> =
sum_i a_i . b_i` (mathematically equivalent to running scalar CG on the x- and y-components
independently).

- Fluid-only path, `solve_implicit_viscosity` (`pbf.rs:423`): plain CG, relative-residual stopping
  tolerance `VISCOSITY_CG_TOLERANCE = 1e-3` (`|r| <= tol * |b|`, `pbf.rs:96`), warm-started from
  `v` itself, iteration ceiling `VISCOSITY_CG_MAX_ITERS = 50` (`pbf.rs:100`). Used when no ball
  ghost points and no squeeze links exist this sub-step.
- Joint path, `solve_viscosity_with_balls` (`pbf.rs:641`): **Jacobi-preconditioned** CG (inverse of
  the operator's diagonal, `viscosity_diag` `pbf.rs:388` plus the boundary and link diagonals),
  residual measured in the preconditioner norm against the same `1e-3`, ceiling `4 x 50 = 200`
  iterations (the ball unknowns couple distant fluid through a stiff rigid-body mode), and
  **warm-started** with the incoming velocity plus the previous sub-step's viscous correction
  (`viscous_guess`; the right-hand side still uses the incoming velocity alone). See ss6.2.

`A = I + dt*L` is SPD whenever `dt > 0`, so CG converges in exact arithmetic; the iteration
ceilings are generous practical bounds, not tuned-tight ones. This replaced an earlier explicit
XSPH mixing scheme whose effective kinematic viscosity saturated around `~h^2/dt` regardless of the
input value, so it could not represent a genuinely low-viscosity (near-water) slurry; the implicit
solve is unconditionally stable for any `viscosity_pa_s >= 0` at the fixed sub-step and uses the
value directly in Pa·s.

**Joint ball–fluid solve and ball-surface no-slip boundary.** Implemented and used by the shipped
simulator; its equations (unknowns `v_i`, `V_b`, `W_b`; `BoundaryGroup`, `BALL_BC_GAIN = 0.5`) are
in ss6.2, because it is also how viscous drag on a ball is obtained.

**Bingham yield stress** (`bingham_viscosity`, `pbf.rs:854`; `shear_rates`, `pbf.rs:821`). When
`slurry.rheology == Bingham` and `slurry.yield_stress_pa = tau_y > 0` (`pbf.rs:2148`), the
per-particle shear rate is the SPH velocity-gradient estimate

```
grad_v_i = sum_j (m_f/rho_j) (v_j - v_i) (x) grad_W_ij,   D = sym(grad_v),   gamma_dot = sqrt(2 D:D)
```

evaluated from the velocities entering the viscosity step, and the effective viscosity is the
Papanastasiou-regularised Bingham law (exponential regularisation of the yield surface)

```
mu_eff(gamma_dot) = mu + tau_y (1 - exp(-gamma_dot/gamma0)) / gamma_dot,
gamma0 = tau_y / (BINGHAM_PLASTIC_RATIO * mu),   BINGHAM_PLASTIC_RATIO = 100   (pbf.rs:490)
```

with the limit `gamma_dot -> 0` equal to `mu (1 + BINGHAM_PLASTIC_RATIO)`, i.e. the unyielded
viscosity is capped at `101 mu`. `mu_eff` replaces `mu` in `c_ij` and in the ball-surface boundary
weights. Newtonian is the default (`yield_stress_pa = 0`). Only the Bingham case is implemented in
the PBF solver (no Herschel–Bulkley power law; that exists only in the experimental grid track), and
against the grid reference its response is only qualitatively right (docs/VERIFICATION.md, "PBF:
Bingham and coarse-graining k").

`mean_shear_rate` (`pbf.rs:800`): the mean over particles of the same `gamma_dot`, surfaced via
`FluidStepStats` for `metrics::Metrics::mean_shear_rate_per_s`.

`nu = mu/rho` uses each particle's own *measured* SPH density, not `slurry.density_kg_m3` directly
— these agree closely in a well-relaxed interior but diverge near a low-density boundary (free
surface, sparse neighbourhood), where the resulting locally elevated `nu` is itself physically
reasonable (a thin/rarefied film resists shearing more, relative to its own mass, than the bulk).

---

## 6. Two-way ball<->fluid coupling (`coupling.rs`, `pbf.rs` steps 3.5/3.6/6.5/6.6/7/8)

`coupling::step_recorded` (`coupling.rs:71`) orchestrates one shared sub-step: the fluid solves first (`step_coupled`), using
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
OVERLAP_PUSH_MARGIN * h` (`OVERLAP_PUSH_MARGIN = 0.05`, previously `0.25`) of a ball centre is
pushed out along the connecting normal by `push_mag = min(contact_radius - dist, 0.5 *
balls.radius)` — capped at half the ball's radius so a particle deeply embedded in the contact zone
cannot produce an outsized single-substep *position* correction. The grid used to find a ball's
nearby fluid particles here is sized to cover `2 * balls.radius` or `contact_radius`, whichever is
larger (regression test: `pbf::tests::coupling_finds_ball_neighbours_beyond_two_radii`).

**The momentum exchanged is bounded separately from the position correction.** The full push
implies a velocity of `push_mag / dt` once step 5 reconstructs velocity from position — unbounded
in momentum as `dt` shrinks or as particles pile up against a ball, since the position cap alone
doesn't limit it. Instead, the actual momentum exchanged is (`pbf.rs:1404`-`pbf.rs:1431`)

```
impulse_mag = min( m_f * push_mag / dt,  m_red * arrest_speed ),   m_red = m_f * m_b / (m_f + m_b)
arrest_speed = max(0, -(v_i - v_surface) . n_hat)
```

where `arrest_speed` is how much of the particle's velocity *toward* the ball's surface
(`v_surface = balls.v + balls.omega x lever`) this contact actually needs to arrest — zero for a
particle already moving away, or at rest relative to the ball. The bound uses the **reduced mass**
`m_red` of the (fluid particle, ball) pair rather than the fluid particle's mass alone:
`m_f * arrest_speed` is the impulse that stops the particle against an immovable ball, correct only
while the ball is much heavier than one fluid particle. When coarse-graining is off and the true
ball is small (a 2 mm disc is lighter than a fluid particle at a coarse fluid resolution), a single
fluid particle can outweigh the ball and the one-sided `m_f` bound handed the ball an oversized
impulse; `m_red` reduces to `m_f` when `m_b >> m_f` and caps the ball's own reaction near
`m_b * arrest_speed` otherwise. The ball's Newton's-third-law reaction is `-impulse_mag * n_hat` (angular component via
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

Removing this step's ball-side reaction entirely (treating a ball as a purely passive boundary,
like the drum wall) was tried during the 2026-09-22 low-fill investigation (ss9) and reverted: it
measurably regressed `metrics::tests::compression_error_stays_bounded_under_violent_lifter_
cataracting` (worst compression error 0.123 vs. the previously-passing ~0.091 bound) — this step's
direct per-contact reaction does real stabilising work under violent cataracting that ss6.2/6.3's
differently-time-scaled exchange does not replace.

### 6.1a Fluid-fluid cohesion and fluid-boundary adhesion (pbf.rs step 3.6)

Replaced entirely (2026, "replace the fluid coupling with the Akinci formulation" milestone) by an
Akinci-style pairwise-force pair, adapted to 2D: fluid-fluid **cohesion** (step 3.6a,
`slurry.surface_tension_n_m`) and fluid-boundary **adhesion** (step 3.6b, `slurry.wettability`).
Together these give the slurry a genuine contact angle -- cohesion holds a film together against
gravity, adhesion pulls it onto media, and which one wins where is an emergent competition between
two symmetric pairwise forces, not the single ball<->fluid-only knob (no fluid-fluid counterpart)
the previous shell-based mechanism was. Source: Akinci, Akinci & Teschner, "Versatile Surface
Tension and Adhesion for SPH Fluids", 2013.

**Kernels.** Both terms share `cohesion_kernel(r, h)`, the paper's eq. 2 cubic shape (zero at
`r = 0` and `r = h`, peaking near `r ~= 0.6h`, slightly negative for `r < h/2` so two very close
particles get a small locally-repulsive contribution instead of collapsing to coincidence),
independently re-normalised for 2D (`2*pi * integral_0^h kernel(r) r dr == 1`, matching this
module's other kernels' convention) rather than reusing the paper's 3D coefficient. Adhesion reuses
just the outer half of that same shape (`r` in `(h/2, h)`), separately re-normalised over that
restricted domain, instead of the paper's own fractional-power adhesion formula (`0.007 h^-3.25 *
(-4r^2/h + 6r - 2h)^0.25`) -- that formula's non-integer exponents have no clean closed-form 2D
re-derivation, and this kernel needs no property the simpler cubic shape lacks for this solver's
purpose.

**3.6a Cohesion.** For each fluid particle, over its existing fluid-fluid neighbour list (step 2,
no new broad-phase pass):

```
pull_i        = sum_j (particle_mass * cohesion_kernel(r_ij, h) / rest_density) * n_hat_ij
                                                          (n_hat_ij points from i toward j)
cohesion_accel = (surface_tension_n_m / REFERENCE_SURFACE_TENSION_N_M)
                 * COHESION_ACCEL_FACTOR * |GRAVITY|      (REFERENCE = 0.072 N/m, real water;
                                                            COHESION_ACCEL_FACTOR = 2.0)
v_i          += -cohesion_accel * pull_i * dt
```

`particle_mass * kernel / rest_density` is dimensionless (kernel is normalised to `1/area`, so the
product is a density, and dividing by `rest_density` cancels units) -- roughly 1 for a particle deep
in bulk fluid and tapering to 0 near a free surface, which is what makes `cohesion_accel` a genuine
peak-acceleration scale rather than needing its own separate normalisation. A fully interior
particle has a (near-)symmetric neighbourhood, so the pairwise pulls cancel and it feels no net
force; only a free-surface/thin-film particle, with an asymmetric neighbourhood, feels a net inward
pull -- ordinary surface tension. `surface_tension_n_m` is a **calibrated proportionality to the
physical value, not a first-principles unit conversion** -- the mapping from an SPH cohesion
coefficient to a macroscopic N/m value is inherently resolution-dependent (a coarser or finer fluid
lattice needs a different coefficient for the same emergent surface tension), a property of the
method itself, not a shortcut this project's implementation takes. Regression test:
`pbf::tests::cohesion_pulls_two_isolated_fluid_particles_together_only_when_surface_tension_is_positive`.

**3.6b Adhesion.** Each ball's circumference is sampled into boundary particles
(`ball_boundary_particles`: roughly `dx` arc-length spacing, at least 6 points, recomputed fresh
every sub-step from the ball's current centre -- cheap, no state carried across sub-steps). The
kernel's support radius is capped, `h_adhesion = h.min(2 * balls.radius)`, regardless of the
fluid's own `h`: at this project's coupling resolution `h` can be several times the ball's own
radius, and using it directly (as the shell mechanism this replaces did, before its own separate
ball-radius-relative cap) let a distant fluid particle read as a floating clump rather than a
clinging film -- see below for that mechanism's own history.

Bounded per ball (not per contributing particle), mirroring ss6.2's drag closure for the identical
reason (an earlier version applied the peak acceleration independently to every particle in range
and let the ball absorb the sum, amplifying its response by the contributing-particle count):

```
weight_sum          = sum over (fluid j, boundary point k of ball b) of
                       (particle_mass * adhesion_kernel(r_jk, h_adhesion) / rest_density)
fluid_pull_dir_sum   = weighted sum of n_hat_jk (points from j toward k)
ball_pull_dir        = -normalize(fluid_pull_dir_sum)     (opposite: toward the fluid, not the
                                                             ball's own surface)
coverage             = min(1, weight_sum)
m_contrib            = particle_mass * weight_sum
attach               = m_contrib / (balls.mass + m_contrib)   <= 1, -> 0 as m_contrib -> 0
adhesion_accel        = wettability * ADHESION_ACCEL_FACTOR * |GRAVITY|   (ADHESION_ACCEL_FACTOR = 2.0)
dv_b                  = ball_pull_dir * (adhesion_accel * coverage * dt * attach)
```

`coupling.impulses[b] += balls.mass * dv_b` (angular term via a single mean `lever = ball_pull_dir *
balls.radius` across every contributing boundary point, the same level of aggregate-lever
simplification the drag closure and this mechanism's predecessor both already used); each
contributing fluid particle `j` absorbs `dv_j = -(balls.mass * frac_j / mass) * dv_b` (`frac_j` its
share of `weight_sum`), so the fluid side's total momentum change is exactly `-balls.mass * dv_b`
regardless of how many boundary points or fluid particles were involved. Balls are processed
sequentially (Gauss-Seidel), same reason as ss6.2: a fluid particle can be near more than one ball's
boundary at once. Both cohesion and adhesion are pure velocity kicks (`accel * dt`, no position
change), so -- like ss6.1's push excess -- they are accumulated separately and added back into
`self.v` right after step 5's `v = (x - x0) / dt` reconstruction, which would otherwise silently
discard them. Regression tests:
`coupling::tests::adhesion_pulls_a_ball_and_nearby_fluid_together_only_when_wettability_is_positive`,
`coupling::tests::adhesion_range_is_capped_at_twice_the_ball_radius_regardless_of_fluid_resolution`.

**History: the shell-based predecessor this replaced.** Before this pass, `slurry.wettability`
pulled fluid in a purely geometric shell `(shell_start, adhesion_radius]` beyond `contact_radius`
(ss6.1's overlap-push standoff), gated by each contributing particle's own SPH density
(`adhesion_wetness`) reading "wet" versus "an isolated airborne droplet". That gate existed to fix
a review-reported defect: ball+slurry clumps floating indefinitely through the air, caused by (1)
the shell reaching roughly 2.6 ball radii out from the surface at this project's typical coupling
resolution, wide enough that a distant fluid particle read as a separate floating clump rather than
a clinging film, and (2) nothing testing whether there was any actual bulk liquid there, so a
single stray droplet pulled at full strength identical to one touching a real pool. The
density-based wetness gate suppressed exactly the particles that would form a genuine thin wetting
film (necessarily low-density, being sub-resolution relative to `h`), and had no fluid-fluid
cohesion counterpart to hold a film together once pulled on -- a structural limitation the current
formulation's real contact angle and real geometric range cap (rather than a shell width tied to
the fluid's own resolution) do not share.

### 6.2 Ball drag: viscous (implicit no-slip boundary, step 7) and form (pbf.rs step 6.5)

Drag on a ball has two parts in the current closure, with two different mechanisms. The Stokes-regime
relaxation that used to carry the viscous part is **disabled** (`inv_tau_stokes = 0.0`,
`pbf.rs:1807`); its derivation and the investigations built on it are kept, clearly marked, in
§6.2.1.

**Viscous drag: ball-surface no-slip boundary in the joint implicit solve** (`pbf.rs` step 7,
`solve_viscosity_with_balls`, `pbf.rs:641`; `BoundaryGroup`, `pbf.rs:498`). The ball surface is a
no-slip (Dirichlet) boundary of the viscosity operator, and the fluid and ball velocities are
solved for *together*. Ghost points `k` on the circumference of ball `b` (spacing `dx = h/2`, at
least 6 points, `ball_boundary_particles`, `pbf.rs:190`) carry the rigid-body surface velocity
`u_k = V_b + W_b perp(l_k)`, `l_k` the offset from the ball centre. Each fluid particle `i` within
`h` of a ghost point couples with a Morris-type, isotropic (not central) weight (`pbf.rs:2166`-
`pbf.rs:2230`)

```
c_ik = G * m_f * 2 mu_i / (rho_i * rho0) * |delta_ik . grad_W_ik| / (|delta_ik|^2 + eta^2),
G = BALL_BC_GAIN = 0.5                                                      (pbf.rs:487)
```

(ghost points have mass `m_f` and density `rho0`; `mu_i` is `mu` or the Bingham `mu_eff`; the gain
`G` compensates for a single ghost layer being thinner than a true wall half-space and was
calibrated on the grid reference, docs/VERIFICATION.md "PBF ball no-slip boundary"). All ghost
points of one (fluid, ball) pair are merged into three sums, `c = sum c_ik`, `sum c_ik l_k`,
`sum c_ik |l_k|^2`. The unknowns are `v_i`, `V_b`, `W_b`, and the symmetric positive definite system
is

```
m_f v_i  + dt m_f (L v)_i + dt m_f sum_k c_ik (v_i - u_k)                          = m_f v*_i
m_b V_b  - dt m_f sum_{(i,k) in b} c_ik (v_i - u_k) + dt sum_{links of b} c n n^T (V_b - V_b') = m_b V*_b
I_b W_b  - dt m_f sum_{(i,k) in b} c_ik l_k x (v_i - u_k)                          = I_b W*_b
```

where `v*`, `V*`, `W*` are the values after the earlier steps of the sub-step (overlap push, form
drag, buoyancy: `ball_v0 = v + impulse/m_b`, `pbf.rs:2237`), so the solve does not double-count a
relative velocity those steps already relaxed. The second equation's last term is the squeeze-film
link of §6.5. The system is solved by Jacobi-preconditioned CG (tolerance `1e-3`, up to `4 x 50`
iterations, warm-started with the previous sub-step's viscous correction; §5.3). The ball's
impulse is exactly what the wall term took from the fluid: `impulse += m_b (V_b - V_b0)`,
`angular_impulse += I_b (W_b - W_b0)` (`pbf.rs:2290`-`pbf.rs:2298`).

*Conservation (what is and is not exact).* A uniform translation of all unknowns is a null mode of
every coupling term (each depends only on a velocity *difference*), so for the converged solution
the total linear momentum `m_f sum v_i + m_b sum V_b` equals its value before the solve; the ball
receives exactly the linear impulse the wall term removed from the fluid, and the torque
`l_k x F_k` of the same forces at the ghost points (so each ball's own angular momentum balance is
consistent by construction). With the CG tolerance the violation is bounded by the projection of
the final residual on uniform translation. The *total* angular momentum about the drum origin is
conserved only up to the offset between the fluid particle and the ghost point where the pair
force acts (the force on the fluid acts at `x_i`, the matching reaction on the ball is referred to
`x_k`, giving a residual `sum (x_i - x_k) x F`); the code comment on `solve_viscosity_with_balls`
(`pbf.rs:635`) says the same ("total angular momentum only up to the ghost-point offset"). Treating
the ball explicitly instead of in this joint solve was unstable once the slurry acts as a rigid
lump (`mu` of about 10 Pa·s and above).

**Form drag: centre-of-mass relaxation** (`pbf.rs` step 6.5, `pbf.rs:1703`-`pbf.rs:2013`). Each
ball relaxes toward the mean velocity of its locally entrained fluid mass, a closure bounded by
construction and not a per-particle blend. Every fluid particle `i` within `h_c = balls.radius + h`
of a ball centre contributes with a smooth Poly6 taper weight `phi_i = poly6(|r_i|^2, h_c) /
poly6(0, h_c)` (`r_i = x_i - x_b`), giving an entrained mass `m_ent = sum_i m_f phi_i` and a
weighted mean fluid velocity `v_bar = (sum_i m_f phi_i v_i) / m_ent`:

```
v_rel      = v_bar - v_b
1/tau_form = (1/2) * rho0 * C_D * (2 r_true) * |v_rel| / m_true     (0 if |v_rel| <= 1e-6)
beta       = slurry.ball_no_slip * (1 - exp(-dt / tau_form))
a_lin      = beta * m_ent / (m_b + m_ent)                           <= 1 for every dt
dv_b       = a_lin * v_rel
```

with `C_D = BALL_FORM_DRAG_COEFFICIENT = 1.0` (`pbf.rs:116`, a circular cylinder's high-`Re` drag
coefficient), `rho0 = rest_density` and the **true** (un-coarse-grained) particle's `r_true =
Balls::true_radius` (`dem.rs:308`), `m_true = rho_ball pi r_true^2`, `rho_ball = m_b / (pi
r_b^2)` (the material density, recovered exactly since coarse-graining scales mass and area
together). The rate is `F = (1/2) rho0 C_D d |v| v` in the usual convention: the factor `1/2` was
added on 2026-10-19 (`pbf.rs:108`-`pbf.rs:116`, `pbf.rs:1872`-`pbf.rs:1880`); before it the closure
acted as `C_D = 2`. Everything else in the step — the entrainment footprint `h_c`, `m_ent`, the
`a_lin` mass ratio and all momentum bookkeeping — uses the coarse-grained `m_b`, `r_b`, because
`m_ent` is sampled over the coarse ball's own footprint and the ratio `m_ent/(m_b + m_ent)` is
already scale-invariant (both grow with `k^2`); this keeps `simulation.max_balls` and the
coarse-graining mode from changing the drag *rate* (history: §6.2.1). The step runs only if
`mu > 1e-6` (a slurry of zero viscosity has no ball drag at all) and `ball_no_slip > 0`.

The reaction is distributed back across the same weighted fluid neighbours (`dv_i = -(m_b / m_ent)
dv_b phi_i`), so linear momentum is conserved exactly (`sum_i m_f dv_i == -m_b dv_b`). The angular
part mirrors this: `i_ent = sum_i m_f phi_i |r_i|^2`, `omega_bar = sum_i m_f phi_i r_i x v_i /
i_ent`, `a_rot = beta i_ent / (I_b + i_ent)`, a candidate `dw_b = a_rot (omega_bar - omega_b)`, and
a tangential reaction field `c_rot phi_i (perp(r_i) - p_bar)`, `p_bar = perp(r_bar)`, with zero net
linear momentum, plus a `m_b cross2(r_bar, dv_b)` term in the ball's angular impulse — the reaction
torque from spreading the linear reaction over an off-centre weighted mean position `r_bar`, needed
for exact angular-momentum conservation. The coefficient is `c_rot = -I_b dw_b / s`,
`s = sum_i m_f phi_i^2 cross2(r_i, perp(r_i) - p_bar)`, a geometric normalizer that can be
comparable to its positive-definite part `sum_i m_f phi_i^2 |r_i|^2` or many orders of magnitude
smaller whenever few fluid particles are entrained (2–4 per ball is routine at the shipped
resolutions). Two guards (`pbf.rs:1934`-`pbf.rs:1989`):

1. the rotational exchange (both `dw_b` and its fluid-side redistribution) is **skipped** for that
   ball this sub-step if `|s| <= 0.1 * sum_i m_f phi_i^2 |r_i|^2` — a *relative* threshold, since a
   fixed absolute epsilon (`1e-9`) let a modest torque be divided by an almost-cancelled
   denominator and amplified by orders of magnitude;
2. `c_rot` is **capped** so no fluid particle's velocity kick exceeds `0.2 * v_max`, `v_max =
   fluid_speed_clamp_v_max(drum)` (§5.2 item 11): `|c_rot| <= 0.2 v_max / (2 h_c)`, and `dw_b` is then
   re-derived from the applied `c_rot` (`dw_b = -c_rot s / I_b`) so the ball's torque and the
   fluid's reaction stay exactly consistent even when the cap engages.

`a_lin <= 1` by construction, so the ball can never be driven past `v_bar` regardless of `dt` or
the entrained/ball mass ratio (`pbf::tests::ball_drag_never_overshoots_the_local_fluid_velocity`;
`pbf::tests::ball_drag_gains_form_drag_at_high_relative_speed_and_low_viscosity` checks that form
drag dominates the viscous estimate at water-like viscosity). The test
`pbf::tests::ball_drag_matches_two_dimensional_stokes_scaling` keeps its historical name; it checks
that the single-ball impulse matches `4 pi mu (v_fluid - v_ball) dt` within 50 % at 10 and
100 Pa·s, which with `inv_tau_stokes = 0` now exercises the step-7 ball-surface boundary (the
isolated-disc calibration of `G`, docs/VERIFICATION.md).

**Balls are processed sequentially (Gauss-Seidel), not Jacobi.** `self.v` is updated immediately
after each ball's contribution, not accumulated into a separate array and applied once after every
ball has been visited. At this project's coupling resolution `h_c` easily spans several
neighbouring balls' worth of a packed charge, so a fluid particle commonly sits within more than
one ball's entrainment radius at once; computing every ball's `v_bar` against the same stale fluid
state and only summing the results afterwards let several individually-bounded (`a_lin <= 1`)
corrections stack on the same fluid particle well beyond what any single ball's own bound was meant
to allow. Sequential application means each ball after the first already sees the fluid's
up-to-date state, so its own relaxation is bounded against reality rather than against a snapshot
every other overlapping ball is also independently correcting. The sub-grid squeeze-film links of
§6.5 enter the same step-7 solve.

Regression test for the whole coupled system: `coupling::tests::a_coupled_charge_in_a_still_drum_settles_with_slurry_on`,
`coupling::tests::coupled_charge_never_gains_more_energy_than_the_wall_supplies`.

#### 6.2.1 History of the drag closure (superseded mechanisms and measured investigations)

Everything below describes mechanisms and numbers as of the dates given. The Stokes term it is built
on was removed on 2026-10-05, the form-drag factor changed on 2026-10-19, and the configurations
quoted (50 Pa·s, `max_balls` 150/600, 10 mm balls) are those of the investigation, not today's
defaults. The `c_rot` normalizer and cap (above) and the Gauss-Seidel ordering are still in the
code; the rest is superseded.

**Stokes-era closure (until 2026-10-05).** The relaxation time was the series combination of a
Stokes-regime and a form-drag term:

```
tau_stokes = rho_ball * r_true^2 / (4 * mu)                       (Stokes-regime disc relaxation time)
tau_form   = mass_true / (rest_density * C_D * (2*r_true) * |v_rel|)   (2D-cylinder form-drag time)
1/tau_eff  = 1/tau_stokes + 1/tau_form
beta       = slurry.ball_no_slip * (1 - exp(-dt / tau_eff))
a_lin      = beta * m_ent / (balls.mass + m_ent)            <= 1 for every mu, dt
dv_b       = a_lin * (v_bar - v_b)
```

`rho_ball = balls.mass / (pi * balls.radius^2)` (the *material* density, recovered exactly
regardless of coarse-graining, since coarse-graining scales mass and area together); `mu =
slurry.viscosity_pa_s` is the *physical* slurry viscosity; `v_rel = v_bar - v_b`; `C_D =
BALL_FORM_DRAG_COEFFICIENT = 1.0` (a circular cylinder's standard high-`Re` drag coefficient; the
`tau_form` line above lacks the factor 1/2 that was added on 2026-10-19, so the closure then acted
as `C_D = 2`).
`r_true`/`mass_true = rho_ball * pi * r_true^2` are the *true* (uncoarsened) individual media
particle's own radius/mass (`Balls::true_radius`, `crates/mill-core/src/dem.rs`) -- **not** the
coarse-grained `balls.radius`/`balls.mass` used everywhere else in this section (`h_c`, `m_ent`,
the momentum bookkeeping in `dv_i`/`coupling.impulses`/`coupling.angular_impulses` below). See
"Coarse-graining independence" below for why the two drag-rate terms specifically need the true
size while everything else correctly keeps using the coarse one. `tau_stokes` alone is the correct
regime at this
project's default of the time, 50 Pa*s (ball `Re ~ 1`), but at water-like viscosity a ball's speed (1-5 m/s)
reaches `Re ~ 10^2-10^4`, where drag is form- (not viscosity-) dominated and `tau_stokes` alone
under-predicts it by orders of magnitude -- an external review's finding. Combining the two
regimes as reciprocal relaxation times means the smaller (faster-relaxing) of the two dominates,
matching how the two drag laws' relative magnitudes actually compare, while keeping the same
`beta`-based closure and hence the same `a_lin <= 1` / `a_rot <= 1` conservation bounds as before.
`tau_form` depends on the ball's own relative speed, so `beta` (and hence `tau_eff`) is now
computed per ball rather than hoisted out as a single scalar for the whole population. The reaction
is distributed back across the same weighted fluid neighbours
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

Capping `m_ent` at the ball's own physical 2D added mass was tried during the 2026-09-22 low-fill
investigation (ss9) and reverted: `pbf::tests::ball_drag_matches_two_dimensional_stokes_scaling`
(and its sibling) specifically calibrate this closure's *large*-`m_ent` limit against real 2D
Stokes drag — `m_ent >> balls.mass` (an effectively-infinite fluid reservoir) is the physically
correct regime there, not a discretisation artefact to suppress.

**Coarse-graining independence (2026-09-26).** Follow-up to a "balls floating/launching at high
`max_balls`" investigation (see ss9): the tunnelling-risk hypothesis that investigation's previous
round left open was directly measured and retracted (zero ball-ball position-order-flip or deep-
overlap events even at `max_substep_displacement_over_diameter` ~= 1, well above the level the
user's own screenshots showed) -- but a genuine, separate bug was found instead, in this section.
Before this fix, `tau_stokes`/`tau_form` above used the *coarse-grained* `balls.radius`/`balls.mass`
(the same ones DEM contact mechanics correctly uses), which meant `simulation.max_balls` -- a pure
performance/resolution knob that is only supposed to trade off simulated ball count against
accuracy of the DEM contact network, never the *physics* itself -- silently changed the physical
strength of ball<->fluid drag. Measured directly (a controlled single-ball-in-a-uniform-fluid-patch
harness, this project's default media/viscosity): the one-sub-step drag relaxation fraction
(`|dv_b| / |v_rel|`) was **0.82 at the true 10mm particle diameter, but only 0.63 at `max_balls=600`
(20.25mm coarse diameter) and 0.26 at the 150-ball preset's default (40.5mm coarse diameter)** --
the *default*, most-used configuration was under-predicting real drag by roughly 3x, not the
high-`max_balls` configuration over-predicting it, the opposite of what the "floating at high
`max_balls`" symptom's surface appearance suggested.

The fix (`Balls::true_radius`, `crates/mill-core/src/dem.rs`; `rho_ball`/`tau_stokes`/`tau_form`
above) uses the true particle size/mass for the drag-*rate* terms only. Everything else in this
step -- the entrainment kernel's own geometric footprint `h_c = balls.radius + h`, and every
momentum-bookkeeping quantity (`dv_i`'s `balls.mass`, `coupling.impulses`'s `balls.mass * dv_b`,
`coupling.angular_impulses`'s `i_b`) -- deliberately keeps using the coarse-grained values. A first
attempt also swapped `a_lin`/`a_rot`'s own ball-side mass/inertia term (`balls.mass + m_ent` ->
`mass_true + m_ent`) to the true values and was wrong: `m_ent` is sampled over the *coarse* ball's
actual (larger) footprint, so pairing it with the tiny, now-fixed true mass over-corrected `a_lin`
toward saturation *harder* at larger `max_balls` (measured 0.82 -> 0.88 -> 0.95 -> 0.98 at
increasing coarse radius) -- backwards from the intent. Reverting just that one term (keeping
`a_lin = beta * m_ent / (balls.mass + m_ent)` on the coarse mass, since a locally-uniform fluid
field means this specific ratio is already scale-invariant under coarse-graining -- both `m_ent`
and `balls.mass` grow by the same `scale_factor^2` together) gives the intended result: **0.8177,
0.8177, 0.8178 at `scale_factor` 1x/2x/4x respectively**, matching the true-particle baseline
(0.8179) to within measurement noise. Regression test:
`pbf::tests::ball_drag_relaxation_fraction_is_independent_of_coarse_graining` (removed 2026-10-05 together with the Stokes closure; no test of that name exists in the crate today; see docs/VERIFICATION.md "Ball drag measurement").

One existing regression tightened its bound against the *old*, physically-weaker drag and needed
updating, not the fix reverting: `metrics::tests::compression_error_stays_bounded_under_violent_
lifter_cataracting` (`max_balls=150`, lifters, violent cataracting) measured worst compression
error rising from ~0.091 to ~0.183 with the corrected (genuinely stronger) coupling -- still
comfortably inside `docs/METRICS.md`'s own documented **10-20% healthy range for active splashing/
impact events**, so the bound moved `0.12 -> 0.20` (the top of that documented range) rather than
treating the new reading as a regression to chase.

**Angular resonance bug and fix (2026-09-26, same day).** Rebuilding the app after the fix above and
replaying the user's exact default config (Realtime preset of that date, `max_balls=150`, `resolution=15`)
visually still showed a dramatic scattering/launching charge, and `dissipated_power_w` reading ~3x
`power_draw_w` in steady state -- contradicting this file's own documented expectation (ss4.2) that
slurry drag should make `dissipated_power_w < power_draw_w`, the same signature as the two earlier
"gas-inflation" bugs (ss9). An earlier hypothesis that this was merely a `dissipated_energy_j`
accounting-scope gap (ball-contact dissipation legitimately exceeding wall-only power draw once the
fluid routes wall-sourced energy through a much-stronger drag) was retracted after a decisive,
per-sub-step full-system (ball + fluid) mechanical-energy-vs-cumulative-wall-work invariant check
(mirroring ss4.3's dry-only version) found a genuine violation: **worst single-sub-step margin
-8262.8 J**, against a pre-fix (`HEAD~1`) baseline noise floor of only -34.8 J over the same harness
-- a clean A/B, not measurement noise.

Root cause: `beta` saturating toward 1 (the fix above) meant step 6.2's angular-relaxation reaction
(`dw_b`, distributed back to the fluid via `c_rot = -i_b * dw_b / s`) now attempts much larger
corrections per sub-step. `s` (the exact-conservation normalizer)'s existing degeneracy gate
(`|s| > 0.1 * s_scale`, ss6.2 above) still let a small-but-passing `|s|` through, and `c_rot` -- a
plain division, with no bound of its own -- reached into the thousands, injecting a 10+ m/s
single-sub-step velocity kick into the handful of contributing fluid particles. Those particles'
now-huge velocity was then read back, next sub-step, as a neighbouring ball's own local
`omega_bar`, closing a resonant feedback loop that drove individual balls' spin to 1000+ rad/s in
testing -- unphysical rotation that then dissipated via ordinary DEM friction/restitution, exactly
matching the observed `dissipated_power_w` excess.

Fix: cap `c_rot` (`crates/mill-core/src/pbf.rs`, step 6.2) so the resulting per-particle velocity
kick can never exceed `0.2 * FLUID_SPEED_SAFETY_FACTOR * (wall speed + free-fall-across-diameter
speed)` -- the fluid speed clamp's own physically-anchored scale, at a fraction of it (a
per-*reaction* bound, well below the whole-fluid ceiling it shares a formula with). This is the
same accepted trade-off as the fluid speed clamp (ss5.2) and the ball-side impulse clamp (ss6.4):
the exact-conservation identity `c_rot` computes is deliberately broken, on the rare sub-step where
it would otherwise blow up, rather than ever injecting an unphysical kick. Measured effect at the
default config: worst single-sub-step margin -8262.8 J -> **-14.7 J**, back to the pre-fix noise
floor; `dissipated_power_w`/`power_draw_w` back to the expected `<` relationship. Regression test:
`coupling::tests::coupled_charge_never_gains_more_energy_than_the_wall_supplies`.

A second regression appeared to surface at a different config (`media.fill_fraction=0.10`, weaker
coarse-graining, `coupling::tests::low_fill_cataracting_charge_does_not_gain_energy_from_the_fluid`),
where step 6.2's correction is usually *skipped* rather than merely bounded (`s`'s degeneracy gate
failing outright, `applied = false`): that test's *ball-only* mechanical-energy-slope metric read
+26.6 W/m (failing its `<5.0` bound). **Root-caused (2026-09-27) as a false alarm, not a physics
bug.** A ball-only slope cannot distinguish "the coupling created energy" from "the fluid
legitimately handed the ball some of the energy the wall gave the fluid" -- exactly the ball<->fluid
drag this project intends. Checked with the same full-system (ball+fluid) per-sub-step invariant
this section's fix already validates against, at this exact config, across 2 seeds, over a 20s
window: the invariant held at every sub-step both times (worst margins +1.78 J and +7.80 J), while
the ball-only slope swung sign run-to-run and window-to-window (+13.8/+1.0/+4.3 W/m at 5/10/20s one
run, -1.2/+3.9 W/m at 5/20s the other) -- noise, not a steady drift. The regression test now uses
the full-system invariant directly; no code in `pbf.rs`/`dem.rs` changed as a result of this
finding. An angular ball-speed backstop added while chasing this false alarm was, together with the
pre-existing linear one, ablation-tested (2026-09-27, ss9) and removed once shown not load-bearing.

**Update 2026-10-05.** The `tau_stokes` term was removed: free-space 2D Stokes drag has no basis and the
closure was saturated, so the ball merely followed its entrained fluid. Viscous drag now comes from a
no-slip ball-surface boundary in the implicit viscosity solve (pbf.rs step 7, balls solved together with
the fluid); step 6.5 keeps only the form-drag relaxation. See docs/VERIFICATION.md "PBF ball no-slip
boundary".

### 6.3 Buoyancy (pbf.rs step 6.6)

Balls are typically sub-resolution relative to the fluid spacing and do not contribute to the
density-constraint sum, so buoyancy is modelled directly rather than emerging from the PBF pressure
field. Each ball samples the local fluid density by reusing the same Poly6 kernel and the fluid's
own step-2 neighbour grid at the ball's position:

```
rho_local = sum_j m_j * W_poly6(|x_ball - x_j|^2, h)   (over nearby fluid particles j)
rho_eff = min(rho_local, rest_density)
omega_local = (sum_j w_j * cross2(x_j, v_j)/|x_j|^2) / sum_j w_j, clamped to [-|omega|, |omega|]
g_eff = (0, GRAVITY) + omega_local^2 * x_ball  (apparent gravity in the *fluid's own* rotating frame)
impulse_on_ball = -rho_eff * (pi * r^2) * g_eff * dt   (acts through the centroid: no torque)
```

(`w_j = m_j * W_poly6(...)`, the same weights used for `rho_local`; `omega` is `drum.omega`.)
Clamping to `rest_density` avoids over-buoyancy from a locally compacted pocket and tapers smoothly
to zero as a ball nears the free surface (lower sampled density there) instead of an on/off cutoff.
`g_eff` is the general Archimedes formula `F = -rho * V * g_eff` applied with the *apparent* gravity
a slurry parcel actually feels in its own (possibly accelerating) rest frame: a parcel in solid-body
rotation has real acceleration `-omega_local^2 * x_ball` (`x_ball` measured from the drum's centre,
the world origin -- see `Drum::wall_velocity`), so `g_eff = g_vec - a_fluid = g_vec + omega_local^2
* x_ball`. `omega_local` is the *measured* kernel-weighted mean angular rate of the sampled fluid
neighbourhood, not `drum.omega` itself: an earlier version used `drum.omega` unconditionally, which
was correct for slurry genuinely centrifuged against the wall but gave every ball in a bottom pool
held mostly by gravity (not in solid-body rotation with the drum -- an external review's finding) a
large fictitious inward pull, since `omega_local` there is far below `drum.omega`. The `|drum.omega|`
clamp is a physical ceiling (a wall-driven pool cannot out-rotate the wall in steady state) that also
makes this term strictly non-increasing relative to the pre-fix `drum.omega^2 * x_ball` version, so
the fix can only remove a fictitious force, never add a new one. At `omega_local = 0` (a still pool,
or `drum.omega = 0`) this reduces to the plain world-vertical expression; a genuinely centrifuged
pool still gets the full inward buoyancy this mechanism was added for -- see this section's earlier
history: an even-earlier version omitted the centripetal term entirely, buoying every ball only
against gravity regardless of rotation speed. The reaction is applied immediately as a velocity
change split across the contributing fluid
particles, weighted by each one's share `w_j / rho_local` of the sampled local density — this must
run after step 6's velocity reconstruction, not before, or the reconstruction would overwrite it.
This mechanism did not exist in the original coupling design (`docs/PLAN.md` ss3.4 describes
buoyancy as emergent from the density-constraint push alone, with an Akinci-style boundary-particle
extension noted as future work); it is now a distinct, directly-modelled step.

Because the density is sampled at the ball *centre* with the kernel radius `h` and step 3.5 keeps
fluid particles at least `r_b + 0.05 h` away from that centre, a ball with `r_b` larger than about
`0.95 h` samples no fluid (`rho_local = 0`, the ball is skipped, `pbf.rs:2066`): for such balls
buoyancy arises only from the overlap-projection pressure of step 3.5. At the shipped defaults
(`r_b = 1 mm`, `h = 2 R / resolution = 1.6 mm` at `resolution = 40`) balls are within the regime
this step serves.

### 6.4 Stability clamp (pbf.rs step 8)

The accumulated per-ball impulse this sub-step is clamped to
`impulse_clamp(mass, dt) = F_CLAMP_G_MULTIPLE * mass * 9.81 * dt`, `F_CLAMP_G_MULTIPLE = 60.0`
(angular impulse clamped to that magnitude times `balls.radius`). `pbf.rs`'s doc comment on
`F_CLAMP_G_MULTIPLE` records that this constant used to be `3.0`, tuned down purely to survive an
earlier 3D-sphere-vs-2D-disc mass mismatch between balls and fluid particles, then `20.0` once both
were unit-depth discs with dimensionally consistent masses. It is `60.0` (`pbf.rs:957`) because a
*velocity*-matching ball<->fluid relaxation (no-slip drag) has a natural impulse scale that does not
shrink with `dt` the way a position-capped correction (the overlap push, ss6.1) does, so a clamp
tuned only against the latter would clip physically legitimate no-slip corrections (achieving
no-slip within one sub-step is what "highly viscous" means). The comment on the constant
(`pbf.rs:939`-`pbf.rs:956`) and the numbers below were written for the Stokes-era closure (§6.2.1),
whose relaxation factor saturated above roughly 10 Pa·s; today the viscous part of that exchange
arrives through the step-7 joint solve (its ball impulse is added to the same clamped total,
`pbf.rs:2290`-`pbf.rs:2298`) and step 6.5 contributes only form drag, but the clamp value was
kept. Lowering the clamp to match ss6.1's now-much-smaller contribution was tried and confirmed
wrong (it broke `pbf::tests::ball_drag_never_overshoots_the_local_fluid_velocity`, a scenario with
no overlap-push contribution at all). History (measured with the Stokes-era closure at the then
sub-step rate): residual clamp rate ~5.5-6.5% of ball-substeps at 50 and 200 Pa*s, guarded by
`coupling::tests::cascading_charge_keeps_coupling_clamp_hits_rare_once_settled_at_50_pa_s` and
`..._at_200_pa_s` (bound `<10%`, tightened from an earlier `<15%` once ss6.1's contribution became
negligible). Persistent clamping still indicates a real problem, just a rarer one than before either
fix. `clamp_hits` (count of
balls clamped this sub-step) and `fluid_momentum_change` (sum of fluid momentum change from every
mechanism above, for a Newton's-third-law cross-check against `sum(impulses)` when no clamp fired)
are both tracked in `CouplingImpulses` purely as diagnostics.

This impulse clamp is a distinct backstop from the fluid speed clamp (ss5.2 step 11). A separate
ball speed clamp (linear, `dem::BALL_SPEED_SAFETY_FACTOR`) plus an angular counterpart were added
during the 2026-09-22..26 investigation below to backstop a still-unidentified energy-injection
source; both were removed 2026-09-27 once an ablation (disable, re-run the full suite and the
full-system energy invariant at three configs) showed neither was load-bearing once the real
sources (ss6.2's `c_rot` cap, and the broad-phase margin in ss4.2 step 2) were fixed. See ss9's
2026-09-27 follow-up.

### 6.5 Sub-grid squeeze-film links (`squeeze_links`, `pbf.rs:544`; model `mac/lubrication.rs`)

When two balls (or a ball and the circular drum wall) approach to a gap that the fluid lattice does
not resolve, the particle fluid underestimates the lubrication resistance of the squeezed film and
cannot be refined enough to carry it. The shipped solver therefore adds a sub-grid film force
between such pairs, taken from the experimental grid track's calibrated model
(`crate::mac::lubrication`, imported at `pbf.rs:552` — the only place the shipped simulator uses
`mac/`). It enters the joint implicit solve of §6.2, so it is implicit and conserves momentum
between the two balls.

**Film coefficient** (2D Reynolds squeeze film between two curved surfaces, with correction terms
fitted to the exact bipolar-coordinate solutions of a disc near the drum wall and of two equal
discs in unbounded fluid; `lubrication.rs:1`-`lubrication.rs:14`). For surfaces of signed curvature
`kappa_1`, `kappa_2` (`+1/a` for a disc of radius `a`, `-1/b` for the concave drum wall of radius
`b`: `disc_curvature`, `wall_curvature`, `lubrication.rs:23`, `lubrication.rs:28`), reduced
radius `1/R_c = kappa_1 + kappa_2` and gap `h`:

```
c(h) = 3 sqrt(2) pi mu R_c^(3/2) h^(-3/2) * (1 + c1 eps + c2 eps^2),     eps = h / R_c
Q    = R_c^3 (kappa_1^3 + kappa_2^3)
c1   = 0.788 + 0.265 Q,            c2 = max(0, 0.006 + 0.057 (Q - 0.25))
```

(`SqueezeFilm::new`, `lubrication.rs:53`; `coefficient`, `lubrication.rs:67`). For two equal discs
`R_c = a/2` and `Q = 1/4`; for a disc against a flat wall `R_c = a` and `Q = 1`. `c` has units of
Pa·s (force per unit depth per unit approach speed) when `mu` is the dynamic viscosity, which is
what `squeeze_links` passes (`pbf.rs:559`). The leading `h^(-3/2)` is the 2D (cylinder) analogue of
the sphere's `h^(-1)` divergence; the corrections reproduce the exact solutions to about 1 % for
`eps <= 0.4` (`lubrication.rs:12`-`lubrication.rs:13`).

**Blend with the resolved fluid.** The film is used only for the part the lattice does not already
deliver. `resolved_ratio(h/dx_f)` (`lubrication.rs:95`) is the fraction of the exact coefficient
that the resolved flow provides at a gap of `h/dx_f` fluid cells, tabulated from grid measurements:
1 for gaps of 8 cells or more, 0.99 at 4, 0.90 at 2, 0.25 at 1, 0.06 at 0.5 and 0 for touching
surfaces (linear interpolation). With `dx_f = h_kernel/2` the fluid spacing and gap
`g = |x_a - x_b| - 2 r_b` (ball–wall: `g = R - r_b - |x_a|`), `squeeze_links` builds a link when
`g < reach = H_START * dx_f` with `H_START = 8` (`lubrication.rs:19`), with

```
s       = 1 - resolved_ratio(g / dx_f)                          (weight of the sub-grid film)
h_f     = max(g, MIN_FILM_GAP_RATIO * r_b),   MIN_FILM_GAP_RATIO = 0.02     (pbf.rs:533)
c_link  = min( SQUEEZE_GAIN * s * c(h_f),  c_max ),   SQUEEZE_GAIN = 1.0    (pbf.rs:524)
c_max   = SQUEEZE_STIFFNESS_CAP * m_b / dt,           SQUEEZE_STIFFNESS_CAP = 8   (pbf.rs:529)
```

`MIN_FILM_GAP_RATIO` evaluates the film at no less than 2 % of the ball radius (surface roughness:
touching or slightly overlapping balls would otherwise get an unbounded coefficient). The stiffness
cap bounds `dt c_link / m_b <= 8`: beyond a few, the film already locks the relative normal motion
within the sub-step, and a larger value would only make the joint CG system ill-conditioned (a stiff
chain through a dense, settled charge).

**Closing-only filter.** A link is kept only if the partners are approaching,
`n . (V_a - V_b) > 0` with `n` the unit vector from ball `a` toward its partner (outward for the
wall), evaluated with the ball velocities at the start of the solve (`ball_v0`, which already
include this sub-step's earlier impulses; `pbf.rs:2237`-`pbf.rs:2241`). A film resists squeezing
only (it cavitates when the gap opens), and a symmetric film kept a settled charge jittering; the
filter also spares about half of the links their film model. The active set is fixed for the whole
solve, so the system stays linear.

**How it enters the joint solve.** The film force on ball `a` is
`-c_link [n . (V_a - V_b)] n`, and the opposite on `b` (the wall is immovable, `V_b = 0`). In the
implicit system of §6.2 this is the term `dt c_link n n^T (V_a - V_b)` in the ball rows (and in the
Jacobi preconditioner's diagonal) — the same operator, applied to ball unknowns only, with no
fluid unknowns attached (`pbf.rs:724`-`pbf.rs:736`). Because both partners' rows carry the matching
`+/-` terms it is symmetric positive semi-definite and conserves the pair's linear momentum; wall
links remove momentum to the wall. Limitations: only the circular drum wall is modelled (not
lifters), and only the normal film (the tangential coefficient in `lubrication.rs:79` is used by the
grid track only). Regression tests: `pbf::tests::squeeze_links_exist_only_for_gaps_below_the_resolved_range`,
`pbf::tests::squeeze_link_damps_the_approach_and_conserves_momentum`.

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

## 8. Metrics & grinding instrumentation (`metrics.rs`, `stress.rs`, `lib.rs`)

`metrics::compute(balls, fluid, drum, drum_angle)` is a pure function of the current populations and
derives:

- **Toe/shoulder angles** (`charge_toe_shoulder`): balls within `WALL_MARGIN_BALL_RADII (2.5) *
  balls.radius` of the wall, per `drum.sdf_world` (so a ball riding a lifter's tip counts as
  wall-adjacent, not just one near the plain circular wall between lifters -- an earlier version
  compared against the bare undecorated-circle distance and missed lifter-borne balls; identical
  to that version when `lifters.count == 0`), approximate the charge's outer layer; if there are at
  least
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
- **Lacey mixing index** (`mixing_index`, `metrics.rs:505`): fluid particles (the dye tracer) are
  binned into a `MIXING_GRID_SIZE (16) x 16` grid over the drum bounding box;
  `M = (S0^2 - S^2) / (S0^2 - Sr^2)` with, over occupied cells `c` (particle count `n_c`, mean dye
  `cbar_c`, overall mean dye fraction `p`):

  ```
  S^2  = sum_c n_c (cbar_c - p)^2 / sum_c n_c              (particle-count-weighted variance)
  S0^2 = p (1 - p)                                          (fully segregated)
  Sr^2 = S0^2 / n_bar,     n_bar = sum_c n_c^2 / sum_c n_c  (weighted mean occupancy)
  ```

  The variance and the occupancy are both **particle-count weighted** (the weighted-variance
  "effective sample size"), so a sparse splash cell holding one stray particle contributes as little
  as its one particle warrants; an earlier unweighted version let such cells dominate and understated
  the index. Clamped to `[0,1]`; if `S0^2 - Sr^2` is ~0 (including on average one particle per
  occupied cell) the index is degenerate and returns `1.0`; `None` for fewer than 2 occupied cells.
  Regression: `metrics::tests::mixing_index_weights_by_cell_particle_count_not_cell_count`.
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
zero/empty-defaults these fields; `Simulation::metrics()` (`lib.rs:370`) fills them in from its own
accumulated state after calling `compute`:

- `power_draw_w`, `torque_nm`, `dissipated_power_w`: EMA-smoothed (time constant
  `GRINDING_STATS_EMA_TAU_S = 1.0` s, `lib.rs:48`, i.e. `alpha = 1 - exp(-sub_dt/tau)` applied every
  sub-step in `Simulation::update_grinding_stats`, `lib.rs:516`). `power_draw_w =
  (dem_stats.wall_work_j + fluid_stats.wall_work_j) / sub_dt` (`lib.rs:529`) -- the drum's total
  work rate against the ball charge (friction, and the normal force on a lifter face) **and**
  against the slurry (the no-slip blend's wall work, ss5.2 step 7), not ball-wall friction alone;
  `torque_nm = power_draw_w
  / omega` (0 when `|omega| <= 1e-6`, since torque is undefined rather than infinite at zero
  rotation); `dissipated_power_w = max(dissipated_energy_j, 0)/sub_dt` (DEM-side only -- see ss9's
  note on why it excludes fluid-side dissipation).
- `collision_rate_per_s` and `impact_energy_histogram` (`ImpactEnergyHistogram { bin_edges_j,
  counts_per_s }`) come from the **recorded-impact statistics** of `stress.rs` (below), not from
  `DemStepStats::collision_count` / the 12-bin 2D-disc histogram of ss4.2 (those stay internal to
  the DEM step): `collision_rate_per_s = StressAccum::total_rate()` (debounced impacts per second,
  whole simulated population; `lib.rs:552`), and the histogram is the pooled 28-bin stress-energy
  histogram with `bin_edges_j = se_bin_edges() * (m_3 g d)` in joules and `counts_per_s` the
  EMA-smoothed rate per bin (`lib.rs:391`-`lib.rs:402`).
- `stress`: `stress::StressMetrics`, below.
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

### 8.1 Energy budget (`EnergyBudget`, `lib.rs:61`)

A cumulative, read-only account of the whole system's mechanical energy over a window, all in J per
metre of mill depth, available through `Simulation::energy_budget()` (`lib.rs:218`; restart the
window with `reset_energy_budget`, `lib.rs:211`; deliberately not part of `Metrics`). It is
accumulated every sub-step (`lib.rs:492`-`lib.rs:504`):

```
shaft_work              += dem_wall_work + fluid_wall_work          (the motor's work, = what power_draw integrates)
ball_contact_dissipation += DemStepStats::dissipated_energy_j
fluid_wall_slip          += FluidStepStats::wall_work_j - ke_delta_noslip_j
fluid_viscous            -= ke_delta_viscosity_j
fluid_clamp_removed      -= ke_delta_clamp_j
interface_created        += DemStepStats::coupling_work_j + ke_delta_drag_j + ke_delta_buoyancy_j
unattributed = delta_mechanical - ( shaft_work - ball_contact_dissipation - fluid_wall_slip
                                     - fluid_viscous - fluid_clamp_removed + interface_created )
```

(`delta_mechanical` is the change of `dem::mechanical_energy_f64 + FluidParticles::mechanical_energy_j`
since the window start.) `interface_created` is the net energy *created* at the ball-fluid
interface (the kinetic-energy change the balls receive from the applied fluid impulses plus the
fluid-side change of the drag and buoyancy exchanges), physically `<= 0` since drag dissipates.
`unattributed` is the closure error: for PBF the position-projection numerical dissipation and the
unmetered cohesion/overlap-push steps land there; it is zero in a perfectly closed budget.
`fluid_wall_angular_impulse` accumulates the fluid's angular impulse on the drum (mean torque
`-fluid_wall_angular_impulse / elapsed_s`).

### 8.2 3D-converted stress statistics (`stress.rs`, `lib.rs`)

Grinding itself is not simulated. What is reported are the quantities mechanistic grinding models
are driven by: the distribution of **stress energies** (`SE`, the relative kinetic energy of an
impact) and the frequency of stressing events, plus the power drawn (module comment,
`stress.rs:1`-`stress.rs:21`). The 2D disc mass overstates `SE` by a ball-size-dependent factor
(about 750 for a 2 mm ball), so `SE` uses the mass of the **3D sphere** a disc stands for
(`Simulation::stress_scales`, `lib.rs:233`):

```
m_3  = rho * (4/3) pi r^3 = m_disc * (4/3) r         (unit-depth disc mass rho pi r^2 times 4r/3)
SE   = (1/2) m_red (v_n^2 + v_s^2),     m_red = m_3/2 (ball-ball),  m_3 (ball-wall)      (classify, stress.rs:68)
E^   = SE / (m_3 g d)                                  (dimensionless; d = 2r)
```

`v_n` is the pre-solve normal approach speed and `v_s` the tangential slip speed at the contact
point including spin, both from the `ImpactRecord` (ss4.2). Recording only happens for approach
speeds above the size-relative threshold `restitution_velocity_threshold(r)` (ss4.2). An impact is
**drop** type if `v_s <= v_n` and **shear** type otherwise (`vs > vn`), giving four classes
(`CLASS_NAMES`, `stress.rs:33`: ball-ball and ball-wall, drop and shear).

- **Histogram.** `E^` is binned into `SE_BINS = 28` log-spaced bins starting at `10^-3` with 4 bins
  per decade (`stress.rs:28`-`stress.rs:30`; edges `se_bin_edges`, `stress.rs:51`), values outside
  the range going to the end bins.
- **Debounce** (`DEBOUNCE_FACTOR = 0.25`, `stress.rs:48`). Re-contacts of the same pair (or of a
  ball with the wall) closer in time than `0.25 sqrt(d/g)` count as one impact episode; without it
  the count keeps rising as the sub-step shrinks (finer steps resolve more micro-bounces of a
  resting contact). The window is about 10x the sub-step of the validated accuracy settings.
- **Thresholds.** The rate and cumulative energy carried by impacts above
  `E^ in {0.03, 0.1, 0.3, 1, 3, 10}` (`THRESHOLDS_REL`, `stress.rs:42`) are tracked as the
  "effective" stress energy; the time two conditions need to reach the same cumulative effective
  energy is a relative grinding-time proxy.
- **Reported** (`StressMetrics`, `stress.rs:243`): impacts per ball per second per class and in
  total; the mean `E^` per class, for drop and shear pooled, and over all impacts (the
  stress-intensity measure `stress_intensity_scaled`); wall power per ball; specific power; the
  effective power per ball per threshold; and cumulative wall / specific / effective energies.
- **Per-ball power and the depth basis.** "Per ball" power uses a one-ball-thick slab (depth equal
  to the ball diameter): `P_ball = P_2D d / N` (`StressAccum::metrics`, `stress.rs:186`), specific
  power `P_ball / m_3`, cumulative specific energy likewise. This is equivalent to an implied 3D
  solid fraction: a sphere of volume `(4/3) pi r^3` is represented by a cylinder of volume
  `pi r^2 d = 2 pi r^3`, so `phi_3D = (2/3) phi_2D` (0.547 at the default `phi_2D = 0.82`); ratios
  between conditions do not depend on this depth basis.

The results are meant for **relative** comparison between conditions run with the same code and the
same numerical settings; with coarse-graining (`k > 1`) the simulated ball is not the real ball and
these statistics must not be used (`stress.rs:19`-`stress.rs:21`; ss9, ss3.2).

---

## 9. Known limitations

> **Status as of 2026-10-20.** The simulator runs the PBF solver described above (the DFSPH and
> grid-based rebuild tracks of 2026-10-01..05 are reference instruments, not the shipped solver;
> docs/VERIFICATION.md, "Simulation back on PBF (2026-10-05)"). Bulk quantities (power, torque) are
> *not* invariant under coarse-graining or fluid resolution: coarse-graining changes the power draw
> by +17...39 % at large `D/d` and does not preserve per-impact statistics, and the fluid resolution
> sets the converged power to a few percent only at the finer tiers (docs/VERIFICATION.md, "What
> coarse-graining preserves"; the quality presets' notes in `web/src/params/presets.ts`). Against
> the strict verification gates the shipped PBF solver passes 6 of 22 and fails 16 (hydrostatic
> compression, wall pressure, spin-up, Taylor–Couette ball and wall torque and others), re-measured
> on 2026-10-19: see docs/VERIFICATION.md, "2026-10-19: strict gates re-measured, ball drag
> measured", for the table and for the unmet gates' causes. Items below that are marked *History*
> describe a defect or mechanism that has since been fixed or replaced; they are kept for the
> measurements and the reasoning, with the defaults of their date.

- **Rolling resistance (ss4.2 step 7) damps world `omega` toward zero, not the relative spin at
  each contact.** It applies no reaction torque to a ball-ball contact's partner (so it does not
  conserve angular momentum) and, for a ball-wall contact, targets zero rather than `drum.omega`
  (so it fights, rather than assists, a ball correctly rolling without slipping on a rotating
  drum). A per-contact, momentum-conserving, wall-frame-aware rewrite was tried (an external
  review finding, symmetric ball-ball correction plus a `drum.omega`-referenced wall correction)
  and reverted: it destabilised this crate's core cataracting/cascading-density regression tests
  (`dem::tests::a_cascading_charge_stays_dense`,
  `metrics::tests::compression_error_stays_bounded_under_violent_lifter_cataracting`), which guard
  against the much more severe fluidised-charge energy-injection bug this solver was extensively
  tuned against (`MAX_RECOVERY_FRACTION`, step 6's bidirectional restitution). At this project's
  default `rolling_friction = 0.01` the simplification's practical effect is small (friction, step
  5, is what actually drives a resting/rolling ball toward the wall's tangential speed; rolling
  resistance only trims the residual spin on top), but revisiting it needs dedicated re-tuning
  against those regressions, not a drive-by correctness fix.
- **Continuous collision detection covers the wall/lifters, not ball-ball.**
  `Drum::toi_swept` (ss2) clamps a ball's predicted advance to its geometric time-of-impact against
  the wall/lifters before the non-penetration solve runs, gated to engage only when the predicted
  end-of-step overlap exceeds `MAX_RECOVERY_FRACTION`'s one-pass recovery budget -- so it acts as a
  backstop for genuinely severe events, not a routine per-sub-step intervention (an unconditional
  version measurably regressed `max_ball_wall_overlap_fraction` on a packed, lifters-enabled
  cataracting charge, by repeatedly truncating ordinary tangential sliding along the wall's curved
  boundary -- see `dem.rs` step 2.5's own doc comment for the full account). A ball-ball CCD pass
  was tried and reverted for a related but distinct reason: in a dense granular bed a ball routinely
  has several simultaneous near-touching neighbours (ordinary jostling, not a tunnelling risk), and
  clamping to the *minimum* time-of-impact across all of them froze far more motion than the
  discrete solver's own bounded-recovery pass already handles safely (the mechanism validated by
  the fluidised-charge energy-injection fix, `MAX_RECOVERY_FRACTION`) -- balls piling up unable to
  settle pushed *more* of them into the wall, not fewer. `DemStepStats::
  max_substep_displacement_over_diameter = |v| * dt / (2*radius)` remains the live diagnostic for
  ball-ball tunnelling risk specifically; at this project's current defaults it stays comfortably
  under the XPBD `< 1` criterion (`docs/METRICS.md`), and a genuine, isolated ball-ball tunnelling
  event has not been observed in this crate's regression tests (`dem::tests::
  two_balls_head_on_collision_conserves_momentum`'s anti-tunnelling assertion).

  **History -- 2026-09-26 follow-up: this margin was widened (current formula: ss4.2 step 2), and
  the `< 1` criterion is not actually comfortable at every configuration.** Investigating a user report of balls suddenly launching
  "with no visible cause" and moving in coherent clumps at a high `simulation.max_balls` found two
  things. First, a real gap in step 2's broad-phase margin: `ball_ball_pairs` is built once (from
  predict-step positions) and reused unchanged for every one of step 3's `iterations` depenetration
  passes plus steps 5-7 -- but the margin only budgeted for the predict step's own displacement
  (`max_speed * dt`), not for the *additional* motion a ball can pick up during the solve itself (up
  to `MAX_RECOVERY_FRACTION * 2r` per pass, from a contact resolving). A third ball just outside the
  predict-only margin could be pushed into a genuine new overlap by that solve-phase motion and
  never be recorded as a contact at all this sub-step (absent from `ball_ball_pairs`/`ContactBook`
  entirely, not merely under-resolved). Fixed by adding `iterations.max(1) * MAX_RECOVERY_FRACTION *
  2r` to the margin (`dem.rs` step 2) -- a widening of the existing discrete search radius, not a
  new swept/CCD subsystem, so it does not reintroduce the ball-ball CCD regression described above.
  Measured across 4 seeds at `max_balls = 1500`: this consistently (4/4 seeds) reduced the rate of
  large (`>0.5 m/s`) single-sub-step speed jumps by roughly 5-14%, with no regression across the
  full 124-test suite.

  Second, and more fundamentally: **`max_substep_displacement_over_diameter` does *not* stay under
  the `< 1` criterion at every configuration** -- at `max_balls = 1500` (small coarse-grained balls,
  `d_eff ~= 12.8 mm` vs. the Realtime default's `~40.5 mm`), the same margin fix measured this ratio
  at 0.64-1.07 across seeds (both *before and after* the margin widening above -- the widening does
  not, and structurally cannot, change this ratio, which is a pure kinematic quantity computed
  before the broad-phase even runs), i.e. genuinely crossing the documented safety threshold in
  roughly half of tested seeds. This is not a bug the margin fix (or any broad-phase change) can
  close: it is a direct consequence of `max_substep_displacement_over_diameter = |v| * dt /
  (2*radius)` scaling inversely with ball radius at a fixed `dt` and fixed typical peak collision
  speed (~2-3 m/s here, itself unremarkable) -- shrinking `d_eff` by raising `max_balls` shrinks the
  denominator directly. The occasional resulting "sudden unexplained pop"/clump-launch visual
  symptom is most plausibly a real (if rare per-run) discrete-collision-detection near-miss at this
  ball size, not a coding defect elsewhere.

  **History -- 2026-09-27 follow-up: implemented and measured** (current definition: ss3.3).
  `Params::effective_substeps` (`params.rs:690`, right after `Params::effective_media`) mirrors `Params::effective_fluid_resolution`'s auto-raise
  pattern exactly: it computes a reference impact speed `v_ref = sqrt(2 * g * D)` (free fall across
  the full drum diameter, a deliberately generous upper bound), requires
  `v_ref * dt <= TARGET_RATIO * d_eff` with `TARGET_RATIO = 0.5` (half the `< 1` criterion) and
  `dt = 1 / (60 * substeps)`, and returned `simulation.substeps.max(min_substeps).min(16)` at that
  date (today `substeps.max(min_substeps.min(16))`, so an explicit request above 16 is honoured) --
  raising the sub-step rate, which divides `dt` and so this ratio directly, whenever a coarse-grained
  `d_eff` would otherwise push it too close to 1. `Simulation::step` (`lib.rs`) calls this method
  instead of reading `simulation.substeps` directly when sizing `fixed_sub_dt()`. `TARGET_RATIO`
  was chosen so the then `Params::default()` (63 mm ball in a 1 m drum) and every quality preset
  of the time at the 0.30 fill fraction were unchanged (`substeps` stayed at 8; the test that
  checked this, `effective_substeps_matches_requested_at_every_quality_preset`, no longer exists
  -- it was replaced by `effective_substeps_matches_requested_for_the_legacy_1m_drum` and
  `effective_substeps_is_capped_at_every_quality_preset_at_the_shipped_default` when the default
  became the 2 mm / 63 mm bench configuration, which needs 16); at the exact
  `max_balls = 1500` / Realtime `resolution = 15` repro above, it raised `substeps` from 8 to 12.
  Measured at this exact config, a fresh 3-simulated-second window: `max_substep_displacement_over_
  diameter`'s worst value over the window went from 0.64-1.07 across seeds (before this fix) to
  **0.3697** (after) -- comfortably back under 1. Like `effective_fluid_resolution`, this is
  capped at 16 for the auto-raise (`params::tests::effective_substeps_is_capped_at_the_validated_maximum`;
  the validated range of the requested `simulation.substeps` itself is 1-64) so a pathological
  combination (`max_balls` far above `N_true` with a very small `media.ball_diameter_m`) cannot
  make a frame's DEM/PBF cost unboundedly large; in that residual regime the ratio may still
  approach or exceed 1, a real disclosed limit of the auto-raise, not a bug. Surfaced (never
  applied silently) in the parameters panel's derived-values block as "Effective sub-steps"
  (`web/src/ui/paramsPanel.ts`, mirrored for instant UI feedback by
  `web/src/params/derived.ts`'s `effectiveSubsteps`) alongside "Effective slurry resolution".

  **History -- 2026-09-27 follow-up: the widened margin is load-bearing, confirmed by ablation.** While
  removing other stability backstops added during this investigation that turned out not to earn
  their place (see the "Unphysical scatter" entry's final follow-up below), disabling this margin
  alone (`solve_margin = 0`) was tested the same way: the full suite still passes, but
  `coupling::tests::coupled_charge_never_gains_more_energy_than_the_wall_supplies` (the full-system
  energy invariant, at the browser's default config) fails with a genuine **-22.3 J** single-sub-
  step violation -- a real missed-contact defect, not a theoretical nicety. This margin stays.
- **`max_ball_wall_overlap_fraction` reads elevated under `lifters.count > 0`** (measured in 2026-09
  at the then-default 63 mm media: ~0.87-0.90 worst-case over a 20 s cataracting run at
  `lifters.count = 8`, versus a much smaller figure with no lifters) -- a real, XPBD contact-convergence residual (cataracting off a lifter reaches higher
  peak ball speeds than no-lifter cascading), separate from the CCD backstop above, which targets
  only the pathological >=100% "ball centre fully passed through solid" case and leaves this
  residual untouched by design. *History:* an earlier hypothesis attributed part of this reading to
  the lifter SDF's near-corner approximation (the min over each edge's infinite supporting line,
  exact for interior points but under-estimating the distance to an exterior convex vertex), which
  biased contact to register *earlier*, not later -- the opposite of a tunnelling risk. That
  approximation no longer exists: the SDF is an exact convex-polygon signed distance (ss2,
  `geometry.rs:247`). Reducing this residual further is a `dem_iterations`/performance-budget
  question, not a geometry bug.
- **Rheology: only a regularised Bingham option, qualitatively validated.** The PBF solver
  implements Newtonian and Papanastasiou-regularised Bingham viscosity (ss5.3; `Rheology::Bingham`
  with `yield_stress_pa > 0`), with the unyielded viscosity capped at `101 mu`
  (`BINGHAM_PLASTIC_RATIO = 100`); there is no Herschel–Bulkley power law (that exists only in the
  experimental grid track), no thixotropy and no slurry-borne fine solids. Against the grid
  reference the Bingham response is only qualitatively right (docs/VERIFICATION.md, "PBF: Bingham
  and coarse-graining k").
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
- **History: the broad-phase displacement margin (ss4.2) covers the largest single-ball speed, not
  the largest pair *closing* speed** (concern raised when the margin was `max(2r*1.05, 2r +
  max_speed*dt)`; the margin today also carries the `n_it * 0.2 * 2r` solve term, ss4.2 step 2).
  The margin uses `max_speed`,
  the fastest any one ball moves this sub-step; a head-on pair (each moving toward the other at
  close to `max_speed`) closes at up to twice that rate, so the margin can in principle be half of
  what step 2's own rationale intends for that specific pair.

  **History -- 2026-09-27 follow-up: tested, not adopted.** Widening the margin's `max_speed * dt` term to
  `2.0 * max_speed * dt` (to cover a head-on pair's closing speed rather than one ball's own speed)
  was measured via `coupling::tests::settle_and_measure_full_system_energy_invariant` at
  `max_balls = 1500` / `resolution = 15` (1.0 s settle + 2.0 s measure, single seed):
  `worst_margin_j` went from **-40.49** (baseline) to **-39.17** (2x margin) -- no meaningful
  improvement, well within this chaotic many-body system's run-to-run noise. For calibration,
  raising `substeps` alone from 8 to 16 on the *same* baseline config (unrelated to this margin)
  shifted the measured margin from -40.49 J to -90.81 J, showing this single-seed energy metric is
  far too noisy across any `dt`/margin change to detect a small real effect this way. **Decision:
  not adopted, reverted, `dem.rs` unchanged.** The theoretical concern turned out not to have a
  clean single-sub-step failure mode in practice: the broad-phase grid is built from this
  sub-step's *post-predict* positions (already reflecting the full actual relative motion for this
  `dt`, not a pre-motion estimate), so a head-on pair's closing speed is already fully reflected in
  their current separation by the time the margin is evaluated -- a genuine full tunnel-through
  within one sub-step is the already-documented, separate "no ball-ball continuous collision
  detection" limitation above, not this margin.
- **`dissipated_power_w` excludes fluid-side viscous dissipation.** It is a purely DEM-side
  accounting (ss4.2's `dissipated_energy_j`, net of the wall's own work input); the slurry's
  viscosity (ss5.3) also removes real mechanical energy from the balls (ss6.2's drag) that never
  appears in this metric. A coupled steady-state run legitimately shows
  `dissipated_power_w < power_draw_w`; this is not a missing-energy bug.
- **Cohesion/adhesion accelerations are calibrated against `|g|`, not a first-principles surface-
  tension unit conversion.** `COHESION_ACCEL_FACTOR`/`ADHESION_ACCEL_FACTOR` (pbf.rs ss6.1a)
  express each mechanism's peak acceleration as a multiple of gravitational acceleration, so a
  Bond-number similarity argument (`Bo = rho*g*L^2/gamma`) would not survive an arbitrary change of
  `g`. `GRAVITY` is a compile-time constant in this crate (ss1), so no inconsistency can arise in
  practice; this is a calibration choice, not a bug, but it means these two knobs are not
  independently portable to a different gravity without re-deriving their scale.
- **The coupling impulse clamp (ss6.4, `pbf.rs` step 8) is not momentum-conserving.** A ball's
  impulse is clamped for stability *after* the fluid has already received the full (unclamped)
  reaction, so each clamp hit loses momentum from the system as a whole -- the same trade-off the
  fluid speed clamp (ss5.2 item 11) already makes deliberately. `coupling_clamp_hits` is surfaced as a
  metric precisely so a run where this matters (persistently nonzero) is visible; a healthy run
  keeps it rare (`coupling::tests::cascading_charge_keeps_coupling_clamp_hits_rare_once_settled_at_*`).
- **`coupling::tests::settle_and_measure_full_system_energy_invariant`'s tolerance floor does not
  scale with ball count (found 2026-09-27, still open).** Its tolerance formula
  (`tolerance_j_per_step = (1e-2 * dem.balls.mass * (omega*radius_m)^2.max(1.0)).max(15.0)`,
  `coupling.rs:1083`; the floor was 5 J when this was measured and was raised to 15 J on 2026-09-27
  for the 63 mm media of the day) uses *per-ball* mass and a fixed floor, both effectively
  independent of ball *count* -- it was only ever calibrated/validated at the browser's 150-ball
  default of that date. At `max_balls = 1500` the 5 J floor was measurably too tight even at
  baseline, unrelated to the `effective_substeps`/broad-phase-margin
  work above: a fresh measurement (1.0 s settle + 2.0 s measure, single seed) found
  `worst_margin_j = -40.49` J at the current `substeps = 8`, and `-90.81` J at `substeps = 16` --
  both most likely aggregate per-particle floating-point/contact-resolution noise summed across
  many more (1500 vs. 150) independent balls, not a genuine "energy created from nothing" solver
  defect. This is a test-calibration gap (the tolerance should probably scale with ball count),
  left open for a future session -- the *displacement-ratio* criterion (`max_substep_displacement_
  over_diameter < 1`, docs/METRICS.md) is the oracle actually used to validate the `effective_
  substeps` fix above, not this energy invariant.
- **2D areal packing (`packing_fraction_2d`, default `0.82`) is not 3D voidage.** Random close
  packing of equal discs in 2D is denser (~18% void fraction) than random close packing of equal
  spheres in 3D (~36-40% void fraction), so the same `slurry.fill_fraction` corresponds to a much
  higher *interstitial* filling `U = slurry_area / (media.fill_fraction * (1 - packing_fraction_2d)
  * drum_area)` than the equivalent 3D mill, and the simulated slurry pool reads correspondingly
  deeper than a real one at the same nominal fill fraction. This is an unavoidable consequence of
  the project's 2D cross-section scope (see the entry above), not a parameter-tuning bug; the params
  panel surfaces `U` directly (docs/PARAMETERS.md) so its magnitude is visible rather than hidden
  inside the fill-fraction number.
- **Ball drag in water is overestimated (known limit, measured 2026-10-19).** At the shipped resolution
  (fluid spacing 1/2 to 1/4 of the ball diameter) a ball moving through water meets 5-17 times the cylinder drag
  at Re 1e3-3e3, and about 0.9-1.9 times at Re 1e4, with a factor-2 dependence on the sub-step count there. About
  60-95 % of it is the step 3.5 momentum exchange that stands in for a pressure on the ball; the form-drag closure
  and the viscous solve are negligible. Scaling the exchange, adding the balls to the density constraint, and a
  quadratic exchange rule were tried and not adopted: each either leaves the excess or degrades compression and
  the energy balance of the whole mill (docs/VERIFICATION.md, 2026-10-19). Consequence: motion of balls through
  fast water (settling, impact into a pool, a ball thrown into the slurry) is damped too strongly; bulk
  quantities (power, toe/shoulder angles) are not affected beyond the other documented errors. A real fix needs a
  pressure-based solver, which is 2-6x slower than the PBF solver here. This is one of the causes of the
  unmet strict gates: the full 22-gate table (6 pass, 16 fail), the attribution of the excess drag to the
  step 3.5 exchange and the rejected experiments are in docs/VERIFICATION.md, "2026-10-19: strict gates
  re-measured, ball drag measured", with raw logs under `docs/measurements/2026-10-19/`.
- **Coarse-graining distorts per-impact statistics.** The collision rate, the impact-energy
  histogram and all stress statistics (ss4.2, ss8.2) are derived from *simulated* impacts, whose
  count scales as `1/k^2` and whose energy changes with the ball size relative to the true
  (uncoarsened) population (ss3.2's mass-preserving substitution). Total charge mass is preserved
  regardless of `k`, but that does **not** make bulk quantities invariant: power draw changes by
  +17...39 % at large `D/d` (docs/VERIFICATION.md, "What coarse-graining preserves"; the
  2026-09-era claim that power, torque and toe/shoulder are unaffected was wrong). Coarse-graining
  is Off by default (ss3.1). The UI marks the affected readings as unavailable ("Coarse-grained",
  `web/src/metrics/specs.ts`) whenever `coarse_graining_factor > 1` rather than showing a number
  that cannot be compared to a real mill.
- **Restitution (ss4.2 step 6) can apply a tensile (pulling) normal impulse to a ball-ball or
  ball-wall pair.** An external review read this as violating Signorini's non-tensile contact
  condition (`F_n >= 0`). Step 6 drives the *relative normal velocity* to
  `target = max(-e*v_n_pre, 0, v_n_pre)` -- in both directions, not just topped up when it falls
  short. When step 3's bounded depenetration has already pushed a pair apart faster than `target`
  justifies (`v_n_now > target`), `delta_v_n = target - v_n_now < 0` and the impulse pulls them back
  together -- but the *result* of that impulse is the pair's relative normal velocity landing
  exactly at `target >= 0`, never negative. No pair is ever left approaching after this step, so no
  interpenetration or adhesion follows from it; this is the deliberate energy cap
  `MAX_RECOVERY_FRACTION` and this step's own doc comments describe (closing the fluidised-charge
  energy-injection bug), not a Signorini violation at the velocity level a downstream step could
  observe.
- **`RESTITUTION_VELOCITY_THRESHOLD` (0.02 m/s) looks smaller than one sub-step's gravity increment
  (`g*dt ~= 0.041 m/s` at 240 Hz, the sub-step rate when this was analysed; `0.0102 m/s` at today's
  `1/960 s`), suggesting a settled bed would misread as colliding every sub-step.** It does not: `v_pre` (the approach speed step 6 tests against the threshold) is
  snapshotted *before* step 1 applies gravity, and step 6 itself drives every contact's actual
  relative normal velocity to (for a resting contact) exactly `0`, not merely toward it. A settled
  bed therefore enters each new sub-step already at `v_n_pre ~= 0`, not `-g*dt`, and is correctly
  read as resting rather than as a fresh impact.
- **`morris_weights` (ss5.3) evaluates the Morris viscosity Laplacian with `spiky_grad`, a kernel
  gradient that does not vanish at the origin.** Morris's (1997) second-order argument assumes a
  kernel whose gradient does go to zero there; using one that doesn't is a theoretical mismatch an
  external review correctly flagged. In this discretisation it stays benign: the coefficient
  `c_ij = m_j*2*mu/(rho_i*rho_j) * |x_ij . grad_W_ij| / (|x_ij|^2 + eta^2)` (applied today in the
  central-force form of ss5.3, which uses the same coefficients) has a numerator that is
  `O(|x_ij|)` (Spiky's gradient magnitude itself is `O(h - |x_ij|)`, finite at the origin, but the
  dot product `x_ij . grad_W_ij` still carries one explicit factor of `|x_ij|`), so `c_ij -> 0` as
  `|x_ij| -> 0` despite the non-vanishing gradient, and `c_ij >= 0` for every pair regardless --
  `L` stays symmetric positive semi-definite, the only property `solve_implicit_viscosity`'s
  conjugate-gradient solve relies on. Swapping in an origin-vanishing kernel gradient is possible
  but would change the discretisation's effective viscosity scale, so it is a recalibration, not a
  drop-in correctness fix.
- **History -- unphysical scatter at low fill fraction / coarse fluid resolution (2026-09-22
  investigation; configurations are those of that date, not today's defaults).**
  A user report of media balls "flying" unphysically out of a cascading charge was reproduced at a
  specific configuration: Realtime quality preset, `max_balls = 150`, `simulation.resolution = 15`,
  `media.fill_fraction = 0.10` (well below this project's default `0.30`). Root operating-envelope
  mismatch: at low `media.fill_fraction`, coarse-graining (ss3) weakens, so the effective simulated
  ball diameter (`d_eff`) shrinks toward or below the fluid lattice spacing `dx = drum_radius_m /
  resolution` — the params panel's "Fluid spacing vs ball diameter" (`dx/d_eff`) warning
  (`web/src/ui/paramsPanel.ts`) already flags exactly this ratio exceeding 1. Under-resolved
  coupling of this kind lets the ball<->fluid exchange (ss6) see badly-aliased local fluid state.
  Two changes landed initially, validated against the entire existing suite (`cargo test -p
  mill-core`, 120/120 passing, no regression): a fluid speed clamp running immediately after
  velocity reconstruction, before the ball<->fluid exchange reads it, and a ball speed safety
  backstop mirroring the fluid's own. **Both were later removed (2026-09-27, see the final
  follow-up below) once ablation showed neither was load-bearing** once the real energy-injection
  sources (ss6.2's `c_rot` cap; the broad-phase margin, ss4.2 step 2/ss9's 2026-09-26 entry) were
  fixed.

  Two more aggressive architectural changes were also tried and reverted after they regressed this
  crate's own existing, tuned regressions (see ss6.1/6.2 for each in place): removing ss6.1's
  ball-side overlap reaction entirely (regressed `metrics::tests::
  compression_error_stays_bounded_under_violent_lifter_cataracting`, 0.123 vs. the previously-
  passing ~0.091 bound), and capping ss6.2's entrained fluid mass `m_ent` at the ball's own physical
  added mass (broke `pbf::tests::ball_drag_matches_two_dimensional_stokes_scaling`, which
  specifically calibrates against the large-`m_ent` limit as the physically correct one). Regression
  test: `coupling::tests::low_fill_cataracting_charge_does_not_gain_energy_from_the_fluid`.

  **Measured effect** at the exact reported configuration: the ball population's total mechanical
  energy slope over a 5s measurement window went from **+16.9 W/m** (net energy gain — the "gas"
  failure mode) before this fix to **-11.1 W/m** (net dissipation, the physically expected sign)
  after. **Honest limitation:** peak instantaneous ball speed did *not* measurably improve in this
  same measurement (stayed in the ~4-5 m/s range either way; the ball speed clamp backstop above
  never engaged in this run) — this is a partial, validated mitigation of the energy-injection
  symptom specifically, not a complete fix of the visually reported "balls flying" symptom. The
  most promising next step, spot-measured to also reduce (though not eliminate) the energy-
  injection number but **not implemented as code in this pass**: auto-raising
  `simulation.resolution` so the fluid lattice spacing never exceeds the coarse-grained effective
  ball diameter (keep `dx/d_eff <= 1`).

  **Follow-up (2026-09-23): the auto-raise is now implemented.** `Params::effective_fluid_resolution`
  (`params.rs`, right after `Params::effective_media`) computes `min_resolution =
  ceil(mill.radius_m() / effective_media().diameter_m)` and returns
  `simulation.resolution.max(min_resolution).min(200)`; `Simulation::new` (`lib.rs`) now calls this
  method instead of reading `simulation.resolution` directly when seeding the fluid lattice
  (`pbf::FluidParticles::seed_lattice`). This was motivated by a second, related repro: raising
  `simulation.max_balls` to 1500 while leaving `simulation.resolution` at the Realtime preset's
  default of 15 (everything else at `Params::default()`) shrinks `d_eff` via coarse-graining until
  `dx/d_eff = 2.60`, well past 1; `effective_fluid_resolution` auto-raises the fluid lattice
  resolution actually used to seed the population to 40 in this case, bringing the ratio back down to
  `~0.98`. Measured at this exact repro: `coupling_clamp_hits` drops from 7944 to 13 over the same
  measurement window (~600x), and peak ball speed stops climbing (capped rather than still rising
  past 7.8 m/s). This was a no-op for all three official quality presets of that date (Realtime 150/15,
  Balanced 300/25, Accuracy 600/40; the presets today are 300/25, 300/30, 600/50 and 600/25
  "Impact statistics", `web/src/params/presets.ts`) since each already satisfies `dx <= d_eff` at
  its own defaults (docs/PERF.md's browser table) — it only engages once a custom
  `max_balls`/`simulation.resolution` combination (or a smaller `media.ball_diameter_m`) pushes the
  ratio past 1, and when it does engage it costs real extra fluid-solver time (`~resolution^2`
  particles), a disclosed correctness-over-a-previously-silent-broken-performance-assumption
  trade-off. Regression test:
  `tests::simulation_new_auto_raises_fluid_resolution_when_coarse_graining_shrinks_balls`.

  **Follow-up (2026-09-26): apparent regression, later found to be a false alarm.** After the
  coarse-graining-independence drag fix and the `c_rot` cap (ss6.2) resolved the *default*
  config's own, more severe energy-injection bug, this test's own config (`fill_fraction=0.10`)
  appeared to regress further under the same drag fix (ball-only energy_slope_w +16.9 ->
  **+26.6 W/m**, worse than the original "gas" baseline). An angular ball-speed backstop was
  added and brought the reading to 15.2 W/m (still over the `<5.0` bound), and tightening it
  further made it *worse* (19.5 W/m at half the ceiling) -- non-monotonic, the first sign this
  was noise, not a real defect.

  **Follow-up (2026-09-27): root-caused as a test-design bug, not a physics bug, and both ball
  speed backstops removed.** A ball-only mechanical-energy slope cannot distinguish "the coupling
  created energy" from "the fluid legitimately handed the ball some of the energy the wall gave
  the fluid" -- exactly the ball<->fluid drag this project intends, and step 6.2's correction
  being usually skipped here (degeneracy gate failing outright) doesn't change that. Checked with
  the same full-system (ball+fluid) per-sub-step energy invariant ss6.2's fix validates against,
  at this exact config, 2 seeds, 20s window: the invariant held at every sub-step both times
  (worst margins +1.78 J and +7.80 J), while the ball-only slope swung sign run-to-run and
  window-to-window (+13.8/+1.0/+4.3 W/m at 5/10/20s one run, -1.2/+3.9 W/m at 5/20s the other) --
  ordinary noise. This test now asserts the full-system invariant directly (`coupling::tests::
  low_fill_cataracting_charge_does_not_gain_energy_from_the_fluid`); no code in `pbf.rs`/`dem.rs`
  changed as a result.

  With no real defect behind either the linear or angular ball speed backstop (the linear one had
  never actually engaged in the measurements that motivated it either, see "Measured effect"
  above), both were ablation-tested: disabling each individually, then together, still passed the
  full 126-test suite and held the full-system invariant (worst margin still positive) at the
  default config, this config, and a 1500-ball high-count config. Both were removed outright
  (`dem::BALL_SPEED_SAFETY_FACTOR`, `ball_speed_clamp_hits`, and the pre-coupling fluid speed
  clamp duplicate call tried alongside them, ss5.2) rather than kept as unused insurance --
  restoring the pre-investigation NaN-sanitization-only behaviour (ss4.2 step 0) for ball
  velocity/spin.
- **History -- periodic media-charge oscillation ("surging"/"slumping") -- already reproducible,
  not a gap (2026-09-27 investigation; the 50 Pa*s, `max_balls` 150/600 and 1 m drum quoted in this
  entry are the defaults of that date; today's slurry default is water, 1 mPa*s).** A user report of real wet mills' charge visibly oscillating as a
  block at low speed (tens of `%Nc`), which this simulator was not observed to reproduce, was
  hypothesised to stem from `media.friction_ball_ball`/`friction_ball_wall` having no separate
  static-vs-kinetic coefficients (ss4.2 step 5's Coulomb clamp uses one `mu` for both the stick and
  slide regimes -- a non-anchored model with no cross-substep stiction memory, docs/PLAN.md step
  5). That hypothesis did not hold up: a new diagnostic binary,
  `crates/mill-core/examples/oscillation_probe.rs` (60 Hz sampling of charge centroid angle,
  toe/shoulder, and wall-slip ratio after a settle period; periodicity is assessed on the
  *linearly-detrended first difference* of the centroid angle, requiring a genuine
  trough-then-rebound autocorrelation shape -- not just a raw threshold crossing, which produces
  false positives from the charge's own smooth settle-in drift), found a clear, repeatable periodic
  oscillation **already present with `friction_ball_ball == friction_ball_wall`'s existing
  single-coefficient model, no code change**, at `mill.speed_value` in 20-30 `%Nc`, no lifters:
  - Dry, `max_balls = 150`: period 1.43 s. Dry, `max_balls = 600` (`Params::default`'s own count):
    period 1.35 s. Wet (`slurry.viscosity_pa_s = 50`, the project default) at `max_balls = 150`:
    period 1.32-1.42 s across 2 of 3 seeds (the third borderline). Measured period is essentially
    independent of `%Nc` (10/20/30 tested) and of wet vs. dry, and tracks a physical-pendulum
    estimate `2*pi*sqrt(R/g)` (R = the charge centroid's mean distance from the drum axis, ~0.30 m
    here) of ~1.10-1.11 s reasonably well (observed periods run ~20-30% longer, consistent with a
    damped, not undamped, pendulum). This points to the mechanism being a gravity-driven bulk
    "sloshing" mode of the charge's centroid about the drum's low point -- an emergent granular
    effect of the existing multi-contact XPBD solve, not something requiring a static/kinetic
    friction split to exist at all.
  - **Why it doesn't show up at this project's actual UI defaults** (`max_balls = 600`,
    `slurry.viscosity_pa_s = 50`): that specific combination sits in an over-damped regime. At
    `max_balls = 600` wet, the same probe found no clear periodic signal at the default 50 Pa*s
    (even extending the measurement window to 25 s), but lowering `slurry.viscosity_pa_s` alone (to
    5, 1, or 0.5 Pa*s, everything else at `Params::default()`) restored a clear, strong periodic
    signal at every one of those three values (periods 1.45-1.75 s). A secondary, smaller effect
    compounds this at low `max_balls`: `pbf.rs`'s ball<->fluid drag law is calibrated against the
    *true* (uncoarsened) ball radius (see the "true-radius drag fix" this crate's git history
    documents), so a coarse-grained population's larger, heavier simulated balls receive
    proportionally less drag relative to their own inertia than the true population would -- which
    is why `max_balls = 150`/`300` still oscillated at the full 50 Pa*s default while `max_balls =
    600` did not. Both are consequences of already-deliberate, documented modelling choices (a
    thick default slurry viscosity; drag calibrated to the true, not coarse-grained, ball size),
    not defects, so neither `dem.rs`/`pbf.rs`/`params.rs` nor any default value was changed as a
    result of this investigation.
  - **How to observe it in the running app**, no code change needed: set `slurry.viscosity_pa_s` to
    5 Pa*s or lower (or toggle `slurry.enabled` off), leave `mill.speed_value` around 20-30 `%Nc`
    with no lifters, and watch the existing "Toe angle" / "Shoulder angle" metrics-panel sparklines
    (`web/src/metrics/specs.ts`, already `sparkline: true`) -- both already surface the same
    vertical-degree convention `oscillation_probe.rs` samples.
  - **Large-amplitude recipe (2026-09-28 follow-up, still no code change).** The prior entry only
    established that a periodic mode *exists*; it did not target amplitude, and the default ball
    diameter has since moved to 63 mm (no coarse-graining) so period/amplitude numbers needed
    re-measuring. `oscillation_probe.rs` gained `--rolling-friction`/`--fill`/`--slurry-fill`/
    `--restitution-wall` flags and a "median swing per detected period" amplitude metric (distinct
    from the pre-existing whole-window `p2p`, which a single outlier swing or settle-in drift can
    inflate). A ~50-run sweep (wet, no lifters, default 63 mm balls) found: `slurry.viscosity_pa_s`
    dominates by far (only <=~7 Pa*s ever shows a periodic swing at all; the 50 Pa*s default is
    almost always flat), and within the periodic regime raising `media.friction_ball_ball` above its
    0.25 default consistently increases amplitude (more ball-ball grip makes the charge move more as
    a rigid block rather than shedding energy through internal rolling), while
    `media.friction_ball_wall` needs enough headroom above 0.35 to keep the mode from breaking back
    into non-periodic sloshing at low `%Nc`. Recipe periodic across all 3 tested seeds (1/2/3),
    amplitude (median swing per period) 6.1-9.5 deg, whole-window centroid p2p 12.2-15.9 deg, period
    ~1.3-1.6 s: `slurry.viscosity_pa_s = 5`, `mill.speed_value = 20` (`%Nc` mode), `lifters.count =
    0`, `media.friction_ball_wall = 0.6`, `media.friction_ball_ball = 0.5`,
    `media.rolling_friction` left at its 0.01 default. Reproduce with:
    `cargo run -p mill-core --release --example oscillation_probe -- --percent-critical 20 --slurry
    on --viscosity 5 --friction-ball-wall 0.6 --friction-ball-ball 0.5 --seed <1|2|3>`. The
    amplitude landscape is narrow and jaggy near its edges -- e.g. `%Nc = 40` with the same
    frictions gave a larger single-seed amplitude (12.5 deg) but lost periodicity on a third seed --
    so this `%Nc = 20` combination was chosen over higher-peak alternatives specifically for being
    periodic on every seed tested, not for the single highest amplitude observed. No
    `dem.rs`/`pbf.rs`/`params.rs` change was needed to reach this; static-vs-kinetic wall friction
    (which would add a mechanism for even more pronounced stick-slip and a wider stable parameter
    band) remains a possible future addition if a wider/more-robust surging band is ever needed, not
    a blocker for reproducing the effect today.
  - **Static/kinetic friction split added (2026-09-28), but a symmetric both-side "rocks through
    the drum's low point" swing was still not reproduced.** A user report described a real wet
    mill's charge swinging left and right by roughly equal angles (i.e. spending comparable time,
    and reaching comparable extents, on *both* sides of the drum's vertical low point), distinct
    from the one-sided wobble measured above (every centroid-angle trace so far sits at a
    steady-state offset toward the ascending/ "up" side, wobbling *within* that one side, never
    crossing back past vertical to the other). Two changes followed:
    1. `MediaParams` gained `friction_ball_ball_static`/`friction_ball_wall_static` (defaulting to
       their existing kinetic counterparts -- an old client's JSON deserializes to the exact prior
       single-coefficient behaviour) and `friction_velocity_scale_m_s`. `dem.rs`'s new
       `effective_friction(mu_kinetic, mu_static, v_t, velocity_scale)` blends the two with
       `mu_kinetic + (mu_static - mu_kinetic) * exp(-(v_t/velocity_scale)^2)` -- smooth (an even,
       everywhere-differentiable function of the contact's relative tangential speed `v_t`) rather
       than a hard Karnopp-style if/else switch, specifically chosen over a hard switch (an earlier
       version of this plan) to avoid injecting a velocity discontinuity into step 5's fixed-
       substep, per-contact-then-resync solve, which can make a contact chatter between "just
       stuck" and "just slipping" every sub-step at the switch boundary. Verified bit-exact for
       every shipped default via `examples/perf_probe.rs`'s hash oracle (identical hashes,
       `git stash`-compared same-session) -- this is a strictly additive, opt-in capability.
    2. `oscillation_probe.rs` gained `--friction-ball-wall-static`/`--friction-ball-ball-static`/
       `--friction-velocity-scale`, `--ball-diameter-mm`/`--drum-diameter-mm` (to test other mill
       scales), and a "vertical-crossing symmetry" report line (extents of the raw centroid-angle
       offset from the 180 deg low point on each side, and whether the trace crosses it at all).
    Swept broadly (all wet, no lifters): `%Nc` 1-80 at the shipped 63 mm ball / 1 m drum scale,
    `%Nc` 10-50 at a 2 mm ball / 63 mm drum scale (the user's own earlier small-scale experiment --
    note this scale needs `simulation.resolution` capped low, e.g. 10, or the auto-raised fluid
    resolution makes a CLI sweep impractically slow), viscosity from thick (50 Pa*s) down to real
    water's (~0.001 Pa*s), and -- with the new split -- `friction_ball_wall_static` up to 2.0 (a
    physically extreme value; most real material pairs' static coefficient is below ~0.8) crossed
    with `friction_ball_ball_static` up to 1.5. In every one of these runs (~100 total), the
    centroid-angle trace's `side_a_extent`/`side_b_extent` symmetry metric came back one-sided
    (`crosses_bottom=false`): the mean position always sits 10-35 deg to one side of vertical, and
    the oscillation on top of that mean stays confined to a 3-15 deg window without ever crossing
    back through the low point to the other side. The likely reason: a *continuously rotating* wall
    inherently carries the charge toward the ascending side (that carry is exactly what drives
    cascading at all), so the charge's time-averaged position cannot sit at vertical -- the
    static/kinetic split changes how abruptly the charge grips and releases the wall, and how large
    the resulting wobble is, but doesn't remove that one-sided bias, because the bias's cause (net
    angular momentum injection from a one-directional wall) is unrelated to which friction law is
    used. Tested down to `%Nc = 1` (near-zero net rotation) with no different outcome -- the charge
    just settles into a small, still-one-sided jitter rather than swinging. No default values were
    changed by this investigation; the static/kinetic friction capability is available for future
    use (e.g. a UI control, or lifter-driven surging).
  - **Free-pendulum release tried next (still 2026-09-28): decisively overdamped, no crossing at
    all.** `oscillation_probe.rs` gained `--stop-after-settle`: runs the usual driven settle phase
    (lifting the charge up the ascending side as normal), then calls `Simulation::set_params` to
    set `mill.speed_value = 0` (wall stops dead, mid-run, no reseed) and measures what the charge
    does purely under gravity + its own momentum from then on -- a genuine release, unlike every
    other measurement in this section which keeps the wall actively driving the charge throughout.
    Tried from both a modest initial displacement (settled at 30 `%Nc`, released ~15 deg off
    vertical) and a large one (settled at 70 `%Nc`, released ~28 deg off vertical), both at real
    water's viscosity. In both cases the centroid angle **crept monotonically back toward vertical
    over several seconds and never overshot past it even once** (e.g. the 70 `%Nc` case: released
    at 152.6 deg, still only at 177.3 deg -- short of 180 -- after a further 8 s) -- not a damped
    oscillation with decaying overshoot, a plain overdamped relaxation with no overshoot at all.
    At the project's default friction (0.35 wall / 0.25 ball) this looked like it ruled out "just
    release it" entirely -- but see the next entry, which isolates *why* and finds this is a
    friction-magnitude threshold, not a dead end.
  - **Root cause isolated (still 2026-09-28): kinetic friction alone controls whether a release
    swing crosses the low point -- and it does, well below this project's shipped defaults.**
    `oscillation_probe.rs` gained `--release-friction-ball-wall`/`--release-friction-ball-ball`/
    `--release-rolling-friction`/`--release-restitution-wall`/`--release-restitution-ball`:
    overrides applied at the same instant as `--stop-after-settle`'s wall stop, so the *settle*
    phase still uses normal friction (the wall needs real friction to lift/displace the charge at
    all -- see below) while the *release* itself can use different coefficients, letting the
    "what displaces the charge" and "what damps its swing" questions be answered independently.
    Sweeping the release-phase Coulomb friction coefficient (`friction_ball_wall`/
    `friction_ball_ball`, applied equally, settle phase left at the 0.35/0.25 defaults, `slurry
    off`, released from ~27 deg at 70 `%Nc`, seed 1) found a clean threshold: **`crosses_bottom =
    true` for release friction `<= ~0.1`, `false` for `>= ~0.15`-`0.2`** -- e.g. at `mu = 0.05`:
    `symmetry_ratio = 0.55`, p2p 45 deg, a clear multi-cycle decaying swing (153 -> 186 -> 172 ->
    198 -> 160 deg ... settling near 180 deg by ~14 s, confirmed by inspecting the raw trace); at
    `mu = 0.005` (near-frictionless): `symmetry_ratio = 0.88`, amplitude (first-period) 35.7 deg.
    **Restitution and rolling friction turned out not to matter**: re-running `mu = 0.05` with
    release restitution swept 0.7 (the project's own default) through 0.99 gave near-identical
    results every time (`symmetry_ratio` 0.52-0.55 throughout), and leaving `rolling_friction` at
    its 0.01 default instead of forcing it to 0 changed nothing either -- Coulomb (kinetic) friction
    is the sole controlling variable for this mode's damping ratio in this solver, confirming the
    "effective damping ratio `>= 1`" read above was specifically about friction, not restitution.
    **But this crossing behaviour is release-only, not a driven-rotation phenomenon**: re-running
    the driven-rotation sweep (no stop, continuously rotating wall) at this same low friction
    (`mu = 0.05`, `%Nc` 5-15, water viscosity) still came back `crosses_bottom = false` every time
    -- a continuously rotating wall injects net one-directional angular momentum regardless of how
    low the friction coefficient is (friction only sets *how efficiently* that momentum couples
    in, not its direction), so low friction alone cannot turn the driven steady state symmetric --
    it only unlocks the underdamped free-release mode. **Net picture**: a genuine, large,
    left-right-symmetric swing through the drum's low point is achievable in this model, but only
    as a *decaying transient right after the drum stops* (e.g. an operator cutting power, not
    steady-state running), and only with ball-wall/ball-ball kinetic friction well below this
    project's ceramic-media defaults (roughly `<= 0.1`, vs. the shipped 0.35/0.25) -- consistent
    with a much smoother/harder media-and-liner pairing than the default YSZ-on-steel model.
    Reproduce: `cargo run -p mill-core --release --example oscillation_probe --
    --percent-critical 70 --slurry off --settle-s 3 --measure-s 10 --stop-after-settle
    --release-friction-ball-wall 0.05 --release-friction-ball-ball 0.05 --seed 1`.
  - **User rejected the release-only finding (2026-09-29): needs the crossing swing DURING
    continuous rotation specifically, not as a stop-transient.** An extensive further search of
    driven-rotation (never-stopped) mechanisms all failed to produce a robust `crosses_bottom =
    true`: high static + low kinetic wall friction (the classic "climbs while stuck, breaks free,
    slides, overshoots" mental model -- still one-sided, and *smaller* amplitude than the Phase-1
    recipe); low ball-ball friction with normal/high wall friction (a detached avalanche layer
    flowing freely); fill fraction 0.5-0.7 (rigid-block precession); a 2:1 parametric-resonance
    check (`%Nc` 55-75, where the rotation period is close to twice the pendulum period -- no
    resonance effect, bias angle stayed pinned at ~163.5 deg across the whole range); and a
    combined wall-stick-slip + low-ball-ball-friction avalanche. A structural explanation for why:
    a continuously rotating wall's "stuck" state ties the ball to the *wall's own nonzero,
    constantly-refreshed speed*, not to rest -- every time the slip speed re-approaches zero the
    ball re-locks near wall speed rather than continuing to build independent momentum, capping
    achievable swing energy near `~(wall speed)^2` regardless of the friction law's shape (this is
    different from the release case, where `v_t = 0` means the wall itself, and thus the "target"
    speed to re-lock to, is at rest). Two alternative model extensions were proposed (a
    torque-limited/compliant-motor dynamic-drum model; an azimuthally wall-friction-textured
    liner) and both declined by the user as too specialized a modelling assumption.
  - **Literature search (2026-09-29, WebSearch): the classic single-species "slumping regime" is
    exactly this model's one-sided result; a genuinely symmetric driven oscillation in the
    literature needs a second granular fraction this simulator doesn't have.** The rotating-drum
    granular-physics literature describes the classic low-Froude-number "slumping regime" (e.g.
    Rajchenbach-style avalanche dynamics) as the free-surface slope oscillating *between the angle
    of repose and the angle of marginal stability* -- i.e. the same one-sided sawtooth this
    project's DEM already reproduces (Phase 1's recipe); it is not described as crossing the
    drum's low point either. A real, literature-documented *symmetric* "self-oscillating"/
    "auto-excited" pulsation during continuous rotation exists (Deineka et al.,
    Eastern-European Journal of Enterprise Technologies, "Establishing the rotation speed
    variation range limits for auto-excitation of self-oscillating grinding in a tumbling mill";
    "Revealing the mechanism of stability loss of a two-fraction granular flow in a rotating
    drum"), but its driving mechanism is **dilatancy between two granular fractions** (coarse
    grinding media + fine ground material/ore), not wall friction at all -- this simulator has no
    second solid granular phase (only media balls + a continuum SPH liquid). The user's own
    reference, however, was water + media only (no ore/fines), which this specific mechanism does
    not explain either. A related, independently-sourced term from wet-milling literature --
    "sloshing" (the charge moves as a liquid rather than tumbling as a coherent granular mass when
    slurry solids content is below ~70 wt%) -- prompted a follow-up test: low media fill
    (0.05-0.15) with high slurry fill (0.6-0.8), tracking the *slurry population's own centroid*
    (`oscillation_probe.rs`'s new `fluid_centroid`, independent of the sparse ball population).
    Result: weak and inconsistent (`symmetry_ratio` only 0.07-0.39 across the configs that crossed
    at all, vs. 0.55-0.88 for the genuine release finding above) -- not a robust phenomenon, though
    not conclusively ruled out either (only the fluid's mass-weighted centroid was checked, not its
    free-surface extent, which a real sloshing wave crest might show more clearly). **Status: open,
    unresolved** -- no driven-rotation mechanism tried so far reproduces a robust symmetric
    low-point-crossing swing; see the project memory entry for this investigation for the
    prioritized list of what to try next.
- **Slurry does not seep into a settled/stationary ball bed, even at low viscosity with the drum
  stopped (2026-09-27 investigation, UI note added, no solver change; measured at the 10 mm /
  `max_balls = 150` / `resolution = 15` configuration of that date).** A user report: after
  lowering `slurry.viscosity_pa_s` and setting `mill.speed_value = 0`, the slurry pool visible above
  the settled charge does not drain down into the charge's own interstitial gaps over time, even
  waiting tens of seconds. A new diagnostic binary, `crates/mill-core/examples/infiltration_probe.rs`,
  measures the fraction of fluid particles that end up "covered" (a ball sits directly above them,
  i.e. they are physically underneath/between balls rather than exposed to the open pool) over a
  long idle period after the drum stops. Two compounding, already-documented modelling limitations
  explain the observation; **no bug was found and no defaults were changed.**
  - **Numerical resolution: `dx` is comparable to the ball's own size, not to the much smaller gaps
    between packed balls.** At a repro matching the report (`max_balls = 150`, `resolution = 15`,
    `slurry.viscosity_pa_s = 1`, static settle from `t = 0`, no rotation), the fluid lattice spacing
    is `dx = 0.0333` m against an effective (coarse-grained) ball diameter `d_eff = 0.0405` m --
    `dx/d_eff = 0.82` -- comfortably inside `Params::effective_fluid_resolution`'s only guarantee
    (`dx <= d_eff`, sized for ball<->fluid coupling *stability*, ss3/ss6, not for resolving pore
    throats). A single fluid particle is nearly as wide as an entire ball, which cannot fit through
    a gap between two touching balls (a gap that is at most a small fraction of the ball's own
    radius). Measured "covered" fraction plateaus at **41%** within about 15-20 s of settling and
    does not move further given another 10 s. Re-running the same repro at 4x finer resolution
    (`resolution = 60`, `dx/d_eff = 0.21`, `dx` now well under half the ball radius) raises the
    plateau to **54%** (reached by `t = 10` s, essentially flat through `t = 30` s) -- confirming
    resolution is a real, measurable factor, not a red herring.
  - **2D cross-section topology and voidage are the second, resolution-independent factor.** Even at
    4x finer resolution the plateau stayed well under "fully filled". This is consistent with ss3's
    "not the true interstitial void structure" note and ss9's own "2D areal packing is not 3D
    voidage" entry above: a 2D packed-circle bed has both a lower void fraction (`1 -
    packing_fraction_2d`, ~18% at the project's default `0.82`, vs. ~36-40% for a random 3D sphere
    packing) and far fewer alternate percolation paths around any one near-contact than a real 3D
    bed has (in 2D, two touching circles' point of tangency is the *only* route between the void
    pockets on either side of it; in 3D the same local near-contact is bypassed from many more
    directions). Raising `media.packing_fraction_2d` down toward the 3D range does **not** fix this:
    it only rescales `Params::true_ball_count`'s solid-disc count for a given `fill_fraction` (and
    hence the "Interstitial filling (U)" volume bookkeeping, ss3), not the local gap geometry between
    whichever balls actually end up touching after settling -- those still pack as densely as XPBD
    contact resolution allows, regardless of how many of them there are.
  - **Conclusion:** given the slurry's own `fill_fraction` is typically comparable to or larger than
    the charge bed's own void volume, most slurry is expected to remain visible as a pool above the
    settled bed in this simulator, independent of viscosity or whether the drum is turning -- this is
    consistent, expected behaviour of a coarse-grained 2D model, not a solver defect. A UI note was
    added (`web/src/ui/paramsPanel.ts`'s Slurry group) explaining this so it is not mistaken for a
    bug; no `dem.rs`/`pbf.rs`/`params.rs` change or default-value change was made.

### 2 mm / 63 mm shipped default: measured side effects (2026-09-29)

The default became a 2 mm ball in a 63 mm drum (118 rpm ~ 70% Nc at that date, `N_true` ~ 244,
Realtime `max_balls` 150 -> 300; today's default speed is 100 rpm ~ 59% Nc, slurry is water at
1 mPa·s, and the quality presets are listed in ss3.1). Measured consequences on that date, disclosed
rather than tuned away:

- `Params::effective_substeps` now returns 16 (the cap) at every quality tier: `D / d = 31.5`
  needs ~19 sub-steps for the 0.5 displacement-ratio target, so the shipped default costs ~2x the
  former 8 sub-steps per frame.
- Transient ball-ball overlap at `dem_iterations = 2` is 0.3-0.5 of a radius (dry and static drum
  alike, so it is DEM contact resolution of a dense pile); `dem_iterations = 6` gives ~0.12-0.15.
- Wet: the maximum PBF compression error is 6-12 (sustained over 150 frames, not growing) --
  single fluid particles pinched between the tiny balls; see the infiltration note above.
  Energy-invariant tests still pass.

---

## Notes for the maintainer

Points found while reconciling this document with the source (2026-10-20). They are flagged, not
decided.

1. **PRNG constants (`rng.rs:6`-`rng.rs:34`).** `Rng::next_u64` applies Marsaglia's xorshift64
   shifts `(13, 7, 17)` and then multiplies the new state by `0x2545F4914F6CDD1D`, the multiplier of
   Vigna's xorshift64\*. The published xorshift64\* generator uses the shift triple `(12, 25, 27)`,
   not `(13, 7, 17)`, so the code is a hybrid that is not xorshift64\* as published; its doc comment
   says exactly this ("xorshift64 steps (13, 7, 17) followed by the xorshift64\* multiplier"), while
   other documents (`docs/PLAN.md`, `docs/TECHNOLOGY.md`) call it a xorshift generator /
   xorshift-with-multiplier. In the crate it is used only for the ball-lattice jitter
   (`Balls::seed_lattice`, `dem.rs:361`; the fluid lattice is deterministic), so statistical quality
   is not load-bearing, but changing the constants would change every seeded initial condition and
   invalidate all bit-exact baselines
   (`examples/perf_probe.rs` hashes, tolerance-probe references). Decision needed: keep and document
   as is, or switch to `(12, 25, 27)` and re-baseline.
2. **Angular momentum in the joint ball–fluid solve (`pbf.rs:627`-`pbf.rs:639`, ss6.2).** Verified
   from the equations: the converged solve conserves total linear momentum (uniform translation is a
   null mode of every coupling term) and gives each ball exactly the torque `l_k x F_k` of the forces
   at its ghost points; the *total* angular momentum about the drum origin is conserved only up to
   the offset between the fluid particle (where the force on the fluid acts) and the ghost point
   (to which the ball's reaction is referred). The doc comment on `solve_viscosity_with_balls` now
   states this ("each ball's torque exactly, total angular momentum only up to the ghost-point
   offset"); this document says the same. Whether that offset error matters at the shipped
   resolutions was not measured here.
3. **Stale code comments about the removed Stokes term (resolved 2026-10-20, comments only).** The comment on
   `F_CLAMP_G_MULTIPLE` (`pbf.rs`) and the comments of `cascading_charge_keeps_coupling_clamp_hits_rare_once_settled_at_*`
   (`coupling.rs`) now say that the value `60` and the 10 % bounds date from the Stokes-era closure and were not re-derived;
   they are kept as measured bounds. Still open (decision: keep): re-deriving `60` for the current closure.
4. **`RESTITUTION_VELOCITY_THRESHOLD`'s vestigial `allow(dead_code)` (resolved 2026-10-20)** was removed (`dem.rs`).
5. **`Params::default()` versus the UI default.** `Params::default()` has `simulation.resolution =
   40`, while the web's default quality preset ("Fast (default)") sets `resolution = 25` and
   `max_balls = 300` (`web/src/params/presets.ts`). Numbers in this document that depend on the
   fluid resolution use `Params::default()` unless stated.

**Could not be verified (and why).** Every dated measurement in ss6.2.1 and ss9 (drag fractions,
energy-invariant margins, oscillation periods, infiltration plateaus, clamp rates) is quoted from the
earlier text or from `docs/VERIFICATION.md` and was not re-run; `cargo` was not run, so whether the
tests named here (for example `ball_drag_matches_two_dimensional_stokes_scaling`, whose name is
historical, ss6.2) pass was not checked, only that each named test exists; the figures "+17...39 %
power change under coarse-graining" and "6 of 22 strict gates pass" are taken from the comment in
`params.rs:38`-`params.rs:40` and `docs/VERIFICATION.md` respectively; the statement that the UI
marks per-impact readings unavailable under coarse-graining was spot-checked in
`web/src/metrics/specs.ts` only.
