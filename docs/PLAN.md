# MillDynamics3 — Ball Mill Media/Slurry Simulator (Rust → WASM, real-time 2D)

## Context

Goal: an interactive browser app that simulates a tumbling ball mill cross-section — grinding media (balls) and a
viscous slurry — and animates it in near real time. The two things the user must be able to see are
(1) where the balls are and (2) where the slurry free surface is. Mixing of the slurry should also be visible.
Optional lifters (wall protrusions) must be supported, but the default drum wall is perfectly smooth.
Parameters are edited in a modal window. Everything runs client-side in WebAssembly.

Decisions already made with the user:
- **2D cross-section** (plane perpendicular to the mill axis), not 3D.
- **Rust → WASM** core (Rust is not installed yet; Node 24 / npm 11 / git / MinGW gcc / VS present).
- **Custom solver core** (no rapier/salva): DEM for balls + Position Based Fluids (PBF) for slurry.
- **Vite + TypeScript, no framework**; Canvas 2D rendering; native `<dialog>` modal.
- Conversation with the user in Japanese; *all* code, comments, docs, UI text, commits in English.
- Fable/Opus only for planning and hard parts; Sonnet/Haiku subagents for the rest (recorded in CLAUDE.md).

### Why DEM + PBF instead of pure coarse-grained DEM (or SPH / LBM / MPM)

| Method for slurry | Verdict |
|---|---|
| Coarse-grained DEM "liquid" particles with damping/cohesion | Single solver, but no incompressibility → slurry compacts/sinks, viscosity is uncontrolled, free surface is ill-defined. Rejected. |
| WCSPH (weakly compressible SPH) | Physically faithful viscosity, but dt ≈ 1e-4 s → ~170 substeps per 60 Hz frame → ~10× too slow for real time in single-thread WASM at useful resolution. Kept as a possible high-fidelity mode later. |
| **PBF (Position Based Fluids)** | Unconditionally stable at dt = 1/240 s, 3–4 iterations; ~10–20× cheaper than WCSPH; clean free surface; well-proven for real-time (NVIDIA Flex). Viscosity is qualitative (XSPH) in v1, upgradeable to implicit viscosity. **Chosen.** |
| LBM free-surface | Grid-based, free-surface tracking and moving lifters are complex. Rejected. |
| MPM | Great for non-Newtonian, but grid+particles transfer cost is too high for real time here. Rejected. |

### Why position-based (XPBD) rigid discs instead of explicit soft-sphere DEM, for the media

The original plan (this section, earlier revision) called for an explicit spring–dashpot soft-sphere
DEM for the balls, matching common DEM practice. While deriving the M1 contact constants this was
found to be **incompatible with real time even after coarse-graining**: an explicit spring–dashpot
contact's stable time step scales as `t_c ~ pi * (allowed overlap) / (impact speed)`. Keeping the
overlap within a realistic 1% of the ball radius at the drum's actual wall speed (~3 m/s at the
defaults) gives `dt_dem` on the order of **microseconds** — tens of thousands of sub-steps per
rendered frame — regardless of how much coarse-graining (ss3.2) enlarges the effective ball
diameter, since bigger, lighter effective balls don't relax that overlap/speed ratio enough at
achievable `max_balls` values. This is the same instability-vs-stiffness trade-off that ruled out
WCSPH for the fluid above, so the same fix applies:

| Method for grinding media | Verdict |
|---|---|
| Explicit soft-sphere DEM (linear spring–dashpot, Cundall–Strack) | Physically standard, but stable `dt` is set by contact stiffness/impact speed and lands in the microsecond range at real mill wall speeds — tens of thousands of sub-steps/frame, independent of coarse-graining. Rejected. |
| **Position-based rigid discs (XPBD-style)** | Non-penetration solved as a rigid (zero-compliance) geometric constraint, projected iteratively — unconditionally stable at any sub-step size, including the same 1/240 s / `substeps`-per-frame rate already used for PBF. Friction, restitution, and rolling resistance are added as position/velocity corrections on top of the converged contact solve. **Chosen.** |

Balls are therefore solved with the **same fixed sub-step as the fluid** (`crate::FIXED_DT` =
1/240 s, `simulation.substeps` per rendered frame), removing the separate DEM time-step/stiffness
concept entirely and simplifying the M4 coupling loop (ss3.4) to a single shared sub-step for both
solvers. Ball charge behaviour (cascading/cataracting/centrifuging, toe/shoulder) is still the
primary output and this method reproduces it: it is the same class of solver used for real-time
rigid body contact in production physics engines (e.g. NVIDIA PhysX, Rapier), applied here to 2D
discs.

---

## 1. Repository layout

```
MillDynamics3/
  CLAUDE.md                     # project rules (language + model-delegation policy) — see §2
  README.md                     # how to build/run
  .gitignore
  docs/
    PLAN.md                     # copy of this plan
    PHYSICS.md                  # equations, units, parameter derivations (written with M1/M3)
    PARAMETERS.md               # every UI parameter: symbol, unit, range, default, hot-swappable?
  Cargo.toml                    # workspace
  crates/
    mill-core/                  # pure Rust simulation, no wasm deps; unit tests + benches
      Cargo.toml
      src/
        lib.rs                  # pub API: Simulation::new(&Params), step(dt_wall), views
        params.rs               # Params struct (serde), Default impl, validation, derived values
        geometry.rs             # drum SDF (circle + optional lifters), wall velocity, rotating frame
        grid.rs                 # uniform grid spatial hash (cell list) shared by DEM & PBF
        dem.rs                  # balls: state, contact model, integration, wall/lifter contacts
        pbf.rs                  # fluid: kernels, density constraint, XSPH viscosity, boundaries
        coupling.rs             # two-way ball↔fluid exchange
        surface.rs              # density grid splat + marching squares → free-surface polylines
        metrics.rs              # toe/shoulder angles, slurry level, mixing index, energy, power
        rng.rs                  # small xorshift RNG (deterministic seeds)
      benches/step.rs           # criterion benchmark (500/1000/2000 balls × 2k/4k/8k fluid)
    mill-wasm/                  # thin wasm-bindgen wrapper exposing mill-core
      Cargo.toml
      src/lib.rs
  web/
    package.json  vite.config.ts  tsconfig.json  index.html  styles.css
    src/
      main.ts                   # boot: load wasm worker, wire UI, rAF render loop
      worker.ts                 # owns the wasm Simulation; fixed-step accumulator; posts frames
      protocol.ts               # typed messages main↔worker (Init, SetParams, Pause, Frame, Stats)
      params/schema.ts          # parameter definitions (id, label, unit, min/max/step, default, group, hot)
      params/presets.ts         # named presets (Lab mill 0.3 m, Pilot 0.6 m, Lifters ×8, High viscosity…)
      ui/paramsModal.ts         # builds <dialog> form from schema; validate; Apply/Reset/Cancel
      ui/hud.ts                 # sim time, rpm, %Nc, fps, ms/step, achieved time scale, metrics
      ui/toolbar.ts             # Play/Pause/Step/Reset/Parameters/Presets/Screenshot buttons
      render/canvas.ts          # Canvas 2D renderer (drum, lifters, surface, fluid, balls, dye)
      render/colors.ts
      state.ts                  # app state (params, running, latest frame)
    tests/
      unit/*.test.ts            # vitest: schema validation, protocol, surface polygon helpers
      e2e/*.spec.ts             # playwright: modal opens, edit rpm, apply, canvas updates, lifters
  scripts/
    build-wasm.ps1 / build-wasm.sh   # wasm-pack build --target web --release → web/src/wasm/
```

Rust workspace uses only: `glam` (Vec2 math), `serde`/`serde_json` (params), `wasm-bindgen`, `js-sys`,
`criterion` (dev). No physics libraries.

---

## 2. CLAUDE.md (to be created verbatim in step M0)

```markdown
# MillDynamics3 — Project Rules

## Language policy
- Conversation with the user: **Japanese**.
- Everything else — code, comments, identifiers, commit messages, docs, UI strings, file names,
  plan files, subagent prompts, test names — **English only**.

## Model / cost policy
Fable and Opus are very expensive. Use them only for:
- planning and architecture decisions,
- the hard numerical core (DEM contact model, PBF solver, ball–fluid coupling, stability/perf bugs),
- reviewing changes to `crates/mill-core/src/{dem,pbf,coupling}.rs`.

Delegate everything else to subagents via the Agent tool with an explicit `model`:
- `sonnet`: UI (web/), rendering, build/tooling config, tests, docs, benches, wasm bindings,
  geometry/surface/metrics modules, refactors, moderate bug fixes.
- `haiku`: boilerplate, file scaffolding, formatting, running commands, simple lookups, small edits,
  writing/updating docs from existing content.
When in doubt, start with `sonnet`; escalate to the main (Fable/Opus) session only if the subagent
fails twice or the problem is numerical/physical.

## Project summary
2D cross-section tumbling ball mill simulator: DEM balls + PBF slurry, Rust → WASM, Vite + TS
frontend, Canvas 2D, native <dialog> modal for parameters. Default drum wall has **no lifters**.
See docs/PLAN.md, docs/PHYSICS.md, docs/PARAMETERS.md.

## Repository
- Remote: https://github.com/ToshihiroIguchi/MillDynamics3 (branch `main`). Commit per milestone; push only when the user asks.

## Conventions
- SI units everywhere in the core (m, kg, s, Pa·s). UI may show mm / rpm and convert in `schema.ts`.
- Deterministic: every run is reproducible from `Params` + `seed`.
- `cargo fmt`, `cargo clippy -D warnings`, `npm run lint` must pass before finishing a task.
- Do not commit `web/src/wasm/` build output or `target/`.
```

---

## 3. Physics specification (2D)

Units SI. World frame: drum center at origin, gravity −y, drum rotates counter-clockwise at ω = 2π·rpm/60.

### 3.1 Drum geometry (`geometry.rs`)
- Radius `R` (from diameter `D`). Critical speed `N_c = 42.3/√D` rpm (D in m); UI accepts `rpm` or `% N_c`.
- **Signed distance function in the drum's rotating frame**: `sdf(p_local) = min(R − |p|, sdf_lifters(p_local))`
  where positive = inside free space. Lifters (default **count = 0** → circle only): `n` identical bars evenly
  spaced, each a convex quad defined by `height`, `base_width`, `top_width` (trapezoid, face angle derived),
  optional `phase` angle. `sdf_lifters` = −(max over lifter quads of convex-polygon SDF) (i.e., subtract lifters from free space).
- `wall_velocity(p) = ω × p` (rigid rotation), used for friction/no-slip.
- Contact query: transform particle into drum frame (rotate by −θ_drum), evaluate `sdf` and its numeric gradient
  (central differences, ε = 1e-4 R) → normal; rotate back. Only evaluated for particles with `|p| > R − r − h_margin`
  or inside the lifter annulus `|p| > R − lifter_height − r − margin` to keep cost O(surface particles).

### 3.2 Balls — position-based rigid discs (`dem.rs`)

**Coarse-graining (particle scaling), implemented in `params.rs` (`Params::effective_media`).**
At the current defaults (D = 1 m, d = 2 mm, J = 0.30) the true 2D media population is ~75,000
balls, far above what a single-threaded WASM solver can step in real time even with the
unconditionally-stable XPBD approach below (the cost is per-contact, not per-substep). Rather than
exposing this as a raw performance cliff, the solver is always seeded from an **effective** media
population instead of the raw UI values: if the true ball count
`N_real = J·(π R²) / (π r_true²)` exceeds `simulation.max_balls` (default 2000, chosen to meet the
M6 "≥1.0x real time" target), it is replaced with `N_sim ≈ max_balls` larger, lighter balls using
scale factor `k = sqrt(N_real / max_balls)`:
- `d_eff = d_true · k` — preserves total footprint area (`N_sim·π r_eff² ≈ N_real·π r_true²`), i.e.
  the fill fraction the user set.
- `ρ_eff = ρ_true / k` — preserves total charge mass (a 3D-sphere mass `∝ r³` grows by `k³` while
  the count shrinks by `k²`, net `k`, canceled by dividing density by `k`).

This is the standard coarse-grained DEM approximation (see e.g. Sakai & Koshizuka 2009; Bierwisch
et al. 2009): it reproduces bulk charge behaviour (cascading/cataracting/centrifuging, toe/shoulder
angles) but not the true interstitial void structure or single-collision statistics at the real
particle size. `d_eff`/`ρ_eff`/`N_sim`/`k` are always shown in the parameters modal's derived-values
panel (ss4.3) so the approximation is never silent; `k = 1.0` (no coarse-graining) whenever
`N_real <= max_balls`. Mass/radius/inertia below are the *effective* (post coarse-graining) values.

**Solver.** No contact stiffness or DEM-specific time step exists in this design (see the rationale
above): balls advance on the same fixed sub-step as the fluid (`FIXED_DT` = 1/240 s,
`simulation.substeps` per rendered frame, `simulation.dem_iterations` constraint-solver iterations
per sub-step, analogous to PBF's `pbf_iterations`). Per sub-step `dt`:

1. **Predict**: `v += g·dt`; `x += v·dt`; `θ += ω·dt` (save `x0`, `θ0` for the velocity
   reconstruction in step 4). State per ball: `x, v (Vec2), θ, ω (Vec<f32> SoA)`, plus the uniform
   effective `r, m, I = 2/5·m·r²` from `effective_media` (size distributions are future work).
2. **Broad-phase**: `grid.rs` uniform grid over predicted `x`, cell = `2·r_eff`; gather candidate
   ball–ball pairs once (reused across iterations).
3. **Solve non-penetration** (`dem_iterations` Gauss–Seidel passes, zero compliance = rigid):
   - Ball–ball: `C = |x_i − x_j| − (r_i + r_j)`; if `C < 0`, `Δλ = −C / (w_i + w_j)` (`w = 1/m`),
     apply `Δx_i = w_i·Δλ·n̂`, `Δx_j = −w_j·Δλ·n̂`, accumulate `λ_n` for this contact this sub-step.
   - Ball–wall/lifter: `C = sdf(x_i) − r_i` against [`crate::geometry::Drum`]; same projection with
     `w_wall = 0` (infinite wall mass), accumulating `λ_n` per ball.
4. **Reconstruct velocities**: `v = (x − x0)/dt`, `ω = (θ − θ0)/dt`.
5. **Friction** (one pass, Coulomb-clamped position correction, not a persistent tangential
   spring): for each contact with `λ_n > 0`, `v_t = (v_i−v_j)·t̂ − r_iω_i − r_jω_j` (wall case: `v_j`
   replaced by the wall's rigid velocity at the contact point, no `r_jω_j` term); the correction
   that would fully cancel sliding this sub-step is clamped to `μ·λ_n` before being applied to
   position (linear) and orientation (`∝ r·w_rot`), then velocities are reconstructed again from
   the (now friction-corrected) positions. This is a simplified, non-anchored Coulomb model (no
   cross-substep stiction memory); acceptable for the bulk charge-motion behaviour this project
   targets, revisited if validation shows drift.
6. **Restitution** (one pass): for contacts with `λ_n > 0` whose *pre-solve* normal velocity was
   approaching faster than a small threshold (avoids resting-contact buzz), add back an extra
   normal impulse so the *post-solve* separating velocity matches `e · v_n,pre` (per-contact `e`:
   `restitution_ball_ball` or `restitution_ball_wall`).
7. **Rolling resistance**: `Δω = −sign(ω)·min(μ_r·(λ_n,total/dt)·r/I·dt, |ω|)` per ball, where
   `λ_n,total` sums this ball's accumulated normal impulses (ball–ball + wall) this sub-step —
   capped so it cannot reverse the sign of `ω` in one sub-step.

This is the same family of solver used for real-time rigid-body contact in production physics
engines (XPBD / "small steps" style, e.g. Macklin, Müller & Chentanez 2016; Müller et al. 2020),
applied here to 2D discs sharing the fluid's uniform grid ([`crate::grid`]) and fixed sub-step.

### 3.3 Slurry — PBF (`pbf.rs`)
- Particle spacing `dx` from `resolution` param (`dx = R / res`, default res = 40 → 80 particles across the drum);
  kernel radius `h = 2·dx`; rest density `ρ0 = ρ_slurry`; particle mass `m_f = ρ0·dx²` (2D, unit depth).
- Fill: slurry volume given as **fraction of drum area** `U_s` (or "% of charge voids" as an alternative input) → initial
  particles on a hexagonal lattice in the bottom of the drum (interstitial with balls; balls initialised first, fluid particles overlapping a ball are removed).
- Kernels: 2D Poly6 (`4/(π h⁸)`) for density, 2D Spiky gradient (`−30/(π h⁵)`) for ∇W.
- Substep (dt = 1/240 s, `substeps` param): predict `x* = x + dt·(v + dt·g)` → neighbour grid (cell = h) → `iters` (default 3,
  Jacobi-style: all λ/Δp computed from the current positions, then applied together) of density constraint
  `C_i = ρ_i/ρ0 − 1`, `λ_i = −C_i/(Σ|∇C|² + ε)` (ε = 200; plan default), `Δp_i = 1/ρ0·Σ(λ_i+λ_j+s_corr)∇W` → boundary
  projection → update `v = (x*−x)/dt` → no-slip blending → viscosity → `x = x*`.
- **Artificial pressure `s_corr`, disabled in v1.** The literature-default `s_corr = −k·(W(r)/W(Δq))⁴` (k = 0.1, Δq =
  0.2h) is implemented but set to `k = 0`: at this project's real-SI-unit scale (ρ0 ~1000-2000 kg/m³, small h), k = 0.1
  produced a correction 10-30x larger than the density-constraint terms it's meant to supplement, causing runaway
  dispersal instead of preventing clustering. Disabling it gives a stable settled puddle (~1% mean density error over
  2s simulated); several smaller `k` values tried did not find a working point in the time available. Revisit only if
  visual clustering artifacts appear in practice.
- Boundaries: SDF projection against drum/lifters with wall velocity blending
  `v ← (1−β)·v + β·v_wall` for particles that were projected (touching the wall) this substep (β = no-slip factor,
  default 1 → no-slip).
- **Viscosity (v1, implemented)**: XSPH `v_i += c·Σ (m_j/ρ_j)(v_j − v_i) W_ij` with a smooth saturating coefficient
  `c = sqrt(μ)/(sqrt(μ)+sqrt(μ_half))`, `μ_half = 15 Pa·s` (the viscosity at which `c = 0.5`; `c → 1` asymptotically,
  never exactly reached, for any `μ`). An earlier version used a hard-clamped `c = clamp(μ/2, 0, 1)`, which saturated
  at 2 Pa·s — every viscosity from 2 Pa·s upward produced identical output, silently capping the UI's usable range.
  The current curve is monotone across the full 0–200 Pa·s the UI exposes but is still qualitative (ordinal, not
  quantitatively calibrated to real Pa·s: XSPH's effective kinematic viscosity is bounded by `~h²/dt` regardless of
  `c`) — consistent with this project's project-wide 2D/qualitative caveat. See `pbf.rs`'s `xsph_coefficient` doc
  comment for the full derivation.
- **Viscosity (v2, phase M7)**: implicit viscosity (Weiler et al. 2018) with conjugate gradient; Bingham/Herschel–Bulkley via
  Papanastasiou-regularised effective viscosity `μ_eff = K·γ̇^(n−1) + τ_y(1−e^{−m γ̇})/γ̇`.
- Dye: each fluid particle carries `dye ∈ [0,1]` (initial: left half 0 / right half 1, or top/bottom). Pure Lagrangian tracer (no diffusion) → mixing index computed in `metrics.rs`.

### 3.4 Coupling (`coupling.rs`) — one shared sub-step for both solvers, implemented
Balls (ss3.2) and fluid (ss3.3) run on the *same* fixed sub-step (`FIXED_DT`, `simulation.substeps`
per frame) rather than a separate DEM/fluid time-step ratio; `coupling::step` orchestrates the
staggered exchange each sub-step (fluid solves first, using last sub-step's ball positions/
velocities as a fixed boundary; balls then advance using the resulting impulses, so their new state
becomes the fluid's boundary for the *next* sub-step). The exchange is expressed as **impulses**,
not forces (`CouplingImpulses`, applied as `Δv_ball = impulse * inv_mass` directly, no extra
`* dt`) -- a force-based version (impulse/dt, re-integrated with another `* dt` on the ball side)
was tried first and is more `dt`-sensitive than necessary; working in impulses throughout matches
how the rest of the solver already applies position/velocity corrections directly (ss3.2/3.3).
1. During `FluidParticles::step_coupled`, after its density-constraint solve, fluid particles inside a ball
   (dist < r_b + 0.5·dx) are projected to the ball surface by a position correction `push`, **capped at
   half the ball's radius** per sub-step (bounds the correction itself, not just its downstream impulse,
   for particles that start or persistently end up deep inside the contact zone). Since a position
   correction becomes a velocity via `Δv = push/dt` (ss3.3 step 5's reconstruction), the fluid's own
   momentum change (impulse) is `m_f · push / dt`; the ball's Newton's-third-law reaction impulse is the
   exact opposite, `-(m_f · push) / dt`.
2. Viscous no-slip at ball surface: fluid particles within `h` of a ball surface blend towards the ball's surface velocity
   `v_surf = v_b + ω_b × (x_i − x_b)` with factor `β_b·c_eff`; here `Δv` is already a velocity, so the
   ball's reaction impulse is simply `-(m_f · Δv)` (no `/dt` at all, unlike step 1's position-based `push`).
3. Both contributions' angular impulses use the lever arm from the ball's centre to the contact point/fluid
   particle. The accumulated per-ball impulse is clamped to `|impulse| ≤ 3·m_b·g·dt` (angular impulse
   clamped to that times `r_ball`) for stability, then applied as a direct velocity change in the ball
   solver's predict step (`DemState::step_with_external_forces`, ss3.2 step 1). The `3x`-gravity multiple
   (down from an initially-planned `20x`) was reached empirically: at this project's default scale, where
   the coarse-grained ball radius can end up *smaller* than the fluid's own particle spacing (e.g.
   defaults: ~6 mm balls vs. ~12.5 mm fluid spacing at `resolution = 40`), `20x` was still visibly ejecting
   balls from the charge over several seconds of simulated time; `3x`, combined with the push cap in step 1,
   stayed visually stable over 30+ seconds while still giving a meaningful buoyancy/drag effect.
   (Buoyancy emerges from the density-constraint push; steel in slurry is a minor effect so v1 does not add ball
   boundary particles to the density sum. v2 option: sample ball perimeter as boundary particles contributing to
   `ρ_i` — Akinci-style — for correct buoyancy.)
4. Balls advance their one sub-step; new positions/velocities become the moving boundary for the fluid's next sub-step.

### 3.5 Free surface & metrics (`surface.rs`, `metrics.rs`)
- Scalar field `φ` on a `G×G` grid (G = 128, spanning the drum bbox): splat each fluid particle with a smooth kernel of
  radius 1.5 h; marching squares at `φ = 0.5·φ_full` → closed polylines; Chaikin smoothing ×1. `φ_full` is computed
  *analytically*, not read off the grid's measured maximum: `φ_full = n·∫K dA`, with number density `n = rest_density /
  particle_mass` (exact, since `particle_mass = rest_density·dx²` at seeding) and `∫K dA = π·radius²/4` in closed form
  for the cubic-falloff splat kernel. Reading `φ_full` off the grid instead is unsafe -- at high rotation speed the
  slurry can centrifuge into a thin wall-hugging film whose field values are depressed everywhere, while a transient
  PBF-solver over-compacted pocket elsewhere can inflate the measured maximum well past the true bulk value, starving
  the film's threshold and making the free surface disappear. As a graceful fallback for degenerate cases where the
  analytic threshold yields no contour at all, extraction retries once with the grid's measured-peak-based threshold
  (the original approach) instead. Exported as `Vec<f32>` [n_polys, len_0, x,y,…]. Optional: the polygon is filled in
  the renderer, particles hidden.
  The field is masked by the drum wall SDF before marching squares (no bulk value survives outside the wall or inside a
  lifter bar), and after smoothing every contour point is projected back onto the wall/lifter surface along its normal
  if it still landed outside -- together these prevent the rendered surface from bulging through the wall or lifters.
- Slurry-level metrics: pool angular extent along the wall (fluid particles within `1.5 h` of the wall → min/max angle),
  free-surface mean line (fit to boundary particles not near the wall: angle and offset), pool depth at bottom.
- Charge metrics: toe/shoulder angles from angular histogram of balls near the wall; charge center of mass;
  power draw estimate `P = Σ (F_wall × v_wall)` (torque on drum × ω) — smoothed.
- Mixing index (Lacey): dye variance over occupied grid cells normalised by initial variance, 0 → 1.
- Energy checks (debug): total KE, DEM overlap max, fluid density error max.

---

## 4. Frontend specification

### 4.1 Runtime architecture
- `worker.ts` instantiates the wasm module and `Simulation`. Loop: on each `requestFrame` message (sent by main
  thread rAF) the worker runs fixed 1/240 s substeps until `sim_time` catches up with wall time × `time_scale`, but
  never more than `budget_ms` (default 12 ms) per frame. It then posts `Frame { ballsXYRT: Float32Array, fluidXYD: Float32Array,
  surface: Float32Array, drumAngle, simTime, stats }` with **transferable** buffers (double-buffered). Achieved
  time scale is reported so the HUD shows e.g. "0.6× real time" instead of freezing → "near real time" by construction.
- SharedArrayBuffer + wasm threads are *not* required (Vite dev server would need COOP/COEP); listed as future work.
- Params: `SetParams { params, reset: boolean }`. Hot-swappable (no reset): rpm, viscosity, time scale, substeps/iters,
  no-slip factors, display options. Everything else (geometry, fill, ball sizes, resolution, seed, lifters) triggers `reset`.

### 4.2 Rendering (`render/canvas.ts`, Canvas 2D, DPR-aware)
Layers per frame: background → drum interior disc → lifters (rotated by drumAngle) → slurry surface polygon (fill, alpha 0.55)
→ fluid particles (optional; colour = dye, blue→orange) → balls (fill grey, spin indicator line) → drum wall ring +
rotation marker → overlays (toe/shoulder rays, free-surface line) → HUD text. 2000 balls + 8000 dots at 60 fps is fine on
Canvas 2D (dots drawn as `fillRect`; balls as `arc`). WebGL2 instancing is a documented fallback if profiling shows
render > 6 ms.

### 4.3 Parameters modal (`ui/paramsModal.ts`) — native `<dialog showModal()>`
Tabs (buttons toggling sections): **Mill**, **Media**, **Slurry**, **Lifters**, **Simulation**, **Display**.
Each field generated from `schema.ts` entries `{ id, group, label, unit, type: number|select|boolean, min, max, step,
default, hot, help }` with live validation (red outline + message, Apply disabled). Footer: `Presets ▾`, `Reset to defaults`,
`Cancel`, `Apply` (Apply asks "Restart simulation?" only if a non-hot param changed). Esc / backdrop click = Cancel.
Derived read-only values shown in the modal: critical speed, %Nc, true vs. simulated ball count,
coarse-graining factor `k` and effective media diameter (see ss3.2; `k = 1.0`/"no coarse-graining"
when the true count is already <= `simulation.max_balls`), fluid particle count, shared sub-step
rate (`FIXED_DT` x `substeps`), estimated cost.

| Group | Parameters (defaults) |
|---|---|
| Mill | diameter D = 1.0 m; speed mode = rpm (or %Nc); speed = 30 rpm (~70 %Nc for D = 1 m); rotation direction CCW |
| Media | ball diameter 10 mm (+ optional distribution rows); ball fill J = 0.30 (fraction of drum area incl. voids, packing 0.6 in 2D); density 6000 (ZrO2/YSZ ceramic); restitution 0.7 (ball–ball) / 0.5 (ball–wall); friction μ 0.25 / 0.35; rolling μ_r 0.01 |
| Slurry | enabled = true; fill U_s = 0.15 of drum area; density 1800 kg/m³; viscosity 0.5 Pa·s; rheology = Newtonian (Bingham in M7); wall no-slip β = 1.0; ball no-slip β_b = 1.0; dye pattern = left/right |
| Lifters | count = **0** (default, smooth wall); height 20 mm; base width 30 mm; top width 20 mm; phase 0° |
| Simulation | resolution 40 (particles across R); substeps 4; PBF iterations 3; DEM (XPBD) iterations 4; max balls 2000; time scale 1.0; frame budget 12 ms; seed 1 |
| Display | show fluid particles / surface polygon / dye / toe-shoulder / free-surface line / spin marker / velocity vectors; ball colour by speed |

### 4.4 Toolbar & HUD
Play/Pause (Space), Step (one frame), Reset, Parameters (opens modal), Presets menu, Screenshot (PNG via `canvas.toBlob`).
HUD: sim time, rpm and %Nc, balls/fluid counts, fps, sim ms/frame, achieved time scale, toe/shoulder, slurry pool angles, mixing index.

---

## 5. Milestones (each ends with a runnable, verified state)

Model assignment per CLAUDE.md: **[F]** = Fable (main session), **[S]** = sonnet subagent, **[H]** = haiku subagent.

### M0 — Bootstrap (½ day) [H/S]
1. Write `CLAUDE.md` (§2), `.gitignore`, `README.md`; `git init -b main`, add remote
   `https://github.com/ToshihiroIguchi/MillDynamics3` (repo exists, public, currently empty — verified with `gh`).
   Commit at the end of every milestone; push to `main` (user confirms first push).
2. Install Rust: `winget install Rustlang.Rustup` (or rustup-init.exe), `rustup target add wasm32-unknown-unknown`,
   `cargo install wasm-pack` (or `wasm-bindgen-cli`). Verify host toolchain links (MSVC from the installed VS 18, fall back to
   `stable-x86_64-pc-windows-gnu` with rtools gcc if MSVC C++ tools are absent). `cargo test` on a hello crate must pass.
3. Scaffold workspace (§1), `mill-core` with `Params` + empty `Simulation`, `mill-wasm` exposing `new/step/ptrs`, Vite app
   with worker + canvas drawing a rotating empty drum from wasm-provided `drumAngle`.
   Verify: `scripts/build-wasm` + `npm run dev` → rotating drum with rotation marker at correct rpm.

### M1 — Balls (position-based rigid discs) in rotating drum (1–2 days) [F core, S tests/UI]
1. [F] `grid.rs`, `geometry.rs` (circle-only SDF + wall velocity), `dem.rs` XPBD contact solver, ball init lattice.
2. [S] Unit tests: single ball drop restitution (bounce height ratio ≈ e²), two-ball head-on momentum/energy, ball resting on wall
   (overlap < 1 % r, no drift), rotating drum centrifuging test (at 120 % Nc all balls stay within 1.5 r of wall after 5 s),
   at 70 % Nc a cascading charge with toe ~ 200–230° and shoulder ~ 40–60° (measured from the vertical, loose tolerance).
3. [S] Renderer for balls; params modal v1 (Mill + Media + Simulation groups, schema-driven); HUD basics; Play/Pause/Reset.
4. [S] criterion bench: 500/1000/2000 balls step cost; record in docs.
   Verify: visually cascading at 70 %, cataracting at 85–90 %, centrifuging > 100 %.

### M2 — Lifters (½ day) [S, F review of SDF]
1. Lifter quads in `geometry.rs` (drum frame), numeric gradient, restricted evaluation region; `Lifters` modal tab; render.
2. Tests: SDF sign/gradient sanity at sampled points; ball resting on a lifter face is stable; default `count = 0` produces the
   identical trajectory to M1 (regression: hash of positions after 2 s equals circle-only run).

### M3 — Slurry PBF alone (2 days) [F core, S surface/metrics/UI]
1. [F] `pbf.rs`: kernels, neighbour grid, density constraint, s_corr, boundary projection with wall velocity, XSPH viscosity, dye.
2. [S] `surface.rs` marching squares + export; `metrics.rs` pool angles / free-surface line; renderer surface polygon + dye colouring;
   Slurry + Display tabs.
3. Tests: hydrostatic column at rest (max density error < 2 % after 1 s, no drift); volume conservation (particle count constant,
   bulk area within 5 %); rotating drum with viscous fluid — steady-state free-surface tilt increases monotonically with viscosity
   (3 viscosities); marching-squares unit test on a synthetic disc field (area within 3 %).
4. [F] Calibrate `μ → c_eff` mapping; document in PHYSICS.md.

### M4 — Two-way coupling (2 days) [F]
1. `coupling.rs` as §3.4; substep orchestration in `lib.rs` (`step_frame(wall_dt)` with accumulator and budget).
2. Tests: ball dropped into a pool decelerates (terminal velocity lower than dry drop); fluid particle count inside balls == 0 after
   each substep; momentum exchange symmetric (Σ impulses fluid = −Σ impulses balls within 1e-3 relative); 30 s run at
   1000 balls + 4000 fluid with no NaN, max overlap < 3 % r, density error < 5 %.
3. Visual: slurry is dragged up by the charge, pool forms at the toe, dye mixes.

### M5 — UI completion (1 day) [S/H]
Presets, full validation UX, derived values panel, Step button, Screenshot, mixing-index sparkline, keyboard shortcuts,
responsive layout (canvas fits viewport, modal scrolls on small screens), `docs/PARAMETERS.md` generated from `schema.ts` [H].
Playwright e2e: open modal → change rpm → Apply → HUD rpm updates without reset; change lifter count → Apply → reset confirm
→ lifters visible; invalid input blocks Apply.

### M6 — Performance (1–2 days) [F profiling/SIMD, S tooling]
Targets on a mid laptop, single thread: **≥ 1.0× real time at 500 balls + 2000 fluid; ≥ 0.5× at 1000 + 4000**.
Actions in order until targets met: `-C target-feature=+simd128` + `opt-level=3` + `lto`; SoA layouts and f32; neighbour-list reuse
across PBF iterations; DEM ball–ball broadphase only for balls (fluid grid separate, cell = h); avoid allocation per step;
`wasm-opt -O3`. Add an "auto resolution" option that lowers `resolution` when achieved time scale < 0.5 for 3 s.
Deliver a perf table in README.

### M7 — Fidelity extensions (optional, after user review) [F]
Implicit viscosity + Bingham/Herschel–Bulkley; Akinci boundary particles on balls (proper buoyancy); WCSPH "accurate" mode
(non-real-time); size distribution for media; WebGL2 renderer; SharedArrayBuffer/wasm-threads; data export (CSV of metrics).

---

## 6. Verification (end-to-end)

- `cargo fmt --check && cargo clippy --all-targets -- -D warnings && cargo test -p mill-core` — all physics tests in §5.
- `cargo bench -p mill-core` — step cost table (balls × fluid).
- `scripts/build-wasm` → `web/src/wasm/` produced; `npm run build` succeeds; `npm run test` (vitest) and `npm run e2e` (playwright, Chromium)
  pass.
- Manual/visual acceptance (use the `run` skill / Playwright screenshots):
  1. Default params (no lifters, 70 % Nc, 0.5 Pa·s): cascading charge, slurry pool at toe, free-surface polygon drawn, ≥ 0.8× real time.
  2. Set lifters = 8: balls lifted higher, visible cataracting.
  3. Viscosity 5 Pa·s vs 0.05 Pa·s: free surface tilts more with higher viscosity; dye mixes slower.
  4. Speed 110 % Nc: centrifuging; slurry film on wall.
  5. Modal: every field validated; hot params apply without reset; others prompt for reset.
- Determinism: same params + seed → identical position hash after 5 s (test in `mill-core`).

## 7. Risks and mitigations

| Risk | Mitigation |
|---|---|
| PBF viscosity only qualitative | Calibrated mapping in v1 documented as such; implicit viscosity + Bingham in M7 |
| Coupling instability (light fluid pushing heavy balls is fine, but many overlapping corrections) | Impulse clamp, shared fixed sub-step for both solvers (ss3.2/3.4, no separate DEM/fluid time-step ratio to get wrong), symmetric momentum test, per-substep overlap resolution |
| WASM single-thread too slow at high resolution | Frame budget + achieved-time-scale HUD, auto resolution, SIMD; targets defined in M6 |
| 2D slice vs real 3D mill quantitatively different | State clearly in README/HUD ("2D cross-section, qualitative") |
| Rust toolchain on Windows (MSVC linker) | Check in M0; fall back to GNU toolchain with rtools gcc; wasm target needs no linker |
| Lifter SDF gradient numeric noise at corners | Round corners with small radius in SDF; test resting stability |
| Default media (2 mm) at default mill diameter (1 m) with J = 0.30 implies ~75,000 balls in 2D — far above the M6 real-time budget (500–2000 balls) | **Resolved, implemented** (`Params::effective_media`, ss3.2): coarse-graining is applied automatically whenever the true ball count exceeds `simulation.max_balls` (default 2000), substituting fewer, larger, lighter balls that preserve total media mass and fill fraction. The effective diameter, coarse-graining factor, and simulated vs. true ball count are always shown explicitly in the parameters modal's derived-values panel (ss4.3) — never applied silently. |
