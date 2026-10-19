# Simulation Parameters Reference

This reference documents every user-facing simulation parameter. The table data is generated from `web/src/params/schema.ts` (UI-facing ranges, units, and display scaling) and `crates/mill-core/src/params.rs` (SI defaults, validation ranges, and physical descriptions). SI units are used internally throughout `mill-core`; the UI converts via the `displayScale` factor noted below. Regenerate this table after changing parameter defaults or ranges in either source file.

## Panel layout: Advanced sections and hidden fields

In the parameters panel, rarely-touched fields (`advanced: true` in `schema.ts`) are folded into a collapsed **Advanced** sub-section per group: Media (packing fraction, restitution, friction, rolling friction), Lifters (base width, top width, phase), Slurry (wall no-slip, wettability, surface tension, dye pattern) and Numerical accuracy (sub-steps, ball solver iterations, seed). Fields with no effect under the current settings are hidden instead of greyed out: lifter geometry when Lifter count is 0; Max balls when Coarse-graining is Off (it stays visible but greyed out in Manual); Manual k unless Coarse-graining is Manual; all other Slurry rows when Slurry is disabled. A hidden field still keeps and sends its value, and a field that fails validation is always revealed (with its group and Advanced section opened). This is UI-only: `Params`, presets and reproducibility are unchanged.

## Mill Parameters

| Path | Label | Unit (UI) | Default (SI) | Min–Max (UI) | Step (UI) | Notes |
|------|-------|-----------|--------------|--------------|-----------|-------|
| `mill.diameter_m` | Drum diameter | mm | 0.063 (63 mm) | 30–5000 | 1 | Inner drum diameter; the low end covers small bench/lab jar mills (~3 cm). Validation range is (0.02, 20.0) m SI. |
| `mill.speed_mode` | Speed mode | — | `rpm` | `percent_critical` / `rpm` | — | Whether speed is specified as absolute rpm or percentage of critical speed. |
| `mill.speed_value` | Speed | — | 100.0 | 0–300 | 1 | Rotation speed: interpretation depends on `speed_mode`. Default (100 rpm, ~59% Nc for the default 63 mm drum, Nc = 168.5 rpm). |
| `mill.direction` | Direction | — | `counter_clockwise` | `counter_clockwise` / `clockwise` | — | Drum rotation direction viewed from the standard right-handed 2D frame (+x right, +y up). |

## Media Parameters

| Path | Label | Unit (UI) | Default (SI) | Min–Max (UI) | Step (UI) | Notes |
|------|-------|-----------|--------------|--------------|-----------|-------|
| `media.ball_diameter_m` | Ball diameter | mm | 0.002 (2 mm) | 0.5–200 | 0.1 | Diameter of grinding balls. Defaults model yttria-stabilized zirconia (ZrO₂) ceramic media in a lab-scale 63 mm drum; `True ball count` (~244) stays under every quality preset's `Max balls` (no coarse-graining out of the box, `k = 1.0`) at the default mill / fill / packing fraction. Validation range (0.0005, 1.0) m SI. |
| `media.fill_fraction` | Fill fraction (J) | — | 0.30 | 0–0.9 | 0.01 | Fraction of drum cross-sectional area filled with media including voids (conventional "J" ball-loading fraction). |
| `media.packing_fraction_2d` | 2D packing fraction | — | 0.82 | 0.5–0.907 | 0.001 | Areal packing fraction of settled disc charge (solid area / footprint area). Random close packing ≈ 0.82; hexagonal upper bound ≈ 0.907. |
| `media.density_kg_m3` | Media density | kg/m³ | 6000.0 | 100–20000 | 100 | Ball material density. Default (6000 kg/m³) models ZrO₂ ceramic. Balls modelled as unit-depth discs. A "Media material" dropdown in the UI (`web/src/params/materials.ts`) offers representative densities for common grinding media as a shortcut into this field: YSZ (6000, the default), stainless steel (7700), alumina (3600), mullite (3000); an explicit "(custom)" option keeps the raw number editable. |
| `media.restitution_ball_ball` | Restitution (ball-ball) | — | 0.7 | 0–1 | 0.01 | Coefficient of restitution for ball–ball collisions. Higher (≈0.7–1) for ceramic; lower (≈0.3–0.5) for steel. |
| `media.restitution_ball_wall` | Restitution (ball-wall) | — | 0.5 | 0–1 | 0.01 | Coefficient of restitution for ball–wall collisions. |
| `media.friction_ball_ball` | Friction (ball-ball) | — | 0.25 | 0–2 | 0.01 | Coefficient of kinetic friction for ball–ball contact. Validation range [0, ∞). |
| `media.friction_ball_wall` | Friction (ball-wall) | — | 0.35 | 0–2 | 0.01 | Coefficient of kinetic friction for ball–wall contact. |
| `media.rolling_friction` | Rolling friction | — | 0.01 | 0–1 | 0.001 | Rolling-resistance coefficient (dimensionless torque coefficient). |

Not yet UI-exposed (no `schema.ts` entry): `media.friction_ball_ball_static`/`media.friction_ball_wall_static` (static/stick Coulomb friction coefficients, each defaulting to its kinetic counterpart above -- no behavior change until set higher) and `media.friction_velocity_scale_m_s` (default 0.02 m/s, the relative-speed scale over which friction blends from the static to the kinetic coefficient). See `dem.rs`'s `effective_friction` and docs/PHYSICS.md §9's 2026-09-28 entry.

## Lifters Parameters

| Path | Label | Unit (UI) | Default (SI) | Min–Max (UI) | Step (UI) | Notes |
|------|-------|-----------|--------------|--------------|-----------|-------|
| `lifters.count` | Lifter count | — | 0 | 0–64 | 1 | Number of evenly-spaced lifter bars on the drum wall. Zero (the default) means a perfectly smooth wall with no lifters. |
| `lifters.height_m` | Height | m | 0.020 | 0–0.2 | 0.005 | Lifter blade height. Only enforced if `count > 0`. |
| `lifters.base_width_m` | Base width | m | 0.030 | 0–0.3 | 0.005 | Base width of lifter blade. Only enforced if `count > 0`. |
| `lifters.top_width_m` | Top width | m | 0.020 | 0–0.3 | 0.005 | Top width of lifter blade. Only enforced if `count > 0`. |
| `lifters.phase_deg` | Phase offset | deg | 0.0 | −180–180 | 1 | Angular offset of the first lifter from the +x axis (degrees). |

## Slurry Parameters

| Path | Label | Unit (UI) | Default (SI) | Min–Max (UI) | Step (UI) | Notes |
|------|-------|-----------|--------------|--------------|-----------|-------|
| `slurry.enabled` | Enabled | — | true | on / off | — | Enable or disable slurry coupling. |
| `slurry.fill_fraction` | Fill fraction | — | 0.35 | 0–0.9 | 0.01 | Fraction of drum cross-sectional area filled with liquid slurry. Kept above `media.fill_fraction`'s default so the liquid surface clears the settled ball heap. |
| `slurry.density_kg_m3` | Slurry density | kg/m³ | 1000.0 | 100–5000 | 50 | Slurry fluid density. |
| `slurry.viscosity_pa_s` | Viscosity | Pa·s | 0.001 | 0–200 | 0.1 | Dynamic viscosity of slurry (default: water; a thick mineral-processing slurry is tens of Pa·s). Accuracy is lower at water-like viscosity: the PBF result is not converged in the cost knobs (docs/VERIFICATION.md). Bingham yield-stress rheology is implemented in the solver (`bingham_viscosity` in `pbf.rs`) but is not exposed in the UI, so the app always runs Newtonian. `step` is a UI spinner granularity hint only, not an enforced constraint — the parameters panel validates against `min`/`max` in its own JS, not via native HTML step-grid validation. |
| `slurry.wall_no_slip` | Wall no-slip (β) | — | 1.0 | 0–1 | 0.05 | No-slip blend factor at drum wall and lifters, in [0, 1] (1 = full no-slip boundary condition). |
| `slurry.wettability` | Wettability | — | 0.6 | 0–1 | 0.05 | Ball<->slurry wettability, in [0, 1]. `0` = non-wetting (fluid is never pulled toward a ball's surface, only pushed off it — a 180° contact angle). `1` = strongly wetting. Drives an Akinci-style boundary-adhesion force, range-capped at `2×` the ball's own radius regardless of fluid resolution. Together with `surface_tension_n_m` this gives a genuine contact angle (wetting pulls fluid onto media, cohesion holds the resulting film together) rather than a single knob trying to do both jobs. See `docs/PHYSICS.md` §6.1a. |
| `slurry.surface_tension_n_m` | Surface tension | N/m | 0.072 | 0–0.2 | 0.005 | Fluid-fluid surface tension, driving an Akinci-style pairwise cohesion force; `0` disables it. `0.072` is real water's value, used here as a calibrated reference point (like the wider SPH surface-tension literature, this is a proportionality to the physical value, not a first-principles unit conversion — see `docs/PHYSICS.md` §6.1a). `#[serde(default = ...)]` in params.rs deserializes an older client's JSON (predating this field) as `0.072`, not `f32`'s bare `0.0`, so wetting stays physically coherent regardless of client version. |
| `slurry.dye_pattern` | Dye pattern | — | `left_right` | `left_right` / `top_bottom` / `none` | — | Initial tracer dye pattern for visualizing and measuring mixing. |

**Note:** Fields `slurry.rheology` and `slurry.yield_stress_pa` exist in the data model (`crates/mill-core/src/params.rs`) but are not exposed in the UI. Bingham rheology is implemented in the solver (`bingham_viscosity` in `crates/mill-core/src/pbf.rs`), so setting `rheology`/`yield_stress_pa` through the JSON API takes effect; the UI always sends the Newtonian default. The viscosity solve is implicit (conjugate gradient) with a central-force pairwise form (see `docs/PHYSICS.md`).

**Derived: Interstitial filling (U).** The params panel's Derived section also shows `U = slurry.fill_fraction / (media.fill_fraction * (1 - media.packing_fraction_2d))`, the ratio of slurry volume to the ball charge's void volume; see `docs/PHYSICS.md` §9's "2D areal packing is not 3D voidage" entry for why this project's 2D packing makes the same `slurry.fill_fraction` correspond to a much higher `U` (and deeper-looking pool) than an equivalent real 3D mill.

**Note: limited percolation into a settled bed.** The Slurry group in the params panel shows a one-line note (the longer explanation is in the Slurry resolution field's "?" tooltip) that slurry does not meaningfully seep into a stationary, settled ball bed, even at low viscosity with the drum stopped — expect most slurry to stay visible as a pool above the bed. This is an expected consequence of the fluid lattice spacing being comparable to the gaps between packed balls, compounded by the lower, less-connected void fraction of a 2D packed-circle bed versus a real 3D mill; see `docs/PHYSICS.md` §9's "Slurry does not seep into a settled/stationary ball bed" entry for the measurements.

## Numerical accuracy Parameters

The panel group formerly called "Simulation" is now "Numerical accuracy"; the underlying Params paths remain `simulation.*`.

| Path | Label | Unit (UI) | Default (SI) | Min–Max (UI) | Step (UI) | Notes |
|------|-------|-----------|--------------|--------------|-----------|-------|
| `simulation.substeps` | Sub-steps / frame | — | 8 | 1–64 | 1 | Requested fixed sub-steps per rendered frame at 1× time scale; the solver step is `dt = 1 / (60 * effective_substeps)`. `Params::effective_substeps` raises the count actually used above this request (auto-raise capped at 16; an explicit request may go up to 64) when the ball size relative to the expected peak ball speed would otherwise push the XPBD stability criterion `max_substep_displacement_over_diameter` too close to 1. At the shipped defaults (63 mm drum, 2 mm balls, 100 rpm) the requested 8 is raised to 16, i.e. `dt = 1/960 s`. See `docs/PHYSICS.md` §9 for the rationale and measured numbers. |
| `simulation.dem_iterations` | Ball solver iterations | — | 2 | 1–20 | 1 | XPBD non-penetration solver iterations per sub-step for the ball population. Trades contact-resolution accuracy for cost rather than stability. Halved (from 4) alongside doubling `substeps` (from 4) to keep DEM cost per rendered frame roughly flat while approximately halving the tunnelling-risk ratio `max_substep_displacement_over_diameter`; see `docs/PHYSICS.md` ss9. |
| `simulation.max_balls` | Max balls (Auto coarse-graining only) | — | 600 (mill-core `Default`); **300 at a fresh app load** (Fast preset) | 10–50000 | 10 | Used **only when Coarse-graining is Auto**: target upper bound on simulated ball particles; if the true media population exceeds it, balls are enlarged by `k`. Greyed out in the UI when the mode is Off (the default) or Manual. See `docs/PLAN.md` ss3.2 and the Quality presets section below. |
| `simulation.coarse_graining_mode` | Coarse-graining | — | **`off`** (true-size balls) | `off` / `auto` / `manual` | — | Off simulates true-size balls and is the default. Auto (target Max balls) and Manual (fixed `simulation.coarse_graining_k`) enlarge the balls and **change the physics**: kinetic energy, per-impact stress statistics and the shoulder angle are not preserved, and mill power draw is off by 17-39 % at large D/d (`docs/VERIFICATION.md`, *What coarse-graining preserves*). Use it only on explicit choice, for qualitative charge-shape checks. With Off the UI warns above 3000 true balls (much slower than real time) and blocks Apply above 30000. The HUD shows a `COARSE-GRAINED k=…` tag while the live k > 1. |
| `simulation.resolution` | Slurry resolution (particles/radius) | — | 40 (mill-core `Default`); **25 at a fresh app load** (Fast preset) | 4–200 | 1 | Fluid particles spanning the drum radius; sets the PBF lattice spacing `dx = drum_radius_m / resolution` and kernel radius `h = 2*dx`, baked in at seed time (`FluidParticles::seed_lattice`) -- changing it via `Simulation::set_params` alone has no effect on an already-seeded population, a full reset is required (`web/src/params/schema.ts`'s `resetRequired: true`). The resolution actually used to seed the fluid can be raised above this value by `Params::effective_fluid_resolution` when the balls (e.g. coarse-grained ones) would otherwise be narrower than the fluid lattice spacing; the raise is surfaced in the UI, not silent -- see `docs/PHYSICS.md` ss9 for the full rationale. |
| `simulation.time_scale` | Time scale | — | 1.0 | 0.1–5 | 0.1 | Wall-clock-to-simulation-time multiplier requested by the user. Achieved rate is reported back and may be lower depending on frame budget. |
| `simulation.seed` | Random seed | — | 1 | 0–1,000,000,000 | 1 | Random number seed for reproducible initial ball configuration. All runs with identical `Params` and `seed` are deterministic. |

### Quality presets

`web/src/params/presets.ts` offers four presets (three pre-measured `(max_balls, resolution)` pairs plus "Impact statistics", below) as one-click
choices in the parameters panel's Numerical accuracy group ("Quality preset" selector), still requiring
Apply (both fields reset the simulation) like any other field edit. Measured in-browser, WASM
release build, default drum/media/slurry otherwise, `lifters.count = 0`, steady cascading after
~10 s of sim time -- see `docs/PERF.md` for the full methodology. "Fluid spacing / ball diameter"
is this same panel's live derived-values readout (`h / (2 * effective ball diameter)`, a fixed
multiple of the `h/r` wetting-resolution ratio `docs/PHYSICS.md` discusses); every tier keeps it
well above the point where a thin wetting film is visually resolvable -- that trade-off is real
and disclosed here, not hidden, and reducing it needs much more solver speed (the M6 and M7 performance passes have landed and
were not enough) or much more aggressive ball coarse-graining than any of these tiers use.

The fourth preset, "Impact statistics (slow)", additionally sets 32 sub-steps, 12 ball-solver iterations and coarse-graining off (`max_balls` 600, `resolution` 25): it exists because the Grinding proxy (3D) impact statistics were measured to be only reasonably converged (a few percent) at that setting, whereas the app defaults (8 sub-steps, 2 iterations) are not; it runs several times slower than real time, impact counts still drift roughly +5 % per doubling of sub-steps, and runs should be compared only with identical settings.

**Stale `Achieved speed` column.** These figures were measured before `media.ball_diameter_m`'s
default changed to a 2 mm ball in a 63 mm drum (chosen so the shipped config needs no ball coarse-graining
-- see the Media table above). At the old default, each tier coarse-grained a
different simulated ball count; at the new default `N_true` (~244) sits under all three tiers'
`max_balls` (300/300/600), so every tier now simulates the *same* ~244 real 2 mm balls -- only `resolution` (fluid particle count) still differs between tiers. `Fluid
spacing / ball diameter` below is recomputed for the new default (`dx / d_eff` with `d_eff` fixed
at 2 mm across all three tiers; the solver auto-raises the fluid resolution when this exceeds the
warning threshold, `Params::effective_fluid_resolution`); `Achieved speed` has not been re-measured against the new,
much smaller ball population and should be treated as directional only until `docs/PERF.md` is
refreshed.

| Preset | `max_balls` | `resolution` | Achieved speed (stale, see above) | Fluid spacing / ball diameter |
|---|---|---|---|---|
| **Fast** (default at a fresh app load; id `realtime`) | 300 | 25 | not re-measured (native 67 ms/frame) | ~0.63 |
| Balanced | 300 | 30 | not re-measured (native 72 ms/frame) | ~0.5 |
| Accuracy | 600 | 50 | not re-measured (native 202 ms/frame) | ~0.3 |

(Pre-2026-10-08 note, for the old Realtime at resolution 15:) the ratio (~1.05) sat just over the panel's own `> 1.0` "a single fluid particle is wider
than a ball" warning threshold at the new 2 mm default (the solver auto-raises its effective
resolution in that case), Balanced and Accuracy stay under it, and all three are still coarse relative to what a visually smooth *thin film* specifically
needs (a film a few millimetres thick is not resolvable by particles a third-to-half a ball diameter
across) -- the warning threshold and "a thin wetting film looks smooth" are different bars, and none
of these three tiers clears the second one.

**Note:** Fields `simulation.pbf_iterations` (default 3, controls the PBF density-constraint solver) and `simulation.frame_budget_ms` (default 12.0 ms) exist in the data model but are not exposed in the UI. `simulation.resolution` is exposed (as "Slurry resolution", above).
