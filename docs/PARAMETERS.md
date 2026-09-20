# Simulation Parameters Reference

This reference documents every user-facing simulation parameter. The table data is generated from `web/src/params/schema.ts` (UI-facing ranges, units, and display scaling) and `crates/mill-core/src/params.rs` (SI defaults, validation ranges, and physical descriptions). SI units are used internally throughout `mill-core`; the UI converts via the `displayScale` factor noted below. Regenerate this table after changing parameter defaults or ranges in either source file.

## Mill Parameters

| Path | Label | Unit (UI) | Default (SI) | Min–Max (UI) | Step (UI) | Notes |
|------|-------|-----------|--------------|--------------|-----------|-------|
| `mill.diameter_m` | Drum diameter | m | 1.0 | 0.1–5 | 0.05 | Inner drum diameter. Validation range is (0.05, 20.0) m SI. |
| `mill.speed_mode` | Speed mode | — | `rpm` | `percent_critical` / `rpm` | — | Whether speed is specified as absolute rpm or percentage of critical speed. |
| `mill.speed_value` | Speed | — | 30.0 | 0–300 | 1 | Rotation speed: interpretation depends on `speed_mode`. Default (~70% Nc for the default 1.0 m drum). |
| `mill.direction` | Direction | — | `counter_clockwise` | `counter_clockwise` / `clockwise` | — | Drum rotation direction viewed from the standard right-handed 2D frame (+x right, +y up). |

## Media Parameters

| Path | Label | Unit (UI) | Default (SI) | Min–Max (UI) | Step (UI) | Notes |
|------|-------|-----------|--------------|--------------|-----------|-------|
| `media.ball_diameter_m` | Ball diameter | mm | 0.010 (10 mm) | 0.5–200 | 0.1 | Diameter of grinding balls. Defaults model yttria-stabilized zirconia (ZrO₂) ceramic media; 10 mm is typical for lab/bench tumbling mills. Validation range (0.0005, 1.0) m SI. |
| `media.fill_fraction` | Fill fraction (J) | — | 0.30 | 0–0.9 | 0.01 | Fraction of drum cross-sectional area filled with media including voids (conventional "J" ball-loading fraction). |
| `media.packing_fraction_2d` | 2D packing fraction | — | 0.82 | 0.5–0.907 | 0.001 | Areal packing fraction of settled disc charge (solid area / footprint area). Random close packing ≈ 0.82; hexagonal upper bound ≈ 0.907. |
| `media.density_kg_m3` | Media density | kg/m³ | 6000.0 | 100–20000 | 100 | Ball material density. Default (6000 kg/m³) models ZrO₂ ceramic. Balls modelled as unit-depth discs. |
| `media.restitution_ball_ball` | Restitution (ball-ball) | — | 0.7 | 0–1 | 0.01 | Coefficient of restitution for ball–ball collisions. Higher (≈0.7–1) for ceramic; lower (≈0.3–0.5) for steel. |
| `media.restitution_ball_wall` | Restitution (ball-wall) | — | 0.5 | 0–1 | 0.01 | Coefficient of restitution for ball–wall collisions. |
| `media.friction_ball_ball` | Friction (ball-ball) | — | 0.25 | 0–2 | 0.01 | Coefficient of kinetic friction for ball–ball contact. Validation range [0, ∞). |
| `media.friction_ball_wall` | Friction (ball-wall) | — | 0.35 | 0–2 | 0.01 | Coefficient of kinetic friction for ball–wall contact. |
| `media.rolling_friction` | Rolling friction | — | 0.01 | 0–1 | 0.001 | Rolling-resistance coefficient (dimensionless torque coefficient). |

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
| `slurry.fill_fraction` | Fill fraction | — | 0.15 | 0–0.9 | 0.01 | Fraction of drum cross-sectional area filled with liquid slurry. |
| `slurry.density_kg_m3` | Slurry density | kg/m³ | 1800.0 | 100–5000 | 50 | Slurry fluid density. |
| `slurry.viscosity_pa_s` | Viscosity | Pa·s | 50 | 0–200 | 0.1 | Dynamic viscosity of slurry (a real mineral-processing slurry, not near-water). Only Newtonian rheology is implemented (M1–M6); Bingham yield stress (M7) is reserved for future work. `step` is a UI spinner granularity hint only, not an enforced constraint — the parameters panel validates against `min`/`max` in its own JS, not via native HTML step-grid validation. |
| `slurry.wall_no_slip` | Wall no-slip (β) | — | 1.0 | 0–1 | 0.05 | No-slip blend factor at drum wall and lifters, in [0, 1] (1 = full no-slip boundary condition). |
| `slurry.wettability` | Wettability | — | 0.6 | 0–1 | 0.05 | Ball<->slurry wettability, in [0, 1]. `0` = non-wetting (fluid is never pulled toward a ball's surface, only pushed off it — a 180° contact angle). `1` = strongly wetting. Drives an Akinci-style boundary-adhesion force, range-capped at `2×` the ball's own radius regardless of fluid resolution. Together with `surface_tension_n_m` this gives a genuine contact angle (wetting pulls fluid onto media, cohesion holds the resulting film together) rather than a single knob trying to do both jobs. See `docs/PHYSICS.md` §6.1a. |
| `slurry.surface_tension_n_m` | Surface tension | N/m | 0.072 | 0–0.2 (suggested) | 0.005 | Fluid-fluid surface tension, driving an Akinci-style pairwise cohesion force; `0` disables it. `0.072` is real water's value, used here as a calibrated reference point (like the wider SPH surface-tension literature, this is a proportionality to the physical value, not a first-principles unit conversion — see `docs/PHYSICS.md` §6.1a). Not yet exposed in the UI (`web/src/params/schema.ts`); deserializes from an older client's JSON as `0.072` (not `0.0`), not this crate's own `f32` default, so wetting stays physically coherent until the UI control lands. |
| `slurry.dye_pattern` | Dye pattern | — | `left_right` | `left_right` / `top_bottom` / `none` | — | Initial tracer dye pattern for visualizing and measuring mixing. |

**Note:** Fields `slurry.rheology` and `slurry.yield_stress_pa` exist in the data model (`crates/mill-core/src/params.rs`) but are not yet exposed in the UI because Bingham rheology is unimplemented (scheduled for M7 per `docs/PLAN.md`). The solver currently always runs an implicit (conjugate-gradient) Newtonian viscosity solve (Morris SPH Laplacian, see `docs/PHYSICS.md`) regardless of the `rheology` field.

## Simulation Parameters

| Path | Label | Unit (UI) | Default (SI) | Min–Max (UI) | Step (UI) | Notes |
|------|-------|-----------|--------------|--------------|-----------|-------|
| `simulation.substeps` | Sub-steps / frame | — | 8 | 1–16 | 1 | Fixed sub-steps per rendered frame at 1× time scale. Nominal PBF step rate is `substeps × 60` Hz before frame-budget throttling. |
| `simulation.dem_iterations` | Ball solver iterations | — | 2 | 1–20 | 1 | XPBD non-penetration solver iterations per sub-step for the ball population. Trades contact-resolution accuracy for cost rather than stability. Halved (from 4) alongside doubling `substeps` (from 4) to keep DEM cost per rendered frame roughly flat while approximately halving the tunnelling-risk ratio `max_substep_displacement_over_diameter`; see `docs/PHYSICS.md` ss9. |
| `simulation.max_balls` | Max balls (coarse-graining target) | — | 600 | 10–50000 | 10 | Target upper bound on simulated ball particles. If the true media population exceeds this, coarse-graining applies transparent particle scaling to preserve total charge mass and footprint area. See `docs/PLAN.md` ss3.2. |
| `simulation.time_scale` | Time scale | — | 1.0 | 0.1–5 | 0.1 | Wall-clock-to-simulation-time multiplier requested by the user. Achieved rate is reported back and may be lower depending on frame budget. |
| `simulation.seed` | Random seed | — | 1 | 0–1,000,000,000 | 1 | Random number seed for reproducible initial ball configuration. All runs with identical `Params` and `seed` are deterministic. |

**Note:** Fields `simulation.resolution` (default 40, controls PBF spatial resolution) and `simulation.pbf_iterations` (default 3, controls PBF density-constraint solver) exist in the data model but are not exposed in the UI. Frame budget (`simulation.frame_budget_ms`, default 12.0 ms) is also not user-tunable.
