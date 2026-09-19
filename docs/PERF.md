# Performance baseline

Recorded from `cargo bench -p mill-core` (release, single-threaded, native — not WASM) before the
physics-fixes pass described in this plan (unit-depth disc masses, buoyancy/drag coupling closure,
implicit viscosity, tunnelling guard, grinding instrumentation). Re-run and update this table after
that work lands so a regression is visible.

Machine: developer workstation (Windows, native `cargo bench`, not WASM — WASM is typically
slower; see docs/PLAN.md ss5 M6 for the real-time targets these numbers are meant to stay under).

| Benchmark | Meaning | Time |
|---|---|---|
| `dem_step/500_balls` | One ball-solver sub-step (`DemState::step`, `dem_iterations = 4`) at 500 balls | 0.77 ms |
| `dem_step/1000_balls` | Same, 1000 balls | 1.32 ms |
| `dem_step/2000_balls` | Same, 2000 balls | 2.46 ms |
| `drum_only_step` | One `Simulation::step(FIXED_DT)` call (one 1/240 s sub-step: DEM + PBF + coupling) at default `Params` | 12.4 ms |

`drum_only_step` is a single fixed sub-step (`FIXED_DT = 1/240 s`), not a full rendered frame —
`Simulation::step(dt)` at the default `simulation.substeps = 4` and `dt = 1/60 s` (one rendered
frame at 1x time scale) would run four of these per frame, so a frame currently costs roughly
`4 * 12.4 ms ~= 50 ms` native, well above the 16.7 ms (60 fps) budget and above the M6 target of
real time at 500 balls + 2000 fluid. This crate has not yet had the M6 performance pass (SIMD,
`wasm-opt`, neighbour-list reuse) applied; these numbers are a pre-optimization baseline for
regression comparison during the physics-fixes work, not a statement that M6 is met.

## After Phase 1 (unit-depth masses, buoyancy/drag coupling closure)

| Benchmark | Time | vs. baseline |
|---|---|---|
| `dem_step/500_balls` | 0.47 ms | slightly faster (cheaper disc mass formula) |
| `dem_step/1000_balls` | 1.04 ms | slightly faster |
| `dem_step/2000_balls` | 2.27 ms | slightly faster |
| `drum_only_step` | 15.0 ms | +21% (new per-ball buoyancy neighbour query, `crates/mill-core/src/pbf.rs` step 6.6) |

## After Phase 2 (implicit Morris/CG viscosity, deterministic `UniformGrid`)

| Benchmark | Time | vs. Phase 1 |
|---|---|---|
| `dem_step/500_balls` | 0.54 ms | +15% (`grid.rs`'s `BTreeMap` switch, shared by the ball-ball grid; see its module doc comment for why — fixes a real non-determinism bug, not a pure regression) |
| `dem_step/1000_balls` | 1.27 ms | +21% |
| `dem_step/2000_balls` | 2.72 ms | +20% |
| `drum_only_step` | 20.0 ms | +33% (conjugate-gradient viscosity solve, `crates/mill-core/src/pbf.rs`'s `solve_implicit_viscosity`, replacing the O(1)-per-particle XSPH pass) |

Cumulative `drum_only_step` change from the original baseline: 12.4 ms -> 20.0 ms (+61%). Still the
same order of magnitude and not yet alarming, but worth watching: the M6 performance pass (not yet
done) will need to budget for the CG solve's iteration count, not just a fixed per-particle cost.

## After Phase 3b (grinding instrumentation: power, collisions, dissipation)

| Benchmark | Time | vs. Phase 2 |
|---|---|---|
| `dem_step/500_balls` | 0.55 ms | no material change |
| `dem_step/1000_balls` | 1.43 ms | no material change (noisy run-to-run on this machine) |
| `dem_step/2000_balls` | 2.55 ms | no material change |
| `drum_only_step` | 17.7 ms | no material change |

Measured back-to-back with a `git stash`/`stash pop` of the Phase 2 commit in the same session, to
rule out cross-session machine-state noise (an earlier same-machine, different-session run of
*this exact Phase 2 commit* showed `dem_step/500_balls` at 1.18 ms and `drum_only_step` at 16.5 ms
-- more than 2x `dem_step`'s originally-recorded Phase 2 number above -- purely from ambient
system load, not any code change). One real regression *was* found and fixed during this phase:
an initial attempt at deterministic ordering for `dem.rs`'s per-substep contact bookkeeping
(`ContactBook`) used the same `HashMap` -> `BTreeMap` swap as `grid.rs`'s fix, but `ContactBook`'s
map is accumulated into (`.entry(key).or_insert(0.0) += ...`) inside step 3's hot solver loop
(`dem_iterations` times per contact pair every sub-step), where `BTreeMap`'s `O(log n)` lookup
cost is significant -- this alone roughly doubled `dem_step`. Fixed by keeping `ContactBook` a
`HashMap` (that accumulation's result does not depend on iteration order -- see its module doc
comment) and instead collecting + sorting its entries once per sub-step for the passes that *do*
need deterministic order (friction/restitution/rolling-resistance, steps 5-7), which cost nothing
measurable at this bench's ball counts.

## After the fluidised-charge energy-injection fix (DEM/coupling physics + new defaults)

| Benchmark | Time | vs. Phase 3b |
|---|---|---|
| `dem_step/500_balls` | 0.41 ms | faster (noisy run-to-run on this machine; no code path removed that would explain a ~25% drop, treat as machine-state noise rather than a real gain) |
| `dem_step/1000_balls` | 0.93 ms | faster, same caveat |
| `dem_step/2000_balls` | 2.00 ms | faster, same caveat |
| `drum_only_step` | 13.17 ms | not directly comparable to Phase 3b's 17.7 ms -- see below |

This round fixed the DEM depenetration/restitution and ball<->fluid coupling energy-injection bugs
that let the charge inflate into a "gas" instead of cascading, and changed
`SimulationParams::default()` to `max_balls = 600` (from 2000), `substeps = 8` (from 4),
`dem_iterations = 2` (from 4). The `dem_step/*` benchmarks hardcode their own ball counts
(`benches/step.rs`'s `effective_media_with_count` helper) and are unaffected by the defaults
change, so they remain a fair apples-to-apples comparison against every earlier section above --
the ~500 balls case simply got a little faster, most plausibly ambient machine-load noise (see the
Phase 3b section's own note on the same kind of run-to-run variance on this machine) rather than a
result of the physics fix itself, since the fixed depenetration/restitution logic changed magnitude
clamps, not algorithmic cost.

`drum_only_step` calls `Simulation::step(FIXED_DT)` (`FIXED_DT = 1/240 s`, a crate-level constant,
*not* derived from `simulation.substeps`) at *default* `Params`. `Simulation::step(dt)` splits `dt`
into `simulation.substeps` equal-size internal sub-steps (`crates/mill-core/src/lib.rs`'s
`step`/`advance_one_sub_step`), so this one call always advances `FIXED_DT` of sim time total, but
now runs that as **8** internal sub-steps of `FIXED_DT / 8` each (previously 4 sub-steps of
`FIXED_DT / 4` each, back when `FIXED_DT` and the true fixed sub-step size happened to coincide at
the old default `substeps = 4` -- see `Simulation::fixed_sub_dt`'s doc comment). Each of those 8
sub-steps now also does DEM work for only 600 balls at `dem_iterations = 2`, instead of 4 sub-steps
of ~2000 balls at `dem_iterations = 4`. Three things changed at once relative to Phase 3b (ball
count, sub-step count, DEM iteration count), so the new number (13.17 ms) is a defaults-and-
physics-fix snapshot, not a clean isolated regression check against Phase 3b's 17.7 ms -- it
happens to come out lower, consistent with the ball-count and DEM-iteration drops outweighing the
doubled sub-step count, but that's a coarser statement than "X% faster/slower for the same reason"
the other three rows above support.
