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

## Browser (WASM) achieved-speed measurements and quality presets

The table above is native `cargo bench`, not the actual in-browser experience; this section is. A
second external review measured `Achieved speed` at **0.12x** in-browser at (then-)default
parameters, with the HUD reading `speed: 0.12x` -- a real, reproducible reading, but the review's
own framing ("the M6 performance pass has not been done") was incomplete: `web/src/worker.ts`
recomputed `fluidSurface()` (marching squares over a 128x128 grid) and `metricsJson()` (a full
toe/shoulder/pool/mixing-index pass) **unconditionally every rendered frame**, regardless of how
fast the underlying simulation itself was actually advancing -- pure waste, fixable without any
solver-level optimization at all. Throttling both to `SLOW_UPDATE_INTERVAL_MS = 1000/15` (~15 Hz,
well below render rate, and below the metrics panel's own 100 ms DOM-paint throttle) alone raised
the reading from 0.12x to **~0.35x** at the original default parameters (`max_balls = 600`,
`resolution = 40`) -- measured in-browser, same machine, same params, only that one change.

That alone does not reach 1.0x; the remaining gap is genuine `mill-core` solver cost (DEM + PBF +
coupling), which the full M6 SIMD/wasm-opt/neighbour-list-reuse pass this project has still not
done would address directly. Absent that, the only other lever is fewer particles -- lowering
`simulation.max_balls`/`resolution` trades fidelity for speed along the same curve `docs/PERF.md`'s
native benches already show (`dem_step` scales roughly linearly with ball count; the PBF/coupling
side dominates `drum_only_step` and scales with fluid particle count, roughly `resolution^2`).
Measured in-browser (WASM release build, default drum/media/slurry otherwise, `lifters.count = 0`,
steady cascading after ~10 s of sim time), with the metrics/surface throttling fix already applied:

| `max_balls` | `resolution` | Achieved speed | Fluid spacing / ball diameter (`dx / effective_diameter`) |
|---|---|---|---|
| 150 | 15 | ~1.0-1.1x | ~0.82 |
| 200 | 18 | ~0.9-1.0x (borderline) | ~0.73 |
| 300 | 25 | ~0.7-0.75x | ~0.70 |
| 600 | 40 (mill-core's own `Default`) | ~0.3-0.4x | ~0.62 |

The first, third, and fourth rows are exposed as the "Realtime" / "Balanced" / "Accuracy" quality
presets (`web/src/params/presets.ts`, selectable in the parameters panel's Simulation group); a
fresh app load with no explicit params starts at Realtime (`web/src/worker.ts`), not mill-core's
own `SimulationParams::default()` (kept at the `600`/`40` "Accuracy" values since native
tests/benches on this page compare against it directly). All three presets stay under the
parameters panel's own "a fluid particle is wider than a ball" warning (`> 1.0`), but all three
still under-resolve a *thin* wetting film specifically -- a stricter bar than that warning checks,
and one none of these three tiers clears; see `docs/PARAMETERS.md`'s quality-preset table for the
full disclosure of that trade-off.

Achieved-speed readings above vary run-to-run by several points (observed range roughly ±0.1x at a
fixed configuration on this machine) with ambient system load, the same caveat the native bench
numbers above already carry -- read them as orders of magnitude and relative comparisons, not
precise guarantees, and expect a slower end-user machine to read lower across the board.

## 2026-09-21: re-measured after raising the default slurry fill fraction

`bd19299` ("Raise default slurry fill and make PBF density solve adaptive") raised
`SlurryParams::fill_fraction`'s default from 0.15 to 0.35 and made `pbf_iterations` a floor that
the density solve escalates past (up to `pbf::ADAPTIVE_MAX_EXTRA_ITERATIONS` extra passes) when
still unconverged after a violent sub-step. Neither the quality presets
(`web/src/params/presets.ts`) nor `SimulationParams::default()` override `slurry.fill_fraction` --
all three tiers, and mill-core's own default, inherit whatever the global `SlurryParams` default
is. `FluidParticles::seed_lattice` (`crates/mill-core/src/pbf.rs`) seeds
`target_count = fill_fraction * pi * drum_radius_m^2 / cell_area`, with `cell_area` fixed only by
`resolution`, so this change scales the fluid particle count at *every* fixed `resolution` by
`0.35 / 0.15 ~= 2.33x` -- this was expected to cost real PBF time, not just change the free-
surface height as the commit's own rationale focused on, and it was re-measured rather than
assumed.

Native `cargo bench -p mill-core` (release, same machine as every section above), best-of-three
runs after the first run's cold-start outlier settled (a >2x first-run inflation that vanished on
immediate re-run, the same ambient-noise pattern already documented in the Phase 3b section above):

| Benchmark | Time | vs. fluidised-charge-fix section |
|---|---|---|
| `dem_step/500_balls` | 0.42 ms | no material change (was 0.41 ms) |
| `dem_step/1000_balls` | 1.00-1.12 ms | no material change (was 0.93 ms) |
| `dem_step/2000_balls` | 2.16-2.17 ms | no material change (was 2.00 ms) |
| `drum_only_step` | 52.7-53.0 ms | **+~4x** (was 13.17 ms) |

`dem_step/*` hardcodes its own ball counts and is unaffected by `fill_fraction`, exactly as
expected -- confirmed, not just assumed. `drum_only_step` runs at *default* `Params`
(`slurry.fill_fraction = 0.35`, `simulation.resolution = 40` unchanged), so it is not a clean
single-cause comparison: since the last recorded number six commits landed
(`1cc6cc7` wettability/adhesion, `8ecd966` review fixes, `c81899d` wall/lifter CCD, `b7f0c06`
Akinci cohesion+adhesion replacing the old shell model, `74b85d1` hot-apply/presets, `bd19299`
itself), several of which add real per-substep PBF or DEM work. The fluid-particle-count math above
(2.33x more particles at the same resolution, each now also possibly paying extra adaptive PBF
passes) plausibly accounts for most of a ~4x `drum_only_step` increase on its own; this section
does not attempt to isolate `bd19299`'s share from the other five commits' -- that would need a
`git stash` bisection this task didn't budget for, and the browser measurement below is the number
that actually matters for the UI-facing question this re-measurement exists to answer.

Browser (WASM, freshly rebuilt from this commit via `scripts/build-wasm.sh` -- the committed
`web/src/wasm/` bundle predated `db96784` by a few minutes and was not trusted as-is), same
methodology as the table above (default drum/media/slurry otherwise, `lifters.count = 0`, steady
cascading, HUD/metrics-panel `Achieved speed (x)` reading, several samples per preset):

| `max_balls` | `resolution` | Achieved speed (this measurement) | Achieved speed (prior table) | Est. fluid particles |
|---|---|---|---|---|
| 150 (Realtime) | 15 | ~1.0x (samples: 1.02, 1.00, 1.00) | ~1.0-1.1x | 247 |
| 300 (Balanced) | 25 | ~0.25-0.38x (samples: 0.38, 0.25, 0.37, 0.38) | ~0.7-0.75x | ~690 (extrapolated) |
| 600 (Accuracy) | 40 | ~0.12-0.13x (samples: 0.13, 0.13, 0.13, 0.13, 0.12) | ~0.3-0.4x | 1759 |

The verdict: **the fill-fraction bump (plus the other physics work that landed alongside it) did
measurably shift the achieved-speed numbers, and the shift is well outside this doc's own ±0.1x
noise band.** Realtime is unaffected (~1.0x, matching the prior table) because at `resolution = 15`
the absolute fluid particle count stays small (247) regardless of `fill_fraction`, and Realtime's
cost is DEM/render-bound, not PBF-bound. Balanced and Accuracy both regressed by roughly half:
Balanced ~0.7-0.75x -> ~0.3x, Accuracy ~0.3-0.4x -> ~0.13x -- consistent with PBF cost dominating
at higher `resolution` and scaling with the ~2.33x particle-count increase (plus whatever the
adaptive-iteration escalation and the Akinci cohesion/adhesion passes added on top).
`web/src/params/presets.ts`'s own per-preset `note` strings (last written before `bd19299`) are now
stale for Balanced and Accuracy -- both undersell how slow those tiers actually are -- and should
be refreshed together with this table the next time someone touches that file; this task's scope
was measurement only, so the `note` strings were left as-is pending that follow-up. Realtime's
"default" status stays justified (it is still the only tier at >= 1.0x), but a user who deliberately
picks Balanced or Accuracy today gets a slower experience than the UI copy currently promises.
