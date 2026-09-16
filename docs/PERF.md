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
