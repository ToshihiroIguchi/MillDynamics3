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
