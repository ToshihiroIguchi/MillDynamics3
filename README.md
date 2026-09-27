# MillDynamics3

A browser-based simulator of a tumbling ball mill cross-section: grinding media (balls)
modeled as position-based (XPBD) rigid discs, and a viscous slurry modeled with Position Based Fluids (PBF),
coupled two-way with buoyancy, viscous drag, and Akinci-style cohesion/adhesion. The simulation core is written in Rust
and compiled to WebAssembly; the frontend is a framework-free Vite + TypeScript app rendering to Canvas 2D.
Targets real time (>= 1.0x) at the default "Realtime" quality preset on typical desktop hardware; two
coarser-grained alternatives (Balanced, Accuracy) trade that speed for fidelity -- see
`docs/PARAMETERS.md`'s quality-preset table and `docs/PERF.md` for the measured trade-off each one makes.
The HUD always reports the actually-achieved speed, not just the setting, and a "SLOW MOTION" banner
appears unmissably whenever the achieved rate falls behind.

See [`docs/PLAN.md`](docs/PLAN.md) for the full design and milestone plan, and
[`CLAUDE.md`](CLAUDE.md) for project conventions.

## Status

Functional simulator with DEM ball solver (position-based rigid discs, wall/lifter continuous collision
detection), PBF slurry solver with implicit Newtonian viscosity and Akinci-style cohesion/adhesion, two-way
ball–fluid coupling, free-surface rendering, a live metrics panel (power draw, mixing index, collision
energy histogram, CSV export), and a parameters panel with three quality presets. See `docs/PLAN.md`
milestones M0–M7 for roadmap. Most of M6's performance work has landed: WASM SIMD (`+simd128`),
`wasm-opt`, and neighbour-list/broad-phase reuse via a counting-sort CSR spatial grid (replacing a
`BTreeMap`), plus removing per-substep allocations from the PBF/coupling hot path (see the 2026-09-27
section of `docs/PERF.md` for the full before/after). The "Realtime" preset reaches >= 1.0x on a
typical dev machine mainly by throttling the per-frame metrics/free-surface recomputation to well
below render rate and reducing the default ball/fluid particle counts; Balanced/Accuracy are real-time
in neither preset even after M6's algorithmic work, since PBF/coupling cost still scales with fluid
particle count at those tiers' higher resolution (see `docs/PERF.md`). Not done: a preconditioned
(Jacobi) viscosity solve and multi-rate/looser-tolerance options, both left as future, fidelity-
trading work rather than pure-speed wins.

## Live demo

https://toshihiroiguchi.github.io/MillDynamics3/

`main` auto-deploys to GitHub Pages via GitHub Actions (`.github/workflows/deploy.yml`) on every
push, gated on the same checks as `## Development` below.

## Performance

These are native `cargo bench` results (single-threaded, release, not WASM—WASM is typically slower).
The table below is kept as originally recorded, after Phase 3b (grinding instrumentation) and before
the 2026-09-27 M6 solver pass -- absolute numbers vary a lot machine-to-machine (see that section's
own caveat), so it is not directly comparable to a fresh `cargo bench` run today; it is retained here
as the historical entry point into `docs/PERF.md`'s full dated history, which has the actual
before/after numbers for every change including M6's.

| Benchmark | Meaning | Time |
|---|---|---|
| `dem_step/500_balls` | One ball-solver sub-step at 500 balls | 0.55 ms |
| `dem_step/1000_balls` | One ball-solver sub-step at 1000 balls | 1.43 ms |
| `dem_step/2000_balls` | One ball-solver sub-step at 2000 balls | 2.55 ms |
| `drum_only_step` | One full sub-step (DEM + PBF + coupling) at default params | 17.7 ms |

All masses, energies, and power values in the core are per metre of mill axial length (2D unit-depth convention).

## Prerequisites

- [Rust](https://www.rust-lang.org/tools/install) (stable) with the `wasm32-unknown-unknown` target
- [`wasm-pack`](https://rustwasm.github.io/wasm-pack/installer/)
- [Node.js](https://nodejs.org/) 20+ and npm

```sh
rustup target add wasm32-unknown-unknown
cargo install wasm-pack
```

## Build & run

```sh
# Build the WASM core (release)
./scripts/build-wasm.sh      # or .ps1 on Windows

# Run the web app
cd web
npm install
npm run dev
```

## Development

```sh
# Rust: format, lint, test, benchmark
cargo fmt
cargo clippy --all-targets -- -D warnings
cargo test -p mill-core
cargo bench -p mill-core

# Web: lint, unit tests, e2e tests
cd web
npm run lint
npm run test
npm run e2e  # requires scripts/build-wasm.sh (or .ps1 on Windows) to have been run first
```

## Notes

- The simulation is a 2D cross-section of the mill (perpendicular to the rotation axis), not a full
  3D model. Results are qualitative, not a quantitative replacement for a 3D DEM/CFD study.
- Default drum wall is smooth (no lifters); lifters can be enabled and configured in the Parameters
  modal.
