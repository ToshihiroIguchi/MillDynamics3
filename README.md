# MillDynamics3

A real-time, browser-based simulator of a tumbling ball mill cross-section: grinding media (balls)
modeled with soft-sphere DEM and a viscous slurry modeled with Position Based Fluids (PBF), coupled
two-way. The simulation core is written in Rust and compiled to WebAssembly; the frontend is a
framework-free Vite + TypeScript app rendering to Canvas 2D.

See [`docs/PLAN.md`](docs/PLAN.md) for the full design and milestone plan, and
[`CLAUDE.md`](CLAUDE.md) for project conventions.

## Status

Early development (see `docs/PLAN.md` milestones M0–M7). Not yet runnable.

## Live demo

https://toshihiroiguchi.github.io/MillDynamics3/

`main` auto-deploys to GitHub Pages via GitHub Actions (`.github/workflows/deploy.yml`) on every
push, gated on the same checks as `## Development` below.

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
