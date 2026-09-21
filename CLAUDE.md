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
frontend, Canvas 2D, a persistent left-side parameters panel (`<aside class="params-panel">`,
`web/src/ui/paramsPanel.ts`; replaced the earlier `<dialog>` modal). Default drum wall has
**no lifters**.
See docs/PLAN.md, docs/PHYSICS.md, docs/PARAMETERS.md.

## Repository
- Remote: https://github.com/ToshihiroIguchi/MillDynamics3 (branch `main`). Commit per milestone; push only when the user asks.

## Conventions
- SI units everywhere in the core (m, kg, s, Pa·s). UI may show mm / rpm and convert in `schema.ts`.
- Deterministic: every run is reproducible from `Params` + `seed`.
- `cargo fmt`, `cargo clippy -D warnings`, `npm run lint` must pass before finishing a task.
- Do not commit `web/src/wasm/` build output or `target/`.
