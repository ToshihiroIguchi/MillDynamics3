# Reading the metrics panel

A guide to what each number in the live metrics panel (`web/src/ui/metricsPanel.ts`, driven by
`web/src/metrics/specs.ts`) actually measures, what a healthy reading looks like at this project's
current defaults, and — the point of this document — which readings are *expected* to look
alarming to someone unfamiliar with this project's specific solver choices, and why. It exists
because an external review of a browser session misread several of these numbers as evidence of
broken physics; each misread is called out explicitly below with the corrected reading, so the next
reviewer doesn't repeat the same diagnosis. See `docs/PHYSICS.md` for the underlying algorithms and
exact formulas this document summarises in plain terms.

All readings below were measured in-browser at default `Params` (drum 1 m, 30 rpm / 71% Nc, no
lifters unless noted, media fill 0.30, slurry fill 0.15 / viscosity 50 Pa·s / wettability 0.6),
after the initial settling transient (roughly `t > 15` s of simulated time).

---

## Grinding group

### Power draw (W/m), Torque (N·m/m)

The mill's instantaneous mechanical power/torque draw, from the drum wall's own friction work on
the charge. **Healthy range**: a few hundred to a few thousand W/m depending on fill/speed/lifters;
zero (or near it) if centrifuging (see `docs/PHYSICS.md` §8, `centrifuging_draws_less_power_than_cascading`).

### Dissipated power (W/m)

Net mechanical energy the DEM contact solve removes from the ball population, per second — friction,
restitution `e < 1`, and rolling resistance combined, minus (via `wall_work_j`, see
`docs/PHYSICS.md` §4.2) whatever the wall's own motor supplied that sub-step. **This should always
read positive** (up to floating-point slack) and of the same rough order as *Power draw* in a dry
(no-slurry) steady state, since energy in ≈ energy out once the charge's own kinetic + potential
energy stops trending up or down. **A previous version of this metric could read hundreds of watts
negative** on an entirely ordinary cascading run — not a thermodynamics violation, but a metric
definition bug (it omitted the wall's own energy input, and omitted rotational KE and gravitational
PE from the energy balance). Fixed; see `docs/PHYSICS.md` §4.2's `dissipated_energy_j` description
and the regression tests
`dem::tests::dissipated_energy_is_never_materially_negative_while_cascading` /
`tests::dissipated_power_approaches_power_draw_in_a_settled_dry_charge`.

**With slurry enabled, expect `Dissipated power < Power draw`** — this metric is purely DEM-side and
does not include the fluid's own viscous dissipation of the balls' kinetic energy (a real,
separate energy sink; see `docs/PHYSICS.md` §9). Measured at defaults: Power draw ~3300-3700 W/m,
Dissipated power ~4000-4100 W/m (slurry coupling impulses are folded into the DEM predict step, so
this particular reading can occasionally exceed Power draw too — both are within the same order of
magnitude, which is the property that matters).

---

## Media group

### Max ball overlap (of radius)

Largest ball-ball penetration depth as a fraction of the ball **radius**, not diameter —
`(2r - dist) / r`. Two ball centres exactly coincident would read `200%`; a half-diameter overlap
reads `100%`, not `50%`. This is XPBD (zero-compliance position-based), not a spring-dashpot DEM —
overlap here is a converged-solver-quality/contact-accuracy indicator, not a stability indicator
(the solver is unconditionally stable at any overlap; see `docs/PHYSICS.md` §4). **Measured at
defaults: roughly 20-45%** (of radius) at `dem_iterations = 2`, the current default — a packed,
actively-cascading charge does not fully converge its contact network in only 2 Gauss-Seidel
passes per sub-step. This is a real, known trade-off (`dem_iterations` was halved, from 4 to 2,
when `substeps` was doubled, to keep DEM cost roughly flat while improving the tunnelling-risk
ratio below) — raising `dem_iterations` back up reduces this at a real-time performance cost; left
as-is pending the M6 performance pass creating headroom. A companion **Max ball-wall overlap (of
radius)** reads much smaller with no lifters (~1-2%) but **~31-33% with `lifters.count = 8`** — see
`docs/PHYSICS.md` §9 for that specific, separately-noted finding.

---

## Slurry group

### Max/mean compression error

The PBF density constraint is deliberately **one-sided**
(`(rho_i/rho0 - 1).max(0.0)`, `docs/PHYSICS.md` §5.2 step 3): it only ever pushes *over*-dense
particles apart, and never tries to pull an under-dense one back together (that used to be the
solver's only fluid-fluid attraction, and it balled slurry into mid-air blobs at the free surface —
see `docs/PHYSICS.md` §5.2's `c_i` discussion). This pair of metrics reports exactly that one-sided
quantity, over the whole fluid population — the honest convergence readout for "is the
incompressibility constraint actually converged." **Healthy: mean well under 1% at defaults; the
max can spike higher (10-20%) during active splashing/impact, settling back down between events** —
this is not a computed-once static reading, it's Sampled every frame during active motion, so a
transient impact naturally shows up. A settled, undisturbed puddle converges to mean well under
1%, max under a few percent (`metrics::tests::settled_puddle_compression_error_is_small`).

### Density spread, max/mean (incl. free surface)

The *absolute-value* reading, `|rho_i/rho0 - 1|`, over the whole population — kept for CSV/back-
compat and because it is a real "how ragged is the density field" signal, including regions the
solver deliberately does not correct. **Structurally large on a completely healthy run**: every
free-surface and near-wall particle has fewer SPH neighbours than an interior particle by simple
geometry, in *any* SPH-family method, so its measured density reads low — and the one-sided
constraint above leaves that alone by design. **Measured at defaults: max ~72%, mean ~15-19%.**

**This is exactly the pair of numbers an external review read as "the incompressibility constraint
has failed" (72.4% max, 16-26% mean, matching the measured range here almost exactly).** It hadn't
failed — see *Max/mean compression error* above for the number that actually tracks convergence,
which stays under 1% mean over the same run. The absolute-value reading is dominated by the free
surface's inherent neighbour deficiency, not solver divergence.

---

## Solver group

### Substep displacement / diameter

`|v| * dt / (2 * radius)` — the largest fraction of its own diameter any ball moves in one
sub-step, a tunnelling-risk diagnostic (`docs/PHYSICS.md` §9). **The reviewer's expectation that
this "must be < 0.05-0.1" is the criterion for an explicit spring-dashpot DEM contact model,
which this project does not use.** This solver is XPBD (zero-compliance position-based, no
tangible collision "speed" to overshoot) — the documented criterion here is `< 1`, not `< 0.1`;
values well under 1 are healthy, not merely tolerable. **Measured at defaults: ~0.16-0.21 with no
lifters, ~0.30 with `lifters.count = 8`** (cataracting off a lifter reaches higher peak ball
speeds) — both comfortably under 1. A reviewer's `0.402` reading (with lifters enabled) is
consistent with this pattern, not a regression.

A real, separate gap noted in `docs/PHYSICS.md` §9: the broad-phase margin that keeps a fast pair
from being lost between sub-steps is sized off the single fastest ball's speed, not the worst-case
pair's *closing* speed — a head-on pair can in principle close at up to twice that margin's
assumption. Not yet fixed (would cost more candidate pairs per sub-step; gated on the M6
performance pass).

### Coupling clamp hits

Count of balls whose fluid-coupling impulse hit the stability clamp this call. **Healthy: 0, or
very near it, once the charge has settled** — persistent clamping indicates the coupling forces are
pinned at an artificial ceiling rather than reflecting the physical interaction. The Phase 3
adhesion-gating fix (below) only ever *reduces* adhesion strength, so it cannot raise this rate;
verified by `coupling::tests::cascading_charge_keeps_coupling_clamp_hits_rare_once_settled_at_{50,200}_pa_s`.

---

## Not a metrics-panel number, but related: the airborne slurry clump and the wetting film

Not a displayed metric, but two concrete defects external reviews reported, both traced back to
the same root cause: the old shell-based adhesion mechanism had no fluid-fluid cohesion
counterpart, so it could only ever be tuned to trade one defect for the other.

- **A ball flying through the air (cataracting) could keep a small cluster of slurry particles
  glued to it indefinitely.** Caused by the old mechanism's shell reaching too far from the
  surface (up to ~2.6 ball radii at this project's typical coupling resolution) with nothing
  testing whether there was any actual bulk liquid there.
- **After a density-based gate was added to fix the clump defect, wetting itself visibly
  weakened.** The gate suppressed adhesion for any particle whose local SPH density read low --
  which is exactly what a genuine thin wetting film reads as, being sub-resolution relative to the
  fluid's own kernel radius `h` at this project's coupling resolution. The two defects were not
  independently fixable with a single ball<->fluid-only coefficient: strengthening the pull to
  restore wetting reintroduced the floating-clump defect, and weakening it to suppress the clump
  starved the film.

Both are now fixed together by replacing the mechanism entirely with an Akinci-style pairwise
cohesion (`slurry.surface_tension_n_m`, new) and adhesion (`slurry.wettability`) pair
(`docs/PHYSICS.md` §6.1a): adhesion's range is capped at `2 * balls.radius` regardless of `h`
(fixing the clump defect geometrically, not via a density gate that also suppresses real films),
and cohesion gives the slurry a genuine, separately-tunable attraction that holds a pulled-on film
together against gravity instead of relying on adhesion alone to do both jobs. Regression tests:
`coupling::tests::adhesion_range_is_capped_at_twice_the_ball_radius_regardless_of_fluid_resolution`,
`pbf::tests::cohesion_pulls_two_isolated_fluid_particles_together_only_when_surface_tension_is_positive`.

---

## Known, deliberate, out-of-scope items (not defects)

These were also raised by the same review and are documented project scope decisions, not bugs:

- **Monodisperse media** (no ball size distribution) — `docs/PHYSICS.md` §9, reserved for a future
  milestone (M7).
- **Real-time performance** — the M6 SIMD/wasm-opt/neighbour-list-reuse pass has not been done yet;
  `Achieved speed` reading well under `1.0x` (commonly 0.2-0.5x at default parameters) is the known,
  current state, not a regression. See `docs/PERF.md`.
- **2D cross-section only** — every value in this crate is per metre of unit mill depth; no axial
  transport, end-wall friction, or 3D particle shape. A deliberate, project-wide scope decision.
- **Parameter changes reset the simulation** — by design; not a bug to fix.
