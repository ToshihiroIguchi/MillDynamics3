# MillDynamics3 — Physics Reference

Equations and algorithms actually implemented in `crates/mill-core/src/{params,geometry,grid,dem,
pbf,coupling,surface,metrics}.rs`, as of the `wip/physics-fixes-metrics-panel` branch. This
documents the real, current source — not the design narrative in `docs/PLAN.md` ss3, whose ss3.2
and ss3.4 text predates the coarse-graining/mass-model and coupling fixes described below (call
sites are cited so the two can be cross-checked; where they disagree, this document and the source
win).

See `docs/PLAN.md` ss3 for the original design rationale (solver-choice tables, milestone history)
and `docs/PARAMETERS.md` for the UI-facing parameter list.

---

## 1. Units & conventions

- SI throughout the core: metres, kilograms, seconds, Pascal-seconds, radians (`params.rs` module
  doc comment). The UI may display mm/rpm and converts in `web/src/params/schema.ts`; nothing in
  `mill-core` itself uses a non-SI unit.
- **2D cross-section, unit depth-of-1-metre convention.** Every ball and fluid particle is modelled
  as a unit-depth (1 m) slice, not a 3D sphere/voxel. Consequently every mass, energy and power
  value produced anywhere in this crate is *per metre of mill length* (`dem.rs` module doc comment
  and `DemStepStats` doc comment); a caller wanting an absolute value for a real mill multiplies by
  the mill's actual axial length.
- World frame: drum centre at the origin, gravity along `−y` (`GRAVITY = -9.81` m/s^2 in both
  `dem.rs` and `pbf.rs`), drum rotates about the origin.
- Angle convention used consistently by `geometry.rs`, `dem.rs` and `metrics.rs`: `atan2(y, x)` in
  `[0, 2*pi)`, measured counter-clockwise from `+x`; "down" (6 o'clock) is `3*pi/2`. `metrics.rs`
  additionally exposes `to_vertical_degrees` to convert to the mill-literature "0 deg at 12
  o'clock, clockwise-increasing for CCW rotation" convention used by acceptance-test literature
  values.
- Determinism: every run is reproducible from `Params` + `seed` (`rng.rs`'s small deterministic
  PRNG seeds initial lattice jitter; `grid.rs`'s `UniformGrid` uses a `BTreeMap`, not `HashMap`,
  specifically so broad-phase candidate-pair iteration order — and hence f32 summation order in the
  contact/kernel accumulations — is identical run-to-run for identical input; see that module's
  doc comment for the full argument). `crates/mill-core/src/dem.rs`'s `ContactBook` similarly uses
  a `HashMap` for the hot per-substep accumulation (measured ~2x faster than `BTreeMap` at default
  ball counts) but replays results through a sorted `Vec` (`sorted_ball_ball`/`sorted_ball_wall`)
  before any pass that iterates *across* contacts and mutates shared state sequentially (friction,
  restitution, rolling resistance), which is where non-deterministic order would actually change
  the physical result.
- Fixed sub-step: `crate::FIXED_DT = 1/240 s` (`lib.rs`), shared by the ball solver and the PBF
  fluid solver. `Simulation::step(dt)` splits a caller-supplied `dt` into `simulation.substeps`
  equal sub-steps (`sub_dt = dt / substeps`, not literally `FIXED_DT` unless `dt` is exactly
  `substeps / 60`; the "1/240 s at 1x time scale" figure assumes the nominal 60 Hz frame rate the
  UI targets). `simulation.dem_iterations` and `simulation.pbf_iterations` control the respective
  solvers' per-substep Gauss-Seidel/Jacobi iteration counts.

---

## 2. Drum geometry & SDF (`geometry.rs`)

The drum (`Drum { radius_m, omega, lifters }`) is represented in its own rotating frame; callers
rotate a world-space point into that frame by `-drum_angle` before querying it.

- **Base wall.** `dist_to_wall(p_local) = radius_m − |p_local|`; positive = free space, negative =
  inside solid (`Drum::sdf`).
- **Lifters.** `lifters.count == 0` (the default) means a perfectly smooth circle: `sdf = dist_to_wall`
  exactly, with no extra evaluation cost. For `count = n > 0`, each of the `n` evenly-spaced bars is
  a trapezoidal cross-section (`lifter_cross_section_sdf`): in its own local (radial `pr`,
  tangential `pt`) frame it spans `pr` in `[R − height, R]` with tangential half-width `base_width/2`
  at the base (`pr = R`) and `top_width/2` at the tip (`pr = R − height`); a possibly-degenerate
  trapezoid (`top_width == base_width` is a rectangle, `top_width == 0` a triangle). The signed
  distance to the convex quad is the minimum, over its four counter-clockwise-wound edges, of the
  point's distance to that edge's infinite line along the edge's inward normal — exact on-edge,
  an approximation near corners (a standard, documented simplification for a `min`/`max`-combined
  convex SDF). `Drum::sdf_lifters` takes the maximum "inside" value over all `n` bars (nearest
  lifter wins) and negates it to match the wall's positive-outside convention; `Drum::sdf` then
  takes `dist_to_wall.min(sdf_lifters(..))` (nearest solid boundary, of wall or any lifter, wins).
- **Normal.** `Drum::sdf_world(p_world, drum_angle)` rotates into the local frame, evaluates `sdf`,
  and takes a central-difference numeric gradient (`EPS = 1e-4`) to get the local normal, falling
  back to `-p_local.normalize_or_zero()` if the gradient is degenerate (e.g. exactly at the drum
  centre); the normal is then rotated back to world space. Points from the wall towards free space.
  A non-finite `p_world` short-circuits to `(0.0, Vec2::ZERO)` rather than letting the numeric
  gradient manufacture a NaN from `-inf − -inf` (documented in detail on `sdf_world`) — this keeps a
  caller's typical `p += -d * normal` correction a safe no-op instead of injecting NaN.
- **Wall velocity.** `wall_velocity(p) = omega * perp(p) = omega * (-p.y, p.x)`, i.e. `omega x r` in
  2D, evaluated at the world-space point currently coincident with that wall material point.

---

## 3. Media / coarse-graining (`params.rs`)

Balls are seeded from an *effective* media population, not the raw UI ball diameter, whenever the
true population would be too large to step in real time.

- **True disc count** (`Params::true_ball_count`):
  `N_true = fill_fraction * packing_fraction_2d * drum_area / (pi * r_true^2)`
  where `drum_area = pi * radius_m^2` and `r_true = ball_diameter_m / 2`. `fill_fraction` (the
  conventional "J" ball-filling fraction) is the charge footprint *including voids*;
  `packing_fraction_2d` (default `0.82`, valid range `[0.5, 0.907]`) converts that footprint into a
  solid disc count — random close packing of equal discs is ~0.82, `pi / (2*sqrt(3)) ~= 0.907` is
  the hexagonal upper bound.
- **Coarse-graining** (`Params::effective_media`): if `N_true > simulation.max_balls`, the true
  population is replaced by `N_sim ~= max_balls` larger, lighter discs with scale factor
  `k = sqrt(N_true / max_balls)`:
  - `diameter_m = true_diameter_m * k` (preserves total footprint area, `N_sim * pi * r_eff^2 ~=
    N_true * pi * r_true^2`, i.e. the fill fraction the user set).
  - `density_kg_m3` is **left unchanged** at the true media density. Because balls are modelled as
    unit-depth discs (`ball_mass = rho * pi * r^2`, ss4 below), preserving the total solid
    footprint area at constant density automatically preserves total charge mass (`mass ~ r^2`, and
    `N_sim * r_eff^2 == N_true * r_true^2` by construction) — no compensating `1/k` density scaling
    is needed or applied. (`Params::effective_media`'s doc comment states this explicitly, and
    `effective_media_preserves_total_charge_mass`/`effective_media_coarse_grains_when_true_count_is_large`
    in `params.rs`'s test module assert it directly.)
  - `ball_count = round(N_true / k^2)` (`~= max_balls`).
  - When `N_true <= max_balls`, `scale_factor = 1.0` and the true media parameters pass through
    unchanged.
- This is a standard coarse-grained DEM approximation: it preserves bulk charge mass and footprint
  area, not the true interstitial void structure or single-collision statistics at the real
  particle size. `EffectiveMedia { true_diameter_m, diameter_m, density_kg_m3, ball_count,
  scale_factor }` carries both the true and simulated values so the approximation is never silent
  (surfaced in the UI's derived-values panel and in `metrics::Metrics`'s
  `true_ball_count`/`simulated_ball_count`/`coarse_graining_factor`/`effective_ball_diameter_m`).

---

## 4. Ball solver — XPBD rigid discs (`dem.rs`)

Non-penetration (ball-ball and ball-wall/lifter) is solved as a rigid, zero-compliance geometric
position constraint, projected iteratively — unconditionally stable at any sub-step size, unlike an
explicit spring-dashpot contact (whose stable time step would need to be on the order of
microseconds at this project's wall speeds/ball stiffness; see `docs/PLAN.md`'s intro section for
the derivation of why that was rejected). Friction, restitution and rolling resistance are added as
position/velocity corrections on top of the converged contact solve.

### 4.1 Ball properties

- `ball_mass(diameter_m, density_kg_m3) = density_kg_m3 * pi * r^2` (`r = diameter_m / 2`) — a
  unit-depth disc's areal mass, kg per metre of mill length. This is the same 2D-slice convention
  the fluid uses (`particle_mass = rest_density * dx^2`, ss5.1), which is what makes ball<->fluid
  momentum exchange (ss6) dimensionally meaningful. `dem.rs`'s doc comment on `ball_mass` notes an
  earlier version used a 3D-sphere mass here instead, which made a fluid particle hundreds of times
  heavier than a coarse-grained ball.
- `ball_inertia(mass, radius_m) = 0.5 * mass * radius_m^2` — moment of inertia of a **uniform disc**
  about its centre (`I = 1/2 m r^2`), consistent with the unit-depth-disc mass model above. This is
  *not* the 3D-sphere value `2/5 m r^2`.
- All balls in a population currently share one radius/mass/inertia (`Balls` struct; a size
  distribution is future work), seeded on a hexagonal lattice filling the drum's cross-section with
  small random jitter (`Balls::seed_lattice`, deterministic via `seed`).

### 4.2 Per-substep algorithm (`DemState::step_with_external_forces`)

Given a fixed sub-step `dt`, drum pose, media parameters and iteration count:

1. **Predict.** For each ball: `v.y += GRAVITY * dt`; if external (fluid coupling) impulses are
   supplied (`Option<&CouplingImpulses>`), `v += impulse * inv_mass` and `omega += angular_impulse *
   inv_inertia` are applied here too — as a direct velocity change, with **no extra `* dt`** (see
   ss6 for why). Then `x += v * dt`, `theta += omega * dt`. Pre-predict velocity (`v_pre`) and
   pre-step position/orientation (`x0`, `theta0`) are saved for later steps.
2. **Broad-phase.** A `UniformGrid` is built over the predicted ball positions with cell size
   `2 * r * 1.05`; `for_each_candidate_pair` enumerates ball-ball candidate pairs once, reused
   across all solver iterations this sub-step. Ball-wall contact is checked directly per ball
   (no broad-phase needed: one `Drum::sdf_world` query per ball).
3. **Solve non-penetration** (`dem_iterations` Gauss-Seidel passes, zero compliance):
   - Ball-ball: `C = |x_i - x_j| - (r_i + r_j)`; if `C < 0`, `d_lambda = -C / (w_i + w_j)` (`w =
     1/mass`), apply `x_i += w_i * d_lambda * n_hat`, `x_j -= w_j * d_lambda * n_hat`, accumulate
     `lambda_n` for that `(i, j)` pair in `ContactBook`.
   - Ball-wall/lifter: `C = drum.sdf_world(x_i) - r_i`; if `C < 0`, `d_lambda = -C / w_i` (wall has
     `w_wall = 0`, infinite mass), apply `x_i += w_i * d_lambda * n_hat`, accumulate `lambda_n` for
     ball `i`.
4. **Reconstruct velocities**: `v = (x - x0) / dt`, `omega = angle_diff(theta, theta0) / dt`.
5. **Friction** (one pass, Coulomb-clamped position correction — not a persistent tangential
   spring), iterating contacts in a deterministic sorted order:
   - Ball-ball: `v_t = (v_i - v_j).dot(t_hat) - r*omega_i - r*omega_j`; the raw correction that
     would fully cancel `v_t` this sub-step is `raw = -v_t * dt / w_sum_t` (`w_sum_t = w_i + w_j +
     r^2*w_rot_i + r^2*w_rot_j`), clamped to `[-mu*lambda_n, +mu*lambda_n]`
     (`media.friction_ball_ball`) and applied as a position correction (linear on `x`, rotational on
     `theta`, both scaled by the respective inverse mass/inertia).
   - Ball-wall: same form, with the wall's own `wall_velocity` substituted for the "other body"'s
     velocity and no wall-side rotational term; `media.friction_ball_wall` is the clamp coefficient.
     The tangential impulse `d_lambda_t` dotted with the wall's velocity there gives the work done
     against wall friction this contact (`wall_work_j`, summed over all ball-wall contacts) — the
     wall's *normal* impulse does no work since the wall's velocity is purely tangential.
   - Velocities are reconstructed a second time from the friction-corrected positions/orientations.
6. **Restitution** (one pass, using the *pre-solve* approach velocity `v_pre`), applied only to
   contacts whose pre-solve normal approach speed exceeded `RESTITUTION_VELOCITY_THRESHOLD = 0.02`
   m/s (below this, a contact is treated as already-resting so restitution does not re-fire every
   sub-step and cause a resting contact to buzz):
   - Ball-ball: `v_n_pre = (v_pre_i - v_pre_j).dot(n_hat)`; if approaching faster than the
     threshold, the target post-solve separating speed is `-e * v_n_pre`
     (`e = media.restitution_ball_ball`), applied as a normal-velocity impulse split by inverse
     mass. Counted as one collision in `DemStepStats::collision_count`, with impact energy `E = 0.5
     * (mass/2) * v_n_pre^2` (reduced mass `mass/2` for two equal masses) binned into
     `impact_energy_histogram`.
   - Ball-wall: same but with the wall's "infinite mass" (reduced mass is just the ball's own,
     `E = 0.5 * mass * v_n_pre^2`), `e = media.restitution_ball_wall`.
7. **Rolling resistance**: for each ball with nonzero `omega` and nonzero accumulated normal
   impulse (`total_lambda_n`, summed over its ball-ball and ball-wall contacts this sub-step), the
   angular deceleration is capped at `max_delta_omega = media.rolling_friction *
   (total_lambda_n/dt) * r / inertia * dt`, applied toward zero and clamped so it cannot overshoot
   past `omega = 0` (i.e. cannot reverse the sign of `omega` in one sub-step) — a dimensionless
   torque-coefficient model of rolling friction.

`DemStepStats` (per-substep diagnostics, all per-metre-of-mill-depth) returned by this function:
`wall_work_j`, `collision_count`, `impact_energy_histogram` (12 log-spaced bins from `1e-6` J to
`1.0` J, `impact_energy_bin_edges`), `dissipated_energy_j` (total ball KE right after predict, i.e.
including gravity/external forces but before any contact correction, minus final KE — a net,
mechanism-agnostic dissipation estimate that does not separately attribute loss to friction vs.
restitution vs. rolling resistance), and `max_substep_displacement_over_diameter` (see ss9).

---

## 5. Slurry solver — PBF + implicit viscosity (`pbf.rs`)

Position Based Fluids with 2D Poly6 (density) / Spiky-gradient (pressure-gradient) kernels, an
iterative density-constraint (incompressibility) solve, drum-wall boundary projection with
no-slip velocity blending, and an implicit (backward-Euler, conjugate-gradient) Newtonian viscosity
solve. `s_corr` artificial pressure is implemented but disabled.

### 5.1 Seeding & kernels

- `FluidParticles::seed_lattice`: spacing `dx = drum_radius_m / resolution`, kernel radius `h = 2 *
  dx`, `rest_density = slurry.density_kg_m3`, `particle_mass = rest_density * dx^2` (2D unit-depth
  convention, matching `ball_mass`). Particles are placed on a hexagonal lattice filling
  `slurry.fill_fraction` of the drum's cross-sectional area from the bottom upward; lattice sites
  overlapping an already-seeded ball are skipped.
- **Poly6** (density kernel): `W(r^2, h) = 4/(pi*h^8) * (h^2 - r^2)^3` for `r < h`, else `0`.
- **Spiky gradient**: `grad_W(delta, r, h) = -30/(pi*h^5) * (h-r)^2 * (delta/r)` for `0 < r < h`,
  else `0` — points from `x_i` toward its neighbour (Spiky decreases with `r`).
- `FluidParticles::densities()`: `rho_i = sum_j m_j * W_poly6(r_ij, h)`, including the particle's
  own self-contribution at `r = 0`.

### 5.2 Per-substep algorithm (`FluidParticles::step_coupled`; `step` is the ball-free special case)

0. **Sanitize.** Any incoming non-finite `x`/`v` (which can only have entered from outside this
   struct, e.g. a diverged ball position read in steps 3.5/6.5) is reset to `(Vec2::ZERO,
   Vec2::ZERO)` before anything else runs — every correction this solver itself performs is
   individually bounded, so this guard exists purely against externally-introduced NaN/inf.
1. **Predict**: `v.y += GRAVITY * dt`; `x += v * dt`.
2. **Neighbour lists**: a `UniformGrid` over the predicted positions with cell size `h`;
   `build_neighbor_lists` keeps only candidate pairs actually within `h` (mutual: `j` in
   `neighbors[i]` implies `i` in `neighbors[j]`), rebuilt once per sub-step and reused across every
   density-constraint iteration.
3. **Density-constraint solve** (`pbf_iterations` Jacobi-style passes — all deltas computed from
   current positions, then applied together):
   - `C_i = rho_i/rest_density - 1`.
   - `lambda_i = -C_i / (sum_j |grad_W_ij/rest_density|^2 + grad_self^2 + EPSILON_RELAX)`,
     `EPSILON_RELAX = 200.0` (CFM-style relaxation against divide-by-zero for sparse
     neighbourhoods).
   - `delta_p_i = 1/rest_density * sum_j (lambda_i + lambda_j + s_corr_ij) * grad_W_ij`, applied
     directly to `x`.
   - **`s_corr` is implemented but disabled** (`S_CORR_K = 0.0`, formula kept as `-k *
     (W(r)/W(delta_q))^n` with `delta_q = 0.2h`, `n = 4`): at this project's real-SI scale
     (`rest_density` ~1000-2000 kg/m^3, small `h`), the literature-default `k = 0.1` produced a
     correction 10-30x larger than the density-constraint terms it's meant to gently supplement,
     causing runaway dispersal instead of preventing clustering. Disabling it gives a stable
     settled puddle at ~1% mean density error over 2 s of simulated time; several smaller `k`
     values were tried without finding a stable working point. Documented in `pbf.rs` as revisit-if
     visual clustering artefacts are observed, not resolved.
4. **Ball overlap projection** (coupling step 1 — see ss6.1).
5. **Boundary projection** (position only): for any particle with `drum.sdf_world(x) = d < 0`,
   `x += -d * normal`; marks the particle `touched_wall` for the no-slip blend below.
6. **Velocity reconstruction**: `v = (x - x0) / dt` for every particle (`x0` = pre-predict
   position).
7. **No-slip wall blending**: for particles marked `touched_wall`, `v = (1-beta)*v + beta*v_wall`
   (`v_wall = drum.wall_velocity(x)`), `beta = slurry.wall_no_slip` clamped to `[0,1]`.
8. **Ball viscous drag** (coupling step 2 — see ss6.2).
9. **Buoyancy** (coupling step 3 — see ss6.3).
10. **Implicit Newtonian viscosity** (`mu = slurry.viscosity_pa_s`; skipped entirely, at zero
    iteration cost, when `mu == 0`): density is recomputed once more at the final,
    post-boundary-projection positions, `morris_weights` are built, and `solve_implicit_viscosity`
    solves the diffusion system in place. `mean_shear_rate` is computed from the same final
    positions/densities/neighbour lists.
11. **Fluid speed clamp** (stability backstop): any particle speed above `v_max =
    FLUID_SPEED_SAFETY_FACTOR * (|omega|*radius_m + sqrt(4*|GRAVITY|*radius_m))`,
    `FLUID_SPEED_SAFETY_FACTOR = 5.0`, is rescaled down to `v_max`. Deliberately not a CFL-style
    `h/dt` bound (that would be coupled to `simulation.resolution` and could clip real motion at
    high resolution instead of only blow-ups). Momentum removed here is *not* credited back to any
    ball.
12. **Coupling impulse clamp** (see ss6.4).

### 5.3 Implicit viscosity detail

The Newtonian viscous term `(mu/rho) * laplacian(v)` is discretised with the Morris (1997) SPH
viscous Laplacian (`morris_weights`):

```
c_ij = m_j * (mu_i + mu_j) / (rho_i * rho_j) * |x_ij . grad_W_ij| / (|x_ij|^2 + eta^2)
eta^2 = VISCOSITY_ETA_FACTOR * h^2,  VISCOSITY_ETA_FACTOR = 0.01
```

(single-phase fluid here, so `mu_i = mu_j = mu` and `mu_i + mu_j = 2*mu`). `(Lv)_i = sum_j c_ij *
(v_i - v_j)` is the resulting graph Laplacian (`laplacian_apply`) — symmetric and positive
semi-definite since every `c_ij >= 0`. The backward-Euler system `(I + dt*L) v_new = v` is solved
matrix-free by conjugate gradient (`solve_implicit_viscosity`), operating directly on the `Vec2`
field with the real inner product `<a,b> = sum_i a_i . b_i` (mathematically equivalent to running
scalar CG on the x- and y-components independently, without the bookkeeping of splitting them).
Relative-residual stopping tolerance `VISCOSITY_CG_TOLERANCE = 1e-3` (`|r| <= tol * |b|`), warm-started
from `v` itself, iteration ceiling `VISCOSITY_CG_MAX_ITERS = 50`. `A = I + dt*L` is SPD whenever `dt >
0`, so CG is guaranteed to converge in exact arithmetic; the 50-iteration ceiling is a generous
practical bound, not a tuned-tight one. This replaced an earlier explicit XSPH mixing scheme whose
effective kinematic viscosity saturated around `~h^2/dt` regardless of the input value, so it could
not represent a genuinely low-viscosity (near-water) slurry; the implicit solve is unconditionally
stable for any `viscosity_pa_s >= 0` at the fixed sub-step and uses the value directly in Pa*s.

`mean_shear_rate`: `gamma_dot = sqrt(2 * D:D)`, `D = sym(grad v)`, with `grad_v_i = sum_j (m_j/rho_j)
* (v_j - v_i) (x) grad_W_ij` — the standard SPH velocity-gradient estimator (same neighbour-weighting
form the density constraint and viscosity already use). Surfaced via `FluidStepStats` for
`metrics::Metrics::mean_shear_rate_per_s` and reserved for a future Bingham/Herschel-Bulkley
extension (not implemented; `Rheology::Bingham` exists in `params.rs` as a reserved enum variant but
only `Rheology::Newtonian` is wired into the solver).

`nu = mu/rho` uses each particle's own *measured* SPH density, not `slurry.density_kg_m3` directly
— these agree closely in a well-relaxed interior but diverge near a low-density boundary (free
surface, sparse neighbourhood), where the resulting locally elevated `nu` is itself physically
reasonable (a thin/rarefied film resists shearing more, relative to its own mass, than the bulk).

---

## 6. Two-way ball<->fluid coupling (`coupling.rs`, `pbf.rs` steps 3.5/6.5/6.6/7.5/8)

`coupling::step` orchestrates one shared sub-step: the fluid solves first (`step_coupled`), using
*last* sub-step's ball positions/velocities as a fixed boundary and producing this sub-step's
reaction impulses on the balls; the balls then advance using those impulses
(`DemState::step_with_external_forces`), so their new state becomes the fluid's boundary for the
*next* sub-step (a staggered/semi-implicit scheme).

**Impulses, not forces.** `CouplingImpulses { impulses: Vec<Vec2>, angular_impulses: Vec<f32>,
clamp_hits: u32, fluid_momentum_change: Vec2 }` accumulates a linear + angular impulse per ball. The
exchange is expressed and applied as impulses (`delta_v = impulse * inv_mass`, no extra `* dt`)
rather than as forces re-integrated with another `* dt` on the ball side; the doc comment on
`CouplingImpulses` notes an earlier force-based version double-round-tripped through `dt` and was
visibly unstable (balls flung out of the charge) at this project's default sub-step rate.

### 6.1 Overlap projection (pbf.rs step 3.5)

After the density-constraint solve, any fluid particle within `contact_radius = balls.radius + 0.25
* h` of a ball centre is pushed out along the connecting normal by `push_mag = min(contact_radius -
dist, 0.5 * balls.radius)` — capped at half the ball's radius so a particle deeply embedded in the
contact zone cannot produce an outsized single-substep correction. This position correction becomes
a velocity via the later reconstruction (`Δv = push/dt`), so the fluid particle's own momentum
change is `mass * push / dt`; the ball's Newton's-third-law reaction impulse is the exact opposite,
`-(mass * push) / dt`, applied at the contact point (`lever = n_hat * balls.radius`) for the angular
component.

### 6.2 Viscous no-slip drag (pbf.rs step 6.5)

Fluid within `balls.radius + h` of a ball's surface blends toward the ball's local surface velocity
`v_surf = v_b + omega_b x (x_i - x_b)` over this sub-step:

```
tau = rho_ball * r^2 / (4 * mu)          (Stokes-regime disc relaxation-time closure)
beta = slurry.ball_no_slip * (1 - exp(-dt / tau))
dv = beta * (v_surf - v_fluid)
```

`rho_ball = balls.mass / (pi * balls.radius^2)` (the ball's own effective areal density); `mu =
slurry.viscosity_pa_s` is the *physical* slurry viscosity (not the qualitative XSPH coefficient the
now-removed explicit scheme used). A highly viscous fluid (small `tau`) reaches full no-slip within
one sub-step; an inviscid fluid (`mu -> 0`, `tau -> inf`) applies essentially no drag. `dv` is
already a velocity change, so the ball's reaction impulse is `-(mass * dv)` directly, no `/dt`. This
formula replaced an earlier `ball_no_slip * xsph_c` blend (a qualitative XSPH mixing coefficient,
not tied to the physical viscosity value).

### 6.3 Buoyancy (pbf.rs step 6.6)

Balls are typically sub-resolution relative to the fluid spacing and do not contribute to the
density-constraint sum, so buoyancy is modelled directly rather than emerging from the PBF pressure
field. Each ball samples the local fluid density by reusing the same Poly6 kernel and the fluid's
own step-2 neighbour grid at the ball's position:

```
rho_local = sum_j m_j * W_poly6(|x_ball - x_j|^2, h)   (over nearby fluid particles j)
rho_eff = min(rho_local, rest_density)
impulse_on_ball = (0, -rho_eff * (pi * r^2) * GRAVITY * dt)   (acts through the centroid: no torque)
```

Clamping to `rest_density` avoids over-buoyancy from a locally compacted pocket and tapers smoothly
to zero as a ball nears the free surface (lower sampled density there) instead of an on/off cutoff.
The reaction is applied immediately as a velocity change split across the contributing fluid
particles, weighted by each one's share `w_j / rho_local` of the sampled local density — this must
run after step 6's velocity reconstruction, not before, or the reconstruction would overwrite it.
This mechanism did not exist in the original coupling design (`docs/PLAN.md` ss3.4 describes
buoyancy as emergent from the density-constraint push alone, with an Akinci-style boundary-particle
extension noted as future work); it is now a distinct, directly-modelled step.

### 6.4 Stability clamp (pbf.rs step 8)

The accumulated per-ball impulse this sub-step is clamped to
`impulse_clamp(mass, dt) = F_CLAMP_G_MULTIPLE * mass * 9.81 * dt`, `F_CLAMP_G_MULTIPLE = 20.0`
(angular impulse clamped to that magnitude times `balls.radius`). `pbf.rs`'s doc comment on
`F_CLAMP_G_MULTIPLE` records that this constant used to be `3.0`, tuned down purely to survive an
earlier 3D-sphere-vs-2D-disc mass mismatch between balls and fluid particles; now that both are
unit-depth discs with dimensionally consistent masses, `20x` is a pure numerical-stability backstop
against a transient large overlap rather than a value that shapes normal behaviour. `clamp_hits`
(count of balls clamped this sub-step) and `fluid_momentum_change` (sum of fluid momentum change
from every mechanism above, for a Newton's-third-law cross-check against `sum(impulses)` when no
clamp fired) are both tracked in `CouplingImpulses` purely as diagnostics.

---

## 7. Free surface extraction (`surface.rs`)

Produces closed (rarely open) polylines approximating the slurry free surface for rendering, via:

1. **Splat.** Every fluid particle splats a smooth, compactly-supported cubic-falloff kernel
   (`splat_kernel(r2, radius) = (1 - r2/radius^2)^3` for `r2 < radius^2`, else 0; unnormalized —
   only its threshold-relative shape matters) of radius `SPLAT_RADIUS_FACTOR * fluid.h`
   (`SPLAT_RADIUS_FACTOR = 1.5`) onto a `GRID_SIZE x GRID_SIZE` (`GRID_SIZE = 128`) scalar occupancy
   field spanning `[-drum.radius_m, drum.radius_m]` on each axis (`build_field`).
2. **Mask.** Any grid point outside the wall or inside a lifter bar (per `drum.sdf`, evaluated in
   the drum-local frame) is zeroed — this is what stops the marching-squares contour from bulging
   through the wall in the first place.
3. **Threshold.** Marching squares runs at `THRESHOLD_FRACTION * phi_full`
   (`THRESHOLD_FRACTION = 0.5`), where `phi_full` is computed **analytically**
   (`analytic_phi_full`), not read off the grid's measured maximum:
   `phi_full = number_density * integral(K dA) = (rest_density / particle_mass) * (pi *
   splat_radius^2 / 4)`. `number_density = 1/dx^2` is exact at the PBF density constraint's
   equilibrium because `particle_mass = rest_density * dx^2` at seeding; the kernel's integral has
   the closed form `pi * splat_radius^2 / 4` (substitute `u = r^2/splat_radius^2`). By symmetry, the
   field value exactly at a flat bulk/free-surface boundary is exactly half the deep-bulk value, so
   the `0.5` contour lands precisely on the physical surface. Using the grid's own measured maximum
   instead is unsafe: at high rotation speed the slurry can centrifuge into a thin wall-hugging film
   whose field values are depressed everywhere, while a transient over-compacted PBF pocket
   elsewhere can inflate the measured maximum well past the true bulk value and starve the film's
   threshold, making the real free surface invisible to marching squares. If the analytic threshold
   yields no contour at all (a degenerate case), extraction retries once with the grid's measured
   peak as the threshold reference instead, as a graceful fallback (`extract_surface`).
4. **Contour.** Standard marching squares (`cell_edge_pairs`, with the two ambiguous 4-corner cases
   resolved by a fixed, non-adaptive choice) produces segments per cell, which are chained into
   polylines by exact `EdgeId` matching (`extend_chain`).
5. **Smooth & clamp.** One pass of Chaikin corner-cutting (`chaikin_smooth_closed`) softens the
   blocky contour; every resulting point still outside the wall (up to one grid cell of overshoot
   from linear interpolation/smoothing) is projected back onto the wall along its SDF normal
   (`project_inside_wall`).

`flatten_polygons` serialises the result as `[n_polys, len_0, x, y, ..., len_1, x, y, ...]` for the
wasm/frontend boundary.

---

## 8. Metrics & grinding instrumentation (`metrics.rs`, `lib.rs`)

`metrics::compute(balls, fluid, drum, drum_angle)` is a pure function of the current populations and
derives:

- **Toe/shoulder angles** (`charge_toe_shoulder`): balls within `WALL_MARGIN_BALL_RADII (2.5) *
  balls.radius` of the wall approximate the charge's outer layer; if there are at least
  `MIN_WALL_BALLS (6)` such balls and their angular gap exceeds `CENTRIFUGE_GAP_THRESHOLD_RAD (0.3
  * pi)` (i.e. not centrifuged), the two ends of the largest angular gap in that wall-adjacent
  population are the toe/shoulder (assignment to toe vs. shoulder depends on `drum.omega`'s sign).
  `None`/`None` otherwise.
- **Charge centroid** (`charge_centroid`): plain unweighted mean of all ball positions (all balls
  share one mass, so this is the correct centroid).
- **Slurry pool angular extent** (`slurry_pool_angular_extent`): fluid particles within
  `POOL_WALL_MARGIN_H (1.5) * fluid.h` of the wall (per `drum.sdf_world`, so lifters count),
  clustered the same wraparound-aware way as toe/shoulder.
- **Free-surface line fit** (`fluid_free_surface_line`): fluid particles *not* near the wall (the
  pool's own exposed top), fit by the standard 2x2 covariance-matrix principal-eigenvector method;
  returns `(angle_rad mod pi, offset_m)` where `offset_m` is the signed perpendicular distance of
  the fitted line from the drum centre. `None` if fewer than 3 such particles.
- **Pool depth at bottom** (`slurry_pool_depth_at_bottom`): among fluid particles within
  `POOL_DEPTH_BAND_RAD (+/- 5 deg)` of straight-down, `drum.radius_m - min(distance_from_centre)`.
- **Lacey mixing index** (`mixing_index`): particles binned into a `MIXING_GRID_SIZE (16) x 16`
  grid; `M = (S0^2 - S^2) / (S0^2 - Sr^2)` where `S^2` is the variance of per-cell mean dye value
  (unweighted across occupied cells), `S0^2 = p*(1-p)` is the fully-segregated variance (`p` =
  overall mean dye fraction), `Sr^2 = S0^2 / n_bar` is the fully-mixed variance using the average
  occupied-cell particle count as the per-cell sample size. Clamped to `[0,1]`; degenerate cases
  (near-zero denominator) return `1.0`.
- **Total ball kinetic energy** (`total_kinetic_energy_j`): `sum(0.5 * mass * |v|^2)`.
- **Max ball-ball overlap fraction** (`max_ball_overlap_fraction`): `(2r - dist)/r` over
  broad-phase candidate pairs (same `UniformGrid` cell size as the DEM contact solve).
- **Max/mean fluid density error fraction** (`max_fluid_density_error_fraction`,
  `mean_fluid_density_error_fraction`): `|rho_i - rest_density| / rest_density`, reusing
  `FluidParticles::densities()`.

**Grinding/solver diagnostics are not derivable from `(balls, fluid, drum)` alone** — `compute`
zero/empty-defaults these fields; `Simulation::metrics()` fills them in from its own accumulated
state after calling `compute`:

- `power_draw_w`, `torque_nm`, `collision_rate_per_s`, `dissipated_power_w`,
  `impact_energy_histogram` (`ImpactEnergyHistogram { bin_edges_j, counts_per_s }`): EMA-smoothed
  (time constant `GRINDING_STATS_EMA_TAU_S = 1.0` s, i.e. `alpha = 1 - exp(-sub_dt/tau)` applied
  every sub-step in `Simulation::update_grinding_stats`) from `DemStepStats`. `power_draw_w =
  wall_work_j / sub_dt`; `torque_nm = power_draw_w / omega` (0 when `|omega| <= 1e-6`, since torque
  is undefined rather than infinite at zero rotation); `collision_rate_per_s =
  collision_count/sub_dt`; `dissipated_power_w = dissipated_energy_j/sub_dt`; each histogram bin's
  EMA is `count/sub_dt`.
- `coupling_clamp_hits`: sum of `CouplingImpulses::clamp_hits` over every sub-step of the most
  recent `Simulation::step` call — reset each call, **not** EMA-smoothed (meant to read as "did
  this happen just now", not a smoothed rate).
- `effective_ball_diameter_m`, `simulated_ball_count`, `true_ball_count`, `coarse_graining_factor`:
  from `Params::effective_media`/`Params::true_ball_count` (ss3).
- `max_substep_displacement_over_diameter`: the largest value of
  `DemStepStats::max_substep_displacement_over_diameter` seen across every sub-step of the most
  recent `Simulation::step` call — reset each call, **not** EMA-smoothed (a tunnelling-risk spike
  is exactly the kind of transient an average would hide). See ss9.
- `mean_shear_rate_per_s`, `viscosity_solver_iterations`: from the most recently completed
  sub-step's `FluidStepStats`.

---

## 9. Known limitations

- **Tunnelling / continuous-collision guard is diagnostic only.** `DemStepStats::
  max_substep_displacement_over_diameter = |v| * dt / (2*radius)` (the fraction of a ball's own
  diameter it moved in one sub-step) is computed every sub-step and surfaced through
  `Simulation::max_substep_displacement_over_diameter` and `metrics::Metrics::
  max_substep_displacement_over_diameter`, but **there is no actual swept/continuous-collision
  correction implemented**, and the DEM broad-phase cell size (`2 * r * 1.05`, `dem.rs` step 2)
  includes no relative-displacement margin. A fast ball's motion within a single discrete sub-step
  can therefore, in principle, skip past a collision the broad-phase would otherwise have caught.
  This is expected to read non-trivially high (observed ~0.85 at default parameters) during normal
  cascading, not only in pathological cases — it is a live risk indicator, not evidence of a bug,
  and not yet mitigated by any correction in the solver.
- **`s_corr` artificial pressure is implemented but disabled** (`S_CORR_K = 0.0` in `pbf.rs`) — see
  ss5.2. Revisit only if visual clustering artefacts are observed; no working non-zero `k` was found
  for this project's SI parameter scale.
- **Bingham/Herschel-Bulkley rheology is not implemented.** `params::Rheology::Bingham` exists as a
  reserved enum variant (with `yield_stress_pa` already a validated `SlurryParams` field) but only
  `Rheology::Newtonian` is wired into `pbf.rs`'s viscosity solve; `mean_shear_rate` is computed and
  surfaced specifically in anticipation of this future extension.
- **No ball size distribution.** `Balls` (and `EffectiveMedia`) model one shared radius/mass/inertia
  for the whole population; per-particle size variation is not implemented.
- **Buoyancy has no boundary-particle contribution to the fluid's own density field.** Balls do not
  contribute to the PBF density-constraint sum (`step_coupled` step 3), so a ball never "occupies"
  space from the fluid's perspective the way an Akinci-style boundary-particle scheme would; the
  directly-modelled buoyancy impulse (ss6.3) is a substitute for, not an emergent consequence of,
  incompressibility.
- **2D cross-section only.** Every value in this crate is per metre of unit mill depth; no axial
  (out-of-plane) transport, end-wall effects, or 3D particle shape are modelled. This is a
  deliberate, project-wide scope decision (see `docs/PLAN.md`), not a defect, but it means results
  are qualitative relative to a real 3D mill.
