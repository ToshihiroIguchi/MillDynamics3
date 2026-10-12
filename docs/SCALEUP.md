# Grinding proxies of the running simulation (for relative scale-up studies)

Purpose: help carry a small-scale grinding experiment over to a larger mill by comparing
**relative** quantities between runs of the *main simulation*. Grinding itself is not simulated and
absolute values are not trusted; what is meaningful is a *ratio* between two conditions run with the
same code and the same numerical settings. Code: `crates/mill-core/src/stress.rs`; the values are
`Metrics.stress` and the "Grinding proxy (3D)" rows of the metrics panel (docs/METRICS.md), with a
"Set reference" control that shows each row as a ratio to a stored run.

History: an earlier separate "scale-up runner" with its own settings (no slurry, finer steps, its
own panel) was built and then removed, because numbers that differ from the main screen are
confusing. Its validation results (docs/VERIFICATION.md, "Scale-up runner") still describe how the
recorded statistics behave.

## What is measured

Every impact the ball solver sees (`dem::ImpactRecord`: normal approach speed, tangential slip speed
at the contact point with spin included, restitution used) is turned into a **stress energy**
`SE = 1/2 m_red (v_n^2 + v_t^2)` with the **3D sphere mass** `m = rho (4/3) pi r^3` (ball-ball
`m_red = m/2`, ball-wall `m_red = m`). The 2D unit-depth disc mass used inside the solver is *not*
used: it overstates `SE` by a ball-size-dependent factor (about 750x for a 2 mm ball). Only impacts
faster than `0.1428 sqrt(g d)` (0.02 m/s at 2 mm) are recorded; the solver's own restitution logic is
unchanged.

Impacts are classified, for ball-ball and ball-wall separately, as **drop** (`v_n >= v_t`) or
**shear** (`v_t > v_n`): four classes. All values are EMA-smoothed (1 s) except the cumulative ones.

| Quantity | Meaning |
|---|---|
| `specific_power_w_per_kg`, `wall_power_per_ball_w` | Drum power (incl. slurry drag on the wall) per ball, one-ball-thick slab (depth = ball diameter), and per 3D ball mass. Scales as `sqrt(L)` between similar mills. |
| `impacts_per_ball_per_s` | Impact frequency per ball, per class and total. |
| `mean_se_scaled_drop/shear`, `stress_intensity_scaled` | Mean `SE / (m g d)`: dimensionless, comparable across sizes (Kwade-type stress intensity). |
| `effective_power_per_ball_w` | Rate of stress energy carried by impacts above each threshold (`[0.03, 0.1, 0.3, 1, 3, 10] m g d`). |
| `cumulative_*` | Wall work and above-threshold stress energy per ball since the simulation started (monotone; take differences between two times). The time two conditions need to reach the same cumulative effective energy is the relative-grinding-time proxy. |
| `impact_energy_histogram` | Pooled stress-energy histogram (28 log bins, joules, 3D mass). |

## How to use it for scale-up

1. Same %Nc, fill J, ball material; change the quantity under study. Keep **every** numerical
   setting identical between compared runs (sub-steps, solver iterations, coarse-graining), see below.
2. Keep coarse-graining off (`k = 1`): the rows are meaningless when `k > 1` (the simulated ball is
   not the real ball). This limits the size that can be simulated (a few thousand balls).
3. Let the run settle, then use "Repeat runs (seeds)" in the Grinding proxy (3D) group: it runs the
   same simulation with the same applied settings for seeds 1..N over whole revolutions and gives mean
   +- standard error per row, so the ratio to the lab condition ("Use as reference" on its batch) comes
   with a propagated error (ratios above 10 % relative error are marked). The report window with
   "Set reference" remains a single-run alternative.
4. Compare dimensionless values (`SE/(m g d)`) across ball sizes, per-ball frequencies and specific
   power across drum sizes.
5. Use the ratios only if they exceed the uncertainties below.

## What is missing before a number can be trusted for scale-up

Measured results are in docs/VERIFICATION.md ("Scale-up runner"); items marked *(not done)* are
open work.

1. **Numerical settings.** The app default (2 contact iterations) understates the stress energies by
   about 15 % against 12 iterations (12 vs 20: 1-3 %); sub-steps 16 / 32 / 64 move the stress intensity by
   about 4 % and the specific power by < 4 %. Use the "Impact statistics" accuracy preset (32 sub-steps,
   12 iterations, coarse-graining off; about 4x slower than the default) and compare only runs with
   identical settings (the UI warns when the settings differ from the stored reference).
   Measurements: docs/VERIFICATION.md, "Impact statistics of the main simulation".
2. **Impact counts are not converged in the sub-step** (+4 % then +7 % per doubling even with
   re-contacts of the same pair within `0.25 sqrt(d/g)` counted once; the debounce window itself is a
   modelling choice that scales the absolute count). Treat impact frequency as comparable only at
   identical settings, and prefer the mean stress energy and the specific power, which are more stable
   *(a time-step-independent impact definition is not done: counting only contact onsets was tried and is worse, +25 % per doubling, docs/VERIFICATION.md)*.
3. **Statistics.** The charge is chaotic and the EMA is a 1 s average. "Repeat runs (seeds)" in the
   metrics panel averages whole revolutions over several seeds of the same `Simulation` with the same
   applied settings and reports standard errors and propagated ratio errors (docs/METRICS.md,
   "Repeat runs (seeds)"). It is slow with slurry on, and the number of seeds and revolutions is the
   user's choice; the warm-up needed for steady state is flagged through the power-drift indicator.
4. **Reaching a large mill.** Do not extrapolate stress energies in D/d. With true-size balls
   (k = 1) at 59 % Nc the stress intensity is 0.22 for D/d 20-45, then rises: 0.26 (60), 0.28 (100),
   0.36 (150), 0.62 (250, 1 seed); the impact rate per ball rises from about 100 to 200 (D/d 45 -> 150) and
   specific power grows more slowly than `sqrt(D)` at large D/d (docs/VERIFICATION.md; an earlier "flat up to
   D/d 350" result was a coarse-graining artefact, see the CORRECTION there). The simple free-fall /
   shear-rate scaling laws are not supported. For a 1 m drum with 2 mm media (D/d = 500, about 61 500
   discs, above the 50 000-ball limit) no number is available: the trend is steep and unexplained, and
   direct simulation at D/d 500 is out of reach. Coarse-graining cannot substitute: even k = 1.2 changed the
   intensity by 14 % and the shear energy by 25 %. What would be needed: an explanation of why the
   stress energies grow with D/d (e.g. a change of the charge's motion regime with the relative charge
   depth), checked with the existing `Impact statistics` preset and the seed batch, and runs at D/d
   300-400 with steady-state checks.
5. **Contact coefficients are not calibrated.** +-20 % on ball-ball friction changes the shear stress
   energy by +17 % / -8 %; +-20 % on ball-ball restitution changes the impact rate by -11 % / +14 %.
   Ratios within about 10-15 % are inside this uncertainty. Check a conclusion by re-running with
   the coefficients changed.
6. **Model scope.** 2D cross-section (no end walls, different coordination number: relative only),
   monodisperse media, no powder/ore, no breakage, no size distribution; the energy threshold for
   "effective" impacts is an input that has to come from the material, not from this model.
7. **Slurry.** The particle slurry is not converged at water-like viscosity (+-12 % in power, see
   docs/VERIFICATION.md), has no solids-concentration dependence, and the wet impact-cushion effect
   is only what the fluid solver produces. Compare wet runs only with the same fluid; a dry run is
   the cleaner reference for ball dynamics.
