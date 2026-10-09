# Scale-up study: relative grinding proxies

Purpose: help carry a small-scale (lab) grinding experiment over to a larger mill by comparing
**relative** DEM quantities between conditions. Grinding itself is not simulated and absolute
values are not trusted; what is meaningful is a *ratio* between two conditions computed with the
same code and the same contact coefficients. Code: `crates/mill-core/src/scaleup.rs` (pure Rust,
CPU only, runs in WASM and natively with identical results), UI: the "Scale-up study" panel,
validation: `cargo run --release -p mill-core --example scaleup_probe` (results in
`docs/VERIFICATION.md`, section "Scale-up runner").

## What is measured

Every impact the ball solver sees (`dem::ImpactRecord`: normal approach speed, tangential slip speed
at the contact point with spin included, restitution used) is turned into a **stress energy**
`SE = 1/2 m_red (v_n^2 + v_t^2)` using the **3D sphere mass** `m = rho (4/3) pi r^3` (ball-ball
`m_red = m/2`, ball-wall `m_red = m`). The 2D unit-depth disc mass used inside the solver is *not*
used: it overstates `SE` by a factor that depends on the ball size (about 750x for a 2 mm ball) and
would corrupt comparisons between different ball sizes.

Impacts are classified, for ball-ball and ball-wall separately, as **drop** (`v_n >= v_t`) or
**shear** (`v_t > v_n`): four classes (`CLASS_NAMES`). The two types scale differently with drum
size `D` and ball size `d` (drop: `SE/(m g d) ~ D/d`; shear: `~ d/D`), so they are extrapolated
separately.

Outputs (`Summary`, per run; ensemble mean and standard error over seeds):

| Quantity | Meaning |
|---|---|
| `wall_power_per_ball_w`, `dissipated_power_per_ball_w` | Drum power and DEM contact dissipation per ball. Convention: a one-ball-thick slab (depth = ball diameter), i.e. the 2D per-metre power times `d` divided by the ball count. |
| `specific_power_w_per_kg` | `wall_power_per_ball_w / m` (3D ball mass). Scales as `sqrt(L)` between similar mills. |
| `impacts_per_ball_per_s.<class>` | Impact frequency per ball and class. |
| `mean_se_scaled_drop`, `mean_se_scaled_shear`, `stress_intensity_scaled` | Mean `SE / (m g d)`: dimensionless, comparable across sizes (Kwade-type stress intensity). |
| `effective_power_per_ball_w.thr*`, `effective_power_fraction.thr*` | Rate of stress energy carried by impacts above a threshold (default `[0.03, 0.1, 0.3, 1, 3, 10] m g d`, or absolute thresholds in J from a material test), and its fraction of wall power. |
| Time series (`Interval`) | Per recording interval (default 1/8 revolution): powers, impact rates, mean `SE`, **cumulative** wall energy and cumulative effective energy per threshold, and the `SE` histogram per class. `ScaleupResult::time_to_cumulative_effective_energy` gives the time to reach a given cumulative effective energy; the ratio of two conditions' times is the relative-grinding-time proxy. |
| `power_drift_rel_per_rev`, `power_rel_std` | Steady-state indicators for the recorded window. |

Time-series export: CSV (`ScaleupResult::intervals_csv`, scalars) and JSON (everything, including
histograms), per seed.

## Numerical settings are ball-relative

So that conditions of different size are computed equally accurately and a geometrically similar
mill reproduces the same dimensionless results:

- sub-step `dt = dt_ratio sqrt(d/g)` (not a fixed frequency);
- the fresh-impact speed threshold of the restitution pass and the friction speed scale are
  `0.1428 sqrt(g d)` (`media.restitution_threshold_froude`; identical to the historical 0.02 m/s
  at the shipped 2 mm ball; the interactive app keeps the old absolute value, bit-identical);
- no coarse-graining (`CoarseGrainingMode::Off`): the true ball count is simulated (limit 50000);
- impact recording uses the same `0.1428 sqrt(g d)` threshold.

## Wet operation (what is and is not modelled)

Particle slurry is *not* part of this runner. Wet impact damping enters through the
**lubrication-cushion** law `e = e_dry exp(-35 / St)`, `St = rho_ball v_n d / (9 mu)` (Legendre et
al. 2006; `dem::wet_cushion_restitution`, `media.wet_cushion_viscosity_pa_s`). Because `v_n`
grows with `sqrt(g D)`, `St` changes with scale, which is the main size-dependent wet effect.
Not modelled: viscous drag on the charge motion, buoyancy, solids-concentration-dependent
rheology, tangential lubrication, slurry cushioning of ball-liner contact beyond the same law.
The constant 35 is a literature value for sphere-wall impacts and has **not** been calibrated
against this project's resolved-flow solver (`mac/`); treat wet ratios as qualitative.

## Extrapolating to a larger mill without simulating it

Real production mills hold 10^5-10^7 balls, far beyond a browser. The runner therefore offers a
**D/d sweep** (`sweep_dd`): `d`, `% critical` and all dimensionless settings fixed, `D` varied over
a range that is affordable (default `D/d` = 20-90, up to ~2000 balls), then a log-log power-law fit
per quantity (`PowerLawFit`, with the exponent expected from the drop/shear scaling shown next to
the fitted one). A fit with R^2 < 0.95, or extrapolation more than ~1.5x outside the swept range,
must be treated as unreliable. `matching_ball_diameter(d_lab, D_lab, D_target)` gives the ball size
that keeps the drop-type stress energy (`d ~ D^(-1/3)`) as a starting point for a candidate search.

## Sensitivity to the contact coefficients

Restitution and friction are not calibrated. `sensitivity()` re-runs a condition with each of the
four coefficients scaled by 1 +- rel (default 20 %). A *ratio* between two conditions is only
trustworthy if it survives this perturbation; the panel highlights changes above 10 %.

## Known limits (read before using a number)

- 2D cross-section: no end walls or axial transport; coordination number and impact frequency
  differ from 3D. Only relative comparison between 2D runs is meaningful.
- Rigid (XPBD) contacts: no contact duration or Hertz stiffness; fine for hard media, but the
  impact spectrum is only as converged as the numerical settings (see VERIFICATION.md).
- No powder/ore, no size distribution, no breakage; the energy threshold is an input.
- Monodisperse media; lifters are supported, liner wear is not.
- No comparison with measurements has been made (not part of this work).
