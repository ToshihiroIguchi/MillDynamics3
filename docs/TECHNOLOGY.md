# MillDynamics3 - Technology and Innovation Overview

MillDynamics3 simulates the cross-section of a tumbling ball mill in real time. Grinding media are rigid discs advanced by a position-based (XPBD) contact solver, and the slurry is a particle fluid advanced by Position Based Fluids (PBF) with implicit viscosity. The two are coupled two-way every sub-step by impulse exchange, including a joint implicit no-slip solve and a sub-grid squeeze-film model. This document writes out every model that the shipped core (`crates/mill-core`) actually runs, with its equations, discretisation, constants and source location, followed by the research-track solvers that are in the repository but not used by the simulation. The source code is the ground truth; where an older design note disagrees, the code is described.

**Contents**

- 0. Summary
- 1. What is simulated
- 2. Physics and numerical models
  - 2.1 Drum geometry
  - 2.2 Grinding media: rigid discs, XPBD contact
  - 2.3 Coarse-graining and automatic accuracy scaling
  - 2.4 Slurry: Position Based Fluids
  - 2.5 Ball-fluid coupling
  - 2.6 Free surface for rendering
  - 2.7 Metrics and the 2D-to-3D grinding proxy
  - 2.8 Research tracks not used by the shipped simulation
  - 2.9 Determinism
- 3. Software technology
- 4. Verification and quality assurance
- 5. What is innovative, and what is not
- 6. References

Status: written against commit `2d35f15` (2026-10-18). Source locations are `path:line`; line numbers drift as the code changes, so treat the symbol name as the stable anchor. Where `docs/PHYSICS.md` and the source disagree, this document follows the source. Literature references are numbered in order of first citation and listed in Section 6; a limited literature survey underlies Section 5 (see there for its scope).

---

## 0. Summary

| Layer | Technique | Key source | Main reference |
|---|---|---|---|
| Drum | Signed-distance field of a circle plus trapezoidal lifters (exact convex-polygon SDF), central-difference normal, rigid wall velocity $\omega\times r$, swept time-of-impact guard | `crates/mill-core/src/geometry.rs:37`, `crates/mill-core/src/geometry.rs:146` | none (standard) |
| Media contact | XPBD rigid discs, zero compliance, bounded overlap recovery, Coulomb friction as clamped position correction, bidirectional restitution, rolling resistance | `crates/mill-core/src/dem.rs:530` | [1, 2]; contrast [3]; mill DEM context [4, 5] |
| Media scaling | Mass-preserving coarse-graining $k=\sqrt{N_\text{true}/N_\max}$; automatic fluid resolution and sub-step count | `crates/mill-core/src/params.rs:582`, `crates/mill-core/src/params.rs:648`, `crates/mill-core/src/params.rs:690` | context [6, 7, 8] |
| Slurry pressure | PBF: one-sided density constraint, Jacobi projection with adaptive iteration count, 2D Poly6 and Spiky kernels | `crates/mill-core/src/pbf.rs:1241` | [9, 10, 11] |
| Slurry surface forces | Akinci-style cohesion and boundary adhesion with 2D re-normalised kernels | `crates/mill-core/src/pbf.rs:1455`, `crates/mill-core/src/pbf.rs:1493` | [12] |
| Slurry viscosity | Implicit backward-Euler, Morris weights in central-force form, Jacobi-preconditioned CG; Papanastasiou-regularised Bingham option | `crates/mill-core/src/pbf.rs:2122`, `crates/mill-core/src/pbf.rs:640` | [13] |
| Ball-fluid coupling | Impulse exchange; joint implicit solve of fluid and ball velocities with a ball-surface no-slip boundary; form-drag relaxation; rotating-frame buoyancy; impulse clamp | `crates/mill-core/src/pbf.rs:640`, `crates/mill-core/src/pbf.rs:1702`, `crates/mill-core/src/coupling.rs:71` | project-specific; [13] for the weights |
| Near-contact lubrication | Sub-grid 2D Reynolds squeeze-film links between ball pairs and ball-wall, blended with the resolved fluid | `crates/mill-core/src/pbf.rs:543`, `crates/mill-core/src/mac/lubrication.rs:53` | classical [14, 15]; corrections in resolved solvers [16, 17, 18]; in DEM [19, 20] |
| Rendering | Particle splat to a $128^2$ field, analytic threshold, marching squares, Chaikin smoothing | `crates/mill-core/src/surface.rs:233` | [21, 22] |
| Outputs | Toe/shoulder, Lacey mixing index, wall-work power and torque, energy budget, 3D-converted stress-energy statistics | `crates/mill-core/src/metrics.rs:195`, `crates/mill-core/src/stress.rs:68` | [23, 24, 25] |
| Determinism | Counting-sort CSR neighbour grid, fixed iteration order, xorshift-multiply PRNG | `crates/mill-core/src/grid.rs:101`, `crates/mill-core/src/rng.rs:7` | [26] |
| Research tracks (not shipped) | DFSPH with tabulated smooth-solid kernel integrals (`fluid/`); MAC cut-cell Navier-Stokes with variational viscosity, multigrid-CG, WENO5/TVD-RK3 level set, Newton body coupling, Herschel-Bulkley by ALG2 (`mac/`) | `crates/mill-core/src/fluid/mod.rs:1`, `crates/mill-core/src/mac/mod.rs:1` | [27, 28]; [29, 30, 31, 32, 33, 34, 35] |

What is new in this combination, and what is not, is discussed in Section 5.

---

## 1. What is simulated

**System.** A horizontal tumbling mill seen as a 2D cross-section perpendicular to its axis. A circular drum of radius $R$, optionally with $n_L$ trapezoidal lifter bars (default: none, `crates/mill-core/src/params.rs:417`), rotates at a prescribed angular speed $\omega$ about the world origin. Inside are $N$ rigid circular grinding balls and a viscous liquid slurry. Gravity acts along $-y$ with $g=9.81\ \mathrm{m/s^2}$ (`GRAVITY = -9.81`, `crates/mill-core/src/dem.rs:20`, `crates/mill-core/src/pbf.rs:45`).

**Unit-depth convention.** Every body is a slice of unit axial depth (1 m). A ball has areal mass $m=\rho\pi r^2$ (kg per metre of mill length), a fluid particle has areal mass $\rho_0\,\Delta A$. All masses, energies and powers reported by the core are therefore per metre of mill length. A real mill's value is obtained by multiplying by its axial length. All quantities are SI inside the core; the UI converts mm and rpm.

**Determinism.** A run is a pure function of `Params` and `seed` (Section 2.9).

**Shipped solvers.** `Simulation` (`crates/mill-core/src/lib.rs:82`) runs exactly two solvers, the DEM/XPBD ball solver (`dem.rs`) and the PBF slurry solver (`pbf.rs`), through `coupling::step_recorded` (called at `crates/mill-core/src/lib.rs:479`, defined at `crates/mill-core/src/coupling.rs:71`). The `fluid/` (DFSPH) and `mac/` (grid Navier-Stokes) modules are compiled into the crate (`crates/mill-core/src/lib.rs:22`, `crates/mill-core/src/lib.rs:25`) but are not referenced by `Simulation`. One function of `mac/` is used by the shipped code: the squeeze-film model `mac::lubrication` is imported by `squeeze_links` (`crates/mill-core/src/pbf.rs:551`). Everything else in `fluid/` and `mac/` is verification and calibration infrastructure (Section 2.8).

**Sub-step size.** One rendered frame at 60 Hz is divided into $n_s$ fixed sub-steps,

$$
\Delta t = \frac{1}{60\,n_s},\qquad n_s=\texttt{effective\_substeps}(\text{params}),
$$

(`Simulation::fixed_sub_dt`, `crates/mill-core/src/lib.rs:426`). The requested value is `simulation.substeps`, default 8 (`crates/mill-core/src/params.rs:488`); the effective value is raised automatically for small balls (Section 2.3) and is 16 at the shipped 2 mm / 63 mm defaults, i.e. $\Delta t = 1/960$ s. The constant `FIXED_DT = 1/240` s (`crates/mill-core/src/lib.rs:41`) corresponds to $n_s=4$ and is used only by the benchmark harness, not by the solver. `Simulation::step(dt)` splits an arbitrary `dt` into $n_s$ equal sub-steps (`crates/mill-core/src/lib.rs:456`); `Simulation::step_fixed` advances one sub-step (`crates/mill-core/src/lib.rs:438`) for the web worker's wall-clock accumulator loop.

**Per-sub-step pipeline** (`advance_one_sub_step`, `crates/mill-core/src/lib.rs:470`; the fluid advances first, using the balls' old state as a moving boundary, then the balls advance using the fluid's reaction impulses, so their new state is the boundary of the next fluid step - a staggered, semi-implicit scheme):

1. **Fluid, `FluidParticles::step_coupled`** (`crates/mill-core/src/pbf.rs:1173`; called at `crates/mill-core/src/coupling.rs:83`):
   1. Sanitise non-finite state (`pbf.rs:1189`).
   2. Predict: $v_y \mathrel{+}= -g\,\Delta t$, $x \mathrel{+}= v\,\Delta t$ (`pbf.rs:1223`).
   3. Neighbour lists on a counting-sort grid of cell size $h$ (`pbf.rs:1236`).
   4. Density-constraint (PBF pressure) projection, adaptive iteration count (`pbf.rs:1241`).
   5. Ball overlap projection with bounded momentum exchange (step 3.5, `pbf.rs:1337`).
   6. Cohesion and adhesion velocity kicks, computed here and applied after step 9 (step 3.6, `pbf.rs:1446`).
   7. Drum-wall boundary projection at the end-of-step drum angle (step 4, `pbf.rs:1633`).
   8. Velocity reconstruction $v=(x-x_0)/\Delta t$, then add back the deferred cohesion/adhesion kicks and remove the geometric excess of the overlap push (step 5, `pbf.rs:1654`).
   9. Wall no-slip blending (step 6, `pbf.rs:1675`).
   10. Ball form-drag relaxation (step 6.5, `pbf.rs:1702`).
   11. Rotating-frame buoyancy (step 6.6, `pbf.rs:2012`).
   12. Implicit viscosity, solved jointly with the ball velocities when balls are near, plus squeeze-film links (step 7, `pbf.rs:2122`).
   13. Fluid speed clamp (step 7.5, `pbf.rs:2303`) and per-ball impulse clamp (step 8, `pbf.rs:2317`).
2. **Balls, `DemState::step_recorded`** (`crates/mill-core/src/dem.rs:530`; called at `crates/mill-core/src/coupling.rs:85`), applying the fluid impulses in the predict step:
   1. Predict with gravity and fluid impulses (`dem.rs:565`).
   2. Broad phase on a counting-sort grid (`dem.rs:597`).
   3. Swept time-of-impact clamp against wall and lifters (step 2.5, `dem.rs:635`).
   4. Non-penetration projection, Gauss-Seidel (`dem.rs:710`).
   5. Velocity reconstruction (`dem.rs:752`).
   6. Friction (`dem.rs:765`), restitution (`dem.rs:892`), rolling resistance (`dem.rs:1021`).
3. **Bookkeeping** (`crates/mill-core/src/lib.rs:491`): accumulate the energy budget and EMA-smoothed power, torque and impact statistics (`update_grinding_stats`, `lib.rs:516`), then rotate the drum, $\theta_D \leftarrow (\theta_D+\omega\Delta t) \bmod 2\pi$ (`lib.rs:507`).

Both solvers evaluate the wall at the drum angle at the end of the sub-step, $\theta_D+\omega\Delta t$ (`drum_angle_next`, `dem.rs:674`, `pbf.rs:1234`), because positions have already been predicted to $t+\Delta t$.

---

## 2. Physics and numerical models

### 2.1 Drum geometry

Source: `crates/mill-core/src/geometry.rs`.

The drum frame rotates with the drum. A world point $p$ is rotated into the drum frame by $-\theta_D$ before evaluating the signed distance (`Drum::sdf_world`, `geometry.rs:87`). Convention: $d>0$ in free space, $d<0$ inside solid.

**Wall.** $d_\text{wall}(p)=R-|p|$ (`Drum::sdf`, `geometry.rs:37`).

**Lifters.** For $n_L>0$ bars at angles $\theta_k=\varphi+2\pi k/n_L$ ($\varphi$ = phase), a point has radial and tangential coordinates in bar $k$'s frame

$$
p_r = p\cdot(\cos\theta_k,\sin\theta_k),\qquad p_t = p\cdot(-\sin\theta_k,\cos\theta_k).
$$

The bar cross-section is the convex quadrilateral with counter-clockwise vertices $(R,-b/2)$, $(R,b/2)$, $(R-H,t/2)$, $(R-H,-t/2)$ in $(p_r,p_t)$, where $H$ is the lifter height, $b$ the base width and $t$ the top width (`lifter_cross_section_sdf`, `geometry.rs:211`). The signed distance to a convex polygon with vertices $a_i$ is exact (`convex_polygon_sdf`, `geometry.rs:247`):

$$
\phi_\text{poly}(p) = \sigma\,\min_i \operatorname{dist}\bigl(p,\ [a_i,a_{i+1}]\bigr),\qquad \sigma=\begin{cases}+1 & p \text{ left of every edge (inside)}\\ -1 & \text{otherwise,}\end{cases}
$$

with the distance to each edge taken to the closed segment (the projection parameter is clamped to $[0,1]$), so vertex-nearest exterior points have the true Euclidean distance. The drum SDF is

$$
d(p)=\min\Bigl(R-|p|,\ -\max_k \phi_{\text{poly},k}(p)\Bigr)
$$

(`Drum::sdf_lifters`, `geometry.rs:50`). With $n_L=0$ it reduces exactly to the circle.

**Normal.** The unit normal (pointing into free space) is the central-difference gradient of $d$ with step $\epsilon=10^{-4}$ m,

$$
\hat n_\text{local}=\frac{\nabla d}{|\nabla d|},\quad \partial_x d\approx\frac{d(p+\epsilon\hat x)-d(p-\epsilon\hat x)}{2\epsilon},
$$

falling back to $-\hat p$ if the gradient is degenerate, then rotated back to the world frame (`geometry.rs:94`-`geometry.rs:107`). A non-finite input returns $(0,\mathbf 0)$ so a caller's correction $p\mathrel{+}=-d\hat n$ is a no-op rather than a NaN (`geometry.rs:88`).

**Wall velocity.** The rigid-body velocity of the wall material point currently at $p$ is $v_\text{wall}(p)=\omega\,(-p_y,p_x)$ (`Drum::wall_velocity`, `geometry.rs:112`).

**Swept time of impact.** `Drum::toi_swept` (`geometry.rs:146`) returns the first fraction $\tau\in[0,1]$ of a sub-step at which a disc of radius $r$, whose centre moves on a straight world-space line $x_0\to x_1$ while the drum turns from $\theta_0$ to $\theta_1$, touches the solid. It sphere-traces the static drum-frame SDF along the path expressed in the instantaneous drum frame, advancing

$$
\Delta\tau = \max\!\Bigl(\frac{d(\tau)-r}{|x_1-x_0| + |\theta_1-\theta_0|\,(|x_0|+|x_1-x_0|)},\ 10^{-6}\Bigr),
$$

where the denominator is a rigorous upper bound on the speed of the point in the drum frame, so the step never overshoots (the SDF is 1-Lipschitz). It stops at the first $\tau$ with $d-r\le 0$, returns `None` for a ball already touching at $\tau=0$ (recovering an existing overlap is the discrete solver's job) and after 32 iterations conservatively reports contact. How the solver uses it is in Section 2.2 (step 2.5).

### 2.2 Grinding media: rigid discs, XPBD contact

Source: `crates/mill-core/src/dem.rs`. Entry point: `DemState::step_recorded` (`dem.rs:530`).

**Contrast with spring-dashpot DEM.** Classical DEM [3], and the mill DEM built on it [4, 5], computes contact forces from overlap-dependent springs and dashpots; stability then requires a time step smaller than the contact duration, microseconds for stiff media at wall speeds of order 1 m/s. Here contact is a rigid geometric constraint solved by projection, in the lineage of position based dynamics [1] and its compliance-aware extension XPBD [2] with compliance $\alpha=0$. It is stable at any sub-step; the price is that contact forces are reconstructed from position corrections, not integrated, and that contact mechanics are not those of a calibrated spring-dashpot law. Position-based dynamics originates in computer graphics [1, 2, 9, 11] and has precedents in engineering use (outside the scope of this document).

**Body properties.** Balls are uniform discs of unit depth (`ball_mass`, `dem.rs:230`; `ball_inertia`, `dem.rs:237`):

$$
m=\rho\,\pi r^2,\qquad I=\tfrac12 m r^2,\qquad w=m^{-1},\quad w_\text{rot}=I^{-1}.
$$

All balls of a population share one radius, mass and inertia. Mechanical energy per metre is $E=\sum_i\bigl(\tfrac12 m|v_i|^2+\tfrac12 I\omega_i^2+m g\,y_i\bigr)$ (`mechanical_energy_j`, `dem.rs:253`). Balls are seeded on a hexagonal lattice of spacing $2r\cdot1.02$ with random jitter of $\pm0.05r$ from the seeded PRNG (`Balls::seed_lattice`, `dem.rs:355`).

**Per-sub-step algorithm** (given $\Delta t$, $n_\text{it}$ = `dem_iterations`, default 2, `params.rs:490`):

0. **Sanitise.** Non-finite $x$ or $v$ are reset to zero (`dem.rs:552`).

1. **Predict** (`dem.rs:565`). For each ball, with fluid impulses $J_i$ (linear) and $\Gamma_i$ (angular) from the coupling:

$$
v_i \leftarrow v_i + (0,-g)\Delta t + w\,J_i,\qquad \omega_i\leftarrow\omega_i+w_\text{rot}\Gamma_i,\qquad x_i\leftarrow x_i+v_i\Delta t,\quad \theta_i\leftarrow\theta_i+\omega_i\Delta t .
$$

   The impulses enter as a direct velocity change with no extra factor $\Delta t$. The pre-step state $(x^0,\theta^0)$, the pre-predict velocity $v^\text{pre}$ and $E_\text{pred}=E$ after this step are saved. The tunnelling diagnostic $\max_i |v_i|\Delta t/(2r)$ is recorded here (`dem.rs:591`).

2. **Broad phase** (`dem.rs:597`). A counting-sort `UniformGrid` (Section 2.9) is built on the predicted positions with cell size

$$
s=\max\Bigl(2r\cdot1.05,\ 2r+v_{\max}\Delta t+n_\text{it}\cdot0.2\cdot 2r\Bigr),
$$

   (`dem.rs:626`-`dem.rs:630`), i.e. the largest predicted displacement of the sub-step plus the maximum motion the solve itself can add: each of the $n_\text{it}$ passes can move a ball by up to `MAX_RECOVERY_FRACTION` $\times 2r$. The candidate pair list is built once and reused for all passes and for friction, restitution and rolling resistance. Removing the solve-margin term caused a measured full-system energy-invariant violation (comment at `dem.rs:618`-`dem.rs:625`). Ball-wall contact needs no broad phase: one SDF query per ball.

3. **Swept wall clamp** (step 2.5, `dem.rs:635`-`dem.rs:706`). For each ball, the end-of-step constraint value $C_\text{end}=d(x_i;\theta_D+\omega\Delta t)-r$ is evaluated. Only if the predicted overlap is deeper than one bounded-recovery pass can absorb, $C_\text{end}<-0.2\cdot 2r$, `Drum::toi_swept` (Section 2.1) is called; if it reports contact at fraction $\tau$, the predicted advance is truncated to

$$
x_i \leftarrow x^0_i + (x_i-x^0_i)\min\Bigl(\tau+\frac{0.01\,r}{|x_i-x^0_i|},1\Bigr),
$$

   a small physical skin (`CCD_CONTACT_SKIN_FRACTION = 0.01`, `dem.rs:72`) past the exact crossing so the contact is a real overlap and is seen by friction and restitution. Swept collision detection is deliberately applied to wall and lifters only: a ball-ball version froze motion in dense beds and made packed charges worse (comment at `dem.rs:645`-`dem.rs:657`), and an unconditional wall version repeatedly truncated ordinary tangential sliding along the curved wall.

4. **Non-penetration projection** (`dem.rs:710`-`dem.rs:750`), $n_\text{it}$ Gauss-Seidel passes, each over all candidate pairs and then all ball-wall contacts. Ball-ball constraint and bounded correction:

$$
C_{ij}=|x_i-x_j|-(r_i+r_j),\qquad \hat n=\frac{x_i-x_j}{|x_i-x_j|},
$$

$$
\text{if } C_{ij}<0:\quad C_\text{eff}=\max\bigl(C_{ij},\,-f_\text{rec}\cdot 2r\bigr),\quad f_\text{rec}=\texttt{MAX\_RECOVERY\_FRACTION}=0.2,
$$

$$
\Delta\lambda=\frac{-C_\text{eff}}{w_i+w_j},\qquad x_i\mathrel{+}=w_i\,\Delta\lambda\,\hat n,\quad x_j\mathrel{-}=w_j\,\Delta\lambda\,\hat n,\quad \lambda^n_{ij}\mathrel{+}=\Delta\lambda .
$$

   This is the XPBD multiplier update with zero compliance [2]. Ball-wall: $C_i=d(x_i;\theta_D+\omega\Delta t)-r$, $\hat n$ the SDF normal, wall inverse mass 0, so $\Delta\lambda=-C_\text{eff}/w_i$ and only the ball moves. The cap `MAX_RECOVERY_FRACTION` (`dem.rs:58`) prevents a deep overlap from being recovered in one pass, which step 5 would turn into a large separation velocity $\Delta x/\Delta t$; the accumulated multiplier $\lambda^n$ is computed from the capped correction, so a partially recovered contact carries proportionally less normal load.

5. **Velocity reconstruction** (`dem.rs:752`): $v_i=(x_i-x^0_i)/\Delta t$, $\omega_i=\operatorname{angdiff}(\theta_i,\theta^0_i)/\Delta t$.

6. **Friction** (`dem.rs:765`-`dem.rs:890`). Contacts are processed in a deterministic $(i,j)$-sorted order (`ContactBook::sorted_ball_ball`, `dem.rs:459`). With tangent $\hat t=(-\hat n_y,\hat n_x)$, the relative tangential contact-point speed is

$$
v_t^{bb}=(v_i-v_j)\cdot\hat t-r\omega_i-r\omega_j,\qquad v_t^{bw}=(v_i-v_\text{wall}(x_c))\cdot\hat t-r\omega_i,
$$

   with the wall velocity sampled at the contact point $x_c=x_i-r\hat n$ (not the ball centre). The correction that would cancel the slip in this sub-step is clamped to the Coulomb cone and applied as a position correction:

$$
w_t=w_i+w_j+r^2 w_\text{rot,i}+r^2 w_\text{rot,j}\ \ (\text{ball-wall: } w_i+r^2 w_\text{rot,i}),\qquad
\Delta\lambda_t=\operatorname{clamp}\!\Bigl(\frac{-v_t\,\Delta t}{w_t},\ \pm\mu_\text{eff}\,\lambda^n\Bigr),
$$

$$
x_i\mathrel{+}=w_i\Delta\lambda_t\hat t,\quad x_j\mathrel{-}=w_j\Delta\lambda_t\hat t,\quad \theta_{i,j}\mathrel{-}=r\,w_\text{rot}\Delta\lambda_t .
$$

   The velocity and spin are updated in the same instant, $\Delta v=w\,\Delta\lambda_t\hat t/\Delta t$, $\Delta\omega=-r\,w_\text{rot}\Delta\lambda_t/\Delta t$, so a later contact of the same ball sees its true current state; each single-contact correction then only reduces relative slip and cannot add energy. The effective coefficient blends static and kinetic friction smoothly in the slip speed (`effective_friction`, `dem.rs:128`):

$$
\mu_\text{eff}(v_t)=\mu_k+(\mu_s-\mu_k)\,e^{-(v_t/v_0)^2},\qquad v_0=\texttt{friction\_velocity\_scale\_m\_s}=0.02\ \mathrm{m/s}.
$$

   At the defaults $\mu_s=\mu_k$ (ball-ball 0.25, ball-wall 0.35, `params.rs:216`-`params.rs:219`) and $\mu_\text{eff}$ is the constant $\mu_k$. The wall's work this sub-step, which becomes the motor power, is the tangential and the normal wall impulse dotted with the wall velocity at the contact point:

$$
W_\text{wall}=\sum_{\text{ball-wall}}\Bigl[\frac{\Delta\lambda_t}{\Delta t}\,\hat t\cdot v_\text{wall}+\frac{\lambda^n}{\Delta t}\,\hat n\cdot v_\text{wall}\Bigr]
$$

   (`dem.rs:870`, `dem.rs:881`). On a smooth wall $\hat n\cdot v_\text{wall}=0$, so only the lifter faces contribute the second term. Velocities are re-reconstructed from positions after the pass (`dem.rs:887`) as a self-healing resync.

7. **Restitution** (`dem.rs:892`-`dem.rs:1019`), bidirectional. Using the approach speed before the solve, $v_n^\text{pre}=(v^\text{pre}_i-v^\text{pre}_j)\cdot\hat n$ (ball-wall: relative to $v_\text{wall}$), a contact is a fresh impact if $v_n^\text{pre}<-v_\text{thr}$, $v_\text{thr}=0.02$ m/s (`RESTITUTION_VELOCITY_THRESHOLD`, `dem.rs:32`). The coefficient is $e$ for a fresh impact and 0 otherwise, and the actual relative normal speed is driven to a target in both directions:

$$
v_n^\text{target}=\max\bigl(0,\ -e\,v_n^\text{pre},\ v_n^\text{pre}\bigr),\qquad
\Delta v_n=v_n^\text{target}-v_n^\text{now},
$$

   applied as an impulse $\Delta v_n/(w_i+w_j)$ split by inverse mass (ball-wall: $\Delta v_n\hat n$ on the ball) (`dem.rs:961`, `dem.rs:1012`). Driving it in both directions removes separation velocity manufactured by step 4's depenetration, which is what bounds the energy; the third term floors the target at the pair's pre-solve separation speed so a pair that was already separating is never pulled back together (the one-sided Signorini condition). Impact kinetic energy for statistics is $E=\tfrac12 (m/2)\,{v_n^\text{pre}}^2$ for ball-ball and $\tfrac12 m\,{v_n^\text{pre}}^2$ for ball-wall (`dem.rs:924`, `dem.rs:986`). For the recorded impact log used by the stress statistics, the threshold is size-relative, $v_\text{thr}(r)=0.1428\sqrt{2gr}$ (`dem.rs:44`; equal to 0.02 m/s at $d=2$ mm).

8. **Rolling resistance** (`dem.rs:1021`-`dem.rs:1055`). With $\Lambda_i$ the sum of the ball's normal multipliers over its contacts,

$$
\Delta\omega_\text{max}=\mu_r\,\frac{\Lambda_i}{\Delta t}\,\frac{r}{I}\,\Delta t,\qquad
\omega_i\leftarrow\omega_i-\operatorname{sign}(\omega_i)\min(\Delta\omega_\text{max},|\omega_i|),
$$

   $\mu_r=0.01$ by default. This damps the world-frame spin toward zero, not the relative spin at a contact (known simplification: it exerts no reaction torque on a partner and, at a wall contact, opposes rather than assists rolling without slip on a rotating drum).

**Energy accounting and the no-energy-creation property.** Net dissipation over the contact solve is

$$
E_\text{diss}=E_\text{pred}+W_\text{wall}-E_\text{final}
$$

(`dem.rs:1062`). Every stage after the predict step is non-increasing in mechanical energy except the work done by the moving wall: bounded depenetration, per-contact friction clamped to the cone, bidirectional restitution with $e\le1$, and rolling resistance. The property is tested every sub-step, not only at the end: $E_\text{final}\le E_\text{pred}+W_\text{wall}$ for the dry charge (`crates/mill-core/src/lib.rs:759`, `dem.rs:1784`), and for the whole ball plus fluid system, $\Delta E_\text{total}\le W^\text{DEM}_\text{wall}+W^\text{fluid}_\text{wall}+\text{tolerance}$ (`crates/mill-core/src/coupling.rs:1149`). $E_\text{diss}$ excludes fluid-side viscous dissipation, so with slurry on it is expected to be below the power draw. The step-3 cap and the step-7 bidirectional restitution together removed an earlier energy-injection failure in which a fluidised charge inflated like a gas.

### 2.3 Coarse-graining and automatic accuracy scaling

Source: `crates/mill-core/src/params.rs`.

**True population** (`Params::true_ball_count`, `params.rs:556`). With fill fraction $J$ (charge footprint including voids, default 0.30), 2D areal packing fraction $\phi_2$ (default 0.82; random close packing of equal discs $\approx0.82$, hexagonal limit $\pi/(2\sqrt3)\approx0.907$), drum area $A=\pi R^2$ and true ball radius $r_t$:

$$
N_\text{true}=\frac{J\,\phi_2\,\pi R^2}{\pi r_t^2}.
$$

At the shipped 2 mm balls in a 63 mm drum this is $\approx244$.

**Coarse-graining** (`Params::effective_media`, `params.rs:582`). The mode is explicit, `CoarseGrainingMode` (`params.rs:35`): `Off` (default, $k=1$), `Auto` ($k=\sqrt{N_\text{true}/N_\max}$ when $N_\text{true}>N_\max$, `simulation.max_balls`, else 1) and `Manual` ($k=$ `coarse_graining_k`). For $k>1$ a smaller number of larger discs replaces the true population:

$$
d_\text{eff}=k\,d_\text{true},\qquad \rho_\text{eff}=\rho_\text{true},\qquad N_\text{sim}=\operatorname{round}(N_\text{true}/k^2).
$$

Because the unit-depth disc mass scales as $r^2$, $N_\text{sim}\,\pi r_\text{eff}^2=N_\text{true}\,\pi r_t^2$ preserves both the footprint area and the total charge mass at unchanged density. It does not preserve the interstitial void structure or per-impact statistics at the true particle size; `EffectiveMedia` (`params.rs:709`) carries both sizes so the approximation is never silent. The idea of replacing many small particles by fewer larger ones is established in granular simulation [6, 7] and in wet bead milling [8] (context: the first two are not mill studies; the third uses an enlarged-particle model for wet bead milling). The true radius is kept in `Balls::true_radius` (`dem.rs:308`) and used only by the form-drag rate (Section 2.5).

**Automatic fluid resolution** (`Params::effective_fluid_resolution`, `params.rs:648`). The fluid lattice spacing is $\Delta x = R/\text{res}$. To prevent one fluid particle from spanning several balls (which saturated the coupling clamps and injected energy), the resolution is raised so that $\Delta x\le d_\text{eff}$:

$$
\text{res}_\text{eff}=\min\bigl(200,\ \max(\text{res},\ \lceil R/d_\text{eff}\rceil)\bigr).
$$

**Automatic sub-step count** (`Params::effective_substeps`, `params.rs:690`). The XPBD tunnelling criterion is $\max|v|\Delta t/(2r) < 1$, which worsens as $d_\text{eff}$ shrinks. Using a generous free-fall reference speed across the drum diameter $D$, $v_\text{ref}=\sqrt{2gD}$, and requiring $v_\text{ref}\Delta t\le \texttt{TARGET\_RATIO}\cdot d_\text{eff}$ with $\Delta t=1/(60 n_s)$ and `TARGET_RATIO = 0.5` (`params.rs:692`):

$$
n_{s,\min}=\Bigl\lceil \frac{v_\text{ref}}{60\cdot0.5\cdot d_\text{eff}}\Bigr\rceil,\qquad n_s=\max\bigl(\texttt{substeps},\ \min(n_{s,\min},16)\bigr).
$$

An explicit request above 16 (up to 64) is honoured as is (`params.rs:702`). At the defaults, $v_\text{ref}=1.11$ m/s, $n_{s,\min}=19$, so $n_s=16$.

**Tunnelling diagnostic.** `max_substep_displacement_over_diameter` $=\max_i|v_i|\Delta t/(2r)$ per sub-step (`dem.rs:591`) is reported as a metric; the maximum over a frame is kept by `Simulation` (`lib.rs:562`). Values approaching 1 indicate that a ball moves about its own diameter in one sub-step.

### 2.4 Slurry: Position Based Fluids

Source: `crates/mill-core/src/pbf.rs`. Entry point: `FluidParticles::step_coupled` (`pbf.rs:1173`).

The slurry is a set of particles of equal areal mass $m_f$, kernel radius $h=2\Delta x$ and rest density $\rho_0$ (`slurry.density_kg_m3`, default 1000 kg/m³). Position Based Fluids [9] enforce incompressibility as a density constraint solved by position projection, in the position-based-dynamics framework [1, 11].

**Seeding** (`FluidParticles::seed_lattice`, `pbf.rs:1043`). Spacing $\Delta x=R/\text{res}_\text{eff}$, $h=2\Delta x$. Particles are placed on a hexagonal lattice, columns $\Delta x$ apart and rows $\Delta x\sqrt3/2$ apart, filling $\texttt{slurry.fill\_fraction}$ (default 0.35) of the drum area from the bottom upwards, skipping sites inside a ball. The areal mass is the hexagonal cell area times density,

$$
m_f=\rho_0\,\Delta x\cdot\frac{\sqrt3}{2}\Delta x,
$$

so the seeded Poly6 density sum already equals $\rho_0$ to about 0.1% and the run starts without a compression transient. Initial velocities are zero; a passive dye tracer in $\{0,1\}$ (left/right, top/bottom or none) is carried by each particle for mixing statistics.

**Kernels** (2D, each normalised to integrate to 1 over the disc of radius $h$). Poly6 for density and Spiky for the pressure gradient are the kernels of [10]; the 2D constants are re-derived here (`pbf.rs:118`, `pbf.rs:130`):

$$
W(r,h)=\frac{4}{\pi h^{8}}\,(h^2-r^2)^3\ \ (r<h),\qquad
\nabla W(\boldsymbol\delta,r,h)=-\frac{30}{\pi h^{5}}\,(h-r)^2\,\frac{\boldsymbol\delta}{r}\ \ (0<r<h),
$$

with $\boldsymbol\delta=x_i-x_j$. The gradient points from $x_i$ toward $x_j$ because the kernel decreases with $r$.

**Density.** $\rho_i=m_f\sum_j W(|x_i-x_j|,h)$ including the self term $m_f W(0,h)$.

**One-sided density constraint** (`pbf.rs:1241`-`pbf.rs:1333`). Per particle,

$$
C_i=\max\!\Bigl(0,\ \frac{\rho_i}{\rho_0}-1\Bigr),\qquad
\mathbf g_{ij}=\frac{m_f}{\rho_0}\nabla W_{ij},\qquad \mathbf g_{ii}=\sum_j \mathbf g_{ij},
$$

$$
\lambda_i=\frac{-C_i}{\sum_j|\mathbf g_{ij}|^2+|\mathbf g_{ii}|^2+\varepsilon},\qquad
\Delta x_i=\frac{m_f}{\rho_0}\sum_j(\lambda_i+\lambda_j)\nabla W_{ij},
$$

where $\varepsilon=$ `EPSILON_RELAX` $=200$ (`pbf.rs:48`) is the constraint-force-mixing relaxation of [9]. Two differences from the published formulation: the particle mass $m_f/\rho_0$ appears in the gradient because the constraint is written with $\rho_i=m_f\sum W$ (omitting it made the correction too weak by a factor of order $1/m_f$ at small drum sizes), and the constraint is clamped to $C_i\ge0$, so only over-dense particles are projected apart. The published method leaves $C_i$ unclamped and suppresses the resulting tensile instability with an artificial-pressure term; here the clamp removes the instability at its source, and fluid-fluid attraction is supplied only by the explicit cohesion term below. All $\Delta x_i$ are computed from the current positions and then applied together (Jacobi style).

**Adaptive iteration count.** `simulation.pbf_iterations` (default 3) is a floor. After each pass the loop stops if at least that many passes have run and the largest compression error measured at the start of the pass, $\max_i C_i$, is $\le$ `COMPRESSION_CONVERGENCE_TOLERANCE` $=0.02$ (`pbf.rs:61`); otherwise it continues up to `ADAPTIVE_MAX_EXTRA_ITERATIONS` $=6$ further passes (`pbf.rs:57`, `pbf.rs:1326`-`pbf.rs:1332`). A violent impact that leaves a large residual therefore costs extra passes only when needed.

**Ball overlap projection (step 3.5).** Described in Section 2.5 because it exchanges momentum with the balls.

**Cohesion (step 3.6a)** [12] (`pbf.rs:1455`-`pbf.rs:1490`). For each particle, over the same neighbour list,

$$
\mathbf p_i=\sum_j \frac{m_f\,K(r_{ij},h)}{\rho_0}\,\hat{\mathbf n}_{ij},\qquad
a_c=\frac{\sigma}{\sigma_\text{ref}}\,A_c\,g,\qquad
\Delta v_i=-a_c\,\mathbf p_i\,\Delta t ,
$$

with $\hat{\mathbf n}_{ij}=(x_i-x_j)/r_{ij}$ so that the net kick points toward the neighbours; $\sigma$ = `slurry.surface_tension_n_m` (default 0.072 N/m), $\sigma_\text{ref}=0.072$ N/m (`REFERENCE_SURFACE_TENSION_N_M`, `pbf.rs:80`), $A_c=2$ (`COHESION_ACCEL_FACTOR`, `pbf.rs:85`). The kernel is the cubic shape of Akinci et al. (their eq. 2), zero at $r=0$ and $r=h$, peaking near $0.6h$ and slightly negative for $r<h/2$, re-normalised for 2D (`cohesion_kernel`, `pbf.rs:152`):

$$
K(r,h)=\frac{35840}{209\pi h^{8}}\times\begin{cases}2\,(h-r)^3 r^3-\dfrac{h^6}{64}, & r\le h/2\\[4pt] (h-r)^3 r^3, & h/2<r<h .\end{cases}
$$

The term $m_f K/\rho_0$ is dimensionless, about 1 deep in the bulk and tapering at a free surface, so $a_c$ is a genuine peak-acceleration scale. A particle with a symmetric neighbourhood feels no net force; only free-surface particles are pulled inward. The mapping from $\sigma$ to the SPH coefficient is a calibrated proportionality, not a first-principles unit conversion (as in the original method, where it is resolution dependent). The kick is applied to the velocity after the position-based velocity reconstruction (step 5), which would otherwise discard it.

**Adhesion (step 3.6b)** is a ball-fluid force and is described in Section 2.5.

**Boundary projection** (step 4, `pbf.rs:1633`). A particle with $d(x;\theta_D+\omega\Delta t)<0$ is moved to the surface, $x\mathrel{+}=-d\,\hat n$, and flagged. Unlike the ball solver, the wall here is a position projection only; no-slip is applied afterwards at velocity level.

**Velocity reconstruction** (step 5, `pbf.rs:1654`). $v_i=(x_i-x_i^0)/\Delta t$ where $x^0$ is the position before the predict step. The geometric excess of the overlap push is subtracted and the deferred cohesion and adhesion kicks are added.

**Wall no-slip blending** (step 6, `pbf.rs:1675`). For flagged particles,

$$
v_i\leftarrow(1-\beta_w)\,v_i+\beta_w\,v_\text{wall}(x_i),\qquad\beta_w=\texttt{slurry.wall\_no\_slip}\in[0,1]\ (\text{default }1),
$$

and the work $m_f\,\Delta v\cdot v_\text{wall}$ is summed into the fluid wall work, which enters power and torque (Section 2.7).

**Implicit viscosity** (step 7, `pbf.rs:2122`-`pbf.rs:2293`). Density is recomputed at the final positions. The viscous term is discretised with the Morris weights [13] (`morris_weights`, `pbf.rs:311`),

$$
c_{ij}=\frac{m_f\,(\mu_i+\mu_j)}{\rho_i\rho_j}\,\frac{|\boldsymbol\delta_{ij}\cdot\nabla W_{ij}|}{|\boldsymbol\delta_{ij}|^2+\eta^2},\qquad \eta^2=0.01\,h^2
$$

(`VISCOSITY_ETA_FACTOR`, `pbf.rs:106`; $\rho_i$ is each particle's own measured SPH density), but applied in a central-force form that acts along the line of centres only:

$$
(L\mathbf v)_i=\sum_j 4\,c_{ij}\ \hat{\mathbf e}_{ij}\,\bigl[\hat{\mathbf e}_{ij}\cdot(\mathbf v_i-\mathbf v_j)\bigr],\qquad \hat{\mathbf e}_{ij}=\frac{x_i-x_j}{|x_i-x_j|}
$$

(`laplacian_apply_into`, `pbf.rs:358`; factor `CENTRAL_VISCOSITY_FACTOR = 4`, `pbf.rs:383`, which is the 2D angular-average factor that recovers the same Laplacian for incompressible flow). The plain Morris form $c_{ij}(\mathbf v_i-\mathbf v_j)$ exerts a spurious force on a rigidly rotating fluid wherever the kernel sum is one-sided (free surface, wall); the central form vanishes exactly for rigid rotation and conserves angular momentum (unit test `pbf.rs:2421`). Implicit (backward Euler) time stepping makes it unconditionally stable for any $\mu\ge0$:

$$
(I+\Delta t\,L)\,\mathbf v^{\text{new}}=\mathbf v^{*}.
$$

$A=I+\Delta t L$ is symmetric positive definite for $\Delta t>0$ since $L$ is symmetric positive semi-definite ($c_{ij}\ge0$, symmetric in $i,j$). It is solved matrix-free by conjugate gradient on the $2n$ real unknowns (`solve_implicit_viscosity`, `pbf.rs:422`), warm-started from the incoming velocity, with relative residual tolerance $|r|\le10^{-3}|b|$ (`VISCOSITY_CG_TOLERANCE`, `pbf.rs:96`) and at most 50 iterations (`VISCOSITY_CG_MAX_ITERS`, `pbf.rs:100`). It is skipped when $\mu=0$. When balls are near the fluid, the system is instead solved jointly with the ball velocities (Section 2.5).

**Bingham rheology** (`bingham_viscosity`, `pbf.rs:853`; used if `slurry.rheology == Bingham` and the yield stress $\tau_y>0$, `pbf.rs:2141`). The per-particle shear rate is the standard SPH velocity-gradient estimate (`shear_rates`, `pbf.rs:820`),

$$
\nabla\mathbf v_i=\sum_j\frac{m_f}{\rho_j}(\mathbf v_j-\mathbf v_i)\otimes\nabla W_{ij},\qquad D=\operatorname{sym}(\nabla\mathbf v),\qquad \dot\gamma=\sqrt{2\,D\!:\!D},
$$

evaluated from the velocities entering the viscosity step. The effective viscosity is the Papanastasiou-regularised Bingham law (exponential regularisation of the yield surface; no key for the original paper is held in the reference list):

$$
\mu_\text{eff}(\dot\gamma)=\mu+\tau_y\,\frac{1-e^{-\dot\gamma/\dot\gamma_0}}{\dot\gamma},\qquad \dot\gamma_0=\frac{\tau_y}{R_p\,\mu},\quad R_p=\texttt{BINGHAM\_PLASTIC\_RATIO}=100\ (\texttt{pbf.rs:489}),
$$

with the limit $\dot\gamma\to0$ giving $\mu(1+R_p)$, i.e. the unyielded viscosity is capped at 101 $\mu$. $\mu_\text{eff}$ replaces $\mu$ in $c_{ij}$ and in the ball-surface boundary weights. The mean of $\dot\gamma$ over all particles is reported as `mean_shear_rate_per_s` (`pbf.rs:799`). Against the grid reference the Bingham response of the PBF solver is documented as only qualitatively right (`docs/VERIFICATION.md`, "PBF: Bingham and coarse-graining k").

**Fluid speed clamp** (step 7.5, `pbf.rs:2303`; `clamp_fluid_speed`, `pbf.rs:986`). Any particle speed above

$$
v_{\max}=5\Bigl(|\omega|R+\sqrt{4gR}\Bigr)
$$

(`FLUID_SPEED_SAFETY_FACTOR = 5`, `pbf.rs:973`; `fluid_speed_clamp_v_max`, `pbf.rs:978`) is rescaled to $v_{\max}$, and non-finite velocities are zeroed. It is the wall speed plus the free-fall speed across the drum diameter, times 5: a backstop against position-corrections accumulating into unphysical speeds, not a model parameter. Momentum removed here is not returned to any ball.

### 2.5 Ball-fluid coupling

Sources: `crates/mill-core/src/coupling.rs`, `crates/mill-core/src/pbf.rs`.

**Impulse, not force.** Each fluid mechanism accumulates a linear impulse $J_b$ and an angular impulse $\Gamma_b$ on ball $b$ (`CouplingImpulses`, `coupling.rs:37`), applied in the ball predict step as $\Delta v=J_b/m$, $\Delta\omega=\Gamma_b/I$ with no extra factor $\Delta t$. Expressing the exchange as force and re-integrating with another $\Delta t$ made the system sensitive to the sub-step and unstable. The matching reaction is applied to the fluid in the same call, so a Newton's-third-law cross-check against the summed impulses is possible.

**Overlap projection with bounded momentum (step 3.5, `pbf.rs:1337`-`pbf.rs:1443`).** After the density solve, a fluid particle within

$$
r_c=r_b+0.05\,h\qquad (\texttt{OVERLAP\_PUSH\_MARGIN}=0.05,\ \texttt{pbf.rs:70})
$$

of a ball centre is pushed out along the normal $\hat n=(x_i-x_b)/|x_i-x_b|$ by $\Delta=\min(r_c-|x_i-x_b|,\ r_b/2)$. The position push is a geometric constraint; the momentum it may exchange is bounded separately by what the contact actually needs to arrest. With the ball's surface velocity $v_s=v_b+\omega_b\,\mathrm{perp}(\hat n\,r_b)$,

$$
u_\text{arrest}=\max\bigl(0,\ -(v_i-v_s)\cdot\hat n\bigr),\qquad
\mu_r=\frac{m_f m_b}{m_f+m_b},\qquad
J=\min\Bigl(\frac{m_f\Delta}{\Delta t},\ \mu_r\,u_\text{arrest}\Bigr),
$$

the ball receives $-J\hat n$ (and angular impulse $(r_b\hat n)\times(-J\hat n)$), and the part of the push beyond $J$ is subtracted from the particle's reconstructed velocity so it never appears as implicit fluid speed. A resting overlap therefore exchanges no momentum, only a position correction; the reduced mass keeps a heavy fluid particle from handing an oversized impulse to a light ball.

**Adhesion (step 3.6b, `pbf.rs:1493`-`pbf.rs:1630`)** [12] (shape only). Each ball's circumference is sampled into boundary points at about $\Delta x$ arc spacing (at least 6, `ball_boundary_particles`, `pbf.rs:189`). With $h_a=\min(h,2r_b)$ and the adhesion kernel, the outer half of the cohesion shape re-normalised on $(h/2,h)$ (`adhesion_kernel`, `pbf.rs:174`; Akinci's own fractional-power adhesion formula is not used):

$$
K_a(r,h)=\frac{35840}{163\pi h^{8}}(h-r)^3r^3\ \ (h/2<r<h),
$$

contributions of all boundary points of ball $b$ to each fluid particle $j$ are merged, giving weights $w_j=m_f K_a/\rho_0$ and directions $\hat{\mathbf n}_j$ toward the ball. The per-ball response is bounded, not per particle:

$$
W=\sum_j w_j,\quad \text{cov}=\min(1,W),\quad m_c=m_f W,\quad \text{att}=\frac{m_c}{m_b+m_c},\quad
a_a=\text{wet}\cdot A_a\,g\ (A_a=2,\ \texttt{pbf.rs:91}),
$$

$$
\Delta v_b=\hat{\mathbf d}\ a_a\ \text{cov}\ \Delta t\ \text{att},\qquad \hat{\mathbf d}=-\frac{\sum_j w_j\hat{\mathbf n}_j}{|\sum_j w_j\hat{\mathbf n}_j|}\ ,
$$

the ball receives $m_b\Delta v_b$ and each contributing particle $j$ receives $-m_b (w_j/W)\Delta v_b/m_f$, so the fluid's total momentum change is exactly $-m_b\Delta v_b$. With cohesion this gives the slurry an emergent, empirically calibrated contact angle; wet = `slurry.wettability`, default 0.6.

**Form drag (step 6.5, `pbf.rs:1702`-`pbf.rs:2006`).** Each ball relaxes toward the mean velocity of its entrained fluid, a centre-of-mass closure bounded by construction. Fluid particles within $h_c=r_b+h$ of the ball centre contribute with the smooth Poly6 taper $\phi_i=W(|r_i|^2,h_c)/W(0,h_c)$, $\mathbf r_i=x_i-x_b$:

$$
m_\text{ent}=\sum_i m_f\phi_i,\qquad \bar{\mathbf v}=\frac{\sum_i m_f\phi_i\mathbf v_i}{m_\text{ent}},\qquad \mathbf v_\text{rel}=\bar{\mathbf v}-\mathbf v_b,
$$

$$
\frac1{\tau_\text{form}}=\frac{\rho_0\,C_D\,(2r_t)\,|\mathbf v_\text{rel}|}{m_t},\qquad
\beta=\beta_b\bigl(1-e^{-\Delta t/\tau_\text{form}}\bigr),\qquad
a_\text{lin}=\beta\,\frac{m_\text{ent}}{m_b+m_\text{ent}}\le1,\qquad
\Delta\mathbf v_b=a_\text{lin}\,\mathbf v_\text{rel},
$$

with $C_D=1$ (`BALL_FORM_DRAG_COEFFICIENT`, `pbf.rs:115`, the high-Reynolds-number circular-cylinder value), $\beta_b=$ `slurry.ball_no_slip`, and $r_t$, $m_t=\rho_b\pi r_t^2$ the true (un-coarse-grained) ball radius and mass, so the drag rate does not depend on the coarse-graining knob; the entrained-mass ratio and the kernel footprint use the coarse-grained $m_b$, $r_b$ because both scale together. The ball receives $m_b\Delta\mathbf v_b$ and each entrained particle $i$ the reaction $-(m_b/m_\text{ent})\phi_i\Delta\mathbf v_b$, so linear momentum is conserved exactly. The angular part is the mirror image with the inertia $I_b$, a moment-weighted mean angular velocity $\bar\omega=\sum_i m_f\phi_i\,\mathbf r_i\times\mathbf v_i/\sum_i m_f\phi_i|\mathbf r_i|^2$ and a tangential reaction field with zero net linear momentum, scaled by a coefficient $c_\text{rot}=-I_b\,\Delta\omega_b/s$ with $s=\sum_i m_f\phi_i^2\,\mathbf r_i\times(\mathrm{perp}(\mathbf r_i)-\bar{\mathbf p})$; the rotational exchange is skipped when $|s|\le0.1\sum_i m_f\phi_i^2|\mathbf r_i|^2$ (near-degenerate geometry) and $c_\text{rot}$ is capped so that no particle's velocity kick exceeds $0.2\,v_{\max}$ (with $v_{\max}$ the fluid speed clamp above), re-deriving $\Delta\omega_b$ from the capped $c_\text{rot}$ so the pair stays consistent. Balls are processed sequentially (Gauss-Seidel) because a fluid particle can lie within the entrainment radius of several balls. The step is gated on $\mu>10^{-6}$ (`pbf.rs:1774`). The Stokes term of the older closure is disabled (`inv_tau_stokes = 0.0`, `pbf.rs:1800`): viscous drag now comes from the no-slip boundary of the implicit solve below, and only the inertial form-drag relaxation remains here.

**Buoyancy in the rotating frame (step 6.6, `pbf.rs:2012`-`pbf.rs:2115`).** A ball sampling fluid at its centre with the Poly6 kernel gets $\rho_\text{loc}=\sum_j m_f W(|x_b-x_j|,h)$ and $\rho_\text{eff}=\min(\rho_\text{loc},\rho_0)$. The apparent gravity in the fluid's own rest frame, for a neighbourhood rotating at the locally measured rate $\omega_\ell$ about the drum centre, is

$$
\mathbf g_\text{eff}=(0,-g)+\omega_\ell^2\,\mathbf x_b,\qquad \omega_\ell=\operatorname{clamp}\Bigl(\frac{\sum_j w_j\,(\mathbf x_j\times\mathbf v_j)/|\mathbf x_j|^2}{\sum_j w_j},\ \pm|\omega|\Bigr),\quad w_j=m_fW_j,
$$

and the impulse on the ball is the Archimedes expression $-\rho_\text{eff}\,(\pi r_b^2)\,\mathbf g_\text{eff}\,\Delta t$, acting through the centroid (no torque), with the reaction distributed over the sampled particles in proportion to $w_j/\rho_\text{loc}$. Using the measured $\omega_\ell$ rather than $\omega$ avoids a fictitious inward force on a still bottom pool; for a freely rotating, centrifuged pool it recovers the full centripetal term. Because the fluid density is sampled at the ball centre, this term is inactive for balls larger than $h$ (the density there is zero); buoyancy then arises from the overlap-projection pressure.

**Joint implicit no-slip solve (step 7, `solve_viscosity_with_balls`, `pbf.rs:640`).** Viscous drag on a ball is obtained by treating the ball surface as a no-slip (Dirichlet) boundary of the viscosity operator and solving for the fluid and ball velocities together. The ghost points $k$ on the circumference of ball $b$ (spacing $\Delta x=h/2$) carry the rigid-body surface velocity $\mathbf u_k=\mathbf V_b+W_b\,\mathrm{perp}(\mathbf l_k)$, with $\mathbf l_k$ the offset from the ball centre. Each fluid particle $i$ within $h$ of a ghost point is coupled with the Morris-type weight (isotropic, not central, `pbf.rs:2166`-`pbf.rs:2202`)

$$
c_{ik}=G\,\frac{m_f\cdot2\mu_i}{\rho_i\rho_0}\ \frac{\bigl|\boldsymbol\delta_{ik}\cdot\nabla W_{ik}\bigr|}{|\boldsymbol\delta_{ik}|^{2}+\eta^2},\qquad G=\texttt{BALL\_BC\_GAIN}=0.5\ (\texttt{pbf.rs:486}),
$$

(ghost points have mass $m_f$ and density $\rho_0$; the gain $G$ compensates for the single ghost layer being thinner than a true wall half-space and was calibrated on the grid reference, `docs/VERIFICATION.md`, "PBF ball no-slip boundary"). All pairs of one (fluid, ball) couple are merged into three sums $c=\sum c_{ik}$, $\sum c_{ik}\mathbf l_k$, $\sum c_{ik}|\mathbf l_k|^2$ (`BoundaryGroup`, `pbf.rs:497`). The unknowns are $\mathbf v_i$, $\mathbf V_b$, $W_b$ and the system is

$$
\begin{aligned}
m_f\mathbf v_i+\Delta t\,m_f(L\mathbf v)_i+\Delta t\,m_f\sum_k c_{ik}(\mathbf v_i-\mathbf u_k)&=m_f\mathbf v_i^{*},\\
m_b\mathbf V_b-\Delta t\,m_f\sum_{(i,k)\in b}c_{ik}(\mathbf v_i-\mathbf u_k)+\Delta t\sum_{\text{links}\ni b}c\,\hat{\mathbf n}\hat{\mathbf n}^{T}(\mathbf V_b-\mathbf V_{b'})&=m_b\mathbf V_b^{*},\\
I_bW_b-\Delta t\,m_f\sum_{(i,k)\in b}c_{ik}\,\mathbf l_k\times(\mathbf v_i-\mathbf u_k)&=I_bW_b^{*},
\end{aligned}
$$

where $\mathbf v^{*}$, $\mathbf V^{*}$, $W^{*}$ are the values after the earlier steps of the sub-step (overlap push, drag, buoyancy), so the solve does not double-count the same relative velocity. The last term of the second equation is the squeeze-film link below. The matrix is symmetric positive definite and is solved by Jacobi-preconditioned CG with the same relative tolerance $10^{-3}$ and up to $4\times50$ iterations (the ball unknowns couple distant fluid through a stiff rigid-body mode), warm-started with the previous sub-step's viscous correction (`pbf.rs:660`-`pbf.rs:790`). Momentum property: a uniform translation of all unknowns is a null mode of every coupling term, so for the converged solution the total linear momentum $m_f\sum\mathbf v_i+m_b\sum\mathbf V_b$ equals its value before the solve; the ball receives exactly the linear impulse the wall term removed from the fluid and the torque $\mathbf l_k\times\mathbf F_k$ of the same forces at the ghost points. With the CG tolerance the violation is bounded by the projection of the final residual on uniform translation. (Angular momentum is exact for ball spin about its own centre by the construction above; the fluid-side force acts at the fluid particle, not the ghost point, so the total angular momentum about the origin is conserved only to that offset.) Treating the ball explicitly was unstable once the slurry acts as a rigid lump (about 10 Pa·s and above). The change in $\mathbf V_b$ and $W_b$ is added to the ball's impulse (`pbf.rs:2283`-`pbf.rs:2291`).

**Sub-grid squeeze-film links (`squeeze_links`, `pbf.rs:543`; model `crates/mill-core/src/mac/lubrication.rs`).** When two balls (or a ball and the circular drum wall) approach to a gap below what the fluid lattice resolves, the resolved flow underestimates the lubrication resistance and the particle fluid cannot be refined to resolve it. Classical lubrication theory for a sphere approaching a plane gives a resistance that diverges like $h^{-1}$ as the gap $h\to0$ [14, 15]; the 2D (cylinder) analogue used here diverges like $h^{-3/2}$ (the classical papers are context for the idea, not the source of the 2D formula); resolved-fluid solvers add a lubrication correction for the unresolved gap [16, 17, 18], and DEM models for wet milling use lubrication as a force law [19, 20]. The squeeze-film coefficient used here is the 2D Reynolds result for two curved surfaces, with correction terms fitted to exact bipolar-coordinate solutions of a disc near the drum wall and of two discs in unbounded fluid (`crates/mill-core/src/mac/lubrication.rs:1`-`lubrication.rs:14`):

$$
c(h)=3\sqrt2\,\pi\,\mu\,R_c^{3/2}\,h^{-3/2}\bigl(1+c_1\varepsilon+c_2\varepsilon^2\bigr),\qquad \varepsilon=\frac h{R_c},\qquad \frac1{R_c}=\kappa_1+\kappa_2,
$$

$$
Q=R_c^3(\kappa_1^3+\kappa_2^3),\qquad c_1=0.788+0.265\,Q,\qquad c_2=\max\bigl(0,\ 0.006+0.057\,(Q-0.25)\bigr)
$$

(`SqueezeFilm::new`, `lubrication.rs:53`; `coefficient`, `lubrication.rs:67`). $\kappa$ is the signed curvature: $+1/a$ for a disc of radius $a$ and $-1/b$ for the concave drum wall (`disc_curvature`, `wall_curvature`, `lubrication.rs:23`, `lubrication.rs:28`). For two equal discs $R_c=a/2$, $Q=1/4$; for a disc against a flat wall $Q=1$. $c$ has units of Pa·s here (force per unit depth per unit approach speed), with $\mu$ the dynamic viscosity (`pbf.rs:558`). The model is used only for gaps below the resolved range, blended with the resolved fluid by the fraction of the film that the lattice already delivers, tabulated from grid measurements (`resolved_ratio`, `lubrication.rs:95`; ratio 1 for gaps of 4 cells and more, 0.90 at 2, 0.25 at 1, 0 at touching). In `squeeze_links`, with $\Delta x_f=h/2$ the fluid spacing and gap $g_{ab}=|x_a-x_b|-2r_b$ (ball-wall: $R-r_b-|x_a|$):

- A link exists for $g_{ab}<8\Delta x_f$ (`H_START = 8`, `lubrication.rs:19`).
- Weight: $s=1-\text{resolved\_ratio}(g_{ab}/\Delta x_f)$, applied to the sub-grid film, $c_\text{link}=\min\bigl(G_s\, s\, c(h_f),\ c_\text{max}\bigr)$ with $h_f=\max(g_{ab},\,0.02\,r_b)$ (`MIN_FILM_GAP_RATIO = 0.02`, `pbf.rs:532`: surface roughness; touching balls would otherwise get an unbounded coefficient), the gain $G_s=1$ (`SQUEEZE_GAIN`, `pbf.rs:523`) and the stiffness cap $c_\text{max}=8\,m_b/\Delta t$ (`SQUEEZE_STIFFNESS_CAP = 8`, `pbf.rs:528`), so that $\Delta t\,c/m_b\le8$ and the joint system stays well conditioned in a dense settled charge.
- Closing-only filter: a link is kept only if the partners approach, $\hat{\mathbf n}\cdot(\mathbf V_a-\mathbf V_b)>0$ at the start of the solve (a film cavitates when the gap opens and resists squeezing only; a symmetric film kept a settled charge jittering). The active set is fixed during the solve so the system stays linear.
- The film force on ball $a$ is $-c_\text{link}\,[\hat{\mathbf n}\cdot(\mathbf V_a-\mathbf V_b)]\hat{\mathbf n}$ (opposite on $b$), entering the ball rows of the joint system as the $\Delta t\,c\,\hat{\mathbf n}\hat{\mathbf n}^T$ term. Only the circular wall is modelled, not the lifters, and only the normal film (the tangential coefficient in `lubrication.rs:79` is used by the grid track only).

**Stability clamp on the impulse (step 8, `pbf.rs:2317`).** The accumulated impulse on each ball is clamped to

$$
|J_b|\le F_g\,m_b\,g\,\Delta t,\qquad |\Gamma_b|\le F_g\,m_b\,g\,\Delta t\,r_b,\qquad F_g=\texttt{F\_CLAMP\_G\_MULTIPLE}=60
$$

(`impulse_clamp`, `pbf.rs:958`; `pbf.rs:956`). It is a backstop: the fluid-side reaction has already been applied unclamped, so a clamp event is a deliberate, counted (`coupling_clamp_hits`) conservation violation; in a settled charge it fires on a few percent of ball-sub-steps at most, and persistent clamping indicates a real problem.

### 2.6 Free surface for rendering

Source: `crates/mill-core/src/surface.rs` (visualisation only; it does not feed back into the physics).

1. **Splat** (`build_field`, `surface.rs:84`). Every particle adds the compactly supported kernel $K(r)=(1-r^2/\rho^2)^3$ ($r<\rho$) with radius $\rho=1.5\,h$ (`SPLAT_RADIUS_FACTOR`, `surface.rs:29`; `splat_kernel`, `surface.rs:39`) to a $128\times128$ scalar field spanning $[-R,R]^2$ (`GRID_SIZE`, `surface.rs:27`).
2. **Mask.** Grid points outside the wall or inside a lifter (SDF $<0$) are set to zero.
3. **Threshold.** The iso-level is $0.5\,\phi_\text{full}$ (`THRESHOLD_FRACTION`, `surface.rs:31`), where the bulk field value is computed analytically, not read off the grid (`analytic_phi_full`, `surface.rs:71`):

$$
\phi_\text{full}=\frac{\rho_0}{m_f}\cdot\frac{\pi\rho^2}{4},
$$

   since the kernel integrates to $\int_0^\rho(1-r^2/\rho^2)^3 2\pi r\,dr=\pi\rho^2/4$ and the hexagonal lattice has number density $\rho_0/m_f$. By symmetry the field value at a flat bulk/free-surface boundary is half the bulk value. A measured-maximum fallback is used if the analytic threshold yields no contour (`extract_surface`, `surface.rs:233`).
4. **Contour.** Marching squares, the 2D analogue of marching cubes [21], with segments chained into polylines by exact edge identifiers (`polygons_at_threshold`, `surface.rs:279`).
5. **Smooth and clamp.** One Chaikin corner-cutting pass [22] (new points at $0.25$ and $0.75$ along each edge, `chaikin_smooth_closed`, `surface.rs:403`), then contour points still outside the wall are projected back along the SDF normal (`project_inside_wall`, `surface.rs:359`).

### 2.7 Metrics and the 2D-to-3D grinding proxy

Sources: `crates/mill-core/src/metrics.rs`, `crates/mill-core/src/stress.rs`, `crates/mill-core/src/lib.rs`.

**Charge shape** (`charge_toe_shoulder`, `metrics.rs:327`). Balls within $2.5\,r_b$ of the wall by the drum SDF (so lifter-borne balls count) are collected (`WALL_MARGIN_BALL_RADII`, `metrics.rs:44`); if there are at least 6 (`MIN_WALL_BALLS`) and the largest angular gap between them exceeds $0.3\pi$ (`CENTRIFUGE_GAP_THRESHOLD_RAD`, `metrics.rs:52`), the two ends of the cluster on the far side of that gap are toe and shoulder, assigned by the sign of $\omega$; otherwise the charge is centrifuged and both are undefined. Angles are $\operatorname{atan2}(y,x)\in[0,2\pi)$; `to_vertical_degrees` (`metrics.rs:275`) converts to the mill-literature convention (0 at 12 o'clock, increasing in the direction of rotation).

**Lacey mixing index** [23] (`mixing_index`, `metrics.rs:505`). Fluid particles are binned into a $16\times16$ grid over the drum bounding box (`MIXING_GRID_SIZE`, `metrics.rs:62`) with the dye value $c\in\{0,1\}$ as the tracer:

$$
M=\frac{S_0^2-S^2}{S_0^2-S_r^2},\qquad S_0^2=\bar p(1-\bar p),\qquad S_r^2=\frac{S_0^2}{\bar n},
$$

$$
S^2=\frac{\sum_c n_c(\bar c_c-\bar p)^2}{\sum_c n_c},\qquad \bar n=\frac{\sum_c n_c^2}{\sum_c n_c},
$$

where $c$ runs over occupied cells, $n_c$ is the particle count and $\bar c_c$ the mean dye of cell $c$, $\bar p$ the overall mean dye fraction; $S^2$ is particle-count-weighted so a sparse splash cell does not dominate, and $\bar n$ is the matching weighted occupancy. The result is clamped to $[0,1]$; if the denominator is below $10^{-9}$ it returns 1. Other readouts: slurry pool angular extent, free-surface line by a 2×2 covariance principal axis, pool depth at the bottom within $\pm5^\circ$ of straight down (`metrics.rs:377`-`metrics.rs:477`), maximum ball-ball and ball-wall overlap fractions (relative to the ball radius), and the maximum and mean fluid density and compression errors.

**Power and torque** (`update_grinding_stats`, `lib.rs:516`). The motor power per metre of mill is the work the wall delivers to balls and fluid per sub-step, smoothed by an exponential moving average with $\alpha=1-e^{-\Delta t/\tau}$, $\tau=1$ s (`GRINDING_STATS_EMA_TAU_S`, `lib.rs:48`):

$$
P_\text{inst}=\frac{W^\text{DEM}_\text{wall}+W^\text{fluid}_\text{wall}}{\Delta t},\qquad T=\frac{P}{\omega}\ \ (T=0\text{ for }|\omega|\le10^{-6}).
$$

$W^\text{DEM}_\text{wall}$ is the ball-wall friction and lifter normal work (Section 2.2); $W^\text{fluid}_\text{wall}$ is $\sum m_f\,\Delta\mathbf v\cdot\mathbf v_\text{wall}$ over wall-flagged fluid particles in the no-slip blend. The net dissipated power $\max(E_\text{diss},0)/\Delta t$ (DEM side only) and the impact rate are smoothed the same way.

**Energy budget** (`EnergyBudget`, `lib.rs:61`; accumulated at `lib.rs:492`-`lib.rs:504`). A cumulative read-only account of the whole system's mechanical energy: shaft work in, ball-contact dissipation, fluid wall slip, fluid viscous dissipation, speed-clamp removal and the net energy created at the ball-fluid interface (drag, buoyancy and coupling work), with the unattributed remainder $\Delta E-(W_\text{shaft}-D_\text{contact}-D_\text{slip}-D_\text{visc}-D_\text{clamp}+C_\text{interface})$ reported as the closure error (`Simulation::energy_budget`, `lib.rs:218`).

**3D-converted stress statistics** (`stress.rs`; module comment `stress.rs:1`-`stress.rs:21`). Grinding itself is not simulated; what is reported are quantities that mechanistic grinding models are driven by: the distribution of stress energies and the frequency of stressing events, the quantities introduced as stress energy and stress number in the stirred-media-mill literature [24, 25]. The 2D disc mass overstates stress energy by a size-dependent factor (about 750 for 2 mm balls), so energies use the mass of the 3D sphere that a disc stands for (`Simulation::stress_scales`, `lib.rs:233`):

$$
m_3=\rho\,\tfrac43\pi r^3=m_\text{disc}\cdot\tfrac43 r,\qquad \text{SE}=\tfrac12 m_\text{red}\bigl(v_n^2+v_s^2\bigr),\qquad m_\text{red}=\begin{cases}m_3/2&\text{ball-ball}\\ m_3&\text{ball-wall}\end{cases}
$$

where $v_n$ is the pre-solve normal approach speed and $v_s$ the tangential slip speed at the contact point including spin (`classify`, `stress.rs:68`). An impact is of drop type if $v_n\ge v_s$ and of shear type otherwise, giving four classes (ball-ball and ball-wall, drop and shear; `CLASS_NAMES`, `stress.rs:33`). Energies are made dimensionless by $m_3 g d$: $\hat E=\text{SE}/(m_3 g d)$, binned in 28 log-spaced bins from $10^{-3}$ upward at 4 bins per decade (`stress.rs:28`-`stress.rs:30`). Re-contacts of the same pair closer in time than $0.25\sqrt{d/g}$ count as one episode (`DEBOUNCE_FACTOR`, `stress.rs:48`), so the count does not grow as the sub-step shrinks. Reported: impacts per ball per second per class, the mean $\hat E$ (the stress-intensity measure, mean of $\text{SE}/(m_3 g d)$ over all impacts), and the rate and cumulative energy carried by impacts above thresholds $\hat E\in\{0.03,0.1,0.3,1,3,10\}$ (`THRESHOLDS_REL`, `stress.rs:42`); the time two conditions need to reach the same cumulative effective energy is a relative grinding-time proxy.

**Per-ball power and the depth basis.** "Per ball" power uses a one-ball-thick slab: $P_\text{ball}=P_{2D}\,d/N$, specific power $P_\text{ball}/m_3$, cumulative specific energy likewise (`StressAccum::metrics`, `stress.rs:186`). This is equivalent to an implied 3D solid fraction: a sphere of volume $\tfrac43\pi r^3$ is represented by a cylinder of volume $\pi r^2 d=2\pi r^3$, hence

$$
\phi_{3D}=\tfrac23\,\phi_{2D}\quad(0.547\text{ at the default }\phi_{2D}=0.82).
$$

Results are meant for relative comparison between conditions run with the same code and settings. With coarse-graining ($k>1$) the simulated ball is not the real ball and these statistics must not be used.

### 2.8 Research tracks not used by the shipped simulation

Two further fluid solvers are in the repository as reference and verification instruments. `Simulation` does not call either (verified by search; the only cross-use is `pbf.rs` importing the squeeze-film model from `mac`).

**DFSPH with smooth solids (`crates/mill-core/src/fluid/`).** A divergence-free SPH solver [27] intended as a replacement for PBF (`fluid/mod.rs:1`). Fluid density is $\rho_i=m\sum_j W_{ij}+\rho_0\sum_S\Gamma_S$, with the 2D cubic-spline kernel (support $H$, $\sigma=40/(7\pi H^2)$, `fluid/kernel.rs:6`). The drum wall, lifters and balls are not particle layers but smooth solids represented by tabulated kernel integrals in the sense of density maps [28]: $\Gamma(x)=\int_S W(|x-x'|)\,dA'$ (solid volume fraction seen by the kernel) and the viscous coupling tensor $T(x)=\int_S K(|x-x'|)\,\hat e\hat e^T\,dA'$ with $K=|W'|r/(r^2+\eta^2)$, stored in 1-D tables per body radius (`fluid/solid.rs:1`). A constant-density pressure solve gives pressure acceleration $a_i=-\sum_j m(\kappa_i+\kappa_j)\nabla W_{ij}-\kappa_i\rho_0\nabla\Gamma$ with $\kappa=p/\rho^2$ accumulated and clamped at $\ge0$ (no suction at a free surface), iterated Jacobi-style (`fluid/pressure.rs:1`). Viscosity is implicit with pairwise-central coefficients $C_{ij}=8\mu m^2/(\rho_i\rho_j)\,|W'|r/(r^2+\eta^2)$ and the ball velocities as unknowns of the same CG system (`fluid/viscosity.rs:1`). It conserves linear and angular momentum exactly and was built to fix PBF's lack of solid-wall density, energy created or destroyed by position-based velocity reconstruction and wall-touching-only no-slip. It currently serves as a comparison solver; it did not reach its verification gates (`docs/VERIFICATION.md`, "Smooth solid boundary").

**MAC cut-cell Navier-Stokes (`crates/mill-core/src/mac/`).** A fixed-grid incompressible solver in the marker-and-cell tradition [29] with a staggered layout (face-normal velocities as unknowns, `mac/staggered.rs:1`) and cut cells for the circular drum and for discs (face apertures, `mac/embedded.rs:1`). Its components are:
- Variational (symmetric) discretisation of the viscous and yield stress, where the stress force is exactly the transpose of the strain-rate operator so the implicit operator and the explicit stress terms are adjoint [30] (`mac/variational.rs:1`); Dirichlet Helmholtz on irregular domains by the symmetric second-order scheme of [31] (`mac/viscous.rs:1`).
- Aperture-weighted Neumann pressure Poisson system solved by a multigrid-preconditioned CG, with Galerkin coarse operators from 2×2 aggregation (`mac/multigrid.rs:1`).
- A cell-centred level set for the free surface advected with WENO5 upwind derivatives [32] and TVD-RK3 time stepping, re-initialised to a signed distance by the PDE method with the sub-cell fix of [33] (`mac/levelset.rs:1`).
- Rigid discs in the fluid with a Newton iteration on all body velocities at once (resolved pressure plus wall traction, contact, lubrication and friction; `mac/bodies.rs:995`, `mac/newton.rs:1`), and the squeeze-film model (`mac/lubrication.rs`).
- Yield-stress fluids by the augmented-Lagrangian method (ALG2) [34] for a Herschel-Bulkley law [35], $\tau=\tau_y+k\dot\gamma^n$ (`mac/rheology.rs:1`).

**Role today.** The grid track is the reference solution against which the shipped PBF is calibrated and checked: the ball no-slip gain `BALL_BC_GAIN` (`pbf.rs:486`) and the squeeze-film gain and fit (`mac/lubrication.rs:1`; `SQUEEZE_GAIN`, `pbf.rs:523`) were tuned on it, and the verification gates (Taylor-Couette torque, disc settling, drum torque balance, a dense disc pile, Bingham paste) are recorded in `docs/VERIFICATION.md`. They are not part of the browser build's step loop.

### 2.9 Determinism

Every run is reproducible from `Params` and `seed`. Three mechanisms ensure that the floating-point summation order is identical run to run:

- **Counting-sort CSR grid** (`UniformGrid`, `crates/mill-core/src/grid.rs:79`; `build`, `grid.rs:101`). Points are bucketed by cell with a two-pass counting sort (histogram, prefix sum, stable fill) into flat arrays, so candidate pairs and neighbour lists are enumerated in a fixed $(c_x,c_y)$-lexicographic cell order and point-index order within a cell, with no hash map whose iteration order depends on a per-process random seed. A sparse fallback with a stable sort gives the identical order for pathologically spread bounding boxes. The ball solver sorts the contact list by $(i,j)$ before the passes that mutate shared per-ball state (`dem.rs:459`). Non-finite points are skipped.
- **Fixed iteration order and fixed neighbour order.** Gauss-Seidel loops run over balls and contacts in index order, and the fluid's CSR neighbour lists preserve the order in which pairs are enumerated (`build_neighbor_lists`, `pbf.rs:261`).
- **Seeded PRNG** (`Rng`, `rng.rs:9`). The state is advanced by a Marsaglia xorshift64 step (shifts 13, 7, 17: `x ^= x << 13; x ^= x >> 7; x ^= x << 17`, `rng.rs:27-31`) and the returned value is the new state multiplied by `0x2545F4914F6CDD1D` (`rng.rs:33`). It supplies the lattice jitter and is seeded from `simulation.seed` (0 is remapped to a fixed non-zero constant). The code comments call it xorshift64*; the multiplier is the one of the xorshift64* generator analysed in [26], whose published shift triple is (12, 25, 27), not the (13, 7, 17) used here. The generator serves only non-critical jitter; its statistical quality is not relied upon.

---

## 3. Software technology

This section describes how the model of Sections 1-2 is turned into a program that runs unmodified in a web browser, how its speed and reproducibility are engineered, and how the surrounding experiment and distribution tooling is built. All source locations refer to the repository at the commit this document was written against; where a project document and the source disagree, the source is taken as correct.

### 3.1 Workspace layout and the Rust/WASM boundary

The repository is a Cargo workspace with two crates (`Cargo.toml:1-3`), a Vite/TypeScript front end in `web/`, and a standard-library-only Python launcher in `python/`.

- `crates/mill-core` is a plain Rust library. Its manifest declares only `glam` (with the `libm` feature and default features off), `serde` and `serde_json` as dependencies (`crates/mill-core/Cargo.toml:8-11`) and describes itself as having "no WASM dependencies" (`crates/mill-core/Cargo.toml:6`). Everything physical lives here: the ball solver (`dem.rs`), the slurry solver (`pbf.rs`), the ball-fluid coupling (`coupling.rs` and the coupling steps of `pbf.rs`), geometry (`geometry.rs`), spatial grid (`grid.rs`), free-surface extraction (`surface.rs`), metrics (`metrics.rs`), the grinding-proxy statistics (`stress.rs`), the headless seed runner (`batch.rs`) and the deterministic random number generator (`rng.rs`) (module list at `crates/mill-core/src/lib.rs:18-31`). Because the crate has no browser dependency, the same code is compiled natively for unit tests, the verification suite, benchmarks and the diagnostic probes of Section 4.
- `crates/mill-wasm` is a thin `wasm-bindgen` wrapper (`crates/mill-wasm/src/lib.rs:1-5`) built as `cdylib` and `rlib` (`crates/mill-wasm/Cargo.toml:9`). It owns all JavaScript-facing naming, JSON marshalling and the panic hook.
- `crates/mill-core/src/fluid/` (a DFSPH solver) and `crates/mill-core/src/mac/` (a fixed-grid Navier-Stokes solver with resolved ball coupling) are experimental solvers used only for verification and as calibration references (Section 5.2, item 3). The shipped `Simulation` runs the ball solver and the position-based slurry solver only (`docs/PHYSICS.md:93-100`, consistent with the module use in `crates/mill-core/src/lib.rs`).

**Exported API.** The class `Simulation` exposes (`crates/mill-wasm/src/lib.rs:15-167`): a constructor taking an optional JSON `Params` string (`:24-34`); `setParams(json)` (`:39-46`); `step(dt)` (`:49`); `fixedSubDt()` and `stepFixed()` for the fixed-sub-step accumulator (`:55-66`); `resetFrameStats()` (`:71`); `drumAngle()`, `simTime()`, `ballCount()`, `ballRadiusM()` (`:76-98`); the array accessors `ballPositions()`, `ballOrientations()`, `fluidPositions()`, `fluidDye()` and `fluidSurface()` (`:103-150`); and `paramsJson()` and `metricsJson()` (`:156-166`). Free functions are `defaultParamsJson()` (`:171`), `setPanicHook()` (`:178`) and `stressRunJson(paramsJson, seed, warmupRevolutions, measureRevolutions)` (`:187-204`), which runs the same `Simulation` headlessly for one seed (Section 3.5).

**Data crossing.** Parameters and metrics cross the boundary as JSON strings (`serde_json`). Particle arrays are returned as `Vec<f32>`, which `wasm-bindgen` marshals into a fresh `Float32Array` per call: the data are copied, not shared with WASM memory (`crates/mill-wasm/src/lib.rs:7-9`, `:100-102`). The worker then moves each copy to the main thread as a transferable buffer (Section 3.4), so there is one copy out of WASM memory per array per frame and no copy on the worker-to-main hop. A zero-copy view is noted in the source as a possible future change; no profile in `docs/PERF.md` attributes a measured cost to the copy.

**Build configuration.** The release profile is `opt-level = 3`, `lto = true`, `codegen-units = 1`, `panic = "abort"` (`Cargo.toml:10-14`). For the `wasm32-unknown-unknown` target only, `.cargo/config.toml:10-11` sets `-C target-feature=+simd128`; native builds are unaffected. The file's own comment states that this changes which instructions the compiler may emit, not floating-point semantics, and that WASM SIMD is shipped by the major browsers (`.cargo/config.toml:1-9`). `docs/PERF.md` records that SIMD was kept because it is free and non-regressive and "not because it was isolated and measured on its own". `scripts/build-wasm.sh:8-13` runs `wasm-pack build crates/mill-wasm --target web --release` and writes the bindings to the git-ignored `web/src/wasm/`. An optional Cargo feature `stage-timing` (`crates/mill-core/Cargo.toml:14-15`) compiles per-stage wall-clock marks into the solver for native profiling; with the feature off the `stage_mark!` macro expands to nothing (`crates/mill-core/src/stage_timing.rs:1-6`).

### 3.2 Performance engineering

All performance changes recorded below were made under a constraint stated in `docs/PERF.md`: a pass meant to be "mechanical" must reproduce a recorded state hash bit for bit (Section 3.3), and approximate speed-ups must pass a physical-metric tolerance gate (Section 3.3). The speed-up figures are quoted exactly as recorded in `docs/PERF.md`; they are same-session A/B measurements (`git stash` or interleaved builds) on the author's machine, and the document itself warns that absolute milliseconds vary by roughly 10 % between runs and by up to a factor of 10 between sessions (`docs/PERF.md:289-299`, `:420-422`). Only the relative changes should be read.

1. **Counting-sort CSR spatial grid** (`crates/mill-core/src/grid.rs`). The earlier `BTreeMap<(i32,i32),Vec<u32>>` grid, rebuilt four times per sub-step, was replaced by a compressed-sparse-row structure built with a histogram, a prefix sum and a stable fill (`grid.rs:155-198`). The `Dense` layout addresses cell $(c_x,c_y)$ by the linear index $k=(c_x-c_{x,\min})\,H+(c_y-c_{y,\min})$ with $H$ the number of cell rows, giving $O(1)$ lookup and one flat `Vec<u32>` for the whole point set. A `Sparse` layout (sorted keys, binary search) is used when the bounding box is large relative to the point count (`grid.rs:44-48`, `:101-153`, `:204-223`). Both layouts visit cells in row-major order and points in ascending index order, which reproduces the iteration order of the earlier map exactly (`grid.rs:83-84`). The dense-cell budget constant was tuned from 8 to 2 points per cell after a measurement; `docs/PERF.md` records that this still leaves the isolated ball-only micro-benchmark 7-19 % slower at 1000-2000 balls, which was accepted because the full coupled step improved (`docs/PERF.md:273-287`).
2. **CSR neighbour lists and a fused density pass** (`crates/mill-core/src/pbf.rs:228-258`, `:1244-1262`). Fluid neighbour lists are one flat array with offsets rather than one `Vec` per particle, built from a directed edge list by a stable counting sort so that each particle sees its neighbours in the same order as before. The density and Lagrange-multiplier loops are fused into one pass over each particle's neighbours, and each directed edge's raw kernel gradient is cached in a parallel buffer and reused by the position-correction pass (`pbf.rs:1259`, `:1278`, `:1318`).
3. **Hoisted scratch buffers.** Per-particle and per-ball scratch vectors in the overlap, cohesion, drag and buoyancy steps were hoisted out of inner loops and cleared per use, and the conjugate-gradient viscosity solve preallocates its working buffers once per call (`docs/PERF.md:230-237`, `:250-257`). This removes allocations from the innermost loops; it does not make a sub-step allocation-free, since per-sub-step buffers such as `density`, `lambda`, `delta_p` and `grad_cache` are still allocated (`pbf.rs:1256-1259`). The `stage-timing` build installs a counting allocator in the probe to measure allocation counts (`crates/mill-core/examples/perf_probe.rs:24-52`).
4. **Merged ball-boundary groups and a warm-started viscosity solve** (M7 pass, `docs/PERF.md:360-400`). Ghost-point pairs between one fluid particle and one ball are merged into a single `BoundaryGroup` that stores $c=\sum_k c_k$, $\sum_k c_k\boldsymbol l_k$ and $\sum_k c_k|\boldsymbol l_k|^2$ (`pbf.rs:491-505`); the document records about 20.7k pairs reducing to about 6.5k groups for the same linear operator. The conjugate-gradient iteration is warm-started from the previous sub-step's viscous velocity correction (`pbf.rs:1018-1023`), which `docs/PERF.md:384` records as lowering the mean iteration count at 50 Pa s from 102.5 to 81.4. The document also lists approximate changes that were tried and rejected (a looser CG tolerance, a shorter squeeze-link range, fewer adaptive density iterations) because they failed the tolerance gate (`docs/PERF.md:402-412`).

Recorded results (all from `docs/PERF.md`, same session, author's machine):

| Pass | Measurement | Before | After | Change |
|---|---|---|---|---|
| M6 (grid, CSR, hoisting), native `drum_only_step`, default wet parameters | ms per step | 379.7 | 216.5 | -43 % (`docs/PERF.md:309`) |
| M6, browser, Balanced tier | achieved speed | ~0.04-0.06x | ~0.11-0.12x | roughly 2x (`docs/PERF.md:320`) |
| M7, `ship_fast`, water | ms per frame | 65.09 | 49.39 | -24.1 % (`docs/PERF.md:395`) |
| M7, `ship_fast`, 50 Pa s | ms per frame | 175.57 | 94.92 | -45.9 % (`docs/PERF.md:397`) |
| M7, `ship_acc`, 50 Pa s | ms per frame | 498.60 | 291.66 | -41.5 % (`docs/PERF.md:400`) |

The browser-speed figures of `docs/PERF.md` predate later model changes (water default, central-force viscosity, squeeze films, higher sub-step counts) and are described there as stale (`docs/PERF.md:350-358`); the browser speed of the current quality tiers has not been measured (`docs/PERF.md:356`). The remaining lever identified in the document is sub-cycling the fluid solve, because about 90 % of a frame is fluid work while the sub-step count is dictated by ball-solver stability; this is a physics change that has not been attempted (`docs/PERF.md:414-418`).

### 3.3 Determinism and reproducibility

A run is a function of `Params` (which includes `simulation.seed`) and nothing else. The design choices that support this are:

- **Seeded generator.** `rng.rs` implements a small xorshift-with-multiplier generator (Section 2.9; the xorshift family, cf. [26]) instead of depending on `rand`, "so that the simulation is bit-reproducible across platforms/builds given the same seed" (`crates/mill-core/src/rng.rs:1-5`). A zero seed is remapped to a fixed constant (`rng.rs:15-24`). The statement of cross-platform bit-reproducibility is a design intent; the repository verifies bit-exactness between two builds of the same target (see below), not between native and WASM executions.
- **Fixed sub-step.** The sub-step is $\Delta t = 1/(60\,n_\text{sub})$ with $n_\text{sub}$ the effective sub-step count, independent of wall-clock time and of `time_scale` (`crates/mill-core/src/lib.rs:426-428`; the front-end loop is described in Section 3.4).
- **Ordered iteration.** Containers with unspecified iteration order are avoided in the solver; the grid, the neighbour lists and the contact sets are built so that cells, particles and pairs are visited in a fixed order, and the Gauss-Seidel contact passes sort contacts by pair index (`grid.rs:83-84`, `pbf.rs:241-258`, `docs/PERF.md:257-260`). Accumulations over neighbours sum in the same order as before every performance change.
- **Hash oracle.** `crates/mill-core/examples/perf_probe.rs` is a standalone probe, not a test. For each named configuration (quality presets, a dry run, a lifter run, and the shipped-default configurations `ship_fast`, `ship_bal`, `ship_acc` with 63 mm drum and 2 mm balls) it settles the charge, times a fixed number of `Simulation::step(1/60)` calls, and prints a 64-bit FNV-1a hash over the raw `f32` bits of every ball $x$, $v$, $\omega$ (and, in the extended variant, $\theta$) and every fluid particle $x$, $v$ (and dye), taken in a fixed field and index order (`perf_probe.rs:222-262`; header `perf_probe.rs:1-22`; shipped configurations `docs/PERF.md:362-367`). A refactor that is meant to change nothing physical must reproduce the identical hash for every configuration; the M6 pass did so for four successive changes (`docs/PERF.md:262-271`). The hash is sensitive to any change in floating-point summation order, not only to a change of the physical outcome (`perf_probe.rs:222-230`).
- **Tolerance gate for approximate changes.** `crates/mill-core/examples/tolerance_probe.rs` compares physical metrics between two builds over seeds $1..N$ after 2 settle and 3 measured revolutions (`tolerance_probe.rs:1-20`). Gates (`tolerance_probe.rs:11-27`, `docs/VERIFICATION.md:1672-1690`): power, torque, specific power and impact rate pass if $|\Delta|\le\max\!\big(2\sqrt{\mathrm{se}_a^2+\mathrm{se}_b^2},\,2\,\%\,|\bar x_\text{base}|\big)$; toe and shoulder angles use the same rule with a 2 degree floor; the unattributed-energy fraction passes if $|f_\text{new}|\le 1.10\,|f_\text{base}|+0.005$ (relative, because the reference solver's own closure is already off by 0.17 to 0.38 depending on configuration); mean and worst-5 % compression pass if $x_\text{new}\le 1.10\,x_\text{base}+10^{-4}$. A methodological note recorded with the M7 result: three seeds underestimated the noise (apparent +1.3 to +2.4 % power shifts vanished at eight seeds), so at least eight seeds are used near the 2 % floor (`docs/VERIFICATION.md:1706-1708`).

### 3.4 Web runtime

The front end is TypeScript on Vite with no UI framework. The simulation never runs on the main thread.

**Threads and protocol.** `web/src/main.ts` creates one module Web Worker running `web/src/worker.ts` (`main.ts:35`) and exchanges typed messages defined in `web/src/protocol.ts`: from main to worker `init`, `setParams`, `play`, `pause`, `step`, `requestFrame`; from worker to main `ready`, `frame`, `frameSkipped`, `error` (`protocol.ts:13-125`). A second, independent worker (`web/src/batch/batchWorker.ts`) runs the headless seed runs of Section 3.5, so a blocking multi-minute WASM call never shares a thread with the interactive simulation (`batchWorker.ts:1-4`).

**Fixed sub-step accumulator with a wall-clock budget.** The main thread's `requestAnimationFrame` loop (`main.ts:290-319`) asks for a frame; the worker keeps an accumulator $A$ of simulated time owed (`worker.ts:44-47`). On each `requestFrame` the wall-clock delta is clamped, $\Delta t_\text{wall}\leftarrow\min(\Delta t_\text{wall},1/15\,\mathrm{s})$ (`worker.ts:49-53`, `:256`), and $A \leftarrow A + s\,\Delta t_\text{wall}$ with $s$ the requested time scale (`worker.ts:257`). `drainPendingSimTime` then calls `stepFixed()` while $A\ge\Delta t_\text{sub}$, subtracting $\Delta t_\text{sub}$ each time, and stops when the budget `frame_budget_ms` (default 12 ms, `crates/mill-core/src/params.rs:492`) of wall time is exhausted, checked once per sub-step (`worker.ts:137-146`). If $A$ still exceeds $\max(8\,\Delta t_\text{sub},\,1\,\mathrm{s})$ the remainder is discarded, so a stretch of slow frames cannot later produce a destabilising catch-up burst (`worker.ts:147-150`). The reply carries a measured achieved time scale $=\text{steps run}\cdot\Delta t_\text{sub}/\Delta t_\text{wall}$ (`worker.ts:259-262`), which the head-up display shows. Because only the number of sub-steps per frame varies and never their size, a slow machine produces slow motion, not a larger, less stable step (`worker.ts:8-18`).

**Credit-based frame requests.** Left alone, a faster main thread would queue unbounded `requestFrame` messages behind a slow worker and starve `pause` and `setParams`. `FrameRequestGate` allows at most `DEFAULT_MAX_IN_FLIGHT = 2` unanswered requests (`web/src/frameGate.ts:9`, `:26-34`); `main.ts:296-301` sends a request only when the gate issues a credit, and `main.ts:93-97` settles the credit when a `frame` or `frameSkipped` reply arrives. The gate's clock is reset on pause and on reset so the paused interval is not folded into the next `wallDt` (`frameGate.ts:46-52`, `main.ts:108`, `:127`).

**Transferable buffers and throttled extras.** Each `frame` message carries the four particle arrays plus, on throttled frames, the free-surface array; their underlying `ArrayBuffer`s are listed in the transfer list so that ownership moves to the main thread without a further copy (`worker.ts:101-117`, transfer list at `:109-116`). The free-surface contour (marching squares over a regular grid, see `crates/mill-core/src/surface.rs:5-12`) and the metrics JSON are recomputed at most every $1000/15\approx 67$ ms (`SLOW_UPDATE_INTERVAL_MS`, `worker.ts:59-69`, applied at `:93-97`); the main thread keeps the last received value on frames where they are absent (`main.ts:76-80`). The project does not use `SharedArrayBuffer`, WASM threads or `OffscreenCanvas`: a search of `web/src` finds none of them, and `docs/PLAN.md:395` records that `SharedArrayBuffer` and WASM threads were deliberately not required, because they would need cross-origin isolation (COOP/COEP) headers on the development server, and lists them as future work.

**Hot-apply versus reset.** Editing a parameter does not always restart the run. `paramsChangeRequiresReset(old, new)` walks the schema and returns the label of the first changed field that is marked `resetRequired` (`web/src/params/schema.ts:184-192`). Fields baked into the seeded populations (drum and media geometry, `simulation.max_balls`, `resolution`, `coarse_graining_*`, `seed`, initial tracer colouring) are marked as requiring a reset (`schema.ts:143`, `:152-171`); `substeps`, `dem_iterations` and `time_scale` are read fresh every sub-step and are hot-applied (`schema.ts:168-170`). On Apply, `main.ts:123-138` either sends `init` (new `Simulation`) or `setParams` (the worker calls `Simulation::set_params`, keeping drum angle, time and populations, `worker.ts:218-227`). The worker caches `time_scale` and `frame_budget_ms` after both paths (`worker.ts:159-163`).

### 3.5 Rendering and user interface

**Rendering.** `web/src/render/canvas.ts` draws with the Canvas 2D API only: the drum wall and a radial rotation marker fixed to the wall (`canvas.ts:65-88`), lifters as trapezoids rotating with the drum if present (`canvas.ts:90-110`), the slurry free surface as translucent polygons from the flattened contour array `[n_polys, len_0, x, y, ...]` (`canvas.ts:136-142`, array layout `crates/mill-core/src/surface.rs:420`), fluid particles coloured by a dye tracer (`canvas.ts:168`), and balls with a spin indicator showing their orientation (`canvas.ts` ball loop, using `ballOrientations`). The surface extraction uses marching squares [21] followed by one Chaikin corner-cutting pass [22] (`surface.rs:5-12`). World coordinates are $y$-up and are flipped for the canvas. Rendering is decoupled from the solver: the canvas redraws every animation frame from the most recent state, whatever fraction of real time the solver achieved.

**Parameters panel.** `web/src/ui/paramsPanel.ts` is a persistent left-side panel generated from the declarative schema in `web/src/params/schema.ts`, which defines each field's path, group, label, type, limits, whether it requires a reset, and a help string. Groups are Mill, Media, Slurry, Lifters and Numerical accuracy. Units are converted at the UI boundary (mm, rpm; the core is SI). Display-only mirrors of core formulas live in `web/src/params/derived.ts` (critical speed $N_c=42.3/\sqrt{D}$ rpm, with $D$ in metres, true ball count, effective resolution and sub-steps, fluid-spacing-to-ball-diameter ratio) and are labelled as mirrors of `params.rs`, which remains authoritative (`derived.ts:1-8`, `:15-17`). The panel shows the auto-raised effective resolution and sub-steps alongside the requested values, so that the automatic adjustments of Section 5.2 item 4 are never silent.

**Quality tiers.** `web/src/params/presets.ts` defines three tiers, "Fast" (default), "Balanced" and "Accuracy", as (max balls, fluid resolution) pairs 300/25, 300/30 and 600/50, plus an "Impact statistics (slow)" preset (`presets.ts` `QUALITY_PRESETS`). The note strings quote native milliseconds per frame and the power error against a converged reference; they state that the tiers' power errors are +4.6 %, +2.5 % and -1.0 % (`docs/PERF.md:358`, `docs/VERIFICATION.md:1344`).

**Metrics panel and HUD.** `web/src/ui/metricsPanel.ts` presents the metric groups Drum, Media, Grinding, a 3D-basis proxy group, Slurry and Solver (`web/src/metrics/specs.ts:390`), with sparklines, an impact-energy histogram, CSV export and the report controls of Section 3.6. The speed badge in `web/src/ui/hud.ts` is fed through `SpeedSmoother` (`web/src/ui/speedSmoother.ts`): an exponential moving average with time constant 1000 ms, $\alpha = 1-\exp(-\Delta t/\tau)$ (`speedSmoother.ts:53`), published at most every 500 ms, and a hysteresis switch between "REAL TIME" and "SLOW MOTION" that enters slow motion below 0.98 and returns above 0.99 of real-time speed (`speedSmoother.ts:37-40`, `:58-66`). This avoids a flickering badge, since the raw achieved rate is an integer number of fixed sub-steps per frame.

### 3.6 Experiment tooling

The experiment tooling is built into the application so that a result and its provenance leave the browser together.

**Batch runner.** `web/src/batch/runner.ts` owns the lazily created batch worker and runs one seed per request; cancellation terminates the worker because "a blocking wasm call cannot be interrupted any other way" (`runner.ts:1-3`, `:49-52`, `:59-66`). The worker calls `stressRunJson`, which runs `mill_core::batch::run_seed`: the given `Params` with only `simulation.seed` replaced, a warm-up of $w$ revolutions, then frame-by-frame averaging of the proxy statistics over $m$ whole revolutions (the charge oscillates once per revolution, so the window must be a whole number, `crates/mill-core/src/batch.rs:36-39`, `:40-120`). Cumulative energies are used for the window-mean specific power and effective powers so the window means are exact rather than smoothed (`batch.rs:79-93`). The plan layer enforces at most 50 seeds (`web/src/batch/plan.ts:13`), uses the literal seeds $1..N$ and never the main simulation's seed (`plan.ts:28-31`), and refuses to start if coarse-graining is active ($k>1$, because the ball is then not the real ball) or the drum speed is zero (`plan.ts:79-92`). A duration estimate is $T\approx (\text{runs})\cdot(w+m)\,T_\text{rev}\cdot 60\cdot n_\text{sub}/\dot n_\text{achieved}$ with $T_\text{rev}=60/\text{rpm}$ and $\dot n_\text{achieved}$ the main simulation's achieved sub-steps per second (`plan.ts:96-107`).

**Statistics.** For each scalar the per-seed values $x_1..x_n$ give the mean $\bar x$ and the standard error of the mean
$$\mathrm{se}=\sqrt{\frac{1}{n(n-1)}\sum_{i=1}^{n}(x_i-\bar x)^2},$$
with $\mathrm{se}=0$ for $n=1$ (`web/src/batch/stats.ts:47-72`, mirrored by `batch::combine`, `batch.rs:162-195`). The relative error is $|\mathrm{se}/x|$ (`stats.ts:76-79`). A ratio of the current run $c$ to a stored reference run $r$ carries a first-order propagated relative error
$$\frac{\delta(c/r)}{c/r}\approx\sqrt{\left(\frac{\mathrm{se}_c}{c}\right)^2+\left(\frac{\mathrm{se}_r}{r}\right)^2},$$
which treats the two runs as independent and neglects covariance (`stats.ts:82-106`); a ratio whose relative error exceeds 10 % is flagged (`stats.ts:108`). A steady-state indicator is the least-squares slope of the wall power over the window divided by its mean, per revolution (`batch.rs:142-160`); the interface warns if the mean absolute drift exceeds 5 % per revolution (`stats.ts:110-127`).

**Sensitivity.** `web/src/batch/sensitivity.ts` runs a base condition plus each of four contact coefficients (ball-ball and ball-wall friction, ball-ball and ball-wall restitution) scaled by 0.8 and 1.2, one seed per condition (`sensitivity.ts:21-26`, `:56-64`); static friction is scaled with its kinetic value and kept at least as large, and restitution is clamped to 1 (`sensitivity.ts:36-52`). Each row reports the fractional change $(x-x_\text{base})/|x_\text{base}|$ and highlights changes above 10 % (`sensitivity.ts:9`, `:67-69`, `:93-95`). This is a one-at-a-time sensitivity at a single seed; it is not a variance-based analysis, and a single seed's noise is not separated from the effect.

**Report aggregator and PDF.** For live runs, `web/src/report/aggregator.ts` accumulates per-metric statistics online with Welford's algorithm, $\delta = x-\bar x_{k-1}$, $\bar x_k=\bar x_{k-1}+\delta/k$, $M_{2,k}=M_{2,k-1}+\delta\,(x-\bar x_k)$, giving the sample standard deviation $\sqrt{M_2/(k-1)}$ in $O(1)$ memory per metric (`aggregator.ts:108-113`, `:187`). Samples outside a configurable window (warm-up, end time) are ignored (`aggregator.ts:94-95`), and the time series kept for plots is decimated by doubling a stride whenever 1000 stored points are exceeded (`aggregator.ts:26`, `:121-132`). `web/src/report/pdf.ts` builds an A4 report with jsPDF [declared dependency `jspdf ^4.2.1`] entirely in the browser (`pdf.ts:9`, `:416-418`), including the run configuration, per-metric statistics, a canvas snapshot and the model-basis statements.

**Basis statements.** Every export (PDF, live CSV and batch CSV) carries the same text from `web/src/report/assumptions.ts:14-27`: the model is a 2D cross-section, quantities marked per metre are per metre of drum length, stress energy uses the mass of the 3D sphere each disc represents, per-ball powers assume a one-ball-thick slab with an implied 3D solid fraction $\varphi_{3D}=\tfrac23\,\varphi_{2D}$ (0.547 for the default 0.82), ratios between conditions do not depend on this basis, impact statistics require $k=1$, and contact coefficients are uncalibrated. The batch CSV also contains the exact parameter JSON, the seeds and the software commit (`web/src/batch/exporters.ts:40-49`).

**Video.** `web/src/video/recorder.ts` encodes the canvas with WebCodecs through the `mediabunny` library [declared dependency `mediabunny ^1.61.0`], preferring MP4/H.264 and falling back to WebM with VP9, VP8 or AV1 (`recorder.ts:36-44`, bitrate 6 Mbit/s and 30 frames per simulated second, `:29`, `:131-132`). Frames are stamped with simulated time, $t_\text{video}=t_\text{sim}-t_{\text{sim},0}$, rather than wall-clock time, and are captured when at least $1/30$ s of simulated time has elapsed since the previous frame (`recorder.ts:1-14`, `:48-51`, `:146-147`); consequently a run that computes at 0.1x real time still yields a video whose duration equals the simulated duration, and paused intervals add no frames. The feature disables itself if `VideoEncoder` is unavailable (`recorder.ts:101`).

### 3.7 Distribution

**Toolchain.** The web app is TypeScript ^5.9 on Vite ^8.3 with Vitest ^5.0 for unit tests and Playwright ^1.63 for browser tests (`web/package.json:14-19`); runtime dependencies are `jspdf` and `mediabunny` (`package.json:20-23`). Scripts are `dev`, `build`, `test` (`vitest run`), `e2e` (`playwright test`) and `lint` (`tsc --noEmit`) (`package.json:6-13`). The production build roots asset URLs at `/MillDynamics3/` for GitHub Pages and stamps the short git commit into exported outputs (`web/vite.config.ts:14-27`, `:12`).

**Continuous integration.** `.github/workflows/deploy.yml` runs on every push to `main` with four jobs (`deploy.yml:15`, `:29`, `:47`, `:92`, `:107`):

1. `test` (`:29-43`): `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`.
2. `build` (`:47-90`, after `test`): builds the WASM core with `scripts/build-wasm.sh`, then `npm ci`, `npm run test`, `npm run lint`, `npm run build`, and uploads the bundle both as a Pages artifact and as a plain artifact.
3. `deploy` (`:92-105`, after `build`): publishes to GitHub Pages with Pages-write permission scoped to this job only.
4. `python-release` (`:107-`, after `build`): builds the Python wheel and a static zip, runs the Python unit tests, and replaces the rolling prerelease `python-latest`.

Playwright end-to-end tests are not run in this workflow; they are run locally (`npm run e2e`, `web/playwright.config.ts`, which starts the dev server, Chromium only). The fixed numbers of tests are given in Section 4.

**Python-only launcher.** A user without Node or Rust can install the pure-Python package `milldynamics` (`pyproject.toml:5-19`; `requires-python >=3.8`, no dependencies). The package is standard-library only: `server.py` is an `http.server` static file server that serves the Vite bundle under the `/MillDynamics3/` base path, with the `application/wasm` MIME type registered, and binds to `127.0.0.1` on port 8765 or a free port (`python/milldynamics/server.py:11-30`, `:50-87`); `fetch.py` downloads and unpacks `milldynamics-static.zip` from the `python-latest` release on first use, with a path-traversal check on extraction (`fetch.py:9-28`, `:31-57`); `__main__.py` is the `milldynamics` command with `--host`, `--port`, `--no-browser` (`__main__.py:6-20`). `scripts/bundle_python.py` is the maintainer-only local alternative: it runs `wasm-pack`, `npm run build`, and copies `web/dist` into `python/milldynamics/static/` (`scripts/bundle_python.py:1-24`), which is git-ignored. The simulation itself runs in the user's browser; Python only serves static files.

## 4. Verification and quality assurance

### 4.1 What is gated and what the gates are

`docs/VERIFICATION.md` was created to judge the fluid solver and the ball-fluid coupling "by numbers, not by feel", independently of the solver implementation (`docs/VERIFICATION.md:3-6`). The gated tests live in `crates/mill-core/tests/verification.rs`, with the solver adapter, analytic helpers and all tolerances in `crates/mill-core/tests/common/verification_common.rs` (`gates` module at `:14-38`). A solver swap edits only the adapter. The cases (`docs/VERIFICATION.md:24-36`):

1. **Hydrostatic** (63 mm drum, fill 0.35). The worst-5 % per-particle compression error must be at most 0.02, the mean at most 0.001, the residual kinetic energy $\mathrm{KE}/(M g\,\Delta x)\le 0.05$, and the wall-pressure profile (wall normal impulse per polar-angle bin against $\rho g\,(y_\text{surface}-y_\text{wall})\times\text{arc}$) must have relative $L_2$ error at most 0.05 (`verification_common.rs:15-18`).
2. **Impulsive spin-up.** A full disc of radius $R$ with no-slip wall started impulsively at angular speed $\omega_0$, fluid initially at rest. With $s=r/R$, $T=\nu t/R^2$ and $\lambda_n$ the positive zeros of $J_1$, the exact azimuthal velocity is
$$\frac{u_\theta}{\omega_0 R}=s+2\sum_{n}\frac{J_1(\lambda_n s)}{\lambda_n J_0(\lambda_n)}\,e^{-\lambda_n^2 T}$$
(`verification_common.rs:121-139`), and the angular-momentum fraction is
$$\frac{L(T)}{L_\infty}=1-8\sum_n\frac{e^{-\lambda_n^2 T}}{\lambda_n^2}$$
(`verification_common.rs:155-162`). The gate requires $|L-L_\text{exact}|/L_\infty\le 0.03$ at $T=0.02,\,0.05,\,0.1,\,0.2$ (`SPINUP_FRACTIONS`, `:37`) with $R^2/\nu=0.25$ s and $\Delta t=1/1920$ s, and the wall angular impulse must match $dL/dt$ within 0.03 (including the lattice's exact gravity torque).
3. **Taylor-Couette torque.** A ball of radius $R_1=R_2/2$ rotating at $\omega_1$ inside the drum of radius $R_2$, outer wall at rest. The torque per unit depth is [36]
$$T=\frac{4\pi\mu\,\omega_1R_1^2R_2^2}{R_2^2-R_1^2},$$
(`verification_common.rs:164-167`); both the ball torque and the wall torque must be within 3 % (`TC_TORQUE_REL`, `:21`).
4. **Energy closure.** With browser defaults at 60 % of critical speed and slurry on, over at least one revolution, $|E_\text{unattributed}|/W_\text{shaft}\le 0.02$, where
$$E_\text{unattributed}=\Delta E_\text{mech}-\big(W_\text{shaft}-D_\text{ball contact}-D_\text{fluid wall slip}-D_\text{fluid viscous}-E_\text{clamp removed}+E_\text{interface created}\big)$$
is the net energy change that no instrumented term explains (`crates/mill-core/src/lib.rs:56-60`, ledger fields `:61-78`).
5. **Dry sub-step convergence.** Mean wall power with 4 and 8 sub-steps within 3 % of that with 16 (two seeds); the validated maximum of `substeps` for these tests is 16, so 4/8/16 stand in for 8/16/32 (`verification.rs:445-472`).

There are further gates in the file for a rotating partially filled drum (power and torque change at most 3 % per doubling of resolution; wall power by work equals $\omega\times$ torque impulse within 1 %; fluid-only energy closure within 2 %; `verification_common.rs:24-35`). The convergence goal behind them is that mill power and torque must not change by more than 3 % when resolution, coarse-graining factor or sub-steps change (`docs/VERIFICATION.md:6-8`).

### 4.2 Test inventory

Counts below are from `grep` of the source at the time of writing, not from a test-runner report.

| Suite | Location | Count |
|---|---|---|
| Rust unit tests (in `mill-core/src`, including `mac/tests.rs` 12 and `fluid/tests.rs` 15) | `crates/mill-core/src/**` | 188 `#[test]`, 1 of them `#[ignore]` (`fluid/tests.rs:385`) |
| Rust verification suite | `crates/mill-core/tests/verification.rs` | 26 `#[test]`, of which 21 are `#[ignore]` with the measured failure in the reason string, so 5 run by default |
| Vitest unit tests | `web/tests/unit/*.test.ts` (16 files) | about 150 `it`/`test` call sites |
| Playwright end-to-end tests | `web/tests/e2e/*.spec.ts` (8 files) | 15 `test(` call sites; run locally, not in CI |
| Python launcher tests | `python/tests/test_server.py` | 7 `def test` |

The five verification tests that run by default are an analytic-helper self-consistency check, the hydrostatic kinetic-energy gate at two resolutions, the DFSPH power-equals-impulse identity, and the dry sub-step test (`verification.rs:73`, `:147-153`, `:345`, `:449`). The 21 ignored tests are listed with measured numbers in their `#[ignore = "..."]` reasons (`verification.rs:137-432`) and are run with `cargo test --release -p mill-core --test verification -- --ignored` (`docs/VERIFICATION.md:20-22`). CI runs `cargo test --workspace` (default, non-ignored tests), `npm run test` (Vitest), `npm run lint` and the Python tests; the Playwright suite and the ignored verification tests are not run in CI.

### 4.3 Full-system energy-invariant tests

The cheapest check that a position-based solver is not creating energy is to compare its total mechanical energy against the work the driving wall supplies. Four tests do so (the first three are asserted invariants; the fourth guards the diagnostic definition):

- `dem::tests::total_energy_never_increases_in_a_still_drum` (`crates/mill-core/src/dem.rs:1784`): with the drum stationary the wall supplies no energy, so ball kinetic plus gravitational energy must never increase.
- `tests::steady_cascading_charge_never_gains_more_energy_than_the_wall_supplies` (`crates/mill-core/src/lib.rs:759`): for a dry rotating drum, at every sub-step $\Delta E\le W_\text{wall}+\varepsilon$, with a fresh per-sub-step tolerance $\varepsilon$ (never a single budget accumulated over many sub-steps, `coupling.rs:1002-1008`). The test records that before the 2026-09 fix this was violated "by orders of magnitude on essentially every sub-step" (`lib.rs:774-776`).
- `coupling::tests::coupled_charge_never_gains_more_energy_than_the_wall_supplies` (`crates/mill-core/src/coupling.rs:1149`) and `coupling::tests::low_fill_cataracting_charge_does_not_gain_energy_from_the_fluid` (`coupling.rs:930`): the same invariant for the full ball-plus-fluid system (`total_mechanical_energy_j`, `lib.rs:206-208`) at the browser configuration and at a low-fill cataracting configuration. The helper `settle_and_measure_full_system_energy_invariant` (`coupling.rs:1013`) returns the worst per-sub-step margin $W_\text{wall}+\varepsilon-\Delta E$, the maximum ball speed and a centroid-distance outlier ratio.
- `dem::tests::dissipated_energy_is_never_materially_negative_while_cascading` (`dem.rs:2130`): the reported dissipated power, defined as the residual $E_\text{after predict}+W_\text{wall}-E_\text{final}$ of the same invariant, stays non-negative up to rounding.

`docs/PHYSICS.md:1163-1166` records an ablation of a broad-phase safety margin (`solve_margin`, `dem.rs:621-629`): with it set to zero the whole suite still passes except the full-system invariant, which fails with a single-sub-step violation of -22.3 J. The invariant tests therefore detect a real defect that the other tests do not.

### 4.4 Diagnostic probes

`crates/mill-core/examples/` holds standalone diagnostics, none of which is part of the test suite (each is described by its header comment):

| Probe | Purpose |
|---|---|
| `verification_probe` | Prints the verification-suite metrics at several fluid resolutions, optionally for PBF and DFSPH |
| `perf_probe` | Times configurations and prints the state hash (Section 3.3); `PERF_STAGES=1` gives a per-stage table |
| `tolerance_probe` | Physical-metric A/B gate for solver changes (Section 3.3) |
| `coarse_graining_probe` | Which output metrics are independent of the coarse-graining factor $k$; also the converged-reference power comparison |
| `stress_probe` | Convergence of the 3D-converted impact statistics against sub-steps, contact iterations and the impact-episode debounce |
| `oscillation_probe` | Reproduces and classifies periodic charge oscillation at low speed with a smooth wall |
| `infiltration_probe` | Measures whether slurry seeps into a settled ball bed |
| `escape_probe` | Containment: how far balls poke through the wall and how fast they get |
| `compression_probe`, `pbf_resolution_probe`, `still_pool_probe` | Fluid compression-error distribution; fluid-only PBF compression against resolution; a still pool's residual motion |
| `pbf_bench`, `grid_probe` | PBF side and fixed-grid side of the PBF-versus-grid benchmark |

### 4.5 Project rules

`CLAUDE.md` requires that `cargo fmt`, `cargo clippy -D warnings` and `npm run lint` pass before a task is finished; the first two are enforced in CI (`deploy.yml:41-42`) and the third in the build job (`:75`). Core quantities are SI, and runs must be reproducible from `Params` plus seed.

### 4.6 Recorded status: what passes, what fails, what is open

The following is the state recorded in `docs/VERIFICATION.md` and in the `#[ignore]` reasons; it is a negative result for the shipped slurry solver on the strict analytic criteria and should be read as such.

- **Re-measured 2026-10-19** (22 gate tests, commit f8c2d16): 6 pass, 16 fail. Changes since 2026-10-01: spin-up wall torque now passes; Taylor-Couette ball torque error fell from 99.9 % to 9-16 %; energy-closure error 0.174 to 0.092; rotating-drum 25-to-50 change 0.181 to 0.055. Everything else still fails by similar margins (docs/VERIFICATION.md, 2026-10-19). The paragraph below gives the older 2026-10-01 values.
- **Shipped position-based fluid (PBF) solver against the strict analytic gates.** The 2026-10-01 baseline states "Only the dry sub-step test and the hydrostatic KE gate pass; everything else fails" (`docs/VERIFICATION.md:44`). The recorded ignore reasons give the measured numbers: hydrostatic worst-5 % compression 0.270 at resolution 25 and 0.914 at 50 against a gate of 0.02 (mean 0.047 and 0.138 against 0.001); wall-pressure $L_2$ error 0.229 and 0.197 against 0.05; spin-up angular-momentum error 0.52/0.63/0.44/0.20 at $T=0.02/0.05/0.1/0.2$ against 0.03; spin-up wall-torque mismatch 0.468; Taylor-Couette ball torque about $10^{-3}$ to $10^{-4}$ of the analytic value and wall torque 3 % to 10 % of it; energy closure $|E_\text{unattributed}|/W_\text{shaft}=0.174$ against 0.02 (`verification.rs:137-432`). The unattributed energy grows with resolution (0.20, 0.17, 0.80, 4.35, 8.40 at resolutions 15, 25, 40, 60, 100, `docs/VERIFICATION.md:55`). These numbers date from 2026-10-01; the ball-surface no-slip boundary (which acts on the Taylor-Couette ball-torque item) and the central-force viscosity were added afterwards (`docs/VERIFICATION.md:1156-1189`, `:1259-1275`) and the ignore reasons were not necessarily re-measured after each. In addition, the harness is pinned to 1800 kg/m$^3$ and 50 Pa s, and a still pool at the shipped water viscosity is noisier (residual $\mathrm{KE}/(Mg\Delta x)=3.85$ against 0.05, `docs/VERIFICATION.md:1333-1337`).
- **Fluid motion in the shipped solver is not at rest in a still pool.** At resolution 50 the RMS speed is about 0.015 m/s and constant for 10 s; it is particle noise, not a decaying slosh. A global pressure solve (preconditioned CG) was tried and rejected because it increased the noise (`docs/VERIFICATION.md:1357-1382`).
- **Convergence of mill power with resolution (app geometry, water).** Power against a converged reference (resolution 80, 32 sub-steps, 3.026 W per metre): +19.4 %, +10.5 %, +4.6 %, +2.5 %, -4.9 %, -1.0 %, +1.7 % at resolutions 15, 20, 25, 30, 40, 50, 60; "the +-3 % gate is met only at res 50 / 60" (`docs/VERIFICATION.md:1342-1346`). The default Fast tier (resolution 25) is therefore about 5 % off the reference.
- **Coarse-graining.** A matched-mass pile test gives about -3 % at $k=2$ and +9 % at $k=3$ relative to $k=1$, with ensemble noise of about 8-9 %; at the app geometry, $k=2$ gives +20 % at resolution 25 and +12 % at 40. These are outside the 3 % gate and no correction is applied (`docs/VERIFICATION.md:1384-1402`, `:1347-1349`).
- **Peak ball overlap** at the converged reference is 13.4 % of the radius against a target below 10 % (`docs/VERIFICATION.md:1349-1350`).
- **Comparison with the resolved grid solver (10 heavy discs, $\nu=0.1$).** Converged PBF peak gravity moment is about -1.2 against the grid's -1.30 (-7 %) at $\nu=0.1$, and +10 to +14 % above the grid at $\nu=1$; PBF's dependence on viscosity is much weaker than the grid's. At resolution $\ge$ 60 and $\ge$ 32 sub-steps the PBF result is within about 3 % of its own converged value, at roughly 9 times the cost of resolution 40 with 8 sub-steps (`docs/VERIFICATION.md:1275-1313`).
- **Bingham (yield-stress) slurry.** A Papanastasiou-regularised Bingham viscosity is implemented in the shipped solver (`pbf.rs:850-853`, applied at `:2141-2143`; `BINGHAM_PLASTIC_RATIO` at `:489`). The recorded comparison with the grid solver shows the target is not met: the pile's peak does not fall with yield stress as the grid's does (`docs/VERIFICATION.md:1404-1422`). The cause is not found.
- **Passed gates of the experimental solvers.** For the fixed-grid track, the rigid-body-in-fluid gates (E4), the normal squeeze-film lubrication gate (E5) and the ball-size-versus-grid study (E7) are recorded as met; the pure-slurry free-surface and rotating-drum gates (E2, E3) are not met; the full-system pile agrees between grid resolutions to within 2 % for $\nu\ge 0.1$ but produces NaN at water-like viscosity (`docs/VERIFICATION.md:1063-1071`). The DFSPH solver passes the hydrostatic and spin-up gates to within the stated numbers when the original boundary was used but is recorded as not converging in the rotating-drum case (+38 % for 25 to 50) and is not used by the simulator (`docs/VERIFICATION.md:91-125`, `:1073-1080`).
- **Open items recorded in the document:** the dominant sub-step-dependent error term of the PBF ball-fluid coupling at $\nu\ge 1$ (torque balance not converged in the sub-step, "cause not found"), the Bingham response, and the energy-closure gap (`docs/VERIFICATION.md:1286-1292`, `:1329`).

In summary, the analytic gates were introduced to give honest numbers, and for the shipped PBF solver they show that it satisfies the rotating-drum power at the resolutions of the upper quality tiers and the energy-creation invariants of Section 4.3, while failing the strict hydrostatic, spin-up, Taylor-Couette and energy-closure gates. Section 5 states the innovation claims with this in view.

## 5. What is innovative, and what is not

Throughout this section, "to our knowledge" and "in the sources we surveyed" mean that a limited literature and web search was made (recorded in the survey that this section is based on) and that an absence of a result is not proof of absence. Each claim in 5.2 is backed by the implementation (with source location) and by the prior-work comparison of 5.1.

### 5.1 Positioning and lineage

**What is established.** Discrete element simulation of tumbling mills is long established [4, 5, 37], reviewed in [38] and [39]; in the sources we surveyed these tools run offline and the classical mill studies are dry. Coupled DEM and particle-fluid (SPH) simulation of wet mills exists, for example for SAG mills [40, 41] and for a laboratory ball mill [42]; these are offline, three-dimensional studies in which the slurry coupling is either one-way, porous-media or not verified by us to be two-way. Stirred-media-mill DEM is also established [43, 44]. The position-based family originates in computer graphics [1, 2, 9, 11, 45] and has precedents in engineering applications [46, 47, 48].

**Comparison table.** Entries are limited to what the survey supports; items the survey could not verify are marked UNVERIFIED.

| Prior work | What it does | Offline / real-time; dimension; coupling | Difference from this work |
|---|---|---|---|
| [4, 5, 37, 49] | DEM of ball-mill charge motion, power, sensitivity to speed and liner | Offline; dry; no slurry | Here: slurry included, in-browser, interactive; no breakage in either |
| [40, 41] | DEM + SPH for slurry transport in SAG mills; slurry grinding in a 3D pilot SAG mill | Offline; the 2020 study is 3D; coupling described as DEM-based porous media for 2006 (two-way status not verified by us) | Here: 2D, real-time, position-based particle fluid with an implicit joint viscosity solve; no external validation yet |
| [42] | Digital twin of a laboratory ball mill, DEM-SPH calibrated with computer vision | Offline; 3D | Here: no calibration against a physical mill |
| [43, 44] | DEM of stirred media mills (media motion, energy, collisional environment) | Offline; 3D | Different mill type; no slurry in the cited work as surveyed |
| [19, 20] | Lubrication-approximation force models inside DEM for wet milling | DEM-only force models; no particle fluid solver | Here: the same idea (Reynolds squeeze film, cf. [14, 15, 16, 17, 18]) is a sub-grid link inside a particle fluid solver, calibrated on a resolved grid solver |
| [6, 7, 8] | Coarse-grained or enlarged-particle DEM (pneumatic conveying; die filling; wet bead milling) | Offline; mostly 3D | Here: mass-preserving rule tied to fluid resolution and sub-step count in a tumbling mill |
| [46] | Position-based granular simulation (excavator-type training) | Real-time; 3D; no fluid | Establishes game-physics techniques in engineering; not claimed as new here |
| [48] | Real-time multiscale terrain model for earthmoving; agrees with a coupled DEM-multibody reference to 10-25 % | Real-time; 3D | Same; also an example of cross-checking a real-time model against DEM |
| [50] | Distributed particle simulation for mine-to-mill material tracking (digital twin) | Details not verified by us | UNVERIFIED beyond the title; not compared |
| [47] | Position-based fluid for robotic crack sealing, tuned to a laboratory experiment (1.26 mm mean fluid-level error) | Engineering use of PBF; validated to one experiment | Same lineage; not claimed as new |
| [51] | Adaptive sampling for interactive granular simulation | Interactive; graphics | Same lineage |
| Commercial DEM packages (e.g. Rocky DEM, EDEM, Millsoft) | Desktop DEM, in some cases coupled to CFD/SPH | Desktop, offline (Rocky CFD-DEM/SPH capability from a vendor snippet) | EDEM and Millsoft capabilities UNVERIFIED; omitted from comparison |
| [24, 25] | Stress intensity and stress number as grinding descriptors in stirred media mills | Analytical / experimental framework | Here: the framework's stress-energy idea is used as a relative proxy from 2D collisions (item 5) |

The computer-graphics and engineering lineage is the following: the position-based method of this simulator belongs to the PBD/XPBD/PBF family and its use outside graphics is not new.

### 5.2 Contributions we consider new

For each item the claim, the implementation, the comparison and a limit are stated.

**(1) A real-time, in-browser, deterministic, two-way-coupled ball-plus-slurry tumbling-mill model that needs no installation.**
- *Implementation.* The whole model (Section 3.1) compiles to WASM and runs in a Web Worker under the accumulator of Section 3.4 (`worker.ts:137-152`, `:256-262`); the ball solver and the slurry solver exchange momentum both ways in every sub-step (`crates/mill-core/src/coupling.rs`; the joint solve of item 3); runs are reproducible from `Params` and seed (Section 3.3); the Python launcher serves the same bundle from the standard library (`python/milldynamics/server.py`).
- *Prior work.* In the sources we surveyed, mill DEM and DEM-SPH are offline, desktop or HPC tools [40, 41, 42]; we found no browser or WASM mill or DEM simulator and no interactive real-time grinding-mill process simulator for design or training. The search was limited.
- *Limits of the claim.* "Real-time" is conditional. The recorded browser speeds predate later model changes and are stale (`docs/PERF.md:350-358`); the browser speed of the current tiers has not been measured. The native cost of the default Fast tier is about 49-67 ms per 1/60 s frame, i.e. more than one frame of wall time per frame of simulated time on the author's machine at that date (absolute figures are machine- and session-dependent, `docs/PERF.md:420`; `docs/PERF.md:352`, `:395`), so slow-motion operation is expected and is displayed by the speed badge. The model is two-dimensional (item below and 5.3). Existing offline tools are far more detailed.

**(2) Bounded-recovery, bidirectional-restitution position-based contact, verified by full-system energy-invariant tests.**
- *Mechanism.* Ball-ball and ball-wall non-penetration is a rigid (zero-compliance) projection in the position-based family [1, 2] (`dem.rs:1-9`). Two measures keep it energy-bounded. (a) Each iteration recovers at most $0.2$ of a ball diameter of overlap per contact, $c_\text{eff}=\max(c,\,-0.2\cdot 2r)$ (`MAX_RECOVERY_FRACTION`, `dem.rs:58`, `:730`, `:745`), so a deep overlap is recovered over several sub-steps instead of being converted by the velocity update into a large separation speed. (b) Restitution is applied in both directions: the relative normal velocity is driven to $v_n^\star=\max\!\big(-e\,v_n^\text{pre},\,0,\,v_n^\text{pre}\big)$, so velocity that depenetration over-manufactured is removed as well as velocity that fell short being added (`dem.rs:892-965`, target at `:961`). An earlier version only added velocity and could inject "essentially unbounded" energy (`docs/PHYSICS.md:336-346`). A swept collision test against the wall bounds tunnelling, and the broad-phase search distance is widened by the recovery budget (`dem.rs:621-629`, `:635-690`).
- *Evidence.* The invariant tests of Section 4.3 assert, at every sub-step, that the ball solver (dry and in a stationary drum) and the full ball-plus-slurry system never gain more energy than the wall supplies; before the fix the dry invariant was violated by orders of magnitude on essentially every sub-step (`lib.rs:774-776`), and an ablation of the broad-phase margin produced a -22.3 J violation that only the full-system test caught (`docs/PHYSICS.md:1163-1166`). A related coupling defect (an unbounded angular reaction coefficient, observed "in the hundreds to ~13000") was fixed by a cap and is guarded by the same test (`pbf.rs:1776-1782`, `coupling.rs:1141-1148`).
- *Prior work.* The position-based formulations are from [1, 2]. We did not find, in the sources we surveyed, a statement of an energy-conservation (or energy-boundedness) test of this kind for a position-based contact solver coupled to a particle fluid; this is a statement about our survey, not about the field.
- *Limits of the claim.* The invariant is an inequality ("no energy creation"), not a statement that dissipation is correct; restitution and friction coefficients are uncalibrated. Rolling resistance does not conserve angular momentum and a momentum-conserving rewrite was reverted because it destabilised the regression tests (`docs/PHYSICS.md:1067-1081`). Ball-ball continuous collision detection was tried and reverted. The fluid side of the full-system energy balance does not close (Section 4.6): the invariant is one-sided.

**(3) A joint implicit viscosity solve with balls as unknowns, plus calibrated sub-grid squeeze-film links, calibrated against a resolved grid solver in the same code base.**
- *Implementation.* Fluid velocities, ball translational velocities and ball spins are unknowns of one symmetric positive-definite system solved by Jacobi-preconditioned conjugate gradients, with the ball surface treated as a no-slip (Dirichlet) boundary through ghost points carrying the rigid-body velocity; linear and angular momentum are exchanged exactly between fluid and ball. The system is
$$m\,\boldsymbol v_i+\Delta t\,m\,(L\boldsymbol v)_i+\Delta t\,m\sum_k c_{ik}(\boldsymbol v_i-\boldsymbol u_k)=m\boldsymbol v_i^\ast,\qquad \boldsymbol u_k=\boldsymbol V_b+W_b\,\boldsymbol l_k^\perp,$$
with the corresponding ball momentum and angular-momentum rows (`pbf.rs:625-640`). The source notes that treating the ball explicitly was unstable at viscosities above about 10 Pa s (`pbf.rs:632-636`). The viscous operator is of central-force form, which is exactly zero for rigid rotation (`docs/VERIFICATION.md:1267-1273`). For gaps smaller than the fluid particle spacing resolves, a normal squeeze-film link between two balls, or a ball and the drum wall, adds the 2D Reynolds film force $F=c(h)\,u_n$ with
$$c=3\sqrt2\,\pi\,\nu\,R^{3/2}h^{-3/2}\,(1+c_1\varepsilon+c_2\varepsilon^2),\quad \varepsilon=\frac hR,\quad \frac1R=k_1+k_2,$$
$c_1=0.788+0.265\,Q$, $c_2=\max(0.006+0.057(Q-0.25),0)$, $Q=R^3(k_1^3+k_2^3)$, blended in with a smoothstep between gaps of 8 and 3 grid cells and with a correction for the fraction of the film that the resolved flow already delivers (`crates/mill-core/src/mac/lubrication.rs:3-15`, `:33-43`, `:86-118`; links built in `pbf.rs:543-`, gain `SQUEEZE_GAIN = 1` at `:523`, stiffness cap at `:528`; links built only for approaching pairs, `pbf.rs:536-542`). The constants were fitted to exact bipolar-coordinate Stokes solutions (disc against a drum wall; two equal discs), reproduced to 1 % for $\varepsilon\le0.4$ and 2.5 % at 0.8 (`lubrication.rs:10-12`), and the resolved-flow ratio table was measured on the grid solver (`lubrication.rs:90-93`, `docs/VERIFICATION.md:588-620`). The ball no-slip gain `BALL_BC_GAIN = 0.5` was calibrated against the grid solver (`pbf.rs:486`, `docs/VERIFICATION.md:1166-1167`).
- *Prior work.* Lubrication in DEM for wet milling exists as DEM-only force models [19, 20]; the underlying hydrodynamics is classical [14, 15] with corrections used in particle-resolved simulations [16, 17, 18]. In the sources we surveyed we did not find the same idea combined with a particle (position-based) fluid solver and calibrated on a resolved grid solver within one code base, i.e. a multi-fidelity calibration chain (exact solution to grid solver to particle solver).
- *Limits of the claim.* The squeeze film is the normal (approach) component only; tangential and rotational lubrication are open (`docs/VERIFICATION.md:650-656`). In the PBF benchmark the film "does not move the level" of the full-system result at $\nu=0.1$ (-1148 against -1152, `docs/VERIFICATION.md:1301-1302`), and it matters mostly for coarse balls. The dependence of the 2D film on curvature is calibrated, not derived for a general 2D geometry. The calibration targets the grid track, which is itself offline and does not reach water-like viscosity (NaN at $\nu=0.01$, `docs/VERIFICATION.md:1070`). Agreement with the grid at the benchmark is about 7 % to 14 % (Section 4.6), not the 3 % goal.

**(4) Mass-preserving coarse-graining with drag on the true radius and automatic raising of fluid resolution and sub-steps.**
- *Implementation.* When requested (it is off by default, `params.rs:42`), the ball diameter is enlarged by $k=\sqrt{N_\text{true}/N_\text{max}}$ so that the number of discs falls by $k^2$; because discs have unit-depth mass proportional to $r^2$, the total solid area and hence charge mass are preserved and the density is unchanged (`params.rs:571-580`, `:582-613`). The fluid drag relaxation rates use the true particle's radius, mass and inertia, not the enlarged ones (`pbf.rs:1706-1722`, `:1788-1795`; the true radius is stored per population at `dem.rs:308` and set at `:407`). Two automatic raises keep the enlarged population numerically admissible: the fluid resolution is raised to at least $\lceil R_\text{drum}/d_\text{eff}\rceil$ (capped at 200) so the fluid spacing does not exceed the effective ball diameter (`params.rs:648-656`), and the sub-step count is raised to at least $\lceil \sqrt{2gD}/(60\cdot 0.5\,d_\text{eff})\rceil$, capped at 16, so that the sub-step displacement stays near half a ball diameter (`params.rs:690-704`). Both raises are reported in the interface (`derived.ts`, `paramsPanel.ts`).
- *Evidence recorded.* For a charge raised to 1500 balls at coarse resolution, the first raise reduced the ball-fluid coupling's clamp hits from 7944 to 13 in the same window and capped the peak ball speed (`docs/PHYSICS.md:1341-1343`); the second reduced the worst sub-step displacement ratio from 0.64-1.07 to 0.37 (`docs/PHYSICS.md:1148`). Using the true radius in the drag fixed an energy-creation bug (`docs/PHYSICS.md:742-750`).
- *Prior work.* Particle-scaling rules are established [6, 7] and an enlarged-particle model was applied to wet bead milling [8]. In the sources we surveyed we did not find a mass-preserving rule combined with drag on the true radius and with automatic, resolution- and sub-step-dependent raising for a tumbling mill.
- *Limits of the claim.* Coarse-graining changes the physics. The project's own measurements show power draw off by 17-39 % at large diameter ratio, collision count scaling as $1/k^2$ and impact energy as $k^2$, and coarse-grained power and torque outside the 3 % gate (Section 4.6); per-impact statistics are therefore blocked in the interface when $k>1$ (`docs/PHYSICS.md:1255-1262`, `schema.ts:164`, `plan.ts:80-86`). Sub-step and resolution raises are heuristics, not convergence proofs, and are capped; in a pathological regime the ratio can still exceed 1 (`params.rs:690-704`).

**(5) Relative grinding proxies from the 2D run, converted to a 3D basis.**
- *Implementation.* Every contact above a speed threshold of $0.1428\sqrt{g\,d}$ (equal to 0.02 m/s for the 2 mm ball) is recorded with normal and slip speeds (`dem.rs:34-46`, `:550`, `:934-946`). `stress.rs` converts each to a stress energy using the mass of the 3D sphere the disc stands for, $m_3=\rho\cdot\tfrac43\pi r^3$, with reduced mass $m_3/2$ for ball-ball and $m_3$ for ball-wall:
$$\mathrm{SE}=\tfrac12\,m_\text{red}\,(v_n^2+v_s^2),$$
(`stress.rs:68-80`), because the 2D unit-depth disc mass would overstate SE by a size-dependent factor (about 750 for a 2 mm ball, `stress.rs:6-8`). Impacts are classed as drop ($v_n\ge v_s$) or shear and as ball-ball or ball-wall (four classes, `stress.rs:33-38`, `:68-80`), made dimensionless by $\mathrm{SE}/(m_3 g d)$ and binned log-spaced (28 bins, four per decade, `stress.rs:28-29`, `:51-58`); re-contacts closer than $0.25\sqrt{d/g}$ count as one impact episode so that the count does not grow with sub-step (`stress.rs:41-48`). Per-ball powers use a one-ball-thick slab, implying $\varphi_{3D}=\tfrac23\varphi_{2D}$ (`stress.rs:14-17`; Section 3.6).
- *Prior work.* The stress-energy and stress-number framework for stirred media mills is established [24, 25]. We use its ideas as relative descriptors from a 2D tumbling-mill run; we are not proposing a new grinding model.
- *Limits of the claim.* No breakage or grinding kinetics is simulated (Section 5.3). The proxies are for relative comparison between conditions run with the same code and settings, not absolute prediction (`stress.rs:18-21`); the contact coefficients are uncalibrated (`assumptions.ts:23`); the basis is a convention, ratios being independent of it. The proxies depend on sub-steps, contact iterations and the debounce, so `stress_probe` was built to measure the dependence (`docs/VERIFICATION.md:1485-1523`), and impact counting by contact onset was tried and rejected (`docs/VERIFICATION.md:1550`). The validity of SE as a predictor of grinding in a tumbling mill, as opposed to a stirred mill, is not shown here.

**(6) Experiment tooling built into the simulator.**
- *Implementation.* Seeded repeat runs with standard errors, whole-revolution windows and a steady-state check; one-at-a-time contact-coefficient sensitivity at $\pm 20\,\%$; reference-run ratios with propagated error; online-statistics PDF reports; CSV export carrying exact parameters, seeds and commit; simulated-time-stamped video (Sections 3.5-3.6: `batch.rs:40-195`, `stats.ts:55-127`, `sensitivity.ts:28-94`, `aggregator.ts:108-132`, `exporters.ts:40-49`, `recorder.ts:1-14`).
- *Prior work.* None of the individual techniques is new. We did not survey simulator front ends systematically; the only claim is that this combination is built into the same in-browser tool as the simulation, so that provenance travels with each output.
- *Limits of the claim.* The ratio error neglects covariance and assumes independent runs; the sensitivity is a single-seed one-at-a-time scan; seed-to-seed variance reflects only the random initial placement; and the results inherit every limit of the model.

### 5.3 What this is not

- **It is a two-dimensional cross-section model.** Quantities are per metre of drum length; there is no axial flow, no 3D packing, and the 2D void fraction is much lower and less connected than in 3D, so slurry percolates into a settled ball bed far less than in a real mill (`docs/PHYSICS.md:86-91`; the infiltration investigation recorded 41 % "covered" at the default resolution and 54 % at four times finer, and a static note was added to the interface, `schema.ts:167`).
- **No breakage, no grinding kinetics, no product size distribution.** Only collision energies, rates and power are measured (`docs/PHYSICS.md:86-91`, `stress.rs:1-8`).
- **No external validation.** The model has not been compared with a physical mill or with an independent 3D DEM code. In the sources we surveyed, validation of position-based granular flow against DEM and of PBF against analytic viscous flow was not found anywhere; the analytic gates of Section 4 fill part of the second gap but, as Section 4.6 shows, the shipped solver does not meet the strict ones. Validation against DEM or experiments would fill the first gap and is the main missing step. The results are qualitative and for relative comparison.
- **The slurry solver is not grid-converged and is approximate.** It fails the hydrostatic, spin-up, Taylor-Couette and energy-closure gates (Section 4.6); a still pool is noisy; resolution matters (power within 3 % of a converged reference only at resolution 50-60); the first-order sub-step dependence of the ball-fluid coupling is not fully understood; the full-system energy ledger does not close (`docs/VERIFICATION.md:1286-1292`).
- **Surface tension and wettability are resolution-dependent calibrations**, not first-principles properties: the cohesion and adhesion terms are set to reproduce a physical value of surface tension at a given lattice spacing (`params.rs:323-324`, `docs/PHYSICS.md` section 6.1a), and switching them off did not change the benchmarks recorded (`docs/VERIFICATION.md:1314`, `:1322-1324`). A thin wetting film is not resolved at any shipped tier (`presets.ts` header comment).
- **Viscosity regime.** The default slurry is water (0.001 Pa s), where the ball Reynolds number is large and the model uses an inertial form drag on the coupling in addition to the no-slip boundary; the resolved grid reference does not reach this regime (NaN at $\nu=0.01$), so there is no resolved reference for the shipped default. Behaviour at 50 Pa s is better characterised (the earlier default).
- **Coarse-graining distorts impact statistics and, in the project's own measurements, bulk power** (Section 5.2, item 4). It is off by default and per-impact statistics are blocked when it is on.
- **Bingham rheology is implemented in the shipped solver but not reproduced against the grid reference** (Section 4.6); it is not exposed in the UI.
- **Ball drag in water is too large (measured 2026-10-19).** The form-drag closure had no factor $\tfrac12$ (it acted as $C_D=2$); it now reads $1/\tau_\text{form}=\tfrac12\rho_0 C_D(2r_t)|\mathbf v_\text{rel}|/m_t$, $C_D=1$ (`pbf.rs`). This changed no gated physical metric (8 seeds, 27 metrics). But a ball pushed at constant speed through a quiescent water pool meets 5-17 times the cylinder drag at Re $10^3$-$3\cdot10^3$ at the shipped resolution ($d/\Delta x$ 2-4), independent of the closure strength, resolution and sub-steps; at Re $10^4$ it is 0.9-1.9 times and depends on sub-steps. A per-mechanism recording attributes 60-95 % of it to the overlap push of step 3.5 (the fluid-to-ball momentum exchange that stands in for a pressure on the ball), 38 % more to adhesion at Re $10^3$, and almost none to the form-drag closure or the viscous solve (docs/VERIFICATION.md, 2026-10-19); scaling that exchange confirmed it sets the drag, and a prototype that adds the balls to the density constraint as smooth solids (with a two-sided constraint) only halved the excess and was unstable with more iterations, so the shipped model keeps the exchange (docs/VERIFICATION.md, 2026-10-19).
- **Contact coefficients (friction, restitution, rolling resistance) are uncalibrated**, and rolling resistance does not conserve angular momentum (`docs/PHYSICS.md:1067-1081`).
- **Speed.** At the default tier the model runs slower than real time on the author's machine, and the stated speed-ups are same-session A/B figures; they are not a performance guarantee on other hardware.

---

## 6. References

Numbered in order of first citation. Entries marked with a DOI were confirmed against the Crossref record during the survey; where volume or page data could not be confirmed they are omitted rather than guessed.

1. M. Müller, B. Heidelberger, M. Hennix, J. Ratcliff. Position based dynamics. J. Vis. Commun. Image Represent. 18 (2007) 109-118. [doi:10.1016/j.jvcir.2007.01.005](https://doi.org/10.1016/j.jvcir.2007.01.005)
2. M. Macklin, M. Müller, N. Chentanez. XPBD: position-based simulation of compliant constrained dynamics. Proc. 9th Int. Conf. Motion in Games (MIG 2016) 49-54. [doi:10.1145/2994258.2994272](https://doi.org/10.1145/2994258.2994272)
3. P.A. Cundall, O.D.L. Strack. A discrete numerical model for granular assemblies. Géotechnique 29 (1979) 47-65. [doi:10.1680/geot.1979.29.1.47](https://doi.org/10.1680/geot.1979.29.1.47)
4. B.K. Mishra, R.K. Rajamani. The discrete element method for the simulation of ball mills. Appl. Math. Modelling 16 (1992) 598-604. [doi:10.1016/0307-904X(92)90035-2](https://doi.org/10.1016/0307-904X(92)90035-2)
5. P.W. Cleary. Predicting charge motion, power draw, segregation and wear in ball mills using discrete element methods. Minerals Eng. 11 (1998) 1061-1080. [doi:10.1016/S0892-6875(98)00093-4](https://doi.org/10.1016/S0892-6875(98)00093-4)
6. M. Sakai, S. Koshizuka. Large-scale discrete element modeling in pneumatic conveying. Chem. Eng. Sci. 64(3) (2009) 533-539. [doi:10.1016/j.ces.2008.10.003](https://doi.org/10.1016/j.ces.2008.10.003)
7. C. Bierwisch, T. Kraft, H. Riedel, M. Moseler. Three-dimensional discrete element models for the granular statics and dynamics of powders in cavity filling. J. Mech. Phys. Solids 57 (2009) 10-31. [doi:10.1016/j.jmps.2008.10.006](https://doi.org/10.1016/j.jmps.2008.10.006)
8. Y. Yamamoto, R. Soda, J. Kano, F. Saito. DEM simulation of bead motion during wet bead milling using an enlarged particle model. Int. J. Miner. Process. 114-117 (2012) 93-99. [doi:10.1016/j.minpro.2012.10.001](https://doi.org/10.1016/j.minpro.2012.10.001)
9. M. Macklin, M. Müller. Position based fluids. ACM Trans. Graph. 32(4) (2013) art. 104. [doi:10.1145/2461912.2461984](https://doi.org/10.1145/2461912.2461984)
10. M. Müller, D. Charypar, M. Gross. Particle-based fluid simulation for interactive applications. Proc. ACM SIGGRAPH/Eurographics Symp. Computer Animation (SCA 2003) 154-159. [doi:10.2312/SCA03/154-159](https://doi.org/10.2312/SCA03/154-159)
11. J. Bender, M. Müller, M.A. Otaduy, M. Teschner, M. Macklin. A survey on position-based simulation methods in computer graphics. Comput. Graph. Forum 33(6) (2014) 228-251. [doi:10.1111/cgf.12346](https://doi.org/10.1111/cgf.12346)
12. N. Akinci, G. Akinci, M. Teschner. Versatile surface tension and adhesion for SPH fluids. ACM Trans. Graph. 32(6) (2013) art. 182. [doi:10.1145/2508363.2508395](https://doi.org/10.1145/2508363.2508395)
13. J.P. Morris, P.J. Fox, Y. Zhu. Modeling low Reynolds number incompressible flows using SPH. J. Comput. Phys. 136 (1997) 214-226. [doi:10.1006/jcph.1997.5776](https://doi.org/10.1006/jcph.1997.5776)
14. H. Brenner. The slow motion of a sphere through a viscous fluid towards a plane surface. Chem. Eng. Sci. 16 (1961) 242-251. [doi:10.1016/0009-2509(61)80035-3](https://doi.org/10.1016/0009-2509(61)80035-3)
15. R.G. Cox, H. Brenner. The slow motion of a sphere through a viscous fluid towards a plane surface - II. Small gap widths, including inertial effects. Chem. Eng. Sci. 22 (1967) 1753-1777. [doi:10.1016/0009-2509(67)80208-2](https://doi.org/10.1016/0009-2509(67)80208-2)
16. N.-Q. Nguyen, A.J.C. Ladd. Lubrication corrections for lattice-Boltzmann simulations of particle suspensions. Phys. Rev. E 66 (2002) 046708. [doi:10.1103/PhysRevE.66.046708](https://doi.org/10.1103/PhysRevE.66.046708)
17. T. Kempe, J. Fröhlich. Collision modelling for the interface-resolved simulation of spherical particles in viscous fluids. J. Fluid Mech. 709 (2012) 445-489. [doi:10.1017/jfm.2012.343](https://doi.org/10.1017/jfm.2012.343)
18. P.I. Costa, B.J. Boersma, J. Westerweel, W.-P. Breugem. Collision model for fully resolved simulations of flows laden with finite-size particles. Phys. Rev. E 92 (2015) 053012. [doi:10.1103/PhysRevE.92.053012](https://doi.org/10.1103/PhysRevE.92.053012)
19. R. Cabiscol, T. Jansen, M. Marigo, C. Ness. Application of hydrodynamic lubrication in DEM simulations of wet bead milling chambers. Powder Technol. 384 (2021) 542-553. [doi:10.1016/j.powtec.2021.01.071](https://doi.org/10.1016/j.powtec.2021.01.071)
20. D.N. de Klerk, I. Govender, T.L. Moodley, A.N. Mainza. Using the lubrication approximation to model the effects of viscosity in DEM simulations for complex flow regimes. KONA Powder Part. J. (2025/2026). [doi:10.14356/kona.2026007](https://doi.org/10.14356/kona.2026007)
21. W.E. Lorensen, H.E. Cline. Marching cubes: a high resolution 3D surface construction algorithm. Proc. SIGGRAPH '87 (Comput. Graph. 21(4)) 163-169. [doi:10.1145/37401.37422](https://doi.org/10.1145/37401.37422)
22. G.M. Chaikin. An algorithm for high-speed curve generation. Comput. Graph. Image Process. 3 (1974) 346-349. [doi:10.1016/0146-664X(74)90028-8](https://doi.org/10.1016/0146-664X(74)90028-8)
23. P.M.C. Lacey. Developments in the theory of particle mixing. J. Appl. Chem. 4 (1954) 257-268. [doi:10.1002/jctb.5010040504](https://doi.org/10.1002/jctb.5010040504)
24. A. Kwade. Determination of the most important grinding mechanism in stirred media mills by calculating stress intensity and stress number. Powder Technol. 105 (1999) 382-388. [doi:10.1016/S0032-5910(99)00162-X](https://doi.org/10.1016/S0032-5910(99)00162-X)
25. A. Kwade, J. Schwedes. Breaking characteristics of different materials and their effect on stress intensity and stress number in stirred media mills. Powder Technol. 122 (2002) 109-121. [doi:10.1016/S0032-5910(01)00406-5](https://doi.org/10.1016/S0032-5910(01)00406-5)
26. S. Vigna. An experimental exploration of Marsaglia's xorshift generators, scrambled. ACM Trans. Math. Softw. 42(4) (2016) art. 30. [doi:10.1145/2845077](https://doi.org/10.1145/2845077)
27. J. Bender, D. Koschier. Divergence-free smoothed particle hydrodynamics. Proc. SCA 2015, 147-155. [doi:10.1145/2786784.2786796](https://doi.org/10.1145/2786784.2786796)
28. D. Koschier, J. Bender. Density maps for improved SPH boundary handling. Proc. SCA 2017, art. 1. [doi:10.1145/3099564.3099565](https://doi.org/10.1145/3099564.3099565)
29. F.H. Harlow, J.E. Welch. Numerical calculation of time-dependent viscous incompressible flow of fluid with free surface. Phys. Fluids 8 (1965) 2182-2189. [doi:10.1063/1.1761178](https://doi.org/10.1063/1.1761178)
30. C. Batty, F. Bertails, R. Bridson. A fast variational framework for accurate solid-fluid coupling. ACM Trans. Graph. 26(3) (2007) art. 100. [doi:10.1145/1239451.1239551](https://doi.org/10.1145/1239451.1239551)
31. F. Gibou, R.P. Fedkiw, L.-T. Cheng, M. Kang. A second-order-accurate symmetric discretization of the Poisson equation on irregular domains. J. Comput. Phys. 176 (2002) 205-227. [doi:10.1006/jcph.2001.6977](https://doi.org/10.1006/jcph.2001.6977)
32. G.-S. Jiang, D. Peng. Weighted ENO schemes for Hamilton-Jacobi equations. SIAM J. Sci. Comput. 21 (2000) 2126-2143. [doi:10.1137/S106482759732455X](https://doi.org/10.1137/S106482759732455X)
33. G. Russo, P. Smereka. A remark on computing distance functions. J. Comput. Phys. 163 (2000) 51-67. [doi:10.1006/jcph.2000.6553](https://doi.org/10.1006/jcph.2000.6553)
34. R. Glowinski, P. Le Tallec. Augmented Lagrangian and operator-splitting methods in nonlinear mechanics. SIAM, 1989. [doi:10.1137/1.9781611970838](https://doi.org/10.1137/1.9781611970838)
35. W.H. Herschel, R. Bulkley. Konsistenzmessungen von Gummi-Benzollösungen. Kolloid-Z. 39 (1926) 291-300. [doi:10.1007/BF01432034](https://doi.org/10.1007/BF01432034)
36. G.K. Batchelor. An Introduction to Fluid Dynamics. Cambridge University Press, 1967 (no DOI for the 1967 edition).
37. P.W. Cleary. Charge behaviour and power consumption in ball mills: sensitivity to mill operating conditions, liner geometry and charge composition. Int. J. Miner. Process. 63 (2001) 79-114. [doi:10.1016/S0301-7516(01)00037-0](https://doi.org/10.1016/S0301-7516(01)00037-0)
38. N. Weerasekara, M. Powell, P. Cleary, L. Tavares, M. Evertsson, R. Morrison, et al. The contribution of DEM to the science of comminution. Powder Technol. 248 (2013) 3-24. [doi:10.1016/j.powtec.2013.05.032](https://doi.org/10.1016/j.powtec.2013.05.032)
39. L.M. Tavares. A review of advanced ball mill modelling. KONA Powder Part. J. (2017). [doi:10.14356/kona.2017015](https://doi.org/10.14356/kona.2017015)
40. P.W. Cleary, M. Sinnott, R. Morrison. Prediction of slurry transport in SAG mills using SPH fluid flow in a dynamic DEM based porous media. Minerals Eng. 19 (2006) 1517-1527. [doi:10.1016/j.mineng.2006.08.018](https://doi.org/10.1016/j.mineng.2006.08.018)
41. P.W. Cleary, R. Morrison, M. Sinnott. Prediction of slurry grinding due to media and coarse rock interactions in a 3D pilot SAG mill using a coupled DEM + SPH model. Minerals Eng. 159 (2020) 106614. [doi:10.1016/j.mineng.2020.106614](https://doi.org/10.1016/j.mineng.2020.106614)
42. B. Doroszuk, P. Bortnowski, M. Ozdoba, R. Król. Calibrating the digital twin of a laboratory ball mill for copper ore milling: integrating computer vision and DEM-SPH simulations. Minerals 14 (2024) 407. [doi:10.3390/min14040407](https://doi.org/10.3390/min14040407)
43. M. Sinnott, P. Cleary, R. Morrison. Analysis of stirred mill performance using DEM simulation: Part 1 - Media motion, energy consumption and collisional environment. Minerals Eng. 19 (2006) 1537-1550. [doi:10.1016/j.mineng.2006.08.012](https://doi.org/10.1016/j.mineng.2006.08.012)
44. C. Jayasundara, R. Yang, A. Yu, D. Curry. Discrete particle simulation of particle flow in IsaMill - effect of grinding medium properties. Chem. Eng. J. 135 (2008) 103-112. [doi:10.1016/j.cej.2007.04.001](https://doi.org/10.1016/j.cej.2007.04.001)
45. M. Macklin, M. Müller, N. Chentanez, T.-Y. Kim. Unified particle physics for real-time applications. ACM Trans. Graph. 33(4) (2014) art. 153. [doi:10.1145/2601097.2601152](https://doi.org/10.1145/2601097.2601152)
46. D. Holz. Parallel Particles (P2): a parallel position based approach for fast and stable simulation of granular materials. VRIPHYS 2014. [doi:10.2312/vriphys.20141232](https://doi.org/10.2312/vriphys.20141232)
47. S. Schaefer, J. Xu, D. Palin, A. Al-Tabbaa. Position-based fluid simulation for robotic injection sealing of pavement cracks. J. Field Robotics (2024). [doi:10.1002/rob.22339](https://doi.org/10.1002/rob.22339)
48. M. Servin, T. Berglund, S. Nystedt. A multiscale model of terrain dynamics for real-time earthmoving simulation. Adv. Model. Simul. Eng. Sci. (2021). [doi:10.1186/s40323-021-00196-3](https://doi.org/10.1186/s40323-021-00196-3)
49. B.K. Mishra, R.K. Rajamani. Simulation of charge motion in ball mills. Part 1: experimental verifications. Int. J. Miner. Process. 40 (1994) 171-186. [doi:10.1016/0301-7516(94)90042-6](https://doi.org/10.1016/0301-7516(94)90042-6)
50. M. Servin, F. Vesterlund, E. Wallin. Digital twins with distributed particle simulation for mine-to-mill material tracking. Minerals 11 (2021) 524. [doi:10.3390/min11050524](https://doi.org/10.3390/min11050524)
51. S. Gupta, J. Keyser. Adaptive sampling for interactive simulation of granular material. Comput. Animation Virtual Worlds (2025). [doi:10.1002/cav.70062](https://doi.org/10.1002/cav.70062)
