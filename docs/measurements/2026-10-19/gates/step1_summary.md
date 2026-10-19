# Step 1: strict gates re-measured 2026-10-19 on commit f8c2d16 (shipped PBF, unchanged code)
Command: cargo test --release -p mill-core --test verification -- --include-ignored --skip dfsph_rotating_drum --skip dfsph_energy_closure --nocapture --test-threads=3  (78.5 s test time; 6 passed, 16 failed)
Harness inputs as shipped in tests/common/verification_common.rs: hydrostatic + rotating drum pinned 1800 kg/m3 / 50 Pa s; spin-up and Taylor-Couette use Params::default() (density 1000, viscosity set from tau); energy closure = Params::default() (water 0.001 Pa s), res 25.

| gate (PBF) | gate | 2026-10-01 | 2026-10-19 | pass |
|---|---|---|---|---|
| hydrostatic worst-5% compression res25 / res50 | 0.02 | 0.270 / 0.914 | 0.266 / 0.919 | no |
| hydrostatic mean compression res25 / res50 | 0.001 | 0.047 / 0.138 | 0.0472 / 0.138 | no |
| hydrostatic KE/(M g dx) res25 / res50 | 0.05 | pass | pass | yes |
| hydrostatic wall-pressure L2 res25 / res50 | 0.05 | 0.229 / 0.197 | 0.224 / 0.200 | no |
| spin-up |L-Lex|/Linf res25 at t/tau 0.02/0.05/0.1/0.2 | 0.03 | 0.52/0.63/0.44/0.20 | 0.514/0.607/0.419/0.174 | no |
| spin-up |L-Lex|/Linf res50 | 0.03 | 0.52/0.60/0.47/0.23 | 0.520/0.587/0.440/0.183 | no |
| spin-up wall torque vs dL/dt res25 / res50 | 0.03 | 0.468 / 0.451 | PASS / PASS (<0.03) | **yes (changed)** |
| Taylor-Couette ball torque rel err res25 / res50 | 0.03 | 0.999 / 1.000 | 0.0928 / 0.158 (|T|=1.4965e-2 / 1.3894e-2 vs 1.6497e-2 N m) | no (big improvement) |
| Taylor-Couette wall torque rel err res25 / res50 | 0.03 | 0.966 / 0.895 | 0.571 / 0.929 (-7.08e-3 / 1.16e-3 vs 1.6497e-2) | no |
| energy closure |unattributed|/shaft work (full Simulation res25) | 0.02 | 0.174 | 0.0920 (-0.371 J of 4.037 J over 1.183 s; water default now) | no |
| rotating drum 25->50 power/torque change (pinned 1800/50) | 0.03 | 0.181 | 0.0553 (2.156 -> 2.275 W; 0.2036 -> 0.2149 N m) | no (improved) |
| dry sub-step convergence | 0.03 | pass | pass | yes |
DFSPH (not shipped; numbers consistent with their ignore strings): hydro L2 0.104/0.058; spin-up 0.056,0.035 (res25), 0.038 (res50). dfsph_rotating_drum_* and dfsph_energy_closure not re-run (5-20 min).
Notes: TC ball torque quoted as |T_sim| (sign negative = reaction). Counts: 6 pass / 16 fail (of 22 run; includes 4 DFSPH).
