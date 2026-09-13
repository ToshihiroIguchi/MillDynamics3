//! Two-way ball<->fluid momentum exchange.
//!
//! Per fluid sub-step: (1) fluid particles overlapping a ball are projected to its surface during
//! the PBF density-constraint iterations, with the resulting position corrections accumulated as
//! an impulse on that ball; (2) no-slip blending at the ball surface exchanges viscous momentum;
//! (3) the accumulated fluid force/torque is applied (clamped) to the ball for the DEM sub-steps
//! that follow, after which the balls' updated positions/velocities become the moving boundary for
//! the next fluid sub-step. See docs/PLAN.md ss3.4.
//!
//! Implemented starting milestone M4.
