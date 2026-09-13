//! Derived simulation metrics for the HUD and validation tests.
//!
//! Toe/shoulder angles and charge centroid (from the ball population near the wall), slurry pool
//! angular extent and free-surface line, a Lacey mixing index from the dye tracer field, estimated
//! power draw from wall torque, and debug energy/overlap/density-error checks. See docs/PLAN.md
//! ss3.5.
//!
//! Implemented starting milestone M3 (slurry-only metrics) and M1 (charge metrics may land
//! earlier alongside the DEM solver if useful for its own validation tests).
