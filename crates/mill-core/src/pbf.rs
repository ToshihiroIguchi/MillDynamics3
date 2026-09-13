//! Position Based Fluids (PBF) solver for the slurry.
//!
//! 2D Poly6/Spiky SIPH-style kernels, density-constraint projection with artificial-pressure
//! anti-clustering, drum-wall boundary projection with wall-velocity blending, and XSPH-based
//! viscosity. Dye tracer field for mixing visualization/metrics. See docs/PLAN.md ss3.3 for the
//! full per-substep algorithm and the rationale for choosing PBF over WCSPH/LBM/MPM.
//!
//! Implemented starting milestone M3.
