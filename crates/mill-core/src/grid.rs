//! Uniform-grid spatial hash for neighbour search.
//!
//! Shared broadphase structure used by both the DEM ball solver ([`crate::dem`], M1) and the PBF
//! fluid solver ([`crate::pbf`], M3): a dense grid over the drum's bounding box, cell size chosen
//! per-solver (max ball diameter for DEM, kernel radius `h` for PBF), rebuilt each (sub-)step.
//! See docs/PLAN.md ss3.2/3.3 for how each solver uses it.
//!
//! Implemented starting milestone M1.
