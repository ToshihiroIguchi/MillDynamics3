//! Free-surface extraction for rendering.
//!
//! Splats fluid particles onto a scalar density-like field on a grid spanning the drum's bounding
//! box, then runs marching squares at a threshold to produce closed polylines approximating the
//! slurry free surface, lightly smoothed (Chaikin) for a less blocky outline. See docs/PLAN.md
//! ss3.5.
//!
//! Implemented starting milestone M3.
