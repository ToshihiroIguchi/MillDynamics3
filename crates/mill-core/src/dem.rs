//! Soft-sphere DEM for grinding media (balls).
//!
//! Spring-dashpot normal contact + Coulomb-clamped tangential spring + rolling resistance, against
//! both other balls and the drum wall/lifters ([`crate::geometry::Drum`]). Ball population is
//! seeded from [`crate::params::Params::effective_media`] (post coarse-graining, if any) rather
//! than the raw UI media diameter. See docs/PLAN.md ss3.2 for the full contact model and
//! integration scheme.
//!
//! Implemented starting milestone M1.
