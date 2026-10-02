//! Fixed-grid incompressible-flow infrastructure (experiment track E0..E9, see the plan in
//! `docs/VERIFICATION.md`). Independent of the particle solvers in `fluid/` and `pbf/`.
//!
//! Cell-centred pressure, face-based velocity, cut-cell geometry through face apertures, and a
//! multigrid-preconditioned CG solver for the weighted (variational) Poisson problem.

pub mod embedded;
pub mod multigrid;
pub mod verify;

pub use embedded::CircleDomain;
pub use multigrid::{Level, PoissonSolver, SolveStats};

#[cfg(test)]
mod tests;
