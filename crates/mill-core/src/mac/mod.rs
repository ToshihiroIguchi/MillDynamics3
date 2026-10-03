//! Fixed-grid incompressible-flow infrastructure (experiment track E0..E9, see the plan in
//! `docs/VERIFICATION.md`). Independent of the particle solvers in `fluid/` and `pbf/`.
//!
//! Cell-centred pressure, face-based velocity, cut-cell geometry through face apertures, and a
//! multigrid-preconditioned CG solver for the weighted (variational) Poisson problem.

pub mod embedded;
pub mod flow;
pub mod multigrid;
pub mod verify;
pub mod viscous;

pub use embedded::CircleDomain;
pub use flow::{Flow, Mesh};
pub use multigrid::{Level, PoissonSolver, SolveStats};
pub use viscous::{rigid_wall_torque_correction, FluidGrid, Helmholtz, Link, WallLoad};

#[cfg(test)]
mod tests;
