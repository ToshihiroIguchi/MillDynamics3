//! Read-only view of a fluid particle population, implemented by both the DFSPH solver
//! ([`crate::fluid::Fluid`]) and the legacy PBF one ([`crate::pbf::FluidParticles`]) so metrics
//! and the free-surface extraction do not depend on which solver produced the particles.

use glam::Vec2;

pub trait FluidView {
    fn positions(&self) -> &[Vec2];
    fn dyes(&self) -> &[f32];
    /// Kernel support radius.
    fn kernel_h(&self) -> f32;
    fn rest_density(&self) -> f32;
    fn particle_mass(&self) -> f32;
    /// Per-particle SPH density.
    fn densities(&self) -> Vec<f32>;
    fn len(&self) -> usize {
        self.positions().len()
    }
    fn is_empty(&self) -> bool {
        self.positions().is_empty()
    }
}

impl FluidView for crate::fluid::Fluid {
    fn positions(&self) -> &[Vec2] {
        &self.x
    }
    fn dyes(&self) -> &[f32] {
        &self.dye
    }
    fn kernel_h(&self) -> f32 {
        self.h
    }
    fn rest_density(&self) -> f32 {
        self.rest_density
    }
    fn particle_mass(&self) -> f32 {
        self.particle_mass
    }
    fn densities(&self) -> Vec<f32> {
        crate::fluid::Fluid::densities(self)
    }
}

impl FluidView for crate::pbf::FluidParticles {
    fn positions(&self) -> &[Vec2] {
        &self.x
    }
    fn dyes(&self) -> &[f32] {
        &self.dye
    }
    fn kernel_h(&self) -> f32 {
        self.h
    }
    fn rest_density(&self) -> f32 {
        self.rest_density
    }
    fn particle_mass(&self) -> f32 {
        self.particle_mass
    }
    fn densities(&self) -> Vec<f32> {
        crate::pbf::FluidParticles::densities(self)
    }
}
