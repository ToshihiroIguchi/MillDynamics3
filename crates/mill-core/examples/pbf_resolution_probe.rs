//! Fluid-only (no balls) PBF compression error vs. lattice resolution and iteration floor, in a
//! still 63 mm drum at the default slurry. Isolates the density-constraint solver from the
//! ball<->fluid coupling.
//! `cargo run --release -p mill-core --example pbf_resolution_probe -- [settle_s]`

use mill_core::geometry::Drum;
use mill_core::params::{LiftersParams, Params};
use mill_core::pbf::FluidParticles;

fn main() {
    let settle_s: f32 = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(1.0);
    let params = Params::default();
    let radius_m = params.mill.radius_m();
    let slurry = params.slurry;
    let drum = Drum::new(
        radius_m,
        0.0,
        LiftersParams {
            count: 0,
            ..LiftersParams::default()
        },
    );
    let dt = 1.0 / 480.0;
    println!("res  iters  n_particles  mean_c   max_c   p99_c");
    let variant = std::env::var("PROBE_VARIANT").unwrap_or_default();
    let mut slurry = slurry;
    match variant.as_str() {
        "no_st" => slurry.surface_tension_n_m = 0.0,
        "no_visc" => slurry.viscosity_pa_s = 0.0,
        "no_st_no_visc" => {
            slurry.surface_tension_n_m = 0.0;
            slurry.viscosity_pa_s = 0.0;
        }
        _ => {}
    }
    println!("variant={variant:?}");
    let res_list: Vec<u32> = std::env::var("PROBE_RES")
        .map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect())
        .unwrap_or_else(|_| vec![15, 25, 40, 60, 100]);
    for res in res_list {
        for iters in [3u32, 20] {
            let mut fluid = FluidParticles::seed_lattice(&slurry, radius_m, res, &[], 0.0);
            for _ in 0..(settle_s / dt) as usize {
                fluid.step(&drum, 0.0, &slurry, iters, dt);
            }
            let rho = fluid.densities();
            // Where is the over-dense population? Bin the worst 5 % by distance to the wall.
            let mut idx: Vec<usize> = (0..rho.len()).collect();
            idx.sort_by(|&a, &b| rho[b].partial_cmp(&rho[a]).unwrap());
            let top = &idx[..(idx.len() / 20).max(1)];
            let mean_wall_gap = top
                .iter()
                .map(|&i| (radius_m - fluid.x[i].length()) / (radius_m / res as f32))
                .sum::<f32>()
                / top.len() as f32;
            let ymin = fluid.x.iter().map(|p| p.y).fold(f32::MAX, f32::min);
            let mean_depth = top
                .iter()
                .map(|&i| (fluid.x[i].y - ymin) / (radius_m / res as f32))
                .sum::<f32>()
                / top.len() as f32;
            println!(
                "   worst5%: mean wall gap = {mean_wall_gap:.2} dx, mean height above bottom = {mean_depth:.1} dx"
            );
            let mut c: Vec<f32> = rho
                .iter()
                .map(|r| (r / fluid.rest_density - 1.0).max(0.0))
                .collect();
            c.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mean = c.iter().sum::<f32>() / c.len() as f32;
            println!(
                "{res:>3}  {iters:>5}  {:>11}  {mean:.4}  {:.4}  {:.4}",
                c.len(),
                c[c.len() - 1],
                c[(c.len() as f32 * 0.99) as usize]
            );
        }
    }
}
