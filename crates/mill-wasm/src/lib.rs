//! Thin `wasm-bindgen` wrapper around `mill_core`.
//!
//! All WASM-specific glue (JS-facing naming, JSON marshalling, panic hook) lives here so
//! `mill-core` itself stays a plain, natively-testable/benchmarkable Rust crate (see
//! docs/PLAN.md ss1 for the crate split rationale).
//!
//! Exposes drum kinematics, the ball (media) population, the fluid (slurry) population, and its
//! free-surface contour. `Vec<f32>` returns marshal to a fresh `Float32Array` per call (copied,
//! not zero-copy); true zero-copy typed-array views are a possible M6 performance follow-up.

use mill_core::{Params, Simulation as CoreSimulation};
use wasm_bindgen::prelude::*;

/// `wasm-bindgen` entry point wrapping [`mill_core::Simulation`].
#[wasm_bindgen]
pub struct Simulation {
    inner: CoreSimulation,
}

#[wasm_bindgen]
impl Simulation {
    /// Creates a new simulation from a JSON-encoded [`Params`] (or `undefined`/omitted for
    /// mill-core's defaults).
    #[wasm_bindgen(constructor)]
    pub fn new(params_json: Option<String>) -> Result<Simulation, JsValue> {
        let params: Params = match params_json {
            Some(json) => {
                serde_json::from_str(&json).map_err(|e| JsValue::from_str(&e.to_string()))?
            }
            None => Params::default(),
        };
        let inner = CoreSimulation::new(params).map_err(|e| JsValue::from_str(&e))?;
        Ok(Simulation { inner })
    }

    /// Replaces the parameters (JSON-encoded), without resetting drum angle / sim time. Callers
    /// decide whether a given change instead needs a full reset (a new `Simulation`), see
    /// docs/PLAN.md ss4.1.
    #[wasm_bindgen(js_name = setParams)]
    pub fn set_params(&mut self, params_json: &str) -> Result<(), JsValue> {
        let params: Params =
            serde_json::from_str(params_json).map_err(|e| JsValue::from_str(&e.to_string()))?;
        self.inner
            .set_params(params)
            .map_err(|e| JsValue::from_str(&e))
    }

    /// Advances the simulation by `dt` seconds of simulation time.
    pub fn step(&mut self, dt: f32) {
        self.inner.step(dt);
    }

    /// The fixed sub-step size (s) this simulation's `simulation.substeps` implies at the
    /// project's nominal 60 Hz target frame rate. See `mill_core::Simulation::fixed_sub_dt`.
    #[wasm_bindgen(js_name = fixedSubDt)]
    pub fn fixed_sub_dt(&self) -> f32 {
        self.inner.fixed_sub_dt()
    }

    /// Advances the simulation by exactly one fixed sub-step (`fixedSubDt()`). See
    /// `mill_core::Simulation::step_fixed` for the accumulator-loop calling convention this is
    /// meant for (docs/PLAN.md ss4.1) and why `resetFrameStats` is a separate call.
    #[wasm_bindgen(js_name = stepFixed)]
    pub fn step_fixed(&mut self) {
        self.inner.step_fixed();
    }

    /// Resets the per-frame diagnostics driving a `stepFixed` accumulator loop should call once
    /// per rendered frame, before that frame's `stepFixed` calls. See
    /// `mill_core::Simulation::reset_frame_stats`.
    #[wasm_bindgen(js_name = resetFrameStats)]
    pub fn reset_frame_stats(&mut self) {
        self.inner.reset_frame_stats();
    }

    #[wasm_bindgen(js_name = drumAngle)]
    pub fn drum_angle(&self) -> f32 {
        self.inner.drum_angle()
    }

    #[wasm_bindgen(js_name = simTime)]
    pub fn sim_time(&self) -> f64 {
        self.inner.sim_time()
    }

    /// Number of ball (grinding media) particles currently simulated (post coarse-graining, if
    /// any -- see `mill_core::Params::effective_media`).
    #[wasm_bindgen(js_name = ballCount)]
    pub fn ball_count(&self) -> u32 {
        self.inner.balls().len() as u32
    }

    /// The uniform ball radius (m). All balls currently share one effective radius (a size
    /// distribution is future work).
    #[wasm_bindgen(js_name = ballRadiusM)]
    pub fn ball_radius_m(&self) -> f32 {
        self.inner.balls().radius
    }

    /// Ball center positions, flattened as `[x0, y0, x1, y1, ...]` (m). Copied into a fresh
    /// `Float32Array` each call; zero-copy typed-array views are a possible M6 performance
    /// follow-up once profiling calls for it.
    #[wasm_bindgen(js_name = ballPositions)]
    pub fn ball_positions(&self) -> Vec<f32> {
        let balls = self.inner.balls();
        let mut out = Vec::with_capacity(balls.len() * 2);
        for p in &balls.x {
            out.push(p.x);
            out.push(p.y);
        }
        out
    }

    /// Ball orientations (radians), one per ball, in the same order as [`Simulation::ball_positions`].
    /// Useful for rendering a spin indicator confirming rolling behaviour.
    #[wasm_bindgen(js_name = ballOrientations)]
    pub fn ball_orientations(&self) -> Vec<f32> {
        self.inner.balls().theta.clone()
    }

    /// Number of fluid (slurry) particles currently simulated. Zero if `slurry.enabled` was
    /// `false` when the simulation was (re)created.
    #[wasm_bindgen(js_name = fluidCount)]
    pub fn fluid_count(&self) -> u32 {
        self.inner.fluid().len() as u32
    }

    /// Fluid particle positions, flattened as `[x0, y0, x1, y1, ...]` (m).
    #[wasm_bindgen(js_name = fluidPositions)]
    pub fn fluid_positions(&self) -> Vec<f32> {
        let fluid = self.inner.fluid();
        let mut out = Vec::with_capacity(fluid.len() * 2);
        for p in &fluid.x {
            out.push(p.x);
            out.push(p.y);
        }
        out
    }

    /// Fluid dye tracer values (`[0, 1]`), one per particle, same order as
    /// [`Simulation::fluid_positions`]. Used to colour-code mixing.
    #[wasm_bindgen(js_name = fluidDye)]
    pub fn fluid_dye(&self) -> Vec<f32> {
        self.inner.fluid().dye.clone()
    }

    /// Free-surface contour(s) at the current fluid state, flattened as
    /// `[n_polys, len_0, x, y, ..., len_1, ...]` (docs/PLAN.md ss3.5). Recomputed on demand.
    #[wasm_bindgen(js_name = fluidSurface)]
    pub fn fluid_surface(&self) -> Vec<f32> {
        mill_core::surface::flatten_polygons(&self.inner.fluid_surface())
    }

    /// Returns the current parameters, JSON-encoded (e.g. so the UI can read back derived values
    /// via `mill_core::Params::effective_media`).
    #[wasm_bindgen(js_name = paramsJson)]
    pub fn params_json(&self) -> Result<String, JsValue> {
        serde_json::to_string(self.inner.params()).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    /// Derived metrics (toe/shoulder, slurry pool extent, mixing index, debug checks, docs/
    /// PLAN.md ss3.5) at the current state, JSON-encoded. Recomputed on demand.
    #[wasm_bindgen(js_name = metricsJson)]
    pub fn metrics_json(&self) -> Result<String, JsValue> {
        serde_json::to_string(&self.inner.metrics()).map_err(|e| JsValue::from_str(&e.to_string()))
    }
}

/// Returns mill-core's default parameters, JSON-encoded, so the UI can seed the parameters modal
/// (e.g. "Reset to defaults") without duplicating default values in TypeScript.
#[wasm_bindgen(js_name = defaultParamsJson)]
pub fn default_params_json() -> Result<String, JsValue> {
    serde_json::to_string(&Params::default()).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Installs a panic hook that forwards Rust panics to the browser console with a real stack
/// trace, instead of an opaque "unreachable executed" WASM trap. Call once at worker startup.
#[wasm_bindgen(js_name = setPanicHook)]
pub fn set_panic_hook() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}
