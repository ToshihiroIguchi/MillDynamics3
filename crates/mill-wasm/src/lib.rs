//! Thin `wasm-bindgen` wrapper around `mill_core`.
//!
//! All WASM-specific glue (JS-facing naming, JSON marshalling, panic hook) lives here so
//! `mill-core` itself stays a plain, natively-testable/benchmarkable Rust crate (see
//! docs/PLAN.md ss1 for the crate split rationale).
//!
//! Through M0 the exposed surface is just `new`/`step`/angle-and-time getters, matching
//! [`mill_core::Simulation`] at this milestone. Zero-copy typed-array views ("ptrs") for ball and
//! fluid particle data are added in M1/M3 once there is per-particle data to expose.

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

    /// Returns the current parameters, JSON-encoded (e.g. so the UI can read back derived values
    /// via `mill_core::Params::effective_media`).
    #[wasm_bindgen(js_name = paramsJson)]
    pub fn params_json(&self) -> Result<String, JsValue> {
        serde_json::to_string(self.inner.params()).map_err(|e| JsValue::from_str(&e.to_string()))
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
