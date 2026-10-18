//! Opt-in per-stage wall-clock timing for native profiling (cargo feature `stage-timing`).
//!
//! With the feature off (the default) this module contains only the `stage_mark!` macro, which
//! expands to nothing, so the solver pays no cost and its arithmetic is untouched. With the
//! feature on, `stage_mark!(Name)` attributes the time elapsed since the previous mark (or since
//! `stage_reset!`) to stage `Name`. Used by `examples/perf_probe.rs` (`PERF_STAGES=1`).

#[cfg(feature = "stage-timing")]
mod imp {
    use std::cell::Cell;
    use std::time::Instant;

    /// Solver stages timed by [`mark`].
    #[derive(Clone, Copy, Debug)]
    pub enum Stage {
        PredictGrid,
        Density,
        PushCohesion,
        WallProj,
        DragBuoy,
        DensityRecompute,
        ViscWeights,
        ViscBoundary,
        ViscLinks,
        ViscosityCg,
        ShearRate,
        FluidOther,
        Dem,
    }

    /// Number of [`Stage`] variants.
    pub const STAGE_COUNT: usize = 13;

    /// Display names, indexed by `Stage as usize`.
    pub const STAGE_NAMES: [&str; STAGE_COUNT] = [
        "PBF predict + grid build",
        "density iterations",
        "push/cohesion/adhesion (3.5-3.6)",
        "wall projection (4-6)",
        "ball drag + buoyancy (6.5-6.6)",
        "density recompute",
        "viscosity: Morris weights (+ Bingham rates)",
        "viscosity: ball boundary pairs",
        "viscosity: squeeze links",
        "viscosity: CG solve (+ ball bookkeeping)",
        "mean_shear_rate",
        "fluid clamps/other (7.5-8)",
        "DEM step",
    ];

    thread_local! {
        static LAST: Cell<Option<Instant>> = const { Cell::new(None) };
        static TOTALS_NS: Cell<[u64; STAGE_COUNT]> = const { Cell::new([0; STAGE_COUNT]) };
    }

    /// Starts a new interval at "now" without attributing the elapsed time to any stage.
    pub fn reset_clock() {
        LAST.with(|l| l.set(Some(Instant::now())));
    }

    /// Attributes the time since the previous mark/reset to `stage`.
    pub fn mark(stage: Stage) {
        let now = Instant::now();
        LAST.with(|l| {
            if let Some(prev) = l.get() {
                TOTALS_NS.with(|t| {
                    let mut a = t.get();
                    a[stage as usize] += now.duration_since(prev).as_nanos() as u64;
                    t.set(a);
                });
            }
            l.set(Some(now));
        });
    }

    /// Zeroes the accumulated per-stage totals.
    pub fn clear_totals() {
        TOTALS_NS.with(|t| t.set([0; STAGE_COUNT]));
    }

    /// Accumulated per-stage totals in nanoseconds.
    pub fn totals_ns() -> [u64; STAGE_COUNT] {
        TOTALS_NS.with(|t| t.get())
    }
}

#[cfg(feature = "stage-timing")]
pub use imp::*;

/// Attributes the time since the previous mark to the named [`Stage`]; no-op without the feature.
#[cfg(feature = "stage-timing")]
macro_rules! stage_mark {
    ($s:ident) => {
        $crate::stage_timing::mark($crate::stage_timing::Stage::$s)
    };
}
#[cfg(not(feature = "stage-timing"))]
macro_rules! stage_mark {
    ($s:ident) => {};
}

/// Restarts the stage clock; no-op without the feature.
#[cfg(feature = "stage-timing")]
macro_rules! stage_reset {
    () => {
        $crate::stage_timing::reset_clock()
    };
}
#[cfg(not(feature = "stage-timing"))]
macro_rules! stage_reset {
    () => {};
}
