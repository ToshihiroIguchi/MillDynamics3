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
        PushOverlap,
        Cohesion,
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
    pub const STAGE_COUNT: usize = 15;

    /// Display names, indexed by `Stage as usize`.
    pub const STAGE_NAMES: [&str; STAGE_COUNT] = [
        "PBF predict + grid build",
        "density iterations",
        "ball overlap push (3.5)",
        "fluid cohesion (3.6a)",
        "ball adhesion (3.6b)",
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

    /// Solver event counters (diagnostics only; see [`counters`]).
    #[derive(Clone, Copy, Default, Debug)]
    pub struct Counters {
        /// Sub-steps run.
        pub substeps: u64,
        /// Histogram of the final density-loop iteration count (index = count, clamped to 12).
        pub density_hist: [u64; 13],
        /// Sum of step-7 viscosity CG iterations.
        pub cg_sum: u64,
        /// Histogram of CG iterations: 0, 1, 2-5, 6-10, 11-25, 26-50, 51-100, 101+.
        pub cg_hist: [u64; 8],
        /// Sub-steps where the with-balls viscosity solve ran.
        pub ball_solves: u64,
        /// Boundary pairs built in those sub-steps.
        pub pairs: u64,
        /// Pair c relative to the max c of the same fluid particle: <1%, 1-5%, 5-20%, >=20%.
        pub pair_c_hist: [u64; 4],
        /// Squeeze links built.
        pub links_built: u64,
        /// Squeeze links surviving the active-set filter.
        pub links_active: u64,
    }

    thread_local! {
        static COUNTERS: Cell<Counters> = Cell::new(Counters::default());
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
        COUNTERS.with(|c| c.set(Counters::default()));
        TOTALS_NS.with(|t| t.set([0; STAGE_COUNT]));
    }

    fn upd(f: impl FnOnce(&mut Counters)) {
        COUNTERS.with(|c| {
            let mut v = c.get();
            f(&mut v);
            c.set(v);
        });
    }

    /// Records a sub-step's final density-loop iteration count.
    pub fn density_iters(n: u32) {
        upd(|c| {
            c.substeps += 1;
            c.density_hist[(n as usize).min(12)] += 1;
        });
    }

    /// Records a sub-step's viscosity CG iteration count (0 when no solve ran).
    pub fn cg_iters(n: u32) {
        upd(|c| {
            c.cg_sum += n as u64;
            let b = match n {
                0 => 0,
                1 => 1,
                2..=5 => 2,
                6..=10 => 3,
                11..=25 => 4,
                26..=50 => 5,
                51..=100 => 6,
                _ => 7,
            };
            c.cg_hist[b] += 1;
        });
    }

    /// Records the boundary pairs `(fluid index, c)` of one with-balls solve; `n` = fluid count.
    pub fn ball_pairs(pairs: impl Iterator<Item = (u32, f32)> + Clone, n: usize) {
        let mut max = vec![0.0f32; n];
        let mut count = 0u64;
        for (f, c) in pairs.clone() {
            count += 1;
            let m = &mut max[f as usize];
            *m = m.max(c);
        }
        upd(|k| {
            k.ball_solves += 1;
            k.pairs += count;
            for (f, c) in pairs {
                let r = c / max[f as usize];
                let b = if r < 0.01 {
                    0
                } else if r < 0.05 {
                    1
                } else if r < 0.20 {
                    2
                } else {
                    3
                };
                k.pair_c_hist[b] += 1;
            }
        });
    }

    /// Records squeeze links built.
    pub fn links_built(n: usize) {
        upd(|c| c.links_built += n as u64);
    }

    /// Records squeeze links surviving the active-set filter.
    pub fn links_active(n: usize) {
        upd(|c| c.links_active += n as u64);
    }

    /// Accumulated event counters.
    pub fn counters() -> Counters {
        COUNTERS.with(|c| c.get())
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

/// Records a solver event via `stage_timing::$f(args)`; no-op (args unevaluated) without the feature.
#[cfg(feature = "stage-timing")]
macro_rules! count {
    ($f:ident, $($a:expr),*) => {
        $crate::stage_timing::$f($($a),*)
    };
}
#[cfg(not(feature = "stage-timing"))]
macro_rules! count {
    ($f:ident, $($a:expr),*) => {};
}
