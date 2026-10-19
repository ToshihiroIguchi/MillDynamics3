import io

def load(p):
    return open(p, encoding="utf-8").read()

def save(p, s):
    open(p, "w", encoding="utf-8", newline="\n").write(s)

PBF = r"C:\Users\toshi\python\MillDynamics3\crates\mill-core\src\pbf.rs"
PRB = r"C:\Users\toshi\python\MillDynamics3\crates\mill-core\examples\drag_probe.rs"
s = load(PBF)

def rep(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:70])
    s = s.replace(a, b)

TABLE = '''const BALL_FORM_DRAG_COEFFICIENT: f32 = 1.0;

/// PROTOTYPE switch (examples/drag_probe.rs `--bd`): 1 = balls enter the PBF density constraint as
/// smooth solids (density `rho0 * Gamma(r)` plus its gradient) and step 3.5 stops exchanging momentum.
pub static BALL_DENSITY: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);

/// TEMPORARY diagnostic: x-impulse on ball 0 per mechanism (density term in slot 5).
pub static DRAG_SPLIT: std::sync::Mutex<[f64; 6]> = std::sync::Mutex::new([0.0; 6]);
fn split_add(i: usize, v: f32) {
    if let Ok(mut g) = DRAG_SPLIT.lock() {
        g[i] += v as f64;
    }
}

/// Kernel integral of a disc: `Gamma(r)` = integral over the disc of radius `a` (centre at the
/// origin) of the 2D Poly6 kernel centred at a point at distance `r`; `dgamma` is `dGamma/dr`.
/// Tabulated on `[0, a + h]` (Koschier & Bender's density-map idea, analytic geometry for a disc).
struct GammaTable {
    a: f32,
    h: f32,
    dr: f32,
    gamma: Vec<f32>,
    dgamma: Vec<f32>,
}

impl GammaTable {
    const N: usize = 256;
    fn build(a: f32, h: f32) -> Self {
        let (a64, h64) = (a as f64, h as f64);
        let rmax = a64 + h64;
        let dr = rmax / (Self::N - 1) as f64;
        let k = 4.0 / (std::f64::consts::PI * h64.powi(8));
        let ns = 400usize;
        let gamma_at = |r: f64| -> f64 {
            let mut sum = 0.0;
            for j in 0..ns {
                let sx = (j as f64 + 0.5) / ns as f64 * h64;
                let w = k * (h64 * h64 - sx * sx).powi(3);
                let theta = if r < 1e-12 {
                    if sx < a64 {
                        2.0 * std::f64::consts::PI
                    } else {
                        0.0
                    }
                } else {
                    let q = (a64 * a64 - r * r - sx * sx) / (2.0 * r * sx);
                    if q >= 1.0 {
                        2.0 * std::f64::consts::PI
                    } else if q <= -1.0 {
                        0.0
                    } else {
                        2.0 * std::f64::consts::PI - 2.0 * q.acos()
                    }
                };
                sum += w * sx * theta;
            }
            sum * h64 / ns as f64
        };
        let g: Vec<f64> = (0..Self::N).map(|i| gamma_at(i as f64 * dr)).collect();
        let dg: Vec<f64> = (0..Self::N)
            .map(|i| {
                let (lo, hi) = (i.saturating_sub(1), (i + 1).min(Self::N - 1));
                (g[hi] - g[lo]) / ((hi - lo) as f64 * dr)
            })
            .collect();
        Self {
            a,
            h,
            dr: dr as f32,
            gamma: g.iter().map(|&x| x as f32).collect(),
            dgamma: dg.iter().map(|&x| x as f32).collect(),
        }
    }
    /// `(Gamma, dGamma/dr)` at distance `r` (zero beyond `a + h`).
    fn eval(&self, r: f32) -> (f32, f32) {
        let t = r / self.dr;
        let i = t as usize;
        if i + 1 >= Self::N {
            return (0.0, 0.0);
        }
        let f = t - i as f32;
        (
            self.gamma[i] * (1.0 - f) + self.gamma[i + 1] * f,
            self.dgamma[i] * (1.0 - f) + self.dgamma[i + 1] * f,
        )
    }
}

thread_local! {
    static GAMMA_TABLE: std::cell::RefCell<Option<GammaTable>> = const { std::cell::RefCell::new(None) };
}
'''
rep("const BALL_FORM_DRAG_COEFFICIENT: f32 = 1.0;\n", TABLE)

rep("        let base_iterations = iterations.max(1);\n        let mut iter_count = 0u32;\n        loop {\n",
'''        let base_iterations = iterations.max(1);
        let mut iter_count = 0u32;
        let ball_density = BALL_DENSITY.load(std::sync::atomic::Ordering::Relaxed) == 1
            && !balls.is_empty()
            && balls.radius > 0.0;
        let gamma_range = balls.radius + h;
        let bd_grid = if ball_density {
            GAMMA_TABLE.with(|t| {
                let mut t = t.borrow_mut();
                let stale = t
                    .as_ref()
                    .is_none_or(|g| (g.a - balls.radius).abs() > 1e-9 || (g.h - h).abs() > 1e-9);
                if stale {
                    *t = Some(GammaTable::build(balls.radius, h));
                }
            });
            Some(UniformGrid::build(&balls.x, gamma_range.max(1e-6)))
        } else {
            None
        };
        let mut bd_near: Vec<u32> = Vec::new();
        loop {
''')

rep('''                density[i] = rho;

                // Pressure-only constraint''', '''                if let Some(bg) = &bd_grid {
                    bd_near.clear();
                    bg.for_each_near(self.x[i], |b| bd_near.push(b));
                    GAMMA_TABLE.with(|t| {
                        let t = t.borrow();
                        let t = t.as_ref().unwrap();
                        for &b in &bd_near {
                            let d = self.x[i] - balls.x[b as usize];
                            let r = d.length();
                            if r >= gamma_range || r < 1e-9 {
                                continue;
                            }
                            let (g, dg) = t.eval(r);
                            rho += rest_density * g;
                            grad_self += (dg / r) * d;
                        }
                    });
                }
                density[i] = rho;

                // Pressure-only constraint''')

rep('''                delta_p[i] = sum * (mass / rest_density);
            }
''', '''                delta_p[i] = sum * (mass / rest_density);
                if let Some(bg) = &bd_grid {
                    bd_near.clear();
                    bg.for_each_near(self.x[i], |b| bd_near.push(b));
                    GAMMA_TABLE.with(|t| {
                        let t = t.borrow();
                        let t = t.as_ref().unwrap();
                        for &b in &bd_near {
                            let d = self.x[i] - balls.x[b as usize];
                            let r = d.length();
                            if r >= gamma_range || r < 1e-9 {
                                continue;
                            }
                            let (_, dg) = t.eval(r);
                            // Boundary term of this particle's own constraint: displacement along
                            // the Gamma gradient; its momentum reaction goes to the ball.
                            let dx = lambda[i] * (dg / r) * d;
                            delta_p[i] += dx;
                            coupling.impulses[b as usize] -= mass * dx / dt;
                            if b == 0 {
                                split_add(5, -(mass * dx.x / dt));
                            }
                        }
                    });
                }
            }
''')

rep("                    let impulse_mag = full_impulse_mag.min(reduced_mass * arrest_speed);\n",
'''                    let impulse_mag = if ball_density {
                        0.0
                    } else {
                        full_impulse_mag.min(reduced_mass * arrest_speed)
                    };
''')
save(PBF, s)

s = load(PRB)
s = s.replace("use mill_core::pbf::FluidParticles;", "use mill_core::pbf::{FluidParticles, BALL_DENSITY};\nuse std::sync::atomic::Ordering;", 1)
a = '    let (wet, st) = (opt("--wet"), opt("--st"));\n'
assert a in s
s = s.replace(a, a + '    BALL_DENSITY.store(one("--bd", 0.0) as u32, Ordering::Relaxed);\n', 1)
save(PRB, s)
print("patched")
