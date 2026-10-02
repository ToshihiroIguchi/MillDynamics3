//! Cut-cell geometry of a circular domain on a uniform square grid.
//!
//! The grid has `n x n` cells over `[-half, half]^2`. A face aperture is the exact fraction of the
//! face segment inside the disc; weights are stored per face (`wx`: vertical faces, `wy`:
//! horizontal faces) and are dimensionless, so the weighted Laplacian is independent of `dx`.

/// Circular fluid domain centred at the origin.
#[derive(Clone, Debug)]
pub struct CircleDomain {
    pub radius: f64,
    pub n: usize,
    pub half: f64,
    /// Vertical faces: index `i + (n + 1) * j`, face `i` lies at `x = -half + i dx`, spans cell row `j`.
    pub wx: Vec<f64>,
    /// Horizontal faces: index `i + n * j` (`j` in `0..=n`), spans cell column `i`.
    pub wy: Vec<f64>,
}

impl CircleDomain {
    /// `n` cells across, domain half-width `half = radius * (1 + margin)`.
    pub fn new(radius: f64, n: usize, margin: f64) -> Self {
        let half = radius * (1.0 + margin);
        let dx = 2.0 * half / n as f64;
        let mut wx = vec![0.0; (n + 1) * n];
        let mut wy = vec![0.0; n * (n + 1)];
        for j in 0..n {
            let (y0, y1) = (-half + j as f64 * dx, -half + (j + 1) as f64 * dx);
            for i in 0..=n {
                let x = -half + i as f64 * dx;
                wx[i + (n + 1) * j] = chord_overlap(radius, x, y0, y1) / dx;
            }
        }
        for j in 0..=n {
            let y = -half + j as f64 * dx;
            for i in 0..n {
                let (x0, x1) = (-half + i as f64 * dx, -half + (i + 1) as f64 * dx);
                wy[i + n * j] = chord_overlap(radius, y, x0, x1) / dx;
            }
        }
        // Faces on the domain border carry no flux.
        for j in 0..n {
            wx[(n + 1) * j] = 0.0;
            wx[n + (n + 1) * j] = 0.0;
        }
        for i in 0..n {
            wy[i] = 0.0;
            wy[i + n * n] = 0.0;
        }
        Self {
            radius,
            n,
            half,
            wx,
            wy,
        }
    }

    pub fn dx(&self) -> f64 {
        2.0 * self.half / self.n as f64
    }

    pub fn cell_center(&self, i: usize, j: usize) -> (f64, f64) {
        let dx = self.dx();
        (
            -self.half + (i as f64 + 0.5) * dx,
            -self.half + (j as f64 + 0.5) * dx,
        )
    }

    /// Integral of `f(x, y)` over the part of cell `(i, j)` inside the disc (tensor Gauss rule on
    /// the exact chord range, 4 x-subintervals).
    pub fn integrate_cell(&self, i: usize, j: usize, f: impl Fn(f64, f64) -> f64) -> f64 {
        const GX: [f64; 4] = [
            -0.861_136_311_594_052_6,
            -0.339_981_043_584_856_3,
            0.339_981_043_584_856_3,
            0.861_136_311_594_052_6,
        ];
        const GW: [f64; 4] = [
            0.347_854_845_137_453_85,
            0.652_145_154_862_546_1,
            0.652_145_154_862_546_1,
            0.347_854_845_137_453_85,
        ];
        let dx = self.dx();
        let (x0, y0) = (-self.half + i as f64 * dx, -self.half + j as f64 * dx);
        let y1 = y0 + dx;
        let mut total = 0.0;
        let sub = 4usize;
        let hs = dx / sub as f64;
        for s in 0..sub {
            let xa = x0 + s as f64 * hs;
            for (gx, gw) in GX.iter().zip(&GW) {
                let x = xa + 0.5 * hs * (1.0 + gx);
                let h2 = self.radius * self.radius - x * x;
                if h2 <= 0.0 {
                    continue;
                }
                let h = h2.sqrt();
                let (lo, hi) = (y0.max(-h), y1.min(h));
                if hi <= lo {
                    continue;
                }
                let mut col = 0.0;
                for (gy, gwy) in GX.iter().zip(&GW) {
                    let y = lo + 0.5 * (hi - lo) * (1.0 + gy);
                    col += 0.5 * (hi - lo) * gwy * f(x, y);
                }
                total += 0.5 * hs * gw * col;
            }
        }
        total
    }
}

/// Length of `[a, b]` inside the disc along the line at coordinate `c` (the other coordinate
/// ranges over `[a, b]`).
fn chord_overlap(radius: f64, c: f64, a: f64, b: f64) -> f64 {
    let h2 = radius * radius - c * c;
    if h2 <= 0.0 {
        return 0.0;
    }
    let h = h2.sqrt();
    (b.min(h) - a.max(-h)).max(0.0)
}
