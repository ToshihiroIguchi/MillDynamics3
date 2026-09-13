//! A tiny, deterministic pseudo-random number generator.
//!
//! We intentionally avoid pulling in the `rand` crate: the simulation must be bit-reproducible
//! across platforms/builds given the same `seed`, and a small xorshift64* generator is more than
//! enough for particle-lattice jitter and other non-critical randomness.

/// xorshift64* PRNG. Not cryptographically secure; deterministic given the same seed.
#[derive(Debug, Clone)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Creates a new generator from a seed. A seed of `0` is remapped to a fixed non-zero value,
    /// since xorshift generators cannot recover from an all-zero state.
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x9E37_79B9_7F4A_7C15 } else { seed },
        }
    }

    /// Returns the next raw `u64` in the sequence.
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Returns a `f32` uniformly distributed in `[0, 1)`.
    pub fn next_f32(&mut self) -> f32 {
        // Take the top 24 bits for a clean f32 mantissa's worth of entropy.
        ((self.next_u64() >> 40) as f32) / (1u32 << 24) as f32
    }

    /// Returns a `f32` uniformly distributed in `[lo, hi)`.
    pub fn range_f32(&mut self, lo: f32, hi: f32) -> f32 {
        lo + self.next_f32() * (hi - lo)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_for_same_seed() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn zero_seed_is_remapped() {
        let mut rng = Rng::new(0);
        // Should not panic or get stuck at zero.
        assert_ne!(rng.next_u64(), 0);
    }

    #[test]
    fn f32_in_unit_range() {
        let mut rng = Rng::new(7);
        for _ in 0..1000 {
            let v = rng.next_f32();
            assert!((0.0..1.0).contains(&v));
        }
    }
}
