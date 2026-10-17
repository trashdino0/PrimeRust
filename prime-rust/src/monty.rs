//! Hand-written Montgomery arithmetic for the Miller-Rabin hot path.
//!
//! Profile (§10): one 6300-bit `modpow` costs 130.2 ms — 99.9% of a
//! candidate. num-bigint's portable Montgomery (base 2^32, generic carry
//! code LLVM won't turn into `adc`/`mulx`) is the ceiling; this module owns
//! the loop instead, limb by limb, with x86 carry intrinsics on the fast
//! path and a portable fallback elsewhere.
//!
//! Notation: `R = 2^(64n)`, `n0 = -m^-1 mod 2^64`. `mont(x) = x·R mod m`.
//! All moduli here are odd (prime-search candidates); even inputs fall back
//! to num-bigint (asserted, never silently wrong).
//!
//! Roadmap is in `progress.md`. Phase 0 (this file so far): exact,
//! oracle-tested scaffolding — [`n0inv`] and limb conversions. Kernels land
//! in Phase 1 behind `#[target_feature]` dispatch.

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::One;

/// Bits per limb. Fixed at 64: matches `_addcarry_u64` / `mulx` and keeps
/// one limb per register on x86-64.
#[cfg(test)]
pub(crate) const LIMB_BITS: u32 = 64;

/// Number of u64 limbs covering `bits` bits.
#[cfg(test)]
pub(crate) fn nlimbs(bits: u32) -> usize {
    bits.div_ceil(LIMB_BITS) as usize
}

/// Montgomery `n0 = -m0^-1 mod 2^64` for the least-significant limb `m0` of
/// an odd modulus (Newton iteration, 6 doublings: 1→64 bits).
///
/// Pure wrapping arithmetic — no division, no branches. `m0` must be odd
/// (debug-asserted); the result satisfies `m0·n0 ≡ -1 (mod 2^64)` exactly.
pub(crate) fn n0inv(m0: u64) -> u64 {
    debug_assert!(m0 & 1 == 1, "Montgomery needs an odd modulus");
    let mut inv = 1u64;
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(m0.wrapping_mul(inv)));
    }
    inv.wrapping_neg()
}

/// BigUint → little-endian u64 limbs, left-padded with zeros to `n` limbs.
///
/// Byte-based on purpose: O(n) memcpy against the O(n²) kernel, and robust
/// to num-bigint's internal digit width. Panics if the value needs more
/// than `n` limbs (caller bug — sizes derive from the candidate's bits).
pub(crate) fn to_limbs_padded(x: &BigUint, n: usize) -> Vec<u64> {
    let bytes = x.to_bytes_be();
    assert!(
        bytes.len() <= n * 8,
        "value needs more than {n} limbs ({} bytes)",
        bytes.len()
    );
    let mut padded = vec![0u8; n * 8];
    padded[n * 8 - bytes.len()..].copy_from_slice(&bytes);
    let (chunks, _) = padded.as_chunks::<8>();
    chunks
        .iter()
        .rev()
        .map(|c| u64::from_be_bytes(*c))
        .collect()
}

/// Little-endian u64 limbs → BigUint (inverse of [`to_limbs_padded`]).
pub(crate) fn from_limbs(limbs: &[u64]) -> BigUint {
    let mut bytes = Vec::with_capacity(limbs.len() * 8);
    for &l in limbs.iter().rev() {
        bytes.extend_from_slice(&l.to_be_bytes());
    }
    BigUint::from_bytes_be(&bytes)
}

/// Fixed-width Montgomery context for one odd modulus (< 2^6399).
///
/// Holds the modulus limbs, `n0`, `R² mod m`, and `mont(1)` so repeated
/// exponentiations against the same modulus (Miller-Rabin runs dozens per
/// candidate) pay setup once. All values crossing the boundary are
/// `< m`; conversions are O(n) byte copies against the O(n²) kernel.
pub(crate) struct MontyCtx100 {
    m: [u64; 100],
    n0: u64,
    r2: [u64; 100],
    one_m: [u64; 100],
}

impl MontyCtx100 {
    /// Build for `m`. Returns `None` for even moduli or `m >= 2^6399`
    /// (no spare top limb) — callers fall back to num-bigint.
    pub(crate) fn new(m: &BigUint) -> Option<Self> {
        if m.is_even() || m.bits() > 6399 {
            return None;
        }
        let ma = limbs_100(m);
        let n0 = n0inv(ma[0]);
        let r2big = (BigUint::one() << 12800u64) % m;
        let r2 = limbs_100(&r2big);
        let one_m = monty_cios_100(&limbs_100(&BigUint::one()), &r2, &ma, n0);
        Some(Self {
            m: ma,
            n0,
            r2,
            one_m,
        })
    }

    /// `x` (BigUint, `< m`) → Montgomery domain (`x·R mod m`).
    pub(crate) fn enter(&self, x: &BigUint) -> [u64; 100] {
        let xa = limbs_100(x);
        debug_assert!(lt_limbs(&xa, &self.m), "enter requires x < m");
        monty_cios_100(&xa, &self.r2, &self.m, self.n0)
    }

    /// Montgomery multiply / square (stays in-domain).
    pub(crate) fn mul(&self, a: &[u64; 100], b: &[u64; 100]) -> [u64; 100] {
        monty_cios_100(a, b, &self.m, self.n0)
    }

    pub(crate) fn square(&self, a: &[u64; 100]) -> [u64; 100] {
        monty_square_100(a, &self.m, self.n0)
    }

    /// `mont(1)` for Montgomery-domain comparisons (MR tail).
    pub(crate) fn mont_one(&self) -> &[u64; 100] {
        &self.one_m
    }

    /// Montgomery domain → standard value (`X·R^-1 mod m`). Test/oracle
    /// use; the MR path never leaves the domain.
    #[allow(dead_code)]
    pub(crate) fn exit(&self, a: &[u64; 100]) -> BigUint {
        let mut one = [0u64; 100];
        one[0] = 1;
        from_limbs(&monty_cios_100(a, &one, &self.m, self.n0))
    }

    /// In-domain sliding-window (`w = 5`) exponentiation: `base_m^exp`
    /// (both in-domain, exponent as 100 limbs). ~6300 squares + ~1050 muls
    /// per 6300-bit round; table of 16 odd powers (12.8 KB) stays L1-hot.
    pub(crate) fn pow_mont(&self, base_m: &[u64; 100], exp: &[u64; 100]) -> [u64; 100] {
        let base2 = self.square(base_m);
        let mut table = [*base_m; 32];
        table[0] = self.one_m;
        for v in (3..32).step_by(2) {
            table[v] = self.mul(&table[v - 2], &base2);
        }
        // Top set bit of the exponent (d is odd and nonzero here).
        let mut top = 0usize;
        for i in (0..100).rev() {
            if exp[i] != 0 {
                top = i * 64 + (64 - exp[i].leading_zeros() as usize) - 1;
                break;
            }
        }
        let bit = |k: usize| (exp[k / 64] >> (k % 64)) & 1 == 1;
        let mut acc = self.one_m;
        let mut i = top as i64;
        while i >= 0 {
            if !bit(i as usize) {
                acc = self.square(&acc);
                i -= 1;
            } else {
                // Longest window (≤ 5 bits) ending in a 1 bit.
                let mut l = 5.min(i + 1);
                while !bit((i - l + 1) as usize) {
                    l -= 1;
                }
                let mut v = 0usize;
                for k in (i - l + 1..=i).rev() {
                    acc = self.square(&acc);
                    v = (v << 1) | usize::from(bit(k as usize));
                }
                acc = self.mul(&acc, &table[v]);
                i -= l;
            }
        }
        acc
    }
}

/// Strict `<` on 100-limb little-endian values (for `enter` contracts).
fn lt_limbs(x: &[u64; 100], m: &[u64; 100]) -> bool {
    for i in (0..100).rev() {
        if x[i] != m[i] {
            return x[i] < m[i];
        }
    }
    false
}

/// BigUint (`< 2^6400`) → fixed 100 limbs. Panics on overflow (caller
/// bug — candidates are `< 2^6399` by construction).
pub(crate) fn limbs_100(x: &BigUint) -> [u64; 100] {
    to_limbs_padded(x, 100)
        .as_slice()
        .try_into()
        .expect("100 limbs")
}

/// Fixed-width Montgomery squaring: `out = x·x·R^-1 mod m` at 100 limbs.
///
/// Same contract as [`monty_cios_100`] restricted to squaring, at ~half the
/// limb-muls: accumulate the upper triangle `S = Σ_{i≤j} x[i]x[j]·2^64(i+j)`
/// (n²/2 muls), then answer `2S − D` with `D = Σ_i x[i]²·2^128i` (the
/// diagonal, subtracted back after doubling). Reduction + finale identical
/// to SOS. This matters because ~2/3 of a square-multiply exponentiation
/// is squarings — and num-bigint runs those through the same full
/// multiplier as products, so a dedicated square beats it structurally.
///
/// Exactness: all accumulators are `u128` (never narrowed except by
/// low-word-take with the high half explicitly threaded); diagonal
/// subtract is an exact borrow chain. Cross-checked against
/// [`monty_cios_100`]`(x, x)` as well as the division oracle.
pub(crate) fn monty_square_100(x: &[u64; 100], m: &[u64; 100], n0: u64) -> [u64; 100] {
    debug_assert!(m[0] & 1 == 1, "Montgomery needs an odd modulus");
    debug_assert!(
        m[99] < (1u64 << 63),
        "spare top limb required: m must be < 2^(64n-1)"
    );
    // Upper triangle S (fits 200 limbs: S < x² < 2^12800). Carry `c` stays
    // u128 throughout — exact, never narrowed.
    let mut s = [0u64; 200];
    for i in 0..100 {
        let mut c: u128 = 0;
        for j in i..100 {
            let uv = s[i + j] as u128 + x[i] as u128 * x[j] as u128 + c;
            s[i + j] = uv as u64;
            c = uv >> 64;
        }
        let mut k = i + 100;
        while c != 0 {
            debug_assert!(k < 200, "square carry overflow");
            let nv = s[k] as u128 + c;
            s[k] = nv as u64;
            c = nv >> 64;
            k += 1;
        }
    }
    // Answer x² = 2S − D, D = Σ_i x[i]²·2^128i. Double S keeping bit 12800
    // in `hic`, then subtract D with the borrow chain extended into `hic`;
    // hic == 0 at the end is exactly x² < 2^12800.
    let mut hic = 0u64;
    for v in s.iter_mut() {
        let nv = (*v as u128) * 2 + hic as u128;
        *v = nv as u64;
        hic = (nv >> 64) as u64;
    }
    let mut borrow = 0u64;
    for i in 0..100 {
        let sq = x[i] as u128 * x[i] as u128;
        let (d1, o1) = s[2 * i].overflowing_sub(sq as u64);
        let (d2, o2) = d1.overflowing_sub(borrow);
        s[2 * i] = d2;
        let b1 = u64::from(o1 | o2);
        let (e1, p1) = s[2 * i + 1].overflowing_sub((sq >> 64) as u64);
        let (e2, p2) = e1.overflowing_sub(b1);
        s[2 * i + 1] = e2;
        borrow = u64::from(p1 | p2);
    }
    hic = hic.wrapping_sub(borrow);
    debug_assert!(hic == 0, "square exceeds 200 limbs?");
    // Montgomery reduction of the 200-limb square. Same loop as SOS but on
    // fixed arrays; t gets the SOS guard limb (201) — copy the shape that
    // already passed every oracle.
    let mut t = [0u64; 201];
    t[..200].copy_from_slice(&s);
    for i in 0..100 {
        let q = t[i].wrapping_mul(n0);
        let mut carry = 0u64;
        for j in 0..100 {
            let uv = t[i + j] as u128 + q as u128 * m[j] as u128 + carry as u128;
            t[i + j] = uv as u64;
            carry = (uv >> 64) as u64;
        }
        let mut k = i + 100;
        while carry != 0 {
            debug_assert!(k < 201, "reduction carry overflow");
            let (nv, ov) = t[k].overflowing_add(carry);
            t[k] = nv;
            carry = u64::from(ov);
            k += 1;
        }
    }
    // Result in t[100..200] (< 1.5m): conditional subtract, mask select.
    debug_assert!(t[200] == 0, "result exceeds 100 limbs (m too large?)");
    let mut sub = [0u64; 100];
    let mut borrow = 0u64;
    for i in 0..100 {
        let (d1, o1) = t[100 + i].overflowing_sub(m[i]);
        let (d2, o2) = d1.overflowing_sub(borrow);
        sub[i] = d2;
        borrow = u64::from(o1 | o2);
    }
    let mask = borrow.wrapping_neg();
    let mut out = [0u64; 100];
    for (i, o) in out.iter_mut().enumerate() {
        *o = (t[100 + i] & mask) | (sub[i] & !mask);
    }
    out
}

/// Fixed-width CIOS for the production size: 6300-bit moduli in 100 limbs
/// (99 value + spare top limb). Identical algorithm to [`monty_mul_cios`];
/// concrete array sizes let LLVM unroll and keep `t` on the stack (no
/// allocs). Correctness rides the same differential tests (dedicated
/// 100-limb cases below); any divergence from the slice version is a bug.
pub(crate) fn monty_cios_100(
    a: &[u64; 100],
    b: &[u64; 100],
    m: &[u64; 100],
    n0: u64,
) -> [u64; 100] {
    debug_assert!(m[0] & 1 == 1, "Montgomery needs an odd modulus");
    debug_assert!(
        m[99] < (1u64 << 63),
        "spare top limb required: m must be < 2^(64n-1)"
    );
    let mut t = [0u64; 101];
    for &bi in b.iter() {
        let mut carry = 0u64;
        for j in 0..100 {
            let uv = t[j] as u128 + a[j] as u128 * bi as u128 + carry as u128;
            t[j] = uv as u64;
            carry = (uv >> 64) as u64;
        }
        let spill_sum = t[100] as u128 + carry as u128;
        t[100] = spill_sum as u64;
        let spill = (spill_sum >> 64) as u64;
        let q = t[0].wrapping_mul(n0);
        debug_assert!(
            t[0].wrapping_add(q.wrapping_mul(m[0])) == 0,
            "low word must cancel"
        );
        let mut c = (t[0] as u128 + q as u128 * m[0] as u128) >> 64;
        for j in 1..100 {
            let uv = t[j] as u128 + q as u128 * m[j] as u128 + c;
            t[j - 1] = uv as u64;
            c = uv >> 64;
        }
        let uv = t[100] as u128 + (spill as u128) * (1u128 << 64) + c;
        t[99] = uv as u64;
        t[100] = (uv >> 64) as u64;
    }
    let mut sub = [0u64; 101];
    let mut borrow = 0u64;
    for i in 0..=100 {
        let mi = if i < 100 { m[i] } else { 0 };
        let (d1, o1) = t[i].overflowing_sub(mi);
        let (d2, o2) = d1.overflowing_sub(borrow);
        sub[i] = d2;
        borrow = u64::from(o1 | o2);
    }
    let mask = borrow.wrapping_neg();
    let mut out = [0u64; 100];
    for (i, o) in out.iter_mut().enumerate() {
        *o = (t[i] & mask) | (sub[i] & !mask);
    }
    // Only the SELECTED top matters: in the keep-t branch sub[100] is
    // wrapping garbage (0 - 0 - borrow), which is fine — it is discarded.
    debug_assert!(
        (t[100] & mask) | (sub[100] & !mask) == 0,
        "unreduced output"
    );
    out
}

/// CIOS Montgomery multiplication: same contract as [`monty_mul`].
///
/// Coarsely Integrated Operand Scanning: per outer step, one multiply row
/// (`t += a·b[i]`) immediately followed by one reduction row
/// (`t += q·m`, low word discarded). Single pass over a rolling `n+1`
/// window instead of SOS's full `2n` product plus second pass.
///
/// Exactness (no truncation anywhere — audit): every `u64×u64` sum lives in
/// `u128`; the multiply-row spill `t[n] + carry < 2^65` rides in a `u128`
/// into the reduction row (`spill` is 0/1, threaded exactly, never
/// asserted away). The finale borrow chain covers all `n+1` window limbs;
/// only `out[n] == 0` (true value `< m`) is asserted — and every
/// differential case re-verifies reducedness through the BigUint oracle,
/// including all-ones moduli (maximal carry propagation).
#[cfg(test)]
pub(crate) fn monty_mul_cios(a: &[u64], b: &[u64], m: &[u64], n0: u64) -> Vec<u64> {
    let n = m.len();
    debug_assert!(a.len() == n && b.len() == n, "limb count mismatch");
    debug_assert!(m[0] & 1 == 1, "Montgomery needs an odd modulus");
    debug_assert!(
        n == 0 || m[n - 1] < (1u64 << 63),
        "spare top limb required: m must be < 2^(64n-1)"
    );
    if n == 0 {
        return Vec::new();
    }
    let mut t = vec![0u64; n + 1];
    for &bi in b {
        // Multiply row: t[0..=n] += a[0..n]·bi; spill bit exact via u128.
        let mut carry = 0u64;
        for j in 0..n {
            let uv = t[j] as u128 + a[j] as u128 * bi as u128 + carry as u128;
            t[j] = uv as u64;
            carry = (uv >> 64) as u64;
        }
        let spill_sum = t[n] as u128 + carry as u128;
        t[n] = spill_sum as u64;
        let spill = (spill_sum >> 64) as u64; // 0 or 1, exact
                                              // Reduction row: q zeroes t[0]; add q·m, shift window down,
                                              // consuming the spill bit exactly.
        let q = t[0].wrapping_mul(n0);
        debug_assert!(
            t[0].wrapping_add(q.wrapping_mul(m[0])) == 0,
            "low word must cancel"
        );
        let mut c = (t[0] as u128 + q as u128 * m[0] as u128) >> 64;
        for j in 1..n {
            let uv = t[j] as u128 + q as u128 * m[j] as u128 + c;
            t[j - 1] = uv as u64;
            c = uv >> 64;
        }
        let uv = t[n] as u128 + (spill as u128) * (1u128 << 64) + c;
        t[n - 1] = uv as u64;
        t[n] = (uv >> 64) as u64;
    }
    // Conditional subtract over all n+1 window limbs, branchless select.
    let mut sub = vec![0u64; n + 1];
    let mut borrow = 0u64;
    for i in 0..=n {
        let mi = if i < n { m[i] } else { 0 };
        let (d1, o1) = t[i].overflowing_sub(mi);
        let (d2, o2) = d1.overflowing_sub(borrow);
        sub[i] = d2;
        borrow = u64::from(o1 | o2);
    }
    let mask = borrow.wrapping_neg(); // borrow=1 → window was < m, keep t
    let mut out: Vec<u64> = sub
        .iter()
        .enumerate()
        .map(|(i, &s)| (t[i] & mask) | (s & !mask))
        .collect();
    debug_assert!(out[n] == 0, "unreduced output (value >= m?)");
    out.truncate(n);
    out
}

/// SOS Montgomery multiplication: `out = a·b·R^-1 mod m`, `R = 2^(64n)`.
///
/// All inputs are `n` little-endian limbs with `a, b < m` and
/// `m < 2^(64n-1)` — i.e. the top limb of `m` must be zero (callers pass
/// one spare limb: `n = nlimbs(bits) + 1`). Proof sketch: with `m < R/2`
/// and `a,b < m`, the reduced value is `< m²/R + m < 1.5m`, so one
/// conditional subtract yields `< m` in `n` limbs; the `2n+1`-limb
/// accumulator provably holds every intermediate (product `< m²`,
/// reduction add `< m·R`). Without the spare limb the `< 2m` output of a
/// top-bit-set modulus would not fit `n` limbs — same reason mbedTLS-style
/// code keeps guard limbs.
///
/// `n0 = n0inv(m[0])`. Returns `n` limbs `< m`.
///
/// Separated Operand Scanning: full schoolbook product (u128 accumulators —
// one `mul` per inner step on x86-64) then Montgomery reduction. O(n²)
// limb-muls like num-bigint's Montgomery but in native u64 limbs instead
// of its base-2^32 form (4x fewer limb ops). Safe Rust throughout —
// intrinsics only if this disappoints (Phase 1 microbench decides).
#[cfg(test)]
pub(crate) fn monty_mul(a: &[u64], b: &[u64], m: &[u64], n0: u64) -> Vec<u64> {
    let n = m.len();
    debug_assert!(a.len() == n && b.len() == n, "limb count mismatch");
    debug_assert!(m[0] & 1 == 1, "Montgomery needs an odd modulus");
    debug_assert!(
        n == 0 || m[n - 1] < (1u64 << 63),
        "spare top limb required: m must be < 2^(64n-1)"
    );
    if n == 0 {
        return Vec::new();
    }
    // Product t = a·b (< m², fits 2n limbs; array has one guard limb).
    let mut t = vec![0u64; 2 * n + 1];
    for (i, &bi) in b.iter().enumerate() {
        let mut carry = 0u64;
        for (j, &aj) in a.iter().enumerate() {
            let uv = t[i + j] as u128 + aj as u128 * bi as u128 + carry as u128;
            t[i + j] = uv as u64;
            carry = (uv >> 64) as u64;
        }
        t[i + n] = carry;
    }
    // Reduction: fold out the low n limbs, one Montgomery step each.
    for i in 0..n {
        let q = t[i].wrapping_mul(n0);
        let mut carry = 0u64;
        for j in 0..n {
            let uv = t[i + j] as u128 + q as u128 * m[j] as u128 + carry as u128;
            t[i + j] = uv as u64;
            carry = (uv >> 64) as u64;
        }
        // Carry chain into the guard region (value bounds guarantee it
        // terminates inside the array; the assert pins that).
        let mut k = i + n;
        while carry != 0 {
            debug_assert!(k < 2 * n + 1, "reduction carry overflow");
            let (nv, ov) = t[k].overflowing_add(carry);
            t[k] = nv;
            carry = u64::from(ov);
            k += 1;
        }
    }
    // Result in t[n..2n+1] is < 1.5m < 2^64n: top guard limb is zero.
    debug_assert!(t[2 * n] == 0, "result exceeds n limbs (m too large?)");
    // Conditional subtract, branchless mask select.
    let mut sub = vec![0u64; n];
    let mut borrow = 0u64;
    for i in 0..n {
        let (d1, o1) = t[n + i].overflowing_sub(m[i]);
        let (d2, o2) = d1.overflowing_sub(borrow);
        sub[i] = d2;
        borrow = u64::from(o1 | o2);
    }
    let mask = borrow.wrapping_neg(); // borrow=1 → t was < m, keep t
    sub.iter()
        .enumerate()
        .map(|(i, &s)| (t[n + i] & mask) | (s & !mask))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::RandBigInt;
    use num_traits::One;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    #[test]
    fn n0inv_satisfies_defining_identity() {
        // m0·n0 ≡ -1 (mod 2^64) for odd limbs: the whole contract.
        let mut rng = StdRng::seed_from_u64(0x5EED_0001);
        for _ in 0..200 {
            let mut m0 = rng.gen::<u64>();
            m0 |= 1;
            assert_eq!(m0.wrapping_mul(n0inv(m0)), u64::MAX, "m0={m0}");
        }
        // Edges: 1, all-ones, large prime limb.
        for m0 in [1u64, u64::MAX, 0xFFFFFFFF_FFFFFFFFu64 - 12, 3u64] {
            assert_eq!(m0.wrapping_mul(n0inv(m0)), u64::MAX, "m0={m0}");
        }
    }

    #[test]
    fn limb_conversions_roundtrip() {
        let mut rng = StdRng::seed_from_u64(0x5EED_0002);
        for bits in [2u32, 3, 63, 64, 65, 511, 512, 513, 6299, 6300] {
            for _ in 0..5 {
                let x = rng.gen_biguint(bits.into());
                let n = nlimbs(bits.max(1));
                let limbs = to_limbs_padded(&x, n);
                assert_eq!(limbs.len(), n);
                assert_eq!(from_limbs(&limbs), x, "bits={bits}");
            }
        }
        // Exact layout: 0x0102030405060708 → limbs [0x0102.., ...] little-endian.
        let x = BigUint::from(0x0102_0304_0506_0708u64);
        assert_eq!(to_limbs_padded(&x, 2), [0x0102_0304_0506_0708, 0]);
        assert_eq!(from_limbs(&[0x0102_0304_0506_0708, 0]), x);
        assert_eq!(from_limbs(&[]), BigUint::from(0u32));
    }

    #[test]
    fn monty_mul_matches_division_oracle() {
        // Spec check without any inverse: kernel output G (n limbs, < m)
        // must satisfy G·R ≡ a·b (mod m), R = 2^(64n) — using only BigUint
        // mul+mod as the oracle.
        let mut rng = StdRng::seed_from_u64(0x5EED_0003);
        let zero = BigUint::from(0u32);
        for limbs in [1usize, 2, 3, 7, 8, 9, 16, 33, 64, 98, 99, 100] {
            let bits = limbs as u32 * 64;
            // One spare zero top limb: the kernel requires m < 2^(64n-1).
            let n = limbs + 1;
            for _ in 0..15 {
                let mut m = rng.gen_biguint(bits.into());
                m.set_bit(bits as u64 - 1, true);
                m.set_bit(0, true);
                let mv = to_limbs_padded(&m, n);
                let n0 = n0inv(mv[0]);
                let a = rng.gen_biguint_range(&zero, &m);
                let b = rng.gen_biguint_range(&zero, &m);
                let av = to_limbs_padded(&a, n);
                let bv = to_limbs_padded(&b, n);
                for got in [
                    monty_mul(&av, &bv, &mv, n0),
                    monty_mul_cios(&av, &bv, &mv, n0),
                ] {
                    assert_eq!(got.len(), n);
                    let g = from_limbs(&got);
                    assert!(g < m, "output must be reduced");
                    let r = BigUint::one() << (n as u64 * 64);
                    assert_eq!((g * r) % &m, (&a * &b) % &m, "limbs={limbs}");
                }
            }
        }
        // Edges: zeros, m-1, tiny and all-ones moduli (max carry stress).
        let all_ones = BigUint::from(u64::MAX);
        // Wide all-ones modulus (100 limbs: maximal carry propagation at
        // production width) with max inputs.
        let wide_ones = (BigUint::one() << (100u64 * 64)) - BigUint::from(1u32);
        let edges: Vec<(BigUint, BigUint, BigUint)> = vec![
            (
                "0".parse().unwrap(),
                "0".parse().unwrap(),
                "3".parse().unwrap(),
            ),
            (
                "12".parse().unwrap(),
                "10".parse().unwrap(),
                "13".parse().unwrap(),
            ),
            (
                all_ones.clone() - BigUint::from(1u32),
                all_ones.clone() - BigUint::from(1u32),
                all_ones,
            ),
            (
                wide_ones.clone() - BigUint::from(1u32),
                wide_ones.clone() - BigUint::from(1u32),
                wide_ones,
            ),
        ];
        for (a, b, m) in &edges {
            let n = nlimbs(m.bits() as u32) + 1;
            let mv = to_limbs_padded(m, n);
            let av = to_limbs_padded(a, n);
            let bv = to_limbs_padded(b, n);
            let n0 = n0inv(mv[0]);
            for got in [
                monty_mul(&av, &bv, &mv, n0),
                monty_mul_cios(&av, &bv, &mv, n0),
            ] {
                let g = from_limbs(&got);
                assert!(g < *m, "output must be reduced");
                let r = BigUint::one() << (n as u64 * 64);
                assert_eq!((g * r) % m, (a * b) % m, "a={a} b={b} m={m}");
            }
        }
    }

    #[test]
    fn monty_cios_100_matches_oracle() {
        // Production width (6300-bit m in 100 limbs): 30 random cases plus
        // the all-ones crucible, all against the division oracle.
        let mut rng = StdRng::seed_from_u64(0x5EED_0100);
        let zero = BigUint::from(0u32);
        for _ in 0..30 {
            let mut m = rng.gen_biguint(6300u64);
            m.set_bit(6299, true);
            m.set_bit(0, true);
            let mv = to_limbs_padded(&m, 100);
            let ma: &[u64; 100] = mv.as_slice().try_into().expect("100 limbs");
            let n0 = n0inv(ma[0]);
            let a = rng.gen_biguint_range(&zero, &m);
            let b = rng.gen_biguint_range(&zero, &m);
            let avv = to_limbs_padded(&a, 100);
            let aa: &[u64; 100] = avv.as_slice().try_into().expect("100");
            let bvv = to_limbs_padded(&b, 100);
            let bb: &[u64; 100] = bvv.as_slice().try_into().expect("100");
            let got = monty_cios_100(aa, bb, ma, n0);
            let g = from_limbs(&got);
            assert!(g < m, "output must be reduced");
            let r = BigUint::one() << (100u64 * 64);
            assert_eq!((g * r) % &m, (&a * &b) % &m);
        }
        let ones = (BigUint::one() << (100u64 * 64)) - BigUint::from(1u32);
        // Largest odd below 2^6399: maximal carry propagation that still
        // honors the spare top limb (m < 2^(64n-1)).
        let m = ones >> 1;
        let mvv = to_limbs_padded(&m, 100);
        let ma: &[u64; 100] = mvv.as_slice().try_into().expect("100");
        let a = m.clone() - BigUint::from(2u32);
        let avv = to_limbs_padded(&a, 100);
        let aa: &[u64; 100] = avv.as_slice().try_into().expect("100");
        let got = monty_cios_100(aa, aa, ma, n0inv(ma[0]));
        let g = from_limbs(&got);
        assert!(g < m, "output must be reduced");
        let r = BigUint::one() << (100u64 * 64);
        assert_eq!((g * r) % &m, (&a * &a) % &m);
    }

    #[test]
    fn monty_square_100_matches_oracle_and_mul() {
        // Square kernel vs division oracle AND vs monty_cios_100(x, x):
        // two independent implementations must agree with each other and
        // with G·R ≡ x² (mod m).
        let mut rng = StdRng::seed_from_u64(0x5EED_5000);
        let zero = BigUint::from(0u32);
        for _ in 0..30 {
            let mut m = rng.gen_biguint(6300u64);
            m.set_bit(6299, true);
            m.set_bit(0, true);
            let mvv = to_limbs_padded(&m, 100);
            let ma: &[u64; 100] = mvv.as_slice().try_into().expect("100");
            let n0 = n0inv(ma[0]);
            let x = rng.gen_biguint_range(&zero, &m);
            let xvv = to_limbs_padded(&x, 100);
            let xx: &[u64; 100] = xvv.as_slice().try_into().expect("100");
            let sq = monty_square_100(xx, ma, n0);
            let via_mul = monty_cios_100(xx, xx, ma, n0);
            assert_eq!(sq, via_mul, "square/mul cross-check");
            let g = from_limbs(&sq);
            assert!(g < m, "output must be reduced");
            let r = BigUint::one() << (100u64 * 64);
            assert_eq!((g * r) % &m, (&x * &x) % &m);
        }
        // Maximal-carry edge: x = m-2 under the 2^6399-1 modulus.
        let m = (BigUint::one() << 6399u32) - BigUint::from(1u32);
        let mvv = to_limbs_padded(&m, 100);
        let ma: &[u64; 100] = mvv.as_slice().try_into().expect("100");
        let x = m.clone() - BigUint::from(2u32);
        let xvv = to_limbs_padded(&x, 100);
        let xx: &[u64; 100] = xvv.as_slice().try_into().expect("100");
        let sq = monty_square_100(xx, ma, n0inv(ma[0]));
        assert_eq!(sq, monty_cios_100(xx, xx, ma, n0inv(ma[0])));
        let g = from_limbs(&sq);
        assert!(g < m, "output must be reduced");
        let r = BigUint::one() << (100u64 * 64);
        assert_eq!((g * r) % &m, (&x * &x) % &m);
    }

    #[test]
    fn monty_modpow_matches_biguint() {
        // Small moduli (debug-fast): random + edges vs BigUint::modpow.
        let mut rng = StdRng::seed_from_u64(0x5EED_6000);
        for bits in [64u32, 129, 255, 256, 511, 512] {
            for _ in 0..5 {
                let mut m = rng.gen_biguint(bits.into());
                m.set_bit(bits as u64 - 1, true);
                m.set_bit(0, true);
                let ctx = MontyCtx100::new(&m).expect("odd, fits");
                let a = rng.gen_biguint_range(&BigUint::from(0u32), &m);
                let d = rng.gen_biguint(bits.into());
                let dl = to_limbs_padded(&d, 100);
                let da: &[u64; 100] = dl.as_slice().try_into().unwrap();
                let got = ctx.exit(&ctx.pow_mont(&ctx.enter(&a), da));
                assert_eq!(got, a.modpow(&d, &m), "bits={bits}");
            }
        }
        // Fixed small exponents (window edges: 0, powers of two,
        // window-width boundaries 15/16, zero runs in 256).
        let mut m = rng.gen_biguint(256u64);
        m.set_bit(255, true);
        m.set_bit(0, true);
        let ctx = MontyCtx100::new(&m).expect("odd, fits");
        let a = rng.gen_biguint_range(&BigUint::from(0u32), &m);
        for e in [0u32, 1, 2, 3, 15, 16, 31, 32, 255, 256] {
            let d = BigUint::from(e);
            let dl = to_limbs_padded(&d, 100);
            let da: &[u64; 100] = dl.as_slice().try_into().unwrap();
            let got = ctx.exit(&ctx.pow_mont(&ctx.enter(&a), da));
            assert_eq!(got, a.modpow(&d, &m), "exp={e}");
        }
        // Edges: exp 0/1, base 0/1/m-1, tiny modulus.
        let m3 = BigUint::from(3u32);
        let ctx3 = MontyCtx100::new(&m3).unwrap();
        let one = BigUint::from(1u32);
        let pow = |a: &BigUint, d: &BigUint| {
            let dl = to_limbs_padded(d, 100);
            let da: &[u64; 100] = dl.as_slice().try_into().unwrap();
            ctx3.exit(&ctx3.pow_mont(&ctx3.enter(a), da))
        };
        assert_eq!(pow(&BigUint::from(2u32), &BigUint::from(0u32)), one);
        assert_eq!(
            pow(&BigUint::from(2u32), &one),
            BigUint::from(2u32).modpow(&one, &m3)
        );
        assert_eq!(
            pow(&BigUint::from(0u32), &BigUint::from(5u32)),
            BigUint::from(0u32)
        );
        // Even / oversize moduli fall back (None).
        assert!(MontyCtx100::new(&BigUint::from(4u32)).is_none());
        assert!(MontyCtx100::new(&(BigUint::one() << 6400u64)).is_none());
    }

    #[test]
    #[ignore] // seconds per case in debug; run explicitly in release
    fn monty_modpow_6300_matches() {
        let mut rng = StdRng::seed_from_u64(0x5EED_6300);
        for _ in 0..3 {
            let mut m = rng.gen_biguint(6300u64);
            m.set_bit(6299, true);
            m.set_bit(0, true);
            let ctx = MontyCtx100::new(&m).expect("odd, fits");
            let a = rng.gen_biguint_range(&BigUint::from(0u32), &m);
            let d = rng.gen_biguint(6290u64);
            let dl = to_limbs_padded(&d, 100);
            let da: &[u64; 100] = dl.as_slice().try_into().unwrap();
            let got = ctx.exit(&ctx.pow_mont(&ctx.enter(&a), da));
            assert_eq!(got, a.modpow(&d, &m));
        }
    }

    #[test]
    fn nlimbs_covers_bit_length() {
        assert_eq!(nlimbs(1), 1);
        assert_eq!(nlimbs(64), 1);
        assert_eq!(nlimbs(65), 2);
        assert_eq!(nlimbs(6300), 99);
    }
}
