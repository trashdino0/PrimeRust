//! Prime generation: exact Rust port of `PrimeBenchmark.java`.
//!
//! Three search strategies over the same Miller-Rabin core, mirroring the
//! Java original so benchmark deltas measure implementation, not algorithm:
//!
//! 1. [`method_probable_prime`] — random odd start, sequential `+2` scan
//!    (Java `BigInteger.probablePrime`, certainty 100).
//! 2. [`method_next_prime`] — one random start, `nextProbablePrime` scan
//!    (Java `start.nextProbablePrime()`, certainty 100).
//! 3. [`method_fast_gcd`] — random odds filtered by
//!    `gcd(candidate, SMALL_PRIME_PRODUCT) == 1`, then Miller-Rabin at
//!    certainty 40 (Java `isProbablePrime(40)`).
//!
//! # Faithfulness notes vs Java
//!
//! - Java's `primeToCertainty` runs `rounds = (certainty + 1) / 2`
//!   Miller-Rabin rounds plus a Lucas-Lehmer check. This port runs the same
//!   Miller-Rabin round count ([`rounds_for_certainty`]) but omits the Lucas
//!   test — slightly *less* work per candidate than Java. Any Rust slowdown
//!   measured here is therefore conservative; a win is real.
//! - Java sieves a search window (`largePrime`/`searchSieve`) before
//!   Miller-Rabin. Here the sieve is approximated by trial division inside
//!   [`is_probable_prime`] against a configurable table (default the same 53
//!   odd primes, 3..=251; wider via `--sieve-limit`), plus the explicit GCD
//!   pre-filter in method 3. Same default primes, same order of work,
//!   simpler code.
//! - RNG: Java uses `SecureRandom`; here `rand::thread_rng` (ChaCha-seeded,
//!   CSPRNG). RNG cost is negligible next to 6300-bit modular exponentiation.

use num_bigint::{BigUint, RandBigInt};
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use rand::Rng;
use rayon::prelude::*;
use std::sync::LazyLock;
use thiserror::Error;

use crate::monty::MontyCtx100;

/// Default size from the Java original: 6300 bits ≈ 1897 decimal digits.
pub const DEFAULT_BITS: u32 = 6300;
/// Java `DEFAULT_PRIME_CERTAINTY` used by `probablePrime`/`nextProbablePrime`.
pub const DEFAULT_CERTAINTY: u32 = 100;
/// Java `methodFastGcd` uses `isProbablePrime(40)`.
pub const FAST_GCD_CERTAINTY: u32 = 40;

/// Minimum candidate bit-length for the Montgomery fast path. Measured
/// per-round crossover (identical inputs, release): 4608 → 1.76x loss,
/// 5120 → 1.35x, 5632 → 1.14x, 5888 → 1.03x, 6144 → 0.87x, 6300 → 0.81x.
/// Fixed-100 work scales with exponent bits (linear) while num-bigint
/// scales superlinearly, so the win grows with size; 6144 keeps ≥13%
/// margin. Everything below runs the proven BigUint path untouched.
const MONTY_MIN_BITS: u32 = 6144;

/// The exact 53 odd primes (3..=251) from the Java static initializer.
///
/// Default trial-division table (see `--sieve-limit`): wider tables are built
/// by [`primes_up_to`]. The GCD pre-filter in [`method_fast_gcd`] always uses
/// this 53-prime product — a wider product would slow every `gcd` without
/// changing filter power (trial division already covers the wider range).
pub const SMALL_PRIMES: [u32; 53] = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
];

/// Sieve of Eratosthenes: all odd primes in `[3, limit]`.
///
/// `limit < 3` yields an empty table. `primes_up_to(251)` reproduces
/// [`SMALL_PRIMES`] exactly (asserted by test), so the default benchmark path
/// stays bit-identical to the Java port while wider limits (1000, 10000) are
/// available for head-to-head measurement. Cost is trivial (primes to 10000
/// sieve in ~1 ms) and runs once, outside the timed section.
pub fn primes_up_to(limit: u32) -> Vec<u32> {
    let n = limit as usize;
    if n < 3 {
        return Vec::new();
    }
    let mut sieve = vec![true; n + 1];
    sieve[0] = false;
    sieve[1] = false;
    sieve[2] = false;
    let mut p = 3;
    while p * p <= n {
        if sieve[p] {
            let mut m = p * p;
            while m <= n {
                sieve[m] = false;
                m += 2 * p;
            }
        }
        p += 2;
    }
    (3..=n)
        .step_by(2)
        .filter(|&i| sieve[i])
        .map(|i| i as u32)
        .collect()
}

/// Product of [`SMALL_PRIMES`], computed once and shared by reference.
/// `LazyLock` gives exactly-once initialization across threads (fibonacci-rust
/// used the same pattern for NTT root tables).
pub static SMALL_PRIME_PRODUCT: LazyLock<BigUint> = LazyLock::new(|| {
    SMALL_PRIMES
        .iter()
        .fold(BigUint::one(), |acc, &p| acc * BigUint::from(p))
});

#[derive(Error, Debug, PartialEq, Eq)]
pub enum PrimeError {
    #[error("bit length must be >= 2, got {0}")]
    TooSmall(u32),
}

/// Work counters for one prime search, reported outside the timed section.
///
/// Distinguishes "did more work" (high `candidates`, search luck) from "did
/// work more slowly" (high time per `mr_rounds`). Increments are plain
/// integer adds — negligible next to 6300-bit modular exponentiation.
#[derive(Default, Debug)]
pub struct SearchStats {
    /// Odd candidates examined (drawn or scanned).
    pub candidates: u64,
    /// Candidates rejected by trial division before Miller-Rabin.
    pub trial_rejects: u64,
    /// Candidates that reached Miller-Rabin.
    pub mr_tests: u64,
    /// Total Miller-Rabin rounds executed.
    pub mr_rounds: u64,
}

impl SearchStats {
    /// Fold another search's counters into this one (parallel batch merging).
    ///
    /// Plain integer adds — negligible next to modular exponentiation.
    pub fn accumulate(&mut self, other: &SearchStats) {
        self.candidates += other.candidates;
        self.trial_rejects += other.trial_rejects;
        self.mr_tests += other.mr_tests;
        self.mr_rounds += other.mr_rounds;
    }
}

/// Miller-Rabin rounds for a Java-style certainty value.
///
/// Mirrors OpenJDK `primeToCertainty`: `rounds = (min(certainty, MAX) + 1) / 2`.
/// Certainties here are tiny (≤ 100), so no clamping is needed.
pub fn rounds_for_certainty(certainty: u32) -> u32 {
    certainty.div_ceil(2)
}

/// Random odd `BigUint` with the top bit set (exact `bits` bit-length).
///
/// Borrows `rng` mutably; returns an owned candidate. No cloning involved.
pub fn random_odd_with_bit_length(rng: &mut impl Rng, bits: u32) -> BigUint {
    debug_assert!(bits >= 2, "caller must check bit length first");
    let mut candidate = rng.gen_biguint(bits as u64);
    candidate.set_bit(bits as u64 - 1, true);
    candidate.set_bit(0, true);
    candidate
}

/// Miller-Rabin primality test with Java-equivalent round count.
///
/// - Returns `false` only for definitely-composite inputs (trial division or
///   a Miller-Rabin witness).
/// - Returns `true` for probably-prime inputs; error probability per Java
///   docs is at most `2^-certainty` (this port omits the Lucas second check,
///   so treat the bound as approximate — see module docs).
/// - `certainty == 0` returns `true`, matching `BigInteger.isProbablePrime`.
/// - `small` is the trial-division table (default [`SMALL_PRIMES`]; wider via
///   [`primes_up_to`]). A wider table only rejects more composites before the
///   expensive Miller-Rabin phase.
///
/// Outcome of trial division against 2 and the small-prime table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TrialVerdict {
    /// The candidate IS a small prime (only reachable for tiny inputs).
    Prime,
    /// Definitely composite (even or divisible by a table prime).
    Composite,
    /// No small factor: Miller-Rabin must decide.
    Unknown,
}

/// Trial-divide `candidate` against 2 and `small`.
///
/// Pure function (no stats) shared by the serial path and the parallel
/// fill-stage pre-filter, so both apply byte-identical verdicts. Uses the
/// single-limb `% u32` remainder (measured identical to the `BigUint`
/// form — num-bigint already takes the narrow path) and skips the equality
/// allocation for large candidates (a number over 32 bits cannot equal a
/// table prime). Measured cost: ~1.3 ms per 6300-bit candidate at the 10k
/// table — NOT free (earlier "free" claims rested on a circular metric;
///
/// see docs.md §4 correction), ~1% of a 135 ms MR round.
fn trial_division(candidate: &BigUint, small: &[u32]) -> TrialVerdict {
    if candidate.is_even() {
        if *candidate == BigUint::from(2u32) {
            return TrialVerdict::Prime;
        }
        return TrialVerdict::Composite;
    }
    // Large candidates can never equal a table prime: skip the allocation.
    let tiny = candidate.bits() <= 32;
    for &p in small {
        if tiny && *candidate == BigUint::from(p) {
            return TrialVerdict::Prime;
        }
        if (candidate % p).is_zero() {
            return TrialVerdict::Composite;
        }
    }
    TrialVerdict::Unknown
}

/// Odds per bulk-sieve window in the sequential scanners.
///
/// One `% p` remainder per prime sieves the whole window (see
/// [`sieve_window`); 4096 keeps the transient window small (~3 MB at 6300
/// bits) while making the per-candidate sieving cost negligible.
pub const WINDOW_ODDS: usize = 4096;

/// Bulk trial-division sieve for a contiguous run of ascending odds.
///
/// `window` must hold ascending odds spaced by 2 (exactly what the
/// sequential scanners build). Each odd prime `p` divides every `p`-th
/// window slot, so a single `% p` remainder on `window[0]` plus index
/// arithmetic marks the whole window — ~4000x fewer divisions than
/// per-candidate [`trial_division`] at [`WINDOW_ODDS`].
///
/// Returns `(certain, need_mr, rejects)` in window order: `certain` holds
/// window values that ARE table primes (only possible for tiny inputs; every
/// value `<= limit` is either a table prime or has a table factor, so all
/// certain values precede all `need_mr` values), `need_mr` holds survivors
/// for Miller-Rabin, `rejects` counts composites plus values `< 2` (mirroring
/// [`is_probable_prime`]'s pre-trial `< 2` check, which `trial_division`
/// alone does not do). Verdicts match `trial_division` slot-for-slot
/// (asserted by test); callers keep identical `SearchStats` counting.
fn sieve_window(window: Vec<BigUint>, small: &[u32]) -> (Vec<BigUint>, Vec<BigUint>, u64) {
    let n = window.len();
    if n == 0 {
        return (Vec::new(), Vec::new(), 0);
    }
    // Real searches (bits >= 64) never touch values <= u32::MAX, so every
    // per-slot special case below is skipped there via short-circuit.
    let tiny = window[0].bits() <= 32;
    let mut dead = vec![false; n];
    let mut certain = vec![false; n];
    if tiny {
        for (i, c) in window.iter().enumerate() {
            if c.bits() <= 1 {
                dead[i] = true; // 0 or 1: not prime
            }
        }
    }
    for &p in small {
        let p64 = u64::from(p);
        // One division per prime for the whole window.
        let rem = (&window[0] % p)
            .to_u32()
            .expect("remainder < p always fits u32") as u64;
        // Slot i holds base + 2i; (rem + 2i) % p == 0 gives
        // i == (-rem) * inv2 (mod p), inv2 = (p + 1) / 2 for odd p.
        let inv2 = p64.div_ceil(2);
        let mut idx = ((p64 - rem % p64) % p64 * inv2) % p64;
        while idx < n as u64 {
            let i = idx as usize;
            if tiny && window[i].to_u32() == Some(p) {
                certain[i] = true; // IS a table prime: never dead (nothing else divides it)
            } else {
                dead[i] = true;
            }
            idx += p64;
        }
    }
    let mut primes = Vec::new();
    let mut survivors = Vec::new();
    let mut rejects = 0u64;
    for (i, c) in window.into_iter().enumerate() {
        if certain[i] {
            primes.push(c);
        } else if dead[i] {
            rejects += 1;
        } else {
            survivors.push(c);
        }
    }
    (primes, survivors, rejects)
}

/// Miller-Rabin primality test with Java-equivalent round count.
///
/// - Returns `false` only for definitely-composite inputs (trial division or
///   a Miller-Rabin witness).
/// - Returns `true` for probably-prime inputs; error probability per Java
///   docs is at most `2^-certainty` (this port omits the Lucas second check,
///   so treat the bound as approximate — see module docs).
/// - `certainty == 0` returns `true`, matching `BigInteger.isProbablePrime`.
/// - `small` is the trial-division table (default [`SMALL_PRIMES`]; wider via
///   [`primes_up_to`]). A wider table only rejects more composites before the
///   expensive Miller-Rabin phase.
///
/// Work is recorded in `stats` for out-of-band diagnosis (see [`SearchStats`]).
pub fn is_probable_prime(
    candidate: &BigUint,
    certainty: u32,
    rng: &mut impl Rng,
    stats: &mut SearchStats,
    small: &[u32],
) -> bool {
    if certainty == 0 {
        return true;
    }
    if *candidate < BigUint::from(2u32) {
        return false;
    }
    // Trial division against 2 and the small-prime table. This doubles as the
    // Java sieve approximation and makes tiny inputs exact.
    match trial_division(candidate, small) {
        TrialVerdict::Prime => return true,
        TrialVerdict::Composite => {
            stats.trial_rejects += 1;
            return false;
        }
        TrialVerdict::Unknown => {}
    }
    miller_rabin(candidate, certainty, rng, stats)
}

/// Miller-Rabin witness loop for a trial-division survivor.
///
/// `stats.mr_tests` counts here so pre-sieved (windowed) callers share the
/// same counters as [`is_probable_prime`]. The caller must have excluded
/// inputs `< 2`, even inputs, and `certainty == 0` (all handled before the
/// trial verdict in `is_probable_prime`); window sieves only emit large odds.
fn miller_rabin(
    candidate: &BigUint,
    certainty: u32,
    rng: &mut impl Rng,
    stats: &mut SearchStats,
) -> bool {
    debug_assert!(certainty > 0, "certainty 0 short-circuits before MR");
    debug_assert!(candidate.is_odd(), "evens never reach MR");
    stats.mr_tests += 1;

    // Write n - 1 = 2^s * d with d odd.
    let one = BigUint::one();
    let n_minus_one = candidate - &one;
    let mut d = n_minus_one.clone();
    let mut s: u32 = 0;
    while d.is_even() {
        d >>= 1;
        s += 1;
    }

    let two = BigUint::from(2u32);
    // Random bases in [2, n - 1), i.e. up to n - 2 inclusive.
    let high = candidate - &one;
    let rounds = rounds_for_certainty(certainty);
    // Fast path: one Montgomery context per tested candidate (amortized
    // over the rounds); bases enter the domain per round, tail comparisons
    // stay in-domain — no exit() in the loop. Falls back transparently for
    // even/too-large moduli (constructibility) — and below MONTY_MIN_BITS
    // (measured crossover; fixed-100 work only pays at large sizes).
    if candidate.bits() >= u64::from(MONTY_MIN_BITS) {
        if let Some(ctx) = MontyCtx100::new(candidate) {
            let d_m = crate::monty::limbs_100(&d);
            let nm1_m = ctx.enter(&n_minus_one);
            'witness: for _ in 0..rounds {
                stats.mr_rounds += 1;
                let a = rng.gen_biguint_range(&two, &high);
                let mut x = ctx.pow_mont(&ctx.enter(&a), &d_m);
                if x == *ctx.mont_one() || x == nm1_m {
                    continue;
                }
                for _ in 1..s {
                    x = ctx.square(&x);
                    if x == nm1_m {
                        continue 'witness;
                    }
                }
                return false;
            }
            return true;
        }
    }
    'witness: for _ in 0..rounds {
        stats.mr_rounds += 1;
        let a = rng.gen_biguint_range(&two, &high);
        let mut x = a.modpow(&d, candidate);
        if x == one || x == n_minus_one {
            continue;
        }
        for _ in 1..s {
            x = (&x * &x) % candidate;
            if x == n_minus_one {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// 1. Gold standard: random odd start, sequential `+2` scan at certainty 100.
///
/// Faithful to `BigInteger.probablePrime(BIT_LENGTH, random)`: if the scan
/// overflows the bit length (rare), a fresh random start is drawn, matching
/// Java's regenerate-on-overflow behavior.
pub fn method_probable_prime(
    bits: u32,
    certainty: u32,
    rng: &mut impl Rng,
    stats: &mut SearchStats,
    small: &[u32],
) -> Result<BigUint, PrimeError> {
    if bits < 2 {
        return Err(PrimeError::TooSmall(bits));
    }
    let two = BigUint::from(2u32);
    loop {
        let mut candidate = random_odd_with_bit_length(rng, bits);
        loop {
            // Bulk-sieve one window of contiguous odds in scan order; the
            // survivors are Miller-Rabin-tested in the same order a +2 scan
            // would visit them. A partial window means the scan overflowed
            // the bit length: test its survivors, then draw a fresh random
            // start like Java's regenerate-on-overflow.
            let mut window = Vec::with_capacity(WINDOW_ODDS);
            while window.len() < WINDOW_ODDS && candidate.bits() == bits as u64 {
                stats.candidates += 1;
                // Genuine clone: the window owns its items while the scan
                // position advances (noise next to a millisecond MR round).
                window.push(candidate.clone());
                candidate += &two;
            }
            if window.is_empty() {
                break;
            }
            let full = window.len() == WINDOW_ODDS;
            let (certain, survivors, rejects) = sieve_window(window, small);
            stats.trial_rejects += rejects;
            // Certain values (tiny inputs only) precede every survivor, so
            // the first one is the scan's smallest non-composite: no MR needed.
            if let Some(prime) = certain.into_iter().next() {
                return Ok(prime);
            }
            for survivor in &survivors {
                if miller_rabin(survivor, certainty, rng, stats) {
                    return Ok(survivor.clone());
                }
            }
            if !full {
                break;
            }
        }
    }
}

/// Smallest probably-prime `> start` (odd steps), mirroring
/// `BigInteger.nextProbablePrime()` at Java's default certainty.
///
/// Note the strict inequality: Java returns the first prime *greater than*
/// `this`, so a prime input still advances to the next prime.
pub fn next_probable_prime_from(
    start: &BigUint,
    certainty: u32,
    rng: &mut impl Rng,
    stats: &mut SearchStats,
    small: &[u32],
) -> BigUint {
    let mut candidate = start.clone() + BigUint::one();
    if candidate.is_even() {
        candidate += BigUint::one();
    }
    let two = BigUint::from(2u32);
    loop {
        // Same windowed scan as method_probable_prime: ascending order is
        // preserved, so the first prime found is still the smallest one
        // above the start.
        let mut window = Vec::with_capacity(WINDOW_ODDS);
        while window.len() < WINDOW_ODDS {
            stats.candidates += 1;
            window.push(candidate.clone());
            candidate += &two;
        }
        let (certain, survivors, rejects) = sieve_window(window, small);
        stats.trial_rejects += rejects;
        if let Some(prime) = certain.into_iter().next() {
            return prime;
        }
        for survivor in &survivors {
            if miller_rabin(survivor, certainty, rng, stats) {
                return survivor.clone();
            }
        }
    }
}

/// 2. Linear scanner: one random start (top bit set, parity as drawn), then
///    `nextProbablePrime` scan at certainty 100.
pub fn method_next_prime(
    bits: u32,
    certainty: u32,
    rng: &mut impl Rng,
    stats: &mut SearchStats,
    small: &[u32],
) -> Result<BigUint, PrimeError> {
    if bits < 2 {
        return Err(PrimeError::TooSmall(bits));
    }
    let mut start = rng.gen_biguint(bits as u64);
    start.set_bit(bits as u64 - 1, true);
    Ok(next_probable_prime_from(
        &start, certainty, rng, stats, small,
    ))
}

/// 3. Custom sieve: random odds pre-filtered by
///    `gcd(candidate, SMALL_PRIME_PRODUCT) == 1`, then Miller-Rabin at
///    certainty 40. Exact port of Java `methodFastGcd`.
pub fn method_fast_gcd(
    bits: u32,
    rng: &mut impl Rng,
    stats: &mut SearchStats,
    small: &[u32],
) -> Result<BigUint, PrimeError> {
    if bits < 2 {
        return Err(PrimeError::TooSmall(bits));
    }
    let one = BigUint::one();
    loop {
        let candidate = random_odd_with_bit_length(rng, bits);
        stats.candidates += 1;
        if candidate.gcd(&SMALL_PRIME_PRODUCT) != one {
            continue;
        }
        if is_probable_prime(&candidate, FAST_GCD_CERTAINTY, rng, stats, small) {
            return Ok(candidate);
        }
    }
}

/// Batch oversubscription factor for parallel search: each batch holds
/// `jobs × PAR_BATCH_MULTIPLE` candidates, not `jobs`.
///
/// Measured cause (6300 bits): with one item per thread there is nothing to
/// steal, and the ~40% of items that trial-reject leave their threads idle
/// at the join while others grind full MR tests (2.8× on 6 cores instead of
/// ~5×). Oversubscription feeds work-stealing and amortizes the join.
/// Bounded speculation: at most one batch of extra tests per run's end.
pub const PAR_BATCH_MULTIPLE: usize = 8;

/// Test one batch of pre-sieved candidates on `pool`, returning the first
/// prime in batch order plus honestly merged work counters.
///
/// The batch runs as TWO sequential half-joins: if the first half already
/// holds a prime, the second half never runs (and its work is never
/// counted — bounded speculation shrinks to half a batch). Result is
/// identical to testing the whole batch (a prime in the first half always
/// wins batch order); only the final batch of a search can hit, saving
/// ~a quarter batch per run for one extra join (~µs) per batch.
fn test_batch(
    pool: &rayon::ThreadPool,
    batch: &[BigUint],
    certainty: u32,
) -> (Option<BigUint>, SearchStats) {
    let mid = batch.len().div_ceil(2);
    let (first, rest) = batch.split_at(mid);
    let (found, merged) = test_chunk(pool, first, certainty);
    if found.is_some() {
        return (found, merged);
    }
    let (found_rest, merged_rest) = test_chunk(pool, rest, certainty);
    let mut total = merged;
    total.accumulate(&merged_rest);
    (found_rest, total)
}

/// One join of [`test_batch`]: test every item even after a prime is found
/// elsewhere in the chunk (bounded speculation of one chunk); all work is
/// counted, so parallel `Work:` lines may exceed what a serial scan would
/// do. First-found wins: a valid prime of exact bit length at the stated
/// certainty, but not necessarily the lowest one above the start.
///
/// Items must already have passed trial division, so they go straight to
/// [`miller_rabin`]. Rayon runs at exactly one level here — across
/// candidates — and never inside `modpow` (`learned.md` §3: nested joins
/// hung past 600 s).
fn test_chunk(
    pool: &rayon::ThreadPool,
    batch: &[BigUint],
    certainty: u32,
) -> (Option<BigUint>, SearchStats) {
    // `install` routes this `par_iter` onto `pool`, not the global pool.
    // Each worker draws from its own thread-local RNG; per-item stats merge
    // after the join (no shared mutable state across threads).
    // Batches are deliberately oversubscribed (see PAR_BATCH_MULTIPLE): with
    // one item per thread there is nothing to steal, and threads whose items
    // trial-reject idle at the join while others grind full MR tests.
    // Oversubscription lets rayon balance MR-heavy items via work-stealing.
    let results: Vec<(bool, SearchStats)> = pool.install(|| {
        batch
            .par_iter()
            .map(|candidate| {
                let mut rng = rand::thread_rng();
                let mut stats = SearchStats::default();
                let prime = miller_rabin(candidate, certainty, &mut rng, &mut stats);
                (prime, stats)
            })
            .collect()
    });
    let mut merged = SearchStats::default();
    let mut found = None;
    for (candidate, (prime, stats)) in batch.iter().zip(results) {
        merged.accumulate(&stats);
        if found.is_none() && prime {
            found = Some(candidate.clone());
        }
    }
    (found, merged)
}

/// 1P. Parallel gold standard: sequential `+2` scan in bulk, first-found
/// prime wins.
///
/// Same validity contract as [`method_probable_prime`] (exact bit length,
/// certainty-100 Miller-Rabin); unlike the serial scan the winner is not
/// necessarily the lowest prime above the start. Bit-length overflow draws a
/// fresh random start, mirroring the serial behavior. The scan fills raw
/// odds serially (clone+add are microseconds) and [`sieve_window`] keeps
/// only MR-bound survivors, so test batches run at ~100% density.
pub fn method_probable_prime_parallel(
    bits: u32,
    certainty: u32,
    rng: &mut impl Rng,
    stats: &mut SearchStats,
    small: &[u32],
    pool: &rayon::ThreadPool,
    batch: usize,
) -> Result<BigUint, PrimeError> {
    if bits < 2 {
        return Err(PrimeError::TooSmall(bits));
    }
    // Scan density at the 10k table is ~13%: over-scan to fill the batch.
    let scan = batch.saturating_mul(8);
    let two = BigUint::from(2u32);
    let mut candidate = random_odd_with_bit_length(rng, bits);
    let mut carry: Vec<BigUint> = Vec::new();
    loop {
        let mut bvec = std::mem::take(&mut carry);
        if bvec.len() < batch {
            let mut raw = Vec::with_capacity(scan);
            while raw.len() < scan {
                if candidate.bits() != bits as u64 {
                    candidate = random_odd_with_bit_length(rng, bits);
                    continue;
                }
                stats.candidates += 1;
                // Genuine clone: the batch owns its items for `par_iter`
                // while the scan position advances (~800 B at 6300 bits,
                // noise next to a 135 ms MR round).
                raw.push(candidate.clone());
                candidate += &two;
            }
            let (certain, survivors, rejects) = sieve_window(raw, small);
            stats.trial_rejects += rejects;
            // `certain` is empty at real bit sizes (tiny inputs only); it
            // rides along in order so the batch stays a contiguous scan.
            bvec.extend(certain);
            bvec.extend(survivors);
        }
        if bvec.len() > batch {
            carry = bvec.split_off(batch);
        }
        // Rare shortfall top-up stays serial (scan + shared trial verdicts).
        while bvec.len() < batch {
            if candidate.bits() != bits as u64 {
                candidate = random_odd_with_bit_length(rng, bits);
                continue;
            }
            stats.candidates += 1;
            let next = candidate.clone();
            candidate += &two;
            match trial_division(&next, small) {
                TrialVerdict::Composite => stats.trial_rejects += 1,
                TrialVerdict::Prime | TrialVerdict::Unknown => bvec.push(next),
            }
        }
        let (found, merged) = test_batch(pool, &bvec, certainty);
        stats.accumulate(&merged);
        if let Some(prime) = found {
            return Ok(prime);
        }
    }
}

/// 3P. Parallel custom sieve: GCD pre-filter fills each batch, then
/// the batch is Miller-Rabin-tested in parallel at certainty 40.
///
/// All stages run across candidates (still one parallel level — sequential
/// stages, each internally parallel, never nested): the draw+GCD+trial
/// stage over-draws in parallel because only ~20% of draws pass the GCD
/// filter, and a serial fill would Amdahl-cap the speedup (measured 2.4×
/// on 6 cores before this change). Draw, GCD, and trial fuse in ONE pass
/// (no second sieve join or intermediate Vec) so test batches run at
/// ~100% density. Surplus survivors carry to the next round; a rare shortfall tops up
/// serially. Draws move into the batch (no clone).
pub fn method_fast_gcd_parallel(
    bits: u32,
    rng: &mut impl Rng,
    stats: &mut SearchStats,
    small: &[u32],
    pool: &rayon::ThreadPool,
    batch: usize,
) -> Result<BigUint, PrimeError> {
    if bits < 2 {
        return Err(PrimeError::TooSmall(bits));
    }
    let one = BigUint::one();
    // Draw→survivor rate is ~13% (20% GCD × ~64% trial): over-draw wide so
    // the serial top-up below stays rare.
    let overdraw = batch.saturating_mul(12);
    let mut carry: Vec<BigUint> = Vec::new();
    loop {
        // Surplus survivors from last round's over-draw go first.
        let mut bvec = std::mem::take(&mut carry);
        if bvec.len() < batch {
            // Fused draw+GCD+trial in ONE parallel pass (was: a GCD pass,
            // then a second sieve pass + join + intermediate Vec). Same
            // filters in the same order, same stats: GCD rejects vanish
            // from the counters exactly as before — only trial-division
            // composites increment trial_rejects.
            let verdicts: Vec<(Option<BigUint>, bool)> = pool.install(|| {
                (0..overdraw)
                    .into_par_iter()
                    .map(|_| {
                        let mut rng = rand::thread_rng();
                        let candidate = random_odd_with_bit_length(&mut rng, bits);
                        if candidate.gcd(&SMALL_PRIME_PRODUCT) != one {
                            return (None, false);
                        }
                        match trial_division(&candidate, small) {
                            TrialVerdict::Composite => (None, true),
                            // A `Prime` verdict is only reachable for tiny
                            // inputs; passing it through keeps this exact
                            // for all sizes (test_batch confirms it).
                            TrialVerdict::Prime | TrialVerdict::Unknown => (Some(candidate), false),
                        }
                    })
                    .collect()
            });
            stats.candidates += overdraw as u64;
            for (item, trial_reject) in verdicts {
                match item {
                    Some(candidate) => bvec.push(candidate),
                    None => stats.trial_rejects += u64::from(trial_reject),
                }
            }
        }
        if bvec.len() > batch {
            carry = bvec.split_off(batch);
        }
        // Rare shortfall top-up stays serial (same filters, shared verdicts).
        while bvec.len() < batch {
            let candidate = random_odd_with_bit_length(rng, bits);
            stats.candidates += 1;
            if candidate.gcd(&SMALL_PRIME_PRODUCT) != one {
                continue;
            }
            match trial_division(&candidate, small) {
                TrialVerdict::Composite => stats.trial_rejects += 1,
                TrialVerdict::Prime | TrialVerdict::Unknown => bvec.push(candidate),
            }
        }
        let (found, merged) = test_batch(pool, &bvec, FAST_GCD_CERTAINTY);
        stats.accumulate(&merged);
        if let Some(prime) = found {
            return Ok(prime);
        }
    }
}

/// Cheap structural verification run *outside* the timed section.
///
/// Checks exact bit length, oddness, and freedom from small factors. Full
/// primality re-testing would double the benchmark cost; the timed search
/// itself already ran Miller-Rabin at the stated certainty.
pub fn verify_structure(candidate: &BigUint, bits: u32) -> Result<(), PrimeError> {
    if bits < 2 {
        return Err(PrimeError::TooSmall(bits));
    }
    assert!(
        candidate.bits() == bits as u64,
        "prime bit-length mismatch: got {}, want {bits}",
        candidate.bits()
    );
    assert!(candidate.bit(0), "prime candidate must be odd");
    for &p in &SMALL_PRIMES {
        let prime = BigUint::from(p);
        if candidate != &prime {
            assert!(
                !((candidate % &prime).is_zero()),
                "prime candidate divisible by small prime {p}"
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn seeded_rng() -> StdRng {
        StdRng::seed_from_u64(0xC0FFEE)
    }

    #[test]
    fn rounds_match_java_prime_to_certainty() {
        // Java: rounds = (certainty + 1) / 2.
        let cases = [(40, 20), (100, 50), (1, 1), (2, 1), (0, 0)];
        for (certainty, expected) in cases {
            assert_eq!(
                rounds_for_certainty(certainty),
                expected,
                "certainty {certainty}"
            );
        }
    }

    #[test]
    fn small_prime_product_matches_java_table() {
        // Divisible by every table prime, not by 2.
        for &p in &SMALL_PRIMES {
            assert_eq!(&*SMALL_PRIME_PRODUCT % p, BigUint::zero(), "factor {p}");
        }
        assert!(SMALL_PRIME_PRODUCT.is_odd());
        assert_eq!(SMALL_PRIMES.len(), 53);
        assert_eq!(SMALL_PRIMES[0], 3);
        assert_eq!(SMALL_PRIMES[52], 251);
    }

    #[test]
    fn miller_rabin_accepts_known_primes() {
        let mut rng = seeded_rng();
        // Table-driven: small primes, a Carmichael-adjacent prime, large known primes.
        let primes = ["2", "3", "5", "251", "104729", "1299709", "32416190071"];
        for s in primes {
            let n: BigUint = s.parse().unwrap();
            let mut stats = SearchStats::default();
            assert!(
                is_probable_prime(&n, 40, &mut rng, &mut stats, &SMALL_PRIMES),
                "{s} should test prime"
            );
        }
        // A large prime above the trial-division table must reach Miller-Rabin.
        let big: BigUint = "32416190071".parse().unwrap();
        let mut stats = SearchStats::default();
        assert!(is_probable_prime(
            &big,
            40,
            &mut seeded_rng(),
            &mut stats,
            &SMALL_PRIMES
        ));
        assert_eq!(stats.mr_tests, 1);
        assert_eq!(stats.mr_rounds, 20);
    }

    #[test]
    fn miller_rabin_rejects_known_composites() {
        let mut rng = seeded_rng();
        // Small composites (trial division catches these before Miller-Rabin).
        let composites = [
            "0",
            "1",
            "4",
            "6",
            "8",
            "9",
            "15",
            "25",
            "561",
            "104730",
            "32416190070",
        ];
        for s in composites {
            let n: BigUint = s.parse().unwrap();
            let mut stats = SearchStats::default();
            assert!(
                !is_probable_prime(&n, 40, &mut rng, &mut stats, &SMALL_PRIMES),
                "{s} should test composite"
            );
        }
        // Carmichael number 561 is divisible by 3, so trial division — not
        // Miller-Rabin — rejects it. Assert that honestly.
        let carmichael: BigUint = "561".parse().unwrap();
        let mut stats = SearchStats::default();
        assert!(!is_probable_prime(
            &carmichael,
            40,
            &mut seeded_rng(),
            &mut stats,
            &SMALL_PRIMES
        ));
        assert_eq!(stats.trial_rejects, 1);
        assert_eq!(stats.mr_tests, 0);
        // Constructed semiprime of two large primes: no factor <= 251, so it
        // sails through trial division and Miller-Rabin itself must catch it.
        // (A Fermat-only test would pass this — this assertion is the one
        // with teeth.)
        let p: BigUint = "104729".parse().unwrap();
        let q: BigUint = "1299709".parse().unwrap();
        let semiprime = &p * &q;
        let mut stats = SearchStats::default();
        assert!(!is_probable_prime(
            &semiprime,
            40,
            &mut seeded_rng(),
            &mut stats,
            &SMALL_PRIMES
        ));
        assert_eq!(stats.trial_rejects, 0);
        assert_eq!(stats.mr_tests, 1);
    }

    #[test]
    fn zero_certainty_returns_true_like_java() {
        let mut rng = seeded_rng();
        let composite = BigUint::from(15u32);
        assert!(is_probable_prime(
            &composite,
            0,
            &mut rng,
            &mut SearchStats::default(),
            &SMALL_PRIMES
        ));
    }

    #[test]
    fn random_odd_has_exact_bit_length_and_oddness() {
        let mut rng = seeded_rng();
        for bits in [2, 3, 16, 65, 128] {
            for _ in 0..10 {
                let c = random_odd_with_bit_length(&mut rng, bits);
                assert_eq!(c.bits(), bits as u64, "bits={bits}");
                assert!(c.bit(0), "must be odd");
            }
        }
    }

    #[test]
    fn next_prime_from_known_start() {
        // Strictly greater than start, like Java nextProbablePrime().
        // Windowed scan: every case resolves inside the first window, so
        // candidates == WINDOW_ODDS. Small-prime hits return as certain
        // (mr_tests == 0); 263 is above the table, so Miller-Rabin itself
        // must confirm it (1 test, 20 rounds at certainty 40).
        let cases = [
            ("14", "17", 0, 0),
            ("20", "23", 0, 0),
            ("2", "3", 0, 0),
            ("3", "5", 0, 0),
            ("260", "263", 1, 20),
        ];
        for (start, expected, tests, rounds) in cases {
            let s: BigUint = start.parse().unwrap();
            let e: BigUint = expected.parse().unwrap();
            let mut stats = SearchStats::default();
            assert_eq!(
                next_probable_prime_from(&s, 40, &mut seeded_rng(), &mut stats, &SMALL_PRIMES),
                e
            );
            assert_eq!(stats.candidates, WINDOW_ODDS as u64, "start {start}");
            assert!(stats.trial_rejects > 0, "start {start}");
            assert_eq!(stats.mr_tests, tests, "start {start}");
            assert_eq!(stats.mr_rounds, rounds, "start {start}");
        }
    }

    #[test]
    fn methods_reject_too_small_bit_length() {
        let mut rng = seeded_rng();
        let mut stats = SearchStats::default();
        assert_eq!(
            method_probable_prime(1, DEFAULT_CERTAINTY, &mut rng, &mut stats, &SMALL_PRIMES),
            Err(PrimeError::TooSmall(1))
        );
        assert_eq!(
            method_next_prime(0, DEFAULT_CERTAINTY, &mut rng, &mut stats, &SMALL_PRIMES),
            Err(PrimeError::TooSmall(0))
        );
        assert_eq!(
            method_fast_gcd(1, &mut rng, &mut stats, &SMALL_PRIMES),
            Err(PrimeError::TooSmall(1))
        );
    }

    #[test]
    fn fast_gcd_finds_small_prime_with_structure() {
        let mut rng = seeded_rng();
        let mut stats = SearchStats::default();
        // 64 bits keeps this test fast; certainty path is identical to 6300-bit.
        let p = method_fast_gcd(64, &mut rng, &mut stats, &SMALL_PRIMES).unwrap();
        assert_eq!(p.bits(), 64);
        verify_structure(&p, 64).unwrap();
        // Every Miller-Rabin test at certainty 40 runs *at most* 20 rounds:
        // composites that reveal a witness exit early.
        assert!(stats.candidates > 0);
        assert!(stats.mr_tests > 0);
        assert!(stats.mr_rounds >= stats.mr_tests);
        assert!(stats.mr_rounds <= stats.mr_tests * 20);
        // Independent re-check with fresh randomness at higher certainty.
        assert!(is_probable_prime(
            &p,
            64,
            &mut seeded_rng(),
            &mut SearchStats::default(),
            &SMALL_PRIMES
        ));
    }

    #[test]
    fn sieve_reproduces_exact_port_table() {
        // Default path must stay bit-identical to the Java port's 53 primes.
        assert_eq!(primes_up_to(251), SMALL_PRIMES);
        assert!(primes_up_to(2).is_empty());
        // Known prime-counting values (2 excluded: evenness is checked separately).
        let to_1k = primes_up_to(1000);
        assert_eq!(to_1k.len(), 167);
        assert_eq!(*to_1k.last().unwrap(), 997);
        let to_10k = primes_up_to(10000);
        assert_eq!(to_10k.len(), 1228);
        assert_eq!(*to_10k.last().unwrap(), 9973);
    }

    #[test]
    fn wider_table_rejects_more_before_miller_rabin() {
        // 1009 x 1013: both factors above 251, so the default table lets it
        // through to Miller-Rabin — but a table to 2000 catches factor 1009.
        let semiprime: BigUint = "1022117".parse().unwrap();
        let mut stats = SearchStats::default();
        assert!(!is_probable_prime(
            &semiprime,
            40,
            &mut seeded_rng(),
            &mut stats,
            &SMALL_PRIMES
        ));
        assert_eq!(stats.trial_rejects, 0);
        assert_eq!(stats.mr_tests, 1);
        let wide = primes_up_to(2000);
        let mut stats = SearchStats::default();
        assert!(!is_probable_prime(
            &semiprime,
            40,
            &mut seeded_rng(),
            &mut stats,
            &wide
        ));
        assert_eq!(stats.trial_rejects, 1);
        assert_eq!(stats.mr_tests, 0);
    }

    #[test]
    fn probable_and_next_find_small_primes_with_structure() {
        let mut rng = seeded_rng();
        let mut stats = SearchStats::default();
        let p = method_probable_prime(64, DEFAULT_CERTAINTY, &mut rng, &mut stats, &SMALL_PRIMES)
            .unwrap();
        assert_eq!(p.bits(), 64);
        verify_structure(&p, 64).unwrap();
        // Certainty 100 means at most 50 rounds per Miller-Rabin test
        // (early witness exit uses fewer).
        assert!(stats.candidates > 0);
        assert!(stats.mr_tests > 0);
        assert!(stats.mr_rounds >= stats.mr_tests);
        assert!(stats.mr_rounds <= stats.mr_tests * 50);

        let mut stats = SearchStats::default();
        let q =
            method_next_prime(64, DEFAULT_CERTAINTY, &mut rng, &mut stats, &SMALL_PRIMES).unwrap();
        assert_eq!(q.bits(), 64);
        verify_structure(&q, 64).unwrap();
        assert!(stats.mr_rounds <= stats.mr_tests * 50);
    }

    #[test]
    fn trial_division_verdicts_are_exact() {
        // Tiny inputs resolve exactly; large ones route by small factors.
        let cases = [
            ("2", TrialVerdict::Prime),
            ("3", TrialVerdict::Prime),
            ("251", TrialVerdict::Prime),
            ("4", TrialVerdict::Composite),
            ("9", TrialVerdict::Composite),
            ("561", TrialVerdict::Composite),
            ("263", TrialVerdict::Unknown),
            ("104729", TrialVerdict::Unknown),
            // Semiprime above the table sails through to Miller-Rabin.
            ("136117223861", TrialVerdict::Unknown),
        ];
        for (s, expected) in cases {
            let n: BigUint = s.parse().unwrap();
            assert_eq!(trial_division(&n, &SMALL_PRIMES), expected, "{s}");
        }
        // 104729 * 1299709 (verified by multiplication, not by memory).
        let p: BigUint = "104729".parse().unwrap();
        let q: BigUint = "1299709".parse().unwrap();
        assert_eq!((&p * &q).to_string(), "136117223861");
        // Width matters: 1009x1013 is Unknown at 251, Composite at 2000.
        let mid: BigUint = "1022117".parse().unwrap();
        assert_eq!(trial_division(&mid, &SMALL_PRIMES), TrialVerdict::Unknown);
        assert_eq!(
            trial_division(&mid, &primes_up_to(2000)),
            TrialVerdict::Composite
        );
    }

    #[test]
    fn sieve_window_matches_trial_division() {
        // Slot-for-slot agreement between the bulk sieve and per-candidate
        // trial division, on tiny, sub-2, 64-bit, and wide-table windows.
        let wide = primes_up_to(10000);
        let windows: Vec<(BigUint, usize, Vec<u32>)> = vec![
            ("3".parse().unwrap(), 600, SMALL_PRIMES.to_vec()),
            ("1".parse().unwrap(), 10, SMALL_PRIMES.to_vec()),
            (
                random_odd_with_bit_length(&mut seeded_rng(), 64),
                512,
                SMALL_PRIMES.to_vec(),
            ),
            (
                random_odd_with_bit_length(&mut seeded_rng(), 128),
                300,
                wide.clone(),
            ),
        ];
        let two = BigUint::from(2u32);
        for (base, count, table) in &windows {
            let mut c = base.clone();
            let mut window = Vec::with_capacity(*count);
            for _ in 0..*count {
                window.push(c.clone());
                c += &two;
            }
            let (certain, survivors, rejects) = sieve_window(window.clone(), table);
            // Certain values precede every survivor, both in window order, so
            // their concatenation must equal the alive slots in order. Values
            // < 2 count as rejects (mirroring is_probable_prime's pre-trial
            // check, which trial_division alone does not do).
            let mut alive = certain.clone();
            alive.extend(survivors.clone());
            let mut expect_alive = Vec::new();
            let mut expect_rejects = 0u64;
            for x in &window {
                let verdict = trial_division(x, table);
                if verdict == TrialVerdict::Composite || x.bits() <= 1 {
                    expect_rejects += 1;
                } else {
                    expect_alive.push(x.clone());
                }
            }
            assert_eq!(alive, expect_alive, "base {base}");
            assert_eq!(rejects, expect_rejects, "base {base}");
            // Certain/need_mr split matches Prime/Unknown exactly.
            let expect_certain: Vec<BigUint> = window
                .iter()
                .filter(|x| trial_division(x, table) == TrialVerdict::Prime)
                .cloned()
                .collect();
            let expect_survivors: Vec<BigUint> = window
                .iter()
                .filter(|x| trial_division(x, table) == TrialVerdict::Unknown && x.bits() > 1)
                .cloned()
                .collect();
            assert_eq!(certain, expect_certain, "base {base}");
            assert_eq!(survivors, expect_survivors, "base {base}");
        }
        // Empty window is empty out.
        let (p, s, r) = sieve_window(Vec::new(), &SMALL_PRIMES);
        assert!(p.is_empty() && s.is_empty() && r == 0);
    }

    #[test]
    fn stats_accumulate_merges_all_counters() {
        let mut total = SearchStats::default();
        let part = SearchStats {
            candidates: 10,
            trial_rejects: 7,
            mr_tests: 3,
            mr_rounds: 5,
        };
        total.accumulate(&part);
        total.accumulate(&part);
        assert_eq!(total.candidates, 20);
        assert_eq!(total.trial_rejects, 14);
        assert_eq!(total.mr_tests, 6);
        assert_eq!(total.mr_rounds, 10);
    }

    fn two_thread_pool() -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap()
    }

    #[test]
    fn parallel_methods_reject_too_small_bit_length() {
        let pool = two_thread_pool();
        let mut rng = seeded_rng();
        let mut stats = SearchStats::default();
        assert_eq!(
            method_probable_prime_parallel(
                1,
                DEFAULT_CERTAINTY,
                &mut rng,
                &mut stats,
                &SMALL_PRIMES,
                &pool,
                2
            ),
            Err(PrimeError::TooSmall(1))
        );
        assert_eq!(
            method_fast_gcd_parallel(0, &mut rng, &mut stats, &SMALL_PRIMES, &pool, 2),
            Err(PrimeError::TooSmall(0))
        );
    }

    #[test]
    fn parallel_methods_find_small_primes_with_structure() {
        // 64 bits keeps this fast; batch/pool shape is identical at 6300 bits.
        let pool = two_thread_pool();
        let mut rng = seeded_rng();
        let mut stats = SearchStats::default();
        let p = method_probable_prime_parallel(
            64,
            DEFAULT_CERTAINTY,
            &mut rng,
            &mut stats,
            &SMALL_PRIMES,
            &pool,
            2,
        )
        .unwrap();
        assert_eq!(p.bits(), 64);
        verify_structure(&p, 64).unwrap();
        assert!(stats.candidates > 0);
        assert!(stats.mr_tests > 0);
        assert!(stats.mr_rounds <= stats.mr_tests * 50);
        // Independent serial re-check at higher certainty.
        assert!(is_probable_prime(
            &p,
            64,
            &mut seeded_rng(),
            &mut SearchStats::default(),
            &SMALL_PRIMES
        ));

        let mut stats = SearchStats::default();
        let q =
            method_fast_gcd_parallel(64, &mut rng, &mut stats, &SMALL_PRIMES, &pool, 2).unwrap();
        assert_eq!(q.bits(), 64);
        verify_structure(&q, 64).unwrap();
        assert!(stats.mr_rounds <= stats.mr_tests * 20);
    }
}
