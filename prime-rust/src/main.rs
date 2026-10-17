//! `prime-rust` — 6300-bit prime generation benchmark (Rust port of `PrimeBenchmark.java`).
//!
//! Harness follows fibonacci-rust conventions: warmed release runs, `-a N`
//! mean / 95% CI / stddev / skewness, verification outside the timed section,
//! `anyhow` context at the binary edge.

mod monty;
mod prime;

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use num_bigint::BigUint;
use prime::{PrimeError, SearchStats, DEFAULT_BITS, DEFAULT_CERTAINTY};
use sha2::{Digest, Sha256};
use std::time::Instant;

/// Sharded per-thread heaps: BigUint `modpow` allocates on every multiply,
/// and the default Windows heap lock serializes parallel MR tests.
/// (Kept because the J2 curve diagnosed allocator contention — see docs.md.)
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Copy, Clone, Debug, ValueEnum)]
enum Method {
    /// Java `BigInteger.probablePrime()` (certainty 100).
    ProbablePrime,
    /// Java `start.nextProbablePrime()` (certainty 100).
    NextPrime,
    /// Java custom GCD filter + `isProbablePrime(40)`.
    FastGcd,
}

#[derive(Parser, Debug)]
#[command(
    name = "prime-rust",
    version,
    about = "6300-bit prime generation benchmark (Rust port of PrimeBenchmark.java)"
)]
struct Args {
    /// Bit length of the primes to generate (Java default 6300).
    #[arg(long, default_value_t = DEFAULT_BITS)]
    bits: u32,

    /// Number of timed runs per method for averaging (default 1).
    #[arg(short = 'a', long = "average", default_value_t = 1)]
    runs: u32,

    /// Warm-up rounds per method (timing ignored, default 0 like Java).
    #[arg(short = 'w', long = "warmup", default_value_t = 0)]
    warmup: u32,

    /// Which method(s) to benchmark.
    #[arg(short = 'm', long = "method", value_enum, default_value_t = Method::FastGcd)]
    method: Method,

    /// Benchmark every method instead of just one.
    #[arg(long, default_value_t = false)]
    all: bool,

    /// Trial-division sieve limit (odd primes 3..=limit). Default 100000
    /// (measured ~25% fewer MR tests than 10000; bulk window sieving keeps
    /// sequential filtering near-free, see docs.md). Use 251 for
    /// Java-comparable runs.
    #[arg(long, default_value_t = 100000)]
    sieve_limit: u32,

    /// Parallel candidate-testing threads (rayon, outermost level only).
    /// 1 = serial baseline (exact-port behavior); N > 1 tests batches of
    /// 8*N candidates at once (oversubscribed for work-stealing, see
    /// PAR_BATCH_MULTIPLE), first-found prime wins. next-prime stays serial
    /// (smallest-prime semantics can't parallelize without changing meaning).
    #[arg(long, default_value_t = 1)]
    jobs: usize,
}

/// Two-sided 95% t critical values, df = 1..=30, then normal 1.96.
/// (Same table as fibonacci-rust for comparable statistics.)
fn t_crit_95(df: usize) -> f64 {
    const T: [f64; 30] = [
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160,
        2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056,
        2.052, 2.048, 2.045, 2.042,
    ];
    if df == 0 {
        f64::NAN
    } else if df <= 30 {
        T[df - 1]
    } else {
        1.96
    }
}

fn short_hash(candidate: &BigUint) -> String {
    let mut hasher = Sha256::new();
    hasher.update(candidate.to_bytes_be());
    let digest = hasher.finalize();
    hex_of(&digest[..8])
}

fn hex_of(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn benchmark(
    name: &str,
    bits: u32,
    warmup: u32,
    runs: u32,
    mut f: impl FnMut(&mut SearchStats) -> Result<BigUint, PrimeError>,
) -> Result<()> {
    for _ in 0..warmup {
        let _ = f(&mut SearchStats::default()).context("warm-up generation failed")?;
    }

    let mut times: Vec<f64> = Vec::with_capacity(runs.max(1) as usize);
    let mut total = SearchStats::default();
    let mut last: Option<BigUint> = None;
    for i in 0..runs.max(1) {
        let mut stats = SearchStats::default();
        let start = Instant::now();
        let prime = f(&mut stats).context("prime generation failed")?;
        let secs = start.elapsed().as_secs_f64();
        times.push(secs);
        total.candidates += stats.candidates;
        total.trial_rejects += stats.trial_rejects;
        total.mr_tests += stats.mr_tests;
        total.mr_rounds += stats.mr_rounds;
        last = Some(prime);
        if runs > 1 {
            println!(
                "  Run {}: {:.3} s (candidates: {}, MR tests: {}, MR rounds: {})",
                i + 1,
                secs,
                stats.candidates,
                stats.mr_tests,
                stats.mr_rounds
            );
        }
    }

    let mean = times.iter().sum::<f64>() / times.len() as f64;
    // Verification + reporting outside the timed section (like Java's
    // `prime.toString().length()` after the clock stops).
    let prime = last.context("no prime generated")?;
    prime::verify_structure(&prime, bits).context("generated candidate failed verification")?;
    let digits = prime.to_string().len();

    println!("Method: {name}");
    if times.len() > 1 {
        let n = times.len() as f64;
        let var = times.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / (n - 1.0);
        let sd = var.sqrt();
        let moe = t_crit_95(times.len() - 1) * (sd / n.sqrt());
        let skew = if sd > 0.0 {
            times.iter().map(|t| ((t - mean) / sd).powi(3)).sum::<f64>() / n
        } else {
            0.0
        };
        println!("Average Time: {:.3} ms", mean * 1000.0);
        println!("Runs:     {}", times.len());
        println!("Mean:     {mean:.3} s");
        println!("95% CI:   [{:.3} s, {:.3} s]", mean - moe, mean + moe);
        println!("StdDev:   {sd:.3} s");
        println!("Skewness: {skew:.3}");
    } else {
        println!("Average Time: {:.3} ms", mean * 1000.0);
    }
    println!("Digits: {digits}");
    println!("SHA256(bytes)[..16]: {}", short_hash(&prime));
    println!(
        "Work: {} candidates, {} trial rejects, {} MR tests, {} MR rounds ({:.2} ms/MR round)",
        total.candidates,
        total.trial_rejects,
        total.mr_tests,
        total.mr_rounds,
        if total.mr_rounds > 0 {
            mean * 1000.0 * runs.max(1) as f64 / total.mr_rounds as f64
        } else {
            0.0
        }
    );
    println!("--------------------------------------------------");
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.bits < 2 {
        anyhow::bail!("bit length must be >= 2, got {}", args.bits);
    }
    if args.sieve_limit < 3 {
        anyhow::bail!("sieve limit must be >= 3, got {}", args.sieve_limit);
    }
    if args.jobs < 1 {
        anyhow::bail!("jobs must be >= 1, got {}", args.jobs);
    }
    // Built once, outside every timed section; shared by reference.
    let sieve = prime::primes_up_to(args.sieve_limit);
    println!("Sieve: {} primes (3..={})", sieve.len(), args.sieve_limit);

    let selected: Vec<(Method, &'static str)> = if args.all {
        vec![
            (
                Method::ProbablePrime,
                "BigInteger.probablePrime() [Rust port]",
            ),
            (Method::NextPrime, "nextProbablePrime() [Rust port]"),
            (
                Method::FastGcd,
                "Fast GCD Filter + Miller-Rabin [Rust port]",
            ),
        ]
    } else {
        vec![(
            args.method,
            match args.method {
                Method::ProbablePrime => "BigInteger.probablePrime() [Rust port]",
                Method::NextPrime => "nextProbablePrime() [Rust port]",
                Method::FastGcd => "Fast GCD Filter + Miller-Rabin [Rust port]",
            },
        )]
    };

    if args.warmup > 0 {
        println!("Starting warm-up phase (ignoring these times)...");
    }
    // One pool for the whole run when parallel; serial path never touches rayon.
    // Batches are oversubscribed (jobs x multiple) so work-stealing can
    // balance MR-heavy items (see PAR_BATCH_MULTIPLE).
    let pool = if args.jobs > 1 {
        let batch = args.jobs * prime::PAR_BATCH_MULTIPLE;
        println!("Jobs: {} (parallel batches of {})", args.jobs, batch);
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.jobs)
                .build()
                .context("failed to build rayon thread pool")?,
        )
    } else {
        None
    };
    let mut rng = rand::thread_rng();
    for (method, name) in selected {
        match method {
            Method::ProbablePrime => {
                benchmark(
                    name,
                    args.bits,
                    args.warmup,
                    args.runs,
                    |stats| match &pool {
                        Some(pool) => prime::method_probable_prime_parallel(
                            args.bits,
                            DEFAULT_CERTAINTY,
                            &mut rng,
                            stats,
                            &sieve,
                            pool,
                            args.jobs * prime::PAR_BATCH_MULTIPLE,
                        ),
                        None => prime::method_probable_prime(
                            args.bits,
                            DEFAULT_CERTAINTY,
                            &mut rng,
                            stats,
                            &sieve,
                        ),
                    },
                )?
            }
            Method::NextPrime => {
                if pool.is_some() {
                    eprintln!(
                        "note: next-prime stays serial (smallest-prime semantics); ignoring --jobs"
                    );
                }
                benchmark(name, args.bits, args.warmup, args.runs, |stats| {
                    prime::method_next_prime(args.bits, DEFAULT_CERTAINTY, &mut rng, stats, &sieve)
                })?
            }
            Method::FastGcd => {
                benchmark(
                    name,
                    args.bits,
                    args.warmup,
                    args.runs,
                    |stats| match &pool {
                        Some(pool) => prime::method_fast_gcd_parallel(
                            args.bits,
                            &mut rng,
                            stats,
                            &sieve,
                            pool,
                            args.jobs * prime::PAR_BATCH_MULTIPLE,
                        ),
                        None => prime::method_fast_gcd(args.bits, &mut rng, stats, &sieve),
                    },
                )?
            }
        }
    }
    if args.warmup > 0 {
        println!("Warm-up complete.");
    }
    Ok(())
}
