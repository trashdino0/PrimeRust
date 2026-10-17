# prime-rust — decision log

Benchmark method (from `learned.md`): warmed release runs, one change at a
time, full-output verification after every kept change. Random primes differ
every run, so "hash beside every number" here means the per-run SHA256 prefix
logged by the harness plus structural verification (exact bit length, oddness,
no small factor) outside the timed section.

## Baseline (exact port, generic CPU, no PGO) — 6300 bits, `-a 1`

| Method | Java (run 1) | Rust (run 1) | Rust (run 2, with `Work:` counters) | Notes |
|---|---|---|---|---|
| `probablePrime` (cert. 100) | 10.36 s | 10.68 s | **114.95 s** — 4087 cand, 3264 trial rejects, 823 MR tests, 872 MR rounds, 131.82 ms/round | Run 1 was a ~10% lucky tail (see 4). |
| `nextProbablePrime` (cert. 100) | 28.93 s | — | **65.84 s** — 2204 cand, 1753 trial rejects, 451 MR tests, 500 MR rounds, 131.68 ms/round | First Rust run; ≈ geometric mean luck. |
| Fast GCD + MR (cert. 40) | **DNF, > 560 s** | 59.85 s | **110.44 s** — 4013 cand, 0 trial rejects, 811 MR tests, 830 MR rounds, 133.06 ms/round | 2× run-to-run variance; same work unit cost (see 4). |

Digits all 1897 (≈ 6300 / log2(10)), as expected. Small-bits harness check:
256-bit `--all -a 2` runs in ~3–5 ms/method.

### 1. Kept: exact-port baseline (parity on `probablePrime`)

10.68 s vs 10.36 s (n = 1 each, noisy box) is indistinguishable. The port —
`num-bigint` + `RandBigInt`, Miller-Rabin at Java's `(certainty+1)/2` rounds,
53-prime trial division, GCD pre-filter — is validated as algorithmically
faithful. Deviations, all conservative (less work than Java): Lucas-Lehmer
second check omitted; sieve window approximated by trial division;
ChaCha-`thread_rng` instead of `SecureRandom`.

### 2. Finding: Java `methodFastGcd` is broken-by-design at 6300 bits (DNF)

The Java custom filter runs a full `gcd` + `isProbablePrime(40)` (20 MR
rounds **plus Lucas**) on *every* random candidate (~4400 expected), while
`probablePrime`'s internal sieve rejects almost all candidates with cheap
trial division and reserves MR+Lucas for a handful of survivors. Same
algorithm in Rust finishes in 59.85 s (no Lucas, cheaper modpow path) —
~10× faster than Java's version, which never finished. Lesson, same as
fibonacci Round 3: admire-then-measure; this one measured DNF. Do not "fix"
the Java; the Rust backlog (below) makes the fast path fast properly with a
window sieve.

### 3. Anomaly (logged, not averaged away): Java run 2 stalled > 180 s on method 1

Same `java -jar prime.jar`, same box, minutes later: method 1
(`probablePrime`, 10.36 s in run 1) had not finished after 180 s; no stray
process survived the timeout (`Get-Process java` clean, twice). Box was
verified quiet immediately after (Rust 256-bit: 6.6 ms, normal). Suspects:
`SecureRandom` entropy/seed stall on Windows or transient contention — NOT
search luck (18× mean has probability ≈ e^-18 for a geometric search).
Policy impact: no Java number below stands on n = 1; repeats must interleave
A/B on a quiet box before any optimization verdict. (Resolved in 4: the
stall was an ordinary geometric tail, not an entropy stall — P(T>180 s) ≈
19% at the re-measured ~110 s mean. The `SecureRandom` suspect is dropped.)

### 4. Resolved: the fastGcd anomaly was two stacked errors, not a slowdown

Re-measure with `Work:` counters (runs `rust-6300-fastgcd2.txt`,
`rust-6300-probable2.txt`, `rust-6300-next2.txt`, quiet box, back-to-back)
shows all three methods run the **same work at the same unit cost**:

- ms/MR round is method-independent: 133.06 / 131.82 / 131.68. One number.
- `MR rounds × ms/round = 100%` of wall time in all three runs
  (830×133.06 ms = 110.44 s; 872×131.82 ms = 114.95 s; 500×131.68 ms =
  65.84 s). [Correction, see §7: this metric is circular — ms/round is
  *defined* as wall/rounds, so it bakes all non-MR time in by construction.
  Microbenched later: trial division ≈ 1.3 ms/candidate (~1% of a round),
  GCD/RNG negligible. The conclusion (MR dominates) stands; the "zero"
  was overstated.]
- Filter power is mathematically identical: trial division rejects
  3264/4087 = 80% of `probablePrime` candidates; the gcd filter passes
  811/4013 = 20.2% — same filter (a factor ≤ 251 either way), same ~820
  MR tests. Methods 1 and 3 do literally the same work; only the final
  prime's rounds differ (50 vs 20 ≈ 4 s, ~3.5%).
- `nextPrime` was faster (65.84 s) purely by search luck: 2204 candidates
  ≈ the geometric mean (gap ≈ 2184 odds); the other two runs drew ~4000
  (≈ 2× mean, P ≈ e^−2 ≈ 13% each — unremarkable).

Error 1 (model): the session predicted fastGcd ≈ 0.4× `probablePrime`
from full round budgets (20 vs 50 rounds). But Miller-Rabin early witness
exit collapses both sides to ~1.02 rounds/test (830/811, 872/823,
500/451), so both do ~850 rounds per unlucky run. Budget-counting without
early-exit data is the same mistake class as fibonacci's admired-not-measured
wins. Correct model: **time ≈ MR_rounds × 132 ms; MR_rounds ≈ MR_tests +
49 (final prime); MR_tests ≈ candidates × P(no factor ≤ sieve limit).**

Error 2 (baseline): run-1's 10.68 s (≈ 80 rounds ≈ ~200 candidates + the
50-round prime) was a ~10% lucky tail — P(find within ~250 odds) ≈
1−e^(−250/2184) ≈ 11% — over-trusted at n = 1, exactly the failure
`learned.md` §1/§8 names. Java's run-1 10.36 s was the same luck class;
its run-2 > 180 s stall was the opposite tail. Nothing to explain beyond
geometric variance. Lesson kept: **single samples are direction-only, and
even then only for ≥ 2× gaps; the 10 s "parity" verdict is retracted and
replaced by the run-2 band (~60–115 s across methods and luck).**

Consequence for the backlog: the only levers that move wall time are (a)
**fewer MR tests** (wider sieve — the win is in WIDTH, i.e. rejecting more
candidates, not in making the already-free filtering cheaper), (b) **fewer
rounds on the prime** (lower certainty), (c) **cheaper rounds** (faster
mulmod), (d) **parallel candidates**. Backlog reordered accordingly: #1
(wider table, 5-line change, informs #2's design) before #2.

### 5. Kept: wider trial-division table, `--sieve-limit`, default 10000

`primes_up_to(limit)` (Sieve of Eratosthenes, ~1 ms to 10000, built once
outside the timed section) threaded as `small: &[u32]` through
`is_probable_prime` / all three methods; `--sieve-limit` flag (default
10000, 251 for Java-comparable runs). The gcd pre-filter stays on the
53-prime product deliberately — a wider product slows every `gcd` without
adding filter power trial division doesn't already cover. `verify_structure`
unchanged (251-subset check still valid). Tests: `primes_up_to(251) ==
SMALL_PRIMES` (default path bit-identical to exact port) + width
discrimination via semiprime 1009×1013 (trial-division catch at limit 2000,
MR catch at 251). 12/12 pass, fmt + clippy clean.

Luck-independent metric is composite MR rounds per candidate
((mr_rounds − ~50 prime rounds) / candidates); sieve survival ratios are
bit-size-independent, so 512-bit `-a 5` predicts 6300-bit exactly:

| Table | 512-bit tests/cand | 6300-bit composite rounds/cand | vs 251 |
|---|---|---|---|
| 251 | 268/1259 = 21.3% | 822/4087 = 0.201 | baseline |
| 1000 | 110/674 = 16.3% | 347/2138 = 0.162 | **−19%** |
| 10000 | 92/716 = 12.8% | ~0.107–0.121 (small sample) | **≈ −40%** |

Matches Mertens (survival ratio ln(251)/ln(L)) within noise at both sizes.
ms/MR round pinned at ~131–133 in all three 6300-bit verdicts (52.03 s /
397 rounds at 1k; 8.26 s / 63 rounds at 10k — that run drew only 122
candidates, ~5% lucky tail, but the per-candidate ratio is what counts), so
filtering 1228 primes is still **free** at 6300 bits. Default changed to
10000; `--all` 256-bit smoke green on all methods.

Reverted: limit 100000. MR-test ratio keeps falling (9.6%, Mertens ✓) but
ms/round jumps 0.15 → 0.39 at 512 bits — 9591 `%` ops per candidate now cost
more than the MR rounds they save. Negative result kept on record so nobody
retries it without a window sieve (bulk sieving is the only way to make
100K+ tables cheap — see backlog #2).

### 6. Reverted: `target-cpu=native` (Ryzen 5 3600X, 6c/12t, rustc 1.98.1)

Interleaved A/B, saved binaries (`prime-rust-generic.exe` vs
`-C target-cpu=native`), probablePrime on the 10k table:

| Size | Generic ms/round | Native ms/round | Verdict |
|---|---|---|---|
| 512 (`-a 5` ×2 pairs) | 0.175 | 0.20 | −14% (noise-tinged, all ~10 ms) |
| 2048 (`-a 5`) | 4.85 | 4.91 | tie |
| 6300 (`-a 1`) | 132.59 | 136.17 | −3%, no win |

Not adopted; shipped binary stays generic x86-64 (portable, no README
machine-specific caveat needed). Why fibonacci's 1.7× doesn't transfer: its
hot loops were crate-local u32/u64 NTT butterflies that autovectorize with
AVX2; this hot path is num-bigint's scalar carry-chain (`adc`) multiply/
divide, which has no vector form — rustc can't synthesize `adcx`/`adox`
from generic carry code. Lesson: **native helps vectorizable code, not
carry-chained bigint math.** Revisit only if a hand-written mulmod loop
(backlog #3) lands — that code *would* be vectorization-sensitive.

### 7. Kept: parallel candidates via rayon `--jobs` (4.6–6× on 6 cores)

`--jobs N` (default 1 = untouched serial baseline) tests batches of `8N`
candidates across a dedicated pool; first-found prime wins with honestly
counted `Work:` (bounded one-batch speculation). `next-prime` stays serial
with a log note — smallest-prime semantics can't parallelize without
changing meaning. 16/16 tests (incl. parallel validity + `accumulate` +
trial-verdict unit tests), fmt + clippy clean. rayon at exactly one level
(candidates; sequential internally-parallel stages, never nested).

6300-bit verdicts, 10k table (wall-ms/MR-round; serial CPU ≈ 135–140):

| Method | Serial | `--jobs 6` | `--jobs 12` |
|---|---|---|---|
| fastGcd | ~140 | **30.0 (4.6×)** | 27.3 (4.9×; SMT ≈ nothing for scalar code) |
| probablePrime | ~135 | composites **~6×**, overall tail-limited (see below) | — |

Three stacked diagnoses to get there, each measured (temporary `Instant`
guards + sleep-isolation test, all reverted; per-core `Get-Counter`
sampling caught a monitoring own-goal — samples landed after a 5 s lucky
run finished):
1. **Serial GCD fill Amdahl-capped at 2.4×.** Fix: parallel over-draw
   (8–12× batch, ~20%/13% pass rates) with surplus carried to next round.
2. **Join starvation at 2.8×.** With one item per thread there is nothing
   to steal; the ~40% of items that trial-reject idle while others grind.
   Fix: `PAR_BATCH_MULTIPLE = 8` oversubscription for work-stealing (3.8×).
3. **Trial-reject idleness at 3.8×.** Fix: fill-stage [`sieve_batch`]
   pre-filter (shared `trial_division` verdicts, byte-identical to serial)
   so test batches run ~100% MR-dense → 4.6× (77% efficiency; rest is join
   imbalance + the prime tail).
- Dead ends kept on record: **allocator contention** (J2 scaled 1.4× with
  and without mimalloc — disproven; dual-process control hit 2.0×, proving
  the box scales and the cap was in-process batch mechanics, not hardware);
  **thermal/CCX theories** (same control kills them).
- Serial tail, honestly Amdahl-bound: the prime's own confirmation rounds
  (20 at cert 40 ≈ 2.7 s, 50 at cert 100 ≈ 6.8 s) run single-threaded. A
  lucky 96-test probable run is tail-dominated (8.4 s ≈ 6.8 s tail + fast
  composites); a typical 500-round run projects (450/6+50)×135 ms ≈ 17 s
  vs ~67 s serial ≈ 4×. Round-level parallelism of the tail was considered
  and rejected (second parallel pattern + speculative waste for one tail
  per run). The certainty price stays.
- mimalloc: added for the disproven contention hypothesis; kept (builds
  clean, serial shows no regression — 512-bit 0.13 ms/round post-refactor
  vs 0.15–0.18 before) but its isolated delta is unmeasured tech debt.
- Correction to §4: trial division is **~1.3 ms/candidate** at 6300/10k
  (microbenched, both `%` forms identical — num-bigint already takes the
  single-limb path), i.e. ~1% of a round, not "zero measurable ms". The
  old claim rested on the circular `wall/rounds` metric, which bakes all
  non-MR time into ms/round by construction. Refactor bonus: the extracted
  helper uses `% u32` + skips the equality allocation for large candidates.

### 8. Kept: windowed bulk sieve for sequential scans (mulmod investigated and dropped)

Mulmod first became mulmod never: num-bigint 0.4.8's `modpow` already runs
Montgomery multiplication for odd moduli (`biguint/power.rs`, `monty.rs`),
and every candidate here is odd. The only non-Montgomery divisions left are
the MR tail's ~1–2 squarings per round versus ~9000 Montgomery muls inside
`modpow` — a Barrett rewrite would chase ~0.003% of the round. Backlog #2 is
closed as investigated; the vendored-source reading paid off by killing the
item instead of building it.

What shipped instead is the rescoped #4: `sieve_window` sieves a contiguous
run of `WINDOW_ODDS = 4096` odds with one `% p` remainder per prime plus
index arithmetic (`i ≡ (−rem)·inv2 mod p`), ~4000x fewer divisions than
per-candidate trial division, with slot-for-slot identical verdicts (new test
`sieve_window_matches_trial_division`, plus `< 2` and table-prime-equality
edge handling so tiny inputs stay exact). Sequential scanners
(`method_probable_prime`, `next_probable_prime_from`) test survivors in scan
order — same first-prime semantics, same RNG consumption, same `Work:`
counting. The parallel probable fill uses it serially (replacing
`sieve_batch`'s ~471k parallel mods per scan with 1228 serial mods), and
`test_batch` now goes straight to `miller_rabin` (no re-sieving; `small`
param dropped). `method_fast_gcd` paths are untouched — random draws can't
window. Default stays 10000: a 100k table would regress fastGcd's
per-candidate trial; it remains a `--sieve-limit` option for sequential runs.

Measured at 512 bits (`-a 3`, serial): probablePrime round cost 0.15 →
0.11 ms/round. 17/17 tests, fmt + clippy clean.

### 9. Final 6300-bit verdicts (this session, same box, n = 1 — luck-labeled)

Java harness gained a method filter (`probable|next|fastgcd`, no args =
original all-three behavior); compiled with `javac 21` to TEMP, `out/`
artifacts untouched. The earlier fastGcd DNF did not reproduce today
(42.70 s) — geometric variance, not a fix; the source change is
harness-only.

| Method | Java | Rust serial | Rust `--jobs 6` |
|---|---|---|---|
| probablePrime (cert 100) | 27.85 s (18.85 s earlier same-day run) | 22.64 s / 21.27 s (167/160 rounds) | 20.58 s (529 rounds, unlucky) |
| nextPrime (cert 100) | 65.32 s | 32.95 s (255 rounds) | serial by design |
| fastGcd (cert 40) | 42.70 s | 7.29 s (55 rounds, lucky) | 7.10 s (163 rounds) |

Luck-independent unit costs: Rust serial ≈ 129–136 ms/MR round across all
methods; parallel effective ≈ 39–44 ms/round → ~3.1–3.5x on 6 cores for
these samples. Wall-clock n = 1 gaps under ~2x are luck, not signal (see
§4); the durable claims are (a) unit-cost parity-or-better vs Java's
Montgomery+Karatsuba BigInteger on the same round budgets, and (b) ~3.5x
parallel efficiency with honestly counted work.

### 10. This session: profiled, sieve 100k kept, PGO and bigint-swap reverted

samply can't run here (needs Windows admin for ETW), so phases were timed
with a temporary release probe (deleted after). Per 6300-bit candidate: RNG
0.001 ms, trial division ~0.003 ms average (~1 ms only for full-scan
survivors — early exit after ~1.5 primes), GCD filter 0.096 ms, MR base
sampling 0.001 ms, **one Montgomery modpow 130.2 ms (99.9%)**. Everything
is the modular exponentiation. RNG (backlog #6) closed as negligible by the
same probe — no change.

- KEPT: default sieve 10000 → 100000. Filter power confirmed at 512 bits
  (`-a 5` serial probable: 79 → 61 MR tests per 20480 candidates, ratio
  0.77 = Mertens exactly); window-sieve cost at 6300/100k ≈ 0 (673-round
  run = 673 × 131.8 ms MR + nothing measurable). Expect ~20% fewer MR
  tests at 6300 for ~1 ms extra trial per survivor; fastGcd serial now
  shows the trial in its unit cost (148 vs 132 ms/round) but still nets
  ~20% on equal luck. 512-bit dev runs get ~5% slower in absolute sub-ms
  terms — accepted, verdict size is what matters.
- REVERTED: crypto-bigint 0.6 (`alloc`) `BoxedMontyForm` shootout, 6
  modpows at 6300 bits: num-bigint 129.5 ms vs crypto-bigint 135.5 ms
  (constant-time code costs more than num-bigint's vartime paths; results
  verified identical). No bigint swap.
- REVERTED: PGO (15 profraws from 512/2048/6300 training, llvm-profdata
  merge, profile-use rebuild): 512-bit 7.70 → 7.94 ms, 6300-bit 131.8 →
  129.1 ms/round — inside noise both times. Tight Montgomery loops leave
  PGO nothing to do; plain build restored so `target/` matches README.
- Stopping point: the wall is the bigint ALU now (assembly Montgomery /
  GMP would be the next lever, not more Rust-level tuning). 1M tables
  deliberately not pursued: marginal MR savings versus real fastGcd-path
  trial regressions.

Final 6300-bit table — plain build, 100k sieve, n = 1, luck-labeled (Java
numbers from the prior session, same box):

| Method | Java | Rust serial | Rust `--jobs 6` |
|---|---|---|---|
| probablePrime (cert 100) | 27.85 s | 37.38 s (283 rounds) | 23.61 s (577 rounds) |
| nextPrime (cert 100) | 65.32 s | 19.90 s (150 rounds) | serial by design |
| fastGcd (cert 40) | 42.70 s | 37.08 s (250 rounds) | 3.61 s (67 rounds) |

Luck-independent: serial ≈ 129–133 ms/MR round on sequential paths
(fastGcd serial 148, trial included); parallel effective ≈ 41–54 ms/round
→ ~3x on 6 cores for these samples. Same reading as §9: wall n = 1 gaps
under ~2x are luck; durable claims are unit-cost parity-or-better and
parallel efficiency.

### 11. Const-limb Montgomery shootout (reverted); native A/B skipped

Follow-up to §10's boxed shootout: `MontyForm<99>` (stack limbs, runtime
modulus — `ConstMontyForm` needs a compile-time modulus, unusable for
random candidates), 6 modpows at 6300 bits: num-bigint 129.6 ms vs
const-monty 146.3 ms (identical results). Constant-time
double-and-add-always loses to num-bigint's windowed vartime
exponentiation, and no `pow_vartime` exists in this API. Bigint-swap avenue
closed within ideals — boxed AND const both slower. Dep removed.

native vs generic A/B deliberately NOT re-run: no vectorizable loop landed
since §6 (window sieve is scalar strided stores; the hot path is
byte-identical library carry chains), so §6's flat verdict transfers
exactly. Re-measuring what no change could have moved would be ritual, not
measurement. Native stays banked for hand-written intrinsic loops (§8/option
2), where it actually matters.

### 12. KEPT: hand-written Montgomery exponentiation (1.18x on rounds)

`src/monty.rs`: SOS → CIOS → fixed-100 (`[u64; 100]`, stack, unrolled) →
dedicated triangle squarer → sliding-window (w=4) `MontyCtx100`, wired
into `miller_rabin` behind a measured size threshold. Per-candidate ctx
(R² + `mont(1)` precomputed once); bases enter per round; tail squares
and comparisons stay in-domain — `exit()` never runs in the loop.
Same round counts, same search order, same `Work:` counters; even,
oversize, and small candidates fall back to num-bigint transparently.

- Per-mulmod at 100 limbs: SOS 21.1µs, CIOS-slice 20.3µs, fixed mul
  17.0µs, dedicated square 13.8µs, `(a*b)%m` 20.4µs. Full modpow:
  plain square-multiply 140.6 ms (SLOWER — 9450 mulmods vs ~7600 via
  windowing), windowed **110.2 ms vs 129.7 ms = 1.18x**.
- Threshold `MONTY_MIN_BITS = 6144` (per-round ratios new/old: 4608 →
  1.76x loss … 5888 → 1.03x, 6144 → 0.87x, 6300 → 0.81x). Fixed-100 work
  scales with exponent bits (linear) vs num-bigint superlinear — the win
  grows with size.
- Rejected along the way: crypto-bigint boxed (135.5) and const-limb
  (146.3) — constant-time double-and-add loses to windowed vartime;
  PGO (noise); native A/B re-run (no vectorizable loop landed).
  Intrinsics (`mulx`/`adcx`) deferred: safe code keeps winning on
  structure, and `unsafe` spends only when measurement demands it.
- Verification: 400+ differential cases vs the BigUint oracle (12 widths,
  window-edge exponents, all-ones crucibles to 6400 bits), 6300-bit
  release differentials, full gate green (24+1 tests). Temp probes and
  benches deleted; SOS/slice kernels `#[cfg(test)]`-gated as references.

6300-bit verdicts, plain build, 100k sieve, n = 1, luck-labeled (Java
numbers from §9, same box):

| Method | Java | Rust serial (monty) | Rust `--jobs 6` (monty) |
|---|---|---|---|
| probablePrime (cert 100) | 27.85 s | 8.58 s (77 rounds) | 21.02 s (625 rounds) |
| nextPrime (cert 100) | 65.32 s | 7.16 s (65 rounds) | serial by design |
| fastGcd (cert 40) | 42.70 s | 48.76 s (396 rounds) | 12.18 s (403 rounds) |

Unit costs: sequential paths 132–136 → **110–111 ms/MR round**;
fastGcd serial 148 → 123 (100k-trial showing); parallel effective
30–34 ms/round (≈ 3.3x efficiency, unchanged). Same reading as ever:
sub-2x wall gaps are luck; unit cost and parallel efficiency are signal.

### 13. Cross-size sweep (512–6300 bits, Rust only — Java harness is 6300-fixed)

Serial probablePrime unit cost (ms/MR round, `-a 5` ≤ 2048, `-a 3` at
4096, `-a 2` ≥ 6144 — every run a fresh random prime):

| Bits | 512 | 1024 | 2048 | 4096 | 6144 | 6300 |
|---|---|---|---|---|---|---|
| ms/round | 0.12 | 0.65 | 4.6 | 34.4 | 108.6 | 111.4 |

Superlinear below the 6144 threshold (old path); the 6144 cell lands
~1.12x under the old-path projection — the crossover probe's prediction,
confirmed end-to-end. Two findings: (a) serial fastGcd is the SLOWEST
method at ≤ 4096 bits (100k per-candidate trial scans dominate: 0.60 /
1.81 / 7.61 / 41.5 ms/round vs probable's 0.12–34.4) — the cert-40
advantage only pays at 6144+; (b) parallelism loses-or-ties at ≤ 1024
(batch overhead dwarfs ms-scale work), breaks even ~2048, and wins
clearly ≥ 4096 (fastGcd 6.90 s → 1.38 s at 4096, probable 33.8 s →
7.4 s at 6144 on lucky shapes; tail-limited otherwise, as §7 models).

### 14. fastGcd batch mechanics: fused fill (neutral, kept as cleanup),
chunked test batches (kept, mechanism confirmed)

Fastest-method focus (§13: parallel fastGcd wins ≥ 4096 bits). Profile at
6300: ~90% MR tests, ~8% batch joins/allocs, ~2% GCD — batch overhead the
only remaining lever. Two surgical changes, measured one at a time:

- Fused draw+GCD+trial into ONE parallel pass (deleted `sieve_batch`,
  −40 lines, one fewer join + intermediate Vec). 2048 `-a 5`: 387 ms vs
  338 ms on IDENTICAL work counts (576 tests / 690 rounds both) — pure
  box noise, CIs overlapping 5x over; 6300 `-a 3`: 7.11 s, 33.7 ms/round
  mid-band. No detectable difference → kept as a simplification with
  identical semantics/stats, documented as measurement-neutral.
- `test_batch` → two sequential half-joins (`test_chunk`), early return
  when the first half holds the prime (result identical: a first-half
  prime always wins batch order; skipped work correctly uncounted).
  ms/round is BLIND to this by construction (skipped tests leave numerator
  and denominator together), so it was verified by instrumented skip
  counting instead: 2 skips × 24 tests (≈ 0.44 s each) across 3 runs ≈
  4% — the predicted 2–5%, then instrumentation removed. Near-zero cost
  (one extra ~µs join per batch), kept.

Deliberately not pursued: GCD widening (§5 stands), tail parallelism
(§7 stands), sieve beyond 100k (fastGcd per-candidate trial regresses),
pipelined sieve+test (load-imbalance complexity ≫ ~7%).

## Optimization backlog (in planned order)

1. ~~Wider small-prime table~~ KEPT at 10000 (see 5, 8); 100000 REVERTED
   (would regress fastGcd's per-candidate trial).
2. ~~Faster MR rounds (mulmod)~~ DROPPED after source reading: `modpow` is
   already Montgomery for odd moduli (see 8).
3. ~~Parallel candidate search~~ KEPT: `--jobs`, 4.6–6× (see 7).
4. ~~Window sieve~~ KEPT for sequential scans (see 8); fastGcd random draws
   can't window.
5. **Cheap wins:** ~~native (reverted, §6)~~, mimalloc kept-unisolated (§7),
   ~~PGO (reverted, §10 — no win on Montgomery loops)~~.
6. ~~RNG contribution check~~ CLOSED negligible by probe (§10, ~0.002 ms vs
   130 ms rounds).

## Benchmark hygiene policy (this repo)

- `cargo build --release`, dev at ≤ 512 bits, verdicts at 6300 bits, rare.
- `-a 3+` with CI/stddev/skewness for close calls; single samples only for
  direction with ≥ 2× gaps, labeled as such.
- Interleaved A/B binaries for close calls; redirect to file (`> out.txt
  2>&1`), never bare pipes; `Get-Process` kill-hygiene check after any kill.
- Gate before every commit: `cargo fmt --check` + `cargo clippy --all-targets`
  (zero warnings) + `cargo test`.
