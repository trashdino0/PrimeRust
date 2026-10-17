# monty — hand-written Montgomery kernel project log

Goal: replace num-bigint's Montgomery `modpow` on the MR hot path with our
own CIOS Montgomery kernel (`_addcarry_u64` / `mulx`), keeping the same
Miller-Rabin algorithm, round counts, search order, and `Work:` counters —
the ideals: same work units, our code doing the work.

Baseline (this box, Ryzen 5 3600X): one 6300-bit modpow = 130.2 ms (99.9%
of a candidate). Challengers already beaten: crypto-bigint boxed (135.5),
const-limb (146.3). Target: 1.5–3x on rounds.

## Constraints

- Correctness first: every kernel lands with differential tests vs the
  num-bigint oracle (learned.md §4 style) across sizes incl. edges.
- `unsafe` only inside the limb kernels, justified per site; no `unsafe`
  in search/sieve/RNG logic. Miri can't run x86 intrinsics — differential
  tests + debug asserts are the safety net (recorded, not forgotten).
- Portable binary: fast path via `#[target_feature(bmi2,adx)]` + runtime
  `is_x86_feature_detected` dispatch, generic fallback otherwise. No global
  `target-cpu=native` machine-lock (user approved native, but dispatch keeps
  the shipped binary portable AND fast here).
- u64 limbs; byte-based BigUint conversions at the boundary (O(n) memcpy
  vs O(n²) kernel — negligible, robust to num-bigint internals).

## Phases

- [x] Phase 0: scaffold (`src/monty.rs`): `n0inv`, limb helpers, tests. Gate green.
- [x] Phase 1: SOS + CIOS safe kernels, 360+ differential checks (incl.
      100-limb all-ones crucible), microbench. Gate green (22/22).
- [x] Phase 2: `monty_modpow_100` (sliding w=4 + dedicated square) = 1.18x,
      differentials green incl. ignored 6300-bit release cases.
- [x] Phase 3: integrated (per-candidate ctx, domain tail, threshold
      6144 measured). Gate green (24+1). Temps deleted; SOS/slice
      `#[cfg(test)]`-gated. Verdicts: sequential unit 110–111 ms/round.
- [x] Intrinsics CIOS (`mulx`/`adcx`, dispatch): REVERTED. Correct after
      a carry-chain bug hunt (c_in=0, caught by a 10k-case u128-oracle
      regression over full-range limbs), then measured 23.3µs vs safe
      17.3µs — SLOWER. Root cause, read in stdarch source: `_mulx_u64` /
      `_addcarry_u64` are safe u128-wrapper const fns, not raw
      instructions; the "intrinsics" version was the same ops plus extra
      scalar carry bookkeeping and stack-spilled locals. True MULX/ADCX
      would need `asm!()` — parked as not-worth-it at this prize size.
      Kernel + all adx/temp tests deleted; gate green (24+1, 2.4 s).

## Measurements

- 6300-bit modpow (num-bigint): 130.2 ms (probe, seeded, n=10).
- Per-mulmod at 100 limbs, release: SOS-safe 21.1µs, CIOS-slice 20.3µs,
  fixed-100 mul 17.0µs, **dedicated square 13.8µs** (0.81x of mul —
  reduction is shared so the triangle only saves the product half),
  `(a*b)%m` 20.4µs. Blended exp ≈ 14.9µs avg vs 17.0 → ~1.14x rounds:
  square+CIOS ≈ ties num-bigint's ~14µs internal.
- Integrated verdicts (6300-bit, 100k, n=1): serial probable 8.58 s
  (77r), serial next 7.16 s (65r), serial fastGcd 48.76 s (396r),
  parallel probable 21.02 s (625r, 33.6 ms/r eff.), parallel fastGcd
  12.18 s (403r, 30.2 ms/r eff.). Sequential unit cost 110–111 ms/round.
- Phase 2 LANDED: windowed (w=4) `MontyCtx100::modpow` = **110.2 ms vs
  129.7 ms = 1.18x**. Diagnosis was algorithmic, not instruction-level:
  plain square-multiply measured 140.6 ms (SLOWER — 9450 mulmods vs
  num-bigint's ~7600 via its own windowing); w=4 table (1 square + 7
  muls precomp, ~1300 muls/round) closed it. Small-size differentials +
  fixed exponents (0/1/2/3/15/16/255/256) in debug, 3/3 random 6300-bit
  cases in release (`#[ignore]`d, debug modpow is seconds/case).
  Gate green (26+1).

## Continuity notes

- Phase 0 already paid for itself: the roundtrip test caught a real
  high/low-chunk swap in `to_limbs_padded` (would have poisoned every
  kernel silently). Oracle-style tests stay mandatory per kernel.
- Fixed-100 debugging war story (read before touching the finale): the
  kernel was right three times while asserts/tests were wrong — (1) a
  missing-`i` mirage (caught before compiling: CIOS is a rolling window,
  there IS no `i` offset), (2) an over-strong assert on unselected
  wrapping garbage (release oracle passed throughout; assert fixed to
  check the SELECTED top), (3) an out-of-contract wide-ones edge
  (spare assert correctly refused it; edge re-scoped to 2^6399-1).
  Lesson: when math and experiment disagree, instrument (the TOPNZ print
  settled it in one run) — and debug/release divergence means the bug is
  in a `debug_assert`, not the code.

## Decisions

- CIOS over SOS/FIOS: single loop, fewest passes over `m`, best register
  behavior for 99-limb inputs. SOS kept as the readable reference in tests
  if CIOS ever needs cross-checking (not yet needed — num-bigint is oracle).
- Square dedicated (`square` ≠ `mul(a,a)`): monty squaring saves ~half the
  limb-muls in the exponentiation loop.
- Comparisons stay in Montgomery domain in the MR tail (precompute
  mont(1), mont(n-1) per candidate); single `retrieve()` only for the final
  cross-check in tests, never per round.
