# Session handoff — PrimeRust optimization loop

## Goal for next session
`docs.md` backlog reordered by measured value: **#2 mulmod first** (rounds
are the parallel wall AND the serial prime tail — 1.5× cheaper rounds ≈
1.5× everything; read num-bigint vendored sources first, then integer
Barrett `mulmod`), then #4 rescoped window sieve (~1–1.5 s prize), #5 PGO,
#6 RNG check. Rayon KEPT, native REVERTED — see `docs.md` §6–7.

## Current code state (on disk, gate green: fmt 0, clippy zero warnings, 16/16 tests)
- `src/prime.rs` — prior work (`SearchStats`, `primes_up_to`, `small: &[u32]`)
  + new: `TrialVerdict`/`trial_division` helper (single-limb `% u32`, no
  equality alloc for large candidates; serial path byte-identical via match),
  `sieve_batch` fill-stage pre-filter, `test_batch` (one parallel level),
  `method_{probable_prime,fast_gcd}_parallel` (bulk fill + carry + serial
  top-up), `PAR_BATCH_MULTIPLE = 8`, `SearchStats::accumulate`. GCD product
  stays 53-prime. New tests: `trial_division_verdicts_are_exact` (lesson:
  the 104729×1299709 semiprime is 136117223861, not my misremembered
  136117223 — test caught it), `stats_accumulate…`, parallel validity ×2.
- `src/main.rs` — `--jobs N` (default 1 = serial; batches of 8N,
  first-found wins), dedicated pool, next-prime serial fallback + log note,
  mimalloc global allocator (kept, isolated delta unmeasured — tech debt).
- `Cargo.toml` — rayon 1, mimalloc 0.1. Shipped binary: generic x86-64.
- Run logs this session: `v-6300-fgcd-serial.txt` (61.0 s),
  `v-6300-prob-serial.txt` (11.2 s, lucky), `v-6300-prob-j6.txt` (18.1 s,
  pre-density), `v-6300-fgcd-j6.txt` (26.4 s, pre-density),
  `v2-6300-fgcd-j6.txt` (6.0 s lucky) + `-j2` (19.8 s),
  `x-6300-fgcd-j6.txt` (9.3 s, **4.6×**), `x-6300-prob-j6.txt` (8.4 s,
  tail-limited), `x-6300-fgcd-j12.txt` (5.8 s, 4.9× — SMT ≈ nothing),
  `w-6300-fgcd-j6.txt` (33.7 s pre-density), `ab-*` (native REVERTED:
  −14%/−1%/−3% at 512/2048/6300), `smoke-512-serial.txt` (no refactor
  regression), dual-proc control (2.0× — box scales, cap was in-process).

## Performance model (docs.md §4 + §7 corrections)
- Serial: time ≈ MR_rounds × ~135 ms + candidates × ~1.3 ms trial-div
  (~1% — the old "zero ms" claim was circular-metric artifact).
- Parallel (`--jobs 6`, dense batches): composites ~6×, overall 4.6×
  fastGcd; probable tail-limited by the prime's 50 rounds (6.8 s floor).
  Prime tail is the certainty price — round-parallelism rejected.
- Luck-independent metrics: composite rounds/candidate (tables),
  wall-ms/round same-window (parallelism). Geometric luck still dominates
  n=1 (this session: 252–9216 candidates).

## Backlog (docs.md has detail)
1. ~~Wider table~~ KEPT 10000. 2. **Mulmod (next).** 3. ~~Parallel~~ KEPT.
4. Window sieve rescoped (~1 s prize). 5. PGO (mimalloc unisolated).
6. RNG check (<1% expected).

## Hygiene (from learned.md, bitten twice already)
- Redirect to file, never bare pipes. `Get-Process java/prime-rust` after any kill.
- Gate before every commit: fmt + clippy (zero warnings) + test.
- Dev ≤ 512 bits; 6300-bit verdicts rare, `-a 3+` + interleave for close calls.
- No stray processes (verified). TEMP-DIAG files removed (`tests/diag_*.rs`,
  `diag*.txt`, `mon-*`, `proc-*`). All `eprintln! DIAG` reverted.
