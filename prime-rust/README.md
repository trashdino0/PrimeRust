# prime-rust

6300-bit prime generation benchmark — exact Rust port of
`../prime numbers/src/PrimeBenchmark.java` (untouched original).

Baseline result: **parity with Java** on `probablePrime` at 6300 bits
(Rust 10.68 s vs Java 10.36 s, n = 1, noisy box). See [docs.md](docs.md)
for the full decision log, including two findings: Java's "Fast GCD Filter"
never finishes at 6300 bits (DNF > 560 s; same algorithm in Rust: 59.85 s),
and a stalled Java re-run logged as an anomaly.

## Build, test, bench

```powershell
cargo fmt --check; cargo clippy --all-targets; cargo test   # gate (must be clean)
cargo build --release
.\target\release\prime-rust.exe --bits 256 -a 2 --all       # quick check (~ms)
.\target\release\prime-rust.exe --bits 6300 -a 1 --all      # verdict run (~min)
.\target\release\prime-rust.exe --help                      # -m probable-prime|next-prime|fast-gcd
```

Java baseline (do not modify the original project):

```powershell
java -jar "..\prime numbers\out\artifacts\prime_numbers_jar\prime.jar" > java-run.txt 2>&1
```

Always redirect to a file (killed-pipe output is lost) and check
`Get-Process java` after any killed run. Heavy 6300-bit runs are rare by
policy — dev at ≤ 512 bits.

## Layout

- `src/prime.rs` — the three ported methods + Miller-Rabin core + table-driven tests
- `src/main.rs` — clap CLI, `-a N` mean / 95% CI / stddev / skewness, SHA256 beside every number
- `docs.md` — decision log + optimization backlog (the repeat loop starts there)
