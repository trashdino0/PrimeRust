import java.math.BigInteger;
import java.security.SecureRandom;

public class PrimeBenchmark {

    static final SecureRandom random = new SecureRandom();

    static final int BIT_LENGTH = 6300; // ~1900 digits
    static final int WARMUP_ROUNDS = 0; // Let the JVM JIT compiler warm up
    static final int TEST_ROUNDS = 1;  // Run enough times for a true average

    // Product of the first 53 odd primes (up to 251)
    static final BigInteger SMALL_PRIME_PRODUCT;
    static {
        int[] primes = {
                3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
                151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
                229, 233, 239, 241, 251
        };
        BigInteger prod = BigInteger.ONE;
        for (int p : primes) {
            prod = prod.multiply(BigInteger.valueOf(p));
        }
        SMALL_PRIME_PRODUCT = prod;
    }

    public static void main(String[] args) {
        System.out.println("Starting warm-up phase (ignoring these times)...");
        for (int i = 0; i < WARMUP_ROUNDS; i++) {
            methodProbablePrime();
        }
        System.out.println("Warm-up complete. Starting benchmarks...\n");

        // Optional filter so one method can be timed without the others
        // (the Fast GCD path is extremely slow at 6300 bits and would
        // otherwise block the whole run): no args = all three (original
        // behavior), or pass probable|next|fastgcd.
        boolean all = args.length == 0;
        if (all || args[0].equalsIgnoreCase("probable")) {
            benchmark("BigInteger.probablePrime()", PrimeBenchmark::methodProbablePrime);
        }
        if (all || args[0].equalsIgnoreCase("next")) {
            benchmark("nextProbablePrime()", PrimeBenchmark::methodNextPrime);
        }
        if (all || args[0].equalsIgnoreCase("fastgcd")) {
            benchmark("Fast GCD Filter + Miller-Rabin", PrimeBenchmark::methodFastGcd);
        }
    }

    interface PrimeMethod {
        BigInteger run();
    }

    static void benchmark(String name, PrimeMethod method) {
        long totalMs = 0;
        BigInteger prime = null;

        for (int i = 0; i < TEST_ROUNDS; i++) {
            long start = System.nanoTime();
            prime = method.run();
            long end = System.nanoTime();

            totalMs += (end - start) / 1_000_000;
        }

        double avgMs = (double) totalMs / TEST_ROUNDS;

        System.out.println("Method: " + name);
        System.out.println("Average Time: " + avgMs + " ms");
        System.out.println("Digits: " + prime.toString().length());
        System.out.println("--------------------------------------------------");
    }

    // 1. The Gold Standard
    static BigInteger methodProbablePrime() {
        return BigInteger.probablePrime(BIT_LENGTH, random);
    }

    // 2. The Linear Scanner (Fastest, but slightly biased mathematically)
    static BigInteger methodNextPrime() {
        BigInteger start = new BigInteger(BIT_LENGTH, random)
                .setBit(BIT_LENGTH - 1);
        return start.nextProbablePrime();
    }

    // 3. The Custom Optimized Sieve
    static BigInteger methodFastGcd() {
        while (true) {
            // Generate odd number with exact bit length
            BigInteger candidate = new BigInteger(BIT_LENGTH, random)
                    .setBit(BIT_LENGTH - 1)
                    .setBit(0);

            // Fast filter: If GCD isn't 1, it shares a factor with our small primes
            if (!candidate.gcd(SMALL_PRIME_PRODUCT).equals(BigInteger.ONE)) {
                continue;
            }

            // 40 rounds of Miller-Rabin is industry standard for crypto-safety.
            // 100 rounds is unnecessary for numbers > 1000 bits.
            if (candidate.isProbablePrime(40)) {
                return candidate;
            }
        }
    }
}