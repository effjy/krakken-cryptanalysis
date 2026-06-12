/*
 * sac_test.c  --  Strict Avalanche Criterion (SAC) matrix, completeness,
 *                 and per-output-bit frequency (monobit) for Krakken-2048.
 *
 * These probe the per-bit dependency structure of the permutation, which the
 * mean-avalanche test in krakken_experiments.tex cannot see: a perfect mean of
 * 0.5 is compatible with a badly structured SAC matrix (some input bit always
 * flips a given output bit, or never does). Here we measure, for every
 * (input bit i, output bit j) pair:
 *
 *   SAC[i][j] = Pr[ output bit j flips | input bit i is flipped ]   (ideal 1/2)
 *
 * and derive:
 *   (1) SAC matrix : max |SAC - 1/2| and count of cells beyond 5 sigma.
 *   (2) Completeness: every (i,j) cell strictly between 0 and N
 *       (input i provably influences output j, and not deterministically).
 *   (3) Monobit  : per-output-bit frequency of 1 over random inputs (ideal 1/2).
 *
 * Build:
 *   gcc -O3 -mavx2 sac_test.c krakken_multi.c -o sac_test -lpthread
 * Run:
 *   ./sac_test [rounds] [N]      (default rounds=8, N=4096)
 *   ./sac_test all   [N]         (scan rounds 1..8)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void init_rc_vectors_avx2(void);
void krakken_permute_avx2_rounds(uint64_t state[32], int rounds);

#define NWORDS 32
#define NBITS  2048

static inline uint64_t splitmix64(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* count[i*NBITS + j] : how many of N trials had output bit j flipped when
 * input bit i was toggled. uint32 counters. */
static uint32_t *count;
/* freq[j] : how many random base states had output bit j == 1. */
static uint64_t freq[NBITS];

static void run_sac(int rounds, uint64_t N) {
    memset(count, 0, (size_t)NBITS * NBITS * sizeof(uint32_t));
    memset(freq, 0, sizeof(freq));

    uint64_t seed = 0xC0FFEE000000ULL ^ (uint64_t)rounds;
    uint64_t rng = seed;

    for (int i = 0; i < NBITS; i++) {
        uint32_t *crow = count + (size_t)i * NBITS;
        int iw = i >> 6, ib = i & 63;
        for (uint64_t s = 0; s < N; s++) {
            uint64_t X[NWORDS], Y[NWORDS];
            for (int w = 0; w < NWORDS; w++) X[w] = splitmix64(&rng);
            memcpy(Y, X, sizeof(X));
            Y[iw] ^= (1ULL << ib);

            krakken_permute_avx2_rounds(X, rounds);
            krakken_permute_avx2_rounds(Y, rounds);

            /* frequency: accumulate ones of P(X) (only on first input-bit pass
             * to avoid N*2048 oversampling -- N samples is plenty). */
            if (i == 0) {
                for (int w = 0; w < NWORDS; w++) {
                    uint64_t v = X[w];
                    while (v) {
                        int b = __builtin_ctzll(v);
                        freq[(w << 6) + b]++;
                        v &= v - 1;
                    }
                }
            }

            /* SAC: iterate set bits of the output difference. */
            for (int w = 0; w < NWORDS; w++) {
                uint64_t d = X[w] ^ Y[w];
                while (d) {
                    int b = __builtin_ctzll(d);
                    crow[(w << 6) + b]++;
                    d &= d - 1;
                }
            }
        }
    }

    /* ---- analysis ---- */
    double invN = 1.0 / (double)N;
    double sigma = 0.5 / sqrt((double)N);
    double thr5 = 5.0 * sigma;

    double sac_max_dev = 0.0, sac_sum_dev = 0.0;
    int sac_imax = -1, sac_jmax = -1;
    uint64_t cells_beyond_5sigma = 0;

    uint64_t dead = 0;   /* cell never flipped (==0)  -> incompleteness */
    uint64_t stuck = 0;  /* cell always flipped (==N) -> deterministic link */

    for (int i = 0; i < NBITS; i++) {
        const uint32_t *crow = count + (size_t)i * NBITS;
        for (int j = 0; j < NBITS; j++) {
            uint32_t c = crow[j];
            if (c == 0) dead++;
            else if (c == N) stuck++;
            double p = (double)c * invN;
            double dev = fabs(p - 0.5);
            sac_sum_dev += dev;
            if (dev > thr5) cells_beyond_5sigma++;
            if (dev > sac_max_dev) { sac_max_dev = dev; sac_imax = i; sac_jmax = j; }
        }
    }
    double sac_mean_dev = sac_sum_dev / ((double)NBITS * NBITS);

    /* monobit */
    double freq_max_dev = 0.0; int freq_jmax = -1;
    uint64_t freq_beyond_5sigma = 0;
    for (int j = 0; j < NBITS; j++) {
        double p = (double)freq[j] * invN;
        double dev = fabs(p - 0.5);
        if (dev > thr5) freq_beyond_5sigma++;
        if (dev > freq_max_dev) { freq_max_dev = dev; freq_jmax = j; }
    }

    double expect_fp = (double)NBITS * NBITS * 5.733e-7; /* 2-sided 5 sigma */

    printf("=== rounds=%d  N=%llu  (5 sigma = %.5f) ===\n",
           rounds, (unsigned long long)N, thr5);
    printf("  SAC matrix  : max|p-0.5| = %.5f at (in=%d,out=%d)   mean|p-0.5| = %.5f\n",
           sac_max_dev, sac_imax, sac_jmax, sac_mean_dev);
    printf("                cells > 5sigma = %llu / %d   (random-expected ~ %.1f)\n",
           (unsigned long long)cells_beyond_5sigma, NBITS * NBITS, expect_fp);
    printf("  Completeness: dead cells (never flip) = %llu   stuck cells (always flip) = %llu\n",
           (unsigned long long)dead, (unsigned long long)stuck);
    printf("  Monobit     : max|freq-0.5| = %.5f at out=%d   bits > 5sigma = %llu / %d\n\n",
           freq_max_dev, freq_jmax,
           (unsigned long long)freq_beyond_5sigma, NBITS);
}

int main(int argc, char **argv) {
    init_rc_vectors_avx2();
    count = malloc((size_t)NBITS * NBITS * sizeof(uint32_t));
    if (!count) { fprintf(stderr, "alloc failed\n"); return 1; }

    int scan_all = 0, rounds = 8;
    uint64_t N = 4096;
    if (argc >= 2) {
        if (strcmp(argv[1], "all") == 0) scan_all = 1;
        else rounds = atoi(argv[1]);
    }
    if (argc >= 3) N = strtoull(argv[2], NULL, 10);

    printf("Krakken-2048 SAC / completeness / monobit probe\n");
    printf("================================================\n\n");

    if (scan_all) for (int r = 1; r <= 8; r++) run_sac(r, N);
    else run_sac(rounds, N);

    free(count);
    return 0;
}
