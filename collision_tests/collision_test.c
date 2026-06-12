/*
 * collision_test.c  --  Sponge-hash collision / output-distribution battery
 *                       for Krakken-2048 (krakken_hash_avx2).
 *
 * Probes hash-level properties not covered by the permutation SAC battery or the
 * NIST keystream tests:
 *
 *   1. Birthday collisions : hash N random messages, truncate each digest to t
 *      bits, and count colliding pairs. For an ideal random function the
 *      expected number of colliding pairs is C(N,2)/2^t ~= N^2 / 2^(t+1).
 *      A structural collision weakness shows up as far more collisions than
 *      this at some truncation width.
 *   2. Output diffusion    : for random 1-bit input differences, the full
 *      256-bit digest difference Hamming weight (ideal mean 128). The minimum
 *      over many trials staying far from 0 rules out cheap near-collisions.
 *   3. Digest uniformity   : byte-value frequency over many digests, chi-square
 *      with 255 degrees of freedom (ideal ~= 255 +- sqrt(510)).
 *
 * Build:  make           (gcc -O3 -mavx2 ... -lpthread -lm)
 * Run:    ./collision_test [N]      (default N = 1000000)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void init_rc_vectors_avx2(void);
void krakken_hash_avx2(uint8_t *out, size_t outlen, const uint8_t *in,
                       size_t inlen);

static inline uint64_t splitmix64(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* 24-byte random message -> 64-bit digest truncation, for the birthday test. */
static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static void test_birthday(uint64_t N) {
    printf("--- Test 1: birthday collisions (N = %llu messages) ---\n",
           (unsigned long long)N);
    uint64_t *dig = malloc(N * sizeof(uint64_t));
    if (!dig) { fprintf(stderr, "alloc failed\n"); return; }

    uint64_t rng = 0x1234567890ABCDEFULL;
    for (uint64_t i = 0; i < N; i++) {
        uint8_t msg[24];
        for (int k = 0; k < 3; k++) {
            uint64_t r = splitmix64(&rng);
            memcpy(msg + k * 8, &r, 8);
        }
        uint8_t out[8];
        krakken_hash_avx2(out, 8, msg, sizeof(msg));
        memcpy(&dig[i], out, 8);
    }

    printf("    %-6s %-14s %-14s %-8s\n", "t", "expected", "observed", "obs/exp");
    for (int t = 28; t <= 42; t += 2) {
        uint64_t mask = (t >= 64) ? ~0ULL : ((1ULL << t) - 1);
        uint64_t *tr = malloc(N * sizeof(uint64_t));
        for (uint64_t i = 0; i < N; i++) tr[i] = dig[i] & mask;
        qsort(tr, N, sizeof(uint64_t), cmp_u64);

        /* colliding pairs = sum over equal-value runs of C(run,2) */
        uint64_t pairs = 0, run = 1;
        for (uint64_t i = 1; i < N; i++) {
            if (tr[i] == tr[i - 1]) run++;
            else { pairs += run * (run - 1) / 2; run = 1; }
        }
        pairs += run * (run - 1) / 2;
        free(tr);

        double expect = (double)N * (double)(N - 1) / 2.0 / ldexp(1.0, t);
        printf("    %-6d %-14.3f %-14llu %-8.3f\n",
               t, expect, (unsigned long long)pairs,
               expect > 0 ? pairs / expect : 0.0);
    }
    free(dig);
    printf("\n");
}

static void test_diffusion(uint64_t N) {
    printf("--- Test 2: 1-bit-input digest diffusion (N = %llu pairs, 256-bit out) ---\n",
           (unsigned long long)N);
    uint64_t rng = 0xFEEDFACECAFEBEEFULL;
    long sum = 0, minw = 256, maxw = 0;
    double sumsq = 0;
    for (uint64_t i = 0; i < N; i++) {
        uint8_t msg[32], msg2[32], d1[32], d2[32];
        for (int k = 0; k < 4; k++) {
            uint64_t r = splitmix64(&rng);
            memcpy(msg + k * 8, &r, 8);
        }
        memcpy(msg2, msg, 32);
        int bit = (int)(splitmix64(&rng) % 256);
        msg2[bit >> 3] ^= (uint8_t)(1u << (bit & 7));

        krakken_hash_avx2(d1, 32, msg, 32);
        krakken_hash_avx2(d2, 32, msg2, 32);

        int w = 0;
        for (int b = 0; b < 32; b++) w += __builtin_popcount((unsigned)(d1[b] ^ d2[b]));
        sum += w; sumsq += (double)w * w;
        if (w < minw) minw = w;
        if (w > maxw) maxw = w;
    }
    double mean = (double)sum / N;
    double var = sumsq / N - mean * mean;
    printf("    mean = %.3f (ideal 128.0)   std = %.3f (ideal %.3f)   min = %ld   max = %ld\n\n",
           mean, sqrt(var), sqrt(256 * 0.25), minw, maxw);
}

static void test_uniformity(uint64_t N) {
    printf("--- Test 3: digest byte uniformity, chi-square (N = %llu digests) ---\n",
           (unsigned long long)N);
    uint64_t freq[256] = {0};
    uint64_t rng = 0x0F0F0F0F0F0F0F0FULL;
    for (uint64_t i = 0; i < N; i++) {
        uint8_t msg[16], out[32];
        for (int k = 0; k < 2; k++) {
            uint64_t r = splitmix64(&rng);
            memcpy(msg + k * 8, &r, 8);
        }
        krakken_hash_avx2(out, 32, msg, 16);
        for (int b = 0; b < 32; b++) freq[out[b]]++;
    }
    uint64_t total = N * 32;
    double expect = (double)total / 256.0;
    double chi2 = 0;
    for (int v = 0; v < 256; v++) {
        double d = (double)freq[v] - expect;
        chi2 += d * d / expect;
    }
    /* chi-square, 255 dof: mean 255, std sqrt(2*255)=22.58 */
    double z = (chi2 - 255.0) / sqrt(2.0 * 255.0);
    printf("    total bytes = %llu   chi2 = %.2f  (255 dof, ideal 255 +- 22.6;  z = %+.2f)\n\n",
           (unsigned long long)total, chi2, z);
}

int main(int argc, char **argv) {
    init_rc_vectors_avx2();
    uint64_t N = 1000000;
    if (argc >= 2) N = strtoull(argv[1], NULL, 10);

    printf("Krakken-2048 sponge-hash collision / distribution battery\n");
    printf("=========================================================\n\n");

    test_birthday(N);
    test_diffusion(N < 200000 ? N : 200000);
    test_uniformity(N < 200000 ? N : 200000);
    return 0;
}
