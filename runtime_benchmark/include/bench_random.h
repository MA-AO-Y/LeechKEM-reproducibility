#ifndef BENCH_RANDOM_H
#define BENCH_RANDOM_H

#include <stddef.h>
#include <stdint.h>

/*
 * Deterministic benchmark-only generator.  It deliberately replaces external
 * entropy acquisition so KEM timings measure the cryptographic core.
 * It is NOT suitable for production cryptography.
 */
void bench_random_reset(uint64_t seed);
int randombytes(unsigned char *out, unsigned long long outlen);

#endif
