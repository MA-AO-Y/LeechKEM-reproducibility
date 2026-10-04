#include "bench_random.h"

static uint64_t state = UINT64_C(0x9e3779b97f4a7c15);

void bench_random_reset(uint64_t seed)
{
    state = seed != 0 ? seed : UINT64_C(0x9e3779b97f4a7c15);
}

static uint64_t splitmix64(void)
{
    uint64_t z = (state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

int randombytes(unsigned char *out, unsigned long long outlen)
{
    unsigned long long i = 0;
    while (i < outlen) {
        uint64_t word = splitmix64();
        unsigned j;
        for (j = 0; j < 8 && i < outlen; ++j, ++i) {
            out[i] = (unsigned char)(word >> (8U * j));
        }
    }
    return 0;
}
