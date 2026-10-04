#define _POSIX_C_SOURCE 200809L

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bench_random.h"

#if FRODO_LEVEL == 640
#include "api_frodo640.h"
#define SCHEME_NAME "FrodoKEM-640-SHAKE"
#define KEM_KEYPAIR crypto_kem_keypair_Frodo640
#define KEM_ENCAPS crypto_kem_enc_Frodo640
#define KEM_DECAPS crypto_kem_dec_Frodo640
#elif FRODO_LEVEL == 976
#include "api_frodo976.h"
#define SCHEME_NAME "FrodoKEM-976-SHAKE"
#define KEM_KEYPAIR crypto_kem_keypair_Frodo976
#define KEM_ENCAPS crypto_kem_enc_Frodo976
#define KEM_DECAPS crypto_kem_dec_Frodo976
#elif FRODO_LEVEL == 1344
#include "api_frodo1344.h"
#define SCHEME_NAME "FrodoKEM-1344-SHAKE"
#define KEM_KEYPAIR crypto_kem_keypair_Frodo1344
#define KEM_ENCAPS crypto_kem_enc_Frodo1344
#define KEM_DECAPS crypto_kem_dec_Frodo1344
#else
#error Unsupported FRODO_LEVEL
#endif

typedef struct {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
    uint8_t ct_bad[CRYPTO_CIPHERTEXTBYTES];
    uint8_t ss_enc[CRYPTO_BYTES];
    uint8_t ss_dec[CRYPTO_BYTES];
} kem_case;

typedef struct {
    kem_case *cases;
    size_t count;
    size_t next;
    kem_case *last;
} bench_state;

typedef int (*operation)(bench_state *);

typedef struct {
    double minimum_ms;
    double median_ms;
    double p95_ms;
    uint64_t iterations;
} timing;

static volatile uint64_t checksum_sink = 0;

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1.0e9 + (double)ts.tv_nsec;
}

static kem_case *next_case(bench_state *state)
{
    kem_case *item = &state->cases[state->next];
    state->next = (state->next + 1U) % state->count;
    state->last = item;
    return item;
}

static int op_keygen(bench_state *state)
{
    kem_case *item = next_case(state);
    return KEM_KEYPAIR(item->pk, item->sk) == 0;
}

static int op_encaps(bench_state *state)
{
    kem_case *item = next_case(state);
    return KEM_ENCAPS(item->ct, item->ss_enc, item->pk) == 0;
}

static int op_decaps_valid(bench_state *state)
{
    kem_case *item = next_case(state);
    return KEM_DECAPS(item->ss_dec, item->ct, item->sk) == 0;
}

static int op_decaps_invalid(bench_state *state)
{
    kem_case *item = next_case(state);
    return KEM_DECAPS(item->ss_dec, item->ct_bad, item->sk) == 0;
}

static int op_full_cycle(bench_state *state)
{
    kem_case *item = next_case(state);
    return KEM_KEYPAIR(item->pk, item->sk) == 0 &&
           KEM_ENCAPS(item->ct, item->ss_enc, item->pk) == 0 &&
           KEM_DECAPS(item->ss_dec, item->ct, item->sk) == 0;
}

static int compare_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static void run_batch(bench_state *state, operation op, uint64_t iterations)
{
    uint64_t i;
    for (i = 0; i < iterations; ++i) {
        unsigned j;
        if (!op(state)) {
            fprintf(stderr, "%s operation failed.\n", SCHEME_NAME);
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < CRYPTO_BYTES; ++j) {
            checksum_sink = (checksum_sink << 5) ^ (checksum_sink >> 2) ^
                            state->last->ss_enc[j] ^ state->last->ss_dec[j];
        }
        checksum_sink ^= i;
    }
}

static timing measure(bench_state *state, operation op, unsigned samples,
                      double target_ms)
{
    uint64_t iterations = 1;
    double *values = (double *)malloc(samples * sizeof(double));
    double elapsed;
    unsigned i;
    timing result;
    if (values == NULL) exit(EXIT_FAILURE);

    state->next = 0;
    run_batch(state, op, 1);
    do {
        double start;
        state->next = 0;
        start = now_ns();
        run_batch(state, op, iterations);
        elapsed = now_ns() - start;
        if (elapsed < target_ms * 1.0e6) iterations *= 2U;
    } while (elapsed < target_ms * 1.0e6);

    for (i = 0; i < samples; ++i) {
        double start;
        /* Sample i starts from independently prepared case i. */
        state->next = (size_t)i % state->count;
        start = now_ns();
        run_batch(state, op, iterations);
        values[i] = (now_ns() - start) / (double)iterations / 1.0e6;
    }
    qsort(values, samples, sizeof(double), compare_double);
    result.minimum_ms = values[0];
    result.median_ms = values[samples / 2U];
    result.p95_ms = values[((95U * samples + 99U) / 100U) - 1U];
    result.iterations = iterations;
    free(values);
    return result;
}

static int prepare_cases(bench_state *state, int regenerate_keys)
{
    size_t i;
    for (i = 0; i < state->count; ++i) {
        kem_case *item = &state->cases[i];
        if ((regenerate_keys && KEM_KEYPAIR(item->pk, item->sk) != 0) ||
            KEM_ENCAPS(item->ct, item->ss_enc, item->pk) != 0 ||
            KEM_DECAPS(item->ss_dec, item->ct, item->sk) != 0 ||
            memcmp(item->ss_enc, item->ss_dec, CRYPTO_BYTES) != 0) {
            return 0;
        }
        memcpy(item->ct_bad, item->ct, sizeof(item->ct_bad));
        item->ct_bad[0] ^= 1U;
        if (KEM_DECAPS(item->ss_dec, item->ct_bad, item->sk) != 0 ||
            memcmp(item->ss_enc, item->ss_dec, CRYPTO_BYTES) == 0) {
            return 0;
        }
    }
    state->next = 0;
    return 1;
}

static int verify_full_cycles(bench_state *state)
{
    size_t i;
    for (i = 0; i < state->count; ++i) {
        kem_case *item = &state->cases[i];
        if (KEM_KEYPAIR(item->pk, item->sk) != 0 ||
            KEM_ENCAPS(item->ct, item->ss_enc, item->pk) != 0 ||
            KEM_DECAPS(item->ss_dec, item->ct, item->sk) != 0 ||
            memcmp(item->ss_enc, item->ss_dec, CRYPTO_BYTES) != 0) {
            return 0;
        }
    }
    state->next = 0;
    return 1;
}

static void write_row(FILE *file, const char *operation_name,
                      const timing *value, unsigned samples,
                      const bench_state *state)
{
    fprintf(file, "%s,official-scalar-reference,%s,%.6f,%.6f,%.6f,%" PRIu64
                  ",%u,%zu,%u,%u,%u\n",
            SCHEME_NAME, operation_name, value->median_ms, value->p95_ms,
            value->minimum_ms, value->iterations, samples, state->count,
            (unsigned)CRYPTO_PUBLICKEYBYTES,
            (unsigned)CRYPTO_CIPHERTEXTBYTES,
            (unsigned)CRYPTO_SECRETKEYBYTES);
}

int main(int argc, char **argv)
{
    unsigned samples = 31;
    double target_ms = 20.0;
    const char *output_path = "results/frodo.csv";
    timing keygen, encaps, decaps_valid, decaps_invalid, full_cycle;
    bench_state state;
    FILE *output;
    int arg;

    for (arg = 1; arg < argc; ++arg) {
        if (strcmp(argv[arg], "--samples") == 0 && arg + 1 < argc) {
            samples = (unsigned)strtoul(argv[++arg], NULL, 10);
        } else if (strcmp(argv[arg], "--min-ms") == 0 && arg + 1 < argc) {
            target_ms = strtod(argv[++arg], NULL);
        } else if (strcmp(argv[arg], "--output") == 0 && arg + 1 < argc) {
            output_path = argv[++arg];
        } else {
            fprintf(stderr, "Usage: %s [--samples odd] [--min-ms ms] "
                            "[--output path]\n", argv[0]);
            return EXIT_FAILURE;
        }
    }
    if (samples < 5 || (samples & 1U) == 0 || target_ms <= 0.0) {
        fprintf(stderr, "samples must be odd and >= 5; min-ms must be positive.\n");
        return EXIT_FAILURE;
    }

    memset(&state, 0, sizeof(state));
    state.count = samples;
    state.cases = (kem_case *)calloc(state.count, sizeof(*state.cases));
    if (state.cases == NULL) {
        fprintf(stderr, "Could not allocate %u independent %s cases.\n",
                samples, SCHEME_NAME);
        return EXIT_FAILURE;
    }

    bench_random_reset(UINT64_C(0x46726f646f4b454d) ^ FRODO_LEVEL);
    printf("[%s] self-test and %u-case preparation: ", SCHEME_NAME, samples);
    fflush(stdout);
    if (!prepare_cases(&state, 1)) {
        printf("FAILED\n");
        free(state.cases);
        return EXIT_FAILURE;
    }
    printf("PASS\n");

    keygen = measure(&state, op_keygen, samples, target_ms);
    if (!prepare_cases(&state, 1)) return EXIT_FAILURE;
    encaps = measure(&state, op_encaps, samples, target_ms);
    if (!prepare_cases(&state, 0)) return EXIT_FAILURE;
    decaps_valid = measure(&state, op_decaps_valid, samples, target_ms);
    decaps_invalid = measure(&state, op_decaps_invalid, samples, target_ms);
    if (!verify_full_cycles(&state)) return EXIT_FAILURE;
    full_cycle = measure(&state, op_full_cycle, samples, target_ms);

    printf("  KeyGen:          median %.3f ms, p95 %.3f ms\n",
           keygen.median_ms, keygen.p95_ms);
    printf("  Encaps:          median %.3f ms, p95 %.3f ms\n",
           encaps.median_ms, encaps.p95_ms);
    printf("  Decaps (valid):  median %.3f ms, p95 %.3f ms\n",
           decaps_valid.median_ms, decaps_valid.p95_ms);
    printf("  Decaps (invalid):median %.3f ms, p95 %.3f ms\n",
           decaps_invalid.median_ms, decaps_invalid.p95_ms);
    printf("  Full KEM cycle:  median %.3f ms, p95 %.3f ms\n\n",
           full_cycle.median_ms, full_cycle.p95_ms);

    output = fopen(output_path, "w");
    if (output == NULL) {
        perror(output_path);
        free(state.cases);
        return EXIT_FAILURE;
    }
    fprintf(output, "scheme,implementation,operation,median_ms,p95_ms,min_ms,"
                    "iterations_per_sample,samples,case_pool_size,pk_bytes,"
                    "ct_bytes,sk_bytes\n");
    write_row(output, "keygen", &keygen, samples, &state);
    write_row(output, "encaps", &encaps, samples, &state);
    write_row(output, "decaps_valid", &decaps_valid, samples, &state);
    write_row(output, "decaps_invalid", &decaps_invalid, samples, &state);
    write_row(output, "full_cycle", &full_cycle, samples, &state);
    fclose(output);
    free(state.cases);
    printf("%s results written to %s (checksum 0x%016" PRIx64 ")\n",
           SCHEME_NAME, output_path, checksum_sink);
    return EXIT_SUCCESS;
}
