#define _POSIX_C_SOURCE 200809L

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bench_random.h"
#include "leechkem_ref.h"

typedef struct {
    lk_ctx **cases;
    size_t count;
    size_t next;
    lk_ctx *last;
} lk_pool;

typedef int (*operation)(lk_pool *);

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

static lk_ctx *next_case(lk_pool *pool)
{
    lk_ctx *ctx = pool->cases[pool->next];
    pool->next = (pool->next + 1U) % pool->count;
    pool->last = ctx;
    return ctx;
}

static int op_keygen(lk_pool *pool)
{
    return lk_keygen(next_case(pool));
}

static int op_encaps(lk_pool *pool)
{
    return lk_encaps(next_case(pool));
}

static int op_decaps_valid(lk_pool *pool)
{
    return lk_decaps_valid(next_case(pool));
}

static int op_decaps_invalid(lk_pool *pool)
{
    return lk_decaps_invalid(next_case(pool));
}

static int op_full_cycle(lk_pool *pool)
{
    lk_ctx *ctx = next_case(pool);
    return lk_keygen(ctx) && lk_encaps(ctx) && lk_decaps_valid(ctx);
}

static int compare_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static void run_batch(lk_pool *pool, operation op, uint64_t iterations)
{
    uint64_t i;
    for (i = 0; i < iterations; ++i) {
        if (!op(pool)) {
            fprintf(stderr, "A LeechKEM operation failed.\n");
            exit(EXIT_FAILURE);
        }
        checksum_sink ^= lk_checksum(pool->last) + i;
    }
}

static timing measure(lk_pool *pool, operation op, unsigned samples,
                      double target_ms)
{
    uint64_t iterations = 1;
    double *values = (double *)malloc(samples * sizeof(double));
    double elapsed;
    unsigned i;
    timing result;
    if (values == NULL) exit(EXIT_FAILURE);

    pool->next = 0;
    run_batch(pool, op, 1);
    do {
        double start;
        pool->next = 0;
        start = now_ns();
        run_batch(pool, op, iterations);
        elapsed = now_ns() - start;
        if (elapsed < target_ms * 1.0e6) iterations *= 2U;
    } while (elapsed < target_ms * 1.0e6);

    for (i = 0; i < samples; ++i) {
        double start;
        /* Sample i starts from independently prepared case i. */
        pool->next = (size_t)i % pool->count;
        start = now_ns();
        run_batch(pool, op, iterations);
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

static int create_pool(lk_pool *pool, const lk_params *params, size_t count)
{
    size_t i;
    memset(pool, 0, sizeof(*pool));
    pool->cases = (lk_ctx **)calloc(count, sizeof(*pool->cases));
    if (pool->cases == NULL) return 0;
    pool->count = count;
    for (i = 0; i < count; ++i) {
        if (!lk_ctx_create(&pool->cases[i], params)) return 0;
    }
    return 1;
}

static void destroy_pool(lk_pool *pool)
{
    size_t i;
    if (pool == NULL) return;
    for (i = 0; i < pool->count; ++i) lk_ctx_destroy(pool->cases[i]);
    free(pool->cases);
    memset(pool, 0, sizeof(*pool));
}

static int prepare_cases(lk_pool *pool, int regenerate_keys)
{
    size_t i;
    for (i = 0; i < pool->count; ++i) {
        lk_ctx *ctx = pool->cases[i];
        if ((regenerate_keys && !lk_keygen(ctx)) ||
            !lk_encaps(ctx) || !lk_decaps_valid(ctx) ||
            !lk_shared_secret_matches(ctx)) {
            return 0;
        }
        if (!lk_prepare_invalid_ciphertext(ctx) ||
            !lk_decaps_invalid(ctx) || lk_shared_secret_matches(ctx)) {
            return 0;
        }
    }
    pool->next = 0;
    return 1;
}

static int verify_full_cycles(lk_pool *pool)
{
    size_t i;
    for (i = 0; i < pool->count; ++i) {
        lk_ctx *ctx = pool->cases[i];
        if (!lk_keygen(ctx) || !lk_encaps(ctx) || !lk_decaps_valid(ctx) ||
            !lk_shared_secret_matches(ctx)) {
            return 0;
        }
    }
    pool->next = 0;
    return 1;
}

static void write_header(FILE *file)
{
    fprintf(file, "scheme,implementation,operation,median_ms,p95_ms,min_ms,"
                  "iterations_per_sample,samples,case_pool_size,pk_bytes,"
                  "ct_bytes,sk_bytes\n");
}

static void write_row(FILE *file, const lk_params *p, const char *operation_name,
                      const timing *value, unsigned samples, const lk_pool *pool)
{
    const lk_ctx *ctx = pool->cases[0];
    fprintf(file, "%s,scalar-reference,%s,%.6f,%.6f,%.6f,%" PRIu64
                  ",%u,%zu,%zu,%zu,%zu\n",
            p->name, operation_name, value->median_ms, value->p95_ms,
            value->minimum_ms, value->iterations, samples, pool->count,
            lk_public_key_bytes(ctx), lk_ciphertext_bytes(ctx),
            lk_secret_key_bytes_estimate(ctx));
}

int main(int argc, char **argv)
{
    unsigned samples = 31;
    double target_ms = 20.0;
    const char *output_path = "results/leech_full.csv";
    FILE *output;
    size_t set;
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
    output = fopen(output_path, "w");
    if (output == NULL) {
        perror(output_path);
        return EXIT_FAILURE;
    }
    write_header(output);
    bench_random_reset(UINT64_C(0x4c656563684b454d));

    for (set = 0; set < lk_parameter_count(); ++set) {
        const lk_params *p = lk_parameter_at(set);
        lk_pool pool;
        timing keygen, encaps, decaps_valid, decaps_invalid, full_cycle;
        if (!create_pool(&pool, p, samples)) {
            fprintf(stderr, "Could not allocate %u independent %s cases.\n",
                    samples, p->name);
            destroy_pool(&pool);
            fclose(output);
            return EXIT_FAILURE;
        }
        printf("[%s] self-test and %u-case preparation: ", p->name, samples);
        fflush(stdout);
        if (!lk_selftest(pool.cases[0]) || !prepare_cases(&pool, 1)) {
            printf("FAILED\n");
            destroy_pool(&pool);
            fclose(output);
            return EXIT_FAILURE;
        }
        printf("PASS\n");

        keygen = measure(&pool, op_keygen, samples, target_ms);
        if (!prepare_cases(&pool, 1)) return EXIT_FAILURE;
        encaps = measure(&pool, op_encaps, samples, target_ms);
        if (!prepare_cases(&pool, 0)) return EXIT_FAILURE;
        decaps_valid = measure(&pool, op_decaps_valid, samples, target_ms);
        decaps_invalid = measure(&pool, op_decaps_invalid, samples, target_ms);
        if (!verify_full_cycles(&pool)) return EXIT_FAILURE;
        full_cycle = measure(&pool, op_full_cycle, samples, target_ms);

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

        write_row(output, p, "keygen", &keygen, samples, &pool);
        write_row(output, p, "encaps", &encaps, samples, &pool);
        write_row(output, p, "decaps_valid", &decaps_valid, samples, &pool);
        write_row(output, p, "decaps_invalid", &decaps_invalid, samples, &pool);
        write_row(output, p, "full_cycle", &full_cycle, samples, &pool);
        destroy_pool(&pool);
    }
    fclose(output);
    printf("LeechKEM results written to %s (checksum 0x%016" PRIx64 ")\n",
           output_path, checksum_sink);
    return EXIT_SUCCESS;
}
